import SwiftUI
import Accelerate

// MARK: - 0. LAPACK 直接バインディング (非推奨エラーの回避)
// @_silgen_name により、ヘッダーの非推奨警告をバイパスしてネイティブの LAPACK 関数を直接呼び出します
@_silgen_name("dsygv_")
private func lapack_dsygv(
    _ itype: UnsafeMutablePointer<Int32>,
    _ jobz: UnsafeMutablePointer<CChar>,
    _ uplo: UnsafeMutablePointer<CChar>,
    _ n: UnsafeMutablePointer<Int32>,
    _ a: UnsafeMutablePointer<Double>,
    _ lda: UnsafeMutablePointer<Int32>,
    _ b: UnsafeMutablePointer<Double>,
    _ ldb: UnsafeMutablePointer<Int32>,
    _ w: UnsafeMutablePointer<Double>,
    _ work: UnsafeMutablePointer<Double>,
    _ lwork: UnsafeMutablePointer<Int32>,
    _ info: UnsafeMutablePointer<Int32>
) -> Int32

// MARK: - 1. データ型・モデル定義

enum BoundaryCondition: String, CaseIterable, Identifiable {
    case pinned = "ピン支持 (Pinned)"
    case fixed = "固定 (Fixed)"
    var id: String { self.rawValue }
}

struct FEMParameters {
    // 形状情報
    var lx: Double = 10.0
    var ly: Double = 10.0
    var h: Double = 0.3
    
    // 要素分割数
    var nx: Int = 30
    var ny: Int = 30
    
    // 境界条件
    var bcBottom: BoundaryCondition = .pinned
    var bcTop: BoundaryCondition = .pinned
    var bcLeft: BoundaryCondition = .pinned
    var bcRight: BoundaryCondition = .pinned
    
    // 物性値 (ヤング率は指数入力対応のため文字列として保持)
    var eString: String = "2.5e10"
    var nu: Double = 0.167
    var gamma: Double = 24000.0
    var kappa: Double = 0.833
    
    var E: Double {
        Double(eString) ?? 2.5e10
    }
}

struct ModeResult: Identifiable {
    let id = UUID()
    let modeNumber: Int
    let freq: Double
    let period: Double
    let beta: Double
    let shape: [[Double]] // ny+1 × nx+1 の変位wグリッド (最大振幅1.0)
}

// MARK: - 2. 有限要素法 (FEM) ソルバー

class FEMPlateSolver {
    static func runAnalysis(params: FEMParameters) -> [ModeResult] {
        let lx = params.lx
        let ly = params.ly
        let h = params.h
        let nx = params.nx
        let ny = params.ny
        let E = params.E
        let nu = params.nu
        let gamma = params.gamma
        let kappa = params.kappa
        let g = 9.80665
        let rho = gamma / g
        let G = E / (2.0 * (1.0 + nu))
        
        let nNodesX = nx + 1
        let nNodesY = ny + 1
        let nNodes = nNodesX * nNodesY
        let nElements = nx * ny
        let totalDof = 3 * nNodes
        
        // 節点座標
        var coords = [(x: Double, y: Double)](repeating: (0, 0), count: nNodes)
        for j in 0..<nNodesY {
            for i in 0..<nNodesX {
                let nid = j * nNodesX + i
                coords[nid] = (Double(i) * (lx / Double(nx)), Double(j) * (ly / Double(ny)))
            }
        }
        
        // 要素コネクティビティ
        var elements = [[Int]](repeating: [], count: nElements)
        for j in 0..<ny {
            for i in 0..<nx {
                let eid = j * nx + i
                let n1 = j * nNodesX + i
                let n2 = n1 + 1
                let n3 = (j + 1) * nNodesX + (i + 1)
                let n4 = (j + 1) * nNodesX + i
                elements[eid] = [n1, n2, n3, n4]
            }
        }
        
        // 材料剛性マトリクス
        let Dconst = (E * pow(h, 3)) / (12.0 * (1.0 - pow(nu, 2)))
        let Db: [[Double]] = [
            [Dconst, Dconst * nu, 0.0],
            [Dconst * nu, Dconst, 0.0],
            [0.0, 0.0, Dconst * (1.0 - nu) / 2.0]
        ]
        let Ds: [[Double]] = [
            [kappa * G * h, 0.0],
            [0.0, kappa * G * h]
        ]
        
        // 要素剛性・質量計算用サブルーチン
        func getElementMatrices(xElem: [Double], yElem: [Double]) -> (ke: [[Double]], me: [[Double]]) {
            var ke = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 12)
            var me = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 12)
            
            let gp: [Double] = [-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)]
            let w: [Double] = [1.0, 1.0]
            
            // 1) 曲げ剛性 Ke_b (2x2 ガウス積分)
            for (xi, wxi) in zip(gp, w) {
                for (eta, weta) in zip(gp, w) {
                    let dN_dxi = [-0.25 * (1 - eta),  0.25 * (1 - eta),  0.25 * (1 + eta), -0.25 * (1 + eta)]
                    let dN_deta = [-0.25 * (1 - xi), -0.25 * (1 + xi),  0.25 * (1 + xi),  0.25 * (1 - xi)]
                    
                    var J00 = 0.0, J01 = 0.0, J10 = 0.0, J11 = 0.0
                    for a in 0..<4 {
                        J00 += dN_dxi[a] * xElem[a];  J01 += dN_dxi[a] * yElem[a]
                        J10 += dN_deta[a] * xElem[a]; J11 += dN_deta[a] * yElem[a]
                    }
                    let detJ = J00 * J11 - J01 * J10
                    let invJ00 =  J11 / detJ; let invJ01 = -J01 / detJ
                    let invJ10 = -J10 / detJ; let invJ11 =  J00 / detJ
                    
                    var dN_dx = [Double](repeating: 0, count: 4)
                    var dN_dy = [Double](repeating: 0, count: 4)
                    for a in 0..<4 {
                        dN_dx[a] = invJ00 * dN_dxi[a] + invJ01 * dN_deta[a]
                        dN_dy[a] = invJ10 * dN_dxi[a] + invJ11 * dN_deta[a]
                    }
                    
                    var Bb = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 3)
                    for a in 0..<4 {
                        Bb[0][3 * a + 1] = dN_dx[a]
                        Bb[1][3 * a + 2] = dN_dy[a]
                        Bb[2][3 * a + 1] = dN_dy[a]
                        Bb[2][3 * a + 2] = dN_dx[a]
                    }
                    
                    let factor = detJ * wxi * weta
                    var temp = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 3)
                    for r in 0..<3 {
                        for c in 0..<12 {
                            for k in 0..<3 { temp[r][c] += Db[r][k] * Bb[k][c] }
                        }
                    }
                    for r in 0..<12 {
                        for c in 0..<12 {
                            var sum = 0.0
                            for k in 0..<3 { sum += Bb[k][r] * temp[k][c] }
                            ke[r][c] += sum * factor
                        }
                    }
                }
            }
            
            // 2) せん断剛性 Ke_s (1x1 ガウス積分)
            let dN_dxi_s = [-0.25,  0.25,  0.25, -0.25]
            let dN_deta_s = [-0.25, -0.25,  0.25,  0.25]
            var J00_s = 0.0, J01_s = 0.0, J10_s = 0.0, J11_s = 0.0
            for a in 0..<4 {
                J00_s += dN_dxi_s[a] * xElem[a];  J01_s += dN_dxi_s[a] * yElem[a]
                J10_s += dN_deta_s[a] * xElem[a]; J11_s += dN_deta_s[a] * yElem[a]
            }
            let detJ_s = J00_s * J11_s - J01_s * J10_s
            let invJ00_s =  J11_s / detJ_s; let invJ01_s = -J01_s / detJ_s
            let invJ10_s = -J10_s / detJ_s; let invJ11_s =  J00_s / detJ_s
            
            var Bs = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 2)
            let Ns = [0.25, 0.25, 0.25, 0.25]
            for a in 0..<4 {
                let dN_dx = invJ00_s * dN_dxi_s[a] + invJ01_s * dN_deta_s[a]
                let dN_dy = invJ10_s * dN_dxi_s[a] + invJ11_s * dN_deta_s[a]
                Bs[0][3 * a + 0] = dN_dx
                Bs[0][3 * a + 1] = -Ns[a]
                Bs[1][3 * a + 0] = dN_dy
                Bs[1][3 * a + 2] = -Ns[a]
            }
            
            let factor_s = detJ_s * 4.0
            var temp_s = [[Double]](repeating: [Double](repeating: 0.0, count: 12), count: 2)
            for r in 0..<2 {
                for c in 0..<12 {
                    for k in 0..<2 { temp_s[r][c] += Ds[r][k] * Bs[k][c] }
                }
            }
            for r in 0..<12 {
                for c in 0..<12 {
                    var sum = 0.0
                    for k in 0..<2 { sum += Bs[k][r] * temp_s[k][c] }
                    ke[r][c] += sum * factor_s
                }
            }
            
            // 3) 質量マトリクス Me (2x2 ガウス積分)
            let Irot = rho * pow(h, 3) / 12.0
            let mtrans = rho * h
            for (xi, wxi) in zip(gp, w) {
                for (eta, weta) in zip(gp, w) {
                    let N = [
                        0.25 * (1 - xi) * (1 - eta),
                        0.25 * (1 + xi) * (1 - eta),
                        0.25 * (1 + xi) * (1 + eta),
                        0.25 * (1 - xi) * (1 + eta)
                    ]
                    let dN_dxi = [-0.25 * (1 - eta),  0.25 * (1 - eta),  0.25 * (1 + eta), -0.25 * (1 + eta)]
                    let dN_deta = [-0.25 * (1 - xi), -0.25 * (1 + xi),  0.25 * (1 + xi),  0.25 * (1 - xi)]
                    var J00 = 0.0, J01 = 0.0, J10 = 0.0, J11 = 0.0
                    for a in 0..<4 {
                        J00 += dN_dxi[a] * xElem[a];  J01 += dN_dxi[a] * yElem[a]
                        J10 += dN_deta[a] * xElem[a]; J11 += dN_deta[a] * yElem[a]
                    }
                    let detJ = J00 * J11 - J01 * J10
                    let factor = detJ * wxi * weta
                    
                    for a in 0..<4 {
                        for b in 0..<4 {
                            let valM = mtrans * N[a] * N[b] * factor
                            let valI = Irot * N[a] * N[b] * factor
                            me[3 * a + 0][3 * b + 0] += valM
                            me[3 * a + 1][3 * b + 1] += valI
                            me[3 * a + 2][3 * b + 2] += valI
                        }
                    }
                }
            }
            return (ke, me)
        }
        
        // 全体マトリクスの組み立て
        var K_global = [Double](repeating: 0.0, count: totalDof * totalDof)
        var M_global = [Double](repeating: 0.0, count: totalDof * totalDof)
        
        for elem in elements {
            let xe = elem.map { coords[$0].x }
            let ye = elem.map { coords[$0].y }
            let (ke, me) = getElementMatrices(xElem: xe, yElem: ye)
            
            var edof: [Int] = []
            for n in elem {
                edof.append(contentsOf: [3 * n, 3 * n + 1, 3 * n + 2])
            }
            
            for r in 0..<12 {
                let gr = edof[r]
                for c in 0..<12 {
                    let gc = edof[c]
                    K_global[gr * totalDof + gc] += ke[r][c]
                    M_global[gr * totalDof + gc] += me[r][c]
                }
            }
        }
        
        // 境界条件の適用
        var fixedDofs = Set<Int>()
        let tol = 1e-6
        for (i, coord) in coords.enumerated() {
            // Bottom (y=0)
            if abs(coord.y) < tol {
                if params.bcBottom == .pinned {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1)
                } else {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1); fixedDofs.insert(3 * i + 2)
                }
            }
            // Top (y=Ly)
            if abs(coord.y - ly) < tol {
                if params.bcTop == .pinned {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1)
                } else {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1); fixedDofs.insert(3 * i + 2)
                }
            }
            // Left (x=0)
            if abs(coord.x) < tol {
                if params.bcLeft == .pinned {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 2)
                } else {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1); fixedDofs.insert(3 * i + 2)
                }
            }
            // Right (x=Lx)
            if abs(coord.x - lx) < tol {
                if params.bcRight == .pinned {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 2)
                } else {
                    fixedDofs.insert(3 * i); fixedDofs.insert(3 * i + 1); fixedDofs.insert(3 * i + 2)
                }
            }
        }
        
        var activeDofs: [Int] = []
        for d in 0..<totalDof {
            if !fixedDofs.contains(d) {
                activeDofs.append(d)
            }
        }
        
        let nAct = activeDofs.count
        guard nAct > 5 else { return [] }
        
        // アクティブ部分マトリクスの抽出 (Column-major 配列)
        var K_act = [Double](repeating: 0.0, count: nAct * nAct)
        var M_act = [Double](repeating: 0.0, count: nAct * nAct)
        
        for c in 0..<nAct {
            let gc = activeDofs[c]
            for r in 0..<nAct {
                let gr = activeDofs[r]
                K_act[c * nAct + r] = K_global[gr * totalDof + gc]
                M_act[c * nAct + r] = M_global[gr * totalDof + gc]
            }
        }
        
        // 加振ベクトル (w成分に1.0)
        var r_act = [Double](repeating: 0.0, count: nAct)
        for i in 0..<nAct {
            if activeDofs[i] % 3 == 0 {
                r_act[i] = 1.0
            }
        }
        
        // LAPACK による一般化対称固有値解析 (K * x = lambda * M * x)
        var itype: Int32 = 1
        var jobz: CChar = 86 // 'V': 固有ベクトルを計算
        var uplo: CChar = 85 // 'U': 上三角部を使用
        var n: Int32 = Int32(nAct)
        var lda: Int32 = n
        var ldb: Int32 = n
        var w_eigen = [Double](repeating: 0.0, count: nAct)
        var lwork: Int32 = -1
        var workQuery = [Double](repeating: 0.0, count: 1)
        var info: Int32 = 0
        
        // ワークスペースサイズのクエリ
        _ = lapack_dsygv(&itype, &jobz, &uplo, &n, &K_act, &lda, &M_act, &ldb, &w_eigen, &workQuery, &lwork, &info)
        lwork = Int32(workQuery[0])
        var work = [Double](repeating: 0.0, count: Int(lwork))
        
        // 固有値計算の実行
        _ = lapack_dsygv(&itype, &jobz, &uplo, &n, &K_act, &lda, &M_act, &ldb, &w_eigen, &work, &lwork, &info)
        
        guard info == 0 else {
            print("LAPACK dsygv error: info = \(info)")
            return []
        }
        
        // 元の M_act を再作成（刺激係数計算で利用するため）
        var M_orig = [Double](repeating: 0.0, count: nAct * nAct)
        for c in 0..<nAct {
            let gc = activeDofs[c]
            for r in 0..<nAct {
                let gr = activeDofs[r]
                M_orig[c * nAct + r] = M_global[gr * totalDof + gc]
            }
        }
        
        var results: [ModeResult] = []
        
        for m in 0..<min(5, nAct) {
            let omega2 = w_eigen[m]
            let omega = sqrt(max(0.0, omega2))
            let freq = omega / (2.0 * Double.pi)
            let period = freq > 1e-6 ? (1.0 / freq) : 0.0
            
            // 固有ベクトル抽出 (K_act に列優先で格納されている)
            var phi = [Double](repeating: 0.0, count: nAct)
            for i in 0..<nAct {
                phi[i] = K_act[m * nAct + i]
            }
            
            // M @ phi
            var M_phi = [Double](repeating: 0.0, count: nAct)
            vDSP_mmulD(M_orig, 1, phi, 1, &M_phi, 1, vDSP_Length(nAct), 1, vDSP_Length(nAct))
            
            // モード質量 phi^T @ M @ phi
            var modalMass = 0.0
            vDSP_dotprD(phi, 1, M_phi, 1, &modalMass, vDSP_Length(nAct))
            
            let normFactor = sqrt(modalMass)
            if normFactor > 1e-12 {
                for i in 0..<nAct {
                    phi[i] /= normFactor
                    M_phi[i] /= normFactor
                }
            }
            
            // 刺激係数 beta = phi^T @ M @ r
            var beta = 0.0
            vDSP_dotprD(M_phi, 1, r_act, 1, &beta, vDSP_Length(nAct))
            
            // 全体自由度へのマッピングと変位成分 w の抽出
            var phiFull = [Double](repeating: 0.0, count: totalDof)
            for (idx, dof) in activeDofs.enumerated() {
                phiFull[dof] = phi[idx]
            }
            
            // 2次元グリッド (ny+1, nx+1) への変換
            var shape = [[Double]](repeating: [Double](repeating: 0.0, count: nNodesX), count: nNodesY)
            var maxAmp = 0.0
            for j in 0..<nNodesY {
                for i in 0..<nNodesX {
                    let nid = j * nNodesX + i
                    let wVal = phiFull[3 * nid]
                    shape[j][i] = wVal
                    if abs(wVal) > maxAmp { maxAmp = abs(wVal) }
                }
            }
            
            // 最大振幅を 1.0 にスケーリング
            if maxAmp > 1e-9 {
                for j in 0..<nNodesY {
                    for i in 0..<nNodesX {
                        shape[j][i] /= maxAmp
                    }
                }
            }
            
            results.append(ModeResult(
                modeNumber: m + 1,
                freq: freq,
                period: period,
                beta: beta,
                shape: shape
            ))
        }
        
        return results
    }
}

// MARK: - 3. 描画コンポーネント (メッシュ図 & モード図)

// 要素分割図
struct MeshPlotView: View {
    let nx: Int
    let ny: Int
    let lx: Double
    let ly: Double
    
    var body: some View {
        Canvas { context, size in
            let padding: CGFloat = 20
            let w = size.width - 2 * padding
            let h = size.height - 2 * padding
            
            let scaleX = w / CGFloat(lx)
            let scaleY = h / CGFloat(ly)
            let scale = min(scaleX, scaleY)
            
            let drawW = CGFloat(lx) * scale
            let drawH = CGFloat(ly) * scale
            let originX = padding + (w - drawW) / 2
            let originY = padding + (h - drawH) / 2
            
            // グリッド線
            for i in 0...nx {
                let x = originX + CGFloat(i) * (drawW / CGFloat(nx))
                var path = Path()
                path.move(to: CGPoint(x: x, y: originY))
                path.addLine(to: CGPoint(x: x, y: originY + drawH))
                context.stroke(path, with: .color(.blue.opacity(0.6)), lineWidth: 1)
            }
            for j in 0...ny {
                let y = originY + CGFloat(j) * (drawH / CGFloat(ny))
                var path = Path()
                path.move(to: CGPoint(x: originX, y: y))
                path.addLine(to: CGPoint(x: originX + drawW, y: y))
                context.stroke(path, with: .color(.blue.opacity(0.6)), lineWidth: 1)
            }
            
            // 外枠
            let rect = CGRect(x: originX, y: originY, width: drawW, height: drawH)
            context.stroke(Path(rect), with: .color(.primary), lineWidth: 2)
        }
        .frame(height: 220)
        .background(Color(UIColor.secondarySystemBackground))
        .cornerRadius(12)
    }
}

// モード図 (2.5D等角投影サーフェス)
struct ModeShapePlotView: View {
    let shape: [[Double]]
    
    var body: some View {
        Canvas { context, size in
            let ny = shape.count - 1
            let nx = shape[0].count - 1
            guard ny > 0 && nx > 0 else { return }
            
            let center = CGPoint(x: size.width * 0.5, y: size.height * 0.55)
            let span = min(size.width, size.height) * 0.38
            
            func project(i: Int, j: Int, w: Double) -> CGPoint {
                let u = (Double(i) / Double(nx) - 0.5) * 2.0
                let v = (Double(j) / Double(ny) - 0.5) * 2.0
                let xIso = (u - v) * cos(.pi / 6) * span
                let yIso = (u + v) * sin(.pi / 6) * span - (w * span * 0.5)
                return CGPoint(x: center.x + CGFloat(xIso), y: center.y + CGFloat(yIso))
            }
            
            func getColor(val: Double) -> Color {
                let t = CGFloat((val + 1.0) / 2.0)
                return Color(
                    red: Double(min(max(2.0 * t - 0.5, 0.0), 1.0)),
                    green: Double(1.0 - abs(2.0 * t - 1.0)),
                    blue: Double(min(max(2.0 * (1.0 - t) - 0.5, 0.0), 1.0))
                )
            }
            
            for j in 0..<ny {
                for i in 0..<nx {
                    let p1 = project(i: i, j: j, w: shape[j][i])
                    let p2 = project(i: i + 1, j: j, w: shape[j][i + 1])
                    let p3 = project(i: i + 1, j: j + 1, w: shape[j + 1][i + 1])
                    let p4 = project(i: i, j: j + 1, w: shape[j + 1][i])
                    
                    var poly = Path()
                    poly.move(to: p1)
                    poly.addLine(to: p2)
                    poly.addLine(to: p3)
                    poly.addLine(to: p4)
                    poly.closeSubpath()
                    
                    let avgVal = (shape[j][i] + shape[j][i + 1] + shape[j + 1][i + 1] + shape[j + 1][i]) / 4.0
                    context.fill(poly, with: .color(getColor(val: avgVal).opacity(0.85)))
                    context.stroke(poly, with: .color(.black.opacity(0.3)), lineWidth: 0.5)
                }
            }
        }
        .frame(height: 320)
        .background(Color(UIColor.secondarySystemBackground))
        .cornerRadius(12)
    }
}

// MARK: - 4. 画面ビュー

// 画面１：形状・要素分割数 入力画面 (ContentView)
struct ContentView: View {
    @State private var params = FEMParameters()
    
    var body: some View {
        NavigationStack {
            Form {
                Section(header: Text("形状情報")) {
                    HStack {
                        Text("横の長さ [m]")
                        Spacer()
                        TextField("横の長さ", value: $params.lx, format: .number)
                            .keyboardType(.decimalPad)
                            .multilineTextAlignment(.trailing)
                    }
                    HStack {
                        Text("縦の長さ [m]")
                        Spacer()
                        TextField("縦の長さ", value: $params.ly, format: .number)
                            .keyboardType(.decimalPad)
                            .multilineTextAlignment(.trailing)
                    }
                    HStack {
                        Text("平板の厚さ [m]")
                        Spacer()
                        TextField("厚さ", value: $params.h, format: .number)
                            .keyboardType(.decimalPad)
                            .multilineTextAlignment(.trailing)
                    }
                }
                
                Section(header: Text("要素分割の情報")) {
                    HStack {
                        Text("横の要素分割数 (nx)")
                        Spacer()
                        TextField("nx", value: $params.nx, format: .number)
                            .keyboardType(.numberPad)
                            .multilineTextAlignment(.trailing)
                    }
                    HStack {
                        Text("縦の要素分割数 (ny)")
                        Spacer()
                        TextField("ny", value: $params.ny, format: .number)
                            .keyboardType(.numberPad)
                            .multilineTextAlignment(.trailing)
                    }
                }
                
                Section {
                    NavigationLink(destination: MeshAndMaterialView(params: params)) {
                        Text("メッシュ確認・物性設定へ進む")
                            .frame(maxWidth: .infinity, alignment: .center)
                            .foregroundColor(.blue)
                    }
                }
            }
            .navigationTitle("FEM 平板解析 (1/3)")
        }
    }
}

// 画面２：メッシュ確認 ＆ 境界条件・物性値入力画面
struct MeshAndMaterialView: View {
    @State var params: FEMParameters
    @State private var isCalculating = false
    @State private var results: [ModeResult] = []
    @State private var navigateToResults = false
    
    var body: some View {
        Form {
            Section(header: Text("要素分割図")) {
                MeshPlotView(nx: params.nx, ny: params.ny, lx: params.lx, ly: params.ly)
                    .padding(.vertical, 4)
            }
            
            Section(header: Text("4辺の境界条件")) {
                Picker("下端 (y = 0)", selection: $params.bcBottom) {
                    ForEach(BoundaryCondition.allCases) { Text($0.rawValue).tag($0) }
                }
                Picker("上端 (y = Ly)", selection: $params.bcTop) {
                    ForEach(BoundaryCondition.allCases) { Text($0.rawValue).tag($0) }
                }
                Picker("左端 (x = 0)", selection: $params.bcLeft) {
                    ForEach(BoundaryCondition.allCases) { Text($0.rawValue).tag($0) }
                }
                Picker("右端 (x = Lx)", selection: $params.bcRight) {
                    ForEach(BoundaryCondition.allCases) { Text($0.rawValue).tag($0) }
                }
            }
            
            Section(header: Text("物性の入力項目")) {
                HStack {
                    Text("ヤング率 [N/m²]")
                    Spacer()
                    TextField("例: 2.5e10", text: $params.eString)
                        .keyboardType(.numbersAndPunctuation)
                        .multilineTextAlignment(.trailing)
                }
                HStack {
                    Text("ポアソン比")
                    Spacer()
                    TextField("ν", value: $params.nu, format: .number)
                        .keyboardType(.decimalPad)
                        .multilineTextAlignment(.trailing)
                }
                HStack {
                    Text("単位体積重量 [N/m³]")
                    Spacer()
                    TextField("γ", value: $params.gamma, format: .number)
                        .keyboardType(.decimalPad)
                        .multilineTextAlignment(.trailing)
                }
                HStack {
                    Text("せん断形状係数")
                    Spacer()
                    TextField("κ", value: $params.kappa, format: .number)
                        .keyboardType(.decimalPad)
                        .multilineTextAlignment(.trailing)
                }
            }
            
            Section {
                Button(action: runCalculation) {
                    if isCalculating {
                        HStack {
                            Spacer()
                            ProgressView()
                            Text(" 計算中...")
                            Spacer()
                        }
                    } else {
                        Text("固有値解析を実行")
                            .frame(maxWidth: .infinity, alignment: .center)
                            .foregroundColor(.blue)
                    }
                }
                .disabled(isCalculating)
            }
        }
        .navigationTitle("メッシュ・条件設定 (2/3)")
        .navigationDestination(isPresented: $navigateToResults) {
            ModeResultsSwipeView(results: results)
        }
    }
    
    private func runCalculation() {
        isCalculating = true
        DispatchQueue.global(qos: .userInitiated).async {
            let res = FEMPlateSolver.runAnalysis(params: params)
            DispatchQueue.main.async {
                self.results = res
                self.isCalculating = false
                if !res.isEmpty {
                    self.navigateToResults = true
                }
            }
        }
    }
}

// 画面３：モード図＆固有値出力画面（横スワイプ表示）
struct ModeResultsSwipeView: View {
    let results: [ModeResult]
    @State private var selectedIndex = 0
    
    var body: some View {
        VStack {
            TabView(selection: $selectedIndex) {
                ForEach(Array(results.enumerated()), id: \.offset) { index, mode in
                    VStack(spacing: 20) {
                        Text("Mode \(mode.modeNumber)")
                            .font(.title)
                            .bold()
                        
                        // 独立したモード図
                        ModeShapePlotView(shape: mode.shape)
                            .padding(.horizontal)
                        
                        // 固有振動数、固有周期、刺激係数
                        VStack(spacing: 10) {
                            HStack {
                                Text("固有振動数 (f):")
                                    .bold()
                                Spacer()
                                Text(String(format: "%.3f Hz", mode.freq))
                            }
                            Divider()
                            HStack {
                                Text("固有周期 (T):")
                                    .bold()
                                Spacer()
                                Text(String(format: "%.4f s", mode.period))
                            }
                            Divider()
                            HStack {
                                Text("刺激係数 (β):")
                                    .bold()
                                Spacer()
                                Text(String(format: "%.4e", mode.beta))
                            }
                        }
                        .padding()
                        .background(Color(UIColor.secondarySystemBackground))
                        .cornerRadius(12)
                        .padding(.horizontal)
                        
                        Spacer()
                    }
                    .tag(index)
                }
            }
            .tabViewStyle(.page(indexDisplayMode: .always))
            .indexViewStyle(.page(backgroundDisplayMode: .always))
        }
        .navigationTitle("解析結果 (3/3)")
        .navigationBarTitleDisplayMode(.inline)
    }
}
