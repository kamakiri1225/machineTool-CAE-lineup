"""ref005 の図: 2D トポロジー最適化（SIMP 法、コンプライアンス最小化）の例
   Andreassen et al. "Efficient topology optimization in MATLAB using 88 lines of code" (2011) を Python に移植・簡略化。
   題材: 底面固定のコラムの上端に、横方向の切削反力がかかる（材料を 35 % だけ使って最も硬い形を探す）"""
import os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from figstyle import plt, RED, BLUE

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref005")
os.makedirs(OUT, exist_ok=True)


def lk(nu=0.3):
    k = np.array([1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8])
    return 1/(1-nu**2)*np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]], [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]], [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]], [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]], [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])


def top(nelx, nely, volfrac, fixed, F, penal=3.0, rmin=1.5, iters=80, snap=()):
    E0, Emin = 1.0, 1e-9
    ndof = 2 * (nelx + 1) * (nely + 1)
    KE = lk()
    edof = np.zeros((nelx * nely, 8), dtype=int)
    for elx in range(nelx):
        for ely in range(nely):
            el = ely + elx * nely
            n1 = (nely + 1) * elx + ely; n2 = (nely + 1) * (elx + 1) + ely
            edof[el] = [2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3, 2*n2, 2*n2+1, 2*n1, 2*n1+1]
    iK = np.kron(edof, np.ones((8, 1))).flatten(); jK = np.kron(edof, np.ones((1, 8))).flatten()
    # フィルタ
    iH, jH, sH = [], [], []
    for i in range(nelx):
        for j in range(nely):
            r = i * nely + j
            for k in range(max(i - int(rmin) + 1, 0), min(i + int(rmin), nelx)):
                for l in range(max(j - int(rmin) + 1, 0), min(j + int(rmin), nely)):
                    iH.append(r); jH.append(k * nely + l); sH.append(max(0, rmin - np.hypot(i - k, j - l)))
    H = coo_matrix((sH, (iH, jH)), shape=(nelx*nely, nelx*nely)).tocsc(); Hs = np.array(H.sum(1)).flatten()
    free = np.setdiff1d(np.arange(ndof), fixed)
    x = volfrac * np.ones(nelx * nely); xPhys = x.copy()
    shots = {}
    for it in range(1, iters + 1):
        sK = ((KE.flatten()[np.newaxis]).T * (Emin + xPhys**penal * (E0 - Emin))).flatten(order="F")
        K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
        U = np.zeros(ndof); U[free] = spsolve(K[free][:, free], F[free])
        ce = (np.dot(U[edof], KE) * U[edof]).sum(1)
        dc = -penal * xPhys**(penal - 1) * (E0 - Emin) * ce
        dc = np.asarray(H * (x * dc)[np.newaxis].T / Hs[np.newaxis].T)[:, 0] / np.maximum(1e-3, x)
        l1, l2, move = 0, 1e9, 0.2
        while (l2 - l1) / (l1 + l2) > 1e-3:
            lmid = 0.5 * (l2 + l1)
            xnew = np.maximum(0, np.maximum(x - move, np.minimum(1, np.minimum(x + move, x * np.sqrt(-dc / lmid)))))
            if xnew.sum() > volfrac * nelx * nely: l1 = lmid
            else: l2 = lmid
        x = xnew; xPhys = x
        if it in snap:
            shots[it] = xPhys.reshape((nelx, nely)).T.copy()
    return xPhys.reshape((nelx, nely)).T, shots


if __name__ == "__main__":
    nelx, nely = 60, 120          # 幅 60 × 高さ 120 要素（コラムの縦断面）
    ndof = 2 * (nelx + 1) * (nely + 1)
    nid = lambda ix, iy: (nely + 1) * ix + iy           # iy=0 が上端
    fixed = np.array([d for ix in range(nelx + 1) for d in (2 * nid(ix, nely), 2 * nid(ix, nely) + 1)])
    F = np.zeros(ndof)
    for iy in range(0, 8):                               # 上端の前面（左側）に横向きの力
        F[2 * nid(0, iy)] = 1.0
    xf, shots = top(nelx, nely, 0.35, fixed, F, rmin=3.0, iters=80, snap=(1, 10, 30))
    fig, axs = plt.subplots(1, 4, figsize=(15, 7.2))
    for ax, (k, im) in zip(axs, list(shots.items()) + [(80, xf)]):
        ax.imshow(1 - im, cmap="gray", vmin=0, vmax=1, origin="upper"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("1 回目（灰色 = 中間の密度）" if k == 1 else f"{k} 回目" + ("（最終）" if k == 80 else ""), fontsize=12)
    ax = axs[0]
    ax.annotate("", xy=(0, 4), xytext=(-14, 4), arrowprops=dict(arrowstyle="->", color=RED, lw=3), annotation_clip=False)
    ax.text(-16, -4, "切削反力", color=RED, fontsize=11, ha="left")
    for a in axs:
        a.plot([0, nelx - 1], [nely - 0.5, nely - 0.5], color=BLUE, lw=6)
    axs[0].text(2, nely + 6, "底面を固定（ベッドに取付）", color=BLUE, fontsize=11)
    fig.suptitle("2D トポロジー最適化の例：コラムの縦断面で、材料 35 % で最も硬い形を探す（SIMP 法）", fontsize=13, y=1.03)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "topopt_column_2d.png")); plt.close(fig)
    print("saved")
