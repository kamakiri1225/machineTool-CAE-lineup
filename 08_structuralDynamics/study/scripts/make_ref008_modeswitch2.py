"""ref008 の図: 1 次固有振動数を最大化するトポロジー最適化で、モードの順番が入れ替わる様子（MAC で追跡）
   コラム断面（make_ref005_dyn と同じ）、材料 35 %、一様密度から開始。目的 = 1 次の固有振動数 ω1 を最大（単純な OC 法、追跡なし）
   比較：同じ問題で「最初の 1 次と同じ形のモード」を MAC で追いかけて最大化する（Kim・Kim 2000 の考え方）"""
import os
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import make_ref005_dyn as M
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); NM = 4


def eig(xP):
    K = M.stiff(xP)[M.free][:, M.free]; Md = M.mass(xP)[M.free]
    w2, V = eigsh(K, k=NM, M=diags(Md).tocsc(), sigma=0, which="LM")
    o = np.argsort(w2); V = V[:, o]
    V = V / np.sqrt((V * (Md[:, None] * V)).sum(0))       # 質量正規化
    return w2[o], V, Md


def mac(Va, Vb, Md):
    num = (Va.T @ (Md[:, None] * Vb)) ** 2
    return num / np.outer((Va * (Md[:, None] * Va)).sum(0), (Vb * (Md[:, None] * Vb)).sum(0))


def dlam(xP, lam, v):
    full = np.zeros(M.ndof); full[M.free] = v; Ue = full[M.edof]
    dk = M.PEN * xP**(M.PEN - 1) * (np.dot(Ue, M.KE) * Ue).sum(1)
    _, dm = M.mass_interp(xP); dmm = dm / 4 * (Ue**2).sum(1)
    return dk - lam * dmm


def run(tracking, iters=120, move=0.05):
    x = M.VOL * np.ones(M.nelx * M.nely); x[M.passive] = 1.0
    ref = None; hist = []
    for it in range(iters):
        xP = np.asarray(M.H * x[None].T / M.Hs[None].T)[:, 0]
        lam, V, Md = eig(xP)
        if ref is None: ref = V[:, 0].copy()
        mc = mac(ref[:, None], V, Md)[0]
        k = int(np.argmax(mc)) if tracking else 0          # 最大化するモード
        if tracking: ref = V[:, k].copy()
        hist.append((it, np.sqrt(lam), k, mc.copy()))
        dc = dlam(xP, lam[k], V[:, k])
        dc = np.asarray(M.H * (dc[None].T / M.Hs[None].T))[:, 0]; dc = np.maximum(dc, 1e-12 * np.abs(dc).max())
        l1, l2 = 0, 1e12
        while (l2 - l1) / (l1 + l2) > 1e-4:
            lm = 0.5 * (l1 + l2)
            xn = np.clip(x * np.sqrt(dc / lm), np.maximum(0, x - move), np.minimum(1, x + move)); xn[M.passive] = 1.0
            l1, l2 = (lm, l2) if xn.mean() > M.VOL else (l1, lm)
        x = xn
    return hist, xP


hA, xA = run(False); hB, xB = run(True)
w0 = hA[0][1][0]
fig = plt.figure(figsize=(17, 6.4))
gs = fig.add_gridspec(1, 4, width_ratios=[1.5, 1.5, 0.45, 0.45])
for col, (h, ttl) in enumerate(((hA, "A：いつも「今の 1 次」を上げる（追跡なし）\n→ 1 次と 2 次がくっつき、2 次以上が反復ごとにジグザグ（入れ替わりの兆候）"), (hB, "B：最初の 1 次と同じ形のモードを MAC で追跡して上げる\n→ その形は 4 次まで上がるが、別の形が新しい 1 次になり開始時より低い"))):
    ax = fig.add_subplot(gs[col])
    it = np.array([q[0] for q in h]); om = np.array([q[1] for q in h]) / w0
    for m, c in zip(range(NM), (BLUE, RED, GREEN, ORANGE)):
        ax.plot(it, om[:, m], color=c, lw=1.8, label=f"{m+1} 次（小さい順）")
    tgt = np.array([q[1][q[2]] for q in h]) / w0
    ax.plot(it, tgt, "k.", ms=3, label="最大化しているモード")
    sw = [q[0] for i, q in enumerate(h[1:], 1) if (not col and np.max(q[3]) > 0 and int(np.argmax(q[3])) != 0) or (col and q[2] != h[i - 1][2])]
    if col == 0:
        sw = [q[0] for q in h if int(np.argmax(q[3])) != 0]
    ax.set_xlabel("反復回数"); ax.set_ylabel("固有振動数（開始時の 1 次 = 1）"); ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower right")
    ax.set_title(ttl, loc="left", fontsize=11)
for k, (x, ttl) in enumerate(((xA, "A の形"), (xB, "B の形"))):
    ax = fig.add_subplot(gs[2 + k]); ax.imshow(1 - x.reshape((M.nelx, M.nely)).T, cmap="gray", vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(ttl, fontsize=11)
fig.suptitle("1 次固有振動数の最大化では、1 次と 2 次が近づいて入れ替わる（コラム断面、材料 35 %）", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "mode_switching_eigmax.png"))
for name, h in (("A", hA), ("B", hB)):
    om = np.array([q[1] for q in h]) / w0
    print(name, "final w1..4:", np.round(om[-1], 3), "switch iters:", sum(1 for q in h if int(np.argmax(q[3])) != 0) if name == "A" else sum(1 for i in range(1, len(h)) if h[i][2] != h[i-1][2]))
