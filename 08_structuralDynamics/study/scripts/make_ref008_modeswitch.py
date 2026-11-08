"""ref008 の図: トポロジー最適化の途中でモードの順番が入れ替わるか（MAC で追跡）
   ref008 の多目的最適化（静剛性 → びびりモード重視へ重みを移す継続法）の途中で、低次 4 モードの固有振動数と形を記録する。
   MAC(φ_a, φ_b) = (φ_a^T M φ_b)^2 / ((φ_a^T M φ_a)(φ_b^T M φ_b)) で前の段階のモードと対応付ける"""
import os
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import make_ref005_dyn as M
import make_ref008_multiobj as MO
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008")
NM = 4


def modes(xP):
    K = M.stiff(xP)[M.free][:, M.free]; Md = M.mass(xP)[M.free]
    w2, V = eigsh(K, k=NM, M=diags(Md).tocsc(), sigma=0, which="LM")
    o = np.argsort(w2); return np.sqrt(w2[o]), V[:, o], Md


xs, _, xs_raw = M.optimize(0.0)
ws = M.eig1(xs); we = 0.9 * ws[1]
Cs0 = M.F @ MO.solve(xs, 0.0); Cd0 = abs(M.F @ MO.solve(xs, we))
hist = []; xr = xs_raw.copy(); step = 0
for w in (1.0, 0.75, 0.5, 0.25, 0.0):
    for _ in range(6):                                   # 10 反復ずつ記録
        xP, xr = MO.optimize_w(w, we, Cs0, Cd0, xr, iters=10)
        om, V, Md = modes(xP); hist.append((step, w, om / ws[0], V, Md, xP)); step += 10
# MAC で追跡（前の段階の各モードに、いちばん形が似ているものを対応させる）
track = [list(range(NM))]; snaps = [0]
for k in range(1, len(hist)):
    Vp, Mp = hist[k - 1][3], hist[k - 1][4]; Vc, Mc = hist[k][3], hist[k][4]
    mac = np.array([[(Vp[:, i] @ (Mc * Vc[:, j])) ** 2 / ((Vp[:, i] @ (Mc * Vp[:, i])) * (Vc[:, j] @ (Mc * Vc[:, j]))) for j in range(NM)] for i in range(NM)])
    perm = [int(np.argmax(mac[i])) for i in range(NM)]
    track.append([perm[i] for i in range(NM)])
# 各「同じ形のモード」が、何番目の固有振動数にいるか
ident = np.zeros((len(hist), NM), int); ident[0] = np.arange(NM)
for k in range(1, len(hist)):
    for m in range(NM):
        ident[k, m] = track[k][ident[k - 1, m]]
freq_by_shape = np.array([[hist[k][2][ident[k, m]] for m in range(NM)] for k in range(len(hist))])
order = ident
it = np.array([h[0] for h in hist]); wl = np.array([h[1] for h in hist])
swaps = [(it[k], m) for k in range(1, len(hist)) for m in range(NM) if order[k, m] != order[k - 1, m]]
print("swaps:", swaps[:10])

fig, axs = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw=dict(width_ratios=[1.4, 1]))
ax = axs[0]
cols = [BLUE, RED, GREEN, ORANGE]
for m in range(NM):
    ax.plot(it, freq_by_shape[:, m], "o-", color=cols[m], ms=4, lw=2, label=f"最初の {m+1} 次モードと同じ形のモード")
for k in range(len(hist)):
    for j in range(NM):
        ax.plot(it[k], hist[k][2][j], ".", color="k", ms=3)
for (s, m) in swaps:
    ax.axvline(s, color=GRAY, ls=":", lw=1)
for k in range(1, len(wl)):
    if wl[k] != wl[k - 1]:
        ax.text(it[k], ax.get_ylim()[1] if False else freq_by_shape.max() * 1.02, f"w={wl[k]:.2f}", fontsize=9, color=GRAY, ha="center")
ax.set_xlabel("最適化の反復回数（重み w を 1 → 0 へ移していく）"); ax.set_ylabel("固有振動数（静的最適形状の 1 次 = 1）")
ax.set_title("同じ形のモードを MAC で追いかけると、線が交差する（点線 = 順番が入れ替わった所）", loc="left", fontsize=12)
ax.grid(alpha=0.3); ax.legend(fontsize=9.5)
ax = axs[1]
ax.imshow(order.T + 1, aspect="auto", cmap="Set2", extent=[it[0] - 5, it[-1] + 5, NM + 0.5, 0.5])
for k in range(0, len(hist), 3):
    for m in range(NM):
        ax.text(it[k], m + 1, str(order[k, m] + 1), ha="center", va="center", fontsize=9)
ax.set_yticks(range(1, NM + 1)); ax.set_yticklabels([f"最初の {m+1} 次の形" for m in range(NM)])
ax.set_xlabel("反復回数"); ax.set_title("その形が今「何次」か（数字）", loc="left", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "mode_switching.png")); print("saved", len(swaps))
