"""ref008 の図: 多目的トポロジー最適化（重み付き和）で、静剛性とびびりモードの動的コンプライアンスを同時に下げる
   J(x) = w * Cs(x)/Cs0 + (1-w) * Cd(x)/Cd0
     Cs = f^T u（静的コンプライアンス）, Cd = |f^T u(ω_e)|（加振周波数 ω_e での動的コンプライアンス）
   ω_e は「びびりモード」= 静的最適形状で力によく反応する曲げモードの 0.9 倍。w を 1 → 0 と変えてパレート線を描く。
   題材は ref005 と同じコラム断面（make_ref005_dyn.py の関数を使う）"""
import os
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import make_ref005_dyn as M
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); os.makedirs(OUT, exist_ok=True)


def solve(xP, omega):
    K = M.stiff(xP); Md = M.mass(xP)
    A = (K - diags(omega**2 * Md)).tocsc()[M.free][:, M.free]
    U = np.zeros(M.ndof); U[M.free] = spsolve(A, M.F[M.free]); return U


def sens(xP, U, omega):
    Ue = U[M.edof]; C = M.F @ U
    dk = M.PEN * xP**(M.PEN - 1) * (np.dot(Ue, M.KE) * Ue).sum(1)
    _, dm = M.mass_interp(xP); dmm = dm / 4 * (Ue**2).sum(1)
    return abs(C), -np.sign(C) * (dk - omega**2 * dmm)


def optimize_w(w, we, Cs0, Cd0, x0, iters=60, move=0.08):
    x = x0.copy(); x[M.passive] = 1.0
    for it in range(iters):
        xP = np.asarray(M.H * x[None].T / M.Hs[None].T)[:, 0]
        cs, dcs = sens(xP, solve(xP, 0.0), 0.0)
        cd, dcd = sens(xP, solve(xP, we), we)
        dc = w * dcs / Cs0 + (1 - w) * dcd / Cd0
        dc = np.asarray(M.H * (dc[None].T / M.Hs[None].T))[:, 0]
        dc = np.minimum(dc, -1e-12 * np.abs(dc).max())
        l1, l2 = 0, 1e12
        while (l2 - l1) / (l1 + l2) > 1e-4:
            lm = 0.5 * (l1 + l2)
            xn = np.clip(x * np.sqrt(-dc / lm), np.maximum(0, x - move), np.minimum(1, x + move)); xn[M.passive] = 1.0
            l1, l2 = (lm, l2) if xn.mean() > M.VOL else (l1, lm)
        x = xn
    xP = np.asarray(M.H * x[None].T / M.Hs[None].T)[:, 0]
    return xP, x


if __name__ == "__main__":
    xs, _, xs_raw = M.optimize(0.0)                       # 静的最適
    ws = M.eig1(xs); we = 0.9 * ws[1]
    Cs0 = M.F @ solve(xs, 0.0); Cd0 = abs(M.F @ solve(xs, we))
    res = []; xr = xs_raw
    for w in (1.0, 0.75, 0.5, 0.25, 0.0):                 # 静的寄り → 動的寄りへ少しずつ（継続法）
        xP, xr = optimize_w(w, we, Cs0, Cd0, xr)
        cs = M.F @ solve(xP, 0.0); cd = abs(M.F @ solve(xP, we)); om = M.eig1(xP)
        res.append(dict(w=w, x=xP, cs=cs / Cs0, cd=cd / Cd0, w1=om[0] / ws[0], w2=om[1] / ws[0]))
        print(f"w={w}: Cs={cs/Cs0:.3f} Cd={cd/Cd0:.3f} w1={om[0]/ws[0]:.3f} w2={om[1]/ws[0]:.3f}", flush=True)
    wn = np.linspace(0.05, 2.5, 300) * ws[0]

    fig = plt.figure(figsize=(16, 9.5))
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 0.95])
    for k, r in enumerate(res):
        ax = fig.add_subplot(gs[0, k]); ax.imshow(1 - r["x"].reshape((M.nelx, M.nely)).T, cmap="gray", vmin=0, vmax=1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$w$ = {r['w']:.2f}\n静的 {r['cs']:.2f}・動的 {r['cd']:.2f}", fontsize=11,
                     color=BLUE if r["w"] == 1 else (RED if r["w"] == 0 else "k"))
    fig.axes[0].text(0, M.nely + 6, "← 静剛性を重視", color=BLUE, fontsize=11); fig.axes[4].text(0, M.nely + 6, "びびりモード重視 →", color=RED, fontsize=11)
    ax = fig.add_subplot(gs[1, :2])
    ax.plot([r["cs"] for r in res], [r["cd"] for r in res], "o-", color=GRAY, lw=1.5)
    for r in res:
        ax.plot(r["cs"], r["cd"], "o", ms=12, color=plt.cm.coolwarm(1 - r["w"]))
        ax.annotate(f"$w$={r['w']:.2f}", (r["cs"], r["cd"]), xytext=(8, 4), textcoords="offset points", fontsize=10)
    ax.set_xlabel("静的コンプライアンス（静的最適 = 1、小さいほど硬い）"); ax.set_ylabel("びびり周波数 $\\omega_e$ での動的コンプライアンス（同上）")
    ax.set_title("パレート線：片方を良くすると、もう片方が少し悪くなる", loc="left", fontsize=12); ax.grid(alpha=0.3)
    ax = fig.add_subplot(gs[1, 2:])
    for r in (res[0], res[2], res[4]):
        H = M.frf(r["x"], wn)
        ax.semilogy(wn / ws[0], H, lw=2, color=plt.cm.coolwarm(1 - r["w"]), label=f"$w$ = {r['w']:.2f}")
    ax.axvline(we / ws[0], color=GRAY, ls="--"); ax.text(we / ws[0] + 0.03, 0.92, "びびり周波数 $\\omega_e$", transform=ax.get_xaxis_transform(), color=GRAY)
    ax.set_xlabel("周波数（静的最適形状の 1次を 1）"); ax.set_ylabel("力の位置の $|u/f|$"); ax.grid(alpha=0.3, which="both"); ax.legend()
    ax.set_title("周波数応答（同じ材料量 35 %、減衰 2 %）", loc="left", fontsize=12)
    fig.suptitle("多目的トポロジー最適化の例：$J = w\\,C_s/C_{s0} + (1-w)\\,C_d(\\omega_e)/C_{d0}$（コラム断面、材料 35 %）", fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "multiobj_pareto.png")); plt.close(fig)
    open(os.path.join(OUT, "multiobj_result.md"), "w").write("| w | 静的コンプライアンス | 動的コンプライアンス(ω_e) | 1次固有振動数 | 2次固有振動数 |\n|---|---|---|---|---|\n" +
        "".join(f"| {r['w']:.2f} | {r['cs']:.3f} | {r['cd']:.3f} | {r['w1']:.3f} | {r['w2']:.3f} |\n" for r in res))
    print("saved")
