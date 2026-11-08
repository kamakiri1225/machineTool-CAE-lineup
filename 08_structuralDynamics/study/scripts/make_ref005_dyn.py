"""ref005 の図: 動的コンプライアンス最小化のトポロジー最適化（2D、SIMP 法）
   静的コンプライアンス最小化（make_ref005_figs.py と同じ題材）と比べる。
   題材: 底面固定のコラム。上端に主軸頭の質量（集中質量）があり、前面上端に周期的な切削力 f e^{iωt} がかかる。
   目的: |C(ω)| = |f^T u(ω)|,  (K - ω^2 M) u = f   を最小化（材料 35 %）
   参考: Olhoff & Du, "Topological design of continuum structures subjected to forced vibration" (2005) など"""
import os
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve, eigsh
from figstyle import plt, RED, BLUE, GRAY, ORANGE, GREEN
from make_ref005_figs import lk

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref005")
nelx, nely = 40, 80
VOL, PEN, RMIN = 0.35, 3.0, 2.0
ndof = 2 * (nelx + 1) * (nely + 1)
nid = lambda ix, iy: (nely + 1) * ix + iy             # iy=0 が上端
KE = lk()
edof = np.zeros((nelx * nely, 8), dtype=int)
for elx in range(nelx):
    for ely in range(nely):
        n1 = (nely + 1) * elx + ely; n2 = (nely + 1) * (elx + 1) + ely
        edof[ely + elx * nely] = [2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3, 2*n2, 2*n2+1, 2*n1, 2*n1+1]
iK = np.kron(edof, np.ones((8, 1))).flatten(); jK = np.kron(edof, np.ones((1, 8))).flatten()
fixed = np.array([d for ix in range(nelx + 1) for d in (2 * nid(ix, nely), 2 * nid(ix, nely) + 1)])
free = np.setdiff1d(np.arange(ndof), fixed)
F = np.zeros(ndof)
for iy in range(0, 4):
    F[2 * nid(0, iy)] = 1.0
# 主軸頭の質量（上端の節点に配る）: 構造（全部材料）の 25 % 相当
M_HEAD = 0.25 * nelx * nely
m_head = np.zeros(ndof)
for ix in range(nelx + 1):
    for d in (0, 1):
        m_head[2 * nid(ix, 0) + d] = M_HEAD / (nelx + 1)

# 上端 3 行は主軸頭の取付面（設計しない・常に材料）
passive = np.array([(e % nely) < 3 for e in range(nelx * nely)])

# フィルタ
iH, jH, sH = [], [], []
for i in range(nelx):
    for j in range(nely):
        for k in range(max(i - int(RMIN) + 1, 0), min(i + int(RMIN), nelx)):
            for l in range(max(j - int(RMIN) + 1, 0), min(j + int(RMIN), nely)):
                iH.append(i * nely + j); jH.append(k * nely + l); sH.append(max(0, RMIN - np.hypot(i - k, j - l)))
H = coo_matrix((sH, (iH, jH)), shape=(nelx*nely, nelx*nely)).tocsc(); Hs = np.array(H.sum(1)).flatten()


def stiff(x):
    e = 1e-9 + x**PEN * (1 - 1e-9)
    return coo_matrix(((KE.flatten()[None]).T * e).flatten(order="F"), shape=(ndof, ndof)).tocsc() if False else \
        coo_matrix((np.outer(e, KE.flatten()).flatten(), (np.repeat(edof, 8, axis=1).flatten(), np.tile(edof, (1, 8)).flatten())),
                   shape=(ndof, ndof)).tocsc()


def mass_interp(x):
    """低密度で局所モードが出ないよう、x<0.1 では質量を急に小さくする（Pedersen 2000）"""
    return np.where(x > 0.1, x, 1e5 * x**6), np.where(x > 0.1, 1.0, 6e5 * x**5)


def mass(x):
    m, _ = mass_interp(x)
    md = np.zeros(ndof)
    np.add.at(md, edof, (m / 4)[:, None] * np.ones((1, 8)))     # 集中質量（節点に 1/4 ずつ）
    return md + m_head


def eig1(x, n=3):
    K = stiff(x)[free][:, free]; Md = mass(x)[free]
    w2, _ = eigsh(K, k=n, M=diags(Md).tocsc(), sigma=0, which="LM")
    return np.sqrt(np.sort(w2))


def optimize(omega, iters=120, move=0.1, x0=None):
    x = VOL * np.ones(nelx * nely) if x0 is None else x0.copy(); x[passive] = 1.0; hist = []
    for it in range(iters):
        xP = np.asarray(H * x[None].T / Hs[None].T)[:, 0]
        K = stiff(xP); Md = mass(xP)
        A = (K - diags(omega**2 * Md)).tocsc()[free][:, free]
        U = np.zeros(ndof); U[free] = spsolve(A, F[free])
        C = F @ U; hist.append(abs(C))
        Ue = U[edof]
        dk = PEN * xP**(PEN - 1) * (np.dot(Ue, KE) * Ue).sum(1)
        _, dm = mass_interp(xP); dmm = dm / 4 * (Ue**2).sum(1)
        dc = -np.sign(C) * (dk - omega**2 * dmm)                   # d|C|/dx
        dc = np.asarray(H * (dc[None].T / Hs[None].T))[:, 0]        # 密度フィルタの連鎖律
        dc = np.minimum(dc, -1e-12 * np.abs(dc).max())              # OC 用に負にそろえる（簡易）
        l1, l2 = 0, 1e12
        while (l2 - l1) / (l1 + l2) > 1e-4:
            lm = 0.5 * (l1 + l2)
            xn = np.clip(x * np.sqrt(-dc / lm), np.maximum(0, x - move), np.minimum(1, x + move))
            xn[passive] = 1.0
            l1, l2 = (lm, l2) if xn.mean() > VOL else (l1, lm)
        x = xn; x[passive] = 1.0
    xP = np.asarray(H * x[None].T / Hs[None].T)[:, 0]
    return xP, hist, x


def frf(x, ws, eta=0.02):
    K = stiff(x); Md = mass(x); out = []
    for w in ws:
        A = ((1 + 1j * eta) * K - diags(w**2 * Md)).tocsc()[free][:, free]
        u = np.zeros(ndof, complex); u[free] = spsolve(A, F[free]); out.append(abs(F @ u))
    return np.array(out)


if __name__ == "__main__":
    xs, _, xs_raw = optimize(0.0)                                  # 静的コンプライアンス最小（ω = 0）
    ws = eig1(xs)
    w_exc = 0.9 * ws[1]          # 力でよく揺れるモード（2次・コラムの曲げ）の 0.9 倍で加振（共振の手前）
    xbad, _, _ = optimize(w_exc)                          # いきなり ω_e で解くと「支えを切る」解に落ちる
    print("bad w:", eig1(xbad))
    xr = xs_raw; hist = []
    for fr in (0.3, 0.5, 0.7, 0.85, 1.0):                       # 周波数を少しずつ上げる（静的最適から出発）
        xd, h, xr = optimize(fr * w_exc, iters=60, x0=xr); hist += h
    wd = eig1(xd)
    np.save(os.path.join(os.path.dirname(__file__), "ref005_dyn_bad.npy"), xbad)
    print("hist", np.round(hist[::10], 1)); print("static-opt w:", ws, "dyn-opt w:", wd, "w_exc", w_exc)
    wn = np.linspace(0.05, 3.0, 400) * ws[0]
    Hs_, Hd_ = frf(xs, wn), frf(xd, wn)
    np.savez(os.path.join(os.path.dirname(__file__), "ref005_dyn_result.npz"), xs=xs, xd=xd, ws=ws, wd=wd)

    fig = plt.figure(figsize=(15, 6.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.55, 0.55, 1.5])
    for k, (x, ttl, c) in enumerate([(xs, "静的コンプライアンス最小\n（$\\omega = 0$）", BLUE),
                                     (xd, "動的コンプライアンス最小\n（$\\omega_e = 0.9\\,\\omega_{2,\\mathrm{static}}$）", RED)]):
        ax = fig.add_subplot(gs[k]); ax.imshow(1 - x.reshape((nelx, nely)).T, cmap="gray", vmin=0, vmax=1)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(ttl, fontsize=12, color=c)
        ax.add_patch(plt.Rectangle((-0.5, -2.5), nelx, 2, color=ORANGE, clip_on=False))
        ax.annotate("", xy=(-0.5, 1.5), xytext=(-9, 1.5), arrowprops=dict(arrowstyle="->", color=RED, lw=2.5), annotation_clip=False)
        ax.plot([-0.5, nelx - 0.5], [nely - 0.5, nely - 0.5], color=BLUE, lw=5)
    fig.axes[0].text(0, -4, "主軸頭の質量", color=ORANGE, fontsize=10)
    ax = fig.add_subplot(gs[2])
    ax.semilogy(wn / ws[0], Hs_, color=BLUE, lw=2.3, label="静的コンプライアンス最小の形")
    ax.semilogy(wn / ws[0], Hd_, color=RED, lw=2.3, label="動的コンプライアンス最小の形")
    ax.axvline(w_exc / ws[0], color=GRAY, ls="--")
    ax.text(w_exc / ws[0] + 0.03, 0.93, "加振周波数 $\\omega_e$", color=GRAY, fontsize=11, transform=ax.get_xaxis_transform())
    ax.set_xlabel("周波数 $\\omega / \\omega_{1,\\mathrm{static}}$（静的最適形状の 1次を 1 とする）")
    ax.set_ylabel("力の位置の動コンプライアンス $|u/f|$")
    ax.set_title("周波数応答（同じ材料量 35 %、減衰 $\\eta$ = 2 %）", loc="left", fontsize=12)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=11)
    ie = np.argmin(abs(wn - w_exc)); cs, cd = Hs_[ie], Hd_[ie]
    ax.annotate(f"加振周波数での揺れ\n{cs:.0f} → {cd:.0f}（{(cd/cs-1)*100:+.0f} %）", xy=(w_exc / ws[0], cd), xytext=(1.75, 3e2),
                color=RED, fontsize=11, arrowprops=dict(arrowstyle="->", color=RED))
    ip = np.argmax(Hd_[wn < 0.9 * w_exc]); 
    ax.annotate("狙っていない周波数に\n新しいピークができる", xy=(wn[ip] / ws[0], Hd_[ip]), xytext=(0.15, 4e4),
                color=RED, fontsize=11, arrowprops=dict(arrowstyle="->", color=RED))
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "topopt_dynamic_compliance.png")); plt.close(fig)
    print("saved")

    # 失敗例：いきなり ω_e で解いたとき
    wb = frf(xbad, wn)
    fig, axs = plt.subplots(1, 2, figsize=(11, 6), gridspec_kw=dict(width_ratios=[0.5, 1]))
    axs[0].imshow(1 - xbad.reshape((nelx, nely)).T, cmap="gray", vmin=0, vmax=1); axs[0].set_xticks([]); axs[0].set_yticks([])
    axs[0].plot([-0.5, nelx - 0.5], [nely - 0.5, nely - 0.5], color=BLUE, lw=5)
    axs[0].set_title("いきなり $\\omega_e$ で解いた結果", fontsize=12, color=RED)
    axs[0].text(nelx / 2, nely * 0.8, "ベッドとの\nつながりが無い", ha="center", color=RED, fontsize=12)
    ax = axs[1]
    ax.loglog(wn / ws[0], wb, color=RED, lw=2.3, label="失敗例（支えが切れた形）")
    ax.loglog(wn / ws[0], Hs_, color=BLUE, lw=2.0, label="静的最適の形")
    ax.axvline(w_exc / ws[0], color=GRAY, ls="--")
    ax.set_xlabel("$\\omega / \\omega_{1,\\mathrm{static}}$"); ax.set_ylabel("$|u/f|$"); ax.grid(alpha=0.3, which="both"); ax.legend()
    ax.set_title("共振より上では「浮いた重り」が揺れにくく見える", loc="left", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "topopt_dynamic_pitfall.png")); plt.close(fig)
    print("saved2")
