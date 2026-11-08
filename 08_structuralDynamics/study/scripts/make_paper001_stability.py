"""paper001 の図: Law・Altintas・Phani (2013) の「位置と送り方向で変わるびびり安定限界」を、このプロジェクトの機械モデル（ana002）で再現する。
   - FRF: ana002 の CalculiX 固有モード（30 個、質量正規化）からモード重ね合わせで、テーブル上面の中央と角の XY 2×2 を作る
   - 安定限界: Altintas・Budak のゼロ次近似。送り方向 θ ごとに FRF を回転し、速度によらない絶対安定限界 a_lim,min を求める
       a_lim = -2π Λ_R / (N_t K_t) · (1 + (Λ_I/Λ_R)^2)、Λ = 1/λ0、λ0 = [α0][G(ω_c)] の固有値
   ※ 本来は工具先端とワークの「相対」FRF を使う。ana002 のモデルはテーブル側だけなので、ここでは方法の練習として
     テーブル上面の点の FRF をワーク側の柔らかさとして使う（数値は実機の値ではない）"""
import os, sys
import numpy as np
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
A = os.path.join(os.path.dirname(__file__), "..", "..", "analysis")
sys.path.insert(0, os.path.join(A, "ana002_ccx_fistr_frf"))
from make_inputs import ALPHA, BETA
OUT = os.path.join(os.path.dirname(__file__), "..", "papers", "images", "paper001"); os.makedirs(OUT, exist_ok=True)
MON, CORNER = 20428, 22255


def read_modes(path, nodes, nmax=30):
    fr, phi = [], {n: [] for n in nodes}; disp = False; cnt = 0
    for ln in open(path):
        if ln.startswith("  100CL"):
            freq = float(ln[12:24]); disp = False
        elif ln.startswith(" -4  DISP"):
            disp = True; fr.append(freq); cnt += 1
        elif disp and ln.startswith(" -1") and int(ln[3:13]) in phi:
            phi[int(ln[3:13])].append([float(ln[13:25]), float(ln[25:37]), float(ln[37:49])])
        elif disp and ln.startswith(" -3"):
            disp = False
            if cnt == nmax: break
    return np.array(fr), {n: np.array(v) for n, v in phi.items()}


fr, PH = read_modes(os.path.join(A, "ana002_ccx_fistr_frf", "ccx", "machine.frd"), [MON, CORNER])
wr = 2 * np.pi * fr; zr = ALPHA / (2 * wr) + BETA * wr / 2
f = np.linspace(30, 1200, 6000); w = 2 * np.pi * f


def G_xy(phi):
    """XY の 2×2 FRF [mm/N]（加振点 = 測定点）"""
    p = phi[:, :2]
    den = wr[:, None] ** 2 - w[None, :] ** 2 + 2j * zr[:, None] * wr[:, None] * w[None, :]
    return np.einsum("ri,rj,rw->wij", p, p, 1 / den)


# 切削条件（論文と同じ：AISI4340 の正面フライス、5 枚刃、Kt = 3000 N/mm^2、Kr = 0.24、半径方向かみ込み 80 %、対称）
Nt, Kt, Kr = 5, 3000.0, 0.24
s = np.arcsin(0.8 / 2 * 2 / 2)          # ae/D = 0.8 の対称かみ込み → 中心から ±asin(0.4)
phs, phe = np.pi / 2 - np.arcsin(0.4), np.pi / 2 + np.arcsin(0.4)


def alpha0():
    F = lambda p: np.array([[np.cos(2 * p) - 2 * Kr * p + Kr * np.sin(2 * p), -np.sin(2 * p) - 2 * p + Kr * np.cos(2 * p)],
                            [-np.sin(2 * p) + 2 * p + Kr * np.cos(2 * p), -np.cos(2 * p) - 2 * Kr * p - Kr * np.sin(2 * p)]]) / 2
    return F(phe) - F(phs)


A0 = alpha0()


def alim_min(G, theta):
    """送り方向 θ での速度によらない絶対安定限界 [mm] と、そのときのびびり周波数 [Hz]"""
    c, s_ = np.cos(theta), np.sin(theta); Rm = np.array([[c, s_], [-s_, c]])
    Guv = Rm @ G @ Rm.T
    lam0 = np.linalg.eigvals(A0[None] @ Guv)             # [α0][G] の固有値 λ0
    lam = 1.0 / lam0                                      # det([I] − Λ[α0][G]) = 0 → Λ = 1/λ0
    LR, LI = lam.real, lam.imag
    with np.errstate(divide="ignore", invalid="ignore"):
        a = -2 * np.pi * LR / (Nt * Kt) * (1 + (LI / LR) ** 2)
    a[(LR >= 0) | ~np.isfinite(a)] = np.inf
    k = np.unravel_index(np.argmin(a), a.shape)
    return a[k], f[k[0]]


th = np.radians(np.arange(0, 361, 3))
res = {}
for name, node in (("中央", MON), ("角", CORNER)):
    G = G_xy(PH[node]); out = np.array([alim_min(G, t) for t in th]); res[name] = (G, out)
    i = np.argmin(out[:, 0])
    print(name, f"最小 {out[i,0]:.3f} mm @ {np.degrees(th[i]):.0f}°, びびり {out[i,1]:.0f} Hz; 最大 {out[:,0].max():.3f} mm")

fig = plt.figure(figsize=(16, 6.4))
ax = fig.add_subplot(1, 3, 1)
for (name, (G, out)), c in zip(res.items(), (BLUE, ORANGE)):
    ax.semilogy(f, np.abs(G[:, 0, 0]) * 1e3, color=c, lw=2, label=f"{name}：X-X")
    ax.semilogy(f, np.abs(G[:, 1, 1]) * 1e3, "--", color=c, lw=1.6, label=f"{name}：Y-Y")
    ax.semilogy(f, np.abs(G[:, 0, 1]) * 1e3, ":", color=c, lw=1.6, label=f"{name}：X-Y")
ax.set_xlim(50, 1000); ax.set_ylim(1e-5, 1); ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("コンプライアンス [µm/N]")
ax.set_title("① 位置ごとの FRF（ana002 のモデル）", loc="left", fontsize=12); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8.5, ncol=2)
ax = fig.add_subplot(1, 3, 2, projection="polar")
for (name, (G, out)), c in zip(res.items(), (BLUE, ORANGE)):
    ax.plot(th, out[:, 0], color=c, lw=2.5, label=name)
ax.set_title("② 送り方向ごとの絶対安定限界 [mm]\n（内側が安定。半径 = 削れる深さ）", fontsize=12)
ax.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.12))
ax = fig.add_subplot(1, 3, 3)
for (name, (G, out)), c in zip(res.items(), (BLUE, ORANGE)):
    ax.plot(np.degrees(th), out[:, 1], "o", ms=3.5, color=c, label=name)
for k in range(6):
    ax.axhline(fr[k], color=GRAY, lw=0.6, ls=":"); ax.text(362, fr[k], f"{k+1}次", fontsize=8, color=GRAY, va="center")
ax.set_xlim(0, 360); ax.set_xlabel("送り方向 θ [°]"); ax.set_ylabel("びびり周波数 [Hz]")
ax.set_title("③ どのモードでびびるか", loc="left", fontsize=12); ax.grid(alpha=0.3); ax.legend()
fig.suptitle("位置と送り方向で変わるびびり安定限界（Law ら 2013 の評価方法を ana002 のモデルで再現。練習用で実機の値ではない）", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "position_stability_demo.png")); print("saved")
