"""5.3 補足：ω_c を 5〜3,000 Hz で動かし、k ごとのローブを重ね、最小をつないで線図ができる様子"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
TEETH, KT, KR = 4, 1000.0, 0.3; B0 = TEETH / 4 * np.array([[KR, 1.0], [-1.0, KR]])
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2; d = phi[iw, :2] - phi[it, :2]
f = np.linspace(5, 3000, 60000); w = 2 * np.pi * f
G = np.einsum("ir,jr,wr->wij", d, d, 1 / (wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]))
lam = np.linalg.eigvals(G @ B0); mu = -1 / lam                       # (nw, 2)
A = np.where(mu.real > 0, mu.real * (1 + (mu.imag / mu.real) ** 2) / (2 * KT), np.nan)
eps = np.pi - 2 * np.arctan(mu.imag / mu.real)
fig, axs = plt.subplots(1, 2, figsize=(17, 6.6), gridspec_kw=dict(width_ratios=[0.8, 1.2]))
ax = axs[0]
for j in (0, 1): ax.semilogx(f, A[:, j], color="#555", lw=1.4)
ax.set_xlim(100, 3000); ax.set_ylim(0, 60); ax.grid(alpha=0.3, which="both")
ax.set_xlabel("びびり周波数 $\\omega_c/2\\pi$ [Hz]", fontsize=13); ax.set_ylabel("$a_{\\lim}$ [mm]", fontsize=13)
ax.set_title("① まず周波数ごとに $a_{\\lim}$ を計算する\n（$\\mu_R>0$ の所だけ。山の少し右で小さくなる）", fontsize=13.5, loc="left")
pick = [(185, "#c0392b"), (350, "#2980b9"), (920, "#27ae60")]
for fp, c in pick:
    i = np.argmin(abs(f - fp)); j = np.nanargmin(A[i]); ax.plot(f[i], A[i, j], "o", color=c, ms=10, zorder=5)
    ax.annotate(f"{fp} Hz：{A[i, j]:.1f} mm", (f[i], A[i, j]), xytext=(8, 10), textcoords="offset points", color=c, fontsize=12)
ax = axs[1]
rpm = np.linspace(500, 16000, 1500); env = np.full_like(rpm, np.inf)
cols = plt.cm.viridis(np.linspace(0, 0.9, 8))
for k in range(8):
    for j in (0, 1):
        n = 60 * w / (TEETH * (2 * np.pi * k + eps[:, j])); ok = np.isfinite(A[:, j]) & (n > 400) & (n < 17000)
        ax.plot(n[ok], A[ok, j], ".", ms=1.5, color=cols[k], alpha=0.7)
    ax.plot([], [], "o", color=cols[k], label=f"$k={k}$")
sc = np.load(R / "stability_compare.npz")
ax.plot(sc["rpm"], sc["lim_C"], color="k", lw=2.6, label="いちばん小さい所をつないだ線\n＝安定限界線図")
for fp, c in pick:
    i = np.argmin(abs(f - fp)); j = np.nanargmin(A[i])
    for k in range(8):
        n = 60 * w[i] / (TEETH * (2 * np.pi * k + eps[i, j]))
        if 500 < n < 16000: ax.plot(n, A[i, j], "o", color=c, ms=8, mec="k", zorder=6)
ax.set_xlim(500, 16000); ax.set_ylim(0, 60); ax.grid(alpha=0.3)
ax.set_xlabel("主軸回転数 [min$^{-1}$]", fontsize=13); ax.set_ylabel("軸方向切込み限界 [mm]", fontsize=13)
ax.set_title("② 同じ点を $k = 0, 1, 2, \\dots$ ごとに回転数 $n = 60\\,\\omega_c/(Z(2\\pi k+\\varepsilon))$ へ移す\n（色つきの丸＝①の 1 点が、$k$ ごとに別の回転数に現れる。黒線＝最小をつないだもの）", fontsize=13.5, loc="left")
ax.legend(fontsize=10.5, loc="upper left", ncol=3, markerscale=1.5)
fig.tight_layout(); fig.savefig(R / "lobes_explained.png", dpi=95); print("ok")
