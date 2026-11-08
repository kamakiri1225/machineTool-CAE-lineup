"""5.2 補足：なぜ FRF の実部が負の所（固有振動数の少し上）でびびるのか、ε と k は何か"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt

z = 0.05; r = np.linspace(0.5, 1.6, 2000); G = 1 / (1 - r**2 + 2j * z * r)        # 1 自由度、k=1 で無次元化
fig, axs = plt.subplots(1, 3, figsize=(17, 5.4))
ax = axs[0]
ax.plot(r, G.real, color="#1f618d", lw=2.2, label="実部 Re $G$（力と同じ向きの揺れ）")
ax.plot(r, G.imag, color="#7f8c8d", lw=1.6, ls="--", label="虚部 Im $G$（力より 90° 遅れた揺れ）")
ax.axhline(0, color="k", lw=0.8); ax.axvspan(1, 1.6, color="#f5b7b1", alpha=0.35)
ax.text(1.3, 6, "Re $G<0$\nびびりうる側", color="#c0392b", fontsize=14, ha="center")
ax.text(0.72, 6, "Re $G>0$\nびびらない側", color="#1f618d", fontsize=14, ha="center")
ax.set_xlabel("周波数 / 固有振動数", fontsize=13); ax.set_ylabel("コンプライアンス $G$（$k=1$）", fontsize=13)
ax.set_title("(a) FRF の実部と虚部（1 自由度、$\\zeta=0.05$）", fontsize=14, loc="left"); ax.legend(fontsize=10.5, loc="lower left"); ax.set_ylim(-11, 11)

ax = axs[1]
a = -1 / (2 * G.real); a[G.real >= 0] = np.nan
ax.plot(r, a, color="#c0392b", lw=2.4); i = np.nanargmin(a)
ax.plot(r[i], a[i], "o", color="k", ms=8); ax.annotate(f"最小：周波数 ≈ 固有振動数 ×(1+ζ)\n= {r[i]:.2f} 倍", xy=(r[i], a[i]), xytext=(1.12, 0.75), fontsize=13,
                                                     arrowprops=dict(arrowstyle="->"))
ax.axvline(1, color="k", lw=0.8, ls=":"); ax.text(0.98, 0.9, "固有振動数", ha="right", fontsize=12)
ax.text(0.75, 0.35, "Re $G>0$ では\n$a_{\\lim}<0$ になり\n解がない（びびらない）", fontsize=12.5, color="#1f618d", ha="center")
ax.set_ylim(0, 1.0); ax.set_xlim(0.5, 1.6)
ax.set_xlabel("周波数 / 固有振動数", fontsize=13); ax.set_ylabel("$a_{\\lim} = -1/(2K_t\\,\\mathrm{Re}\\,G)$（$K_t=1$）", fontsize=13)
ax.set_title("(b) 周波数ごとの限界切込み", fontsize=14, loc="left")

ax = axs[2]
k, eps = 3, np.pi / 2; t = np.linspace(0, 1, 800); w = 2 * np.pi * (k + eps / (2 * np.pi))   # τ=1 のあいだに k+ε/2π 波
prev = 0.35 * np.sin(w * t); cur = 0.35 * np.sin(w * t + eps) - 1.0 + 1.0
ax.fill_between(t, prev, cur, color="#f5b041", alpha=0.5, label="削る厚さの変動（今 − 前）")
ax.plot(t, prev, color="#7f8c8d", lw=2, ls="--", label="前の刃が残した波")
ax.plot(t, cur, color="#c0392b", lw=2.2, label="今の刃の揺れ")
ax.set_xlabel("時間（刃 1 枚ぶんの時間 $\\tau$ を 1 とする）", fontsize=12.5); ax.set_yticks([])
ax.set_title("(c) $k$ と $\\varepsilon$：前の波と今の波のずれ", fontsize=14, loc="left")
ax.text(0.02, 0.75, f"$\\tau$ のあいだに波が $k + \\varepsilon/2\\pi$ 個入る\n（この図は $k={k}$、$\\varepsilon=\\pi/2$：3 と 1/4 波）\n$\\varepsilon$ だけ位相がずれるので厚さが変動する", fontsize=12, va="top")
ax.set_ylim(-0.8, 0.95); ax.legend(fontsize=10.5, loc="lower right", ncol=1)
fig.tight_layout(); fig.savefig(R / "phase_explained.png", dpi=100, bbox_inches="tight"); print("ok", r[i])
