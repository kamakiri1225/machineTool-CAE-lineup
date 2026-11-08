"""ana005 4.2：ジャークとは何か、「ジャーク立上げ 40 ms」はどこか"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
import fn_effect as FE

fig, axs = plt.subplots(4, 1, figsize=(13, 12), sharex=True)
for tj, col, lab in ((0.0, RED, "ジャーク立上げ 0 ms（加速度段差）"), (0.04, GREEN, "ジャーク立上げ 40 ms")):
    t, a, tm = FE.profile(tj)
    v = np.cumsum(a) * FE.dt; x = np.cumsum(v) * FE.dt; j = np.gradient(a, FE.dt)
    axs[0].plot(t, x * 1000, color=col, lw=2, label=lab)
    axs[1].plot(t, v * 1000, color=col, lw=2, label=lab)
    axs[2].plot(t, a, color=col, lw=2.4, label=lab)
    if tj > 0: axs[3].plot(t, j, color=col, lw=2, label=lab)
# ジャーク 0 は段差の瞬間に無限大
for tt, sgn in ((0, 1), (0.1, -1), (0.25, -1), (0.35, 1)):
    axs[3].annotate("", xy=(tt, 95 * sgn), xytext=(tt, 0), arrowprops=dict(arrowstyle="-|>", color=RED, lw=2))
axs[3].text(0.0, 102, "段差の瞬間は ∞", color=RED, fontsize=10)
# 40 ms の場所に印
for t0, t1 in ((0, 0.04), (0.10, 0.14), (0.25, 0.29), (0.35, 0.39)):
    for ax in (axs[2], axs[3]):
        ax.axvspan(t0, t1, color=GREEN, alpha=0.12)
    axs[2].annotate("", xy=(t1, 2.35 if t0 < 0.2 else -2.35), xytext=(t0, 2.35 if t0 < 0.2 else -2.35), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.5))
    axs[2].text((t0 + t1) / 2, 2.5 if t0 < 0.2 else -2.85, "40 ms", ha="center", color=GREEN, fontsize=10.5)
axs[0].set_ylabel("位置 [mm]"); axs[1].set_ylabel("速度 [mm/s]"); axs[2].set_ylabel("加速度 [m/s²]"); axs[3].set_ylabel("ジャーク [m/s³]")
axs[0].set_title("① 位置：どちらも 50 mm 動く（40 ms のほうが 0.04 s 遅く着く）", loc="left", fontsize=12)
axs[1].set_title("② 速度：最高 0.2 m/s。40 ms のほうは角が丸い", loc="left", fontsize=12)
axs[2].set_title("③ 加速度：40 ms のほうは 0 → 2 m/s² まで 40 ms かけてかけて上げ下げする（緑の帯が 40 ms、4 か所）", loc="left", fontsize=12)
axs[3].set_title("④ ジャーク = 加速度の変化の速さ = 上昇の傾き。40 ms なら 2 ÷ 0.04 = 50 m/s³。段差は ∞", loc="left", fontsize=12)
axs[2].set_ylim(-3.1, 3.1); axs[3].set_ylim(-120, 120)
for ax in axs: ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=9.5)
axs[3].set_xlabel("時刻 [s]"); axs[3].set_xlim(0, 0.45)
fig.suptitle("ジャークとは：位置 →（微分）速度 →（微分）加速度 →（微分）ジャーク", fontsize=14)
fig.tight_layout(); fig.savefig(D / "results" / "jerk_explained.png"); print("ok")
