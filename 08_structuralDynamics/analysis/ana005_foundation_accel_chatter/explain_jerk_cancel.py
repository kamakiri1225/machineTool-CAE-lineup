"""ana005 4.2b：なぜジャーク立上げを周期と同じ長さにすると揺れが消えるのか
   対象：加速の「始まり」1 回だけ（加速度 0 → 2 m/s²）に対する、機械本体の床に対する X 変位 x_b（4.1 の 1 自由度モデル、説明のため減衰なし）
   縦軸は「加速度 2 m/s² が続いたときに落ち着く変位（静的変位 m_t A / k）」を 1 とした比
   加速度の上昇を 8 個の小さな段差に分け、段差 1 つずつの揺れ（細線）と、その合計（太線）を描く"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
fn = 22.9; P = 1 / fn; n = 8
fig = plt.figure(figsize=(17, 14))
gs = fig.add_gridspec(4, 3, height_ratios=[0.55, 0.75, 0.95, 1.05])
t = np.linspace(-0.01, 0.16, 3000)
cases = ((0.0, "段差（0 ms）", RED), (0.02, "上昇時間 20 ms（周期の半分）", ORANGE), (P, f"上昇時間 {P*1000:.0f} ms（周期と同じ）", GREEN))
for col, (T, title, cc) in enumerate(cases):
    steps = np.array([0.0]) if T == 0 else (np.arange(n) + 0.5) / n * T
    vc = [plt.cm.viridis(k / max(1, len(steps) - 1)) for k in range(len(steps))]
    # ① 入力：加速度指令（加速度の上昇を小さな段差に分けたもの）
    ax = fig.add_subplot(gs[0, col])
    a_true = np.clip(t / T, 0, 1) * 2 if T > 0 else np.where(t >= 0, 2.0, 0)
    ax.plot(t * 1000, a_true, color=cc, lw=2.5, label="加速度指令")
    a_st = sum(np.where(t >= s, 2.0 / len(steps), 0) for s in steps)
    ax.step(t * 1000, a_st, color=GRAY, lw=1, where="post", label=f"{len(steps)} 個の小さな段差で近似")
    for k, s in enumerate(steps): ax.plot(s * 1000, 2.0 * (k + 1) / len(steps), "o", color=vc[k], ms=5)
    ax.set_ylim(-0.3, 2.6); ax.set_xlim(-10, 160); ax.grid(alpha=0.3); ax.axvspan(0, T * 1000, color=GREEN, alpha=0.12)
    ax.set_title(f"{title}\n① 入力：加速度指令 [m/s²]（加速の始まり 1 回分）", loc="left", fontsize=11.5)
    if col == 0: ax.legend(fontsize=9, loc="lower right")
    # ② 段差 1 つずつの揺れ
    ax = fig.add_subplot(gs[1, col]); total = np.zeros_like(t)
    for k, s in enumerate(steps):
        y = np.where(t >= s, 1 - np.cos(2 * np.pi * fn * (t - s)), 0) / len(steps)
        total += y; ax.plot(t * 1000, y, color=vc[k], lw=1.3)
    ax.set_xlim(-10, 160); ax.set_ylim(-0.1, 2.1 if T == 0 else 0.3); ax.grid(alpha=0.3); ax.axvspan(0, T * 1000, color=GREEN, alpha=0.12)
    ax.set_title("② 小さな段差 1 つずつが起こす揺れ（色 = ①の丸と対応）\n" + ("段差 1 つ = 全部" if T == 0 else "始まる時刻が少しずつずれる"), loc="left", fontsize=11)
    # ③ 合計
    ax = fig.add_subplot(gs[2, col])
    ax.plot(t * 1000, total, color=cc, lw=2.8); ax.axhline(1, color=GRAY, ls=":")
    ax.set_xlim(-10, 160); ax.set_ylim(-0.2, 2.2); ax.grid(alpha=0.3); ax.axvspan(0, T * 1000, color=GREEN, alpha=0.12)
    res = 1 if T == 0 else abs(np.sin(np.pi * fn * T) / (np.pi * fn * T))
    ax.text(155, 0.1, f"残る揺れ ±{res:.2f}", ha="right", fontsize=13, color=cc, weight="bold")
    ax.set_title("③ 合計 = 本体の変位 $x_b$（②を全部足したもの）", loc="left", fontsize=11); ax.set_xlabel("時刻 [ms]")
    if col == 0:
        for r_, lab in ((0, "① 加速度"), (1, "② 揺れ（比）"), (2, "③ 変位（比）")):
            fig.axes[r_].set_ylabel(lab)
ax = fig.add_subplot(gs[3, :2])
Tg = np.linspace(0, 0.1, 1000)
for f_, c in ((11.5, RED), (22.9, BLUE), (45.8, GREEN)):
    r = np.where(Tg == 0, 1, np.abs(np.sin(np.pi * f_ * Tg) / (np.pi * f_ * np.maximum(Tg, 1e-12))))
    ax.plot(Tg * 1000, r, color=c, lw=2.2, label=f"$f_n$ = {f_} Hz（周期 {1000/f_:.0f} ms）")
    for k in (1, 2, 3):
        if 1000 * k / f_ <= 100: ax.plot(1000 * k / f_, 0, "o", color=c, ms=8)
ax.axvline(20, color=GRAY, ls=":"); ax.axvline(44, color=GRAY, ls=":")
ax.set_xlabel("加速度の上昇時間 T（ジャーク立上げ時間） [ms]"); ax.set_ylabel("③で残る揺れの大きさ（段差 = 1）")
ax.set_title("④ ③の「残る揺れ」を上昇時間ごとに描いたもの = $|\\sin(\\pi f_n T)/(\\pi f_n T)|$（丸 = 周期の整数倍で 0）", loc="left", fontsize=12)
ax.grid(alpha=0.3); ax.legend(fontsize=10.5)
ax = fig.add_subplot(gs[3, 2]); ax.axis("off")
ax.text(0.02, 0.98, "読み方\n\n① 指令：加速度を 0→2 m/s² に上げる\n   上昇を 8 個の小さな段差で近似\n② 段差 1 つ = 本体を 22.9 Hz で揺らす\n③ 合計 = 実際の本体の動き\n   1 で止まれば揺れなし\n\n・上昇時間 44 ms：②の揺れが 1 周期に\n  均等にずれ、山と谷が打ち消す\n・上昇時間が短いと打ち消しきれない\n④ 固有振動数が高いほど、\n   短い上昇時間で 0 に届く",
        va="top", fontsize=12, linespacing=1.45)
fig.suptitle(f"なぜ加速度の上昇時間を周期 {P*1000:.0f} ms にすると揺れが消えるのか：加速の始まり 1 回に対する本体の変位 $x_b$（$f_n$ = {fn} Hz、減衰なし）", fontsize=13)
fig.tight_layout(); fig.savefig(D / "results" / "jerk_cancel.png"); print("ok")
