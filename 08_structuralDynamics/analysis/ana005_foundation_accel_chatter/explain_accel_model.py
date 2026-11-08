"""ana005 4.1 の説明図：加減速反力の 1 自由度モデルで「何の変位を計算しているか」"""
import os, sys, json
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import Rectangle, FancyArrowPatch
acc = json.loads((D / "results" / "acceleration.json").read_text())
M, mt, k, c, fn, zeta = acc["M_kg"], acc["moving_mass_kg"], acc["k_N_m"], acc["c_Ns_m"], acc["fn_Hz"], acc["zeta"]


def spring(ax, x, y0, y1, w=0.12, n=5, **kw):
    ys = np.linspace(y0 + 0.08, y1 - 0.08, 2 * n + 1)
    xs = [x] + [x + (w if i % 2 else -w) for i in range(1, 2 * n)] + [x]
    ax.plot([x, x], [y0, y0 + 0.08], **kw); ax.plot(xs, ys, **kw); ax.plot([x, x], [y1 - 0.08, y1], **kw)


def damper(ax, x, y0, y1, **kw):
    ym = (y0 + y1) / 2
    ax.plot([x, x], [y0, ym - 0.05], **kw); ax.plot([x - 0.09, x + 0.09], [ym - 0.05, ym - 0.05], **kw)
    ax.plot([x - 0.12, x - 0.12, x + 0.12, x + 0.12], [ym + 0.12, ym - 0.12, ym - 0.12, ym + 0.12], **kw)
    ax.plot([x, x], [ym + 0.0, y1], **kw)


def arr(ax, a, b, c, lw=2.5):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=18, color=c, lw=lw))


fig, axs = plt.subplots(1, 2, figsize=(16, 7.2), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axs[0]; ax.set_xlim(-0.3, 7.2); ax.set_ylim(-0.9, 5.4); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("横から見た図：何がどう動くか（X = 送り方向）", loc="left", fontsize=13)
ax.plot([-0.2, 7.0], [0, 0], "k", lw=2)
for x in np.arange(-0.2, 7.0, 0.25): ax.plot([x, x - 0.15], [0, -0.15], "k", lw=0.8)
ax.text(6.9, -0.4, "床（動かない）", ha="right", fontsize=11)
# 脚：横ばね（X）を表すため、横向きに描く
for xf in (1.0, 5.0):
    ax.add_patch(Rectangle((xf - 0.35, 0.0), 0.7, 0.18, fc="#7f8c8d"))
    ax.plot([xf - 0.35, xf - 0.35], [0.18, 0.62], color="#555")
    ys = np.linspace(0.35, 0.35, 1)
    # 横ばね：床側の柱から本体へ
    xs = np.linspace(xf - 0.35, xf + 0.15, 11); yy = 0.45 + np.array([0, .08, -.08, .08, -.08, .08, -.08, .08, -.08, .08, 0])
    ax.plot(xs, yy, color=BLUE, lw=2)
    ax.text(xf, -0.7, "脚（4 か所）\n横の剛性 k, 減衰 c", ha="center", fontsize=9.5, color=BLUE)
ax.add_patch(Rectangle((0.6, 0.62), 5.2, 1.0, fc="#d5d8dc", ec="k")); ax.text(1.3, 1.12, "ベッド", ha="center", fontsize=11)
ax.add_patch(Rectangle((4.2, 1.62), 1.2, 3.2, fc="#d5d8dc", ec="k")); ax.text(4.8, 3.6, "コラム", ha="center", fontsize=11)
ax.add_patch(Rectangle((3.0, 3.9), 1.2, 0.6, fc="#d5d8dc", ec="k")); ax.text(3.6, 4.2, "主軸頭", ha="center", fontsize=9)
ax.add_patch(Rectangle((1.2, 1.62), 2.2, 0.55, fc="#f5cba7", ec="k")); ax.text(2.3, 1.9, f"テーブル＋ワーク  m_t = {mt:.0f} kg", ha="center", fontsize=9.5)
ax.text(0.65, 3.25, f"機械本体（ベッド＋コラム＋…）\n全質量 M = {M:.0f} kg（テーブルを含む）", fontsize=10.5)
# 力
arr(ax, (2.3, 2.45), (3.3, 2.45), ORANGE); ax.text(2.3, 2.55, "テーブルを押す力 = $m_t\\,\\ddot r$", color=ORANGE, fontsize=10.5)
arr(ax, (3.9, 1.35), (2.9, 1.35), RED); ax.text(3.95, 0.95, "反力 $-m_t\\,\\ddot r$\n（ボールねじ→ベッド）", color=RED, fontsize=10.5)
# 変位
ax.plot([5.8, 5.8], [0.62, 5.0], color=GREEN, ls=":", lw=1.5); arr(ax, (5.8, 4.95), (6.5, 4.95), GREEN, 2)
ax.text(0.2, 5.1, "$x_b$：本体全体の X 変位（床に対して）← 下段のグラフはこれ", color=GREEN, fontsize=11, weight="bold")
ax.plot([1.2, 1.2], [2.17, 2.6], color=ORANGE, ls=":"); arr(ax, (1.2, 2.32), (1.9, 2.32), ORANGE, 1.6)

ax.text(-0.25, 2.55, "$r(t)$：本体に対する\nテーブルの送り（指令）", color=ORANGE, fontsize=10)

ax = axs[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("運動方程式の読み方", loc="left", fontsize=13)
ax.text(0.2, 8.9, "$M\\,\\ddot x_b + c\\,\\dot x_b + k\\,x_b = -m_t\\,\\ddot r(t)$", fontsize=20)
rows = [("$M\\,\\ddot x_b$", "本体全体（テーブル込み）の慣性", GRAY),
        ("$c\\,\\dot x_b$", "脚の減衰（$c = \\alpha M + \\beta k$）", BLUE),
        ("$k\\,x_b$", "脚の横ばね（4 脚 × 10 N/µm = 40 N/µm）", BLUE),
        ("$-m_t\\,\\ddot r$", "テーブルを加速した反力が本体を逆向きに押す", RED)]
for i, (a, b, col) in enumerate(rows):
    y = 7.6 - i * 1.05
    ax.text(0.4, y, a, fontsize=15, color=col); ax.text(2.6, y + 0.05, b, fontsize=12)
ax.text(0.2, 3.1, f"固有振動数 $f_n = \\dfrac{{1}}{{2\\pi}}\\sqrt{{k/M}}$ = {fn:.1f} Hz、減衰比 ζ = {zeta*100:.1f} %", fontsize=12.5)
ax.text(0.2, 2.1, "・テーブルが加速・減速するたびに、本体が逆向きに押され、脚の上で揺れる", fontsize=11.5)
ax.text(0.2, 1.4, "・止まった後も、この 23 Hz の揺れが残る（＝下段のグラフの残留振動）", fontsize=11.5)
ax.text(0.2, 0.5, "※ 刃先とワークの「相対」変位ではない。本体の剛体並進だけの近似（傾き・部品の変形・サーボは含まない）", fontsize=10.5, color=RED)
fig.tight_layout(); fig.savefig(D / "results" / "accel_model_schematic.png"); print("ok", M, mt, k, fn, zeta)
