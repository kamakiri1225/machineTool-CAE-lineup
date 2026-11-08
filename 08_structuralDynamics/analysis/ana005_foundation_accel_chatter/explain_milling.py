"""ana005 5 章の説明図：①切削条件が何を指すか（絵）、②安定限界線図の計算フロー"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Wedge, FancyBboxPatch, Polygon


def arr(ax, a, b, c, lw=2.2, ms=16):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=ms, color=c, lw=lw))


fig = plt.figure(figsize=(20, 6.8))
# (a) 横から：エンドミルと軸方向切込み ap
ax = fig.add_axes([0.0, 0.0, 0.30, 0.95]); ax.set_xlim(-3, 3); ax.set_ylim(-2.2, 4.3); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(a) エンドミルで溝を削る（横から）", loc="left", fontsize=18)
ax.add_patch(Rectangle((-3, -2.2), 6, 2.2, fc="#d5d8dc", ec="k")); ax.text(2.2, -1.9, "ワーク", fontsize=18)
ax.add_patch(Rectangle((-0.5, -0.8), 1.0, 4.0, fc="#aab7b8", ec="k"))
for y in np.linspace(-0.7, 2.9, 9): ax.plot([-0.5, 0.5], [y, y + 0.35], color="k", lw=1)
ax.add_patch(Rectangle((-0.5, -0.8), 1.0, 0.8, fc=RED, alpha=0.25))
ax.annotate("", xy=(0.75, 0), xytext=(0.75, -0.8), arrowprops=dict(arrowstyle="<->", color=RED, lw=2)); ax.text(0.85, -0.5, "$a_p$ 軸方向切込み\n（線図の縦軸）", color=RED, fontsize=18)
ax.annotate("", xy=(-0.65, 3.6), xytext=(-0.65, 2.8), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2, connectionstyle="arc3,rad=-1.2"))
ax.text(-2.95, 2.0, "回転数 $n$\n500〜16,000 min$^{-1}$\n（線図の横軸）", color=BLUE, fontsize=18)
arr(ax, (-2.8, 0.5), (-1.0, 0.5), GREEN); ax.text(-2.9, 0.75, "送り方向", color=GREEN, fontsize=18)
# (b) 上から：フルスロット、4 刃、接触角
ax = fig.add_axes([0.31, 0.0, 0.34, 0.95]); ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.4, 2.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(b) フルスロット・4 刃・接触角 0〜π（上から）", loc="left", fontsize=18)
ax.add_patch(Rectangle((-2.6, -1.0), 2.6, 2.0, fc="#d5d8dc", ec="none")); ax.text(-2.4, 1.15, "削った溝（工具径と同じ幅）", fontsize=15)
ax.add_patch(Rectangle((0, -2.4), 2.6, 1.4, fc="#d5d8dc", ec="none")); ax.add_patch(Rectangle((0, 1.0), 2.6, 1.4, fc="#d5d8dc", ec="none"))
ax.add_patch(Rectangle((0.0, -1.0), 2.6, 2.0, fc="#f2d7d5", ec="none")); ax.text(1.85, -0.3, "これから\n削る所", ha="center", fontsize=15, color=RED)
ax.add_patch(Wedge((0, 0), 1.0, -90, 90, fc=ORANGE, alpha=0.35)); ax.text(0.05, -1.35, "刃が削る範囲 0〜π（半周）", fontsize=15, color=ORANGE)
ax.add_patch(Circle((0, 0), 1.0, fill=False, ec="k", lw=1.5))
for k in range(4):
    th = np.radians(20 + 90 * k); p = np.array([np.cos(th), np.sin(th)])
    ax.plot([0, p[0]], [0, p[1]], color="k", lw=2); ax.plot(*p, "o", color=RED if k == 0 else "k", ms=7)
ax.text(-1.15, -1.55, "刃数 Z = 4（等ピッチ 90°）", fontsize=15.5)
p0 = np.array([np.cos(np.radians(20)), np.sin(np.radians(20))])
arr(ax, p0, p0 + 0.8 * np.array([-p0[1], p0[0]]), RED); ax.text(*(p0 + [-0.45, 1.0]), "$F_t = K_t\\,a_p\\,h$  接線力", color=RED, fontsize=18); ax.text(-2.5, 2.2, "$K_t$：削る断面 1 mm² あたりの力（材料で決まる係数）\n$h$：1 枚の刃が削る厚さ", color=RED, fontsize=14)
arr(ax, p0, p0 + 0.4 * p0, BLUE); ax.text(*(p0 + [0.5, 0.0]), "$F_r = K_r F_t$\n半径方向力", color=BLUE, fontsize=18)
arr(ax, (-2.4, -1.9), (-1.2, -1.9), GREEN); ax.text(-2.4, -1.75, "送り", color=GREEN, fontsize=18)
# (c) 1 枚の刃が削る厚さ h
ax = fig.add_axes([0.66, 0.0, 0.34, 0.95]); ax.set_xlim(-1.9, 2.4); ax.set_ylim(-1.9, 2.2); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(c) 1 枚の刃が削る厚さ $h$（上から）", loc="left", fontsize=18)
fz = 0.45; th = np.linspace(-np.pi / 2, np.pi / 2, 200)
ax.plot(np.cos(th) * 1.3, np.sin(th) * 1.3, color=GRAY, lw=2, ls="--")
ax.plot(np.cos(th) * 1.3 + fz, np.sin(th) * 1.3, color="k", lw=2.5)
xs = np.concatenate([np.cos(th) * 1.3, (np.cos(th) * 1.3 + fz)[::-1]]); ys = np.concatenate([np.sin(th) * 1.3, np.sin(th)[::-1] * 1.3])
ax.fill(xs, ys, color=ORANGE, alpha=0.35)
ax.text(-1.85, 1.75, "点線：前の刃が通った跡\n実線：今の刃が通る跡\n（工具は 1 刃ぶん $f_z$ だけ送られている）", fontsize=14)
phi = np.radians(55); a0 = np.array([np.cos(phi - np.pi / 2), np.sin(phi - np.pi / 2)])
ax.plot([0, a0[0] * 1.3], [0, a0[1] * 1.3], color=GRAY, lw=1)
p1 = a0 * 1.3; p2 = p1 + a0 * fz * np.sin(phi)
ax.annotate("", xy=p2 + a0 * 0.02, xytext=p1, arrowprops=dict(arrowstyle="<->", color=RED, lw=2.5))
ax.text(p2[0] + 0.1, p2[1] - 0.25, "$h = f_z \\sin\\phi$", color=RED, fontsize=19)
ax.annotate("", xy=(fz, -1.65), xytext=(0, -1.65), arrowprops=dict(arrowstyle="<->", color=GREEN, lw=2))
ax.text(fz + 0.08, -1.75, "$f_z$：1 刃あたりの送り", color=GREEN, fontsize=15)
ax.text(0.35, -0.35, "$\\phi$", fontsize=18)
ax.text(-1.85, -1.35, "揺れがあると、$h$ に\n今と前の揺れの差が加わる\n→ 動的切削力（5.2 節）", fontsize=14, color=BLUE)
fig.savefig(D / "results" / "milling_explained.png")
# 条件の一覧（文字を大きく）
rows = [("刃数 $Z$ = 4、等ピッチ", "刃が 4 枚、90° ずつ均等。1 回転で 4 回削る"),
        ("フルスロット（接触角 0〜π）", "工具径と同じ幅の溝を掘る。刃は半周ずっと削る（いちばんびびりやすい条件）"),
        ("$K_t$ = 1,000 N/mm²", "削る断面 1 mm² あたりの接線方向の力。$F_t = K_t\\,a_p\\,h$。材料で決まる（鋼 2,000、アルミ 700 程度）"),
        ("$K_r$ = 0.3", "半径方向の力 ÷ 接線方向の力。$F_r = 0.3\\,F_t$"),
        ("回転数 $n$ = 500〜16,000 min⁻¹", "線図の横軸。刃が通る周波数 = $Zn/60$（6,000 min⁻¹ で 400 Hz）"),
        ("びびり周波数 5〜3,000 Hz", "びびりが起こりうる振動の周波数を探す範囲。FRF をここまで使う")]
fig, ax = plt.subplots(figsize=(17.5, 6.2)); ax.axis("off"); ax.set_xlim(0, 17.5); ax.set_ylim(0, 6.2)
ax.text(0.1, 5.85, "計算に使った切削条件", fontsize=20, weight="bold")
for i, (a_, b_) in enumerate(rows):
    y = 5.0 - i * 0.92
    ax.add_patch(FancyBboxPatch((0.1, y - 0.36), 17.3, 0.74, boxstyle="round,pad=0.02", fc=BLUE if i % 2 == 0 else ORANGE, alpha=0.07, ec="none"))
    ax.text(0.3, y, a_, fontsize=17, color=BLUE, weight="bold", va="center"); ax.text(6.0, y, b_, fontsize=15, va="center")
fig.savefig(D / "results" / "milling_conditions.png"); 
# 2) 計算フロー（縦に並べ、数式を大きく）
steps = [("① FEM で相対 FRF を求める", "固有値（80 モード）から刃先―ワーク相対 FRF の 2×2 行列をつくる",
          ("MAT", r"$G(\omega)=$", [[r"$G_{xx}$", r"$G_{xy}$"], [r"$G_{yx}$", r"$G_{yy}$"]]), BLUE),
         ("② 方向係数（零次近似）", "刃の向きで変わる力の方向を、刃が削る半周で平均する",
          ("MAT", r"$B_0=\dfrac{Z}{4}$", [[r"$K_r$", r"$1$"], [r"$-1$", r"$K_r$"]]), ORANGE),
         ("③ びびり周波数を 1 つ選ぶ", "5〜3,000 Hz の中の $\\omega_c$ で、$G(\\omega_c)B_0$ の固有値 $\\lambda$ を求める",
          r"$1+a_pK_t\left(1-e^{-i\omega_c\tau}\right)\lambda=0$", GREEN),
         ("④ 限界切込みを求める", "$\\mu=-1/\\lambda$ とおく。$\\mu$ の実部が正のときだけ、びびりが起こりうる",
          r"$a_{\lim}=\dfrac{\mu_R}{2K_t}\left(1+\kappa^2\right),\quad \kappa=\dfrac{\mu_I}{\mu_R}$", RED),
         ("⑤ 回転数を求める", "位相の条件から回転数が決まる（k = 0, 1, 2, … がローブの番号）",
          r"$n=\dfrac{60\,\omega_c}{Z\,(2\pi k+\varepsilon)},\quad \varepsilon=\pi-2\arctan\kappa$", GRAY),
         ("⑥ 線図をつくる", "③〜⑤を全ての $\\omega_c$ で繰り返し、各回転数で最小の $a_{\\lim}$ をつなぐ（下段：そのときの $\\omega_c$）",
          r"$a_{\lim}(n)=\min_{\omega_c,\,k}\ a_{\lim}$", BLUE)]
fig, ax = plt.subplots(figsize=(13, 8.6)); ax.set_xlim(0, 13); ax.set_ylim(-0.15, 8.6); ax.axis("off")
ax.text(0.1, 8.3, "安定限界線図の計算フロー", fontsize=17, weight="bold")
for i, (t, dsc, eq, c) in enumerate(steps):
    y = 7.45 - i * 1.38
    ax.add_patch(FancyBboxPatch((0.1, y - 0.55), 12.8, 1.1, boxstyle="round,pad=0.04", fc=c, alpha=0.10, ec=c, lw=1.5))
    ax.text(0.3, y + 0.25, t, fontsize=14, weight="bold", color=c, va="center")
    ax.text(0.3, y - 0.25, dsc, fontsize=11, va="center")
    if isinstance(eq, tuple):
        _, pre, M = eq; tx = ax.text(7.6, y, pre, fontsize=16, va="center")
        fig.canvas.draw(); bb = tx.get_window_extent().transformed(ax.transData.inverted()); x0 = bb.x1 + 0.1
        for r, row in enumerate(M):
            for cidx, v in enumerate(row):
                ax.text(x0 + 0.2 + cidx * 0.85, y + 0.22 - r * 0.44, v, fontsize=15, va="center")
        for xb, d in ((x0, 1), (x0 + 1.85, -1)):
            ax.plot([xb + 0.1 * d, xb, xb, xb + 0.1 * d], [y + 0.45, y + 0.45, y - 0.45, y - 0.45], color="k", lw=1.5)
    else:
        ax.text(7.6, y, eq, fontsize=16, va="center")
    if i < 5: arr(ax, (1.0, y - 0.58), (1.0, y - 0.8), "k", 1.5, 12)
fig.savefig(D / "results" / "milling_flow.png", bbox_inches="tight", dpi=100); print("ok")
