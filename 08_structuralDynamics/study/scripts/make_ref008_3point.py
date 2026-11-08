"""ref008 6 章の図：3 点支持と 4 点支持の違い
   左：床の沈み δ がベッドをねじる量（ねじれモードへの射影による簡易モデル）
   右：上から見た支持点と、可動部が動いたときの重心の動き（マシニングセンタでの注意点）"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); os.makedirs(OUT, exist_ok=True)

fig, ax = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={"width_ratios": [1, 1.15, 1.15]})

# (a) 模式図：4 点では 1 か所の沈みがベッドのねじれになる
a = ax[0]; a.set_aspect("equal"); a.axis("off")
def bed(x0, y0, pts, title, sink, col):
    a.add_patch(plt.Rectangle((x0, y0), 4, 2.4, fc="#dde6f3", ec="k"))
    for (px, py) in pts:
        c = RED if (px, py) == sink else col
        a.add_patch(plt.Circle((x0 + px, y0 + py), 0.18, fc=c, ec="k", zorder=3))
    a.text(x0 + 2, y0 + 2.75, title, ha="center", fontsize=11)
bed(0, 3.2, [(0.3, 0.3), (3.7, 0.3), (0.3, 2.1), (3.7, 2.1)], "4 点支持（不静定）", (3.7, 0.3), GRAY)
bed(0, 0, [(0.3, 0.3), (0.3, 2.1), (3.7, 1.2)], "3 点支持（静定）", (3.7, 1.2), GRAY)
a.text(4.2, 3.4, "床が $\\delta$ 沈む\n→ 3 点では当たらない\n→ ベッドがねじられる", fontsize=9.5, color=RED, va="bottom")
a.text(4.2, 0.9, "床が $\\delta$ 沈む\n→ 全体が少し傾くだけ\n→ ベッドは変形しない", fontsize=9.5, color=GREEN, va="bottom")
a.set_xlim(-0.2, 7.6); a.set_ylim(-0.3, 6.4); a.set_title("(a) 支持点の数と床の沈み（上から見た図）", loc="left", fontsize=11)

# (b) 簡易モデル：ねじれ量 / δ
r = np.logspace(-2, 2, 200)            # k_f / K_w（足の剛性 / ベッドのねじれ剛性）
w4 = r / (1 + r)                         # 4 点：ベッドのねじれ（角の食い違い）/ δ
b = ax[1]
b.semilogx(r, w4, color=RED, lw=2.5, label="4 点支持")
b.semilogx(r, 0 * r, color=GREEN, lw=2.5, label="3 点支持（常に 0）")
b.axvspan(1, 100, color=ORANGE, alpha=0.12); b.text(1.3, 0.12, "硬い足・硬い床\n（精密機の据付）", fontsize=9.5, color=ORANGE)
b.set_xlabel("足の剛性 / ベッドのねじれ剛性  $k_f / K_w$")
b.set_ylabel("ベッドのねじれ / 床の沈み  $w/\\delta$")
b.set_ylim(-0.05, 1.05); b.grid(alpha=0.3); b.legend(loc="upper left")
b.set_title("(b) 床が $\\delta$ 沈んだときのベッドのねじれ\n$w/\\delta = (k_f/K_w)/(1+k_f/K_w)$（ねじれモードへの射影）", loc="left", fontsize=11)

# (c) マシニングセンタの注意点：可動部の重心が三角形の中にあるか
c = ax[2]; c.set_aspect("equal")
P = np.array([[0.3, 0.3], [0.3, 2.1], [3.7, 1.2], [0.3, 0.3]])
c.add_patch(plt.Rectangle((0, 0), 4, 2.4, fc="#dde6f3", ec="k"))
c.plot(P[:, 0], P[:, 1], color=BLUE, lw=2); c.scatter(P[:3, 0], P[:3, 1], s=120, color=GRAY, ec="k", zorder=3)
xs = np.linspace(0.8, 3.4, 30); ys = 1.2 + 0.0 * xs
c.plot(xs, ys, color=RED, lw=2, ls="--"); c.annotate("", xy=(3.4, 1.2), xytext=(3.0, 1.2), arrowprops=dict(arrowstyle="->", color=RED, lw=2))
c.text(0.8, 1.35, "X 軸を動かしたときの全体重心", color=RED, fontsize=9.5)
yy = np.linspace(0.6, 1.8, 2); c.plot([3.2, 3.2], [0.45, 1.95], color=ORANGE, lw=2, ls=":")
c.text(2.3, 2.2, "Y に振ると三角形の\n先端側で余裕が減る", color=ORANGE, fontsize=9.5, ha="center", va="bottom")
c.set_xlim(-0.2, 4.4); c.set_ylim(-0.3, 3.2); c.axis("off")
c.set_title("(c) 可動部の重心は支持三角形の中にあるか\n（転倒余裕・ロッキングの固有振動数が下がる）", loc="left", fontsize=11)

fig.tight_layout(); fig.savefig(os.path.join(OUT, "three_point.png"), dpi=90); print("ok")
