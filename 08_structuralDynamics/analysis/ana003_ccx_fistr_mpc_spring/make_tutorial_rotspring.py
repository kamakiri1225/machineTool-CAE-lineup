"""ブログ用：並進ばねと回転ばねを分ける理由の図（tutorial/ccx と tutorial/ccx_no_rot の結果を使う）"""
import os, re
import numpy as np
from figstyle_local import plt
from matplotlib.patches import Rectangle, FancyArrowPatch
D = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(D, "tutorial"); OUT = os.path.join(D, "results", "tutorial")


def u(case, n):
    txt = open(os.path.join(T, case, "model.dat")).read()
    m = re.search(rf"^\s+{n}\s+(\S+)\s+(\S+)\s+(\S+)", txt, re.M); return np.array([float(m.group(k)) for k in (1, 2, 3)])


def spring(ax, p0, p1, n=6, w=0.12, **kw):
    p0, p1 = np.array(p0, float), np.array(p1, float); d = p1 - p0; L = np.linalg.norm(d); t = d / L; nr = np.array([-t[1], t[0]])
    pts = [p0, p0 + t * L * .2] + [p0 + t * L * s + nr * w * (1 if i % 2 == 0 else -1) for i, s in enumerate(np.linspace(.2, .8, 2 * n + 1)[1:-1])] + [p0 + t * L * .8, p1]
    pts = np.array(pts); ax.plot(pts[:, 0], pts[:, 1], **kw)


F, k, kth = 1.0, 1000.0, 1e6
u13 = u("ccx", 13)[1]; du = u("ccx", 201)[1] - u("ccx", 101)[1]; u13n = u("ccx_no_rot", 13)[1]
fig = plt.figure(figsize=(16, 6.6))
# 左：横から見た図（YZ 面）
ax = fig.add_axes([0.0, 0.04, 0.42, 0.86]); ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.35, 2.5); ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(Rectangle((-0.6, -1.1), 1.2, 1.0, fc="#d5d8dc", ec="k")); ax.text(0, -0.6, "下のブロック", ha="center")
th = 0.2; c, s_ = np.cos(th), np.sin(th); R = lambda x, y: (0.12 + x * c + y * s_, 0.25 - x * s_ + y * c)
ax.add_patch(plt.Polygon(np.array([R(-0.6, 0), R(0.6, 0), R(0.6, 1.0), R(-0.6, 1.0)]), fc="#aed6f1", ec="k"))
ax.text(*R(0, 0.45), "上のブロック", ha="center", va="center")
ax.add_patch(plt.Polygon(np.array([(-0.6, 0.0), (0.6, 0.0), (0.6, 1.0), (-0.6, 1.0)]) + (0, 0.25), fill=False, ec="#7f8c8d", ls="--"))
spring(ax, (0, -0.1), R(0, 0), color="#e74c3c", lw=2.5)
ax.plot(0, -0.1, "o", color="#f4d03f", ms=10, mec="k"); ax.plot(*R(0, 0), "o", color="#f4d03f", ms=10, mec="k")
p13 = R(-0.6, 1.0); ax.plot(*p13, "o", color="#f39c12", ms=11, mec="k")
ax.add_patch(FancyArrowPatch((p13[0] - 0.9, p13[1]), (p13[0] - 0.05, p13[1]), arrowstyle="-|>", mutation_scale=22, color="#f39c12", lw=3))
ax.text(p13[0] - 0.95, p13[1] + 0.12, "節点 13 に $F$ = 1 N", color="#d35400", fontsize=12)
ax.annotate("", xy=(1.05, -0.1), xytext=(1.05, 1.25), arrowprops=dict(arrowstyle="<->")); ax.text(1.12, 0.55, "$h$ = 100 mm", fontsize=12)
ax.text(0, 2.15, "代表点から離れた節点を押すと、力 × 腕の長さ の\nモーメントで上のブロックが傾く（ばねの周りに回る）", ha="center", fontsize=12, color="#c0392b")
ax.text(-2.1, -1.3, "横から見た図（Y が右、Z が上）。点線 = 元の位置", fontsize=10, color="#555")
# 右上：内訳
ax = fig.add_axes([0.5, 0.56, 0.47, 0.32])
tr, rot = F / k, u13 - F / k
ax.barh([0], [tr], color="#e74c3c", label=f"並進ばねの伸び $F/k$ = {tr:.4f} mm")
ax.barh([0], [rot], left=[tr], color="#58d68d", label=f"回転ばねによる傾き $\\theta \\times r$ = {rot:.4f} mm")
ax.set_yticks([]); ax.set_xlim(0, 0.016); ax.legend(fontsize=11, loc="lower right"); ax.grid(alpha=0.3, axis="x")
ax.set_title(f"節点 13 の Y 変位 {u13:.4f} mm の内訳（CalculiX、回転ばねあり）", loc="left", fontsize=12)
ax.set_xlabel("mm")
# 右下：回転ばね有無
ax = fig.add_axes([0.5, 0.08, 0.47, 0.3])
vals = [("回転ばね あり", u13, "#1f618d"), ("回転ばね なし", abs(u13n), "#7f8c8d")]
ax.barh([0, 1], [v[1] for v in vals], color=[v[2] for v in vals]); ax.set_xscale("log"); ax.set_xlim(1e-3, 1e12); ax.invert_yaxis()
ax.set_yticks([0, 1]); ax.set_yticklabels([v[0] for v in vals], fontsize=12)
ax.text(u13 * 1.8, 0, f"{u13:.4f} mm（手計算 0.0135 と一致）", va="center", fontsize=11)
ax.text(3e-3, 1, f"{abs(u13n):.1e} mm ＝ 止まらずに回る（でたらめ）", va="center", fontsize=11, color="white")
ax.set_xlabel("節点 13 の Y 変位 [mm]（対数目盛）"); ax.grid(alpha=0.3, axis="x")
ax.set_title("回転ばねを外すと、上のブロックが回転を止められない", loc="left", fontsize=12)
fig.savefig(os.path.join(OUT, "step5_why_rotspring.png")); print("ok", u13, du, u13n)
