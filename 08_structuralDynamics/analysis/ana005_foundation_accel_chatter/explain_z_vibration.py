"""5.2 補足：「Z に揺れる」とは何か。工具が主軸方向（上下）にずれたとき、削り厚さ h が変わるかどうかを横から見た図"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
from matplotlib.patches import Rectangle, Polygon, Circle

GRAY, ORANGE, BLUE, RED = "#bdc3c7", "#f5b041", "#1f618d", "#c0392b"
dz = 0.45
fig, axs = plt.subplots(1, 3, figsize=(15, 5.6))


def zarrow(ax, x, y):
    ax.annotate("", xy=(x, y - dz), xytext=(x, y), arrowprops=dict(arrowstyle="->", color=RED, lw=2.5))
    ax.text(x + 0.12, y - dz / 2, "Z に\n$\\delta_z$ ずれる", color=RED, fontsize=12, va="center")


# (a) 平エンドミルの側面削り
ax = axs[0]
ax.add_patch(Rectangle((0.6, -0.8), 2.4, 3.6, fc=GRAY, ec="k"))                 # ワーク
ax.add_patch(Rectangle((0.0, -0.8), 0.6, 3.6, fc=ORANGE, ec="k", hatch="//"))    # 削る層
ax.add_patch(Rectangle((-1.6, -0.6), 1.6, 3.4, fc="#5d6d7e", ec="k", alpha=0.9))   # 工具（実線）
ax.add_patch(Rectangle((-1.6, -0.6 - dz), 1.6, 3.4, fc="none", ec=BLUE, ls="--", lw=2))
ax.annotate("", xy=(0.0, 1.0), xytext=(0.6, 1.0), arrowprops=dict(arrowstyle="<->", lw=2)); ax.text(0.3, 1.15, "$h$", ha="center", fontsize=16)
zarrow(ax, -0.8, 2.4)
ax.set_title("(a) 今回：直刃の平エンドミルで溝の側面を削る", fontsize=15, loc="left")
ax.text(-1.7, -1.6, "刃は側面で削る。上下にずれても\n削る厚さ $h$（横の幅）は変わらない\n→ Z は無視できる（2×2 でよい）", fontsize=12.5, va="top", color=BLUE)

# (b) ボールエンドミル：球の下面で削る
ax = axs[1]
xs = np.linspace(-2.2, 2.6, 2)
ax.add_patch(Rectangle((-2.2, -1.2), 4.8, 1.2, fc=GRAY, ec="k"))
ax.add_patch(Rectangle((0.0, -0.45), 2.6, 0.45, fc=ORANGE, ec="k", hatch="//"))       # 削る層（厚さ h）
ax.add_patch(Rectangle((-0.35, 1.0), 0.7, 2.2, fc="#5d6d7e", ec="k"))
ax.add_patch(Circle((0, 0.55), 1.0, fc="#5d6d7e", ec="k")); ax.add_patch(Rectangle((-1.0, 0.55), 2.0, 0.6, fc="#5d6d7e", ec="none"))
ax.add_patch(Circle((0, 0.55 - dz), 1.0, fc="none", ec=BLUE, ls="--", lw=2))
ax.plot([-1.0, 1.0], [-0.45, -0.45], color="k", lw=1)
ax.annotate("", xy=(1.9, -0.45), xytext=(1.9, 0.0), arrowprops=dict(arrowstyle="<->", lw=2)); ax.text(2.0, -0.25, "$h$", fontsize=16, va="center")
ax.annotate("", xy=(1.3, -0.45 - dz), xytext=(1.3, -0.45), arrowprops=dict(arrowstyle="<->", color=RED, lw=2)); ax.text(1.4, -0.7, "$+\\delta_z$", color=RED, fontsize=14)
zarrow(ax, -0.9, 2.6)
ax.set_title("(b) ボールエンドミル（先端の球で削る）", fontsize=15, loc="left")
ax.text(-2.2, -1.6, "球の下側で削っているので、下にずれると\nその分だけ深く削る → $h$ が $\\delta_z$ 変わる\n→ Z の揺れも再生びびりに入る（3×3）", fontsize=12.5, va="top", color=RED)

# (c) 正面フライス：45° の刃で床を削る
ax = axs[2]
ax.add_patch(Rectangle((-2.2, -1.2), 4.8, 1.2, fc=GRAY, ec="k"))
ax.add_patch(Polygon([(0.0, -0.6), (0.6, 0.0), (2.6, 0.0), (2.6, -0.6)], fc=ORANGE, ec="k", hatch="//"))  # 削る層
P0 = np.array([(-2.0, 2.8), (-2.0, -0.6), (0.0, -0.6), (0.9, 0.3), (0.9, 1.0), (-1.0, 1.0), (-1.0, 2.8)])
ax.add_patch(Polygon(P0, fc="#5d6d7e", ec="k"))
ax.add_patch(Polygon(P0 - [0, dz], fc="none", ec=BLUE, ls="--", lw=2))
ax.annotate("", xy=(1.9, -0.6), xytext=(1.9, 0.0), arrowprops=dict(arrowstyle="<->", lw=2)); ax.text(2.0, -0.3, "$a_p$", fontsize=16, va="center")
ax.annotate("", xy=(0.25, -0.6 - dz), xytext=(0.25, -0.6), arrowprops=dict(arrowstyle="<->", color=RED, lw=2)); ax.text(0.35, -0.85, "$+\\delta_z$", color=RED, fontsize=14)
zarrow(ax, -1.6, 2.4)
ax.set_title("(c) 正面フライス（底の斜めの刃で削る）", fontsize=15, loc="left")
ax.text(-2.2, -1.6, "斜めの刃が床を削っているので、下にずれると\n刃が当たる厚さが変わる → $h$ が変わる\n→ Z の揺れも再生びびりに入る（3×3）", fontsize=12.5, va="top", color=RED)

for ax in axs:
    ax.set_xlim(-2.4, 3.1); ax.set_ylim(-2.9, 3.4); ax.set_aspect("equal"); ax.axis("off")
    ax.text(2.9, 3.1, "灰：ワーク\n橙：削る層\n実線：元の工具\n青破線：Z にずれた工具", fontsize=9.5, ha="right", va="top", color="#555")
fig.suptitle("「Z に揺れる」＝工具が主軸方向（上下）にずれること（横から見た図）", fontsize=14, y=1.03)
fig.tight_layout(); fig.savefig(R / "z_vibration.png", dpi=100, bbox_inches="tight"); print("ok")
