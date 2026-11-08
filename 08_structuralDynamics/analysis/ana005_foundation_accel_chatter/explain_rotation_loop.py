"""5.2 補足：正転・逆転で何が変わるか（揺れ → 力 → 揺れ の一周で、符号が 2 回反転して元に戻る）"""
import sys
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
from matplotlib.patches import FancyBboxPatch
fig, axs = plt.subplots(1, 2, figsize=(16, 5.4))


def box(ax, x, y, t, c):
    ax.add_patch(FancyBboxPatch((x - 1.25, y - 0.45), 2.5, 0.9, boxstyle="round,pad=0.05", fc=c, alpha=0.18, ec=c, lw=2)); ax.text(x, y, t, ha="center", va="center", fontsize=13)


def arr(ax, a, b, t, c="k", off=(0, 0.25)):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", lw=2.2, color=c, mutation_scale=20))
    ax.text((a[0] + b[0]) / 2 + off[0], (a[1] + b[1]) / 2 + off[1], t, ha="center", fontsize=12.5, color=c)


for ax, s, title in [(axs[0], "−", "正転"), (axs[1], "＋", "逆転")]:
    s2 = "＋" if s == "−" else "−"
    box(ax, 1.5, 3.2, "X に揺れる", "#1f618d"); box(ax, 6.5, 3.2, f"{s}Y 向きの力", "#c0392b")
    box(ax, 6.5, 0.9, f"{s}Y に揺れる", "#1f618d"); box(ax, 1.5, 0.9, "＋X 向きの力\n（元と同じ向き）", "#c0392b")
    arr(ax, (2.8, 3.2), (5.2, 3.2), f"$B_0$ の非対角：{s}1", "#c0392b")
    arr(ax, (6.5, 2.7), (6.5, 1.4), "$G_{yy}$\n（Y を押すと Y に揺れる）", "#1f618d", off=(1.6, -0.3))
    arr(ax, (5.2, 0.9), (2.8, 0.9), f"$B_0$ の非対角：{s2}1", "#c0392b", off=(0, -0.55))
    arr(ax, (1.5, 1.4), (1.5, 2.7), "$G_{xx}$", "#1f618d", off=(-0.6, -0.1))
    m1, m2 = ("-", "+") if s == "−" else ("+", "-")
    ax.text(4.0, 2.05, f"一周の符号\n$({m1}1)\\times({m2}1) = -1$", ha="center", va="center", fontsize=14, weight="bold")
    ax.set_xlim(-0.5, 9.2); ax.set_ylim(-0.3, 4.1); ax.axis("off"); ax.set_title(f"{title}：X → Y → X と一周すると", fontsize=15, loc="left")
fig.suptitle("正転と逆転で変わるのは「X に揺れたとき Y に出る力の向き」だけ。一周すると符号は 2 回反転するので、$G_{xy}=0$ の機械では結果は同じ", fontsize=13.5)
fig.tight_layout(); fig.savefig(R / "rotation_loop.png", dpi=95, bbox_inches="tight"); print("ok")
