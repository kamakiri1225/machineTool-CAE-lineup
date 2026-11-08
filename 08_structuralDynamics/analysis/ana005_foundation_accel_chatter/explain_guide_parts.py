"""2.1 補足：リニアガイド（レール・キャリッジ）とボールねじ（ねじ軸・ナット）の実物と、モデルの部品の対応"""
import sys
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle, FancyBboxPatch
fig, axs = plt.subplots(1, 3, figsize=(18, 6.2), gridspec_kw=dict(width_ratios=[1, 1, 1.05]))
ax = axs[0]; ax.imshow(mpimg.imread(D / "photos" / "linear_guide_Mono_Rail.jpg")); ax.axis("off")
ax.set_title("(a) リニアガイド：レールの上をキャリッジが滑る", fontsize=14.5, loc="left")
ax.annotate("キャリッジ（ブロック）\n→ テーブルの裏にボルト止め", xy=(0.36, 0.78), xycoords="axes fraction", xytext=(0.02, 0.97), textcoords="axes fraction",
            fontsize=12.5, color="#b9770e", va="top", arrowprops=dict(arrowstyle="->", color="#b9770e", lw=2))
ax.annotate("レール\n→ ベッドの上にボルト止め", xy=(0.62, 0.42), xycoords="axes fraction", xytext=(0.55, 0.12), textcoords="axes fraction",
            fontsize=12.5, color="#2e4053", arrowprops=dict(arrowstyle="->", color="#2e4053", lw=2))
ax.text(0.99, 0.01, "出典: Wikimedia Commons, Rollon91 / CC BY-SA 3.0", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="white", bbox=dict(fc="black", alpha=0.5, lw=0))
ax = axs[1]; ax.imshow(mpimg.imread(D / "photos" / "ball_screw_BallScrews.jpg")); ax.axis("off")
ax.set_title("(b) ボールねじ：ねじ軸が回るとナットが進む", fontsize=14.5, loc="left")
ax.annotate("ナット\n→ ナットブラケットで\n　 テーブルに固定", xy=(0.62, 0.84), xycoords="axes fraction", xytext=(0.66, 0.62), textcoords="axes fraction",
            fontsize=12, color="#b9770e", bbox=dict(fc="white", alpha=0.9, ec="none"), arrowprops=dict(arrowstyle="->", color="#b9770e", lw=2))
ax.annotate("ねじ軸\n→ 両端を軸受でベッドに支持、\n　 モータで回す", xy=(0.3, 0.80), xycoords="axes fraction", xytext=(0.03, 0.62), textcoords="axes fraction",
            fontsize=12, color="#2e4053", bbox=dict(fc="white", alpha=0.9, ec="none"), arrowprops=dict(arrowstyle="->", color="#2e4053", lw=2))
ax.text(0.99, 0.01, "出典: Wikimedia Commons, Graibeard / CC BY-SA 3.0", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="white", bbox=dict(fc="black", alpha=0.5, lw=0))
ax = axs[2]; ax.set_xlim(0, 10); ax.set_ylim(0, 7.5); ax.axis("off"); ax.set_aspect("equal")
ax.set_title("(c) モデルでの置き換え（送り方向から見た断面）", fontsize=14.5, loc="left")
ax.add_patch(Rectangle((0.3, 0.3), 9.4, 2.0, fc="#7b7d7d", ec="k")); ax.text(5, 0.9, "ベッド", ha="center", color="white", fontsize=14)
for x in (1.4, 7.6):
    ax.add_patch(Rectangle((x, 2.3), 1.0, 0.7, fc="#2e4053", ec="k")); ax.text(x + 0.5, 2.55, "レール", ha="center", color="white", fontsize=10)
    ax.add_patch(Rectangle((x - 0.2, 3.25), 1.4, 0.75, fc="#f1c40f", ec="k")); ax.text(x + 0.5, 3.48, "キャリッジ", ha="center", fontsize=10)
    ax.plot([x + 0.5] * 2, [3.0, 3.25], color="#1f618d", lw=4)
ax.add_patch(Rectangle((4.4, 2.3), 1.2, 0.7, fc="#e67e22", ec="k")); ax.text(5, 2.55, "ナット受け台", ha="center", fontsize=9.5)
ax.add_patch(Rectangle((4.4, 3.25), 1.2, 0.75, fc="#ca6f1e", ec="k")); ax.text(5, 3.48, "ブラケット", ha="center", fontsize=9.5)
ax.plot([5, 5], [3.0, 3.25], color="#1f618d", lw=4)
ax.add_patch(Rectangle((0.6, 4.0), 8.8, 1.2, fc="#2471a3", ec="k")); ax.text(5, 4.45, "テーブル", ha="center", color="white", fontsize=14)
ax.text(0.3, 6.9, "青い線＝長さ 0 のばね（レールとキャリッジの間、\n受け台とブラケットの間）。ねじ軸と軸受は描かず、\n「ナット受け台」とばね 1 本でまとめて表している", fontsize=11.5, va="top")
fig.tight_layout(); fig.savefig(R / "guide_parts.png", dpi=95, bbox_inches="tight"); print("ok")
