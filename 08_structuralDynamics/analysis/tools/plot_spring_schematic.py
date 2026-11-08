"""面と面をつなぐばねの割り振り方の模式図"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "study", "scripts"))
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, FancyArrowPatch


def main(path):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
    # (a) 断面: レールとキャリッジ
    ax = axs[0]; ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-6, 6); ax.set_ylim(-1, 11)
    ax.set_title("(a) 横から見た所：レール上面とキャリッジ下面を「ばね」でつなぐ", loc="left", fontsize=12)
    ax.add_patch(Rectangle((-5.5, 0), 11, 2, fc="#b0b7bf", ec="k")); ax.text(0, 1, "ベッド", ha="center", va="center")
    ax.add_patch(Rectangle((-1, 2), 2, 1.5, fc="#b0b7bf", ec="k")); ax.text(0, 2.75, "レール", ha="center", va="center", fontsize=10)
    ax.add_patch(Rectangle((-1, 4.5), 2, 1.5, fc="#aed6f1", ec="k")); ax.text(0, 5.25, "キャリッジ", ha="center", va="center", fontsize=10)
    ax.add_patch(Rectangle((-4.5, 6), 9, 1.5, fc="#aed6f1", ec="k")); ax.text(0, 6.75, "テーブル", ha="center", va="center")
    for x in (-1, -0.5, 0, 0.5, 1):
        zz = np.linspace(3.5, 4.5, 9); xx = x + 0.12 * np.array([0, 1, -1, 1, -1, 1, -1, 1, 0])
        ax.plot(xx, zz, color=BLUE, lw=1.5)
    ax.text(1.3, 4.0, "ばね（実際は長さ 0：\n同じ位置にある別々の節点を結ぶ）", color=BLUE, fontsize=10, va="center")
    ax.text(-5.8, 10.6, "実物ではレールとキャリッジの間に「ころ」が並び、\nその弾性変形がばねになる。モデルではこの層を\nばねに置き換えて、ころ・溝の形はモデル化しない。", fontsize=10.5, va="top")
    ax.text(-5.8, -0.6, "方向ごとに別のばね：Z（上下＝ラジアル）、Y（横）。X（送り方向）は転がるので 0", fontsize=10)

    # (b) 上から見た界面と割り振り
    ax = axs[1]; ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(b) 上から見た界面：面のばね定数 K を節点に割り振る", loc="left", fontsize=12)
    nx, ny, dx, dy = 4, 2, 20, 20
    for i in range(nx):
        for j in range(ny):
            ax.add_patch(Rectangle((i * dx, j * dy), dx, dy, fc="#fdebd0", ec="k", lw=1))
    W = np.zeros((nx + 1, ny + 1))
    for i in range(nx):
        for j in range(ny):
            for a, b in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                W[i + a, j + b] += 1 / 4 / (nx * ny)
    for i in range(nx + 1):
        for j in range(ny + 1):
            w = W[i, j]
            ax.plot(i * dx, j * dy, "o", ms=8 + 400 * w, color=RED, alpha=0.8, zorder=3)
            ax.text(i * dx + 1.2, j * dy + 1.5, f"{w*32:.0f}/32", fontsize=9.5, zorder=4)
    for i in range(nx):
        for j in range(ny + 1):
            ax.plot(i * dx + dx / 2, j * dy, "s", ms=5, color=GRAY, zorder=3)
    for i in range(nx + 1):
        for j in range(ny):
            ax.plot(i * dx, j * dy + dy / 2, "s", ms=5, color=GRAY, zorder=3)
    ax.text(0, -9, "キャリッジ 1個の下面（80 × 40 mm を 20 mm 角に分割）", fontsize=10.5)
    ax.text(0, -16, "赤丸 = ばねを付ける角の節点。数字は面全体の K に対する割合（合計 = 1）\n"
                    "　隅 1/32、辺 2/32、内部 4/32（各四角形の面積の 1/4 ずつを角に配る）\n"
                    "灰色 = 辺の中点の節点（20節点要素の中間節点）。ばねは付けない", fontsize=10.5, va="top")
    ax.set_xlim(-5, 110); ax.set_ylim(-33, 45)
    ax.text(85, 30, "$k_i = w_i \\, K$\n$\\sum_i w_i = 1$", fontsize=15)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1])
