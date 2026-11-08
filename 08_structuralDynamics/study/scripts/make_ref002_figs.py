"""ref002 の図を生成する: python make_ref002_figs.py"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, Polygon, Circle, FancyArrowPatch, Arc

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref002")
os.makedirs(OUT, exist_ok=True)


def arrow(ax, p, q, color=BLUE, style="<|-|>", lw=1.5, ms=12):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, color=color, lw=lw, mutation_scale=ms))


# ---------------------------------------------------------------- 切削条件
def cutting_params():
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), gridspec_kw=dict(width_ratios=[1, 1.1, 1]))
    for ax in axs:
        ax.set_aspect("equal"); ax.axis("off")

    # ① 送り方向から見た図
    ax = axs[0]
    ax.set_xlim(0, 10); ax.set_ylim(-1.5, 10)
    ax.set_title("① 送り方向から見た図：切込み $a_p$, $a_e$", loc="left", fontsize=13)
    ax.add_patch(Polygon([(0.5, 6), (5, 6), (5, 2.5), (9, 2.5), (9, 0), (0.5, 0)], fc="#d5d8dc", ec="#566573"))
    ax.add_patch(Rectangle((5, 2.5), 1.4, 3.5, fc=ORANGE, alpha=0.45, ec=ORANGE, ls="--"))
    ax.add_patch(Rectangle((5, 6), 2.4, 3.3, fc="#eaecee", ec="#566573"))          # シャンク
    ax.add_patch(Rectangle((5, 2.5), 2.4, 3.5, fc="#aab7b8", ec="#566573", alpha=0.9))  # 刃部
    for y in (3.0, 4.1, 5.2):
        ax.plot([5, 7.4], [y, y + 0.6], color="#566573", lw=1)
    ax.plot([6.2, 6.2], [1.8, 9.8], color=GRAY, ls="-.", lw=0.8)
    ax.text(7.6, 8.4, "エンドミル\n（回転 $n$）", fontsize=11, va="center")
    arrow(ax, (4.6, 2.5), (4.6, 6))
    ax.plot([4.4, 5], [6, 6], color=BLUE, ls=":", lw=1)
    ax.text(0.7, 4.2, "$a_p$\n工具の軸方向の\n切込み（深さ）", color=BLUE, fontsize=10.5, va="center")
    ax.plot([5, 5], [-0.2, -0.9], color=BLUE, ls=":"); ax.plot([6.4, 6.4], [2.4, -0.9], color=BLUE, ls=":")
    arrow(ax, (5, -0.7), (6.4, -0.7))
    ax.text(6.6, -0.7, "$a_e$", color=BLUE, fontsize=13, va="center")
    ax.text(0.5, -1.5, "$a_e$：工具の径方向の切込み（かぶり幅）", color=BLUE, fontsize=10.5, va="center")
    ax.text(0.7, 1.0, "オレンジ ＝ いま削っている断面\n（面積 $a_p \\times a_e$）", color="#a04000", fontsize=10.5, va="center")

    # ② 上から見た図
    ax = axs[1]
    ax.set_xlim(0, 11); ax.set_ylim(-1.5, 10)
    ax.set_title("② 上から見た図：送り速度 $V_f$ と除去量 $Q$", loc="left", fontsize=13)
    ax.add_patch(Polygon([(0.5, 5), (6, 5), (6, 6), (10, 6), (10, 0), (0.5, 0)], fc="#d5d8dc", ec="#566573"))
    ax.add_patch(Rectangle((0.5, 5), 5.5, 1, fc=ORANGE, alpha=0.35, ec=ORANGE, ls="--"))
    ax.add_patch(Circle((6, 6.5), 1.5, fc="#aab7b8", ec="#566573", alpha=0.95, zorder=3))
    ax.add_patch(Arc((6, 6.5), 2.2, 2.2, theta1=30, theta2=150, color=RED, lw=2, zorder=4))
    arrow(ax, (6.9, 7.1), (7.0, 6.95), color=RED, style="-|>", ms=14)
    ax.text(5.8, 7.8, "$n$", color=RED, fontsize=14, zorder=5)
    arrow(ax, (8.0, 7.2), (10.6, 7.2), color=RED, style="-|>", lw=3, ms=18)
    ax.text(8.0, 7.6, "$V_f$ [mm/min]\n工具が進む速さ", color=RED, fontsize=11)
    arrow(ax, (0.5, 4.3), (6, 4.3), color=ORANGE)
    ax.text(0.8, 3.4, "1分間に進む長さ ＝ $V_f$", color=ORANGE, fontsize=11)
    arrow(ax, (10.4, 5), (10.4, 6))
    ax.plot([10, 10.6], [5, 5], color=BLUE, ls=":"); ax.text(10.6, 5.3, "$a_e$", color=BLUE, fontsize=13)
    ax.text(0.5, -0.6, "1分間に削れる体積（材料除去率）", fontsize=11)
    ax.text(0.5, -1.5, "$Q = a_p \\times a_e \\times V_f$　＝ オレンジの断面 × 進んだ長さ", fontsize=12)

    # ③ 比切削抵抗
    ax = axs[2]
    ax.set_xlim(0, 10); ax.set_ylim(-1.5, 10)
    ax.set_title("③ 比切削抵抗 $k_c$", loc="left", fontsize=13)
    ax.text(0.3, 9.0, "1枚の刃が削り取る切りくずの断面", fontsize=11)
    ax.add_patch(Rectangle((2, 6.6), 5, 1.2, fc=ORANGE, alpha=0.8, ec="#b9770e"))
    ax.text(4.5, 7.2, "面積 $A = h \\times b$", ha="center", va="center", fontsize=11)
    arrow(ax, (2, 6.2), (7, 6.2)); ax.text(3.2, 5.55, "幅 $b$（$\\approx a_p$）", color=BLUE, fontsize=11)
    arrow(ax, (1.6, 6.6), (1.6, 7.8)); ax.text(0.1, 7.05, "厚さ $h$", color=BLUE, fontsize=11)
    arrow(ax, (8.3, 9.2), (8.3, 7.0), color=RED, style="-|>", lw=3, ms=18)
    ax.text(8.5, 8.0, "切削力\n$F_c$", color=RED, fontsize=12)
    ax.text(0.3, 5.0, "$F_c = k_c \\times A$", fontsize=16, va="top")
    ax.text(0.3, 3.9, "$k_c$ ＝ 断面 1 mm² を削るのに必要な力 [N/mm²]", fontsize=10.5, va="top")
    ax.text(0.3, 3.0, "炭素鋼 $k_c \\approx 2000$ N/mm²\n→ 1 mm² で約 2000 N（約 200 kg の重さ）", fontsize=10.5, va="top")
    ax.text(0.3, 1.4, "材料が硬いほど大きい（目安）\nアルミ ≈ 700　鋳鉄 ≈ 1100〜1500\nステンレス ≈ 2000〜2500", fontsize=10.5, va="top")
    ax.text(0.3, -1.0, "切削動力 $P = F_c \\times$（刃先の速さ）", fontsize=11, va="top")

    fig.savefig(os.path.join(OUT, "cutting_params.png"))
    plt.close(fig)


# ---------------------------------------------------------------- 理論粗さ
def turning_roughness():
    fig = plt.figure(figsize=(15, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35], height_ratios=[1, 1])

    # 左: 原理図（縮尺は誇張）
    ax = fig.add_subplot(gs[:, 0]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("刃先の丸み（ノーズ半径 $r$）が送り $f$ ごとに残す山", loc="left", fontsize=13)
    r, f = 3.0, 4.4
    h = r - np.sqrt(r**2 - (f / 2) ** 2)
    xs = np.linspace(-1.6 * f, 1.6 * f, 800)
    k = np.round(xs / f)
    prof = r - np.sqrt(np.maximum(r**2 - (xs - k * f) ** 2, 0))
    ax.fill_between(xs, -0.6, prof, color="#d5d8dc")
    ax.plot(xs, prof, color="#566573")
    for c in (-f, 0, f):
        ax.add_patch(Circle((c, r), r, fill=False, ec=GRAY, ls="--", lw=0.8))
    ax.add_patch(Circle((0, r), r, fill=False, ec=BLUE, lw=1.8))
    ax.plot([0, r * np.cos(np.radians(-60))], [r, r + r * np.sin(np.radians(-60))], color=BLUE)
    ax.text(0.5, 2.2, "$r$", color=BLUE, fontsize=15)
    arrow(ax, (0, -1.0), (f, -1.0), color=RED); ax.text(f / 2, -1.7, "$f$（1回転で進む量）", color=RED, ha="center", fontsize=12)
    ax.plot([f / 2, f / 2 + 2.2], [h, h], color=GREEN, ls=":"); ax.plot([f, f + 1.6], [0, 0], color=GREEN, ls=":")
    arrow(ax, (f + 1.4, 0), (f + 1.4, h), color=GREEN)
    ax.text(f + 1.6, h / 2 - 0.1, "$R_z$（山の高さ）", color=GREEN, fontsize=12)
    ax.text(-7.0, 9.2, "工作物（下側）を工具の丸い刃先がなぞる。\n刃先は 1回転ごとに $f$ ずつ横に進むので、\n円弧と円弧の間に山が残る。", fontsize=11, va="top")
    ax.text(-7.0, -3.4, "$R_{z,\\mathrm{th}} = r - \\sqrt{r^2 - (f/2)^2} \\;\\approx\\; \\dfrac{f^2}{8r}$", fontsize=15)
    ax.text(-7.0, -4.5, "（縦方向は誇張して描いている）", fontsize=10, color=GRAY)
    ax.set_xlim(-7.2, 9.5); ax.set_ylim(-4.9, 9.4)

    # 右: 実際の形状（同じ縮尺で 2 条件）
    for i, (f, r, lab) in enumerate([(0.1, 0.4, "仕上げ"), (0.3, 0.8, "粗加工")]):
        ax = fig.add_subplot(gs[i, 1])
        x = np.linspace(0, 1.2, 4000)
        k = np.round(x / f)
        z = (r - np.sqrt(r**2 - (x - k * f) ** 2)) * 1000  # μm
        rz = f**2 / (8 * r) * 1000
        y = z - rz   # 山の頂点を 0 にそろえる（刃先の中心が通った所が谷）
        ax.fill_between(x, -25, y, color="#d5d8dc"); ax.plot(x, y, color="#566573")
        ax.axhline(-rz, color=GREEN, ls=":", lw=1)
        ax.set_ylim(-25, 3); ax.set_xlim(0, 1.2)
        ax.set_ylabel("高さ [μm]"); ax.grid(alpha=0.3)
        ax.set_title(f"{lab}：$f = {f}$ mm/rev,  $r = {r}$ mm  →  $R_{{z,\\mathrm{{th}}}} \\approx {rz:.1f}\\ \\mu\\mathrm{{m}}$",
                     loc="left", fontsize=12)
        if i == 1:
            ax.set_xlabel("送り方向の位置 [mm]（横 1 mm に対して縦は μm。縦を約 50 倍に拡大）")
        else:
            ax.tick_params(labelbottom=False)
    fig.savefig(os.path.join(OUT, "turning_roughness.png"))
    plt.close(fig)


if __name__ == "__main__":
    cutting_params()
    turning_roughness()
    print("saved to", os.path.abspath(OUT))
