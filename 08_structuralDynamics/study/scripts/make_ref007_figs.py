"""ref007 の図: 門型平面研削盤とマシニングセンタの違い"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, Polygon

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref007")
os.makedirs(OUT, exist_ok=True)


def arr(ax, p0, p1, c="k", lw=2.5, st="-|>"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=st, color=c, lw=lw, mutation_scale=18))


def R(ax, x, y, w, h, fc, t=None, fs=11, **kw):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec="k", lw=1.4, **kw))
    if t: ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs)


# ================================================================ 1. 構造の違い
fig, axs = plt.subplots(1, 2, figsize=(15, 6.6))
ax = axs[0]; ax.set_xlim(0, 8); ax.set_ylim(0, 8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("立形マシニングセンタ（C 形・片持ち）", loc="left", fontsize=13)
R(ax, 0.5, 0.5, 7, 1.0, "#d5d8dc", "ベッド")
R(ax, 5.3, 1.5, 1.5, 5.8, "#aed6f1", "コラム")
R(ax, 2.6, 5.2, 2.7, 1.3, "#f9e79f", "主軸頭")
R(ax, 3.7, 3.7, 0.4, 1.5, GRAY); ax.text(3.55, 4.3, "工具", ha="right", fontsize=10)
R(ax, 1.2, 1.9, 3.9, 0.6, "#f5cba7", "テーブル（X・Y）", fs=10)
R(ax, 1.3, 1.5, 3.7, 0.4, "#e8daef", "サドル", fs=9)
R(ax, 2.6, 2.5, 1.8, 0.8, "#eaeded", "ワーク", fs=10)
arr(ax, (7.2, 5.0), (7.2, 6.5), RED); ax.text(7.3, 5.7, "Z", color=RED, fontsize=12)
ax.annotate("", xy=(4.6, 6.9), xytext=(5.4, 6.9), arrowprops=dict(arrowstyle="->", color=RED))
ax.text(0.5, 7.4, "主軸頭がコラムから張り出す（片持ち）\n→ 張り出しが曲がりやすい・熱で前に傾く", fontsize=10.5, color=RED)
ax.text(0.5, 0.05, "回転数: 8,000〜20,000 min$^{-1}$、工具径 φ6〜100", fontsize=11)

ax = axs[1]; ax.set_xlim(0, 9); ax.set_ylim(0, 8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("門型平面研削盤（両持ち・正面から見る）", loc="left", fontsize=13)
R(ax, 0.3, 0.5, 8.4, 1.0, "#d5d8dc", "ベッド")
R(ax, 0.7, 1.5, 1.0, 5.2, "#aed6f1", "コラム", fs=10)
R(ax, 7.3, 1.5, 1.0, 5.2, "#aed6f1", "コラム", fs=10)
R(ax, 0.7, 5.7, 7.6, 1.0, "#aed6f1", "クロスレール（梁）")
R(ax, 3.6, 4.2, 1.8, 1.5, "#f9e79f", "砥石頭", fs=10)
ax.add_patch(Circle((4.5, 3.35), 0.85, fc="#bfc9ca", ec="k", lw=1.5)); ax.text(4.5, 3.35, "砥石\nφ300〜500", ha="center", va="center", fontsize=9)
R(ax, 1.9, 1.5, 5.2, 0.55, "#f5cba7", "テーブル（前後に往復）", fs=10)
R(ax, 2.6, 2.05, 3.8, 0.45, "#eaeded", "ワーク（金型・プレート）", fs=9)
arr(ax, (5.6, 4.9), (7.0, 4.9), RED); ax.text(6.1, 5.05, "左右", color=RED, fontsize=10)
ax.add_patch(Circle((1.3, 1.15), 0.0))
ax.text(0.3, 7.3, "梁を両側のコラムで支える（両持ち）、左右対称\n→ 曲がりにくい・熱で傾きにくい。大きなワークを載せられる", fontsize=10.5, color=BLUE)
ax.text(0.3, 0.05, "砥石回転数: 1,000〜3,600 min$^{-1}$、周速 30〜45 m/s", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "structure_compare.png")); plt.close(fig)

# ================================================================ 2. 削り方・力の向きの違い
fig, axs = plt.subplots(1, 2, figsize=(15, 6))
ax = axs[0]; ax.set_xlim(-3.2, 4.2); ax.set_ylim(-2.4, 3.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("フライス（エンドミル・正面フライス）：数枚の刃で厚く削る", loc="left", fontsize=12)
R(ax, -3, -2.2, 7, 1.6, "#eaeded"); ax.text(2.6, -1.8, "ワーク", fontsize=11)
ax.add_patch(Polygon([[-3, -0.6], [0.1, -0.6], [0.1, 0.2], [4, 0.2], [4, -0.6]], closed=False, fill=False, ls="--", ec=GRAY))
R(ax, -1.0, -0.6, 1.1, 3.6, "#bfc9ca"); ax.text(-0.45, 3.2, "工具 φ20〜80", ha="center", fontsize=10)
for y in np.linspace(-0.3, 2.5, 6): ax.plot([-1.0, 0.1], [y, y + 0.35], color="k", lw=1)
ax.annotate("", xy=(0.1, 0.2), xytext=(0.1, -0.6), arrowprops=dict(arrowstyle="<->", color=BLUE))
ax.text(0.25, -0.35, "切込み $a_p$ = 1〜5 mm", color=BLUE, fontsize=11)
arr(ax, (-0.45, 1.2), (-2.6, 1.2), RED); ax.text(-2.9, 1.4, "主分力 $F_t$\n数百 N〜数 kN", color=RED, fontsize=11)
ax.text(0.5, 3.0, "切削速度 $v_c$ = 100〜300 m/min（2〜5 m/s）", fontsize=11)
arr(ax, (1.0, 2.2), (2.6, 2.2), GREEN); ax.text(1.0, 2.4, "送り $v_f$ 0.5〜5 m/min", color=GREEN, fontsize=11)

ax = axs[1]; ax.set_xlim(-3.2, 4.2); ax.set_ylim(-2.4, 3.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("研削：無数の砥粒で薄く削る（押し付ける力が大きい）", loc="left", fontsize=12)
R(ax, -3, -2.2, 7, 1.6, "#eaeded"); ax.text(2.6, -1.8, "ワーク", fontsize=11)
ax.add_patch(Circle((0.3, 1.25), 1.85, fc="#d5dbdb", ec="k", lw=1.5))
rng = np.random.default_rng(1)
for t in rng.uniform(0, 2 * np.pi, 140):
    r = 1.85 - rng.uniform(0, 0.15); ax.plot(0.3 + r * np.cos(t), 1.25 + r * np.sin(t), "k.", ms=2)
ax.text(0.3, 1.5, "砥石\n（砥粒が無数の刃）", ha="center", fontsize=10)
ax.add_patch(Rectangle((-3, -0.6), 3.3, 0.0001))
ax.plot([-3, 0.3], [-0.55, -0.55], color=BLUE, lw=3)
ax.text(-2.9, -0.15, "切込み $a_e$\n0.001〜0.05 mm", color=BLUE, fontsize=10.5)
arr(ax, (0.3, -0.6), (0.3, 1.1), RED, lw=3.5); ax.text(0.45, 0.0, "法線力 $F_n$", color=RED, fontsize=12)
arr(ax, (0.3, -0.6), (1.5, -0.6), ORANGE, lw=2.5); ax.text(1.6, -0.45, "接線力 $F_t$", color=ORANGE, fontsize=11)
ax.text(1.8, 2.7, "$F_n \\approx (2\\sim3)\\,F_t$", color=RED, fontsize=14)
ax.annotate("", xy=(2.25, 1.85), xytext=(1.9, 0.0), arrowprops=dict(arrowstyle="->", lw=2, connectionstyle="arc3,rad=0.4"))
ax.text(2.2, 0.6, "周速 $v_s$\n30〜45 m/s", fontsize=11)
arr(ax, (-1.2, -1.2), (-2.8, -1.2), GREEN); ax.text(-2.9, -1.05, "", fontsize=1)
ax.text(-1.1, -1.35, "← テーブル $v_w$ 10〜30 m/min", color=GREEN, fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "cutting_vs_grinding.png")); plt.close(fig)

# ================================================================ 3. 数値の比較（対数）
items = [("工具の周速 [m/s]", (2, 5), (30, 45)),
         ("主軸回転数 [min$^{-1}$]", (8000, 20000), (1000, 3600)),
         ("1 回の切込み [mm]", (1, 5), (0.001, 0.05)),
         ("加工時の送り速度 [m/min]", (0.5, 5), (10, 30)),
         ("比エネルギー [J/mm$^3$]（1 mm$^3$ 削るのに要る仕事）", (2, 4), (20, 60)),
         ("面粗さ $Ra$ [µm]", (0.8, 3.2), (0.05, 0.4)),
         ("代表的な力 [N]（片方向）", (300, 3000), (100, 500))]
fig, ax = plt.subplots(figsize=(13, 6.2))
for i, (nm, mc, gr) in enumerate(items):
    ax.plot(mc, [i - 0.15] * 2, color=BLUE, lw=10, solid_capstyle="butt")
    ax.plot(gr, [i + 0.15] * 2, color=RED, lw=10, solid_capstyle="butt")
ax.set_xscale("log"); ax.set_xlim(5e-4, 5e4)
ax.set_yticks(range(len(items))); ax.set_yticklabels([n for n, _, _ in items], fontsize=11); ax.invert_yaxis()
ax.plot([], [], color=BLUE, lw=8, label="マシニングセンタ（鋼のフライス加工）")
ax.plot([], [], color=RED, lw=8, label="平面研削盤")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2, fontsize=11); ax.grid(alpha=0.3, axis="x", which="major")
ax.set_title("数値の目安の比較（対数目盛。条件で大きく変わる目安）", loc="left", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "numbers_compare.png")); plt.close(fig)

# ================================================================ 4. 計算例: 同じ動力でも削れる量・力の向きが違う
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
lab = ["フライス", "研削"]
Q = [1667, 167]; P = [3.3, 5.0]; F = [(1000, 400), (167, 420)]
axs[0].bar(lab, Q, color=[BLUE, RED]); axs[0].set_title("除去量 $Q$ [mm$^3$/s]", fontsize=12)
axs[1].bar(lab, P, color=[BLUE, RED]); axs[1].set_title("必要な動力 $P$ [kW]", fontsize=12)
x = np.arange(2); w = 0.35
axs[2].bar(x - w / 2, [F[0][0], F[1][0]], w, color=ORANGE, label="送り・回転方向の力（$F_c$, $F_t$）")
axs[2].bar(x + w / 2, [F[0][1], F[1][1]], w, color=RED, label="面に押し付ける力（背分力・$F_n$）")
axs[2].set_xticks(x); axs[2].set_xticklabels(lab); axs[2].set_title("力 [N]", fontsize=12); axs[2].legend(fontsize=9.5)
for a, vals in zip(axs[:2], [Q, P]):
    for i, v in enumerate(vals): a.text(i, v, f"{v:g}", ha="center", va="bottom", fontsize=11)
for a in axs: a.grid(alpha=0.3, axis="y")
fig.suptitle("計算例（本文 4 章）：研削は 1/10 の量を削るのに、より大きな動力を使い、押し付ける力が主になる", fontsize=12, y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "example_calc.png")); plt.close(fig)
print("ok")
