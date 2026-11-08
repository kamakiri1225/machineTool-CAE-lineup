"""ref005 の図: 工作機械のどこに、何を目的にトポロジー最適化を使うか（剛性ループ）"""
import os
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref005")
fig, ax = plt.subplots(figsize=(15, 8)); ax.set_xlim(0, 15); ax.set_ylim(0, 8); ax.set_aspect("equal"); ax.axis("off")


def part(x, y, w, h, c, name, tx=0.5, ty=0.5):
    ax.add_patch(Rectangle((x, y), w, h, fc=c, ec="k", lw=1.5, alpha=0.55)); ax.text(x + w * tx, y + h * ty, name, ha="center", va="center", fontsize=12)


FIX, MOV = "#85c1e9", "#f5b041"
part(0.5, 0.5, 7.0, 1.2, FIX, "ベッド（固定）", tx=0.2, ty=0.3)
part(5.4, 1.7, 1.6, 5.3, FIX, "コラム\n（固定）", tx=0.3)
part(3.0, 5.2, 2.4, 1.3, MOV, "主軸頭\n（可動）", tx=0.3, ty=0.35)
ax.add_patch(Rectangle((3.9, 3.6), 0.35, 1.6, fc=GRAY, ec="k")); ax.text(3.75, 4.3, "工具", ha="right", fontsize=10)
part(1.2, 2.0, 4.0, 0.8, MOV, "テーブル（可動）", tx=0.25)
part(2.6, 2.8, 1.8, 0.8, "#d5d8dc", "ワーク", tx=0.75)
for x in (1.6, 4.4): ax.add_patch(Rectangle((x, 1.7), 0.4, 0.3, fc=BLUE))
ax.add_patch(Rectangle((5.25, 5.5), 0.15, 0.8, fc=BLUE))
# 力のループ
loop = [(4.1, 3.6), (4.1, 6.1), (6.6, 6.1), (6.6, 1.3), (3.3, 1.3), (3.3, 2.4), (3.6, 3.6)]
for p0, p1 in zip(loop[:-1], loop[1:]):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", color=RED, lw=2.5, mutation_scale=16, alpha=0.85))
ax.text(4.3, 7.4, "切削力が通る「剛性ループ」（工具 → 主軸頭 → コラム → ベッド → テーブル → ワーク）\nこのループ全体の柔らかさが工具とワークの間の揺れになる",
        fontsize=11, color=RED, ha="center")
ax.text(0.3, 4.5, "青い小さな四角 = ガイド\n（ばね。ref004）", fontsize=10, color=BLUE)
# 右側の説明
rows = [(FIX, "固定部（ベッド・コラム）", "目的: 静剛性・動剛性を最大\n制約: 重さ・鋳物で作れる形\n→ リブの配置・肉厚に読み替える"),
        (MOV, "可動部（主軸頭・テーブル・サドル）", "目的: 軽くして固有振動数を上げる\n→ 加減速・サーボゲインが上がる（ref003）\n制約: 剛性・取付面"),
        ("#ffffff", "共通の注意", "・荷重は 1 つではない（加工方向・位置が変わる）\n・ガイド・ボールねじのばねを入れて解く\n・熱変形（左右対称）・中子・砂抜き穴")]
for i, (c, t, s) in enumerate(rows):
    y = 6.1 - i * 2.35
    ax.add_patch(FancyBboxPatch((8.0, y - 1.0), 6.7, 1.95, boxstyle="round,pad=0.05", fc=c, alpha=0.35, ec="k"))
    ax.text(8.2, y + 0.65, t, fontsize=12, weight="bold"); ax.text(8.2, y + 0.3, s, fontsize=11, va="top")
fig.savefig(os.path.join(OUT, "topopt_where_in_machine.png")); plt.close(fig); print("ok")
