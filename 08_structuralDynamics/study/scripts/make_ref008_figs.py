"""ref008 の図: 設計プロセスの変化（従来 vs トポロジー最適化を前に置く）と、多目的の定式化のメニュー"""
import os
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import FancyBboxPatch
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); os.makedirs(OUT, exist_ok=True)


def box(ax, x, y, w, h, t, c, fs=11):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03", fc=c, alpha=0.18, ec=c, lw=2))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs)


def arr(ax, x0, y0, x1, y1, c="k"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8))


fig, ax = plt.subplots(figsize=(16, 7)); ax.set_xlim(0, 16); ax.set_ylim(0, 7); ax.axis("off")
ax.text(0.2, 6.6, "従来：経験で形を決めてから、解析で確かめる（後追い）", fontsize=13, color=GRAY, weight="bold")
st = ["過去機種の形を\n流用・経験で決める", "3D CAD で\n詳細形状", "FEM で\n確認", "問題があれば\nリブを足す", "試作・\n打撃試験"]
for i, t in enumerate(st):
    box(ax, 0.3 + i * 3.1, 4.9, 2.5, 1.2, t, GRAY)
    if i < 4: arr(ax, 2.85 + i * 3.1, 5.5, 3.35 + i * 3.1, 5.5, GRAY)
ax.annotate("", xy=(7.6, 4.85), xytext=(13.7, 4.85), arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.5, connectionstyle="arc3,rad=-0.25"))
ax.text(10.6, 3.85, "手戻り（重くなりがち）", color=RED, ha="center", fontsize=11)
ax.text(0.2, 3.35, "トポロジー最適化を前に置く：要求（荷重・振動・重さ）から形を出し、人が鋳物に読み替える", fontsize=13, color=BLUE, weight="bold")
st = ["要求を数値で決める\n荷重ケース・びびり周波数\n質量・ばね（ガイド）", "設計領域だけ作る\n（外形・取付面・\n中子の制約）", "多目的トポロジー最適化\n静剛性・固有値・\n動コンプライアンス", "鋳物のリブに\n読み替え\n（設計者の判断）", "全機モデルで\n検証（FRF・安定限界）\n→ 試作"]
cols = [ORANGE, GREEN, BLUE, RED, GRAY]
for i, (t, c) in enumerate(zip(st, cols)):
    box(ax, 0.3 + i * 3.1, 1.2, 2.5, 1.6, t, c, fs=10.5)
    if i < 4: arr(ax, 2.85 + i * 3.1, 2.0, 3.35 + i * 3.1, 2.0, BLUE)
ax.text(0.3, 0.45, "価値：① 最初の案から軽く・硬い（試作回数が減る）　② なぜこの形かを説明できる（構造の合理性）　③ 要求仕様が数字になる（設計プロセスの標準化）",
        fontsize=11.5, color=BLUE)
fig.savefig(os.path.join(OUT, "process_change.png")); plt.close(fig); print("ok")
