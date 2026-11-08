"""paper001 の図: Law・Altintas・Phani (2013) の手法の流れ（自作の模式図）と、論文の数値のグラフ化"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import FancyBboxPatch, Rectangle
OUT = os.path.join(os.path.dirname(__file__), "..", "papers", "images", "paper001"); os.makedirs(OUT, exist_ok=True)


def box(ax, x, y, w, h, t, c, fs=11):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.03", fc=c, alpha=0.16, ec=c, lw=2))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs)


def arr(ax, a, b, c="k", rad=0.0):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8, connectionstyle=f"arc3,rad={rad}"))


# ---------------------------------------------------------------- 1. 手法の流れ
fig, ax = plt.subplots(figsize=(16, 8.4)); ax.set_xlim(0, 16); ax.set_ylim(0, 8.4); ax.axis("off")
ax.text(0.2, 8.0, "Law・Altintas・Phani (2013) の手法：「速い全機モデル」で、位置ごとのびびり安定限界を見て、弱い部品だけを直す", fontsize=13, weight="bold")
steps = [("① 部品ごとに縮約", "主軸ハウジング・コラム・ベース\nクロススライド・テーブルを別々に\n「接触面の自由度＋重要なモード」\nだけ残す（CMS の改良版）", BLUE),
         ("② 代表点で組み立て", "接触面の中心に代表点を置き\n面の節点の重み付き平均（MPC）\nで結ぶ。メッシュが合わなくてよい\nガイド・ねじはばね", GREEN),
         ("③ 位置を変える", "工具の位置（Z ストローク）を変える\nと接触する節点が変わる\n→ 拘束式を付け替えるだけ\n1 位置 約 1 分（全体 FEM は 12 時間）", ORANGE),
         ("④ 工具を結合", "工具・ホルダ（はり要素）の FRF を\n主軸端の FRF に受容性結合\n→ 工具先端の FRF", GRAY)]
for i, (t, d, c) in enumerate(steps):
    x = 0.2 + i * 3.95
    box(ax, x, 5.0, 3.6, 2.5, "", c); ax.text(x + 0.15, 7.15, t, fontsize=12.5, weight="bold", color=c)
    ax.text(x + 0.15, 6.75, d, fontsize=10.5, va="top")
    if i < 3: arr(ax, (x + 3.62, 6.25), (x + 3.93, 6.25))
steps2 = [("⑤ 安定限界", "送り方向 0〜360° ごとに\n速度によらない絶対安定限界\n→ 目標 4 mm に届かない\n位置・方向と原因のモード", RED),
          ("⑥ 補強してみる", "原因 = コラムの曲げ・ねじれ\nフランジ・内部リブを追加\n→ 効くが質量 +15 %", GRAY),
          ("⑦ トポロジー最適化", "コラムを設計領域に\n工具先端に 3000 N（X・Y・Z）\n3 つの位置で、体積最小\n（コンプライアンス制約）", BLUE),
          ("⑧ 再び ①〜⑤", "元と同じ質量で、補強と同じ動特性\n動剛性 +20〜30 %\n安定限界 +20〜45 %", GREEN)]
for i, (t, d, c) in enumerate(steps2):
    x = 12.05 - i * 3.95
    box(ax, x, 1.2, 3.6, 2.6, "", c); ax.text(x + 0.15, 3.45, t, fontsize=12.5, weight="bold", color=c)
    ax.text(x + 0.15, 3.05, d, fontsize=10.5, va="top")
    if i < 3: arr(ax, (x - 0.02, 2.5), (x - 0.33, 2.5))
arr(ax, (14.0, 4.98), (14.0, 3.82))
ax.text(0.2, 0.55, "ポイント：トポロジー最適化そのものは「静的（複数位置の荷重）」。動特性・びびりは、速い全機モデルで「評価」に使い、どこを直すかを決める", fontsize=11.5, color=BLUE)
fig.savefig(os.path.join(OUT, "law2013_flow.png")); plt.close(fig)

# ---------------------------------------------------------------- 2. 論文の数値
fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), gridspec_kw=dict(width_ratios=[1.2, 1, 1]))
ax = axs[0]
parts = ["主軸ハウジング", "コラム", "ベース", "クロススライド", "テーブル"]
full = [21992, 45983, 97214, 50201, 15623]; red = [1037 + 27, 1533 + 28, 1945 + 46, 2598 + 27, 1908 + 24]
y = np.arange(5)
ax.barh(y - 0.2, full, 0.4, color=GRAY, label="全体 FEM の自由度"); ax.barh(y + 0.2, red, 0.4, color=BLUE, label="縮約後（接触面＋重要モード）")
ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(parts); ax.invert_yaxis(); ax.legend(fontsize=9.5, loc="lower right")
ax.set_title("① 縮約で 231,013 → 9,173 自由度（約 1/25）\n1 位置 12 時間 → 約 1 分", loc="left", fontsize=11.5); ax.grid(alpha=0.3, axis="x")
ax = axs[1]
pos = ["上", "中", "下"]; before = [1.8, 3.0, 5.6]; after = [2.5, 3.6, 8.0]
x = np.arange(3)
ax.bar(x - 0.2, before, 0.4, color=GRAY, label="元の機械"); ax.bar(x + 0.2, after, 0.4, color=GREEN, label="トポロジー最適化後のコラム")
ax.axhline(4, color=RED, ls="--"); ax.text(-0.45, 4.15, "目標 4 mm", color=RED, ha="left")
for i in range(3):
    ax.text(i + 0.2, after[i] + 0.15, f"+{(after[i]/before[i]-1)*100:.0f} %", ha="center", color=GREEN, fontsize=11)
ax.set_xticks(x); ax.set_xticklabels([f"主軸頭 {p}" for p in pos]); ax.set_ylabel("絶対安定限界 [mm]（全送り方向の最小）")
ax.set_title("② 位置ごとの安定限界", loc="left", fontsize=11.5); ax.legend(fontsize=9.5, loc="upper left", bbox_to_anchor=(0, 0.92)); ax.grid(alpha=0.3, axis="y")
ax = axs[2]
lab = ["コラム質量", "静剛性", "1 次曲げの\n動剛性", "約 90 Hz の\n動剛性"]
v_stiff = [115, None, None, None]; v_opt = [100, 140, 120, 130]
ax.bar(np.arange(4), v_opt, color=[GREEN, BLUE, BLUE, BLUE])
ax.bar([0], [115], color="none", edgecolor=RED, lw=2, ls="--"); ax.text(0, 118, "補強案 115", ha="center", color=RED, fontsize=10)
for i, v in enumerate(v_opt):
    ax.text(i, v / 2, f"{v}", ha="center", color="white", fontsize=12, weight="bold")
ax.axhline(100, color=GRAY, lw=1); ax.set_xticks(range(4)); ax.set_xticklabels(lab, fontsize=10); ax.set_ylim(0, 160)
ax.set_ylabel("元の機械 = 100"); ax.set_title("③ トポロジー最適化後：同じ質量で硬く", loc="left", fontsize=11.5); ax.grid(alpha=0.3, axis="y")
fig.suptitle("Law・Altintas・Phani (2013) の数値（論文本文から）", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "law2013_numbers.png")); plt.close(fig); print("ok")
