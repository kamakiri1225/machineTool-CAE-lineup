"""ref008 7 章の図：社内でトポロジー最適化を使う部品の見極め（目安）
   (a) 部品 × 判断基準のスコア（3 = 向いている、1 = 向いていない）。点数は 7 章の文献をもとにした目安で、社内データで付け直す前提
   (b) 研削盤とマシニングセンタで、トポロジー最適化の前提条件がどう違うか"""
import os
import numpy as np
from figstyle import plt
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); os.makedirs(OUT, exist_ok=True)

CRIT = ["荷重条件を\n決めやすい", "結合部より構造の\n柔らかさが効く", "作り方の制約に\n合わせやすい", "効果が価値に\nつながる", "台数・回収の\n見込み"]
PARTS = [  # 部品, スコア（各基準 1〜3）
    ("ベッド（固定・鋳物）",            [3, 2, 2, 2, 3]),
    ("コラム（固定・鋳物）",            [3, 3, 2, 3, 3]),
    ("サドル・テーブル（可動）",        [2, 1, 2, 3, 3]),
    ("主軸頭・ラム（可動）",            [2, 2, 2, 3, 3]),
    ("旋回テーブル",                    [2, 2, 2, 3, 2]),
    ("ATC アーム・マガジン",            [3, 3, 3, 2, 3]),
    ("ロボットハンド・治具（積層造形）", [3, 3, 3, 2, 2]),
    ("研削盤の砥石頭・コラム",          [3, 3, 2, 3, 2]),
]
S = np.array([p[1] for p in PARTS], float)

fig, ax = plt.subplots(1, 2, figsize=(16, 6.4), gridspec_kw={"width_ratios": [1.45, 1]})
a = ax[0]
im = a.imshow(S, cmap="RdYlGn", vmin=0.6, vmax=3.4, aspect="auto")
for i in range(S.shape[0]):
    for j in range(S.shape[1]):
        a.text(j, i, ["", "不利", "条件付き", "向く"][int(S[i, j])], ha="center", va="center", fontsize=10.5)
a.set_xticks(range(len(CRIT))); a.set_xticklabels(CRIT, fontsize=10); a.xaxis.tick_top()
a.set_yticks(range(len(PARTS))); a.set_yticklabels([f"{p[0]}（計 {int(sum(p[1]))}）" for p in PARTS], fontsize=10.5)
a.set_title("(a) 部品ごとの向き・不向き（目安。社内データで付け直す）", loc="left", fontsize=11.5, pad=42)

b = ax[1]; b.axis("off")
rows = [("", "研削盤", "マシニングセンタ"),
        ("加工の力", "小さく安定\n（実測しやすい）", "大きく変動\n（工具・条件で変わる）"),
        ("軸の位置の変化", "少ない・遅い", "大きい・速い\n（位置ごとに特性が変わる）"),
        ("重視する性能", "静剛性・熱変位・\n減衰（面粗さ）", "動剛性（びびり）・\n移動質量（加減速）"),
        ("台数・単価", "少量・高単価", "多量・価格競争"),
        ("公開事例", "ナガセ（全体構造）、\nATC アーム、\n砥石マガジン", "Heller・Grob・DMG森精機\n（試作）"),
       ]
y = 0.95
for k, (c0, c1, c2) in enumerate(rows):
    fw = "bold" if k == 0 else "normal"
    b.text(0.0, y, c0, fontsize=10.5, va="top", fontweight="bold", transform=b.transAxes)
    b.text(0.30, y, c1, fontsize=10.5, va="top", fontweight=fw, color="#1f6f3f" if k else "k", transform=b.transAxes)
    b.text(0.66, y, c2, fontsize=10.5, va="top", fontweight=fw, color="#8a3b12" if k else "k", transform=b.transAxes)
    y -= 0.155 if k else 0.09
    b.plot([0, 1], [y + 0.02, y + 0.02], color="#bbb", lw=0.6, transform=b.transAxes)
b.set_title("(b) 研削盤とマシニングセンタの前提の違い", loc="left", fontsize=11.5, pad=42)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "screening.png"), dpi=90); print("ok")
