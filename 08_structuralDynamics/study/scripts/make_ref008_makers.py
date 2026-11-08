"""ref008 6 章の図：工作機械メーカーが商品づくりにトポロジー最適化を使っている事例の地図
   横 = トポロジー最適化を使った部位、縦 = 会社・機種、印 = 公開情報での段階（販売機種 / 参考出展 / 周辺機器・工具 / 匿名事例）"""
import os
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref008"); os.makedirs(OUT, exist_ok=True)
COLS = ["ベッド・コラム\n（固定部の鋳物）", "テーブル・サドル・\n旋回台（可動部）", "構造部品全般\n（可動・固定）", "ロボットハンド・\n工具（周辺）"]
LV = {"販売機種（新機種として発表）": (RED, "o"), "参考出展・試作": (ORANGE, "s"), "周辺機器・工具（販売）": (BLUE, "D"), "匿名事例・研究機関": (GRAY, "^")}
ROWS = [  # 会社・機種, 列, 段階, 注記
    ("Heller H 2000（独, 2019）", 0, "販売機種（新機種として発表）", "鋳物のベッドとコラムをトポロジー最適化。軽く・硬く"),
    ("Heller MPC 10000（独）", 2, "販売機種（新機種として発表）", "鋼製の構造部品をトポロジー最適化"),
    ("Grob G150（独, 2021 試験記事）", 1, "販売機種（新機種として発表）", "軸を構成する全ユニットをトポロジー最適化。テーブルは剛性↑"),
    ("ナガセインテグレックス SGX-126/104\n（日, JIMTOF2022）", 0, "販売機種（新機種として発表）", "曲面の本体。静剛性 2 倍・動剛性 1.25 倍"),
    ("ROMI（伯, Siemens 事例）", 2, "販売機種（新機種として発表）", "最新機種の質量 −30 %（前世代比）"),
    ("DMG森精機（日, 2018 / 2019）", 2, "参考出展・試作", "動剛性 250 %・静剛性 120 %・100 kg 以上軽く。鋳物で製作"),
    ("DMG森精機 Robo2Go ハンド", 3, "周辺機器・工具（販売）", "−62 %、部品数 −60 %"),
    ("Sandvik Coromant CoroMill 390 AM", 3, "周辺機器・工具（販売）", "長い突き出しでびびりにくい"),
    ("Fraunhofer IPA（匿名の工作機械メーカー）", 1, "匿名事例・研究機関", "旋回テーブル −63 %、変形 7.5→9.5 µm"),
]
fig, ax = plt.subplots(figsize=(15, 7.6))
for i, (name, c, lv, note) in enumerate(ROWS):
    y = len(ROWS) - 1 - i; col, mk = LV[lv]
    ax.scatter(c, y, s=260, color=col, marker=mk, zorder=3, edgecolor="k")
    ax.text(c + 0.13, y, note, va="center", fontsize=10)
ax.set_yticks(range(len(ROWS))); ax.set_yticklabels([r[0] for r in ROWS][::-1], fontsize=10.5)
ax.set_xticks(range(len(COLS))); ax.set_xticklabels(COLS, fontsize=11); ax.set_xlim(-0.4, 4.4); ax.set_ylim(-0.7, len(ROWS) - 0.3)
ax.xaxis.tick_top(); ax.grid(alpha=0.3)
for lv, (col, mk) in LV.items():
    ax.scatter([], [], s=140, color=col, marker=mk, edgecolor="k", label=lv)
ax.legend(loc="upper right", fontsize=10, title="公開情報で確認できた段階")
ax.set_title("工作機械メーカーが商品づくりにトポロジー最適化を使っている事例（2026-11 時点、公開情報）", loc="left", fontsize=13, pad=48)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "maker_cases.png"), dpi=90); print("ok")
