"""2 章補足：このモデルが実機のどのクラスに当たるか（写真：Wikimedia Commons CC0）"""
import sys
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
import matplotlib.image as mpimg
fig, (a0, a1) = plt.subplots(1, 2, figsize=(16, 5.6), gridspec_kw=dict(width_ratios=[0.8, 1.2]))
a0.imshow(mpimg.imread(D / "photos" / "mazak_econo_center_museum.jpg")); a0.axis("off")
a0.set_title("同じ形式（C 形コラム・テーブル移動）の立形 MC\nヤマザキマザック工作機械博物館の展示機", fontsize=13, loc="left")
a0.text(0.99, 0.01, "出典: Wikimedia Commons, 先従隗始 / CC0", transform=a0.transAxes, ha="right", va="bottom", fontsize=9, color="white", bbox=dict(fc="black", alpha=0.5, lw=0))
rows = [["項目", "このモデル", "小形立形 MC（例：VCN-430 クラス）"],
        ["形式", "C 形コラム＋移動テーブル", "C 形コラム＋移動テーブル"],
        ["ベッドの大きさ", "1,000 × 900 mm", "—"],
        ["テーブル", "500 × 400 mm", "900 × 430 mm"],
        ["コラム高さ", "ベッド上 800 mm", "Z ストローク 510 mm"],
        ["主軸", "回転数 500〜16,000 min⁻¹ で計算", "12,000 min⁻¹（#40）"],
        ["工具・ワーク", "テストバー 60 mm 角\nワーク 100×100×50 mm", "エンドミル・正面フライス\n小〜中物の部品"]]
a1.axis("off"); tb = a1.table(cellText=rows, loc="center", cellLoc="left", colWidths=[0.2, 0.38, 0.42]); tb.auto_set_font_size(False); tb.set_fontsize(12.5); tb.scale(1, 2.6)
for j in range(3): tb[0, j].set_facecolor("#d6eaf8")
a1.set_title("大きさ・回転数の比較（実機の値は代表的な仕様。正確には各社の仕様書で確認）", fontsize=13, loc="left")
fig.tight_layout(); fig.savefig(R / "machine_class.png", dpi=95, bbox_inches="tight"); print("ok")
