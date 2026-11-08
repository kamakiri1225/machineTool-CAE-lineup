"""5.2 補足：Z 方向（3×3）が必要になる工具の写真（Wikimedia Commons、ライセンスは図中と本文に記載）"""
import sys
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
import matplotlib.image as mpimg
P = [("helical_endmill_Fresa.jpg", "ねじれ刃エンドミル", "刃がらせん状なので、切削力に\n軸方向の成分 $F_z$ が出る",
      "Ricce / Public domain"),
     ("ball_endmill_Fraises_spheriques.jpg", "ボールエンドミル", "先端が球面なので、Z に揺れると\n削り厚さ $h$ が変わる",
      "Bricetofly / CC0"),
     ("face_mill_VHM-WSP-ALU-F-45-IH2.jpg", "正面フライス", "底面側の刃でも削るので、Z に揺れると\n削り厚さ $h$ が変わる",
      "Raizy / CC BY-SA 3.0")]
fig, axs = plt.subplots(1, 3, figsize=(13, 5.6))
for ax, (f, name, why, cr) in zip(axs, P):
    ax.imshow(mpimg.imread(D / "photos" / f)); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(name, fontsize=17, pad=8)
    ax.text(0.5, -0.05, why, transform=ax.transAxes, ha="center", va="top", fontsize=14)
    ax.text(0.99, 0.01, "出典: Wikimedia Commons\n" + cr, transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            color="white", bbox=dict(fc="black", alpha=0.45, lw=0))
fig.suptitle("Z 方向まで入れた 3×3 の FRF が必要になる工具", fontsize=18)
fig.tight_layout(); fig.savefig(R / "tool_types_3d.png", dpi=110, bbox_inches="tight"); print("ok")
