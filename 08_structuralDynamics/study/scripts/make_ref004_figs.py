"""ref004 の図: ボールねじ送り系の軸方向剛性とナット位置"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref004")
os.makedirs(OUT, exist_ok=True)

# 仮定（φ40 のボールねじ。実際はカタログ値に置き換える）
E = 206e3            # MPa
d_root = 34.0        # ねじ谷径 [mm]
A = np.pi / 4 * d_root**2
L = 1000.0           # 支持軸受間の距離 [mm]
K_NUT = 1000e3       # ナット剛性 [N/mm]（= 1000 N/µm、予圧ダブルナットの例）
K_BRG = 1000e3       # サポート軸受 1 組の剛性 [N/mm]


def series(*ks):
    return 1 / sum(1 / k for k in ks)


def k_fixed_support(a):
    return series(A * E / a, K_BRG, K_NUT)          # 固定側から a の位置にナット


def k_fixed_fixed(a):
    left = series(A * E / a, K_BRG); right = series(A * E / (L - a), K_BRG)
    return series(left + right, K_NUT)               # 両側が並列、ナットは直列


a = np.linspace(50, 950, 400)
fig, axs = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
ax = axs[0]
ax.plot(a, k_fixed_support(a) / 1e3, color=RED, lw=2.5, label="固定 − 支持（片側だけ軸方向に固定）")
ax.plot(a, k_fixed_fixed(a) / 1e3, color=BLUE, lw=2.5, label="固定 − 固定（両側を軸方向に固定）")
ax.set_xlabel("ナットの位置 $a$ [mm]（モータ側の固定軸受から）"); ax.set_ylabel("送り系の軸方向剛性 [N/µm]")
ax.set_ylim(0, 700); ax.grid(alpha=0.3); ax.legend(loc="upper right")
ax.set_title("ナット（＝テーブル）の位置で剛性が変わる（φ40、軸受間 1000 mm の例）", loc="left", fontsize=12)
ax.annotate("モータから遠いほど\nねじ軸が長く伸びる → 柔らかい", xy=(900, k_fixed_support(900) / 1e3), xytext=(560, 80),
            color=RED, fontsize=11, arrowprops=dict(arrowstyle="->", color=RED))
ax.annotate("中央がいちばん柔らかい", xy=(500, k_fixed_fixed(500) / 1e3), xytext=(380, 560),
            color=BLUE, fontsize=11, arrowprops=dict(arrowstyle="->", color=BLUE))

ax = axs[1]
pos = [100, 500, 900]
parts = {"ねじ軸": lambda x: A * E / x, "サポート軸受": lambda x: K_BRG, "ナット": lambda x: K_NUT}
bottom = np.zeros(3)
cols = [RED, ORANGE, BLUE]
for (nm, fn), c in zip(parts.items(), cols):
    v = np.array([1 / fn(p) for p in pos]) * 1e3 * 1000    # µm/kN
    ax.bar([f"a = {p} mm" for p in pos], v, bottom=bottom, color=c, label=nm)
    bottom += v
ax.set_ylabel("1 kN あたりの縮み [µm]（＝ 1/剛性）")
ax.set_title("固定−支持の場合の内訳（縮みは直列に足し算）", loc="left", fontsize=12)
ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ballscrew_stiffness.png")); plt.close(fig)
print({p: round(k_fixed_support(p) / 1e3) for p in pos}, {p: round(k_fixed_fixed(p) / 1e3) for p in pos})
