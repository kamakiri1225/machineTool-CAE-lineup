"""docs の図で共通に使うスタイル（日本語フォント + TeX風の数式）"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

_CANDIDATES = [
    os.path.expanduser("~/.fonts/meiryo.ttc"),
    "/mnt/c/Windows/Fonts/meiryo.ttc",
    "/mnt/c/Windows/Fonts/YuGothM.ttc",
]
for f in _CANDIDATES:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["mathtext.fontset"] = "cm"   # 数式は Computer Modern（TeX の見た目）
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.facecolor"] = "white"

BLUE, RED, ORANGE, GRAY, GREEN = "#1f618d", "#c0392b", "#e67e22", "#7f8c8d", "#1e8449"
