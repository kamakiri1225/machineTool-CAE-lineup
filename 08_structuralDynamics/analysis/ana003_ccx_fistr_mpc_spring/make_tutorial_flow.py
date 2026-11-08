"""ブログ用：作業の流れと、どのファイルのどのキーワードに書くかの一覧図"""
import os
from figstyle_local import plt
from matplotlib.patches import FancyBboxPatch
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "tutorial")
steps = [("1", "部品を作る", "節点・要素\n（節点は共有しない）", "*NODE / *ELEMENT", "!NODE / !ELEMENT", ".msh"),
         ("2", "面の節点を集める", "節点集合", "*NSET", "!NGROUP", ".msh"),
         ("3", "代表点を足す", "REF・ROT 節点", "*NODE", "!NODE", ".msh"),
         ("4", "面と代表点を結ぶ", "拘束式 MPC", "*EQUATION", "!EQUATION", ".msh"),
         ("5", "ばねでつなぐ", "ばね要素＋ばね定数", "SPRING2 + *SPRING", "511 (.msh) + !SPRING_D (.cnt)", ".msh/.cnt"),
         ("6", "固定・荷重", "境界条件", "*BOUNDARY / *CLOAD", "!BOUNDARY / !CLOAD", ".cnt"),
         ("7", "解いて確かめる", "F/k と比べる", "ccx -i model", "fistr1", "")]
fig, ax = plt.subplots(figsize=(16, 7.4)); ax.set_xlim(0, 16); ax.set_ylim(0, 7.4); ax.axis("off")
cols = ["#1f618d", "#e74c3c", "#d4ac0d", "#8e44ad", "#c0392b", "#16a085", "#2c3e50"]
ax.text(0.2, 6.95, "手順", fontsize=13, weight="bold"); ax.text(3.0, 6.95, "作るもの", fontsize=13)
ax.text(6.3, 6.95, "CalculiX（model.inp 1 つ）", fontsize=13, color="#1f618d")
ax.text(10.6, 6.95, "FrontISTR（model.msh と model.cnt）", fontsize=13, color="#c0392b")
for i, (n, t, what, c, f, fl) in enumerate(steps):
    y = 6.2 - i * 0.92
    ax.add_patch(FancyBboxPatch((0.2, y - 0.32), 2.6, 0.64, boxstyle="round,pad=0.04", fc=cols[i], alpha=0.85, ec="none"))
    ax.text(0.35, y, f"{n}. {t}", color="white", fontsize=12.5, va="center", weight="bold")
    ax.text(3.0, y, what, fontsize=11.5, va="center")
    ax.text(6.3, y, c, fontsize=12, va="center", family="monospace", color="#1f618d")
    ax.text(10.6, y, f, fontsize=12, va="center", family="monospace", color="#c0392b")
    if i < len(steps) - 1:
        ax.annotate("", xy=(1.5, y - 0.6), xytext=(1.5, y - 0.33), arrowprops=dict(arrowstyle="-|>", color="#888"))
ax.text(0.2, 0.05, "※ CalculiX はすべて 1 つの .inp に書く。FrontISTR は「形（.msh）」と「解析条件・ばね定数（.cnt）」に分かれる", fontsize=11.5, color="#555")
fig.savefig(os.path.join(OUT, "flow.png")); print("ok")
