"""ana003 後処理: 分布ばね（ana001）と RBE2 / RBE3 の固有振動数比較"""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from plot_modes import plt

D = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(D, "results"), exist_ok=True)


def f_ccx(path):
    blk = open(path).read().split("E I G E N V A L U E   O U T P U T")[1].split("P A R T I C I P A T I O N")[0]
    return np.array([float(l.split()[3]) for l in blk.splitlines() if len(l.split()) == 5 and l.split()[0].isdigit()])


def f_fis(path):
    out = []
    for ln in open(path):
        p = ln.split()
        if len(p) >= 4 and p[0].isdigit() and int(p[0]) == len(out) + 1 and "E" in p[3]:
            out.append(float(p[3]))
        if len(out) == 15:
            break
    return np.array(out)


A1 = os.path.join(D, "..", "ana001_ccx_fistr_modal")
cases = {
    "分布ばね（ana001）": (f_ccx(os.path.join(A1, "ccx", "machine.dat")), f_fis(os.path.join(A1, "fistr", "0.log"))),
    "代表点＋RBE2（面を剛体）": (f_ccx(os.path.join(D, "rbe2", "ccx", "machine.dat")), f_fis(os.path.join(D, "rbe2", "fistr", "0.log"))),
    "代表点＋RBE3（面の平均）": (f_ccx(os.path.join(D, "rbe3", "ccx", "machine.dat")), f_fis(os.path.join(D, "rbe3", "fistr", "0.log"))),
}
NAMES = ["コラム Y曲げ", "テーブル X（ボールねじ）", "コラム X曲げ", "テーブル Y（横）", "テーブル Zねじれ", "テーブル Z＋曲げ"]
rows = ["| モード | " + " | ".join(f"{k} CalculiX / FrontISTR [Hz]" for k in cases) + " |",
        "|---|" + "---|" * len(cases)]
for i in range(6):
    rows.append(f"| {i+1} {NAMES[i]} | " + " | ".join(f"{c[0][i]:.1f} / {c[1][i]:.1f}" for c in cases.values()) + " |")
open(os.path.join(D, "results", "frequencies.md"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))

fig, ax = plt.subplots(figsize=(11, 4.6))
x = np.arange(6); wdt = 0.27
for j, (k, (fc, ff)) in enumerate(cases.items()):
    ax.bar(x + (j - 1) * wdt, fc[:6], wdt, label=k)
    for i in range(6):
        ax.text(x[i] + (j - 1) * wdt, fc[i] + 8, f"{fc[i]:.0f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([f"{i+1}次\n{n}" for i, n in enumerate(NAMES)], fontsize=9)
ax.set_ylabel("固有振動数 [Hz]（CalculiX）"); ax.legend(); ax.grid(alpha=0.3, axis="y")
ax.set_title("ばねのつなぎ方による違い（合計の並進剛性・回転剛性は同じ）", loc="left")
fig.tight_layout(); fig.savefig(os.path.join(D, "results", "mpc_compare.png")); plt.close(fig)
