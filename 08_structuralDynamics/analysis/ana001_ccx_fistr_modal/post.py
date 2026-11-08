"""ana001 の後処理: 固有振動数の比較表とモード形状図"""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm
from frd import read_frd
from plot_modes import plot_modes, plt

D = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(D, "results"), exist_ok=True)

coords, blocks = read_frd(os.path.join(D, "ccx", "machine.frd"))
f_ccx = [b[1] for b in blocks]
f_fis = []
for ln in open(os.path.join(D, "fistr", "0.log")):
    m = re.match(r"\s+(\d+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s", ln)
    if m and len(f_fis) < len(f_ccx) and int(m.group(1)) == len(f_fis) + 1 and "E" in m.group(4):
        f_fis.append(float(m.group(4)))

ids = np.array(sorted(coords)); X = np.array([coords[i] for i in ids])
col = ids[X[:, 2] > 999]; tab = ids[(X[:, 2] > 349) & (X[:, 2] < 351)]


def describe(d):
    U = np.array([d[i] for i in ids]); mx = np.abs(U).max()
    uc = np.abs(np.array([d[i] for i in col])).max(0) / mx
    ut = np.abs(np.array([d[i] for i in tab])).max(0) / mx
    part = "コラム" if uc.max() > ut.max() else "テーブル"
    u = uc if part == "コラム" else ut
    return f"{part} {'XYZ'[int(np.argmax(u))]}"


rows = ["| モード | CalculiX [Hz] | FrontISTR [Hz] | 差 [%] | 主にゆれる所・方向 |", "|---|---|---|---|---|"]
for k, (fc, ff) in enumerate(zip(f_ccx, f_fis)):
    rows.append(f"| {k+1} | {fc:.1f} | {ff:.1f} | {100*(ff-fc)/fc:+.2f} | {describe(blocks[k][2])} |")
open(os.path.join(D, "results", "frequencies.md"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))

model = mm.build()
NAMES = ["コラム 前後(Y)に曲がる", "テーブル 送り(X)方向", "コラム 左右(X)に曲がる",
         "テーブル 横(Y)方向", "テーブル Z軸まわりにねじれる", "テーブル 上下(Z)＋曲げ"]   # 図を見て命名
titles = [f"{k+1}次 {f_ccx[k]:.0f} Hz（{NAMES[k]}）" for k in range(6)]
plot_modes(model, coords, [blocks[k][2] for k in range(6)], titles, os.path.join(D, "results", "mode_shapes.png"))

fig, ax = plt.subplots(figsize=(8, 4))
k = np.arange(1, len(f_ccx) + 1)
ax.bar(k - 0.2, f_ccx, 0.4, label="CalculiX 2.21（C3D20 + SPRING2）")
ax.bar(k + 0.2, f_fis, 0.4, label="FrontISTR 5.9（362 + 511/!SPRING_D）")
ax.set_xlabel("モード次数"); ax.set_ylabel("固有振動数 [Hz]"); ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.savefig(os.path.join(D, "results", "freq_compare.png")); plt.close(fig)
