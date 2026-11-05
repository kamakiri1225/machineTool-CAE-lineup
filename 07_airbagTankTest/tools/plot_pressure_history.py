#!/usr/bin/env python3
"""Pressure time history: CFD vs the 0-D comparison model with the SAME constant cp as the CFD
(and optionally the 0-D model with temperature-dependent JANAF cp, plus a JANAF CFD case).
usage: plot_pressure_history.py <out.png> <cfdCase> [<cfdCaseJanaf>]
0-D results are read from zeroD_cpCFD and zeroD_janaf next to the case."""
import glob, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17, "axes.titlesize": 19, "axes.labelsize": 17,
                     "legend.fontsize": 15, "xtick.labelsize": 15, "ytick.labelsize": 15})
out, cases = sys.argv[1], sys.argv[2:]
root = os.path.dirname(os.path.abspath(cases[0]))
P0 = 101325.0


def load(case, name):
    fs = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "volFieldValue.dat")), key=lambda f: float(f.split("/")[-2]))
    return np.vstack([np.loadtxt(f, comments="#") for f in fs]) if fs else None


zs = [(os.path.join(root, "zeroD_cpCFD"), "0 次元モデル（比熱一定、CFD と同じ値）", "k", "--")]
if os.path.exists(os.path.join(root, "zeroD_janaf", "zero_d_result.csv")):
    zs.append((os.path.join(root, "zeroD_janaf"), "0 次元モデル（比熱の温度依存あり）", "#6A1B9A", ":"))
labs = ["CFD（比熱一定）", "CFD（比熱の温度依存あり）"]; cols = [("#C62828", "#1565C0"), ("#EF6C00", "#00838F")]
fig, ax = plt.subplots(1, 2, figsize=(17, 6.8))
for zd, lab, c, ls in zs:
    z = np.loadtxt(os.path.join(zd, "zero_d_result.csv"), delimiter=",", skiprows=1)
    ax[0].plot(z[:, 0] * 1e3, z[:, 1] / 1e6, color=c, ls=ls, lw=2.5, label=lab)
    ax[1].plot(z[:, 0] * 1e3, (z[:, 3] - P0) / 1e3, color=c, ls=ls, lw=2.5, label=lab)
for k, case in enumerate(cases):
    ch, tk = load(case, "chamberAvg"), load(case, "tankAvg")
    if ch is not None:
        ax[0].plot(ch[:, 0] * 1e3, ch[:, 1] / 1e6, color=cols[k][0], lw=3, label=labs[k])
    if tk is not None:
        ax[1].plot(tk[:, 0] * 1e3, (tk[:, 1] - P0) / 1e3, color=cols[k][1], lw=3, label=labs[k] + "（タンクの平均）")
ax[0].set_ylabel("燃焼室の圧力 [MPa]"); ax[0].set_title("インフレーターの燃焼室")
ax[1].set_ylabel("タンクのゲージ圧 [kPa]"); ax[1].set_title("60 L タンクの圧力（タンク試験）")
for a in ax:
    a.set_xlabel("時間 [ms]"); a.grid(alpha=0.3); a.set_xlim(0, 60)
ax[0].legend(frameon=False, loc="upper right", fontsize=14); ax[1].legend(frameon=False, loc="lower right", fontsize=14)
fig.suptitle("圧力の時刻歴：CFD と 0 次元モデルの比較（壁は断熱）", fontsize=20)
fig.tight_layout(); fig.savefig(out, dpi=90, facecolor="white"); print(out)
