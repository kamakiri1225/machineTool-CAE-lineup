#!/usr/bin/env python3
"""Initial residual history (first solve of each field per iteration) + inlet mass flow, for LTS cases.
usage: plot_residuals.py <out.png> <case> [...]"""
import os, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 18, "axes.labelsize": 16})
out, cases = sys.argv[1], sys.argv[2:]
fields = ["Ux", "Uy", "Uz", "p", "h", "O2", "k", "epsilon"]
fig, axes = plt.subplots(2, len(cases), figsize=(7.5 * len(cases), 10), squeeze=False, gridspec_kw={"height_ratios": [2, 1]})
for j, c in enumerate(cases):
    it, res, seen = 0, {f: [] for f in fields}, set()
    for line in open(os.path.join(c, "log.reactingFoam")):
        if line.startswith("Time ="):
            it = float(line.split()[2]); seen = set()
        m = re.search(r"Solving for (\w+), Initial residual = ([\d.eE+-]+)", line)
        if m and m.group(1) in res and m.group(1) not in seen:
            seen.add(m.group(1)); res[m.group(1)].append((it, float(m.group(2))))
    ax = axes[0, j]
    for f in fields:
        if res[f]:
            a = np.array(res[f]); ax.semilogy(a[:, 0], a[:, 1], lw=1.2, label=f)
    ax.set_title(os.path.basename(c).replace("run004_kerf3D_N2_", "")); ax.set_ylabel("初期残差"); ax.grid(alpha=0.3)
    ax.legend(ncol=4, fontsize=12, loc="upper right")
    d = np.loadtxt(os.path.join(c, "postProcessing/inletFlow/0/surfaceFieldValue.dat"), comments="#")
    axes[1, j].plot(d[:, 0], -2e3 * d[:, 1], lw=2, color="#C62828")
    axes[1, j].set_xlabel("反復回数（LTS）"); axes[1, j].set_ylabel("ノズルの流量 [g/s]"); axes[1, j].grid(alpha=0.3)
fig.tight_layout(); fig.savefig(out, dpi=90, facecolor="white"); print(out)
