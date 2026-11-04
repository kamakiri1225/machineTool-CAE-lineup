#!/usr/bin/env python3
"""Front view of the symmetry plane (slot centre plane z = 0) of the 3-D kerf model: fluid region + mesh,
with the nozzle and the plate drawn as the empty (solid) parts. usage: render_symplane.py <case> <out.png>"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from matplotlib.collections import PolyCollection
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
P = json.load(open(os.path.join(case, "caseParams.json")))
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("patch/symmetryPlane", 1); rd.Update()
it = rd.GetOutput().NewIterator(); it.InitTraversal(); pd = None
while not it.IsDoneWithTraversal():
    pd = it.GetCurrentDataObject(); it.GoToNextItem()
g = vtk.vtkGeometryFilter(); g.SetInputData(pd); g.Update(); pd = g.GetOutput()
pts = vtk_to_numpy(pd.GetPoints().GetData()) * 1e3
polys = []
for i in range(pd.GetNumberOfCells()):
    ids = pd.GetCell(i).GetPointIds()
    polys.append(pts[[ids.GetId(k) for k in range(ids.GetNumberOfIds())]][:, :2])
fig, axes = plt.subplots(1, 2, figsize=(17, 8.6), gridspec_kw={"width_ratios": [1.15, 1]})
for ax, (xl, yl, title) in zip(axes, (((-4, 4), (-P["T"] - 2, P["H"] + P["LN"]), "対称面（溝の中心面 z = 0）全体"),
                                     ((-1.6, 1.6), (-P["T"] - 0.5, P["H"] + 0.8), "ノズル出口〜切断溝の拡大"))):
    ax.set_facecolor("#9a9a9a")
    ax.add_collection(PolyCollection(polys, facecolors="#d8ecff", edgecolors="#3a6ea5", linewidths=0.25))
    ax.set_xlim(*xl); ax.set_ylim(*yl); ax.set_aspect("equal")
    ax.set_xlabel("x [mm]（切断方向）"); ax.set_ylabel("y [mm]"); ax.set_title(title)
ax = axes[0]
H, T, LN = P["H"], P["T"], P["LN"]
kw = dict(fontsize=17, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888888", alpha=0.92))
ax.annotate("入口（全圧 p₀）", xy=(0, H + LN), xytext=(1.2, H + LN - 0.6), arrowprops=dict(arrowstyle="-|>"), **kw)
ax.text(-3.6, H + 1.5, "ノズル\n（灰色）", **kw)
ax.text(1.8, -T / 2 - 0.15, "鋼板（灰色）", **kw)
ax.text(-3.7, -T / 2 - 0.15, "切断溝（溝の中\nなので流体）", **kw)
ax.text(-3.7, -T - 1.5, "板の下（大気）", **kw)
ax.annotate("", xy=(3.6, -T - 1.2), xytext=(1.6, -T - 1.2), arrowprops=dict(arrowstyle="-|>", color="#1565C0", lw=3))
ax.text(2.6, -T - 1.0, "切断方向", color="#1565C0", ha="center", fontsize=17)
fig.suptitle(f"run004 のモデル：対称面のメッシュ（水色＝流体、灰色＝ノズル・鋼板）　セル数 {pd.GetNumberOfCells()} 面（対称面上）", fontsize=19)
fig.tight_layout(); fig.savefig(out, dpi=100, facecolor="white"); print(out)
