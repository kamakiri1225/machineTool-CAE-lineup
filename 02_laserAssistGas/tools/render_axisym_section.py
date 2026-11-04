#!/usr/bin/env python3
"""Labeled cross-section (mirrored) of the axisymmetric assist-gas model with its mesh.
usage: render_axisym_section.py <case> <out.png>"""
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
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("patch/front", 1); rd.Update()
it = rd.GetOutput().NewIterator(); it.InitTraversal(); pd = None
while not it.IsDoneWithTraversal():
    pd = it.GetCurrentDataObject(); it.GoToNextItem()
g = vtk.vtkGeometryFilter(); g.SetInputData(pd); g.Update(); pd = g.GetOutput()
pts = vtk_to_numpy(pd.GetPoints().GetData()) * 1e3
r = np.hypot(pts[:, 0], pts[:, 2]); y = pts[:, 1]
polys = []
for i in range(pd.GetNumberOfCells()):
    ids = pd.GetCell(i).GetPointIds(); k = [ids.GetId(j) for j in range(ids.GetNumberOfIds())]
    q = np.c_[r[k], y[k]]; polys += [q, q * [-1, 1]]
RK, RN, RIN, RNOZ, H, LN, T = (P[k] for k in ("RK", "RN", "RIN", "RNOZ", "H", "LN", "T"))
fig, ax = plt.subplots(figsize=(12, 10))
ax.set_facecolor("#9a9a9a")
ax.add_collection(PolyCollection(polys, facecolors="#d8ecff", edgecolors="#3a6ea5", linewidths=0.2))
ax.set_xlim(-6, 6); ax.set_ylim(-T - 2.5, H + LN); ax.set_aspect("equal")
ax.set_xlabel("r [mm]（中心軸からの距離、左右に折り返して表示）"); ax.set_ylabel("y [mm]")
kw = dict(fontsize=18, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888888", alpha=0.95))
ar = dict(arrowstyle="-|>", lw=1.6)
gas = "O₂" if P.get("YO2_in", 1) > 0.5 else "空気"
ax.annotate(f"入口：{gas}、全圧 p₀ = {P['p0_gauge']/1e6:g} MPa(G)", xy=(0, H + LN), xytext=(1.6, H + LN - 0.6), arrowprops=ar, **kw)
ax.text(-5.7, H + LN / 2, "ノズル（先細）", **kw)
ax.annotate(f"ノズル出口 径 {2*RN:g} mm", xy=(RN, H), xytext=(2.0, H + 1.1), arrowprops=ar, **kw)
ax.annotate(f"スタンドオフ {H:g} mm", xy=(-3.5, H / 2), xytext=(-5.7, -0.9), arrowprops=ar, **kw)
ax.text(1.8, -T / 2, f"鋼板（板厚 {T:g} mm）", **kw)
ax.annotate(f"切断溝（丸い穴で近似、径 {2*RK:g} mm）", xy=(RK, -T * 0.7), xytext=(1.8, -T * 0.85), arrowprops=ar, **kw)
ax.text(-5.7, -T - 1.9, "板の下（大気）", **kw)
ax.plot([0, 0], [-T - 2.5, H + LN], "k-.", lw=1); ax.text(0.1, -T - 2.3, "中心軸", fontsize=16)
ax.set_title("軸対称モデル（水色＝流体とメッシュ、灰色＝ノズル・鋼板）", fontsize=20)
fig.tight_layout(); fig.savefig(out, dpi=100, facecolor="white", bbox_inches="tight"); print(out)
