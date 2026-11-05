#!/usr/bin/env python3
"""Contours on the slot-centre plane (z = 0) drawn like Fig. 3 of Matsumura et al. (JSME 2017):
banded rainbow colours, velocity 0-750 m/s and static pressure 0-(p0+1) atm, the same window around the nozzle.
usage: compare_contour_jsme2017.py <out.png> <case> [<case> ...]
       [--paper <page3.png>]   optional: put crops of the paper's Fig. 3 above each CFD panel (for internal use)"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri, numpy as np, vtk
from matplotlib.patches import Polygon
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy as n
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 17})
args = sys.argv[1:]; paper = None
if "--paper" in args:
    i = args.index("--paper"); paper = args[i + 1]; del args[i:i + 2]
out, cases = args[0], args[1:]
XL, YL = (-7.0, 6.0), (-5.8, 4.0)
# paper Fig. 3 panel boxes on page 3 rendered at 200 dpi (left, top, right, bottom); key = p0 [MPa]
PBOX = {"V": {0.5: (158, 295, 447, 515), 2.0: (458, 295, 750, 515), 4.0: (760, 295, 1055, 515), 8.0: (1065, 295, 1355, 515)},
        "p": {0.5: (158, 620, 442, 905), 2.0: (453, 620, 747, 905), 4.0: (757, 620, 1052, 905), 8.0: (1062, 620, 1350, 905)}}


def load(case):
    P = json.load(open(os.path.join(case, "caseParams.json")))
    rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
    tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
    pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, 1e-6); pl.SetNormal(0, 0, 1)
    ct = vtk.vtkCutter(); ct.SetInputData(rd.GetOutput().GetBlock(0)); ct.SetCutFunction(pl)
    tr = vtk.vtkTriangleFilter(); tr.SetInputConnection(ct.GetOutputPort()); tr.Update(); pd = tr.GetOutput()
    xy = n(pd.GetPoints().GetData())[:, :2] * 1e3
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1], n(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:])
    g = lambda k: n(pd.GetPointData().GetArray(k))
    return P, t, tri, np.linalg.norm(g("U"), axis=1), g("p") / 101325.0


def solids(ax, P):
    kw = dict(facecolor="white", edgecolor="#333", lw=1.0, zorder=5)
    H, T, DE, RN, LN, DIN = P["H"], P["T"], P["DE"], P["RNOZ"], P["LN"], P["DIN"]
    ax.add_patch(Polygon([[0, 0], [XL[1] + 1, 0], [XL[1] + 1, -T], [0, -T]], **kw))
    for s in (1, -1):
        ax.add_patch(Polygon([[s * DE / 2, H], [s * RN, H], [s * RN, H + LN + 1], [s * DIN / 2, H + LN + 1], [s * DIN / 2, H + LN]], **kw))


data = [load(os.path.abspath(c)) for c in cases]
nc = len(data); nr = 4 if paper else 2
fig, axes = plt.subplots(nr, nc, figsize=(5.2 * nc + 1.5, 4.1 * nr), squeeze=False)
page = Image.open(paper).convert("RGB") if paper else None
for j, (P, t, tri, V, ps) in enumerate(data):
    p0 = P["p0"] / 1e6; pmax = round(p0 * 1e6 / 101325) + 1
    rows = [("V", V, 0, 750, "速度 V [m/s]"), ("p", ps, 0, pmax, "静圧 pₛ [atm]")]
    for i, (key, f, lo, hi, lab) in enumerate(rows):
        if paper:
            ax = axes[2 * i, j]; ax.imshow(page.crop(PBOX[key][p0])); ax.set_axis_off()
            ax.set_title(f"論文 図 3：{lab}（p₀ = {p0:g} MPa）", fontsize=15)
            ax = axes[2 * i + 1, j]
        else:
            ax = axes[i, j]
        lev = np.linspace(lo, hi, 21)
        cs = ax.tricontourf(tri, np.clip(f, lo, hi), levels=lev, cmap="jet")
        solids(ax, P); ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"今回：{lab}（p₀ = {p0:g} MPa）" if paper else f"今回の計算：{lab}　p₀ = {p0:g} MPa", fontsize=15)
        cb = fig.colorbar(cs, ax=ax, fraction=0.04, pad=0.02, ticks=[lo, (lo + hi) / 2, hi]); cb.ax.tick_params(labelsize=14)
fig.suptitle("切断溝の中心面（z = 0）の速度と静圧：論文（松村ら 2017、空気、Fluent）と今回（N₂、OpenFOAM）" if paper else
             "今回の計算（N₂、OpenFOAM）だけの図：切断溝の中心面（z = 0）の速度と静圧\n（論文の図は含まない。論文の図 3 と見比べやすいよう、同じ色・同じ範囲で描いた）", fontsize=19)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(out, dpi=90, facecolor="white", bbox_inches="tight"); print(out)
