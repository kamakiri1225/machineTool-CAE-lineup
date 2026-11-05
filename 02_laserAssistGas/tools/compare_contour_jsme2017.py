#!/usr/bin/env python3
"""Slot-centre plane (z = 0) contours laid out like Fig. 3 of Matsumura et al. (JSME 2017):
 (a) velocity 0-750 m/s, one colour bar at the right;  (b) static pressure 0-(p0+1) atm, a colour bar above each panel;
 p0 under each panel, banded rainbow colours, the same window around the nozzle.
usage: compare_contour_jsme2017.py <out.png> <case> [<case> ...]
       [--paper <page3.png>]   put the paper's Fig. 3 rows above the CFD rows (internal use; contains the paper's figure)
       cases whose latest iteration is below 4000 or that stalled are labelled with their iteration."""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri, numpy as np, vtk
from matplotlib.patches import Polygon
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy as n
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16})
args = sys.argv[1:]; paper = None
if "--paper" in args:
    i = args.index("--paper"); paper = args[i + 1]; del args[i:i + 2]
out, cases = args[0], args[1:]
XL, YL = (-7.0, 6.0), (-5.8, 4.0)
STALLED = {4.0: "反復 4,500 で停止（参考）", 8.0: "反復 4,500 で停止（参考）"}
# paper Fig. 3 panels on page 3 rendered at 200 dpi (left, top, right, bottom)
PBOX = {"V": {0.5: (158, 295, 447, 515), 2.0: (458, 295, 750, 515), 4.0: (760, 295, 1055, 515), 8.0: (1065, 295, 1355, 515)},
        "p": {0.5: (158, 620, 442, 905), 2.0: (453, 620, 747, 905), 4.0: (757, 620, 1052, 905), 8.0: (1062, 620, 1350, 905)}}
PVBAR = (1370, 300, 1460, 500)          # the paper's velocity colour bar


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


data = sorted([load(os.path.abspath(c)) for c in cases], key=lambda d: d[0]["p0"])
nc = len(data)
rows = ["pV", "V", "pp", "p"] if paper else ["V", "p"]
hr = {"pV": 1.0, "V": 1.0, "pp": 1.25, "p": 1.25}
fig = plt.figure(figsize=(4.0 * nc + 1.6, sum(3.45 * hr[r] for r in rows) + 0.9))
gs = fig.add_gridspec(len(rows), nc + 1, width_ratios=[1] * nc + [0.09], height_ratios=[hr[r] for r in rows],
                      hspace=0.42, wspace=0.08)
page = Image.open(paper).convert("RGB") if paper else None
first = {}
for i, r in enumerate(rows):
    for j, (P, t, tri, V, ps) in enumerate(data):
        p0 = P["p0"] / 1e6
        ax = fig.add_subplot(gs[i, j]); ax.set_xticks([]); ax.set_yticks([])
        if j == 0: first[i] = ax
        note = STALLED.get(p0, "") if t < 4000 or p0 in STALLED else ""
        if r in ("pV", "pp"):
            ax.imshow(page.crop(PBOX[r[1].replace("p", "p") if r == "pp" else "V"][p0])); ax.set_axis_off()
            ax.set_title(f"論文　p₀ = {p0:g} MPa", fontsize=15, y=-0.16)
            continue
        f, hi = (V, 750) if r == "V" else (ps, round(p0 * 10) + 1)   # the paper uses p0 [MPa] x 10 + 1 atm
        cs = ax.tricontourf(tri, np.clip(f, 0, hi), levels=np.linspace(0, hi, 21), cmap="jet")
        solids(ax, P); ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect("equal")
        lab = ("今回　" if paper else "") + f"p₀ = {p0:g} MPa"
        ax.set_title(lab + (f"\n{note}" if note else ""), fontsize=15, y=-0.16 if not note else -0.27, color="#B71C1C" if note else "black")
        if r == "p":   # colour bar above each pressure panel, as in the paper
            cax = ax.inset_axes([0.08, 1.04, 0.72, 0.07])
            cb = fig.colorbar(cs, cax=cax, orientation="horizontal", ticks=[0, hi]); cb.ax.xaxis.set_ticks_position("top")
            cb.ax.tick_params(labelsize=14); cax.text(1.04, 0.5, "pₛ [atm]", transform=cax.transAxes, va="center", fontsize=14)
        if r == "V" and j == nc - 1:
            cax = fig.add_subplot(gs[i, nc]); cb = fig.colorbar(cs, cax=cax, ticks=[0, 325, 750]); cb.ax.tick_params(labelsize=14)
            cax.set_title("V [m/s]", fontsize=14)
    if r in rows[:1] or r == ("pp" if paper else "p"):
        a0 = first[i]
        txt = ("(a) 速度" + ("　上段：論文（松村ら 2017、空気、Fluent）／下段：今回（N₂、OpenFOAM）" if paper else "")) if i == 0 else \
              ("(b) 静圧" + ("　上段：論文／下段：今回" if paper else ""))
        a0.text(0.0, (1.06 if paper else 1.10) if i == 0 else (1.06 if paper else 1.42), txt, transform=a0.transAxes, fontsize=18, ha="left", va="bottom")
title = ("切断溝の中心面（z = 0）の速度と静圧：論文の図 3 と今回の計算の比較" if paper else
         "今回の計算（N₂、OpenFOAM）：切断溝の中心面（z = 0）の速度と静圧（論文の図 3 と同じ並び・色・範囲）")
fig.suptitle(title, fontsize=19, y=0.94 if paper else 1.0)
fig.savefig(out, dpi=90, facecolor="white", bbox_inches="tight"); print(out)
