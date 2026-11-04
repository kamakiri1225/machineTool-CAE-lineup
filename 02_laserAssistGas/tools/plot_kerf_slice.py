#!/usr/bin/env python3
"""Slot-centre plane (z = 0) of the 3-D kerf cases: velocity, static pressure, O2 (air) fraction,
laid out like Fig. 3 of Matsumura et al. (JSME 2017). Also the velocity / O2 along the slot centre line.
usage: plot_kerf_slice.py <out.png> <case1> [<case2> ...]   (latest time of each case)
"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np, vtk
from matplotlib.patches import Polygon
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16,
                     "axes.titlesize": 18, "axes.labelsize": 16, "xtick.labelsize": 14, "ytick.labelsize": 14})
out, cases = sys.argv[1], sys.argv[2:]


def load(case):
    P = json.load(open(os.path.join(case, "caseParams.json")))
    rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
    tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOn()
    rd.UpdateTimeStep(t); rd.Update()
    m = rd.GetOutput().GetBlock(0)
    pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, 1e-6); pl.SetNormal(0, 0, 1)
    ct = vtk.vtkCutter(); ct.SetInputData(m); ct.SetCutFunction(pl); ct.Update()
    tr = vtk.vtkTriangleFilter(); tr.SetInputConnection(ct.GetOutputPort()); tr.Update()
    pd = tr.GetOutput()
    xy = vtk_to_numpy(pd.GetPoints().GetData())[:, :2] * 1e3
    tris = vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    g = lambda n: vtk_to_numpy(pd.GetPointData().GetArray(n))
    U = g("U"); p = g("p"); Y = g("O2")
    X = (Y / 32.0) / (Y / 32.0 + (1 - Y) / 28.01) * 100
    return P, t, mtri.Triangulation(xy[:, 0], xy[:, 1], tris), np.linalg.norm(U, axis=1), p / 101325.0, X, xy


data = [load(os.path.abspath(c)) for c in cases]
n = len(data)
if n == 1:      # one case: the three fields side by side (avoids a tall, empty figure)
    fig, ax1 = plt.subplots(1, 3, figsize=(17, 6.2), squeeze=False); axes = ax1.T
else:
    fig, axes = plt.subplots(3, n, figsize=(4.6 * n + 1.5, 12.5), squeeze=False)


def solids(ax, P):
    kw = dict(facecolor="white", edgecolor="#333333", lw=1.0, zorder=5)
    H, T, DE, RN = P["H"], P["T"], P["DE"], P["RNOZ"]
    ax.add_patch(Polygon([[0, 0], [5, 0], [5, -T], [0, -T]], **kw))          # uncut plate (x > 0)
    for s in (1, -1):
        ax.add_patch(Polygon([[s * DE / 2, H], [s * RN, H], [s * RN, H + P["LN"]], [s * P["DIN"] / 2, H + P["LN"]]], **kw))


vmaxU = 750
for j, (P, t, tri, V, pa, X, xy) in enumerate(data):
    pmax = max(1.0, P["p0"] / 1e5 + 1)
    for i, (f, cmap, vmin, vmax, lab) in enumerate(((V, "jet", 0, vmaxU, "速度 V [m/s]"),
                                                   (pa, "jet", 0, pmax, "静圧 pₛ [atm]"),
                                                   (X, "RdYlBu", 0, 21, "O₂ 濃度 [%]（空気 21 %）"))):
        ax = axes[i, j]
        tc = ax.tripcolor(tri, f, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax)
        solids(ax, P)
        ax.set_xlim(-2.5, 2.0); ax.set_ylim(-P["T"] - 1.5, P["H"] + 1.5); ax.set_aspect("equal")
        if i == 0 and n > 1:
            ax.set_title(f"p₀ = {P['p0']/1e6:g} MPa（反復 {t:g}）")
        if n == 1:
            ax.set_title(f"{lab}", fontsize=17); ax.set_xlabel("x [mm]（左が切断済み）")
            if i == 0: ax.set_ylabel("y [mm]")
        else:
            if j == 0: ax.set_ylabel(f"{lab}\ny [mm]")
            if i == 2: ax.set_xlabel("x [mm]（左が切断済み）")
        cb = fig.colorbar(tc, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=13)
fig.suptitle("切断溝の中心面（z = 0）：速度・静圧・周りの空気から入り込んだ O₂（N₂ アシストガス）" + (f"　p₀ = {data[0][0]['p0']/1e6:g} MPa（反復 {data[0][1]:g}）" if n == 1 else ""), fontsize=19)
fig.tight_layout(rect=(0, 0, 1, 0.93 if n == 1 else 0.96))
fig.savefig(out, dpi=100, facecolor="white", bbox_inches="tight")
print(out)

# centre-line profiles: x = -0.15 mm (middle of the slot behind the front) along y
fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 7))
for P, t, tri, V, pa, X, xy in data:
    yl = np.linspace(-P["T"] - 1.5, P["H"], 300); xl = np.full_like(yl, -0.15)
    iv = mtri.LinearTriInterpolator(tri, V)(xl, yl); ix = mtri.LinearTriInterpolator(tri, X)(xl, yl)
    a1.plot(iv, yl, lw=2.5, label=f"p₀ = {P['p0']/1e6:g} MPa"); a2.plot(ix, yl, lw=2.5, label=f"p₀ = {P['p0']/1e6:g} MPa")
for a in (a1, a2):
    for yy, lab in ((P["H"], "ノズル出口"), (0, "板の上面"), (-P["T"], "板の下面")):
        a.axhline(yy, color="#888", lw=1, ls=":")
        a.text(1.02, yy, lab, transform=a.get_yaxis_transform(), fontsize=14, color="#555", va="center", clip_on=False)
    a.set_ylabel("y [mm]"); a.legend(frameon=False, loc="lower right")
a1.set_xlabel("速度 V [m/s]"); a2.set_xlabel("O₂ 濃度 [%]"); a2.set_xlim(-0.5, 21)
a1.set_title("切断溝の中（x = −0.15 mm）の速度"); a2.set_title("切断溝の中の O₂ 濃度（空気の巻き込み）")
fig.tight_layout()
o2 = out.replace(".png", "_profile.png"); fig.savefig(o2, dpi=100, facecolor="white", bbox_inches="tight"); print(o2)
