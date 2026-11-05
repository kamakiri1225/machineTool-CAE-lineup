#!/usr/bin/env python3
"""Quick look of a 2-D (one-cell-thick) case made by make_case_from_yaml.py: Mach number and O2 mole fraction.
usage: plot_slice2d.py <case> <out.png>"""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
m = rd.GetOutput().GetBlock(0); b = m.GetBounds()
pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, 0.5 * (b[4] + b[5])); pl.SetNormal(0, 0, 1)
ct = vtk.vtkCutter(); ct.SetInputData(m); ct.SetCutFunction(pl)
tr = vtk.vtkTriangleFilter(); tr.SetInputConnection(ct.GetOutputPort()); tr.Update(); pd = tr.GetOutput()
xy = vtk_to_numpy(pd.GetPoints().GetData())[:, :2] * 1e3
tri = mtri.Triangulation(xy[:, 0], xy[:, 1], vtk_to_numpy(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:])
g = lambda n: vtk_to_numpy(pd.GetPointData().GetArray(n))
U, T, Y = np.linalg.norm(g("U"), axis=1), g("T"), g("O2")
Ma = U / np.sqrt(1.4 * 8.314 / 0.02885 * T)
X = (Y / 32.0) / (Y / 32.0 + (1 - Y) / 28.01) * 100
fig, ax = plt.subplots(1, 2, figsize=(16, 5.2))
for a, f, cm, lo, hi, lab in ((ax[0], Ma, "jet", 0, 1.6, "マッハ数"), (ax[1], X, "RdYlBu", 0, 100, "O₂ モル分率 [%]")):
    tc = a.tripcolor(tri, f, shading="gouraud", cmap=cm, vmin=lo, vmax=hi); a.set_aspect("equal")
    a.set_title(lab); a.set_xlabel("x [mm]"); fig.colorbar(tc, ax=a, fraction=0.03, pad=0.02)
ax[0].set_ylabel("y [mm]")
fig.suptitle(f"YAML の例題（反復 {t:g} 回）：N₂ の主噴流（下）と O₂ の随伴流（上）", fontsize=19)
fig.tight_layout(); fig.savefig(out, dpi=100, facecolor="white", bbox_inches="tight"); print(out)
