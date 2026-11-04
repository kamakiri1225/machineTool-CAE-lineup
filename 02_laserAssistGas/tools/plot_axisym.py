#!/usr/bin/env python3
"""Plot an axisymmetric assist-gas case: mirrored cross-section (Mach number + streamlines),
with the nozzle and the steel plate drawn in grey.
usage: python3 plot_axisym.py <caseDir> [time]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import vtk
from matplotlib.patches import Polygon
from vtk.util.numpy_support import vtk_to_numpy

plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17,
                     "axes.titlesize": 20, "axes.labelsize": 18, "xtick.labelsize": 15, "ytick.labelsize": 15})

case = os.path.abspath(sys.argv[1])
P = json.load(open(os.path.join(case, "caseParams.json")))
rd = vtk.vtkOpenFOAMReader()
rd.SetFileName(os.path.join(case, "post.foam"))
rd.UpdateInformation()
times = [rd.GetTimeValues().GetValue(i) for i in range(rd.GetTimeValues().GetNumberOfTuples())]
t = float(sys.argv[2]) if len(sys.argv) > 2 else times[-1]
rd.DisableAllPatchArrays()
rd.SetPatchArrayStatus("internalMesh", 1)
rd.CreateCellToPointOff()
rd.UpdateTimeStep(t)
rd.Update()
mesh = rd.GetOutput().GetBlock(0)
cc = vtk.vtkCellCenters()
cc.SetInputData(mesh)
cc.Update()
xyz = vtk_to_numpy(cc.GetOutput().GetPoints().GetData()) * 1e3      # mm
cd = mesh.GetCellData()
Ma = vtk_to_numpy(cd.GetArray("Ma")) if cd.GetArray("Ma") else None
U = vtk_to_numpy(cd.GetArray("U"))
p = vtk_to_numpy(cd.GetArray("p"))
O2 = vtk_to_numpy(cd.GetArray("O2"))
r, y = np.hypot(xyz[:, 0], xyz[:, 2]), xyz[:, 1]
ur, uy = U[:, 0], U[:, 1]
if Ma is None:
    Ma = np.linalg.norm(U, axis=1) / 340.0

# mirror
R = np.concatenate([r, -r]); Y = np.concatenate([y, y])
tri = mtri.Triangulation(R, Y)
# mask triangles that cross solid parts (long edges or centroid inside solid)
RK, RN, RIN, RNOZ, ROUT = P["RK"], P["RN"], P["RIN"], P["RNOZ"], P["ROUT"]
H, LN, T = P["H"], P["LN"], P["T"]


def in_solid(rr, yy):
    a = np.abs(rr)
    plate = (yy < 0) & (yy > -T) & (a > RK)
    cone_r = RN + (RIN - RN) * (yy - H) / LN
    nozzle = (yy > H) & (a > cone_r) & (a < RNOZ)
    return plate | nozzle


xm = R[tri.triangles].mean(axis=1); ym = Y[tri.triangles].mean(axis=1)
edge = np.max(np.hypot(np.diff(R[tri.triangles][:, [0, 1, 2, 0]], axis=1), np.diff(Y[tri.triangles][:, [0, 1, 2, 0]], axis=1)), axis=1)
tri.set_mask(in_solid(xm, ym) | (edge > 0.6))


def solids(ax):
    kw = dict(facecolor="#9a9a9a", edgecolor="#444444", lw=1.2, zorder=5)
    for s in (1, -1):
        ax.add_patch(Polygon([[s * RK, 0], [s * ROUT, 0], [s * ROUT, -T], [s * RK, -T]], **kw))
        ax.add_patch(Polygon([[s * RN, H], [s * RNOZ, H], [s * RNOZ, H + LN], [s * RIN, H + LN]], **kw))


fig, axes = plt.subplots(1, 2, figsize=(16, 8.2), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axes[0]
vmax = max(1.0, float(np.ceil(Ma.max() * 2) / 2))
tc = ax.tripcolor(tri, np.concatenate([Ma, Ma]), shading="gouraud", cmap="turbo", vmin=0, vmax=vmax)
# streamlines on a regular grid
gx = np.linspace(-4, 4, 260); gy = np.linspace(-T - 2.5, H + LN, 300)
GX, GY = np.meshgrid(gx, gy)
interp_r = mtri.LinearTriInterpolator(tri, np.concatenate([ur, -ur]))
interp_y = mtri.LinearTriInterpolator(tri, np.concatenate([uy, uy]))
UR = interp_r(GX, GY).filled(np.nan); UY = interp_y(GX, GY).filled(np.nan)
ax.streamplot(GX, GY, UR, UY, color="white", density=1.6, linewidth=0.8, arrowsize=0.9)
solids(ax)
ax.set_xlim(-4, 4); ax.set_ylim(-T - 2.5, H + LN)
ax.set_aspect("equal")
ax.set_xlabel("r [mm]"); ax.set_ylabel("y [mm]")
ax.set_title("マッハ数と流線（断面を左右に折り返し）")
cb = fig.colorbar(tc, ax=ax, fraction=0.04, pad=0.02); cb.set_label("Mach 数")
ax.text(0, H + LN - 0.3, "ノズル", ha="center", va="top", color="white", fontsize=16, zorder=6)
ax.text(-ROUT * 0.42, -T / 2, "鋼板", ha="center", color="black", fontsize=17, zorder=6)
ax.text(0, -T - 1.2, "切断溝の出口", ha="center", color="white", fontsize=15, zorder=6)

# axis profiles
ax2 = axes[1]
sel = np.abs(r) < 0.5 * RK
o = np.argsort(y[sel])
ax2.plot(Ma[sel][o], y[sel][o], color="#2F6FB0", lw=3, label="Mach 数（中心軸上）")
ax2b = ax2.twiny()
ax2b.plot(p[sel][o] / 1e5, y[sel][o], color="#C62828", lw=2.5, ls="--", label="静圧 [bar]")
for yy, lab in ((H, "ノズル出口"), (0, "板の上面"), (-T, "板の下面")):
    ax2.axhline(yy, color="#888888", lw=1, ls=":")
    ax2.text(1.02, yy, lab, transform=ax2.get_yaxis_transform(), fontsize=15, color="#555555", ha="left", va="center", clip_on=False)
ax2.set_ylim(-T - 2.5, H + LN)
ax2.set_xlabel("Mach 数", color="#2F6FB0"); ax2b.set_xlabel("静圧 [bar]", color="#C62828")
ax2.set_ylabel("y [mm]")
ax2.set_title("中心軸に沿った変化", pad=48)
gas = "O₂" if P.get("YO2_in", 1) > 0.5 else "空気"
fig.suptitle(f"{P.get('title', '')}：{gas} {P['p0_gauge']/1e6:.2f} MPa(G)、ノズル径 {2*RN:.1f} mm、"
             f"スタンドオフ {H:.1f} mm、板厚 {T:.0f} mm、溝幅 {2*RK:.1f} mm（反復 {t:g}）", fontsize=18)
fig.tight_layout()
out = os.path.join(case, "images")
os.makedirs(out, exist_ok=True)
fn = os.path.join(out, f"axisym_Mach_{int(t)}.png")
fig.savefig(fn, dpi=110, bbox_inches="tight", facecolor="white")
print(fn, "Ma max", Ma.max(), "p max", p.max())
