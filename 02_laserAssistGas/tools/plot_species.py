#!/usr/bin/env python3
"""O2 / N2 composition of the assist-gas case (axisymmetric, mirrored cross-section).
Left: O2 mole fraction (purity) field; right: O2 purity along the axis and along the kerf.
usage: python3 plot_species.py <caseDir> [time]
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
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
tv = rd.GetTimeValues(); t = float(sys.argv[2]) if len(sys.argv) > 2 else tv.GetValue(tv.GetNumberOfTuples() - 1)
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOff()
rd.UpdateTimeStep(t); rd.Update()
mesh = rd.GetOutput().GetBlock(0)
cc = vtk.vtkCellCenters(); cc.SetInputData(mesh); cc.Update()
xyz = vtk_to_numpy(cc.GetOutput().GetPoints().GetData()) * 1e3
Y = vtk_to_numpy(mesh.GetCellData().GetArray("O2"))
X = (Y / 31.9988) / (Y / 31.9988 + (1 - Y) / 28.0134) * 100          # O2 mole (volume) fraction [%]
r, y = np.hypot(xyz[:, 0], xyz[:, 2]), xyz[:, 1]
RK, RN, RIN, RNOZ, ROUT, H, LN, T = (P[k] for k in ("RK", "RN", "RIN", "RNOZ", "ROUT", "H", "LN", "T"))
Rm = np.concatenate([r, -r]); Ym = np.concatenate([y, y])
tri = mtri.Triangulation(Rm, Ym)
xm = Rm[tri.triangles].mean(1); ym = Ym[tri.triangles].mean(1); a = np.abs(xm)
solid = ((ym < 0) & (ym > -T) & (a > RK)) | ((ym > H) & (a > RN + (RIN - RN) * (ym - H) / LN) & (a < RNOZ))
edge = np.max(np.hypot(np.diff(Rm[tri.triangles][:, [0, 1, 2, 0]], axis=1), np.diff(Ym[tri.triangles][:, [0, 1, 2, 0]], axis=1)), axis=1)
tri.set_mask(solid | (edge > 0.6))

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16, 8.2), gridspec_kw={"width_ratios": [1.25, 1]})
tc = ax.tripcolor(tri, np.concatenate([X, X]), shading="gouraud", cmap="RdYlBu_r", vmin=20, vmax=100)
ax.tricontour(tri, np.concatenate([X, X]), levels=[50, 90], colors="k", linewidths=1.0)
kw = dict(facecolor="#9a9a9a", edgecolor="#444444", lw=1.2, zorder=5)
for s in (1, -1):
    ax.add_patch(Polygon([[s * RK, 0], [s * ROUT, 0], [s * ROUT, -T], [s * RK, -T]], **kw))
    ax.add_patch(Polygon([[s * RN, H], [s * RNOZ, H], [s * RNOZ, H + LN], [s * RIN, H + LN]], **kw))
ax.set_xlim(-4, 4); ax.set_ylim(-T - 2.5, H + LN); ax.set_aspect("equal")
ax.set_xlabel("r [mm]"); ax.set_ylabel("y [mm]")
ax.set_title("O₂ の体積割合（純度）[%]")
cb = fig.colorbar(tc, ax=ax, fraction=0.04, pad=0.02); cb.set_label("O₂ 体積割合 [%]（残りは N₂）")

def profile(mask):
    yy, xx = y[mask], X[mask]
    bins = np.linspace(-T - 2.5, H + LN, 160)
    idx = np.digitize(yy, bins)
    yc = [yy[idx == k].mean() for k in np.unique(idx)]
    xc = [xx[idx == k].mean() for k in np.unique(idx)]
    return np.array(xc), np.array(yc)


xc, yc = profile(r < 0.5 * RK)
ax2.plot(xc, yc, color="#2F6FB0", lw=3, label="中心軸上")
xc, yc = profile(np.abs(r - 0.85 * RK) < 0.15 * RK)
ax2.plot(xc, yc, color="#C62828", lw=2.5, ls="--", label="溝の壁の近く（r = 0.85×溝半径）")
for yy, lab in ((H, "ノズル出口"), (0, "板の上面"), (-T, "板の下面")):
    ax2.axhline(yy, color="#888888", lw=1, ls=":")
    ax2.text(1.02, yy, lab, transform=ax2.get_yaxis_transform(), fontsize=15, color="#555555", ha="left", va="center", clip_on=False)
ax2.set_ylim(-T - 2.5, H + LN); ax2.set_xlim(20, 101)
ax2.set_xlabel("O₂ 体積割合 [%]"); ax2.set_ylabel("y [mm]")
ax2.set_title("中心軸・溝の中の O₂ 純度")
ax2.legend(loc="lower left", frameon=False)
gas = "O₂" if P.get("YO2_in", 1) > 0.5 else "空気"
fig.suptitle(f"{gas} {P['p0_gauge']/1e6:.2f} MPa(G)：周りの空気（O₂ 21 %、N₂ 79 %）がどれだけ混ざるか（反復 {t:g}）", fontsize=19)
fig.tight_layout()
fn = os.path.join(case, "images", f"axisym_O2purity_{int(t)}.png")
fig.savefig(fn, dpi=110, bbox_inches="tight", facecolor="white")
ent = X[(y < 0.05) & (y > -0.05) & (r < RK)]
ex = X[(y < -T + 0.05) & (y > -T - 0.05) & (r < RK)]
print(fn, "O2 at kerf entrance %.2f-%.2f %%, at kerf exit %.2f-%.2f %%" % (ent.min(), ent.max(), ex.min(), ex.max()))
