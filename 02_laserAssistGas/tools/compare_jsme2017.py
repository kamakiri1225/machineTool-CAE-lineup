#!/usr/bin/env python3
"""Compare the 3-D kerf cases with Matsumura et al. (JSME 2017) Fig. 5 and Fig. 6.
Paper values are read off the printed figures by eye (about +-3 % of full scale).
  Fig. 6: wall shear stress tau on the cutting front vs p0
  Fig. 5: velocity V and static / dynamic / total pressure at the origin (slot centre on the plate top) vs p0
CFD: tau = mean |wallShearStress| on the cutting front inside the slot (tools/kerf_metrics.py definition);
     origin = (x, y, z) = (-W/2, 0, 0): centre of the slot width, on the plate top surface.
usage: compare_jsme2017.py <out.png> <case> [<case> ...]"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy as n
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17,
                     "legend.fontsize": 15})
out, cases = sys.argv[1], sys.argv[2:]
# ---- paper, read from the figures (p0 in MPa)
P0 = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 8.0])
TAU = np.array([1.1, 1.6, 2.8, 3.9, 4.9, 8.5])                  # kPa, Fig. 6
V = np.array([262, 310, 313, 308, 315, 302])                     # m/s, Fig. 5
PT = np.array([0.44, 0.85, 1.68, 2.40, 3.25, 6.06])              # MPa, Fig. 5 (total)
PS = np.array([0.30, 0.50, 1.00, 1.45, 1.90, 3.72])              # MPa, Fig. 5 (static)
PD = np.array([0.14, 0.34, 0.68, 0.92, 1.30, 2.34])              # MPa, Fig. 5 (dynamic)
rows = []
for case in cases:
    case = os.path.abspath(case); C = json.load(open(os.path.join(case, "caseParams.json")))
    rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
    tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.SetPatchArrayStatus("patch/plate", 1)
    rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
    mb = rd.GetOutput(); m = mb.GetBlock(0); it = mb.NewIterator(); it.InitTraversal(); pl = None
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == "plate":
            pl = it.GetCurrentDataObject()
        it.GoToNextItem()
    W, T = C["W"] * 1e-3, C["T"] * 1e-3
    cc = vtk.vtkCellCenters(); cc.SetInputData(pl); cc.Update(); c = n(cc.GetOutput().GetPoints().GetData())
    sel = (np.abs(c[:, 0]) < 2e-5) & (c[:, 1] < 0) & (c[:, 1] > -T) & (c[:, 2] < W / 2)
    tau = np.linalg.norm(n(pl.GetCellData().GetArray("wallShearStress"))[sel], axis=1).mean()
    pts = vtk.vtkPoints(); pts.InsertNextPoint(-W / 2, 1e-6, 1e-6); pd = vtk.vtkPolyData(); pd.SetPoints(pts)
    pr = vtk.vtkProbeFilter(); pr.SetInputData(pd); pr.SetSourceData(m); pr.Update(); q = pr.GetOutput().GetPointData()
    U = np.linalg.norm(n(q.GetArray("U"))[0]); Tt = n(q.GetArray("T"))[0]; p = n(q.GetArray("p"))[0]
    Y = n(q.GetArray("O2"))[0]; R = 8.314462618 * (Y / 0.0319988 + (1 - Y) / 0.0280134); rho = p / (R * Tt)
    g = 1.4; Ma = U / np.sqrt(g * R * Tt)
    pt = p * (1 + 0.5 * (g - 1) * Ma**2) ** (g / (g - 1))
    rows.append(dict(p0=C["p0"] / 1e6, tau=tau / 1e3, V=U, ps=(p - 101325) / 1e6, pd=0.5 * rho * U**2 / 1e6, pt=(pt - 101325) / 1e6))
    print(f"p0 {C['p0']/1e6:.1f} MPa  tau {tau/1e3:.2f} kPa  V {U:.0f} m/s  Ma {Ma:.2f}  ps {rows[-1]['ps']:.2f}  pd {rows[-1]['pd']:.2f}  pt {rows[-1]['pt']:.2f} MPa(gauge)")
cp = np.array([r["p0"] for r in rows])
fig, ax = plt.subplots(1, 3, figsize=(18, 6.2))
a = ax[0]
a.plot(P0, TAU, "o-", color="#C62828", lw=2.5, ms=8, mfc="white", label="論文（空気、Fluent）")
a.plot(cp, [r["tau"] for r in rows], "s", color="#1565C0", ms=13, label="今回（N₂、OpenFOAM）")
a.set_xlabel("ノズル入口圧力 p₀ [MPa]"); a.set_ylabel("せん断応力 τ [kPa]"); a.set_title("(a) 切断の最前面のせん断応力")
a = ax[1]
a.plot(P0, V, "o-", color="#EF6C00", lw=2.5, ms=8, mfc="white", label="論文")
a.plot(cp, [r["V"] for r in rows], "s", color="#EF6C00", ms=13, label="今回")
a.set_ylim(0, 450); a.set_xlabel("ノズル入口圧力 p₀ [MPa]"); a.set_ylabel("速度 V [m/s]"); a.set_title("(b) 原点の速度")
a = ax[2]
for key, lab, col in (("pt", "全圧 pₜ", "#2E7D32"), ("ps", "静圧 pₛ", "#1565C0"), ("pd", "動圧 p_d", "#C62828")):
    a.plot(P0, {"pt": PT, "ps": PS, "pd": PD}[key], "o-", color=col, lw=2.5, ms=8, mfc="white", label=f"{lab}（論文）")
    a.plot(cp, [r[key] for r in rows], "s", color=col, ms=13, label=f"{lab}（今回）")
a.set_xlabel("ノズル入口圧力 p₀ [MPa]"); a.set_ylabel("圧力（ゲージ）[MPa]"); a.set_title("(c) 原点の圧力")
for a in ax:
    a.set_xlim(0, 8.5); a.grid(alpha=0.3); a.legend(frameon=False, loc="upper left")
ax[1].legend(frameon=False, loc="lower right"); ax[2].legend(frameon=False, loc="upper left", fontsize=13, ncol=2)
fig.tight_layout(); fig.savefig(out, dpi=100, facecolor="white", bbox_inches="tight"); print(out)
