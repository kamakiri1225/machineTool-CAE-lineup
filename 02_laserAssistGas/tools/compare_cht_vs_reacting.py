#!/usr/bin/env python3
"""Compare the axisymmetric assist-gas case solved by reactingFoam (plate = fixed-temperature wall) with the conjugate
heat transfer case solved by chtMultiRegionFoam (steel plate solved, cutting front held at TFRONT).
Left half of each picture: reactingFoam, right half: chtMultiRegionFoam (mirror image about the axis).
usage: compare_cht_vs_reacting.py <out_prefix> <reactingCase> <chtCase>
writes <out_prefix>_fields.png, <out_prefix>_axis.png and prints a metrics table"""
import glob, json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.tri as mtri, numpy as np, vtk
from matplotlib.patches import Polygon
from vtk.util.numpy_support import vtk_to_numpy as n
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 18,
                     "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
out, rc, cc = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
P = json.load(open(os.path.join(rc, "caseParams.json")))


def read(case, block, parallel=False):
    rd = (vtk.vtkPOpenFOAMReader() if parallel else vtk.vtkOpenFOAMReader()); rd.SetFileName(os.path.join(case, "post.foam"))
    if parallel: rd.SetCaseType(1)
    rd.UpdateInformation(); tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus(block, 1); rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
    it = rd.GetOutput().NewIterator(); it.InitTraversal(); ap = vtk.vtkAppendFilter()
    while not it.IsDoneWithTraversal():
        ap.AddInputData(it.GetCurrentDataObject()); it.GoToNextItem()
    ap.Update(); m = ap.GetOutput()
    pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, 0); pl.SetNormal(0, 0, 1)
    ct = vtk.vtkCutter(); ct.SetInputData(m); ct.SetCutFunction(pl)
    tr = vtk.vtkTriangleFilter(); tr.SetInputConnection(ct.GetOutputPort()); tr.Update(); pd = tr.GetOutput()
    xy = n(pd.GetPoints().GetData())[:, :2] * 1e3
    tri = n(pd.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    f = {pd.GetPointData().GetArrayName(i): n(pd.GetPointData().GetArray(i)) for i in range(pd.GetPointData().GetNumberOfArrays())}
    return t, xy, tri, f, m


tR, xyR, triR, fR, mR = read(rc, "internalMesh")
tC, xyC, triC, fC, mC = read(cc, "/air/internalMesh")
_, xyS, triS, fS, _ = read(cc, "/plate/internalMesh")


def mach(f):
    if "Ma" in f: return f["Ma"]
    U = np.linalg.norm(f["U"], axis=1); return U / np.sqrt(1.4 * 287.0 * f["T"])


def xo2(f):
    Y = f["O2"]; return (Y / 32.0) / (Y / 32.0 + (1 - Y) / 28.01) * 100


def solids(ax, plate=True):
    kw = dict(facecolor="white", edgecolor="#333", lw=1.0, zorder=5)
    H, LN, RN, RIN, RNOZ, RK, T, ROUT = P["H"], P["LN"], P["RN"], P["RIN"], P["RNOZ"], P["RK"], P["T"], P["ROUT"]
    for s in (1, -1):
        ax.add_patch(Polygon([[s * RN, H], [s * RNOZ, H], [s * RNOZ, H + LN], [s * RIN, H + LN]], **kw))
    if plate:
        ax.add_patch(Polygon([[-RK, 0], [-ROUT, 0], [-ROUT, -T], [-RK, -T]], **kw))


XL, YL = (-4.0, 4.0), (-4.0, 3.8)
fig, axes = plt.subplots(1, 3, figsize=(23, 6.9))
for ax, (lab, fr, fc, lo, hi, cmap) in zip(axes, (("マッハ数", mach(fR), mach(fC), 0, 2.5, "jet"),
                                                   ("温度 [K]", fR["T"], fC["T"], 150, 1800, "inferno"),
                                                   ("O₂ 濃度 [mol%]（周りの空気の巻き込み）", xo2(fR), xo2(fC), 0, 21, "RdYlBu"))):
    a = ax.tripcolor(mtri.Triangulation(-xyR[:, 0], xyR[:, 1], triR), fr, shading="gouraud", cmap=cmap, vmin=lo, vmax=hi)
    ax.tripcolor(mtri.Triangulation(xyC[:, 0], xyC[:, 1], triC), fc, shading="gouraud", cmap=cmap, vmin=lo, vmax=hi)
    if lab.startswith("温度"):
        ax.tripcolor(mtri.Triangulation(xyS[:, 0], xyS[:, 1], triS), fS["T"], shading="gouraud", cmap=cmap, vmin=lo, vmax=hi)
    solids(ax)
    ax.axvline(0, color="w", lw=1.2, ls="--", zorder=6)
    ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_aspect("equal"); ax.set_xlabel("半径 r [mm]（左：reactingFoam／右：CHT）")
    ax.set_title(lab); fig.colorbar(a, ax=ax, fraction=0.045, pad=0.02)
axes[0].set_ylabel("高さ y [mm]（0：板の上面）")
fig.suptitle("左半分：reactingFoam（板は 300 K 一定の壁）　右半分：CHT（OpenFOAM 13 foamMultiRun、板＝鋼を同時に解く、切断面 1,800 K）", fontsize=18)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(out + "_fields.png", dpi=90, facecolor="white", bbox_inches="tight")
print(out + "_fields.png")

# centre-line profiles (r = 0)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 6.8))
for xy, f, lab, st in ((xyR, fR, "reactingFoam", "-"), (xyC, fC, "CHT（foamMultiRun）", "--")):
    k = xy[:, 0] < 1e-3; o = np.argsort(xy[k, 1]); y = xy[k, 1][o]
    a1.plot(mach(f)[k][o], y, st, lw=2.8, label=lab); a2.plot(f["p"][k][o] / 1e5, y, st, lw=2.8, label=lab)
for a in (a1, a2):
    for yy, t in ((P["H"], "ノズル出口"), (0, "板の上面"), (-P["T"], "板の下面")):
        a.axhline(yy, color="#888", lw=1, ls=":"); a.text(1.01, yy, t, transform=a.get_yaxis_transform(), fontsize=14, va="center")
    a.set_ylabel("y [mm]"); a.legend(frameon=False, loc="upper right"); a.grid(alpha=0.3); a.set_ylim(-P["T"] - 4, P["H"] + P["LN"])
a1.set_xlabel("中心軸上のマッハ数"); a2.set_xlabel("中心軸上の静圧 [bar]")
fig.tight_layout(); fig.savefig(out + "_axis.png", dpi=90, facecolor="white", bbox_inches="tight"); print(out + "_axis.png")


def last(path, col=1):
    fs = sorted(glob.glob(path)); d = np.loadtxt(fs[-1], comments="#", usecols=(0, col)); return d[-1, 1]


mR_ = -last(os.path.join(rc, "postProcessing/massFlow/*/surfaceFieldValue.dat")) * 360 / 5 * 1e3
mC_ = -last(os.path.join(cc, "postProcessing/air/massFlow/*/surfaceFieldValue.dat")) * 360 / 5 * 1e3
hf = sorted(glob.glob(os.path.join(cc, "postProcessing/air/heatToGas/*/wallHeatFlux.dat")))
q = np.loadtxt(hf[-1], comments="#", usecols=(4,)) if hf else np.array([np.nan])   # integral column (W, wedge)
print(f"| ノズル流量 [g/s] | {mR_:.2f} | {mC_:.2f} |")
print(f"| 中心軸上の最大マッハ数 | {mach(fR)[xyR[:,0]<1e-3].max():.2f} | {mach(fC)[xyC[:,0]<1e-3].max():.2f} |")
print(f"| 板から気体へ移る熱 [W]（360°） | — | {np.atleast_1d(q)[-1] * 360 / 5:.1f} |")
print(f"| 板の最高温度 / 平均 [K] | 300（固定） | {fS['T'].max():.0f} / {fS['T'].mean():.0f} |")
for xy, f, lab in ((xyR, fR, "reactingFoam"), (xyC, fC, "CHT")):
    k = (xy[:, 0] < 0.15) & (np.abs(xy[:, 1] + P["T"]) < 0.05)
    print(f"| 板の下面（溝の出口）の O2 濃度 [mol%] {lab} | {xo2(f)[k].mean():.2f} | 溝の中の最高温度 {f['T'][(xy[:,1]<0)&(xy[:,1]>-P['T'])].max():.0f} K |")
print("times", tR, tC)
