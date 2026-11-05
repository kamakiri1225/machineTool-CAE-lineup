#!/usr/bin/env python3
"""GIF of the gas flow in the tank test (decomposed case is read directly):
 panel 1: combustion-gas fraction (1 - air) + velocity arrows, mirrored section
 panel 2: temperature, mirrored section
 panel 3: pressure history (CFD vs 0-D) up to the current time
usage: render_tank_flow.py <case> <zeroD_dir>"""
import glob, io, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np, vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "axes.titlesize": 18, "axes.labelsize": 16})
case, zd = os.path.abspath(sys.argv[1]), sys.argv[2]
OUT = os.path.join(case, "images"); os.makedirs(OUT, exist_ok=True)
rd = vtk.vtkPOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam"))
rd.SetCaseType(0 if os.path.isdir(os.path.join(case, "processor0")) else 1); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOff()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples()) if tv.GetValue(i) > 0]
z = np.loadtxt(os.path.join(zd, "zero_d_result.csv"), delimiter=",", skiprows=1)
fs = sorted(glob.glob(os.path.join(case, "postProcessing/tankAvg/*/volFieldValue.dat")), key=lambda f: float(f.split("/")[-2]))
tk = np.vstack([np.loadtxt(f, comments="#") for f in fs])
P0 = 101325.0


def block(mb):
    it = mb.NewIterator(); it.InitTraversal(); out = []
    while not it.IsDoneWithTraversal():
        out.append(it.GetCurrentDataObject()); it.GoToNextItem()
    ap = vtk.vtkAppendFilter()
    for o in out:
        ap.AddInputData(o)
    ap.Update(); return ap.GetOutput()


frames = []
for t in times:
    rd.UpdateTimeStep(t); rd.Update(); m = block(rd.GetOutput())
    cc = vtk.vtkCellCenters(); cc.SetInputData(m); cc.Update(); c = vtk_to_numpy(cc.GetOutput().GetPoints().GetData()) * 1e3
    cd = m.GetCellData()
    g = lambda n: vtk_to_numpy(cd.GetArray(n))
    r, y = np.hypot(c[:, 0], c[:, 2]), c[:, 1]
    gas = 1 - (g("O2") / 0.233)                         # combustion-gas fraction from the air tracer O2
    T = g("T"); U = g("U"); ur, uy = U[:, 0], U[:, 1]
    R = np.r_[r, -r]; Y = np.r_[y, y]; tri = mtri.Triangulation(R, Y)
    xm = R[tri.triangles].mean(1); ym = Y[tri.triangles].mean(1)
    tri.set_mask((ym < 0) & (np.abs(xm) > 16) | (ym < -5) & (np.abs(xm) > 15) | (ym < 0) & (ym > -5) & (np.abs(xm) > 5.6))
    fig = plt.figure(figsize=(20, 9.5))
    a1 = fig.add_axes([0.03, 0.08, 0.27, 0.80]); a2 = fig.add_axes([0.34, 0.08, 0.27, 0.80]); a3 = fig.add_axes([0.70, 0.14, 0.28, 0.70])
    for a, f, cm, lo, hi, lab in ((a1, np.r_[gas, gas], "inferno", 0, 1, "燃焼ガスの割合（0：空気、1：燃焼ガス）"),
                                  (a2, np.r_[T, T] - 273.15, "jet", 20, 730, "温度 [℃]")):
        tc = a.tripcolor(tri, f, shading="gouraud", cmap=cm, vmin=lo, vmax=hi)
        a.set_xlim(-175, 175); a.set_ylim(-80, 665); a.set_aspect("equal"); a.set_facecolor("#888")
        a.set_xlabel("r [mm]"); a.set_title(lab, fontsize=17)
        cb = fig.colorbar(tc, ax=a, fraction=0.05, pad=0.02)
    a1.set_ylabel("y [mm]")
    # velocity arrows on a regular grid
    gx, gy = np.meshgrid(np.linspace(-160, 160, 15), np.linspace(30, 640, 22))
    iu = mtri.LinearTriInterpolator(tri, np.r_[ur, -ur])(gx, gy); iv = mtri.LinearTriInterpolator(tri, np.r_[uy, uy])(gx, gy)
    sp = np.hypot(iu, iv); ok = ~np.ma.getmaskarray(sp)
    a1.quiver(gx[ok], gy[ok], (iu / (sp + 1e-9))[ok], (iv / (sp + 1e-9))[ok], color="#7FDBFF", scale=28, width=0.004)
    a1.text(-170, 630, "水色の矢印：流れの向き", color="white", fontsize=15)
    a3.plot(z[:, 0] * 1e3, (z[:, 3] - P0) / 1e3, "k--", lw=2, label="0 次元モデル（CFD と同じ比熱）")
    k = tk[:, 0] <= t
    a3.plot(tk[k, 0] * 1e3, (tk[k, 1] - P0) / 1e3, color="#1565C0", lw=3, label="CFD（タンクの平均）")
    a3.axvline(t * 1e3, color="#C62828", lw=1.5)
    a3.set_xlim(0, 50); a3.set_ylim(-5, 180); a3.set_xlabel("時間 [ms]"); a3.set_ylabel("タンクのゲージ圧 [kPa]")
    a3.set_title("タンクの圧力"); a3.legend(frameon=False, loc="lower right"); a3.grid(alpha=0.3)
    fig.suptitle(f"t = {t*1e3:5.2f} ms　インフレーターの燃焼ガスがタンクに噴き出して広がる（断面を左右に折り返して表示）", fontsize=20, x=0.02, ha="left")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=80, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "tank_flow_latest.png"))
pal = frames[-1].quantize(colors=200, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
gp = os.path.join(OUT, "tank_flow.gif")
q[0].save(gp, save_all=True, append_images=q[1:], duration=[300] * (len(q) - 1) + [1500], loop=0, optimize=True)
print(gp, len(q))
