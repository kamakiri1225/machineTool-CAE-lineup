#!/usr/bin/env python3
"""GIF of the spindle heat-up: mirrored axisymmetric section (temperature) + probe history + shaft growth.
usage: render_spindle.py <caseDir>"""
import io, json, os, re, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from matplotlib.collections import PolyCollection
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17,
                     "axes.titlesize": 19, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15})
case = os.path.abspath(sys.argv[1]); P = json.load(open(os.path.join(case, "caseParams.json")))
OUT = os.path.join(case, "images"); os.makedirs(OUT, exist_ok=True)
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOff()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())]
pr = np.loadtxt(os.path.join(case, "postProcessing/Tprobe/0/T"), comments="#")
ALPHA = 11.5e-6
frames = []
for t in times:
    rd.UpdateTimeStep(t); rd.Update(); m = rd.GetOutput().GetBlock(0)
    T = vtk_to_numpy(m.GetCellData().GetArray("T")) - 273.15 if m.GetCellData().GetArray("T") else None
    if T is None:
        continue
    pts = vtk_to_numpy(m.GetPoints().GetData()) * 1e3
    polys, vals, cen = [], [], []
    for i in range(m.GetNumberOfCells()):
        c = m.GetCell(i); ids = [c.GetPointIds().GetId(k) for k in range(c.GetNumberOfPoints())]
        q = pts[ids]; q = q[q[:, 2] <= 1e-9]            # front face (z <= 0)
        r = np.hypot(q[:, 0], q[:, 2]); y = q[:, 1]
        poly = np.c_[r, y]
        cc = poly.mean(0); ang = np.arctan2(poly[:, 1] - cc[1], poly[:, 0] - cc[0]); poly = poly[np.argsort(ang)]
        polys += [poly, poly * [-1, 1]]; vals += [T[i], T[i]]; cen.append(cc)
    cen = np.array(cen)
    # shaft growth: centre-line cells (r < 12.5 mm), integrate alpha*(T-T0) dy
    sh = cen[:, 0] < 12.5
    yy = cen[sh, 1]; tt = T[sh]; o = np.argsort(yy); dy = 300.0 / max(1, len(yy))
    dL = ALPHA * np.sum(tt - (P["T0"] - 273.15)) * dy * 1e3       # um
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 8.2), gridspec_kw={"width_ratios": [1.0, 1.15]})
    pc = PolyCollection(polys, array=np.array(vals), cmap="turbo", edgecolors="face", linewidths=0.3)
    pc.set_clim(20, 36); a1.add_collection(pc)
    a1.set_xlim(-100, 100); a1.set_ylim(-10, 310); a1.set_aspect("equal")
    a1.set_xlabel("r [mm]"); a1.set_ylabel("軸方向 y [mm]")
    a1.set_title(f"t = {t/60:4.1f} 分　温度 [℃]")
    cb = fig.colorbar(pc, ax=a1, fraction=0.05, pad=0.02); cb.set_label("温度 [℃]")
    kw = dict(fontsize=15, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.9))
    a1.text(-97, 50, f"前軸受\n{P['Qf']:.0f} W", **kw); a1.text(-97, 240, f"後軸受\n{P['Qr']:.0f} W", **kw)
    a1.text(-20, 285, "主軸", **kw); a1.text(-97, 140, "ハウジング\n（外周を水で冷却）", **kw)
    k = pr[:, 0] <= t
    for col, lab, c in ((1, "主軸の中心（前端）", "#C62828"), (2, "前軸受", "#EF6C00"), (3, "後軸受", "#2E7D32"), (4, "ハウジング外周", "#1565C0")):
        a2.plot(pr[k, 0] / 60, pr[k, col] - 273.15, lw=3, color=c, label=lab)
        a2.plot(pr[:, 0] / 60, pr[:, col] - 273.15, lw=1, color=c, alpha=0.25)
    a2.set_xlim(0, P["tEnd"] / 60); a2.set_ylim(19, 37)
    a2.set_xlabel("時間 [分]"); a2.set_ylabel("温度 [℃]"); a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, frameon=False)
    a2.set_title("各部の温度の上がり方")
    a2.text(0.03, 0.93, f"主軸の熱膨張による伸び（概算）：{dL:5.1f} µm", transform=a2.transAxes, fontsize=18,
            bbox=dict(boxstyle="round,pad=0.3", fc="#FFF3E0", ec="#EF6C00"))
    fig.suptitle("軸受の発熱で主軸が温まり、ハウジングの冷却ジャケットで冷える（軸対称・固体の熱伝導）", fontsize=19)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=90, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "spindle_final.png"))
pal = frames[-1].quantize(colors=160, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
g = os.path.join(OUT, "spindle_heatup.gif")
q[0].save(g, save_all=True, append_images=q[1:], duration=[180] * (len(q) - 1) + [1800], loop=0, optimize=True)
print(g, len(q), os.path.getsize(g) / 1e6, "MB, dL", dL)
