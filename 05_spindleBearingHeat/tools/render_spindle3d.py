#!/usr/bin/env python3
"""3-D oblique GIF of the spindle heat-up: the 5-deg wedge result revolved 270 deg (cut-away shows the section),
coloured by temperature, + probe history panel. usage: render_spindle3d.py <caseDir>"""
import io, json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17,
                     "axes.titlesize": 19, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15})
case = os.path.abspath(sys.argv[1]); P = json.load(open(os.path.join(case, "caseParams.json")))
OUT = os.path.join(case, "images"); os.makedirs(OUT, exist_ok=True)
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOff()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())]
pr = np.loadtxt(os.path.join(case, "postProcessing/Tprobe/0/T"), comments="#")
lut = vtk.vtkLookupTable(); lut.SetHueRange(0.667, 0.0); lut.SetTableRange(20, 36); lut.Build()
frames = []
for t in times:
    rd.UpdateTimeStep(t); rd.Update(); m = rd.GetOutput().GetBlock(0)
    if m.GetCellData().GetArray("T") is None:
        continue
    g = vtk.vtkGeometryFilter(); g.SetInputData(m); g.Update(); s = g.GetOutput()
    tc = numpy_to_vtk(vtk_to_numpy(s.GetCellData().GetArray("T")) - 273.15, deep=True); tc.SetName("Tc"); s.GetCellData().AddArray(tc)
    sc = vtk.vtkTransform(); sc.Scale(1e3, 1e3, 1e3)
    f0 = vtk.vtkTransformPolyDataFilter(); f0.SetInputData(s); f0.SetTransform(sc); f0.Update(); s = f0.GetOutput()
    ap = vtk.vtkAppendPolyData()
    for k in range(54):                       # 54 x 5 deg = 270 deg, gap faces the camera
        tr = vtk.vtkTransform(); tr.RotateY(-45 + 2.5 + 5 * k)
        f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(s); f.SetTransform(tr); f.Update(); ap.AddInputData(f.GetOutput())
    ap.Update()
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
    mp = vtk.vtkPolyDataMapper(); mp.SetInputConnection(ap.GetOutputPort()); mp.SetScalarModeToUseCellFieldData()
    mp.SelectColorArray("Tc"); mp.SetLookupTable(lut); mp.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(mp); a.GetProperty().SetAmbient(0.3); a.GetProperty().SetDiffuse(0.8); ren.AddActor(a)
    cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, 140, 0); cam.SetPosition(330, 420, 520); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
    ren.ResetCameraClippingRange()
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(900, 900); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(900, 900, -1)[::-1]
    fig = plt.figure(figsize=(17, 8.4))
    ax = fig.add_axes([0.0, 0.0, 0.47, 0.92]); ax.imshow(img); ax.set_axis_off()
    kw = dict(fontsize=16, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.92))
    ax.text(15, 870, "中心軸のまわりに 270° 回して表示（手前 90° は切り欠き）", **kw)
    ax.annotate("前軸受（150 W）", xy=(350, 655), xytext=(15, 790), fontsize=16, bbox=kw["bbox"], arrowprops=dict(arrowstyle="-|>", lw=1.5))
    ax.annotate("後軸受（80 W）", xy=(350, 210), xytext=(15, 330), fontsize=16, bbox=kw["bbox"], arrowprops=dict(arrowstyle="-|>", lw=1.5))
    ax.annotate("主軸", xy=(425, 110), xytext=(560, 40), fontsize=16, bbox=kw["bbox"], arrowprops=dict(arrowstyle="-|>", lw=1.5))
    ax.annotate("ハウジング（外周を水冷）", xy=(600, 420), xytext=(560, 790), fontsize=16, bbox=kw["bbox"], arrowprops=dict(arrowstyle="-|>", lw=1.5))
    cax = fig.add_axes([0.47, 0.18, 0.012, 0.6])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(20, 36)), cax=cax); cb.set_label("温度 [℃]")
    a2 = fig.add_axes([0.58, 0.25, 0.40, 0.62])
    k = pr[:, 0] <= t
    for col, lab, c in ((1, "主軸の中心（前端）", "#C62828"), (2, "前軸受", "#EF6C00"), (3, "後軸受", "#2E7D32"), (4, "ハウジング外周", "#1565C0")):
        a2.plot(pr[:, 0] / 60, pr[:, col] - 273.15, lw=1, color=c, alpha=0.25)
        a2.plot(pr[k, 0] / 60, pr[k, col] - 273.15, lw=3, color=c, label=lab)
    a2.set_xlim(0, P["tEnd"] / 60); a2.set_ylim(19, 37); a2.set_xlabel("時間 [分]"); a2.set_ylabel("温度 [℃]")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False); a2.set_title("各部の温度の上がり方")
    fig.suptitle(f"t = {t/60:4.1f} 分　軸受の発熱（前 {P['Qf']:.0f} W・後 {P['Qr']:.0f} W）で主軸が温まり、ハウジングの水冷で冷える", fontsize=19, x=0.02, ha="left")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=85, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "spindle3d_final.png"))
pal = frames[-1].quantize(colors=200, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
g = os.path.join(OUT, "spindle3d_heatup.gif")
q[0].save(g, save_all=True, append_images=q[1:], duration=[180] * (len(q) - 1) + [1800], loop=0, optimize=True)
print(g, len(q), os.path.getsize(g) / 1e6)
