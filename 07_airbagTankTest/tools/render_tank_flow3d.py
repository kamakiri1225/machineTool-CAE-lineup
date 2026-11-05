#!/usr/bin/env python3
"""Oblique 3-D GIF of the tank test: wedge results revolved 270 deg. Cut faces (2) and the tank wall show
 left: combustion-gas fraction, right: temperature.  usage: render_tank_flow3d.py <case>"""
import io, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17})
case = os.path.abspath(sys.argv[1]); OUT = os.path.join(case, "images")
rd = vtk.vtkPOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam"))
rd.SetCaseType(0 if os.path.isdir(os.path.join(case, "processor0")) else 1); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("patch/front", 1); rd.SetPatchArrayStatus("patch/walls", 1)
rd.CreateCellToPointOff()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples()) if tv.GetValue(i) > 0]


def gather(name):
    it = rd.GetOutput().NewIterator(); it.InitTraversal(); ap = vtk.vtkAppendPolyData()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            g = vtk.vtkGeometryFilter(); g.SetInputData(it.GetCurrentDataObject()); g.Update(); ap.AddInputData(g.GetOutput())
        it.GoToNextItem()
    ap.Update(); return ap.GetOutput()


def addfields(pd):
    cd = pd.GetCellData(); O2 = vtk_to_numpy(cd.GetArray("O2")); T = vtk_to_numpy(cd.GetArray("T"))
    a = numpy_to_vtk(np.clip(1 - O2 / 0.233, 0, 1), deep=True); a.SetName("gas"); cd.AddArray(a)
    b = numpy_to_vtk(T - 273.15, deep=True); b.SetName("Tc"); cd.AddArray(b)
    p = vtk_to_numpy(pd.GetPoints().GetData()) * 1e3
    return pd, p


def flat(pd):
    pd, p = addfields(pd)
    q = np.c_[np.hypot(p[:, 0], p[:, 2]), p[:, 1], np.zeros(len(p))]
    pd.GetPoints().SetData(numpy_to_vtk(q, deep=True)); return pd


def mm(pd):
    pd, p = addfields(pd); pd.GetPoints().SetData(numpy_to_vtk(p, deep=True)); return pd


def rotated(pd, ang):
    tr = vtk.vtkTransform(); tr.RotateY(ang)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(tr); f.Update(); return f.GetOutput()


def luts():
    l1 = vtk.vtkLookupTable(); l1.SetTableRange(0, 1)
    cm = plt.get_cmap("inferno"); l1.SetNumberOfTableValues(256)
    for i in range(256): l1.SetTableValue(i, *cm(i / 255)[:3], 1)
    l1.Build()
    l2 = vtk.vtkLookupTable(); l2.SetHueRange(0.667, 0.0); l2.SetTableRange(20, 730); l2.Build()
    return {"gas": l1, "Tc": l2}


LUT = luts()


def render(sec, wall, field):
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
    ap = vtk.vtkAppendPolyData()
    for ang in (45, 135):
        ap.AddInputData(rotated(sec, ang))
    for k in range(54):
        ap.AddInputData(rotated(wall, 135 + 2.5 + 5 * k))
    ap.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ap.GetOutputPort()); m.SetScalarModeToUseCellFieldData()
    m.SelectColorArray(field); m.SetLookupTable(LUT[field]); m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetAmbient(0.55); a.GetProperty().SetDiffuse(0.55); ren.AddActor(a)
    cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, 300, 0); cam.SetPosition(700, 900, -1550); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
    ren.ResetCameraClippingRange()
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(800, 900); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    return vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(900, 800, -1)[::-1]


frames = []
for t in times:
    rd.UpdateTimeStep(t); rd.Update()
    sec = flat(gather("front")); wall = mm(gather("walls"))
    fig = plt.figure(figsize=(17, 10))
    for k, (field, lab, cm, lo, hi) in enumerate((("gas", "燃焼ガスの割合（0：空気、1：燃焼ガス）", "inferno", 0, 1),
                                                  ("Tc", "温度 [℃]", "jet", 20, 730))):
        ax = fig.add_axes([0.0 + 0.5 * k, 0.02, 0.42, 0.88]); ax.imshow(render(sec, wall, field)); ax.set_axis_off()
        ax.set_title(lab, fontsize=19)
        cax = fig.add_axes([0.42 + 0.5 * k, 0.2, 0.012, 0.55])
        fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(lo, hi)), cax=cax)
    fig.suptitle(f"t = {t*1e3:5.2f} ms　60 L タンクを斜め上から（中心軸のまわりに 270° 回して表示、手前 90° は切り欠き）", fontsize=20, x=0.02, ha="left")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=75, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "tank_flow3d_latest.png"))
pal = frames[len(frames) // 2].quantize(colors=220, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
gp = os.path.join(OUT, "tank_flow3d.gif")
q[0].save(gp, save_all=True, append_images=q[1:], duration=[200] * (len(q) - 1) + [1500], loop=0, optimize=True)
print(gp, len(q), os.path.getsize(gp) / 1e6)
