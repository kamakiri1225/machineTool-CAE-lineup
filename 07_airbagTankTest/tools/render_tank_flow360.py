#!/usr/bin/env python3
"""360-deg oblique GIF: the combustion-gas cloud (iso-surface of gas fraction = ISO, coloured by temperature)
revolved 360 deg inside a translucent 60 L tank. usage: render_tank_flow360.py <case> [ISO]"""
import io, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17})
case = os.path.abspath(sys.argv[1]); ISO = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
OUT = os.path.join(case, "images")
rd = vtk.vtkPOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam"))
rd.SetCaseType(0 if os.path.isdir(os.path.join(case, "processor0")) else 1); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.CreateCellToPointOn()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples()) if tv.GetValue(i) > 0]
lut = vtk.vtkLookupTable(); lut.SetHueRange(0.667, 0.0); lut.SetTableRange(20, 730); lut.Build()


def merged():
    it = rd.GetOutput().NewIterator(); it.InitTraversal(); ap = vtk.vtkAppendFilter()
    while not it.IsDoneWithTraversal():
        ap.AddInputData(it.GetCurrentDataObject()); it.GoToNextItem()
    ap.Update(); return ap.GetOutput()


def tank_actors(ren):
    cyl = vtk.vtkCylinderSource(); cyl.SetRadius(170); cyl.SetHeight(660); cyl.SetCenter(0, 330, 0); cyl.SetResolution(90); cyl.CappingOff()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(cyl.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.6, 0.72, 0.9); a.GetProperty().SetOpacity(0.12); ren.AddActor(a)
    for y in (0, 660):
        c = vtk.vtkRegularPolygonSource(); c.SetNumberOfSides(90); c.SetRadius(170); c.SetCenter(0, y, 0); c.SetNormal(0, 1, 0); c.GeneratePolygonOff()
        m = vtk.vtkPolyDataMapper(); m.SetInputConnection(c.GetOutputPort())
        a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.3, 0.4, 0.6); a.GetProperty().SetLineWidth(2); ren.AddActor(a)
    inf = vtk.vtkCylinderSource(); inf.SetRadius(15); inf.SetHeight(70); inf.SetCenter(0, -40, 0); inf.SetResolution(40)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(inf.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.95, 0.55, 0.2); ren.AddActor(a)


frames = []
for t in times:
    rd.UpdateTimeStep(t); rd.Update(); m = merged()
    pd = m.GetPointData(); O2 = vtk_to_numpy(pd.GetArray("O2")); T = vtk_to_numpy(pd.GetArray("T"))
    g = numpy_to_vtk(np.clip(1 - O2 / 0.233, 0, 1), deep=True); g.SetName("gas"); pd.AddArray(g)
    tc = numpy_to_vtk(T - 273.15, deep=True); tc.SetName("Tc"); pd.AddArray(tc)
    pd.SetActiveScalars("gas")
    ct = vtk.vtkContourFilter(); ct.SetInputData(m); ct.SetValue(0, ISO); ct.Update()
    sc = vtk.vtkTransform(); sc.Scale(1e3, 1e3, 1e3)
    f0 = vtk.vtkTransformPolyDataFilter(); f0.SetInputConnection(ct.GetOutputPort()); f0.SetTransform(sc); f0.Update()
    iso = f0.GetOutput()
    ap = vtk.vtkAppendPolyData()
    for k in range(72):
        tr = vtk.vtkTransform(); tr.RotateY(5 * k)
        f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(iso); f.SetTransform(tr); f.Update(); ap.AddInputData(f.GetOutput())
    ap.Update()
    cl = vtk.vtkCleanPolyData(); cl.SetInputConnection(ap.GetOutputPort()); cl.Update()
    nm = vtk.vtkPolyDataNormals(); nm.SetInputConnection(cl.GetOutputPort()); nm.SetFeatureAngle(60); nm.Update()
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1); tank_actors(ren)
    mp = vtk.vtkPolyDataMapper(); mp.SetInputConnection(nm.GetOutputPort()); mp.SetScalarModeToUsePointFieldData(); mp.SelectColorArray("Tc")
    mp.SetLookupTable(lut); mp.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(mp); a.GetProperty().SetSpecular(0.3); a.GetProperty().SetAmbient(0.25); ren.AddActor(a)
    cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, 300, 0); cam.SetPosition(900, 1050, 1300); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
    ren.ResetCameraClippingRange()
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(900, 1000); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(1000, 900, -1)[::-1]
    fig = plt.figure(figsize=(11, 11)); ax = fig.add_axes([0, 0, 0.86, 0.92]); ax.imshow(img); ax.set_axis_off()
    cax = fig.add_axes([0.86, 0.25, 0.025, 0.5])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(20, 730)), cax=cax); cb.set_label("温度 [℃]")
    fig.suptitle(f"t = {t*1e3:5.2f} ms　燃焼ガスの雲（燃焼ガスの割合 {ISO:g} の面、色は温度）を 360° で表示", fontsize=19, x=0.02, ha="left")
    ax.text(10, 990, "薄い青の円筒：60 L タンク　オレンジ：インフレーター", fontsize=16)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=75, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "tank_flow360_latest.png"))
pal = frames[len(frames) // 3].quantize(colors=220, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
gp = os.path.join(OUT, "tank_flow360.gif")
q[0].save(gp, save_all=True, append_images=q[1:], duration=[200] * (len(q) - 1) + [1500], loop=0, optimize=True)
print(gp, len(q), os.path.getsize(gp) / 1e6)
