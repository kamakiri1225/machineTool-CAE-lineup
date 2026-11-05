#!/usr/bin/env python3
"""Oblique 3-D close-up of the inflator: chamber, orifice and the tank bottom revolved 270 deg (cut-away),
the two cut faces show the computational mesh coloured by region. usage: render_inflator_closeup3d.py <case> <out.png>"""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("patch/front", 1); rd.Update()
it = rd.GetOutput().NewIterator(); it.InitTraversal(); f = None
while not it.IsDoneWithTraversal():
    f = it.GetCurrentDataObject(); it.GoToNextItem()
g = vtk.vtkGeometryFilter(); g.SetInputData(f); g.Update(); sec = g.GetOutput()
pts = vtk_to_numpy(sec.GetPoints().GetData()) * 1e3
p2 = np.c_[np.hypot(pts[:, 0], pts[:, 2]), pts[:, 1], np.zeros(len(pts))]
sec.GetPoints().SetData(numpy_to_vtk(p2, deep=True))
cc = vtk.vtkCellCenters(); cc.SetInputData(sec); cc.Update(); cy = vtk_to_numpy(cc.GetOutput().GetPoints().GetData())[:, 1]
reg = numpy_to_vtk(np.where(cy < -5, 0, np.where(cy < 0, 1, 2)).astype(float), deep=True); reg.SetName("reg"); sec.GetCellData().AddArray(reg)
bx = vtk.vtkBox(); bx.SetBounds(-1, 45, -90, 45, -1, 1)
cl = vtk.vtkClipPolyData(); cl.SetInputData(sec); cl.SetClipFunction(bx); cl.InsideOutOn(); cl.Update(); sec = cl.GetOutput()
lut = vtk.vtkLookupTable(); lut.SetNumberOfTableValues(3); lut.SetTableRange(0, 2)
for i, c in enumerate(((0.96, 0.64, 0.38), (0.90, 0.22, 0.27), (0.81, 0.89, 1.0))):
    lut.SetTableValue(i, *c, 1.0)
lut.Build()
ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
for ang in (45, 135):
    tr = vtk.vtkTransform(); tr.RotateY(ang)
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetInputData(sec); tf.SetTransform(tr); tf.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputData(tf.GetOutput()); m.SetScalarModeToUseCellFieldData(); m.SelectColorArray("reg")
    m.SetLookupTable(lut); m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetEdgeVisibility(1); a.GetProperty().SetEdgeColor(0.2, 0.35, 0.6)
    a.GetProperty().SetLineWidth(0.6); a.GetProperty().SetAmbient(0.7); a.GetProperty().SetDiffuse(0.4); ren.AddActor(a)


def solid(profile, color, op=1.0):
    p = vtk.vtkPoints(); cells = vtk.vtkCellArray(); cells.InsertNextCell(len(profile))
    for i, (r, y) in enumerate(profile):
        p.InsertNextPoint(r, 0.0, y); cells.InsertCellPoint(i)
    pd = vtk.vtkPolyData(); pd.SetPoints(p); pd.SetPolys(cells)
    ex = vtk.vtkRotationalExtrusionFilter(); ex.SetInputData(pd); ex.SetAngle(270); ex.SetResolution(108); ex.CappingOn(); ex.Update()
    tr = vtk.vtkTransform(); tr.RotateY(135); tr.RotateX(-90)
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetInputData(ex.GetOutput()); tf.SetTransform(tr); tf.Update()
    nm = vtk.vtkPolyDataNormals(); nm.SetInputData(tf.GetOutput()); nm.SetFeatureAngle(30); nm.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(nm.GetOutputPort()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*color); a.GetProperty().SetOpacity(op); a.GetProperty().SetAmbient(0.3)
    ren.AddActor(a)


solid([(0, -85), (25, -85), (25, -5), (45, -5), (45, 0), (15, 0), (15, -75), (0, -75)], (0.62, 0.63, 0.67))   # housing / tank floor (display)
solid([(0, -75), (15, -75), (15, -5), (5.5, -5), (5.5, 0), (0, 0)], (0.95, 0.55, 0.25), 0.25)                 # gas in chamber (ghost)
cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, -20, 0); cam.SetPosition(150, 120, 210); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1300, 1100
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(r, y, ang):
    a = np.radians(ang); c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(r * np.cos(a), y, -r * np.sin(a))
    d = c.GetComputedDoubleDisplayValue(ren); return d[0], Hh - d[1]


fig = plt.figure(figsize=(14, 12)); ax = fig.add_axes([0, 0, 1, 0.94]); ax.imshow(img); ax.set_axis_off()
kw = dict(fontsize=19, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", alpha=0.95))
for (r, y, ang), lab, xy in (((8, -40, 45), "燃焼室（50 cc）\nここで燃焼ガスが湧き出す", (0.02, 0.62)),
                            ((3, -2.5, 45), "オリフィス（径 11 mm）", (0.02, 0.30)),
                            ((25, 25, 135), "タンクの中（空気）", (0.62, 0.10)),
                            ((20, -40, 300), "インフレーターの外側（表示用、計算しない）", (0.45, 0.88))):
    x, y2 = proj(r, y, ang)
    ax.annotate(lab, xy=(x, y2), xytext=(xy[0] * W, xy[1] * Hh), arrowprops=dict(arrowstyle="-|>", lw=1.6), **kw)
fig.suptitle("インフレーターのまわりを斜め上から（270° 回して表示、切り口の 2 面に計算メッシュ）", fontsize=20)
fig.savefig(out, dpi=90, facecolor="white", bbox_inches="tight"); print(out)
