#!/usr/bin/env python3
"""Oblique 3-D view of an axisymmetric assist-gas result: nozzle and plate revolved 270 deg (grey, cut-away),
Mach number on the two cut planes, and the sonic surface (Ma = 1) revolved as a translucent surface.
usage: render_axisym_result3d.py <case> <out.png>"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
P = json.load(open(os.path.join(case, "caseParams.json")))
RK, RN, RIN, RNOZ, H, LN, T = (P[k] for k in ("RK", "RN", "RIN", "RNOZ", "H", "LN", "T"))
RPL = 4.0
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.SetPatchArrayStatus("patch/front", 1)
rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
mb = rd.GetOutput()


def blk(name):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()


S = vtk.vtkTransform(); S.Scale(1e3, 1e3, 1e3)


def tf(pd, tr):
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(tr); f.Update(); return f.GetOutput()


def geom(ds):
    g = vtk.vtkGeometryFilter(); g.SetInputData(ds); g.Update(); return g.GetOutput()


front = tf(geom(blk("front")), S)                 # wedge face at -2.5 deg
inner = blk("internalMesh")
lut = vtk.vtkLookupTable(); lut.SetHueRange(0.667, 0.0); lut.SetTableRange(0, 2.5); lut.Build()
ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)


def flat(pd):
    """project the wedge face onto its mid plane (z = 0) so it sits exactly on the cut plane"""
    pts = vtk_to_numpy(pd.GetPoints().GetData()).copy(); r = np.hypot(pts[:, 0], pts[:, 2])
    pts[:, 0] = r; pts[:, 2] = 0.0
    q = vtk.vtkPolyData(); q.DeepCopy(pd); q.GetPoints().SetData(numpy_to_vtk(pts, deep=True)); return q


sec = flat(front)
box = vtk.vtkBox(); box.SetBounds(-0.01, RPL, -T - 2.5, H + LN, -1, 1)
cl = vtk.vtkClipPolyData(); cl.SetInputData(sec); cl.SetClipFunction(box); cl.InsideOutOn(); cl.Update(); sec = cl.GetOutput()
for ang in (45, 135):                             # the two faces of the 90-deg cut-away (solids cover 135..405 deg)
    tr = vtk.vtkTransform(); tr.RotateY(ang)
    m = vtk.vtkPolyDataMapper(); m.SetInputData(tf(sec, tr)); m.SetScalarModeToUseCellFieldData(); m.SelectColorArray("Ma")
    m.SetLookupTable(lut); m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetAmbient(0.6); a.GetProperty().SetDiffuse(0.5); ren.AddActor(a)


def solid(profile):
    pts = vtk.vtkPoints(); cells = vtk.vtkCellArray(); cells.InsertNextCell(len(profile))
    for i, (r, y) in enumerate(profile):
        pts.InsertNextPoint(r, 0.0, y); cells.InsertCellPoint(i)
    pd = vtk.vtkPolyData(); pd.SetPoints(pts); pd.SetPolys(cells)
    ex = vtk.vtkRotationalExtrusionFilter(); ex.SetInputData(pd); ex.SetAngle(270); ex.SetResolution(108); ex.CappingOn(); ex.Update()
    tr = vtk.vtkTransform(); tr.RotateY(135); tr.RotateX(-90)
    nm = vtk.vtkPolyDataNormals(); nm.SetInputData(tf(ex.GetOutput(), tr)); nm.SetFeatureAngle(30); nm.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(nm.GetOutputPort()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.72, 0.73, 0.76); a.GetProperty().SetAmbient(0.3)
    ren.AddActor(a)


solid([(RN, H), (RNOZ, H), (RNOZ, H + LN), (RIN, H + LN)])
solid([(RK, -T), (RPL, -T), (RPL, 0), (RK, 0)])
# sonic surface Ma = 1 revolved
inner.GetPointData().SetActiveScalars("Ma")
ct = vtk.vtkContourFilter(); ct.SetInputData(inner); ct.SetValue(0, 1.0); ct.Update()
iso = flat(tf(geom(ct.GetOutput()), S))
c2 = vtk.vtkClipPolyData(); c2.SetInputData(iso); c2.SetClipFunction(box); c2.InsideOutOn(); c2.Update(); iso = c2.GetOutput()
ap = vtk.vtkAppendPolyData()
for k in range(54):
    tr = vtk.vtkTransform(); tr.RotateY(135 + 2.5 + 5 * k)
    ap.AddInputData(tf(iso, tr))
ap.Update()
m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ap.GetOutputPort()); m.ScalarVisibilityOff()
pass  # sonic surface hidden behind the sections; not drawn
cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, -0.5, 0); cam.SetPosition(7.5, 7.0, -14.5); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1300, 1100
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]
fig = plt.figure(figsize=(15, 11)); ax = fig.add_axes([0, 0, 0.86, 0.93]); ax.imshow(img); ax.set_axis_off()
cax = fig.add_axes([0.87, 0.2, 0.018, 0.55])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(0, 2.5)), cax=cax); cb.set_label("マッハ数")
gas = "O₂" if P.get("YO2_in", 1) > 0.5 else "空気"
fig.suptitle(f"{gas} {P['p0_gauge']/1e6:g} MPa(G)：ノズル〜切断溝のマッハ数（斜め上から見た断面）", fontsize=21, x=0.02, ha="left")
ax.text(20, Hh - 30, "灰色：ノズルと鋼板（中心軸のまわりに 270° 回して表示、手前 90° は切り欠き）", fontsize=17,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", alpha=0.95))
fig.savefig(out, dpi=90, facecolor="white", bbox_inches="tight"); print(out)
