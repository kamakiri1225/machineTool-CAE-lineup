#!/usr/bin/env python3
"""Figure explaining the axisymmetric (wedge) tank model:
 left : the real 60 L tank as a translucent cylinder, the 5-deg wedge that is actually computed (with mesh), the axis
 right: the wedge seen from the front = the r-y section, with the regions labelled.
usage: render_tank_axisym.py <case> <out.png>"""
import os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from matplotlib.collections import PolyCollection
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.SetPatchArrayStatus("patch/front", 1); rd.Update()
mb = rd.GetOutput()


def get(name):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()


def mm(ds, sx=1.0):
    g = vtk.vtkGeometryFilter(); g.SetInputData(ds); g.Update()
    tr = vtk.vtkTransform(); tr.Scale(1e3, 1e3, 1e3 * sx)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(g.GetOutput()); f.SetTransform(tr); f.Update(); return f.GetOutput()


ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
# full tank ghost (360 deg)
for prof, col, op in (([(0, 0), (170, 0), (170, 660), (0, 660)], (0.6, 0.75, 0.95), 0.12),
                      ([(0, -75), (15, -75), (15, -5), (5.5, -5), (5.5, 0), (0, 0)], (0.95, 0.5, 0.2), 0.25)):
    pts = vtk.vtkPoints(); cells = vtk.vtkCellArray(); cells.InsertNextCell(len(prof))
    for i, (r, y) in enumerate(prof):
        pts.InsertNextPoint(r, 0.0, y); cells.InsertCellPoint(i)
    pd = vtk.vtkPolyData(); pd.SetPoints(pts); pd.SetPolys(cells)
    ex = vtk.vtkRotationalExtrusionFilter(); ex.SetInputData(pd); ex.SetAngle(360); ex.SetResolution(120); ex.CappingOff(); ex.Update()
    tr = vtk.vtkTransform(); tr.RotateX(-90)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(ex.GetOutput()); f.SetTransform(tr); f.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputData(f.GetOutput()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*col); a.GetProperty().SetOpacity(op); ren.AddActor(a)
    fe = vtk.vtkFeatureEdges(); fe.SetInputData(f.GetOutput()); fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff()
    em = vtk.vtkPolyDataMapper(); em.SetInputConnection(fe.GetOutputPort()); em.ScalarVisibilityOff()
    ea = vtk.vtkActor(); ea.SetMapper(em); ea.GetProperty().SetColor(0.3, 0.4, 0.6); ren.AddActor(ea)
# computed wedge (real 5 deg), solid + mesh lines
w = mm(get("internalMesh"))
m = vtk.vtkPolyDataMapper(); m.SetInputData(w); m.ScalarVisibilityOff()
a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.15, 0.45, 0.85); a.GetProperty().SetEdgeVisibility(1)
a.GetProperty().SetEdgeColor(0.05, 0.15, 0.35); a.GetProperty().SetLineWidth(0.4); ren.AddActor(a)
ln = vtk.vtkLineSource(); ln.SetPoint1(0, -120, 0); ln.SetPoint2(0, 720, 0)
m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ln.GetOutputPort()); a = vtk.vtkActor(); a.SetMapper(m)
a.GetProperty().SetColor(0.85, 0.1, 0.1); a.GetProperty().SetLineWidth(3); ren.AddActor(a)
cam = ren.GetActiveCamera(); cam.SetFocalPoint(60, 300, 0); cam.SetPosition(950, 900, 1300); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1100, 1200
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(x, y, z):
    c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(x, y, z); d = c.GetComputedDoubleDisplayValue(ren)
    return d[0], Hh - d[1]


fig = plt.figure(figsize=(21, 11.5))
ax = fig.add_axes([0, 0, 0.5, 0.93]); ax.imshow(img); ax.set_axis_off()
kw = dict(fontsize=18, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", alpha=0.95))
ar = dict(arrowstyle="-|>", lw=1.6)
for pt, lab, xy in (((120, 400, 5), "実際に計算する範囲：\n中心軸のまわりの 5° のくさび（青、メッシュ付き）", (0.48, 0.18)),
                    ((-120, 500, -100), "本物のタンク（60 L の円筒、薄い青）\n＝くさびを 360° 回したもの", (0.02, 0.06)),
                    ((0, -100, 0), "中心軸（回転の軸、赤い線）", (0.42, 0.86)),
                    ((-10, -40, 10), "インフレーター", (0.02, 0.88))):
    x, y = proj(*pt)
    ax.annotate(lab, xy=(x, y), xytext=(xy[0] * W, xy[1] * Hh + 40), arrowprops=ar, **kw)
# right: section
f0 = get("front"); g = vtk.vtkGeometryFilter(); g.SetInputData(f0); g.Update(); pd = g.GetOutput()
pts = vtk_to_numpy(pd.GetPoints().GetData()) * 1e3; r = np.hypot(pts[:, 0], pts[:, 2]); y = pts[:, 1]
polys, cols = [], []
for i in range(pd.GetNumberOfCells()):
    ids = pd.GetCell(i).GetPointIds(); k = [ids.GetId(j) for j in range(ids.GetNumberOfIds())]
    q = np.c_[r[k], y[k]]; cy = q[:, 1].mean()
    polys.append(q); cols.append("#F4A261" if cy < -5 else ("#E63946" if cy < 0 else "#cfe3ff"))
a2 = fig.add_axes([0.55, 0.08, 0.20, 0.82]); a2.set_facecolor("#9a9a9a")
a2.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="#3a6ea5", linewidths=0.25))
a2.set_xlim(0, 180); a2.set_ylim(-90, 670); a2.set_aspect("equal")
a2.axvline(0, color="#d32f2f", lw=4); a2.set_xlabel("r [mm]"); a2.set_ylabel("y [mm]")
a2.set_title("くさびの断面（赤：中心軸）")
a3 = fig.add_axes([0.78, 0.08, 0.20, 0.82]); a3.set_facecolor("#9a9a9a")
a3.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="#3a6ea5", linewidths=0.4))
a3.set_xlim(0, 40); a3.set_ylim(-85, 40); a3.set_aspect("equal"); a3.axvline(0, color="#d32f2f", lw=4)
a3.set_xlabel("r [mm]"); a3.set_title("インフレーターの拡大")
kw2 = dict(fontsize=16, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.95))
a3.text(17, -45, "燃焼室\n50 cc", **kw2); a3.annotate("オリフィス\n径 11 mm", xy=(4, -2.5), xytext=(12, 15), arrowprops=dict(arrowstyle="-|>"), **kw2)
a2.text(40, 330, "タンクの中\n（空気）", **kw2); fig.suptitle("軸対称モデル：円筒のタンクを、中心軸のまわりの 5° のくさびだけ計算する（10,700 セル）", fontsize=22, x=0.02, ha="left")
fig.savefig(out, dpi=85, facecolor="white", bbox_inches="tight"); print(out)
