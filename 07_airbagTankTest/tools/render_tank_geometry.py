#!/usr/bin/env python3
"""Oblique cut-away (270 deg) view of the inflator tank-test model with labels, built from the wedge mesh.
usage: render_tank_geometry.py <case> <out.png>"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)


def solid(profile, color, op=1.0):
    pts = vtk.vtkPoints(); cells = vtk.vtkCellArray(); cells.InsertNextCell(len(profile))
    for i, (r, y) in enumerate(profile):
        pts.InsertNextPoint(r, 0.0, y); cells.InsertCellPoint(i)
    pd = vtk.vtkPolyData(); pd.SetPoints(pts); pd.SetPolys(cells)
    ex = vtk.vtkRotationalExtrusionFilter(); ex.SetInputData(pd); ex.SetAngle(270); ex.SetResolution(108); ex.CappingOn(); ex.Update()
    tr = vtk.vtkTransform(); tr.RotateY(135); tr.RotateX(-90)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(ex.GetOutput()); f.SetTransform(tr); f.Update()
    nm = vtk.vtkPolyDataNormals(); nm.SetInputData(f.GetOutput()); nm.SetFeatureAngle(30); nm.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(nm.GetOutputPort()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*color); a.GetProperty().SetOpacity(op); a.GetProperty().SetAmbient(0.3)
    ren.AddActor(a)


solid([(0, 0), (170, 0), (170, 660), (0, 660)], (0.55, 0.75, 1.0), 0.45)           # air in the tank (cut-away)
solid([(0, -75), (15, -75), (15, -5), (0, -5)], (0.95, 0.45, 0.15))                   # combustion chamber
solid([(0, -5), (5.5, -5), (5.5, 0), (0, 0)], (0.85, 0.1, 0.1))                       # orifice
cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, 290, 0); cam.SetPosition(900, 1000, -1500); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1200, 1200
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(x, y, z):
    c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(x, y, z); d = c.GetComputedDoubleDisplayValue(ren)
    return d[0], Hh - d[1]


fig = plt.figure(figsize=(20, 11)); ax = fig.add_axes([0, 0, 0.5, 0.93]); ax.imshow(img); ax.set_axis_off()
kw = dict(fontsize=18, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", alpha=0.95))
for pt, lab, xy in (((0, 450, -120), "60 L タンク（内径 340 mm × 高さ 660 mm）\n初期：空気 101.3 kPa、20 ℃、壁は断熱", (0.40, 0.08)),
                    ((-8, -40, 0), "インフレーター（燃焼室＋オリフィス）", (0.45, 0.90))):
    x, y = proj(*pt)
    ax.annotate(lab, xy=(x, y), xytext=(xy[0] * W, xy[1] * Hh), arrowprops=dict(arrowstyle="-|>", lw=1.6), **kw)
# right: close-up section with mesh
from matplotlib.collections import PolyCollection
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("patch/front", 1); rd.Update()
it = rd.GetOutput().NewIterator(); it.InitTraversal(); pd = None
while not it.IsDoneWithTraversal():
    pd = it.GetCurrentDataObject(); it.GoToNextItem()
g = vtk.vtkGeometryFilter(); g.SetInputData(pd); g.Update(); pd = g.GetOutput()
pts = vtk_to_numpy(pd.GetPoints().GetData()) * 1e3; r = np.hypot(pts[:, 0], pts[:, 2]); y = pts[:, 1]
polys, cols = [], []
for i in range(pd.GetNumberOfCells()):
    ids = pd.GetCell(i).GetPointIds(); k = [ids.GetId(j) for j in range(ids.GetNumberOfIds())]
    q = np.c_[r[k], y[k]]; cy = q[:, 1].mean()
    c = "#F4A261" if cy < -5 else ("#E63946" if cy < 0 else "#cfe3ff")
    polys += [q, q * [-1, 1]]; cols += [c, c]
a2 = fig.add_axes([0.53, 0.07, 0.44, 0.83]); a2.set_facecolor("#9a9a9a")
a2.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="#3a6ea5", linewidths=0.3))
a2.set_xlim(-40, 40); a2.set_ylim(-85, 60); a2.set_aspect("equal")
a2.set_xlabel("r [mm]（左右に折り返し）"); a2.set_ylabel("y [mm]"); a2.set_title("インフレーターのまわりの拡大（断面とメッシュ）")
kw2 = dict(fontsize=17, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.95))
a2.annotate("燃焼室（50 cc）\nガス発生剤が燃えて\n1,000 K の燃焼ガス 25 g を出す\n（N₂ 70 %・CO₂ 20 %・H₂O 10 %）", xy=(5, -45), xytext=(17, -70), arrowprops=dict(arrowstyle="-|>"), **kw2)
a2.annotate("オリフィス（径 11 mm）", xy=(3, -2.5), xytext=(10, 15), arrowprops=dict(arrowstyle="-|>"), **kw2)
a2.text(-38, 45, "タンクの中（空気）", **kw2)
fig.suptitle("火薬式インフレーターのタンク試験モデル（軸対称）", fontsize=22, x=0.02, ha="left")
fig.savefig(out, dpi=85, facecolor="white", bbox_inches="tight"); print(out)
