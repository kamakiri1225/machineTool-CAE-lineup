#!/usr/bin/env python3
"""Oblique 3-D view of the CHT regions (which part is fluid and which is solid), ParaView-like:
each region's 5-degree wedge section ('front' patch) is revolved 270 deg about the axis; the two cut faces show the
real computational mesh. Gas = translucent blue volume + mesh on the cut faces, plate = grey, nozzle = orange.
usage: render_cht_regions3d.py <case> <out.png> [probe_csv]"""
import csv, math, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from matplotlib.patches import Patch
from vtk.util.numpy_support import vtk_to_numpy as n, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
probes = []
if len(sys.argv) > 3 and os.path.exists(sys.argv[3]):
    for row in csv.DictReader(open(sys.argv[3])):
        k = {kk.lower(): vv for kk, vv in row.items()}
        probes.append(((k.get("probe") or k.get("name")).split("_")[0], float(k.get("x", k.get("x_m"))), float(k.get("y", k.get("y_m")))))
START, SWEEP = -45.0, 270.0            # revolve from -45 to 225 deg; the open 90 deg faces the camera

rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()


def section(region):
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus(f"/{region}/patch/front", 1); rd.UpdateTimeStep(0); rd.Update()
    it = rd.GetOutput().NewIterator(); it.InitTraversal(); ap = vtk.vtkAppendPolyData()
    while not it.IsDoneWithTraversal():
        g = vtk.vtkGeometryFilter(); g.SetInputData(it.GetCurrentDataObject()); g.Update(); ap.AddInputData(g.GetOutput()); it.GoToNextItem()
    ap.Update(); pd = ap.GetOutput()
    p = n(pd.GetPoints().GetData()) * 1e3                      # mm
    r = np.hypot(p[:, 0], p[:, 2])
    pd.GetPoints().SetData(numpy_to_vtk(np.c_[r, np.zeros_like(r), p[:, 1]], deep=True))   # (r, 0, y): axis = z
    return pd


def rotz(pd, ang):
    t = vtk.vtkTransform(); t.RotateZ(ang)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(t); f.Update(); return f.GetOutput()


def revolve(pd):
    e = vtk.vtkFeatureEdges(); e.SetInputData(pd); e.BoundaryEdgesOn(); e.FeatureEdgesOff(); e.ManifoldEdgesOff()
    e.NonManifoldEdgesOff(); e.Update()                          # outline of the section -> surface of revolution
    x = vtk.vtkRotationalExtrusionFilter(); x.SetInputData(rotz(e.GetOutput(), START)); x.SetAngle(SWEEP)
    x.SetResolution(108); x.CappingOff(); x.Update(); return x.GetOutput()


ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
COL = {"air": (0.55, 0.78, 0.97), "plate": (0.62, 0.65, 0.70), "nozzle": (0.95, 0.52, 0.18)}
for region in ("plate", "nozzle", "air"):
    sec = section(region)
    m = vtk.vtkPolyDataMapper(); m.SetInputData(revolve(sec)); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*COL[region]); a.GetProperty().SetSpecular(0.15)
    if region == "air":                                          # gas: only the outline of its extent (wireframe edges)
        oe = vtk.vtkFeatureEdges(); oe.SetInputData(revolve(sec)); oe.BoundaryEdgesOn(); oe.FeatureEdgesOn()
        oe.SetFeatureAngle(60); oe.ManifoldEdgesOff(); oe.NonManifoldEdgesOff(); oe.Update()
        m.SetInputData(oe.GetOutput()); a.GetProperty().SetColor(0.25, 0.5, 0.85); a.GetProperty().SetLineWidth(1.5)
    ren.AddActor(a)
    for ang in (START, START + SWEEP):                           # the two cut faces with the real mesh
        f = rotz(sec, ang)
        mm = vtk.vtkPolyDataMapper(); mm.SetInputData(f); mm.ScalarVisibilityOff()
        aa = vtk.vtkActor(); aa.SetMapper(mm); aa.GetProperty().SetColor(*COL[region])
        aa.GetProperty().SetOpacity(1.0); aa.GetProperty().SetAmbient(0.45); aa.GetProperty().SetDiffuse(0.6)
        aa.GetProperty().EdgeVisibilityOn(); aa.GetProperty().SetLineWidth(0.6)
        aa.GetProperty().SetEdgeColor(*(np.array(COL[region]) * 0.45)); ren.AddActor(aa)
for name, x, y in probes:                                       # probes on the first cut face
    s = vtk.vtkSphereSource(); r = x * 1e3; th = math.radians(START)
    s.SetCenter(r * math.cos(th), r * math.sin(th), y * 1e3); s.SetRadius(0.09); s.SetThetaResolution(24); s.SetPhiResolution(24)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(s.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.85, 0.0, 0.1); ren.AddActor(a)
cam = ren.GetActiveCamera(); cam.SetFocalPoint(0, 0, 0.3); cam.SetViewUp(0, 0, 1)
cam.SetPosition(10.5, -27.0, 15.5); cam.SetViewAngle(30); ren.ResetCameraClippingRange()
W, H = 1700, 1250
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, H); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = n(wf.GetOutput().GetPointData().GetScalars()).reshape(H, W, -1)[::-1]


def screen(p):                                                 # world (mm) -> image pixel, for labels
    c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*p)
    x, y = c.GetComputedDisplayValue(ren); return x, H - y


fig = plt.figure(figsize=(21, 11.8)); ax = fig.add_axes([0.26, 0, 0.74, 0.94]); ax.imshow(img); ax.set_axis_off()
th = math.radians(START); te = math.radians(START + SWEEP)
P = lambda r, y, t=te: (r * math.cos(t), r * math.sin(t), y)
lab = [("流体：ガス\n（ノズルから出る窒素＋周りの空気）", P(5.5, 2.2), 0.80, "#1565C0"),
       ("固体：ノズル", P(2.2, 3.0), 0.62, "#B85C00"),
       ("固体：加工する板（鋼）", P(5.0, -1.0), 0.47, "#424242"),
       ("流体：切断溝\n（ガスが通る穴、半径 0.15 mm）", P(0.05, -1.2), 0.32, "#1565C0"),
       ("流体：板の下に抜けるガス", P(5.0, -4.0), 0.16, "#1565C0")]
for t, p, fy, c in lab:
    x, y = screen(p)
    ax.annotate(t, xy=(x, y), xycoords="data", xytext=(0.01, fy), textcoords="figure fraction", fontsize=19,
                color=c, fontweight="bold", va="center", arrowprops=dict(arrowstyle="->", color=c, lw=2), zorder=5)
for name, x, y in probes:
    px, py = screen((x * 1e3 * math.cos(th), x * 1e3 * math.sin(th), y * 1e3))
    ax.text(px + 10, py - 8, name, fontsize=15, color="#B71C1C", fontweight="bold")
ax.legend(handles=[Patch(color=COL["air"], label="流体：ガス（N₂/O₂ の混合ガス）"),
                   Patch(color=COL["nozzle"], label="固体：ノズル"),
                   Patch(color=COL["plate"], label="固体：加工する板（鋼）")],
          loc="upper right", fontsize=17, frameon=True, framealpha=0.9)
fig.suptitle("CHT の領域分け：どこが流体でどこが固体か（中心軸のまわりに 270° 回して表示、切り口に計算メッシュ）",
             fontsize=20, y=0.985)
fig.savefig(out, dpi=100, facecolor="white", bbox_inches="tight"); print(out)
