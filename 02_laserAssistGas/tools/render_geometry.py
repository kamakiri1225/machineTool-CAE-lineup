#!/usr/bin/env python3
"""Oblique-view geometry figure with labels for the assist-gas models.
  mode kerf3d : 3-D slot model (half model reflected to the full plate; nozzle shown as a half cut-away)
  mode axisym : axisymmetric wedge model revolved 270 deg (cut-away)
usage: python3 render_geometry.py <caseDir> kerf3d|axisym <out.png>
VTK renders the surfaces, matplotlib adds the call-outs (world points projected to the image).
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, mode, out = os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3]
P = json.load(open(os.path.join(case, "caseParams.json")))
rd = vtk.vtkOpenFOAMReader()
rd.SetFileName(os.path.join(case, "post.foam"))
rd.UpdateInformation()
rd.DisableAllPatchArrays()
for n in ("patch/nozzle", "patch/plate", "patch/inlet"):
    rd.SetPatchArrayStatus(n, 1)
rd.Update()


def get(name):
    it = rd.GetOutput().NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            g = vtk.vtkGeometryFilter(); g.SetInputData(it.GetCurrentDataObject()); g.Update()
            return g.GetOutput()
        it.GoToNextItem()


def scale(pd, s=1e3):
    t = vtk.vtkTransform(); t.Scale(s, s, s)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(t); f.Update()
    return f.GetOutput()


nozzle, plate, inlet = (scale(get(n)) for n in ("nozzle", "plate", "inlet"))


def reflect_z(pd):
    r = vtk.vtkReflectionFilter(); r.SetInputData(pd); r.SetPlaneToZMin(); r.CopyInputOn(); r.Update()
    g = vtk.vtkGeometryFilter(); g.SetInputConnection(r.GetOutputPort()); g.Update()
    return g.GetOutput()


def revolve(pd, deg=270, step=5):
    ap = vtk.vtkAppendPolyData()
    for k in range(int(deg / step)):
        t = vtk.vtkTransform(); t.RotateY(-135 + k * step)
        f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(t); f.Update()
        ap.AddInputData(f.GetOutput())
    ap.Update()
    return ap.GetOutput()


if mode == "kerf3d":
    plate = reflect_z(plate)          # full plate (slot visible); nozzle left as the z >= 0 half (cut-away)
else:
    nozzle, plate, inlet = revolve(nozzle), revolve(plate), revolve(inlet)

ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
for pd, col, op in ((plate, (0.62, 0.64, 0.68), 1.0), (nozzle, (0.85, 0.62, 0.30), 1.0), (inlet, (0.85, 0.15, 0.15), 1.0)):
    m = vtk.vtkPolyDataMapper(); m.SetInputData(pd); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*col); a.GetProperty().SetOpacity(op)
    a.GetProperty().SetAmbient(0.35); a.GetProperty().SetDiffuse(0.7); a.GetProperty().SetSpecular(0.2)
    ren.AddActor(a)
    fe = vtk.vtkFeatureEdges(); fe.SetInputData(pd); fe.BoundaryEdgesOn(); fe.FeatureEdgesOn(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff()
    fe.SetFeatureAngle(40)
    em = vtk.vtkPolyDataMapper(); em.SetInputConnection(fe.GetOutputPort()); em.ScalarVisibilityOff()
    ea = vtk.vtkActor(); ea.SetMapper(em); ea.GetProperty().SetColor(0.2, 0.2, 0.2); ea.GetProperty().SetLineWidth(1.2)
    ren.AddActor(ea)

H, T = P["H"], P["T"]
cam = ren.GetActiveCamera()
if mode == "kerf3d":
    cam.SetFocalPoint(-0.3, 0.3, 0.0); cam.SetPosition(-11.5, 8.5, -15.5)
else:
    cam.SetFocalPoint(0.0, 0.0, 0.0); cam.SetPosition(9.0, 8.0, 14.0)
cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30); ren.ResetCameraClippingRange()
W, Hh = 1700, 1150
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(x, y, z):
    c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(x, y, z)
    d = c.GetComputedDoubleDisplayValue(ren)
    return d[0], Hh - d[1]


fig = plt.figure(figsize=(W / 100, Hh / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.set_axis_off()


def label(pt, text, tx, ty, color="#111111"):
    """tx, ty: label position as fractions of the image (0..1, from the top-left)"""
    x, y = proj(*pt)
    ax.annotate(text, xy=(x, y), xytext=(tx * W, ty * Hh), fontsize=20, color=color, ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888888", alpha=0.95),
                arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.8))


RN = P.get("DE", 2 * P.get("RN", 0.75)) / 2
gasname = {"N2": "N₂", "O2": "O₂"}.get(P.get("gas", ""), "O₂" if P.get("YO2_in", 0) > 0.5 else "空気")
p0 = P.get("p0", P.get("p0_gauge", 0)) / 1e6
top = H + P.get("LN", 3.0)
label((0, top, 0), f"アシストガス入口（{gasname}、p₀ = {p0:g} MPa(G)）", 0.52, 0.04, "#B71C1C")
label((P.get("RNOZ", 3.0) * 0.95, H + 1.6, 0.0), "ノズル（先細、内側が円すい状に絞られる）", 0.03, 0.10)
label((0.0, H * 0.5, 0.0), f"ノズル出口 径 {2*RN:g} mm／スタンドオフ {H:g} mm", 0.66, 0.30)
if mode == "kerf3d":
    W_ = P["W"]
    label((-2.5, 0.0, 0.0), f"切断溝（幅 {W_:g} mm、x < 0 は切断済み）", 0.50, 0.90)
    label((0.0, 0.0, -W_ / 2), "切断の最前面（x = 0）", 0.66, 0.66)
    label((2.5, 0.0, -2.5), f"鋼板（板厚 {T:g} mm、まだ切っていない側）", 0.03, 0.47)
    ax.text(0.03 * W, 0.93 * Hh, "ノズルは断面を見せるため半分だけ表示（解析も対称面で半分）", fontsize=19, color="#333333")
else:
    label((P["ROUT"] * 0.6, 0.0, 0.0), f"鋼板（板厚 {T:g} mm）", 0.70, 0.56)
    label((0.0, -T, 0.0), f"切断溝（軸対称なので丸い穴、径 {2*P['RK']:g} mm）", 0.03, 0.85)
    ax.text(0.03 * W, 0.93 * Hh, "軸対称モデル（5°のくさび）を 270° 回転して表示", fontsize=19, color="#333333")
fig.savefig(out, dpi=100, facecolor="white")
print(out)
