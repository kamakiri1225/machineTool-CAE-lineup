#!/usr/bin/env python3
"""Oblique 3-D cut-away view (270 deg) of the axisymmetric assist-gas model, built from caseParams.json.
usage: render_axisym_3d.py <case> <out.png>"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
P = json.load(open(os.path.join(case, "caseParams.json")))
RK, RN, RIN, RNOZ, H, LN, T = (P[k] for k in ("RK", "RN", "RIN", "RNOZ", "H", "LN", "T"))
RPL = 4.5  # plate radius shown


def solid(profile, color, angle=270):
    """profile: closed polygon in (r, y); revolved about the vertical (z) axis"""
    pts = vtk.vtkPoints(); poly = vtk.vtkCellArray(); poly.InsertNextCell(len(profile))
    for i, (r, y) in enumerate(profile):
        pts.InsertNextPoint(r, 0.0, y); poly.InsertCellPoint(i)
    pd = vtk.vtkPolyData(); pd.SetPoints(pts); pd.SetPolys(poly)
    ex = vtk.vtkRotationalExtrusionFilter(); ex.SetInputData(pd); ex.SetAngle(angle); ex.SetResolution(120); ex.CappingOn()
    tr = vtk.vtkTransform(); tr.RotateZ(-90)
    tf = vtk.vtkTransformPolyDataFilter(); tf.SetInputConnection(ex.GetOutputPort()); tf.SetTransform(tr)
    nm = vtk.vtkPolyDataNormals(); nm.SetInputConnection(tf.GetOutputPort()); nm.SetFeatureAngle(30)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(nm.GetOutputPort())
    a = vtk.vtkActor(); a.SetMapper(m); p = a.GetProperty()
    p.SetColor(*color); p.SetAmbient(0.35); p.SetDiffuse(0.75); p.SetSpecular(0.25)
    fe = vtk.vtkFeatureEdges(); fe.SetInputConnection(tf.GetOutputPort()); fe.SetFeatureAngle(30)
    fe.BoundaryEdgesOn(); fe.FeatureEdgesOn(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff()
    em = vtk.vtkPolyDataMapper(); em.SetInputConnection(fe.GetOutputPort()); em.ScalarVisibilityOff()
    ea = vtk.vtkActor(); ea.SetMapper(em); ea.GetProperty().SetColor(0.15, 0.15, 0.15); ea.GetProperty().SetLineWidth(1.5)
    return [a, ea]


ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
for a in solid([(RN, H), (RNOZ, H), (RNOZ, H + LN), (RIN, H + LN)], (0.85, 0.62, 0.30)):      # nozzle
    ren.AddActor(a)
for a in solid([(RK, -T), (RPL, -T), (RPL, 0), (RK, 0)], (0.62, 0.64, 0.68)):                 # plate
    ren.AddActor(a)
# gas path (translucent): nozzle bore + gap + kerf + exit plume
gas = [(0, H + LN), (RIN, H + LN), (RN, H), (RN, 0.02), (RK, 0.02), (RK, -T), (RK * 1.6, -T - 2.0), (0, -T - 2.0)]
ga = solid(gas, (0.25, 0.55, 0.95), angle=270)[0]; ga.GetProperty().SetOpacity(0.35); ren.AddActor(ga)
cam = ren.GetActiveCamera()
cam.SetFocalPoint(0, 0, -T * 0.35); cam.SetPosition(-17, -17, 13); cam.SetViewUp(0, 0, 1); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1700, 1250
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(r, y, ang=180):
    a = np.radians(ang); c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(r * np.cos(a), r * np.sin(a), y)
    d = c.GetComputedDoubleDisplayValue(ren); return d[0], Hh - d[1]


fig = plt.figure(figsize=(W / 100, Hh / 100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.set_axis_off()


def label(pt, text, tx, ty, color="#111111"):
    x, y = proj(*pt)
    ax.annotate(text, xy=(x, y), xytext=(tx * W, ty * Hh), fontsize=21, color=color, ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888888", alpha=0.95),
                arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.8))


gasn = "O₂" if P.get("YO2_in", 1) > 0.5 else "空気"
label(((RIN + 0) / 2, H + LN), f"アシストガス入口（{gasn}、全圧 p₀ = {P['p0_gauge']/1e6:g} MPa(G)）", 0.40, 0.05, "#B71C1C")
label(((RNOZ + RIN) / 2, H + LN * 0.6), "ノズル（先細、内側が円すい状に絞られる）", 0.02, 0.13)
label((RN, H), f"ノズル出口（径 {2*RN:g} mm）", 0.66, 0.30)
label((RNOZ * 0.8, H / 2), f"スタンドオフ（すき間）{H:g} mm", 0.02, 0.40)
label((RPL * 0.7, -T / 2), f"鋼板（板厚 {T:g} mm）", 0.02, 0.62)
label((RK, -T * 0.6), f"切断溝（軸対称なので丸い穴、径 {2*RK:g} mm）", 0.55, 0.72)
label((RK * 1.3, -T - 1.2), "板の下へ抜けるガス（大気）", 0.55, 0.86)
ax.text(0.02 * W, 0.95 * Hh, "青：ガスの通り道（ノズル内 → すき間 → 切断溝 → 板の下）　中心軸まわりの回転体を 90° 切り欠いて表示", fontsize=19, color="#333333")
fig.savefig(out, dpi=100, facecolor="white"); print(out)
