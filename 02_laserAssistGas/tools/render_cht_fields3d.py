#!/usr/bin/env python3
"""Oblique 3-D contours of a CHT result (same view as render_cht_regions3d.py): each region's wedge section is revolved
270 deg about the axis and coloured by a field on the two cut faces and on the surfaces of revolution.
Panels: velocity magnitude, temperature, Mach number, O2 mole fraction.  Solids are grey where the field does not exist.
usage: render_cht_fields3d.py <case> <out.png> [time]"""
import math, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy as n, numpy_to_vtk
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
START, SWEEP = -45.0, 270.0
DEC = os.path.isdir(os.path.join(case, "processor0"))      # read the decomposed case directly when it exists
rd = vtk.vtkPOpenFOAMReader() if DEC else vtk.vtkOpenFOAMReader()
rd.SetFileName(os.path.join(case, "post.foam"))
if DEC: rd.SetCaseType(0)
rd.UpdateInformation()
tv = rd.GetTimeValues(); T_ = float(sys.argv[3]) if len(sys.argv) > 3 else tv.GetValue(tv.GetNumberOfTuples() - 1)
regions = [r for r in ("air", "plate", "nozzle")
           if any(rd.GetPatchArrayName(i) == f"/{r}/patch/front" for i in range(rd.GetNumberOfPatchArrays()))]


def section(region):
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus(f"/{region}/patch/front", 1); rd.CreateCellToPointOn()
    rd.UpdateTimeStep(T_); rd.Update()
    it = rd.GetOutput().NewIterator(); it.InitTraversal(); ap = vtk.vtkAppendPolyData()
    while not it.IsDoneWithTraversal():
        g = vtk.vtkGeometryFilter(); g.SetInputData(it.GetCurrentDataObject()); g.Update(); ap.AddInputData(g.GetOutput()); it.GoToNextItem()
    ap.Update(); pd = ap.GetOutput(); P = pd.GetPointData()
    p = n(pd.GetPoints().GetData()) * 1e3; r = np.hypot(p[:, 0], p[:, 2])
    pd.GetPoints().SetData(numpy_to_vtk(np.c_[r, np.zeros_like(r), p[:, 1]], deep=True))
    if region == "air":
        U = n(P.GetArray("U")); T = n(P.GetArray("T")); Y = n(P.GetArray("O2"))
        R = 8.314462618 * (Y / 0.0319988 + (1 - Y) / 0.0280134)
        cp = Y * 918 + (1 - Y) * 1040; g = cp / (cp - R)
        for name, v in (("Umag", np.linalg.norm(U, axis=1)), ("Ma", np.linalg.norm(U, axis=1) / np.sqrt(g * R * T)),
                        ("XO2", (Y / 32.0) / (Y / 32.0 + (1 - Y) / 28.01) * 100)):
            a = numpy_to_vtk(v, deep=True); a.SetName(name); P.AddArray(a)
    return pd


def rotz(pd, ang):
    t = vtk.vtkTransform(); t.RotateZ(ang)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(pd); f.SetTransform(t); f.Update(); return f.GetOutput()


def revolve(pd):
    e = vtk.vtkFeatureEdges(); e.SetInputData(pd); e.BoundaryEdgesOn(); e.FeatureEdgesOff(); e.ManifoldEdgesOff()
    e.NonManifoldEdgesOff(); e.Update()
    x = vtk.vtkRotationalExtrusionFilter(); x.SetInputData(rotz(e.GetOutput(), START)); x.SetAngle(SWEEP)
    x.SetResolution(108); x.CappingOff(); x.Update(); return x.GetOutput()


SEC = {r: section(r) for r in regions}
GREY = {"plate": (0.70, 0.72, 0.75), "nozzle": (0.80, 0.70, 0.58)}


def lut(cmap, lo, hi):
    l = vtk.vtkLookupTable(); l.SetNumberOfTableValues(256); l.SetTableRange(lo, hi); cm = plt.get_cmap(cmap)
    for i in range(256): l.SetTableValue(i, *cm(i / 255)[:3], 1)
    l.Build(); return l


def render(field, cmap, lo, hi, solids_coloured):
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1); L = lut(cmap, lo, hi)
    for r, sec in SEC.items():
        coloured = (r == "air") or solids_coloured
        pieces = [rotz(sec, START), rotz(sec, START + SWEEP)] + ([revolve(sec)] if r != "air" else [])
        for k, pdk in enumerate(pieces):
            m = vtk.vtkPolyDataMapper(); m.SetInputData(pdk)
            if coloured:
                m.SetScalarModeToUsePointFieldData(); m.SelectColorArray(field); m.SetLookupTable(L); m.UseLookupTableScalarRangeOn()
            else:
                m.ScalarVisibilityOff()
            a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetAmbient(0.45); a.GetProperty().SetDiffuse(0.6)
            if not coloured: a.GetProperty().SetColor(*GREY[r])
            ren.AddActor(a)
    # outline of the gas region
    oe = vtk.vtkFeatureEdges(); oe.SetInputData(revolve(SEC["air"])); oe.BoundaryEdgesOn(); oe.FeatureEdgesOn(); oe.SetFeatureAngle(60)
    oe.ManifoldEdgesOff(); oe.NonManifoldEdgesOff(); oe.Update()
    m = vtk.vtkPolyDataMapper(); m.SetInputData(oe.GetOutput()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.3, 0.45, 0.7); a.GetProperty().SetLineWidth(1.2); ren.AddActor(a)
    cam = ren.GetActiveCamera(); cam.SetFocalPoint(0.3, 0, -1.0); cam.SetViewUp(0, 0, 1)
    cam.SetPosition(11.5, -28.5, 13.0); cam.SetViewAngle(30); ren.ResetCameraClippingRange()
    W, H = 1300, 980
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, H); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    return n(wf.GetOutput().GetPointData().GetScalars()).reshape(H, W, -1)[::-1]


PANELS = [("Umag", "流速 |U| [m/s]（ガス）", "jet", 0, 600, False),
          ("T", "温度 [K]（ガス・板・ノズル）", "inferno", 150, 1800, True),
          ("Ma", "マッハ数（ガス）", "jet", 0, 2.5, False),
          ("XO2", "O₂ 濃度 [mol%]（周りの空気の巻き込み）", "RdYlBu", 0, 21, False)]
fig = plt.figure(figsize=(20, 13.4))
for k, (f, lab, cm, lo, hi, sc) in enumerate(PANELS):
    x0, y0 = (k % 2) * 0.5, 0.5 - (k // 2) * 0.49
    ax = fig.add_axes([x0, y0, 0.44, 0.465]); ax.imshow(render(f, cm, lo, hi, sc)); ax.set_axis_off()
    ax.set_title(lab, fontsize=20, pad=2)
    cax = fig.add_axes([x0 + 0.445, y0 + 0.08, 0.012, 0.31])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(lo, hi)), cax=cax); cb.ax.tick_params(labelsize=15)
fig.text(0.01, -0.005, "中心軸のまわりに 270° 回して表示（手前 90° は切り欠き）。灰色・肌色は固体（板・ノズル）で、その量がない領域。"
         f"  {os.path.basename(case)}、反復 {T_:g}", fontsize=16)
fig.savefig(out, dpi=95, facecolor="white", bbox_inches="tight"); print(out, "time", T_)
