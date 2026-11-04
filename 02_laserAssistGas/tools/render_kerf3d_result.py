#!/usr/bin/env python3
"""Oblique 3-D view of the kerf result: nozzle (half, cut-away) and plate (reflected to full width) in grey,
Mach number on the slot-centre symmetry plane (z = 0) and on a slot cross-section (x = XS, reflected).
usage: render_kerf3d_result.py <case> <out.png> [XS_mm]"""
import json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
XS = float(sys.argv[3]) if len(sys.argv) > 3 else -0.6
P = json.load(open(os.path.join(case, "caseParams.json")))
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
rd.DisableAllPatchArrays()
for nme in ("internalMesh", "patch/nozzle", "patch/plate"):
    rd.SetPatchArrayStatus(nme, 1)
rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update(); mb = rd.GetOutput()


def get(name):
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()


def scale(ds):
    g = vtk.vtkGeometryFilter(); g.SetInputData(ds); g.Update()
    tr = vtk.vtkTransform(); tr.Scale(1e3, 1e3, 1e3)
    f = vtk.vtkTransformPolyDataFilter(); f.SetInputData(g.GetOutput()); f.SetTransform(tr); f.Update(); return f.GetOutput()


def mirror(pd):
    r = vtk.vtkReflectionFilter(); r.SetInputData(pd); r.SetPlaneToZMin(); r.CopyInputOn(); r.Update()
    g = vtk.vtkGeometryFilter(); g.SetInputConnection(r.GetOutputPort()); g.Update(); return g.GetOutput()


inner = get("internalMesh")
lut = vtk.vtkLookupTable(); lut.SetHueRange(0.667, 0.0); lut.SetTableRange(0, 2.5); lut.Build()
ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)


def cut(origin, normal):
    pl = vtk.vtkPlane(); pl.SetOrigin(*origin); pl.SetNormal(*normal)
    c = vtk.vtkCutter(); c.SetInputData(inner); c.SetCutFunction(pl); c.Update(); return scale(c.GetOutput())


def clipbox(pd, b):
    bx = vtk.vtkBox(); bx.SetBounds(*b)
    c = vtk.vtkClipPolyData(); c.SetInputData(pd); c.SetClipFunction(bx); c.InsideOutOn(); c.Update(); return c.GetOutput()


Tt, Hn = P["T"], P["H"]
def colored(pd, op=1.0):
    pd.GetPointData().SetActiveScalars("Ma")
    m = vtk.vtkPolyDataMapper(); m.SetInputData(pd); m.SetScalarModeToUsePointFieldData(); m.SelectColorArray("Ma")
    m.SetLookupTable(lut); m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetOpacity(op); a.GetProperty().SetAmbient(0.6); a.GetProperty().SetDiffuse(0.5)
    ren.AddActor(a)


colored(clipbox(cut((0, 0, 1e-7), (0, 0, 1)), (-2.5, 2.0, -Tt - 1.5, Hn + 0.8, -1, 1)))      # symmetry plane
xsec = clipbox(cut((XS * 1e-3, 0, 0), (1, 0, 0)), (-5, 5, -Tt - 1.5, Hn + 0.3, -0.1, 1.5))   # slot cross-section (2-D panel)
for name, col, op, mir in (("plate", (0.70, 0.72, 0.76), 1.0, True), ("nozzle", (0.85, 0.62, 0.30), 1.0, False)):
    pd = scale(get(name))
    flip = vtk.vtkTransform(); flip.Scale(1, 1, -1)                 # put the solid half behind the plane (z < 0)
    ff = vtk.vtkTransformPolyDataFilter(); ff.SetInputData(pd); ff.SetTransform(flip); ff.Update(); pd = ff.GetOutput()
    if mir:
        pd = clipbox(pd, (-2.5, 2.0, -Tt - 0.1, 0.1, -1.5, 0.1))
    else:
        pd = clipbox(pd, (-3.5, 3.5, Hn - 0.1, Hn + 1.2, -3.5, 0.1))
    m = vtk.vtkPolyDataMapper(); m.SetInputData(pd); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*col); a.GetProperty().SetOpacity(op); a.GetProperty().SetAmbient(0.3)
    ren.AddActor(a)
ps = vtk.vtkPlaneSource(); ps.SetOrigin(0, -Tt, 0.001); ps.SetPoint1(2.0, -Tt, 0.001); ps.SetPoint2(0, 0, 0.001)   # plate section on z = 0
m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ps.GetOutputPort()); a = vtk.vtkActor(); a.SetMapper(m)
a.GetProperty().SetColor(0.55, 0.56, 0.6); ren.AddActor(a)
cam = ren.GetActiveCamera(); cam.SetFocalPoint(-0.2, -0.1, 0.0); cam.SetPosition(3.3, 3.2, 8.6); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCameraClippingRange()
W, Hh = 1400, 1100
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, Hh); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(Hh, W, -1)[::-1]


def proj(x, y, z):
    c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(x, y, z); d = c.GetComputedDoubleDisplayValue(ren)
    return d[0], Hh - d[1]


fig = plt.figure(figsize=(21, 10))
ax = fig.add_axes([0, 0, 0.56, 0.92]); ax.imshow(img); ax.set_axis_off()
kw = dict(fontsize=17, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", alpha=0.95))
ar = dict(arrowstyle="-|>", lw=1.6)
for pt, lab, xy in (((1.6, Hn + 0.5, -0.5), "ノズル（奥の半分）", (0.70, 0.05)),
                    ((-1.5, -Tt * 0.5, 0), "切断溝の中（切断済み、x < 0）", (0.02, 0.62)),
                    ((1.2, -Tt * 0.5, -0.2), "鋼板（まだ切っていない側、x > 0）", (0.55, 0.88)),
                    ((0.05, 0.05, 0), "切断の最前面の角（噴流が当たる）", (0.45, 0.40))):
    x, y = proj(*pt)
    ax.annotate(lab, xy=(x, y), xytext=(xy[0] * W, xy[1] * Hh), arrowprops=ar, **kw)
ax.text(0.02 * W, 0.96 * Hh, "手前：溝の中心面（z = 0）のマッハ数　奥：ノズルと鋼板の半分", fontsize=17, bbox=kw["bbox"])
# right: slot cross-section x = XS (z-y plane, mirrored)
import matplotlib.tri as mtri
pts = vtk_to_numpy(xsec.GetPoints().GetData()); Ma = vtk_to_numpy(xsec.GetPointData().GetArray("Ma"))
tri = vtk.vtkTriangleFilter(); tri.SetInputData(xsec); tri.Update(); tp = tri.GetOutput()
pts = vtk_to_numpy(tp.GetPoints().GetData()); Ma = vtk_to_numpy(tp.GetPointData().GetArray("Ma"))
cells = vtk_to_numpy(tp.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
z, y = pts[:, 2], pts[:, 1]
a2 = fig.add_axes([0.60, 0.10, 0.30, 0.78]); a2.set_facecolor("#9a9a9a")
for sgn in (1, -1):
    tc = a2.tripcolor(mtri.Triangulation(sgn * z, y, cells), Ma, shading="gouraud", cmap="jet", vmin=0, vmax=2.5)
a2.set_xlim(-1.5, 1.5); a2.set_ylim(-Tt - 1.5, Hn + 0.3); a2.set_aspect("equal")
a2.set_xlabel("z [mm]（溝の幅の方向）"); a2.set_ylabel("y [mm]")
a2.set_title(f"溝を横切る断面（x = {XS:g} mm）")
for yy, lab in ((0, "板の上面"), (-Tt, "板の下面"), (Hn, "ノズル出口の高さ")):
    a2.axhline(yy, color="white", lw=0.8, ls=":"); a2.text(1.55, yy, lab, fontsize=14, va="center", clip_on=False)
cax = fig.add_axes([0.94, 0.2, 0.012, 0.55])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(0, 2.5)), cax=cax); cb.set_label("マッハ数")
fig.suptitle(f"3 次元の切断溝モデル：N₂ {P['p0']/1e6:g} MPa(G)、ノズル径 {P['DE']:g} mm、スタンドオフ {P['H']:g} mm、板厚 {P['T']:g} mm（反復 {t:g}）",
             fontsize=20, x=0.02, ha="left")
fig.savefig(out, dpi=85, facecolor="white", bbox_inches="tight"); print(out)
