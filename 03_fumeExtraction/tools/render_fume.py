#!/usr/bin/env python3
"""Oblique-view animation of fume particles on the steady suction flow.
usage: python3 render_fume.py <runDir>
  <runDir>/flow      : simpleFoam result (streamlines, sheet, suction openings)
  <runDir>/particles : icoUncoupledKinematicParcelFoam result (lagrangian/kinematicCloud)
Writes <runDir>/images/fume_*.png and fume_extraction.gif. Uses python VTK (pvbatch hangs here).
"""
import glob
import os
import re
import sys

import numpy as np
import vtk
from PIL import Image
from vtk.util.numpy_support import numpy_to_vtk

run = os.path.abspath(sys.argv[1])
FONT = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
FONTB = os.path.expanduser("~/.fonts/NotoSansCJKjp-Bold.otf")
OUT = os.path.join(run, "images")
os.makedirs(OUT, exist_ok=True)
W, H = 1600, 1000


def text(ren, s, x, y, size, bold=False, color=(0.1, 0.1, 0.1)):
    t = vtk.vtkTextActor()
    t.SetInput(s)
    p = t.GetTextProperty()
    p.SetFontFamily(vtk.VTK_FONT_FILE); p.SetFontFile(FONTB if bold else FONT)
    p.SetFontSize(size); p.SetColor(*color)
    t.SetDisplayPosition(x, y)
    ren.AddActor2D(t)


# ---------------------------------------------------------------- steady flow (static actors)
fr = vtk.vtkOpenFOAMReader()
fr.SetFileName(os.path.join(run, "flow", "post.foam"))
fr.UpdateInformation()
ft = fr.GetTimeValues()
fr.UpdateTimeStep(ft.GetValue(ft.GetNumberOfTuples() - 1))
fr.CreateCellToPointOn()
fr.EnableAllPatchArrays()
fr.Update()
fmb = fr.GetOutput()


def find(mb, name):
    it = mb.NewIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == name:
            return it.GetCurrentDataObject()
        it.GoToNextItem()
    return None


internal = find(fmb, "internalMesh")
static = []


def surf(ds, col, op):
    g = vtk.vtkGeometryFilter(); g.SetInputData(ds)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(g.GetOutputPort()); m.ScalarVisibilityOff()
    a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*col); a.GetProperty().SetOpacity(op); a.GetProperty().SetAmbient(0.6); a.GetProperty().SetDiffuse(0.5)
    return a


ol = vtk.vtkOutlineFilter(); ol.SetInputData(internal)
m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ol.GetOutputPort())
a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(0.3, 0.3, 0.3); a.GetProperty().SetLineWidth(2)
static.append(a)
for name, col, op in (("sheet_master", (0.62, 0.62, 0.65), 1.0), ("suction", (0.15, 0.15, 0.15), 1.0),
                      ("roofOpening", (0.75, 0.85, 1.0), 0.25)):
    ds = find(fmb, name)
    if ds is not None:
        static.append(surf(ds, col, op))

# streamlines seeded under the roof and around the cut
seeds = vtk.vtkPointSource()
seeds.SetCenter(0.6, 0.25, 0.4); seeds.SetRadius(0.35); seeds.SetNumberOfPoints(220)
st = vtk.vtkStreamTracer()
internal.GetPointData().SetActiveVectors("U")
st.SetInputData(internal); st.SetSourceConnection(seeds.GetOutputPort())
st.SetMaximumPropagation(3.0); st.SetIntegrationDirectionToBoth(); st.SetIntegratorTypeToRungeKutta45()
sm = vtk.vtkPolyDataMapper(); sm.SetInputConnection(st.GetOutputPort()); sm.ScalarVisibilityOff()
sa = vtk.vtkActor(); sa.SetMapper(sm); sa.GetProperty().SetColor(0.45, 0.6, 0.85); sa.GetProperty().SetOpacity(0.35)
static.append(sa)

# ---------------------------------------------------------------- particles per time
pdir = os.path.join(run, "particles")
tdirs = sorted([d for d in os.listdir(pdir) if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0], key=float)


def read_list(path):
    s = open(path).read()
    m = re.search(r"\n(\d+)\s*\n\(", s)
    if not m:
        m2 = re.search(r"\n(\d+)\{([^}]*)\}", s)
        return np.full(int(m2.group(1)), float(m2.group(2))) if m2 else np.array([])
    body = s[m.end():s.rindex(")")]
    return body


def read_positions(path):
    s = open(path).read()
    pts = re.findall(r"\(\(([-\d.eE+ ]+)\)", s)          # barycentric format: ((x y z) cell ...)? fall back below
    if pts:
        return np.array([list(map(float, p.split())) for p in pts])
    return None


lut = vtk.vtkLookupTable(); lut.SetHueRange(0.0, 0.7); lut.SetTableRange(2, 50); lut.Build()
frames = []
nsuck_path = None
for td in tdirs:
    cdir = os.path.join(pdir, td, "lagrangian", "kinematicCloud")
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
    for a in static:
        ren.AddActor(a)
    npart = 0
    if os.path.isdir(cdir):
        rr = vtk.vtkOpenFOAMReader()
        rr.SetFileName(os.path.join(pdir, "post.foam"))
        rr.UpdateInformation()
        rr.DisableAllPatchArrays()
        for i in range(rr.GetNumberOfPatchArrays()):
            if "lagrangian" in rr.GetPatchArrayName(i):
                rr.SetPatchArrayStatus(rr.GetPatchArrayName(i), 1)
        rr.EnableAllLagrangianArrays()
        rr.UpdateTimeStep(float(td))
        rr.Update()
        cloud = None
        it = rr.GetOutput().NewIterator(); it.InitTraversal()
        while not it.IsDoneWithTraversal():
            d = it.GetCurrentDataObject()
            if d.IsA("vtkPolyData") and d.GetNumberOfPoints() > 0:
                cloud = d
            it.GoToNextItem()
        if cloud is not None:
            npart = cloud.GetNumberOfPoints()
            dd = cloud.GetPointData().GetArray("d")
            dmu = numpy_to_vtk(np.array([dd.GetValue(i) * 1e6 for i in range(npart)]), deep=True)
            dmu.SetName("dmu")
            cloud.GetPointData().AddArray(dmu)
            sph = vtk.vtkSphereSource(); sph.SetRadius(0.006); sph.SetThetaResolution(6); sph.SetPhiResolution(6)
            gl = vtk.vtkGlyph3D(); gl.SetInputData(cloud); gl.SetSourceConnection(sph.GetOutputPort())
            gl.SetScaleModeToDataScalingOff(); gl.SetColorModeToColorByScalar()
            cloud.GetPointData().SetActiveScalars("dmu")
            gm = vtk.vtkPolyDataMapper(); gm.SetInputConnection(gl.GetOutputPort())
            gm.SetLookupTable(lut); gm.UseLookupTableScalarRangeOn(); gm.SetScalarModeToUsePointData()
            ga = vtk.vtkActor(); ga.SetMapper(gm)
            ren.AddActor(ga)
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(0.65, 0.12, 0.40); cam.SetPosition(2.35, 1.45, 2.45); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
    ren.ResetCameraClippingRange()
    sb = vtk.vtkScalarBarActor(); sb.SetLookupTable(lut); sb.SetNumberOfLabels(5); sb.SetLabelFormat("%.0f")
    sb.SetPosition(0.89, 0.15); sb.SetWidth(0.05); sb.SetHeight(0.55); sb.UnconstrainedFontSizeOn()
    tp = sb.GetLabelTextProperty(); tp.SetFontFamily(vtk.VTK_FONT_FILE); tp.SetFontFile(FONT); tp.SetFontSize(22)
    tp.SetColor(0.1, 0.1, 0.1); tp.ItalicOff(); tp.ShadowOff(); tp.BoldOff()
    ren.AddActor2D(sb)
    text(ren, "粒子径 [µm]", int(W * 0.86), int(H * 0.73), 24)
    text(ren, f"t = {float(td):4.2f} s　　加工室内に浮遊中の粒子：{npart} 個", 40, H - 70, 32, bold=True)
    text(ren, "灰色：鋼板（中央に切断点）　黒：集塵室の吸引口（左右）　水色の線：定常の流れの流線　上面：外気の取り入れ口", 40, 30, 22)
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(W, H); rw.AddRenderer(ren); rw.SetMultiSamples(8)
    rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    fn = os.path.join(OUT, f"fume_{float(td):06.3f}.png")
    pw = vtk.vtkPNGWriter(); pw.SetFileName(fn); pw.SetInputConnection(wf.GetOutputPort()); pw.Write()
    frames.append(Image.open(fn).convert("RGB"))
    print(fn, npart)

gif = os.path.join(OUT, "fume_extraction.gif")
small = [f.resize((W * 3 // 4, H * 3 // 4), Image.LANCZOS) for f in frames]
pal = small[len(small) // 2].quantize(colors=200, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in small]
q[0].save(gif, save_all=True, append_images=q[1:], duration=[120] * (len(q) - 1) + [1500], loop=0, optimize=True)
print("gif", gif, os.path.getsize(gif) / 1e6, "MB")
