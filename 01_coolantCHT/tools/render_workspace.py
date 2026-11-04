#!/usr/bin/env python3
"""Render the workspace coolant case (chtMultiRegionTwoPhaseEulerFoam) from an oblique view.

  - enclosure: outline + translucent floor, drains (dark) and coolant inlet (red)
  - table + workpiece (solid region): surface coloured by temperature
  - coolant: iso-surface alpha.liquid = ALPHA_ISO (blue)
Writes one PNG per time step into <case>/images/ and an animated GIF.
usage: python3 render_workspace.py <caseDir> [alphaIso]
VTK (python) is used directly because pvbatch hangs in this WSL environment.
"""
import os
import sys

import vtk
from PIL import Image

case = os.path.abspath(sys.argv[1])
ALPHA_ISO = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
FONT = os.path.expanduser("~/.fonts/NotoSansCJKjp-Regular.otf")
FONTB = os.path.expanduser("~/.fonts/NotoSansCJKjp-Bold.otf")
OUT = os.path.join(case, "images")
os.makedirs(OUT, exist_ok=True)
W, H = 1900, 1150
TMIN, TMAX = 20.0, 25.0     # [degC]; values above TMAX are shown in the top colour

rd = vtk.vtkOpenFOAMReader()
rd.SetFileName(os.path.join(case, "post.foam"))
rd.UpdateInformation()
rd.DisableAllPatchArrays()
for n in ("/water/internalMesh", "/solid/internalMesh", "/solid/patch/solid_to_water", "/water/patch/drain", "/water/patch/inlet", "/water/patch/minY"):
    rd.SetPatchArrayStatus(n, 1)
rd.DisableAllCellArrays()
for n in ("alpha.liquid", "T", "T.liquid"):
    rd.SetCellArrayStatus(n, 1)
rd.CreateCellToPointOff()   # cell->point is done below from interior cells only (no boundary values)
times = [rd.GetTimeValues().GetValue(i) for i in range(rd.GetTimeValues().GetNumberOfTuples())]


def block(mb, path):
    """walk the multiblock by names like ['water','internalMesh']"""
    cur = mb
    for name in path:
        nxt = None
        for i in range(cur.GetNumberOfBlocks()):
            if cur.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME()) == name:
                nxt = cur.GetBlock(i)
                break
        if nxt is None:
            return None
        cur = nxt
    return cur


def text(ren, s, x, y, size, bold=False, color=(0.1, 0.1, 0.1)):
    t = vtk.vtkTextActor()
    t.SetInput(s)
    p = t.GetTextProperty()
    p.SetFontFamily(vtk.VTK_FONT_FILE)
    p.SetFontFile(FONTB if bold else FONT)
    p.SetFontSize(size)
    p.SetColor(*color)
    t.SetDisplayPosition(x, y)
    ren.AddActor2D(t)


lut = vtk.vtkLookupTable()
lut.SetHueRange(0.667, 0.0)
lut.SetTableRange(TMIN, TMAX)
lut.Build()

frames = []
for t in times:
    rd.UpdateTimeStep(t)
    rd.Update()
    mb = rd.GetOutput()
    water = block(mb, ["water", "internalMesh"])
    solid = block(mb, ["solid", "internalMesh"])
    c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(water); c2p.PassCellDataOn(); c2p.Update(); water = c2p.GetOutput()
    c2s = vtk.vtkCellDataToPointData(); c2s.SetInputData(solid); c2s.PassCellDataOn(); c2s.Update(); solid = c2s.GetOutput()
    ren = vtk.vtkRenderer()
    ren.SetBackground(1, 1, 1)

    # enclosure outline
    ol = vtk.vtkOutlineFilter()
    ol.SetInputData(water)
    m = vtk.vtkPolyDataMapper()
    m.SetInputConnection(ol.GetOutputPort())
    a = vtk.vtkActor()
    a.SetMapper(m)
    a.GetProperty().SetColor(0.3, 0.3, 0.3)
    a.GetProperty().SetLineWidth(2)
    ren.AddActor(a)
    for path, col, op in ((["water", "patch", "minY"], (0.85, 0.85, 0.85), 0.5),
                          (["water", "patch", "drain"], (0.35, 0.35, 0.35), 0.9),
                          (["water", "patch", "inlet"], (0.9, 0.1, 0.1), 1.0)):
        b = block(mb, path)
        if b is None:
            continue
        g = vtk.vtkGeometryFilter()
        g.SetInputData(b)
        m = vtk.vtkPolyDataMapper()
        m.SetInputConnection(g.GetOutputPort())
        m.ScalarVisibilityOff()
        a = vtk.vtkActor()
        a.SetMapper(m)
        a.GetProperty().SetColor(*col)
        a.GetProperty().SetOpacity(op)
        ren.AddActor(a)

    # solid (table + workpiece): surface coloured by the fluid-solid interface temperature (face values)
    iface = block(mb, ["solid", "patch", "solid_to_water"])
    g = vtk.vtkGeometryFilter()
    g.SetInputData(iface if iface is not None and iface.GetCellData().GetArray("T") is not None else solid)
    g.Update()
    surfpd = g.GetOutput()
    if surfpd.GetCellData().GetArray("T") is not None:
        from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
        tc = numpy_to_vtk(vtk_to_numpy(surfpd.GetCellData().GetArray("T")) - 273.15, deep=True)
        tc.SetName("Tc")
        surfpd.GetCellData().AddArray(tc)
    m = vtk.vtkPolyDataMapper()
    m.SetInputData(surfpd)
    m.SetScalarModeToUseCellFieldData()
    m.SelectColorArray("Tc")
    m.SetLookupTable(lut)
    m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor()
    a.SetMapper(m)
    g.Update()
    if surfpd.GetCellData().GetArray("Tc") is None:      # t = 0: uniform initial temperature
        m.ScalarVisibilityOff()
        c0 = [0.0, 0.0, 0.0]; lut.GetColor(TMIN, c0); a.GetProperty().SetColor(*c0)
    ren.AddActor(a)
    # right panel: surface temperature only (no coolant), close-up
    ren2 = vtk.vtkRenderer(); ren2.SetBackground(1, 1, 1); ren2.SetViewport(0.56, 0.0, 1.0, 1.0)
    m2 = vtk.vtkPolyDataMapper(); m2.SetInputData(surfpd); m2.SetScalarModeToUseCellFieldData(); m2.SelectColorArray("Tc")
    m2.SetLookupTable(lut); m2.UseLookupTableScalarRangeOn()
    a2 = vtk.vtkActor(); a2.SetMapper(m2)
    if surfpd.GetCellData().GetArray("Tc") is None:
        m2.ScalarVisibilityOff(); c0 = [0.0, 0.0, 0.0]; lut.GetColor(TMIN, c0); a2.GetProperty().SetColor(*c0)
    ren2.AddActor(a2)
    ed = vtk.vtkFeatureEdges()
    ed.SetInputData(surfpd)
    ed.BoundaryEdgesOff(); ed.ManifoldEdgesOff(); ed.NonManifoldEdgesOff(); ed.FeatureEdgesOn()
    em = vtk.vtkPolyDataMapper(); em.SetInputConnection(ed.GetOutputPort()); em.ScalarVisibilityOff()
    ea = vtk.vtkActor(); ea.SetMapper(em); ea.GetProperty().SetColor(0.15, 0.15, 0.15); ea.GetProperty().SetLineWidth(1.5)
    ren.AddActor(ea)
    ea2 = vtk.vtkActor(); ea2.SetMapper(em); ea2.GetProperty().SetColor(0.15, 0.15, 0.15); ea2.GetProperty().SetLineWidth(1.5)
    ren2.AddActor(ea2)

    # coolant: two iso-surfaces of the liquid volume fraction
    water.GetPointData().SetActiveScalars("alpha.liquid")
    amax = water.GetPointData().GetArray("alpha.liquid").GetRange()[1]
    for level, col, op in ((ALPHA_ISO, (0.55, 0.78, 1.0), 0.35), (0.5, (0.05, 0.3, 0.85), 0.85)):
        if amax <= level:
            continue
        c = vtk.vtkContourFilter()
        c.SetInputData(water)
        c.SetValue(0, level)
        sm = vtk.vtkWindowedSincPolyDataFilter()
        sm.SetInputConnection(c.GetOutputPort())
        sm.SetNumberOfIterations(15)
        nm = vtk.vtkPolyDataNormals()
        nm.SetInputConnection(sm.GetOutputPort())
        m = vtk.vtkPolyDataMapper()
        m.SetInputConnection(nm.GetOutputPort())
        m.ScalarVisibilityOff()
        a = vtk.vtkActor()
        a.SetMapper(m)
        a.GetProperty().SetColor(*col)
        a.GetProperty().SetOpacity(op)
        a.GetProperty().SetSpecular(0.3)
        ren.AddActor(a)

    # camera: oblique view from front-right-above, y up
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(0.36, 0.22, 0.22)
    cam.SetPosition(1.80, 1.30, 1.80)
    cam.SetViewUp(0, 1, 0)
    cam.SetViewAngle(30)
    ren.ResetCameraClippingRange()
    ren.SetViewport(0.0, 0.0, 0.56, 1.0)
    cam2 = ren2.GetActiveCamera()
    cam2.SetFocalPoint(0.34, 0.29, 0.20); cam2.SetPosition(1.10, 0.85, 1.10); cam2.SetViewUp(0, 1, 0); cam2.SetViewAngle(30)
    ren2.ResetCameraClippingRange()

    sb = vtk.vtkScalarBarActor()
    sb.SetLookupTable(lut)
    sb.SetNumberOfLabels(4)
    sb.SetLabelFormat("%.0f")
    sb.SetNumberOfLabels(6)
    sb.SetPosition(0.89, 0.22)
    sb.SetWidth(0.06)
    sb.SetHeight(0.55)
    sb.SetTitle("")
    for tp in (sb.GetLabelTextProperty(),):
        tp.SetFontFamily(vtk.VTK_FONT_FILE); tp.SetFontFile(FONT); tp.SetFontSize(30); tp.SetColor(0.1, 0.1, 0.1)
        tp.ItalicOff(); tp.ShadowOff(); tp.BoldOff()
    sb.UnconstrainedFontSizeOn()
    ren2.AddActor2D(sb)

    text(ren, f"t = {t:4.2f} s　クーラントとワークの温度", 40, H - 80, 44, bold=True)
    text(ren, "左：クーラント（35 ℃）の流れ", 40, 140, 32, bold=True)
    text(ren, f"薄い水色：液が少しでもある範囲（α ≥ {ALPHA_ISO:g}）", 40, 90, 30)
    text(ren, "濃い青：ほぼ液だけの範囲（α ≥ 0.5）", 40, 40, 30)
    text(ren2, "右：表面温度（クーラントを消して拡大）", int(W * 0.57), 140, 32, bold=True)
    text(ren2, "初期 20 ℃、25 ℃ 以上は赤", int(W * 0.57), 90, 30)
    text(ren2, "表面温度 [℃]", int(W * 0.855), int(H * 0.80), 30)
    rw = vtk.vtkRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetSize(W, H)
    rw.AddRenderer(ren)
    rw.AddRenderer(ren2)
    rw.SetMultiSamples(8)
    rw.Render()
    wf = vtk.vtkWindowToImageFilter()
    wf.SetInput(rw)
    wf.Update()
    fn = os.path.join(OUT, f"frame_{t:06.3f}.png")
    pw = vtk.vtkPNGWriter()
    pw.SetFileName(fn)
    pw.SetInputConnection(wf.GetOutputPort())
    pw.Write()
    frames.append(Image.open(fn).convert("RGB"))
    print("wrote", fn)

gif = os.path.join(OUT, "coolant_workspace.gif")
small = [f.resize((W * 3 // 4, H * 3 // 4), Image.LANCZOS) for f in frames]
pal = small[-1].quantize(colors=200, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in small]
q[0].save(gif, save_all=True, append_images=q[1:], duration=[400] * (len(q) - 1) + [1500], loop=0, optimize=True)
print("gif", gif, os.path.getsize(gif) / 1e6, "MB")
