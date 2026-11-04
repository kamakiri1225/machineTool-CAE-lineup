#!/usr/bin/env python3
"""GIF for the factory case: oblique 3-D view (air temperature on the mid plane + machine surface temperature)
and the TCP (machine centre) temperature history of machine A (bare) and B (insulated cover).
usage: render_factory.py <caseDir>"""
import io, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, vtk
from PIL import Image
from vtk.util.numpy_support import vtk_to_numpy
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17,
                     "axes.titlesize": 19, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15})
case = os.path.abspath(sys.argv[1]); OUT = os.path.join(case, "images"); os.makedirs(OUT, exist_ok=True)
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays()
for n in ("/air/internalMesh", "/machineA/patch/machineA_to_air", "/machineB/patch/machineB_to_air", "/air/patch/door", "/air/patch/vent"):
    rd.SetPatchArrayStatus(n, 1)
rd.CreateCellToPointOn()
tv = rd.GetTimeValues(); times = [tv.GetValue(i) for i in range(tv.GetNumberOfTuples())]


def load_probe(name):
    import glob
    fs = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "*", "T")), key=lambda f: float(f.split("/")[-2]))
    return np.vstack([np.loadtxt(f, comments="#") for f in fs])


pa, pb = load_probe("TCP"), load_probe("TCPB")
lut = vtk.vtkLookupTable(); lut.SetHueRange(0.667, 0.0); lut.SetTableRange(20, 30); lut.Build()


def block(path):
    it = rd.GetOutput().NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        full = "/".join(filter(None, [it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME())]))
        if full == path[-1]:
            parent = it.GetCurrentFlatIndex()
        it.GoToNextItem()


def get(region, kind, name):
    mb = rd.GetOutput()
    for i in range(mb.GetNumberOfBlocks()):
        if mb.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME()) != region:
            continue
        b = mb.GetBlock(i)
        if kind == "internalMesh":
            for k in range(b.GetNumberOfBlocks()):
                if b.GetMetaData(k).Get(vtk.vtkCompositeDataSet.NAME()) == "internalMesh":
                    return b.GetBlock(k)
        for k in range(b.GetNumberOfBlocks()):
            if b.GetMetaData(k).Get(vtk.vtkCompositeDataSet.NAME()) == "boundary":
                bb = b.GetBlock(k)
                for q in range(bb.GetNumberOfBlocks()):
                    if bb.GetMetaData(q).Get(vtk.vtkCompositeDataSet.NAME()) == name:
                        return bb.GetBlock(q)


def actor_c(pd, arr):
    from vtk.util.numpy_support import numpy_to_vtk
    a0 = pd.GetCellData().GetArray(arr)
    if a0 is None:
        return None
    c = numpy_to_vtk(vtk_to_numpy(a0) - 273.15, deep=True); c.SetName("Tc"); pd.GetCellData().AddArray(c)
    m = vtk.vtkPolyDataMapper(); m.SetInputData(pd); m.SetScalarModeToUseCellFieldData(); m.SelectColorArray("Tc")
    m.SetLookupTable(lut); m.UseLookupTableScalarRangeOn()
    a = vtk.vtkActor(); a.SetMapper(m); return a


frames = []
for t in times:
    if t <= 0:
        continue
    rd.UpdateTimeStep(t); rd.Update()
    ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
    air = get("air", "internalMesh", None)
    ol = vtk.vtkOutlineFilter(); ol.SetInputData(air)
    m = vtk.vtkPolyDataMapper(); m.SetInputConnection(ol.GetOutputPort()); a = vtk.vtkActor(); a.SetMapper(m)
    a.GetProperty().SetColor(0.3, 0.3, 0.3); a.GetProperty().SetLineWidth(2); ren.AddActor(a)
    pl = vtk.vtkPlane(); pl.SetOrigin(0, 0, 2.1); pl.SetNormal(0, 0, 1)
    ct = vtk.vtkCutter(); ct.SetInputData(air); ct.SetCutFunction(pl); ct.Update()
    sl = ct.GetOutput()
    from vtk.util.numpy_support import numpy_to_vtk
    tc = numpy_to_vtk(vtk_to_numpy(sl.GetPointData().GetArray("T")) - 273.15, deep=True); tc.SetName("Tc"); sl.GetPointData().AddArray(tc)
    m = vtk.vtkPolyDataMapper(); m.SetInputData(sl); m.SetScalarModeToUsePointFieldData(); m.SelectColorArray("Tc")
    m.SetLookupTable(lut); m.UseLookupTableScalarRangeOn(); a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetOpacity(0.55)
    ren.AddActor(a)
    # velocity vectors on the mid plane (frozen flow -> same in every frame)
    sl.GetPointData().SetActiveVectors("U")
    mk = vtk.vtkMaskPoints(); mk.SetInputData(sl); mk.SetOnRatio(4); mk.RandomModeOff(); mk.GenerateVerticesOn(); mk.Update()
    ar = vtk.vtkArrowSource(); ar.SetTipLength(0.35); ar.SetTipRadius(0.12); ar.SetShaftRadius(0.04)
    gl = vtk.vtkGlyph3D(); gl.SetInputConnection(mk.GetOutputPort()); gl.SetSourceConnection(ar.GetOutputPort())
    gl.SetVectorModeToUseVector(); gl.SetScaleModeToScaleByVector(); gl.SetScaleFactor(0.6); gl.OrientOn(); gl.ClampingOn()
    gl.SetRange(0, 1.5); gl.Update()
    gm = vtk.vtkPolyDataMapper(); gm.SetInputConnection(gl.GetOutputPort()); gm.ScalarVisibilityOff()
    ga = vtk.vtkActor(); ga.SetMapper(gm); ga.GetProperty().SetColor(0.05, 0.05, 0.05); ren.AddActor(ga)
    for reg in ("machineA", "machineB"):
        pd = get(reg, "patch", f"{reg}_to_air")
        if pd is not None:
            g = vtk.vtkGeometryFilter(); g.SetInputData(pd); g.Update(); ac = actor_c(g.GetOutput(), "T")
            if ac:
                ac.GetProperty().SetOpacity(0.45); ren.AddActor(ac)
    PROBES = [((2.25, 0.9, 2.1), "A 中心", (0.78, 0.16, 0.16)), ((2.25, 1.75, 2.1), "A 上面の近く", (0.78, 0.16, 0.16)),
              ((4.05, 0.9, 2.1), "B 中心", (0.08, 0.40, 0.75)), ((4.05, 1.75, 2.1), "B 上面の近く", (0.08, 0.40, 0.75))]
    for pt, _, col in PROBES:
        sp = vtk.vtkSphereSource(); sp.SetCenter(*pt); sp.SetRadius(0.09); sp.SetThetaResolution(20); sp.SetPhiResolution(20)
        m = vtk.vtkPolyDataMapper(); m.SetInputConnection(sp.GetOutputPort()); a = vtk.vtkActor(); a.SetMapper(m)
        a.GetProperty().SetColor(*col); a.GetProperty().SetAmbient(0.5); ren.AddActor(a)
    for nm, col in (("door", (0.9, 0.2, 0.2)), ("vent", (0.2, 0.2, 0.2))):
        pd = get("air", "patch", nm)
        if pd is not None:
            g = vtk.vtkGeometryFilter(); g.SetInputData(pd); m = vtk.vtkPolyDataMapper(); m.SetInputConnection(g.GetOutputPort())
            m.ScalarVisibilityOff(); a = vtk.vtkActor(); a.SetMapper(m); a.GetProperty().SetColor(*col); ren.AddActor(a)
    cam = ren.GetActiveCamera(); cam.SetFocalPoint(3.0, 1.0, 2.1); cam.SetPosition(-2.5, 6.5, 11.5); cam.SetViewUp(0, 1, 0)
    cam.SetViewAngle(30); ren.ResetCameraClippingRange()
    rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(1100, 850); rw.AddRenderer(ren); rw.SetMultiSamples(8); rw.Render()
    wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
    img = vtk_to_numpy(wf.GetOutput().GetPointData().GetScalars()).reshape(850, 1100, -1)[::-1]
    fig = plt.figure(figsize=(19, 8.2))
    ax = fig.add_axes([0.0, 0.0, 0.52, 0.92]); ax.imshow(img); ax.set_axis_off()
    kw = dict(fontsize=16, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#888", alpha=0.9))
    ax.text(20, 60, "赤：扉（30 ℃ の外気 1 m/s）　黒：換気口", **kw)
    for pt, lab, col in PROBES:
        c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*pt); d = c.GetComputedDoubleDisplayValue(ren)
        x, y = d[0], 850 - d[1]
        ax.annotate(f"● {lab}", xy=(x, y), xytext=(x + (-230 if "A" in lab else 150), y + (-70 if "上面" in lab else 60)),
                    fontsize=15, color=col, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=col, alpha=0.95),
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6))
    ax.text(20, 820, "左の機械 A：むき出し　右の機械 B：断熱カバー（30 mm）\n色：空気の温度（中央の断面）と機械の表面温度（半透明）　黒い矢印：流速　●：温度を測っている点", **kw)
    cax = fig.add_axes([0.53, 0.15, 0.012, 0.65])
    sm = plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(20, 30)); cb = fig.colorbar(sm, cax=cax); cb.set_label("温度 [℃]")
    a2 = fig.add_axes([0.63, 0.12, 0.35, 0.72])
    for d, lab, c in ((pa, "機械 A（むき出し）", "#C62828"), (pb, "機械 B（断熱カバー）", "#1565C0")):
        k = d[:, 0] <= t
        a2.plot(d[:, 0] / 3600, d[:, 1] - 273.15, color=c, lw=1, alpha=0.25)
        a2.plot(d[k, 0] / 3600, d[k, 1] - 273.15, color=c, lw=3, label=lab + "：中心")
        a2.plot(d[k, 0] / 3600, d[k, 2] - 273.15, color=c, lw=2, ls="--", label=lab + "：上面の近く")
    a2.set_xlim(0, max(pa[:, 0]) / 3600); a2.set_ylim(19.8, 24.5)
    a2.set_xlabel("時間 [h]"); a2.set_ylabel("機械の温度 [℃]"); a2.legend(frameon=False, fontsize=14, loc="upper left")
    a2.set_title("機械の温度の上がり方（●の点の温度）")
    fig.suptitle(f"t = {t/3600:4.2f} h　30 ℃ の工場に置いた 20 ℃ の機械（床 20 ℃）", fontsize=21, x=0.02, ha="left")
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=80, facecolor="white"); plt.close(fig); buf.seek(0)
    frames.append(Image.open(buf).convert("RGB"))
frames[-1].save(os.path.join(OUT, "factory_final.png"))
pal = frames[-1].quantize(colors=200, method=Image.Quantize.MEDIANCUT)
q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
g = os.path.join(OUT, "factory_ambient.gif")
q[0].save(g, save_all=True, append_images=q[1:], duration=[250] * (len(q) - 1) + [1800], loop=0, optimize=True)
print(g, len(q))
