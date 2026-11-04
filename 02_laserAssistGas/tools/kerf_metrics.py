#!/usr/bin/env python3
"""Paper (JSME 2017) Fig. 6 quantities for the 3-D kerf cases:
 tau  : mean wall shear stress on the cutting front (plate face x = 0, inside the slot)
 mdot : mass flow through the slot cross-section 1.1 mm below the plate top (x < 0, full slot width, x2 for the half model)
usage: kerf_metrics.py <case> [...]  -> prints a table"""
import json, os, sys
import numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy as n
for case in sys.argv[1:]:
    case = os.path.abspath(case); P = json.load(open(os.path.join(case, "caseParams.json")))
    rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
    tv = rd.GetTimeValues(); t = tv.GetValue(tv.GetNumberOfTuples() - 1)
    rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.SetPatchArrayStatus("patch/plate", 1)
    rd.CreateCellToPointOn(); rd.UpdateTimeStep(t); rd.Update()
    mb = rd.GetOutput(); m = mb.GetBlock(0)
    pl = None
    it = mb.NewIterator(); it.InitTraversal()
    while not it.IsDoneWithTraversal():
        if it.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) == "plate":
            pl = it.GetCurrentDataObject()
        it.GoToNextItem()
    W, T = P["W"] * 1e-3, P["T"] * 1e-3
    # wall shear on the cutting front
    cc = vtk.vtkCellCenters(); cc.SetInputData(pl); cc.Update(); c = n(cc.GetOutput().GetPoints().GetData())
    ws = pl.GetCellData().GetArray("wallShearStress")
    rho_w = pl.GetCellData().GetArray("rho")
    sel = (np.abs(c[:, 0]) < 2e-5) & (c[:, 1] < 0) & (c[:, 1] > -T) & (c[:, 2] < W / 2)
    tau = np.linalg.norm(n(ws)[sel], axis=1).mean() if ws is not None and sel.any() else np.nan
    # wallShearStress from reactingFoam is kinematic? (compressible -> Pa) keep as written
    # mass flow through y = -1.1 mm inside the slot
    plane = vtk.vtkPlane(); plane.SetOrigin(0, -1.1e-3, 0); plane.SetNormal(0, 1, 0)
    ct = vtk.vtkCutter(); ct.SetInputData(m); ct.SetCutFunction(plane); ct.Update(); s = ct.GetOutput()
    sz = vtk.vtkMeshQuality(); sz.SetInputData(s); sz.SetTriangleQualityMeasureToArea(); sz.SetQuadQualityMeasureToArea(); sz.Update(); s2 = sz.GetOutput(); s2.GetCellData().GetArray("Quality").SetName("Area")
    cc2 = vtk.vtkCellCenters(); cc2.SetInputData(s2); cc2.Update(); c2 = n(cc2.GetOutput().GetPoints().GetData())
    p2c = vtk.vtkPointDataToCellData(); p2c.SetInputData(s2); p2c.Update(); cd = p2c.GetOutput().GetCellData()
    A = n(s2.GetCellData().GetArray("Area")); U = n(cd.GetArray("U")); rho = n(cd.GetArray("rho")) if cd.GetArray("rho") else None
    if rho is None:
        Tt = n(cd.GetArray("T")); pp = n(cd.GetArray("p")); rho = pp / (296.8 * Tt)
    k = (c2[:, 0] < 0) & (c2[:, 2] < W / 2)
    mdot = -2 * np.sum(rho[k] * U[k, 1] * A[k]) * 1e3
    print(f"p0={P['p0']/1e6:4.1f} MPa  iter={t:g}  tau_front={tau:8.1f} Pa  mdot_slot={mdot:6.3f} g/s")
