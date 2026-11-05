#!/usr/bin/env python3
"""ParaView-like view of the computational wedge (surface with mesh edges + axes triad).
usage: render_wedge_mesh.py <case> <out.png>"""
import os, sys
import vtk
case, out = os.path.abspath(sys.argv[1]), sys.argv[2]
rd = vtk.vtkOpenFOAMReader(); rd.SetFileName(os.path.join(case, "post.foam")); rd.UpdateInformation()
rd.DisableAllPatchArrays(); rd.SetPatchArrayStatus("internalMesh", 1); rd.Update()
g = vtk.vtkGeometryFilter(); g.SetInputData(rd.GetOutput().GetBlock(0)); g.Update()
tr = vtk.vtkTransform(); tr.Scale(1e3, 1e3, 1e3)
f = vtk.vtkTransformPolyDataFilter(); f.SetInputConnection(g.GetOutputPort()); f.SetTransform(tr); f.Update()
ren = vtk.vtkRenderer(); ren.SetBackground(1, 1, 1)
m = vtk.vtkPolyDataMapper(); m.SetInputConnection(f.GetOutputPort()); m.ScalarVisibilityOff()
a = vtk.vtkActor(); a.SetMapper(m); p = a.GetProperty(); p.SetColor(0.85, 0.85, 0.88)
p.SetEdgeVisibility(1); p.SetEdgeColor(0.0, 0.0, 0.55); p.SetLineWidth(0.6); ren.AddActor(a)
cam = ren.GetActiveCamera(); cam.SetFocalPoint(85, 290, 0); cam.SetPosition(-250, 650, 1350); cam.SetViewUp(0, 1, 0); cam.SetViewAngle(30)
ren.ResetCamera(); cam.Zoom(0.88); ren.ResetCameraClippingRange()
rw = vtk.vtkRenderWindow(); rw.SetOffScreenRendering(1); rw.SetSize(600, 1000); rw.AddRenderer(ren); rw.SetMultiSamples(8)
ax = vtk.vtkAxesActor(); ax.SetTotalLength(1, 1, 1); ax.SetShaftTypeToLine()
om = vtk.vtkOrientationMarkerWidget(); om.SetOrientationMarker(ax); om.SetViewport(0.0, 0.0, 0.22, 0.16)
iren = vtk.vtkRenderWindowInteractor(); iren.SetRenderWindow(rw); om.SetInteractor(iren); om.EnabledOn(); om.InteractiveOff()
rw.Render()
wf = vtk.vtkWindowToImageFilter(); wf.SetInput(rw); wf.Update()
w = vtk.vtkPNGWriter(); w.SetFileName(out); w.SetInputConnection(wf.GetOutputPort()); w.Write(); print(out)
