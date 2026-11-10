#!/usr/bin/env python3
"""Render ParaView-exported real mesh with EGL VTK; not a ParaView screenshot."""
from pathlib import Path
import sys,json,math
import vtk
scene=Path(sys.argv[1]);case=Path(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
font='/home/kamakiri/.fonts/NotoSansCJKjp-Regular.otf'
W,H=1700,1450
ren=vtk.vtkRenderer();ren.SetBackground(1,1,1)
rw=vtk.vtkEGLRenderWindow();rw.SetOffScreenRendering(1);rw.SetSize(W,H);rw.AddRenderer(ren)
def actor(file,color,opacity=1,angle=0):
    rd=vtk.vtkXMLPolyDataReader();rd.SetFileName(str(file));rd.Update()
    if angle:
        tr=vtk.vtkTransform();tr.RotateY(angle)
        f=vtk.vtkTransformPolyDataFilter();f.SetInputData(rd.GetOutput());f.SetTransform(tr);f.Update();data=f.GetOutput()
    else: data=rd.GetOutput()
    mapper=vtk.vtkPolyDataMapper();mapper.SetInputData(data);mapper.ScalarVisibilityOff()
    a=vtk.vtkActor();a.SetMapper(mapper);a.GetProperty().SetColor(*color);a.GetProperty().SetOpacity(opacity);a.GetProperty().SetAmbient(.35)
    ren.AddActor(a);return a
actor(scene/'nozzle_cutaway.vtp',(.82,.40,.12))
actor(scene/'plate_cutaway.vtp',(.65,.69,.72))
for angle in [45,135]: actor(scene/'air_section.vtp',(.17,.55,.88),.32,angle)
def text(s,x,y,size=32,color=(.12,.16,.20)):
    a=vtk.vtkTextActor();a.SetInput(s);a.SetDisplayPosition(x,y);p=a.GetTextProperty();p.SetFontFamily(vtk.VTK_FONT_FILE);p.SetFontFile(font);p.SetFontSize(size);p.SetColor(*color);ren.AddActor2D(a);return a
texts=[]
texts.append(text('レーザー加工 CHT：固体・流体と温度プローブ',45,H-65,38))
texts.append(text('橙：ノズル固体（仮定：銅）',45,175,32,(.65,.29,.04)))
texts.append(text('灰：加工板固体（鋼）',45,130,32))
texts.append(text('青：ガス流体（窒素＋周囲空気）',45,85,32,(.08,.35,.65)))
texts.append(text('軸対称の実メッシュを270°回転表示／単位 mm／温度結果ではありません',45,36,24))
params=json.loads((case/'caseParams.json').read_text())
cam=ren.GetActiveCamera();cam.SetPosition(18,16,-30);cam.SetFocalPoint(0,-.5,0);cam.SetViewUp(0,1,0);cam.ParallelProjectionOn();cam.SetParallelScale(7.6)
for i,p in enumerate(params['probes']):
    r,y,z=[x*1000 for x in p['position_m']];a=math.pi/4;pt=[r*math.cos(a),y,-r*math.sin(a)]
    ball=vtk.vtkSphereSource();ball.SetCenter(*pt);ball.SetRadius(.13);ball.SetThetaResolution(24);ball.SetPhiResolution(24)
    mp=vtk.vtkPolyDataMapper();mp.SetInputConnection(ball.GetOutputPort());ac=vtk.vtkActor();ac.SetMapper(mp);ac.GetProperty().SetColor(.8,.02,.10);ren.AddActor(ac)
ren.ResetCameraClippingRange();rw.Render()
def save(name):
    ren.ResetCameraClippingRange();rw.Render();w=vtk.vtkWindowToImageFilter();w.SetInput(rw);w.Update();p=vtk.vtkPNGWriter();p.SetFileName(str(out/name));p.SetInputConnection(w.GetOutputPort());p.Write()
save('01_CHT領域_斜め図.png')
cam.SetFocalPoint(0,2.0,0);cam.SetPosition(12,10,-21);cam.SetParallelScale(3.3)
for t in texts: ren.RemoveActor2D(t)
text('ノズル固体内の温度プローブ（P1〜P4）',45,H-65,38)
text('P1 先端  /  P2 内壁側中央  /  P3 外壁側中央  /  P4 取付側',45,100,28)
text('赤点は同一の軸対称断面へ回転して表示。採取座標は probe_locations.csv を参照。',45,50,23)
ren.ResetCameraClippingRange();rw.Render()
for i,p in enumerate(params['probes']):
    r,y,z=[x*1000 for x in p['position_m']];a=math.pi/4;pt=[r*math.cos(a),y,-r*math.sin(a)]
    ren.SetWorldPoint(*pt,1);ren.WorldToDisplay();x,yy,_=ren.GetDisplayPoint()
    text('P'+str(i+1),int(x)+15,int(yy)+10,34,(.65,0,.05))
save('02_ノズル温度プローブ位置.png')
print('Actual mesh previews rendered with EGL VTK',flush=True)
