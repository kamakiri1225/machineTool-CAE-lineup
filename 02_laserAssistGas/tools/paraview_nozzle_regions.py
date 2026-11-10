#!/usr/bin/env pvpython
"""Actual mesh -> ParaView cut-away scene, reusable pvsm, and geometry exports.
Run with PYTHONNOUSERSITE=1 pvpython --no-mpi script.py CASE OUTDIR.
PNG rendering is attempted separately because the local X-backed ParaView cannot render headlessly.
"""
from paraview.simple import *
from pathlib import Path
import sys,json,math
case=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
view=CreateView('RenderView');view.ViewSize=[1600,1300];view.UseColorPaletteForBackground=0;view.Background=[1,1,1]
colors={'air':[0.22,0.58,0.91],'plate':[0.66,0.69,0.72],'nozzle':[0.85,0.43,0.12]}
for region in ['air','plate','nozzle']:
    reader=OpenFOAMReader(registrationName=region+' actual mesh',FileName=str(case/'post.foam'))
    reader.MeshRegions=['/'+region+'/patch/front'];reader.CellArrays=['T'];reader.UpdatePipeline(time=0)
    merged=MergeBlocks(Input=reader)
    surf=ExtractSurface(Input=merged)
    flat=Transform(Input=surf);flat.Transform.Rotate=[0,-2.5,0];flat.Transform.Scale=[1000,1000,1000]
    flat.UpdatePipeline()
    SaveData(str(out/(region+'_section.vtp')),proxy=flat)
    if region!='air':
        rotation=RotationalExtrusion(Input=flat);rotation.Angle=270;rotation.Resolution=108;rotation.Capping=1;rotation.RotationAxis=[0,1,0]
        shown=Transform(Input=rotation);shown.Transform.Rotate=[0,135,0]
        SaveData(str(out/(region+'_cutaway.vtp')),proxy=shown)
        rep=Show(shown,view);rep.ColorArrayName=['POINTS',''];rep.DiffuseColor=colors[region];rep.AmbientColor=colors[region]
    else:
        for angle in [45,135]:
            t=Transform(Input=flat);t.Transform.Rotate=[0,angle,0]
            d=Show(t,view);d.ColorArrayName=['POINTS',''];d.DiffuseColor=colors[region];d.AmbientColor=colors[region];d.Opacity=.4
    Hide(reader,view)
params=json.loads((case/'caseParams.json').read_text())
for i,p in enumerate(params['probes']):
    r,y,z=[a*1000 for a in p['position_m']]
    a=math.radians(45)
    ball=Sphere(registrationName=p['name']);ball.Center=[r*math.cos(a),y,-r*math.sin(a)];ball.Radius=.12
    d=Show(ball,view);d.ColorArrayName=['POINTS',''];d.DiffuseColor=[.75,.02,.08]
title=Text();title.Text='CHT regions | copper nozzle + steel plate + gas';d=Show(title,view);d.FontSize=23;d.Color=[.12,.15,.18]
view.CameraPosition=[18,16,-30];view.CameraFocalPoint=[0,-.5,0];view.CameraViewUp=[0,1,0];view.CameraParallelProjection=1;view.CameraParallelScale=7.2
SaveState(str(out/'nozzle_regions.pvsm'))
print('ParaView state and actual-mesh geometry exported',flush=True)
