#!/usr/bin/env python3
"""New physical-time Foundation-13 CHT case with an actual solid nozzle.
Never reuse an existing directory. No artificial reduction in solid density.
"""
from pathlib import Path
import json,runpy,sys,re,math
case=Path(sys.argv[1]).resolve()
if case.exists(): raise SystemExit(f'Refusing existing directory: {case}')
user=json.loads(sys.argv[2]) if len(sys.argv)>2 else {}
base=Path(__file__).parent
P={'endTime':0.05,'ramp':0.0002,'YO2_in':0.,'nProc':1,'nozzleRho':8960.,'nozzleCp':385.,'nozzleKappa':400.,'mountT':300.}
P.update(user)
# Build the known fluid/plate geometry, then replace its missing nozzle volume.
geometry={k:v for k,v in P.items() if k in ['RK','RN','RIN','RNOZ','ROUT','H','LN','T','LB','p0_gauge','TFRONT','DFRONT']}
sys.argv=[str(base/'make_laser_cht_case.py'),str(case),json.dumps(geometry)]
geo=runpy.run_path(sys.argv[0],run_name='__main__')
foundation={k:v for k,v in P.items() if k in ['endTime','ramp','YO2_in','nProc']}
foundation.update(geometry)
sys.argv=[str(base/'make_laser_cht_foundation.py'),str(case),json.dumps(foundation)]
env=runpy.run_path(sys.argv[0],run_name='__main__')
w,foam,field=env['w'],env['foam'],env['field']
C=json.loads((case/'caseParams.json').read_text())
B=geo['B'][:]
B=[(n,z,h,c,'4 1 1' if n=='above' else g) for n,z,h,c,g in B]
noz=geo['hexq'](C['RN'],C['RNOZ'],C['RNOZ'],C['RIN'],C['H'],C['H']+C['LN'])
B.append(('solidNozzle','nozzle',noz,(geo['cells']['rnoz'],geo['cells']['noz']),'3 1 1'))
def face(h,side): return '('+' '.join(str(h[i]) for i in geo['FACE'][side])+')'
patches=dict(geo['patches']); patches.pop('nozzle')
patches['nozzleMount']=('wall',[face(noz,'ymax')])
for name in ['front','back']:
    typ,faces=patches[name]; patches[name]=(typ,faces+[face(noz,name)])
verts=geo['verts']
blk='\n'.join(f"    hex ({' '.join(map(str,h))}) {z} ({nr} {ny} 1) simpleGrading ({g})" for n,z,h,(nr,ny),g in B)
pat='\n'.join(f'    {n}\n    {{\n        type {typ};\n        faces\n        (\n'+ '\n'.join('            '+f for f in faces)+'\n        );\n    }' for n,(typ,faces) in patches.items())
w('system/blockMeshDict','scale 1e-3;\nvertices\n(\n'+'\n'.join(f'    ({x:.12g} {y:.12g} {z:.12g})' for x,y,z in verts)+f'\n);\nblocks\n(\n{blk}\n);\nedges ();\nboundary\n(\n{pat}\n);\n')

# The former adiabatic fluid wall is now a mapped thermal interface.
for f in (case/'0.orig/air').iterdir():
    text=f.read_text().replace('(nozzle|air_to_plate)','(air_to_nozzle|air_to_plate)')
    text=re.sub(r'(?m)^    nozzle\s*$', '    air_to_nozzle',text)
    if f.name=='T':
        text=text.replace('air_to_nozzle\n    {\n        type            zeroGradient;',
                          'air_to_nozzle\n    {\n        type            coupledTemperature;\n        value           uniform 300;')
    f.write_text(text)
field('nozzle','T','uniform 300',{'nozzleMount':{'type':'fixedValue','value':f'uniform {P["mountT"]}'},
    'nozzle_to_air':{'type':'coupledTemperature','value':'uniform 300'}},dim='[0 0 0 1 0 0 0]')
w('constant/nozzle/physicalProperties',foam({'thermoType':'constSolidThermo',
    'rho':{'type':'uniform','value':P['nozzleRho']},'Cv':{'type':'uniform','value':P['nozzleCp']},
    'kappa':{'type':'uniform','value':P['nozzleKappa']}}))
for region in ['air','plate','nozzle']:
    if region=='nozzle':
        schemes=(case/'system/plate/fvSchemes').read_text()
        solution=(case/'system/plate/fvSolution').read_text()
    else:
        schemes=(case/f'system/{region}/fvSchemes').read_text()
        solution=(case/f'system/{region}/fvSolution').read_text()
    schemes=schemes.replace('localEuler','Euler').replace('steadyState','Euler')
    if region=='air':
        schemes=schemes.replace('Gauss linearUpwindV grad(U)','Gauss upwind')
        if P.get('fluidMode')=='quasiSteady':
            schemes=schemes.replace('default         Euler;','default         steadyState;')
        solution=solution.replace('solver          GAMG;\n        smoother        GaussSeidel;','solver          PBiCGStab;\n        preconditioner  DILU;')
        solution=solution.replace('relTol          0.1;','relTol          0;').replace('relTol          0.05;','relTol          0.01;')
    else:
        solution=solution.replace('relTol          0.01;','relTol          0;').replace('0.9;','1;')
    (case/f'system/{region}').mkdir(parents=True,exist_ok=True)
    (case/f'system/{region}/fvSchemes').write_text(schemes)
    (case/f'system/{region}/fvSolution').write_text(solution)
w('system/fvSolution',foam({'PIMPLE':{'nOuterCorrectors':2}}))
rn,rin,rout,h,ln=[C[k]*1e-3 for k in ['RN','RIN','RNOZ','H','LN']]
# Points deliberately inside the solid, not boundary or gas samples.
probes=[('P1_tip',rn+0.00015,h+0.00010),
        ('P2_inner_mid',.5*(rn+rin)+0.00015,h+.5*ln),
        ('P3_outer_mid',rout-0.00015,h+.5*ln),
        ('P4_upper',rin+0.0002,h+.9*ln)]
probe_list='(\n'+'\n'.join(f'    ({r:.12g} {y:.12g} 0)' for name,r,y in probes)+'\n)'
functions={
    'nozzleProbes':{'type':'probes','libs':'("libsampling.so")','region':'nozzle','fields':'(T)',
        'interpolationScheme':'cellPoint','writeControl':'runTime','writeInterval':'1e-5','probeLocations':probe_list},
    'nozzleAverage':{'type':'volFieldValue','libs':'("libfieldFunctionObjects.so")','region':'nozzle','cellZone':'all',
        'operation':'volAverage','fields':'(T)','writeFields':'false','writeControl':'runTime','writeInterval':'1e-5'},
    'fieldBounds':{'type':'volFieldValue','libs':'("libfieldFunctionObjects.so")','region':'air','cellZone':'all',
        'operation':'min','fields':'(p T)','writeFields':'false','writeControl':'runTime','writeInterval':'1e-5'},
    'inletFlow':{'type':'surfaceFieldValue','libs':'("libfieldFunctionObjects.so")','region':'air','select':'patch','patch':'inlet',
        'operation':'sum','fields':'(phi)','writeFields':'false','writeControl':'runTime','writeInterval':'1e-5'},
    'heatNozzle':{'type':'wallHeatFlux','libs':'("libfieldFunctionObjects.so")','region':'nozzle','patches':'(nozzle_to_air nozzleMount)',
        'executeControl':'runTime','executeInterval':'1e-5','writeControl':'runTime','writeInterval':'1e-5'}}
w('system/controlDict',foam({'application':'foamMultiRun','regionSolvers':{'air':'multicomponentFluid','plate':'solid','nozzle':'solid'},
    'startFrom':'latestTime','startTime':0,'stopAt':'endTime','endTime':P['endTime'],'deltaT':'1e-9',
    'writeControl':'adjustableRunTime','writeInterval':'0.001','purgeWrite':0,'writeFormat':'binary','writePrecision':10,
    'timeFormat':'general','timePrecision':12,'runTimeModifiable':'true','adjustTimeStep':'yes','maxCo':0.3,'maxDeltaT':'1e-6',
    'functions':functions}))
if P.get('fluidMode')=='quasiSteady':
    control=case/'system/controlDict'
    txt=control.read_text().replace('deltaT          1e-9;',f'deltaT          {P.get("thermalDeltaT",5e-5)};').replace('adjustTimeStep  yes;','adjustTimeStep  no;')
    control.write_text(txt)
    w('system/fvSolution',foam({'PIMPLE':{'nOuterCorrectors':5}}))
# Physical-time inlet pressure ramp; v13 Function1 dictionary syntax.
pfile=case/'0.orig/air/p'
txt=pfile.read_text()
txt=re.sub(r'p0\s+table\s*\(\(0 .*?;',f'p0\n        {{\n            type table;\n            values ((0 101325) ({P["ramp"]} {101325+C["p0_gauge"]}) (1 {101325+C["p0_gauge"]}));\n        }}',txt,count=1)
pfile.write_text(txt)
C.update(P)
C.update(solver='foamMultiRun Foundation13; Euler all regions; actual solid density',
    materialAssumption='Copper nozzle: provisional; not user-confirmed',
    thermalAssumptions='Nozzle mounting face fixed 300 K; kerf-front plate cells fixed 1800 K. Convection/conduction only; no radiation, direct laser absorption, oxidation, or moving cut.',
    probes=[{'name':name,'position_m':[r,y,0]} for name,r,y in probes])
(case/'caseParams.json').write_text(json.dumps(C,indent=2,ensure_ascii=False))
(case/'probe_locations.csv').write_text('name,x_m,y_m,z_m\n'+'\n'.join(f'{name},{r},{y},0' for name,r,y in probes)+'\n')
(case/'Allrun').write_text('''#!/bin/bash
set -e
cd "$(dirname "$0")"
if [ -d constant/air/polyMesh ]; then echo "Already initialized; refusing reset"; exit 1; fi
blockMesh > log.blockMesh 2>&1
splitMeshRegions -cellZones -overwrite > log.splitMeshRegions 2>&1
topoSet -region plate > log.topoSet 2>&1
mkdir -p 0
cp -r 0.orig/. 0/
for region in air plate nozzle; do checkMesh -region "$region" > "log.checkMesh.$region" 2>&1; done
foamMultiRun > log.foamMultiRun 2>&1
''')
(case/'Allrun').chmod(0o755)
print('Nozzle CHT case ready',case,flush=True)
