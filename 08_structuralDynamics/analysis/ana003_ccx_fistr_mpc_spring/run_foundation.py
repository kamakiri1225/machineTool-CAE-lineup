"""Generate and run 4-foundation-spring variant in both actual solvers."""
from pathlib import Path
import sys,shutil,tempfile,json,re
import numpy as np
from foundation_model import build
from frequency_response import execute,ALPHA,BETA,NMODES
import machine_model as mm
D=Path(__file__).resolve().parent;B=D/'foundation';R=B/'results';R.mkdir(parents=True,exist_ok=True);m=build()
(R/'supports.json').write_text(json.dumps(m.foundations,indent=2))
c=B/'ccx';f=B/'fistr';c.mkdir(exist_ok=True);f.mkdir(exist_ok=True)
mm.write_ccx(m,c/'response.inp',nmodes=NMODES)
txt=(c/'response.inp').read_text()
# Suppress the sweep full-field output; FrontISTR supplies full fields for the animations.
for dof in (1,2,3):
 txt+='\n'.join(['*STEP','*STEADY STATE DYNAMICS','5.0, 750.0, 16, 1.0','*MODAL DAMPING, RAYLEIGH',f',,{ALPHA},{BETA}','*CLOAD, OP=NEW',f'MON, {dof}, 1.0','*NODE PRINT, NSET=MON','U','*NODE FILE, FREQUENCY=0','U','*END STEP'])+'\n'
(c/'response.inp').write_text(txt)
execute([Path.home()/'opt/ccx.sh','-i','response'],c,'solve.log')
mats=mm.write_fistr_mesh(m,f/'machine.msh')
# Explicit face node sets: the actual MPC equations still expand individual node IDs.
txt=(f/'machine.msh').read_text();groups=[]
for foot in m.foundations:groups.extend([f"!NGROUP, NGRP={foot['name']}_FACE",','.join(map(str,foot['surface_nodes']))])
(f/'machine.msh').write_text(txt.replace('!END','\n'.join(groups)+'\n!END'))
mm.write_fistr_eigen_cnt(f/'eigen.cnt',mats,nmodes=NMODES);mm.write_fistr_ctrl(f/'hecmw_ctrl.dat','machine.msh','eigen.cnt','eigen.res')
exe=Path.home()/'opt/frontistr59/bin/fistr1';execute([exe],f,'eigen.out');shutil.copyfile(f/'0.log',f/'eigen_0.log')
# Temporary output location keeps per-frequency full fields out of the workspace.
work=Path(tempfile.mkdtemp(prefix='ana003_foundation_'));shutil.copyfile(f/'machine.msh',work/'machine.msh')
selected=[]
for dof in (1,2,3):
 name=f'sweep_{dof}';mm.write_fistr_freq_cnt(work/f'{name}.cnt',mats,5.,750.,745,ALPHA,BETA,m.monitor,NMODES,load_dof=dof)
 txt=(work/f'{name}.cnt').read_text().replace(' eigen_0.log',' '+str(f/'eigen_0.log')).replace(' 1, 1, 1, 1, 1, 1',' 1, 0, 0, 0, 0, 0').replace('!SOLVER, METHOD=DIRECT','!WRITE,RESULT\n!SOLVER, METHOD=DIRECT');(work/f'{name}.cnt').write_text(txt)
 mm.write_fistr_ctrl(work/'hecmw_ctrl.dat','machine.msh',f'{name}.cnt',f'{name}.res',eig_res=str(f/'eigen.res'),dyna=f'{name}_dyna.res')
 execute([exe],work,f'{name}.out');records=[];files=[]
 for path in sorted(work.glob(f'{name}.res.0.*'),key=lambda p:int(p.name.split('.')[-1])):
  t=path.read_text();assert '\ndisplacement\n' in t
  freq=float(t.split('*data')[0].splitlines()[-1]);match=re.search(rf'^{m.monitor}\s*\n([^\n]+)',t,re.M);assert match
  records.append([freq]+list(map(float,match[1].split()[:3])));files.append(path)
 arr=np.array(records);assert len(arr)==745 and np.isfinite(arr).all();np.savetxt(R/f'rbe2_{dof}_fistr_components.csv',arr*np.array([1,1000,1000,1000]),delimiter=',',header='frequency_Hz,Ux_amplitude_um_per_N,Uy_amplitude_um_per_N,Uz_amplitude_um_per_N')
 # Illustrate the main foundation-related low-frequency response for each force direction.
 low=arr[:,0]<=100;freq=float(arr[low,0][np.argmax(arr[low,dof])]);selected.append(freq)
 for path in files:path.unlink()
 for path in work.glob(f'{name}_dyna.res.0.*'):path.unlink()
 for suffix in ['cnt','out']:shutil.copyfile(work/f'{name}.{suffix}',f/f'{name}.{suffix}')
 shutil.copyfile(work/'0.log',f/f'{name}.log')
 name=f'shape_{dof}';mm.write_fistr_freq_cnt(f/f'{name}.cnt',mats,freq-.001,freq,1,ALPHA,BETA,m.monitor,NMODES,load_dof=dof)
 txt=(f/f'{name}.cnt').read_text().replace(f'{freq-.001}, {freq}, 1, {freq-.001}',f'{freq-.001}, {freq}, 1, {freq}').replace(' 0.0, 1.0',f' 0.0, {1/freq}').replace(' 1, 1, 1, 1, 1, 1',' 1, 0, 0, 0, 0, 0').replace('!SOLVER, METHOD=DIRECT','!WRITE,RESULT\n!SOLVER, METHOD=DIRECT');(f/f'{name}.cnt').write_text(txt)
 mm.write_fistr_ctrl(f/'hecmw_ctrl.dat','machine.msh',f'{name}.cnt',f'{name}.res',eig_res='eigen.res',dyna=f'{name}_dyna.res');execute([exe],f,f'{name}.out');shutil.copyfile(f/'0.log',f/f'{name}.log')
 print('Foundation components and shape',dof,freq,flush=True)
(R/'selected.json').write_text(json.dumps(selected));shutil.rmtree(work)
print('Foundation solver runs completed.',flush=True)
