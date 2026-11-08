"""Extract actual FrontISTR component amplitudes from solver result files.
Temporary full-field sweep files are removed after validated monitor CSV extraction.
"""
from pathlib import Path
import re,shutil,tempfile,numpy as np
from frequency_response import execute
import machine_model as mm
D=Path(__file__).resolve().parent;src=D/'frf/rbe2/fistr';archive=D/'frf/rbe2/component_sweep59';archive.mkdir(exist_ok=True);work=Path(tempfile.mkdtemp(prefix='ana003_components_'))
shutil.copyfile(src/'machine.msh',work/'machine.msh')
for dof in (1,2,3):
 name=f'component_{dof}'
 txt=(src/f'sweep_{dof}.cnt').read_text().replace(' eigen_0.log',' '+str(src/'eigen_0.log')).replace('!SOLVER, METHOD=DIRECT','!WRITE,RESULT\n!SOLVER, METHOD=DIRECT')
 (work/f'{name}.cnt').write_text(txt)
 mm.write_fistr_ctrl(work/'hecmw_ctrl.dat','machine.msh',f'{name}.cnt',f'{name}.res',eig_res=str(src/'eigen.res'),dyna=f'{name}_dyna.res')
 execute([Path.home()/'opt/frontistr59/bin/fistr1'],work,f'{name}.out')
 shutil.copyfile(work/'0.log',work/f'{name}.log');records=[];files=[]
 for path in sorted(work.glob(f'{name}.res.0.*'),key=lambda p:int(p.name.split('.')[-1])):
  txt=path.read_text();assert '\ndisplacement\n' in txt, 'Expected v5.9 amplitude result, not complex real-only data';head=txt.split('*data')[0];freq=float(head.splitlines()[-1]);match=re.search(r'^20428\s*\n([^\n]+)',txt,re.M)
  assert match,path
  amp=list(map(float,match[1].split()[:3]));assert len(amp)==3
  records.append([freq]+amp);files.append(path)
 arr=np.array(records);assert len(arr)==350 and np.isfinite(arr).all()
 reference=np.loadtxt(D/f'results/frf/rbe2_{dof}_fistr.csv',delimiter=',')
 assert np.allclose(arr[:,0],reference[:,0])
 assert np.allclose(np.linalg.norm(arr[:,1:]*1000,axis=1),reference[:,1],rtol=2e-4,atol=1e-10)
 np.savetxt(D/f'results/frf/rbe2_{dof}_fistr_components.csv',arr*np.array([1,1000,1000,1000]),delimiter=',',header='frequency_Hz,Ux_amplitude_um_per_N,Uy_amplitude_um_per_N,Uz_amplitude_um_per_N')
 for path in files:path.unlink()
 for path in work.glob(f'{name}_dyna.res.0.*'):path.unlink()
 for suffix in ['cnt','out','log']:shutil.copyfile(work/f'{name}.{suffix}',archive/f'{name}.{suffix}')
 shutil.copyfile(work/'hecmw_ctrl.dat',archive/f'hecmw_ctrl_{dof}.dat')
 print('Saved actual FrontISTR 5.9 3-component sweep',dof,flush=True)
shutil.rmtree(work)
