"""Run actual solver frequency responses; preserve the original modal study.
Usage: ~/opt/venv/bin/python frequency_response.py run
"""
from pathlib import Path
import sys, os, subprocess, shutil, re, json
import numpy as np
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'tools'))
import machine_model as mm
import mpc_springs as ms
NMODES=30
ALPHA=30.1592894745
BETA=8.4882636316e-6
F0,F1,NPT=50.,750.,350
SELECT={'rbe2':[183.,454.,625.], 'rbe3':[181.,442.,611.]}

def execute(args, cwd, logfile):
    env=os.environ.copy(); root=Path.home()/'opt/ccxroot/usr/lib/x86_64-linux-gnu'
    env['LD_LIBRARY_PATH']=':'.join(map(str,[root,root/'blas',root/'lapack']))+':'+env.get('LD_LIBRARY_PATH','')
    env['OMP_NUM_THREADS']='8'; env['CCX_NPROC_EQUATION_SOLVER']='8'
    print('Running',cwd,logfile,flush=True)
    with (cwd/logfile).open('w') as out:
        subprocess.run(list(map(str,args)),cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT,check=True)

def run():
 for kind in SELECT:
    m=ms.build(kind)
    base=D/'frf'/kind; c=base/'ccx'; f=base/'fistr'
    c.mkdir(parents=True,exist_ok=True); f.mkdir(parents=True,exist_ok=True)
    mm.write_ccx(m,c/'response.inp',nmodes=NMODES)
    txt=(c/'response.inp').read_text()
    for dof in (1,2,3):
        txt+='\n'.join(['*STEP','*STEADY STATE DYNAMICS',f'{F0}, {F1}, 12, 1.0','*MODAL DAMPING, RAYLEIGH',f',,{ALPHA},{BETA}', '*CLOAD, OP=NEW', f'MON, {dof}, 1.0','*NODE PRINT, NSET=MON','U','*NODE FILE','U','*END STEP'])+'\n'
    # FRD stores complex displacement for all sweep points; animate the selected nearest frequency.
    (c/'response.inp').write_text(txt)
    execute([Path.home()/'opt/ccx.sh','-i','response'],c,'solve.log')
    mats=mm.write_fistr_mesh(m,f/'machine.msh')
    mm.write_fistr_eigen_cnt(f/'eigen.cnt',mats,nmodes=NMODES)
    mm.write_fistr_ctrl(f/'hecmw_ctrl.dat','machine.msh','eigen.cnt','eigen.res')
    execute([Path.home()/'opt/frontistr59/bin/fistr1'],f,'eigen.out')
    shutil.copyfile(f/'0.log',f/'eigen_0.log')
    for dof, freq in enumerate(SELECT[kind],1):
        # Sweep logs are compact; selected-frequency output requests only one frequency + phase t=0.
        for tag,lo,hi,n,output in [('sweep',F0,F1,NPT,False),('shape',freq-0.001,freq,1,True)]:
            name=f'{tag}_{dof}'
            mm.write_fistr_freq_cnt(f/f'{name}.cnt',mats,lo,hi,n,ALPHA,BETA,m.monitor,NMODES,load_dof=dof)
            t=(f/f'{name}.cnt').read_text().replace(f'{lo}, {hi}, {n}, {lo}',f'{lo}, {hi}, {n}, {freq}')
            t=t.replace(' 0.0, 1.0',f' 0.0, {1/freq}')
            t=t.replace(' 1, 1, 1, 1, 1, 1',' 1, 0, 0, 0, 0, 0')
            if output: t=t.replace('!SOLVER, METHOD=DIRECT','!WRITE,RESULT\n!SOLVER, METHOD=DIRECT')
            (f/f'{name}.cnt').write_text(t)
            mm.write_fistr_ctrl(f/'hecmw_ctrl.dat','machine.msh',f'{name}.cnt',f'{name}.res',eig_res='eigen.res',dyna=f'{name}_dyna.res')
            execute([Path.home()/'opt/frontistr59/bin/fistr1'],f,f'{name}.out')
            shutil.copyfile(f/'0.log',f/f'{name}.log')
 if __name__=='__main__': print('All solver runs completed.',flush=True)

if __name__=='__main__': run()
