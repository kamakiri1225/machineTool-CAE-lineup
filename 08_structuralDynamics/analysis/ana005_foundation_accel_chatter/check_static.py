from pathlib import Path
import re,json,numpy as np,sys
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent/'ana003_ccx_fistr_mpc_spring'))
from frequency_response import execute
from analyze import read_set
summary={}
for name in ['fixed','foundation','foundation_table_stiff']:
 out=D/name;text=(out/'relative.inp').read_text().split('*STEP')[0]
 for index,dof in enumerate([1,2,3]):
  text+='\n'.join(['*STEP','*STATIC','*CLOAD, OP=NEW',f'WORK, {dof}, 1.0',f'TOOL, {dof}, -1.0','*NODE PRINT, NSET=WORK','U','*NODE PRINT, NSET=TOOL','U','*END STEP'])+'\n'
 reuse='--reuse-static' in sys.argv and (out/'static.inp').exists() and (out/'static.inp').read_text()==text
 if not reuse:
  (out/'static.inp').write_text(text);execute([Path.home()/'opt/ccx.sh','-i','static'],out,'static.log')
 lines=(out/'static.dat').read_text().splitlines();values={'WORK':[],'TOOL':[]}
 for i,line in enumerate(lines):
  for key in values:
   if f'for set {key} and time' in line:values[key].append(np.array(list(map(float,lines[i+2].split()[1:4]))))
 static=np.stack(values['WORK'])-np.stack(values['TOOL']);dynamic=[s[1] for s in read_set(out/'relative.dat','WORK')];tool=[s[1] for s in read_set(out/'relative.dat','TOOL')]
 summary[name]={'static_diagonal_mm_N':np.diag(static).tolist(),'dynamic_5Hz_diagonal_mm_N':[float((dynamic[j][0,j]-tool[j][0,j]).real) for j in range(3)]}
 (out/'static_matrix.json').write_text(json.dumps(static.tolist()))
print(json.dumps(summary,indent=2));(D/'results/static_checks.json').write_text(json.dumps(summary,indent=2))
