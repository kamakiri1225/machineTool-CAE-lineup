"""Actual VMC solid-model cutting and actuator force-pair FRFs."""
from pathlib import Path
import sys,json
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent/'ana003_ccx_fistr_mpc_spring'))
from frequency_response import execute,ALPHA,BETA
from vmc_model import build,write_input,CONTACT,MATERIALS
metadata={}
for name in ['fixed','foundation','foundation_table_stiff']:
 m=build(name);out=D/name;out.mkdir(exist_ok=True);write_input(m,out/'relative.inp')
 nut=next(s for s in m.springs if s[0]=='NUT' and s[1]==1)
 metadata[name]={'tool_tip_node':m.contact['TOOL']['ref'],'work_contact_node':m.monitor,'contact_xyz_mm':CONTACT.tolist(),'contact_patches':m.contact,'mass':m.masses,'materials':MATERIALS,'nmodes':80,'alpha':ALPHA,'beta':BETA,'table_connection_factor':10 if name=='foundation_table_stiff' else 1}
 (D/'metadata.json').write_text(json.dumps(metadata,indent=2))
 text=(out/'relative.inp').read_text()
 for index,dof in enumerate([1,2,3,1]):
  loads=[f'WORK, {dof}, 1.0',f'TOOL, {dof}, -1.0'] if index<3 else [f'{nut[4]}, 1, 1.0',f'{nut[3]}, 1, -1.0']
  text+='\n'.join(['*STEP','*STEADY STATE DYNAMICS','5, 3000, 24, 1.0','*MODAL DAMPING, RAYLEIGH',f',,{ALPHA},{BETA}','*CLOAD, OP=NEW']+loads+['*NODE PRINT, NSET=WORK','U','*NODE PRINT, NSET=TOOL','U','*NODE FILE, FREQUENCY=0','U','*END STEP'])+'\n'
 (out/'relative.inp').write_text(text);execute([Path.home()/'opt/ccx.sh','-i','relative'],out,'solve.log')
 print(name,m.nnode,'nodes',m.masses,flush=True)
print('VMC force-pair FRFs completed.',flush=True)
