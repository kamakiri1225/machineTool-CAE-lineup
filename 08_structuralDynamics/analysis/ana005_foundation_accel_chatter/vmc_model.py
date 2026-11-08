"""Simplified vertical machining center: cast frame/head, steel spindle,
20x20x100 mm carbide tool, cast-iron workpiece bonded to the table.
Separate cutting contact patches share coordinates but not node IDs.
"""
from pathlib import Path
import sys,copy,numpy as np
D=Path(__file__).resolve().parent;sys.path[:0]=[str(D.parent/'tools'),str(D.parent/'ana003_ccx_fistr_mpc_spring')]
import machine_model as mm,mpc_springs as ms
from foundation_model import build as foundation_build
HEAD=(380,620,340,750,650,850)
SPINDLE=(450,550,350,450,500,650)
import os
# 工具：VMC_TOOL=bar（既定）なら打撃試験用の太いテストバー（鋼 60×60×100 mm）、slender なら細い超硬工具（20×20×100 mm）
TOOL_KIND=os.environ.get('VMC_TOOL','bar')
TOOL=(470,530,370,430,400,500) if TOOL_KIND=='bar' else (490,510,390,410,400,500)
WORK=(450,550,350,450,350,400)
CONTACT=np.array([500.,400.,400.])
MATERIALS={'FRAME':{'E':120000.,'nu':.27,'rho':7.3e-9},'SPINDLE':{'E':210000.,'nu':.30,'rho':7.8e-9},'TOOL':({'E':210000.,'nu':.30,'rho':7.8e-9} if TOOL_KIND=='bar' else {'E':600000.,'nu':.22,'rho':14.5e-9})}
def material(center):
 x,y,z=center
 if TOOL[0]<x<TOOL[1] and TOOL[2]<y<TOOL[3] and 400<z<500:return 'TOOL'
 if 450<x<550 and 350<y<450 and 500<z<650:return 'SPINDLE'
 return 'FRAME'
def build(case):
 oldbed,oldtable=mm.BED,mm.TABLE
 try:
  mm.BED=oldbed+[HEAD,SPINDLE,TOOL];mm.TABLE=oldtable+[WORK]
  m=ms.build('rbe2') if case=='fixed' else foundation_build()
 finally:mm.BED,mm.TABLE=oldbed,oldtable
 if case=='foundation_table_stiff':m.springs=[(name,d,k*(10 if name.startswith(('CAR','NUT')) else 1),a,b) for name,d,k,a,b in m.springs]
 contacts={}
 for side,part in [('TOOL',m.parts[0]),('WORK',m.parts[1])]:
  ref,rot=m.nnode+1,m.nnode+2;m.nnode+=2
  m.extra_nodes.extend([(ref,*CONTACT),(rot,*CONTACT)])
  nodes=[c for c in part.coords if abs(c[3]-400)<1e-8 and TOOL[0]-1e-8<=c[1]<=TOOL[1]+1e-8 and TOOL[2]-1e-8<=c[2]<=TOOL[3]+1e-8]
  assert len(nodes)>=4,(side,len(nodes))
  for n,x,y,z in nodes:
   rx,ry,rz=np.array([x,y,z])-CONTACT
   m.equations.extend([[(n,1,1.),(ref,1,-1.),(rot,2,-rz),(rot,3,ry)],[(n,2,1.),(ref,2,-1.),(rot,3,-rx),(rot,1,rz)],[(n,3,1.),(ref,3,-1.),(rot,1,-ry),(rot,2,rx)]])
  contacts[side]={'ref':ref,'rot':rot,'nodes':[c[0] for c in nodes]}
 m.monitor=contacts['WORK']['ref'];m.contact=contacts;m.element_groups={name:[] for name in MATERIALS};m.element_parts=[];m.masses={'total_kg':0.,'moving_kg':0.}
 coord={n:np.array([x,y,z]) for p in m.parts for n,x,y,z in p.coords};eid=0
 for ip,part in enumerate(m.parts):
  for e in part.elems:
   eid+=1;pts=np.array([coord[n] for n in e[:8]]);center=pts.mean(axis=0);volume=np.prod(np.ptp(pts,axis=0));mat=material(center) if ip==0 else 'FRAME'
   if ip==0:m.element_groups[mat].append(eid)
   mass=volume*MATERIALS[mat]['rho']*1000;m.masses['total_kg']+=mass
   if ip==1:m.masses['moving_kg']+=mass
   tag=0 if ip==0 else 1
   if ip==0 and mat=='SPINDLE':tag=2
   if ip==0 and mat=='TOOL':tag=3
   if ip==1 and center[2]>350:tag=4
   m.element_parts.append(tag)
 return m

def write_input(m,path):
 mm.write_ccx(m,path,nmodes=80)
 text=Path(path).read_text();extra=[]
 for name,ids in m.element_groups.items():
  extra.append(f'*ELSET, ELSET={name}')
  for i in range(0,len(ids),12):extra.append(','.join(map(str,ids[i:i+12])))
 for name in ['SPINDLE','TOOL']:
  mat=MATERIALS[name];extra.extend([f'*MATERIAL, NAME={name}','*ELASTIC',f"{mat['E']}, {mat['nu']}",'*DENSITY',str(mat['rho'])])
 extra.extend(['*SOLID SECTION, ELSET=FRAME, MATERIAL=CASTIRON','*SOLID SECTION, ELSET=SPINDLE, MATERIAL=SPINDLE','*SOLID SECTION, ELSET=TOOL, MATERIAL=TOOL'])
 text=text.replace('*SOLID SECTION, ELSET=EBED, MATERIAL=CASTIRON','\n'.join(extra))
 text=text.replace('*BOUNDARY',f"*NSET, NSET=WORK\n{m.contact['WORK']['ref']}\n*NSET, NSET=TOOL\n{m.contact['TOOL']['ref']}\n*BOUNDARY",1)
 Path(path).write_text(text)
