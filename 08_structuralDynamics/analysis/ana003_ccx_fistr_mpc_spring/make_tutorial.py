"""Self-contained 2-block RBE2-like MPC + reference-point spring example.
No dependency on the machine model generator. Generates readable solver decks.
"""
from pathlib import Path
import subprocess,os
D=Path(__file__).resolve().parent/'tutorial';D.mkdir(exist_ok=True)
points={1:(0,0,0),2:(100,0,0),3:(100,100,0),4:(0,100,0),5:(0,0,100),6:(100,0,100),7:(100,100,100),8:(0,100,100),9:(0,0,100),10:(100,0,100),11:(100,100,100),12:(0,100,100),13:(0,0,200),14:(100,0,200),15:(100,100,200),16:(0,100,200),101:(50,50,100),102:(50,50,100),201:(50,50,100),202:(50,50,100)}
faces={'FACE_B':[5,6,7,8],'FACE_T':[9,10,11,12]}
eqs=[]
for face,ref,rot in [('FACE_B',101,102),('FACE_T',201,202)]:
 for n in faces[face]:
  x,y,z=points[n];rx=x-50;ry=y-50
  eqs.extend([[(n,1,1),(ref,1,-1),(rot,3,ry)],[(n,2,1),(ref,2,-1),(rot,3,-rx)],[(n,3,1),(ref,3,-1),(rot,1,-ry),(rot,2,rx)]])
springs=[(100+i,101,201,i,1000.) for i in (1,2,3)]+[(200+i,102,202,i,1e6) for i in (1,2,3)]
# Files contain the face sets explicitly even though equations list node IDs individually.
c=['*HEADING','two separate blocks; face sets -> MPC -> reference springs','*NODE']+[f'{n}, {x}, {y}, {z}' for n,(x,y,z) in points.items()]
c+=['*ELEMENT, TYPE=C3D8, ELSET=SOLIDS','1, 1,2,3,4,5,6,7,8','2, 9,10,11,12,13,14,15,16']
f=['!HEADER','two separate blocks; face sets -> MPC -> reference springs','!NODE']+[f'{n}, {x}, {y}, {z}' for n,(x,y,z) in points.items()]
f+=['!ELEMENT, TYPE=361, EGRP=SOLIDS','1, 1,2,3,4,5,6,7,8','2, 9,10,11,12,13,14,15,16']
for name,ns in {**faces,'FIX':[1,2,3,4],'MON':[101,102,201,202,13]}.items():
 c +=[f'*NSET, NSET={name}',','.join(map(str,ns))]
 f +=[f'!NGROUP, NGRP={name}',','.join(map(str,ns))]
for eq in eqs:
 row=', '.join(f'{n}, {d}, {a}' for n,d,a in eq)
 c+=['*EQUATION',str(len(eq)),row]
f+=['!EQUATION']
for eq in eqs:f +=[f'{len(eq)}, 0.0',', '.join(f'{n}, {d}, {a}' for n,d,a in eq)]
mats=[]
for eid,a,b,d,k in springs:
 name=f'SP{eid}'
 c +=[f'*ELEMENT, TYPE=SPRING2, ELSET={name}',f'{eid}, {a}, {b}',f'*SPRING, ELSET={name}',f'{d}, {d}',str(k)]
 f +=[f'!ELEMENT, TYPE=511, EGRP={name}',f'{eid}, {a}, {b}',f'!MATERIAL, NAME={name}, ITEM=1','!ITEM=1, SUBITEM=2','1.0, 0.3',f'!SECTION, TYPE=INTERFACE, EGRP={name}, MATERIAL={name}','1.0']
 mats +=[f'!MATERIAL, NAME={name}','!SPRING_D',f'{d}, {d}, {k}']
c+=['*MATERIAL, NAME=STEEL','*ELASTIC','200000, 0.3','*DENSITY','7.8e-9','*SOLID SECTION, ELSET=SOLIDS, MATERIAL=STEEL','*BOUNDARY','FIX, 1, 3','*STEP','*STATIC','*CLOAD','13, 2, 1.0','*NODE PRINT, NSET=MON','U','*NODE FILE','U','*END STEP']
f+=['!MATERIAL, NAME=STEEL, ITEM=2','!ITEM=1, SUBITEM=2','200000, 0.3','!ITEM=2, SUBITEM=1','7.8e-9','!SECTION, TYPE=SOLID, EGRP=SOLIDS, MATERIAL=STEEL','1.0','!END']
for name in ['ccx','fistr']:(D/name).mkdir(exist_ok=True)
(D/'ccx/model.inp').write_text('\n'.join(c)+'\n');(D/'fistr/model.msh').write_text('\n'.join(f)+'\n')
(D/'fistr/model.cnt').write_text('\n'.join(['!VERSION','3','!SOLUTION, TYPE=STATIC','!BOUNDARY','FIX, 1, 3, 0.0','!CLOAD','13, 2, 1.0']+mats+['!SOLVER, METHOD=DIRECT','!WRITE,RESULT','!END'])+'\n')
(D/'fistr/hecmw_ctrl.dat').write_text('!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE\nmodel.msh\n!CONTROL, NAME=fstrCNT\nmodel.cnt\n!RESULT, NAME=fstrRES, IO=OUT\nmodel.res\n')
env=os.environ.copy();lib=Path.home()/'opt/ccxroot/usr/lib/x86_64-linux-gnu';env['LD_LIBRARY_PATH']=':'.join(map(str,[lib,lib/'blas',lib/'lapack']));env['OMP_NUM_THREADS']='2'
for name,exe,args in [('ccx',Path.home()/'opt/ccx.sh',['-i','model']),('fistr',Path.home()/'opt/frontistr/bin/fistr1',[])]:
 with (D/name/'solve.log').open('w') as out:subprocess.run([str(exe)]+args,cwd=D/name,env=env,stdout=out,stderr=subprocess.STDOUT,check=True)
print('Generated and solved:',D,flush=True)
# Verify physical spring extension, not the absolute displacement of either solid.
import re
text=(D/'ccx/model.dat').read_text()
cc={}  # populated below from the MON print block
for ln in text.splitlines():
 p=ln.split()
 if len(p)==4 and p[0] in ('101','102','201','202','13'): cc[int(p[0])]=[float(v) for v in p[1:]]
lines=(D/'fistr/model.res.0.1').read_text().splitlines();i=lines.index('*data');nn=int(lines[i+1].split()[0]);nc=int(lines[i+2].split()[0]);nd=sum(map(int,lines[i+3].split()));i+=4+nc
front={}
for _ in range(nn):
 n=int(lines[i]);i+=1;v=[]
 while len(v)<nd:v.extend(map(float,lines[i].split()));i+=1
 front[n]=v[:3]
rows=['| ソルバ | 下側101のY変位 mm | 上側201のY変位 mm | 差 mm | 手計算 F/k mm |','|---|---:|---:|---:|---:|']
for name,lo,hi in [('CalculiX',cc[101][1],cc[201][1]),('FrontISTR',front[101][1],front[201][1])]:
 delta=hi-lo
 assert abs(delta-.001)<1e-9,(name,delta)
 rows.append(f'| {name} | {lo:.9g} | {hi:.9g} | {delta:.9g} | 0.001 |')
# 回転の検算：荷重点 13 は代表点から r=(-50,-50,100)。モーメント M = r x F = (-100, 0, -50) N·mm、回転ばね 1e6 N·mm/rad
rows+=['','| ソルバ | θx(202)−θx(102) rad | 手計算 Mx/kθ | θz(202)−θz(102) rad | 手計算 Mz/kθ | 荷重点 13 の Y 変位 mm |','|---|---:|---:|---:|---:|---:|']
for name,r1,r2,u13 in [('CalculiX',cc[102],cc[202],cc[13][1]),('FrontISTR',front[102],front[202],front[13][1])]:
 tx=r2[0]-r1[0];tz=r2[2]-r1[2]
 assert abs(tx+1e-4)<1e-9 and abs(tz+5e-5)<1e-9,(name,tx,tz)
 rows.append(f'| {name} | {tx:.6g} | -0.0001 | {tz:.6g} | -0.00005 | {u13:.6g} |')
(D/'check.md').write_text('\n'.join(rows)+'\n');print('\n'.join(rows))
