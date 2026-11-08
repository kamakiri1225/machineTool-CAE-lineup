"""Export FE modes and independently integrated consistent Hex20 modal masses."""
from pathlib import Path
import json,numpy as np
from vmc_model import build,MATERIALS,material
D=Path(__file__).resolve().parent;R=D/'results';import sys as _s;CASE=_s.argv[1] if len(_s.argv)>1 else 'foundation'
m=build(CASE);nm=80
phi=np.zeros((m.nnode,3,nm));freq=[]
with (D/CASE/'relative.frd').open() as stream:
 for line in stream:
  if not line.startswith('  100CL'):continue
  f=float(line[12:24]);name=next(stream)[5:13].strip()
  if name!='DISP':continue
  j=len(freq);assert j<nm;freq.append(f)
  for row in stream:
   if row.startswith(' -1'):break
  while row.startswith(' -1'):
   n=int(row[3:13]);phi[n-1,:,j]=[float(row[13:25]),float(row[25:37]),float(row[37:49])];row=next(stream)
assert len(freq)==nm
# CalculiX FRD auxiliary representative nodes can be zero even when DAT is valid.
# Recover contact translations/rotations from the actual MPC-constrained face nodes.
coordinates={n:np.array([x,y,z]) for part in m.parts for n,x,y,z in part.coords}
for side,patch in m.contact.items():
 A=[]
 for n in patch['nodes']:
  rx,ry,rz=coordinates[n]-np.array([500.,400.,400.])
  A.extend([[1,0,0,0,rz,-ry],[0,1,0,-rz,0,rx],[0,0,1,ry,-rx,0]])
 A=np.array(A);nodal=phi[np.array(patch['nodes'])-1].reshape(-1,nm);solution=np.linalg.lstsq(A,nodal,rcond=None)[0]
 assert np.linalg.norm(A@solution-nodal)/np.linalg.norm(nodal)<1e-5
 phi[patch['ref']-1]=solution[:3];phi[patch['rot']-1]=solution[3:]

corners=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
edges=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)];mid=np.array([(corners[a]+corners[b])/2 for a,b in edges])
def shape(q):
 N=[np.prod(1+s*q)*(s@q-2)/8 for s in corners]
 for s in mid:
  axis=np.where(s==0)[0][0];others=[a for a in range(3) if a!=axis];N.append((1-q[axis]**2)*np.prod(1+s[others]*q[others])/4)
 return np.array(N)
x,w=np.polynomial.legendre.leggauss(3);C=np.zeros((20,20))
for i in range(3):
 for j in range(3):
  for k in range(3):
   N=shape(np.array([x[i],x[j],x[k]]));assert np.isclose(N.sum(),1);C+=np.outer(N,N)*w[i]*w[j]*w[k]/8
assert np.isclose(C.sum(),1) and np.linalg.eigvalsh(C).min()>0
coords={n:np.array([x,y,z]) for part in m.parts for n,x,y,z in part.coords};idx=[];mass=[]
for ip,part in enumerate(m.parts):
 for e in part.elems:
  points=np.array([coords[n] for n in e[:8]]);vol=np.prod(np.ptp(points,axis=0));mat=material(points.mean(axis=0)) if ip==0 else 'FRAME';idx.append(np.array(e)-1);mass.append(vol*MATERIALS[mat]['rho'])
idx=np.array(idx);mass=np.array(mass);assert np.isclose(mass.sum()*1000,m.masses['total_kg'])
modal=np.zeros(nm)
for start in range(0,len(idx),128):
 q=phi[idx[start:start+128]];mq=np.einsum('ij,ejcm->eicm',C,q,optimize=True);modal+=np.einsum('eicm,eicm,e->m',q,mq,mass[start:start+128],optimize=True)
assert np.all(modal>0)
metadata=json.loads((D/'metadata.json').read_text())[CASE]
np.savez_compressed(R/f'{CASE}_modal_basis.npz',frequency_hz=np.array(freq),modes=phi,modal_mass=modal,work_index=m.contact['WORK']['ref']-1,tool_index=m.contact['TOOL']['ref']-1,alpha=metadata['alpha'],beta=metadata['beta'],node_ids=np.arange(1,m.nnode+1))
(R/(f'modal_mass_checks.json' if CASE=='foundation' else f'modal_mass_checks_{CASE}.json')).write_text(json.dumps({'consistent_solid_mass_kg':float(mass.sum()*1000),'modal_mass_min':float(modal.min()),'modal_mass_max':float(modal.max()),'nmodes':nm,'integration':'Hex20 consistent mass, 3x3x3 Gauss, density per material'},indent=2))
print('Exported',nm,'modes; modal mass range',modal.min(),modal.max(),flush=True)
