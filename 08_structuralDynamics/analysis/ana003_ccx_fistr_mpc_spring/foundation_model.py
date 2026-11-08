"""Four MPC-coupled foundation patches replacing the all-bottom fixed BC."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'tools'))
import mpc_springs as ms
PATCHES=[(50,150,50,150),(850,950,50,150),(50,150,750,850),(850,950,750,850)]
K={1:1e4,2:1e4,3:1e5}  # N/mm per foot; illustrative values, not identified.
def build():
 m=ms.build('rbe2');n=m.nnode+1;m.fixed=[];m.foundations=[]
 for i,(x0,x1,y0,y1) in enumerate(PATCHES,1):
  xyz=np.array([(x0+x1)/2,(y0+y1)/2,0.]);ref,rot,ground,grot=n,n+1,n+2,n+3;n+=4
  m.extra_nodes.extend([(ref,*xyz),(rot,*xyz),(ground,*xyz),(grot,*xyz)])
  nodes=[(j,x,y) for j,x,y,z in m.parts[0].coords if abs(z)<1e-8 and x0-1e-8<=x<=x1+1e-8 and y0-1e-8<=y<=y1+1e-8]
  assert len(nodes)>=4
  for j,x,y in nodes:
   rx,ry=x-xyz[0],y-xyz[1]
   m.equations.extend([[(j,1,1.),(ref,1,-1.),(rot,3,ry)],[(j,2,1.),(ref,2,-1.),(rot,3,-rx)],[(j,3,1.),(ref,3,-1.),(rot,1,-ry),(rot,2,rx)]])
  for d,k in K.items():m.springs.append((f'FOOT{i}',d,k,ref,ground))
  # Uniform 100x100mm patch rotational stiffness about its centroid.
  kr={1:K[3]*100**2/12,2:K[3]*100**2/12,3:(K[1]+K[2])*100**2/12}
  for d,k in kr.items():m.springs.append((f'FOOT{i}_ROT',d,k,rot,grot))
  m.fixed.extend([ground,grot]);m.foundations.append({'name':f'FOOT{i}','ref':ref,'rot':rot,'ground':ground,'ground_rot':grot,'kr_Nmm_per_rad':kr,'xyz':xyz.tolist(),'surface_nodes':[j for j,x,y in nodes],'k_N_per_mm':K})
 m.nnode=n-1
 return m
