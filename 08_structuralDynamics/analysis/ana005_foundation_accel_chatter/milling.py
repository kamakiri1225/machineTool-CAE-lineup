"""Two-direction zero-order full-slot milling stability using complex relative FRF.
Coordinate convention defined through B(phi), not inferred from unsigned amplitudes.
"""
from pathlib import Path
import json,sys,numpy as np
from scipy.interpolate import interp1d
from analyze import boundary
D=Path(__file__).resolve().parent;R=D/'results';sys.path.insert(0,str(D.parent/'tools'))
from plot_modes import plt
TEETH=4;KT=1000.;KR=.3
# Full-slot phi=0..pi: B0=(Z/(2pi))*integral B(phi)dphi.
B0=TEETH/4*np.array([[KR,1.],[-1.,KR]])
phi=np.linspace(0,np.pi,10001);sn,cs=np.sin(phi),np.cos(phi)
B=np.stack([sn*cs+KR*sn**2,cs**2+KR*sn*cs,-sn**2+KR*sn*cs,-sn*cs+KR*cs**2],axis=1).reshape(-1,2,2)
assert np.allclose(TEETH/(2*np.pi)*np.trapezoid(B,phi,axis=0),B0,rtol=1e-8,atol=1e-8)
labels={'fixed':'底面固定','foundation':'基礎ばね','foundation_table_stiff':'基礎ばね・ガイド/ナット剛性10倍'}
rpm=np.linspace(500,16000,600);grid=np.linspace(5,3000,24000);summary={};fig,axs=plt.subplots(2,1,figsize=(11,9),sharex=True)
for name in labels:
 cols=[np.loadtxt(R/f'{name}_relative_{j}.csv',delimiter=',') for j in (1,2)];f=cols[0][:,0];assert np.array_equal(f,cols[1][:,0])
 H=np.stack([v[:,1:3]+1j*v[:,4:6] for v in cols],axis=2);G=interp1d(f,H,axis=0,kind='linear',bounds_error=True)
 def eigen_interps(sample):
  lam=np.linalg.eigvals(G(sample)@B0)
  for i in range(1,len(lam)):
   if np.sum(abs(lam[i]-lam[i-1]))>np.sum(abs(lam[i,::-1]-lam[i-1])):lam[i]=lam[i,::-1]
  return [interp1d(sample,lam[:,j],kind='linear',bounds_error=True) for j in (0,1)]
 interps=eigen_interps(grid);branch=[boundary(e,rpm,grid,Kc=KT,teeth=TEETH) for e in interps]
 values=np.stack([b[0] for b in branch],axis=1);freqs=np.stack([b[1] for b in branch],axis=1);which=np.argmin(np.where(np.isfinite(values),values,np.inf),axis=1);limit=values[np.arange(len(rpm)),which];critical=freqs[np.arange(len(rpm)),which]
 maxres=0.
 for n,ap,fc in zip(rpm,limit,critical):
  if not np.isfinite(ap):continue
  characteristic=np.eye(2)+KT*ap*(G(fc)@B0)*(1-np.exp(-2j*np.pi*fc*60/(TEETH*n)))
  # Smallest singular value detects a singular characteristic matrix robustly.
  residual=np.linalg.svd(characteristic,compute_uv=False)[-1]/max(np.linalg.norm(characteristic),1.)
  maxres=max(maxres,float(residual))
 assert maxres<2e-3,maxres
 # Refine eigenfrequency interpolation and root sampling at 8 independent spindle speeds.
 ids=np.linspace(0,len(rpm)-1,8,dtype=int);fine=np.linspace(5,3000,48000);fine_e=eigen_interps(fine);vv=np.stack([boundary(e,rpm[ids],fine,Kc=KT,teeth=TEETH)[0] for e in fine_e],axis=1);ref=np.min(np.where(np.isfinite(vv),vv,np.inf),axis=1);valid=np.isfinite(ref)&np.isfinite(limit[ids]);error=float(np.max(abs(ref[valid]-limit[ids][valid])/ref[valid])) if valid.any() else 0.;assert error<.02,error
 np.savetxt(R/f'milling_{name}.csv',np.column_stack([rpm,limit,critical]),delimiter=',',header='rpm,axial_depth_limit_mm,chatter_frequency_Hz')
 axs[0].plot(rpm,limit,label=labels[name]);axs[1].plot(rpm,critical,label=labels[name])
 finite=limit[np.isfinite(limit)];summary[name]={'minimum_axial_depth_mm':float(finite.min()) if len(finite) else None,'maximum_characteristic_relative_residual':maxres,'frequency_grid_refinement_relative_error':error,'finite_points':len(finite)}
 for n in (3000,6000,10000):
  j=np.argmin(abs(rpm-n));summary[name][f'at_{n}_rpm']={'evaluated_rpm':float(rpm[j]),'axial_depth_mm':float(limit[j]) if np.isfinite(limit[j]) else None,'chatter_frequency_Hz':float(critical[j]) if np.isfinite(critical[j]) else None}
axs[0].set_ylabel('軸方向切込み限界 ap [mm]');axs[0].set_ylim(bottom=0);axs[1].set_ylabel('境界のびびり周波数 [Hz]');axs[1].set_xlabel('主軸回転数 [min⁻¹]')
for ax in axs:ax.grid(alpha=.3);ax.legend()
fig.suptitle('刃先―ワーク相対コンプライアンス（X・Y とクロス成分）による安定限界\n4刃、Kt=1000 N/mm²、Kr=0.3、5〜3000 Hz帯域、切削試験未実施');fig.tight_layout();fig.savefig(R/'milling_lobes.png');plt.close(fig)
(R/'milling_checks.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print(json.dumps(summary,indent=2))
