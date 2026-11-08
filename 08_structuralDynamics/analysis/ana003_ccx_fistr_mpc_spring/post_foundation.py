"""Verify and plot foundation-spring case against the fixed-bottom reference."""
from pathlib import Path
import json,re,sys,numpy as np
from post_frequency import ccx_dat,read_front
from foundation_model import build
from plot_modes import plt
D=Path(__file__).resolve().parent;R=D/'foundation/results';m=build();checks=[]
cc=ccx_dat(D/'foundation/ccx/response.dat')
selected=json.loads((R/'selected.json').read_text())
fig,axs=plt.subplots(3,1,figsize=(11,10))
for i,(f,u) in enumerate(cc):
 dof=i+1;norm=np.linalg.norm(abs(u),axis=1)*1000
 np.savetxt(R/f'rbe2_{dof}_ccx.csv',np.column_stack([f,u.real,u.imag,norm]),delimiter=',',header='frequency_Hz,Ux_real_mm,Uy_real_mm,Uz_real_mm,Ux_imag_mm,Uy_imag_mm,Uz_imag_mm,norm_um_per_N')
 front=np.loadtxt(R/f'rbe2_{dof}_fistr_components.csv',delimiter=',')
 shape=read_front(D/f'foundation/fistr/shape_{dof}_dyna.res.0.1');amp=abs(shape[m.monitor])*1000;j=np.argmin(abs(front[:,0]-selected[i]));assert np.allclose(amp,front[j,1:],rtol=2e-4,atol=1e-10)
 # Evaluate MPC residuals on the solver's selected-frequency complex response.
 residuals=[abs(sum(a*shape[n][d-1] for n,d,a in eq)) for eq in m.equations]
 relative=max(residuals)/max(np.linalg.norm(abs(v)) for n,v in shape.items());assert relative<1e-4,relative
 fixed=np.loadtxt(D/f'results/frf/rbe2_{dof}_fistr_components.csv',delimiter=',')
 ax=axs[i];ax.plot(fixed[:,0],fixed[:,dof],label='底面固定（FrontISTR、50 Hz以上）',color='#7f8c8d');ax.plot(f,abs(u[:,i])*1000,label='基礎ばね（CalculiX）',color='#24659c');ax.plot(front[:,0],front[:,dof],'--',label='基礎ばね（FrontISTR）',color='#d17423')
 ax.set_title(f'同じMON：+{"XYZ"[i]}方向1 N加振 → {"XYZ"[i]}方向変位',loc='left');ax.set_ylabel('コンプライアンス [µm/N]');ax.set_ylim(bottom=0);ax.grid(alpha=.3);ax.legend();ax.set_xlim(5,750)
 checks.append({'dof':dof,'selected_frequency_Hz':selected[i],'max_normalized_MPC_residual':relative,'selected_components_um_per_N':amp.tolist()})
axs[-1].set_xlabel('加振周波数 [Hz]');fig.tight_layout();fig.savefig(R/'fixed_vs_foundation.png');plt.close(fig)
(R/'checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
