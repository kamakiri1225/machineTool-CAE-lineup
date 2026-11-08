"""Reproducible exploratory analyses, with explicit reduced-model assumptions."""
from pathlib import Path
import sys,re,json
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.interpolate import interp1d
from scipy.optimize import brentq
D=Path(__file__).resolve().parent;R=D/'results';R.mkdir(exist_ok=True)
sys.path[:0]=[str(D.parent/'tools'),str(D.parent/'ana003_ccx_fistr_mpc_spring')]
from plot_modes import plt
import machine_model as mm
from foundation_model import K
from frequency_response import ALPHA,BETA

def read_set(path,name):
 lines=path.read_text().splitlines();rec=[]
 for i,line in enumerate(lines):
  m=re.match(rf'\s*displacements .* for set {name} and time\s+([\d.E+-]+)',line)
  if m:rec.append((float(m[1]),np.array(list(map(float,lines[i+2].split()[1:4])))))
 assert len(rec)%2==0
 steps=[];cur=[]
 for i in range(0,len(rec),2):
  f,u=rec[i];assert abs(f-rec[i+1][0])<1e-5
  if cur and f<cur[-1][0]:steps.append(cur);cur=[]
  cur.append((f,u+1j*rec[i+1][1]))
 steps.append(cur);assert len(steps)==4,len(steps)
 return [(np.array([x[0] for x in st]),np.array([x[1] for x in st])) for st in steps]

def acceleration():
 vol=lambda b:(b[1]-b[0])*(b[3]-b[2])*(b[5]-b[4])
 mass=json.loads((D/'metadata.json').read_text())['foundation']['mass']
 mt=mass['moving_kg'];M=mass['total_kg']
 k=4*K[1]*1000 # N/mm -> N/m
 c=ALPHA*M+BETA*k;fn=np.sqrt(k/M)/(2*np.pi);zeta=c/(2*np.sqrt(k*M))
 L,V=0.05,0.2;tol=1e-6;table=[];hist=[]
 for A in (1.,2.,4.):
  for tj in (0.,.02,.04):
   ta=V/A+tj;tc=L/V-ta;assert tc>=0;tm=2*ta+tc
   def accel(t):
    if t<0 or t>=tm:return 0.
    if t<ta:q=t;sgn=1.
    elif t<ta+tc:return 0.
    else:q=t-ta-tc;sgn=-1.
    if tj==0:return sgn*A
    if q<tj:return sgn*A*q/tj
    if q<ta-tj:return sgn*A
    return sgn*A*(ta-q)/tj
   # Piecewise integration at every ramp/plateau boundary avoids stepping over short pulses.
   breaks=sorted(set([0.,tj,ta-tj,ta,ta+tc,ta+tc+tj,tm-tj,tm,tm+1.2]))
   state=np.zeros(2);times=[];values=[]
   for a,b in zip(breaks[:-1],breaks[1:]):
    if b-a<1e-12:continue
    sol=solve_ivp(lambda t,y:[y[1],(-c*y[1]-k*y[0]-mt*accel(t))/M],(a,b),state,rtol=1e-9,atol=1e-12,max_step=.0002,dense_output=True)
    t=np.linspace(a,b,max(2,int((b-a)/.0002)+1),endpoint=False);times.extend(t);values.extend(sol.sol(t).T);state=sol.y[:,-1]
   t=np.array(times);y=np.array(values);a=np.array([accel(x) for x in t]);after=t>=tm
   late=t[after][abs(y[after,0])>tol];settle=0. if len(late)==0 else late[-1]-tm
   peak=float(abs(y[after,0]).max());distance=quad(lambda q:(tm-q)*accel(q),0,tm,points=[b for b in breaks if 0<b<tm],epsabs=1e-12)[0]
   final_velocity=quad(accel,0,tm,points=[b for b in breaks if 0<b<tm],epsabs=1e-12)[0]
   assert abs(distance-L)<1e-10 and abs(final_velocity)<1e-10,(A,tj,distance,final_velocity)
   base_acc=(-c*y[:,1]-k*y[:,0]-mt*a)/M
   thrust_peak=float(abs(mt*(a+base_acc)).max())
   table.append({'drive_force_peak_N':thrust_peak,'A_m_s2':A,'jerk_ramp_s':tj,'jerk_m_s3':None if tj==0 else A/tj,'move_s':tm,'base_residual_peak_um':peak*1e6,'base_settle_1um_s':settle,'move_plus_base_settle_s':tm+settle,'integrated_distance_m':distance})
   hist.append((A,tj,t,a,y[:,0],tm))
 fig,axs=plt.subplots(2,1,figsize=(11,7))
 for A,tj,t,a,x,tm in hist:
  if A!=2:continue
  label='加速度段差' if tj==0 else f'ジャーク立上げ {tj*1000:.0f} ms'
  axs[0].plot(t,a,label=f'{label}（移動 {tm:.2f}s）');axs[1].plot(t,x*1e6,label=label)
 axs[0].set_ylabel('相対送り指令の加速度 [m/s²]');axs[1].set_ylabel('基礎上の機械本体変位 [µm]');axs[1].set_xlabel('時刻 [s]')
 for ax in axs:ax.set_xlim(0,1);ax.grid(alpha=.3);ax.legend()
 axs[1].axhline(1,color='gray',ls=':');axs[1].axhline(-1,color='gray',ls=':')
 fig.suptitle('加減速反力の簡易1自由度モデル：同じ50 mm移動・最高速度0.2 m/s・最大加速度2 m/s²');fig.tight_layout();fig.savefig(R/'acceleration_base.png');plt.close(fig)
 (R/'acceleration.json').write_text(json.dumps({'M_kg':M,'moving_mass_kg':mt,'k_N_m':k,'c_Ns_m':c,'fn_Hz':fn,'zeta':zeta,'rows':table},indent=2))
 rows=['| 最大加速度 m/s² | ジャーク立上げ ms | 移動 s | 停止後の本体変位最大 µm | 本体変位±1 µm整定 s | 移動＋本体整定 s |','|---:|---:|---:|---:|---:|---:|']
 for v in table:rows.append(f"| {v['A_m_s2']:.0f} | {v['jerk_ramp_s']*1000:.0f} | {v['move_s']:.3f} | {v['base_residual_peak_um']:.3f} | {v['base_settle_1um_s']:.3f} | {v['move_plus_base_settle_s']:.3f} |")
 (R/'acceleration.md').write_text('\n'.join(rows)+'\n')
 return M,mt,fn,zeta

def boundary(interp,rpm,grid,Kc=1000.,teeth=4):
 result=[];critical=[]
 for n in rpm:
  tau=60/(teeth*n);h=interp(grid)*(1-np.exp(-2j*np.pi*grid*tau));idx=np.where(h.imag[:-1]*h.imag[1:]<0)[0];candidates=[]
  for i in idx:
   f=brentq(lambda q: (interp(q)*(1-np.exp(-2j*np.pi*q*tau))).imag,grid[i],grid[i+1],xtol=1e-9)
   v=interp(f)*(1-np.exp(-2j*np.pi*f*tau))
   if v.real < -1e-12:
    b=-1/(Kc*v.real)
    if b<1e6:candidates.append((b,f,abs(1+Kc*b*v)))
  if candidates:
   b,f,e=min(candidates);assert e<1e-5,e;result.append(b);critical.append(f)
  else:result.append(np.nan);critical.append(np.nan)
 return np.array(result),np.array(critical)

def main():
 M,mt,fn,zeta=acceleration()
 names={'fixed':'底面固定','foundation':'基礎ばね','foundation_table_stiff':'基礎ばね・テーブル結合10倍'}
 response={};fig,axs=plt.subplots(2,2,figsize=(12,8));summary={}
 for name in ['fixed','foundation','foundation_table_stiff']:
  work=read_set(D/name/'relative.dat','WORK');tool=read_set(D/name/'relative.dat','TOOL');response[name]=[]
  for j in range(3):
   f,w=work[j];ft,t=tool[j];assert np.array_equal(f,ft);g=w-t
   response[name].append((f,g[:,j]))
   np.savetxt(R/f'{name}_relative_{j+1}.csv',np.column_stack([f,g.real,g.imag]),delimiter=',',header='frequency_Hz,Hx_real_mm_N,Hy_real_mm_N,Hz_real_mm_N,Hx_imag_mm_N,Hy_imag_mm_N,Hz_imag_mm_N')
  matrix=np.stack([work[j][1]-tool[j][1] for j in range(3)],axis=2)
  reciprocity=float(np.max(abs(matrix-matrix.transpose(0,2,1)))/np.max(abs(matrix)))
  assert reciprocity<1e-3,reciprocity
  assert all((response[name][j][1][0].real>0) for j in range(3))
  summary[name+'_FRF_checks']={'reciprocity_relative_error':reciprocity}
  f,w=work[0];_,t=tool[0];g=w-t
  axs[0,0].plot(f,abs(g[:,0])*1000,label=f'{names[name]}：|Hrel,xx|')
  axs[1,0].plot(f,g[:,0].real*1000,label=f'{names[name]}：Re Hrel,xx')
  if name=='foundation':
   axs[0,1].plot(f,abs(w[:,0])*1000,label='ワーク側変位 |Uw/F|')
   axs[0,1].plot(f,abs(t[:,0])*1000,label='刃先変位 |Ut/F|')
   axs[0,1].plot(f,abs(g[:,0])*1000,label='相対変位 |(Uw−Ut)/F|')
   fd,wd=work[3];_,td=tool[3];gd=wd-td
   axs[1,1].plot(fd,abs(gd[:,0])*1000,label='ナット代表点の駆動力対 → 相対X変位')
 for ax in axs.flat:ax.set_xlabel('周波数 [Hz]');ax.set_ylabel('コンプライアンス [µm/N]');ax.grid(alpha=.3);ax.legend(fontsize=8);ax.set_xlim(5,3000)
 axs[0,0].set_title('逆向き切削力対 → 刃先・ワーク相対X');axs[0,1].set_title('基礎ばね：切削力対に対する2点の応答');axs[1,0].set_title('相対Xコンプライアンスの実部');axs[1,1].set_title('ナット駆動力対に対する相対X応答')
 fig.tight_layout();fig.savefig(R/'relative_frf.png');plt.close(fig)
 rpm=np.linspace(500,16000,600);grid=np.linspace(5,3000,24000);fig,axs=plt.subplots(2,1,figsize=(11,9),sharex=True)
 for direction in (0,1):
  for name in ['fixed','foundation','foundation_table_stiff']:
   f,g=response[name][direction];interp=interp1d(f,g,kind='linear',bounds_error=True)
   limit,fc=boundary(interp,rpm,grid)
   # Root-grid refinement at a subset of spindle speeds, independently of plotting density.
   ids=np.linspace(0,len(rpm)-1,12,dtype=int);refined,_=boundary(interp,rpm[ids],np.linspace(5,3000,48000));valid=np.isfinite(refined)&np.isfinite(limit[ids]);err=float(np.max(abs(refined[valid]-limit[ids][valid])/refined[valid])) if valid.any() else 0.
   assert err<.02,err
   np.savetxt(R/f'lobe_{name}_{direction+1}.csv',np.column_stack([rpm,limit,fc]),delimiter=',',header='rpm,boundary_equivalent_width_mm,critical_chatter_frequency_Hz')
   axs[direction].plot(rpm,limit,label=f'{names[name]}：相対{"XY"[direction]}変位')
   summary[f'{name}_{direction+1}']={'root_grid_relative_error':err,'minimum_width_mm':float(np.nanmin(limit)) if np.isfinite(limit).any() else None,'frequency_band_Hz':[5,3000]}
  axs[direction].set_ylim(bottom=0);axs[direction].grid(alpha=.3);axs[direction].legend();axs[direction].set_ylabel('等価切削幅 b [mm]');axs[direction].set_title(f'{"XY"[direction]}方向：1方向・時間不変再生モデル（フライスの実機予測ではない）')
 axs[-1].set_xlabel('仮定した4刃の回転数 [min⁻¹]');fig.suptitle('探索的な安定境界：Kc=1000 N/mm²、遅れτ=60/(4n)、刃先、5〜3000 Hz帯域内');fig.tight_layout();fig.savefig(R/'simplified_lobes.png');plt.close(fig)
 # Independent single-DOF analytical consistency check of the characteristic equation.
 k,zn,wn,Kc=10000.,.02,2*np.pi*100,1000.;w=wn*np.sqrt(1+2*zn);g=1/(k*(1-(w/wn)**2+2j*zn*w/wn));b=-1/(2*Kc*g.real);expected=2*k*zn*(1+zn)/Kc
 assert np.isclose(b,expected,rtol=1e-10)
 e=1+1/(Kc*b*g);tau=((-np.angle(e))%(2*np.pi)+4*np.pi)/w;res=abs(1+Kc*b*g*(1-np.exp(-1j*w*tau)));assert res<1e-10
 summary['single_DOF_check']={'b_mm':b,'analytical_b_mm':expected,'characteristic_residual':res}
 (R/'checks.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'M':M,'mt':mt,'fn':fn,'zeta':zeta,'chatter':summary},indent=2))
if __name__=='__main__':main()
