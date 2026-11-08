"""All nine H_ij: output i, force input j, at the same MON point."""
from pathlib import Path
import sys,numpy as np
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent/'tools'))
from plot_modes import plt
foundation='--foundation' in sys.argv
R=D/'foundation/results' if foundation else D/'results/frf'
fig,axs=plt.subplots(3,3,figsize=(16,11),sharex=True)
rows=['| 加振点・観測点 | 入力 | 出力 | 記号（出力→入力の添字順） | FrontISTR最大振幅 µm/N |','|---|---|---|---|---:|']
for j in range(3):
 cc=np.loadtxt(R/f'rbe2_{j+1}_ccx.csv',delimiter=',');front=np.loadtxt(R/f'rbe2_{j+1}_fistr_components.csv',delimiter=',')
 U=cc[:,1:4]+1j*cc[:,4:7]
 for i in range(3):
  ax=axs[i,j];c=abs(U[:,i])*1000;f=front[:,i+1];name=f'H{"xyz"[i]}{"xyz"[j]}'
  ax.plot(cc[:,0],c,label='CalculiX',lw=1.5);ax.plot(front[:,0],f,'--',label='FrontISTR',lw=1.3)
  ax.set_title(f'{name}：+{"XYZ"[j]}方向加振 → {"XYZ"[i]}方向変位',fontsize=11)
  ax.set_ylim(0,max(c.max(),f.max())*1.12);ax.grid(alpha=.25,which='both');ax.set_xlim(5 if foundation else 50,750)
  if j==0:ax.set_ylabel('振幅 [µm/N]')
  if i==2:ax.set_xlabel('加振周波数 [Hz]')
  ax.legend(fontsize=8)
  rows.append(f'| MON (500,300,350) mm | +{"XYZ"[j]}方向1 N | {"XYZ"[i]}変位振幅 | {name} | {f.max():.6g} |')
fig.suptitle(('基礎ばね支持・' if foundation else '')+'RBE2：テーブル上面中央MONの9成分コンプライアンス\n行＝観測変位方向、列＝加振方向／対角＝同方向、非対角＝クロス方向',fontsize=15)
fig.tight_layout();fig.savefig(R/'component_frf_3x3.png');plt.close(fig)
(R/'component_summary.md').write_text('\n'.join(rows)+'\n')
print('\n'.join(rows))
