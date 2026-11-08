from pathlib import Path
import sys,numpy as np
D=Path(__file__).resolve().parent;R=D/'results';sys.path.insert(0,str(D.parent/'tools'))
from plot_modes import plt
fig,axs=plt.subplots(3,3,figsize=(16,11),sharex=True)
labels={'fixed':'底面固定','foundation':'基礎ばね','foundation_table_stiff':'テーブル結合10倍'}
for name,label in labels.items():
 for j in range(3):
  a=np.loadtxt(R/f'{name}_relative_{j+1}.csv',delimiter=',');u=a[:,1:4]+1j*a[:,4:7]
  for i in range(3):
   ax=axs[i,j];ax.plot(a[:,0],abs(u[:,i])*1000,label=label,lw=1.2);ax.set_title(f'G{"xyz"[i]}{"xyz"[j]}：{"XYZ"[j]}切削力対 → 相対{"XYZ"[i]}変位',fontsize=10);ax.set_ylim(bottom=0);ax.set_xlim(5,3000);ax.grid(alpha=.3)
   if i==2:ax.set_xlabel('周波数 [Hz]')
   if j==0:ax.set_ylabel('相対コンプライアンス [µm/N]')
for ax in axs.flat:ax.legend(fontsize=7)
fig.suptitle('立型MC：切削点の刃先―ワーク相対FRF 9成分\n行＝観測方向、列＝切削力対の方向（同じ位置、逆向きに各1 N）',fontsize=14);fig.tight_layout();fig.savefig(R/'relative_matrix.png');plt.close(fig)
