#!/usr/bin/env python3
"""Plot actual sampled nozzle temperatures, including incomplete physical-time runs."""
from pathlib import Path
import sys,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
case=Path(sys.argv[1]);out=case/'images';out.mkdir(exist_ok=True)
rows={}
for p in (case/'postProcessing/nozzle/nozzleProbes').glob('*/T'):
    for line in p.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        vals=list(map(float,line.split()))
        if len(vals)==5:rows[vals[0]]=vals
if not rows:raise SystemExit('No sampled data')
a=np.array([rows[t] for t in sorted(rows)])
with (out/'nozzle_temperature_history.csv').open('w') as f:
    writer=csv.writer(f);writer.writerow(['time_s','P1_T_K','P2_T_K','P3_T_K','P4_T_K','P1_deltaT_K','P2_deltaT_K','P3_deltaT_K','P4_deltaT_K'])
    for row in a:writer.writerow(list(row)+list(row[1:]-300))
plt.rcParams.update({'font.family':['Noto Sans CJK JP','DejaVu Sans'],'font.size':21})
fig,ax=plt.subplots(figsize=(12,8))
names=['P1 先端','P2 内壁側中央','P3 外壁側中央','P4 取付側']
for i,name in enumerate(names):ax.plot(a[:,0]*1e3,a[:,i+1]-300,'o-',lw=2.5,label=name)
ax.set_xlabel('時刻 [ms]',fontsize=25);ax.set_ylabel('初期300 Kからの温度変化 ΔT [K]',fontsize=25)
ax.set_title(f'ノズル温度の時刻歴（計算済み区間：0〜{a[-1,0]*1e3:.5g} ms）',fontsize=24,pad=18)
ax.legend(fontsize=20);ax.grid(alpha=.25);ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,3))
fig.text(.12,.02,'予備解析：銅を仮定・取付面300 K・放射／直接レーザー吸収なし。50 ms未完了。',fontsize=14)
fig.tight_layout(rect=[0,.06,1,1]);fig.savefig(out/'03_ノズル温度時刻歴_計算済み区間.png',dpi=140)
print('Latest sample s:',a[-1,0],'; deltaT K:',a[-1,1:]-300,flush=True)
