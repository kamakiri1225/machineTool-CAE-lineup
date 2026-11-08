from pathlib import Path
import numpy as np
from post_frequency import selected_frd
D=Path(__file__).resolve().parent
for i,(f,data) in enumerate(selected_frd(D/'frf/rbe2/ccx/response.frd',[183.,454.,625.]),1):
 ids=np.array(sorted(data));np.savez_compressed(D/f'results/frf/ccx_shape_{i}.npz',frequency=f,ids=ids,U=np.array([data[n] for n in ids]))
 print('cached',i,f,flush=True)
