"""各ケースの固有モードから刃先―ワーク相対 FRF と安定限界線図を作る → results/case_<名前>.npz, results/summary.json
CalculiX に 100 モード要求すると、固有値が負の偽モードが 1 本混ざることがある（12 モード要求では出ない）。
固有値 ≤ 0 のモードは捨て、その本数を summary に残す。"""
import sys, json, re
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; sys.path.insert(0, str(D)); R = D / "results"
import params as P
from lib.modal import contact_vectors
from lib.chatter import lobes, mode_share, G_of, RPM

def eigenvalues(dat):
    t = Path(dat).read_text(); blk = t[t.index("E I G E N V A L U E   O U T P U T"):]
    return np.array([float(m[1]) for m in re.finditer(r"^\s+\d+\s+(\S+)\s+\S+\s+\S+\s+\S+\s*$", blk.split("P A R T I C")[0], re.M)])

summary = json.loads((R / "summary.json").read_text()) if (R / "summary.json").exists() else {}
names = sys.argv[1:] or [n for n in P.CASES if (D / "cases" / n / "model.frd").exists()]
for name in names:
    cd = D / "cases" / name; ev = eigenvalues(cd / "model.dat")
    fr, d, sides = contact_vectors(cd); assert len(ev) == len(fr)
    ok = ev > 0; fr, d = fr[ok], d[:, ok]
    lim, fc, Gg = lobes(fr, d); k = int(np.nanargmin(lim)); sh = mode_share(fc[k], fr, d); top = np.argsort(sh)[::-1][:5]
    G0 = G_of([1e-3], fr, d)[0].real
    np.savez(R / f"case_{name}.npz", fr=fr, d=d, tool=sides["tool"][:, ok], work=sides["work"][:, ok], lim=lim, fc=fc, rpm=RPM,
             Gxx=Gg[:, 0, 0], Gyy=Gg[:, 1, 1], Gxy=Gg[:, 0, 1], mode_index=np.where(ok)[0] + 1)
    dy, dz = P.CASES[name]
    summary[name] = dict(dy=dy, dz=dz, dropped_nonpositive_modes=int((~ok).sum()), fit_residual_max={s: float(sides[s + '_fit_residual'][ok].max()) for s in ('tool', 'work')}, n_modes=int(ok.sum()), f_max_Hz=float(fr.max()),
                         amin_mm=float(lim[k]), rpm_at_min=float(RPM[k]), fc_at_min_Hz=float(fc[k]),
                         static_um_per_N={"xx": G0[0, 0] * 1e3, "yy": G0[1, 1] * 1e3, "zz": G0[2, 2] * 1e3},
                         top_modes=[dict(mode=int(np.where(ok)[0][r] + 1), f_Hz=float(fr[r]), share=float(sh[r])) for r in top])
    print(name, f"min a_lim {lim[k]:.2f} mm @ {RPM[k]:.0f} rpm fc {fc[k]:.0f} Hz dropped {(~ok).sum()}",
          [(s['mode'], round(s['f_Hz']), round(s['share'] * 100)) for s in summary[name]['top_modes'][:3]], flush=True)
(R / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
