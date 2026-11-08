"""全ケースのモデルを作り（params.CASES）、CalculiX で固有値解析を回す。
使い方:  ~/opt/venv/bin/python 01_build_and_run.py [ケース名 ...]   （省略で全ケース。--build-only でモデルだけ作る）"""
import sys, json, subprocess, time, os
from pathlib import Path
D = Path(__file__).resolve().parent; sys.path.insert(0, str(D))
import params as P
from lib.fem import build, write_ccx

args = [a for a in sys.argv[1:] if not a.startswith("--")]
names = args or list(P.CASES)
for name in names:
    dy, dz = P.CASES[name]; out = D / "cases" / name; out.mkdir(parents=True, exist_ok=True)
    m = build(dy, dz); write_ccx(m, out / "model.inp")
    (out / "model_info.json").write_text(json.dumps(m.info, indent=1, ensure_ascii=False))
    print(f"{name}: dy={dy} dz={dz} nodes={m.nnode} elems={len(m.elems)} springs={len(m.springs)} mass={ {k: round(v) for k, v in m.info['mass_kg'].items()} }", flush=True)
    if "--build-only" in sys.argv:
        continue
    t0 = time.time(); env = dict(os.environ, OMP_NUM_THREADS="16", CCX_NPROC_EQUATION_SOLVER="16")
    r = subprocess.run([str(Path.home() / "opt/ccx.sh"), "-i", "model"], cwd=out, env=env,
                       stdout=open(out / "ccx.log", "w"), stderr=subprocess.STDOUT)
    print(f"   ccx rc={r.returncode} {time.time() - t0:.0f} s", flush=True)
    assert r.returncode == 0
