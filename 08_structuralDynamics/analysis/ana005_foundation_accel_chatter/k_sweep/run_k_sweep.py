"""ana005 4.2c：脚の横剛性を 1/4・1・4 倍にした全体 FEM（CalculiX）で、本当に固有振動数が半分・2 倍になるか確かめる
   ・固有値（20 モード）と、加減速反力（ナットのばねの両端に +F・−F、X 方向）に対する周波数応答（1〜300 Hz）
   ・観測点：ベッド側面、コラム上端、刃先（TOOL）、ワーク（WORK）"""
import sys, re, subprocess, json
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; P = D.parent
sys.path[:0] = [str(P), str(P.parent / "ana003_ccx_fistr_mpc_spring"), str(P.parent / "tools")]
import foundation_model as fm
from vmc_model import build, write_input
from frequency_response import ALPHA, BETA
K0 = dict(fm.K)
CASES = {"k025": 0.25, "k1": 1.0, "k4": 4.0}
PTS = {"BEDPT": (0.0, 450.0, 100.0), "COLTOP": (500.0, 800.0, 1000.0)}
info = {}
for name, s in CASES.items():
    fm.K.clear(); fm.K.update({1: K0[1] * s, 2: K0[2] * s, 3: K0[3]})
    m = build("foundation"); out = D / name; out.mkdir(exist_ok=True); inp = out / "model.inp"
    write_input(m, inp)
    coord = {n: np.array((x, y, z)) for p in m.parts for n, x, y, z in p.coords}
    ids = np.array(list(coord)); X = np.array(list(coord.values()))
    pick = {k: int(ids[np.argmin(np.linalg.norm(X - np.array(v), axis=1))]) for k, v in PTS.items()}
    nut = next(sp for sp in m.springs if sp[0] == "NUT" and sp[1] == 1)
    t = inp.read_text()
    t = re.sub(r"(\*FREQUENCY[^\n]*\n)\s*80", r"\g<1>20, 1.0", t)
    t = t.replace("*BOUNDARY", "".join(f"*NSET, NSET={k}\n{v}\n" for k, v in pick.items()) + "*BOUNDARY", 1)
    t += "\n".join(["*STEP", "*STEADY STATE DYNAMICS", "1, 300, 40, 1.0", "*MODAL DAMPING, RAYLEIGH", f",,{ALPHA},{BETA}",
                    "*CLOAD, OP=NEW", f"{nut[4]}, 1, 1.0", f"{nut[3]}, 1, -1.0"] +
                   sum([[f"*NODE PRINT, NSET={k}", "U"] for k in ("BEDPT", "COLTOP", "WORK", "TOOL")], []) + ["*END STEP"]) + "\n"
    inp.write_text(t)
    r = subprocess.run([str(Path.home() / "opt/ccx.sh"), "-i", "model"], cwd=out, stdout=open(out / "ccx.log", "w"), stderr=subprocess.STDOUT)
    info[name] = dict(scale=s, k_h_N_per_um_per_foot=K0[1] * s / 1000, nodes=pick, nut=list(nut[3:5]), M=m.masses, rc=r.returncode)
    print(name, r.returncode, flush=True)
fm.K.clear(); fm.K.update(K0)
(D / "info.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
