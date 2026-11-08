"""rep_spring.py の使い方の例：立方体 2 個を、代表点どうしの 6 成分ばね 1 組でつなぐ（手順書の例題と同じ条件）。
   CalculiX と FrontISTR の両方で解き、手書きの例題（../ccx, ../fistr）と同じ結果になることを確かめる"""
import sys, os, subprocess, re
from pathlib import Path
D = Path(__file__).resolve().parent; sys.path.insert(0, str(D.parents[2] / "tools"))
from rep_spring import RepSpring
P = {1:(0,0,0),2:(100,0,0),3:(100,100,0),4:(0,100,0),5:(0,0,100),6:(100,0,100),7:(100,100,100),8:(0,100,100),
     9:(0,0,100),10:(100,0,100),11:(100,100,100),12:(0,100,100),13:(0,0,200),14:(100,0,200),15:(100,100,200),16:(0,100,200)}
rs = RepSpring(first_node=101, first_elem=101)
c = rs.add("JOINT", [5, 6, 7, 8], [9, 10, 11, 12], center=(50, 50, 100), k=(1000, 1000, 1000, 1e6, 1e6, 1e6), coords=P)
print("代表点:", c)
# ---- CalculiX
(D / "ccx").mkdir(exist_ok=True)
L = ["*NODE"] + [f"{n}, {x}, {y}, {z}" for n, (x, y, z) in P.items()]
L += ["*ELEMENT, TYPE=C3D8, ELSET=SOLIDS", "1, 1,2,3,4,5,6,7,8", "2, 9,10,11,12,13,14,15,16"] + rs.ccx()
L += ["*NSET, NSET=FIX", "1,2,3,4", "*NSET, NSET=MON", f"{c['refA']},{c['refB']},13",
      "*MATERIAL, NAME=STEEL", "*ELASTIC", "200000, 0.3", "*SOLID SECTION, ELSET=SOLIDS, MATERIAL=STEEL",
      "*BOUNDARY", "FIX, 1, 3", "*STEP", "*STATIC", "*CLOAD", "13, 2, 1.0", "*NODE PRINT, NSET=MON", "U", "*END STEP"]
(D / "ccx" / "model.inp").write_text("\n".join(L) + "\n")
subprocess.run([str(Path.home() / "opt/ccx.sh"), "-i", "model"], cwd=D / "ccx", stdout=open(D / "ccx" / "ccx.log", "w"), stderr=subprocess.STDOUT)
# ---- FrontISTR
(D / "fistr").mkdir(exist_ok=True); M, C = rs.fistr()
msh = ["!HEADER", " rep_spring demo", "!NODE"] + [f"{n}, {x}, {y}, {z}" for n, (x, y, z) in P.items()]
msh += ["!ELEMENT, TYPE=361, EGRP=SOLIDS", "1, 1,2,3,4,5,6,7,8", "2, 9,10,11,12,13,14,15,16"] + M
msh += ["!MATERIAL, NAME=STEEL, ITEM=1", "!ITEM=1, SUBITEM=2", "200000, 0.3", "!SECTION, TYPE=SOLID, EGRP=SOLIDS, MATERIAL=STEEL", "1.0",
        "!NGROUP, NGRP=FIX", "1,2,3,4", "!END"]
(D / "fistr" / "model.msh").write_text("\n".join(msh) + "\n")
(D / "fistr" / "model.cnt").write_text("\n".join(["!VERSION", "3", "!SOLUTION, TYPE=STATIC", "!BOUNDARY", "FIX, 1, 3, 0.0", "!CLOAD", "13, 2, 1.0"] + C +
                                                ["!SOLVER, METHOD=DIRECT", "!WRITE,RESULT", "!END"]) + "\n")
(D / "fistr" / "hecmw_ctrl.dat").write_text("!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE\nmodel.msh\n!CONTROL, NAME=fstrCNT\nmodel.cnt\n!RESULT, NAME=fstrRES, IO=OUT\nmodel.res\n")
lib = Path.home() / "opt/ccxroot/usr/lib/x86_64-linux-gnu"
subprocess.run([str(Path.home() / "opt/frontistr/bin/fistr1")], cwd=D / "fistr", env=dict(os.environ, LD_LIBRARY_PATH=f"{lib}/blas:{lib}/lapack"),
               stdout=open(D / "fistr" / "fistr.log", "w"), stderr=subprocess.STDOUT)
# ---- 結果
cc = {int(p[0]): float(p[2]) for p in (l.split() for l in (D / "ccx" / "model.dat").read_text().splitlines()) if len(p) == 4 and p[0].isdigit()}
lines = (D / "fistr" / "model.res.0.1").read_text().splitlines(); i = lines.index("*data"); nn = int(lines[i + 1].split()[0]); nc = int(lines[i + 2].split()[0])
nd = sum(map(int, lines[i + 3].split())); i += 4 + nc; ff = {}
for _ in range(nn):
    n = int(lines[i]); i += 1; v = []
    while len(v) < nd: v += list(map(float, lines[i].split())); i += 1
    ff[n] = v[1]
print(f"CalculiX : ばねの伸び {cc[c['refB']]-cc[c['refA']]:.6f} mm、節点 13 の Y {cc[13]:.6f} mm")
print(f"FrontISTR: ばねの伸び {ff[c['refB']]-ff[c['refA']]:.6f} mm、節点 13 の Y {ff[13]:.6f} mm")
print("手書きの例題: ばねの伸び 0.001000、節点 13 の Y 0.013502")
