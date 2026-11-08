"""ana004: ばねの種類ごとの最小モデル（静解析）を CalculiX と FrontISTR で解き、手計算と比べる。
  A: 地面ばね      CalculiX SPRING1 + *SPRING      / FrontISTR !SPRING（cnt）
  B: 2節点・方向指定 CalculiX SPRING2 + *SPRING      / FrontISTR 511 + !SPRING_D   （長さ 0、向かい合う節点）
  C: 2節点・軸方向   CalculiX SPRINGA + *SPRING      / FrontISTR 511 + !SPRING_A   （長さ 50 mm、節点を結ぶ線の向き）
  実行: ~/opt/venv/bin/python make_and_run.py"""
import os, re, subprocess, glob
import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
E, NU, K = 210000.0, 0.3, 1000.0          # 鋼 [MPa]、ばね定数 1000 N/mm（= 1 N/µm）
FTOT = 1000.0                               # 全荷重 [N]
CASES = {}


def cube(n0, z0, a=100.0):
    """1 辺 a の立方体（8 節点）。節点番号 n0..n0+7。下面 4 → 上面 4（CalculiX・FrontISTR・VTK 共通の順番）"""
    xy = [(0, 0), (a, 0), (a, a), (0, a)]
    nodes = [(n0 + i, x, y, z0) for i, (x, y) in enumerate(xy)] + [(n0 + 4 + i, x, y, z0 + a) for i, (x, y) in enumerate(xy)]
    return nodes, [n0 + i for i in range(8)]


# ---------------------------------------------------------------- モデル定義
# A: 立方体 1 個を、下面 4 節点で地面ばね（X・Y・Z）に載せる。上面に -Z の荷重
nA, eA = cube(1, 0.0)
CASES["A_ground"] = dict(nodes=nA, elems=[eA], kind="ground", spr_nodes=[1, 2, 3, 4], fixed=[], load_nodes=[5, 6, 7, 8],
                         load_dof=3, mon=5, expect=-FTOT / (4 * K))
# B: 下の立方体（固定）と上の立方体を、同じ位置（z=100）の節点どうし 4 組の SPRING2 でつなぐ。上の立方体の下面に +X の荷重
nB1, eB1 = cube(1, 0.0); nB2, eB2 = cube(11, 100.0)
CASES["B_spring2"] = dict(nodes=nB1 + nB2, elems=[eB1, eB2], kind="dir", pairs=[(5, 11), (6, 12), (7, 13), (8, 14)],
                          fixed=[1, 2, 3, 4, 5, 6, 7, 8], load_nodes=[11, 12, 13, 14], load_dof=1, mon=11, expect=FTOT / (4 * K))
# C: 下の立方体（固定）と、50 mm 離れた上の立方体を、縦の SPRINGA 4 本でつなぐ。上の立方体は X・Y を拘束。下面に -Z の荷重
nC1, eC1 = cube(1, 0.0); nC2, eC2 = cube(11, 150.0)
CASES["C_springa"] = dict(nodes=nC1 + nC2, elems=[eC1, eC2], kind="axial", pairs=[(5, 11), (6, 12), (7, 13), (8, 14)],
                          fixed=[1, 2, 3, 4, 5, 6, 7, 8], fix12=list(range(11, 19)), load_nodes=[11, 12, 13, 14], load_dof=3,
                          mon=11, expect=-FTOT / (4 * K))


def write_ccx(c, path):
    L = ["*HEADING", "ana004 spring types", "*NODE, NSET=NALL"] + [f"{n}, {x}, {y}, {z}" for n, x, y, z in c["nodes"]]
    L += ["*ELEMENT, TYPE=C3D8, ELSET=EALL"] + [f"{100 + i}, " + ", ".join(map(str, e)) for i, e in enumerate(c["elems"])]
    L += ["*MATERIAL, NAME=STEEL", "*ELASTIC", f"{E}, {NU}", "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL"]
    eid = 200
    if c["kind"] == "ground":
        # SPRING1: 1 節点の要素。*SPRING の 2 行目は「自由度」1 つ。値ごと・自由度ごとに ELSET を分ける
        for dof in (1, 2, 3):
            L.append(f"*ELEMENT, TYPE=SPRING1, ELSET=GS{dof}")
            for n in c["spr_nodes"]:
                L.append(f"{eid}, {n}"); eid += 1
            L += [f"*SPRING, ELSET=GS{dof}", f"{dof}", f"{K}"]
    elif c["kind"] == "dir":
        # SPRING2: 2 節点の要素。*SPRING の 2 行目は「節点1の自由度, 節点2の自由度」
        for dof in (1, 2, 3):
            L.append(f"*ELEMENT, TYPE=SPRING2, ELSET=S2_{dof}")
            for a, b in c["pairs"]:
                L.append(f"{eid}, {a}, {b}"); eid += 1
            L += [f"*SPRING, ELSET=S2_{dof}", f"{dof}, {dof}", f"{K}"]
    else:
        # SPRINGA: 2 節点を結ぶ線の向きに効く。*SPRING の 2 行目は空行
        L.append("*ELEMENT, TYPE=SPRINGA, ELSET=SA")
        for a, b in c["pairs"]:
            L.append(f"{eid}, {a}, {b}"); eid += 1
        L += ["*SPRING, ELSET=SA", "", f"{K}"]
    L += ["*NSET, NSET=LOADN"] + [", ".join(map(str, c["load_nodes"]))]
    L += ["*STEP", "*STATIC", "*BOUNDARY"]
    L += [f"{n}, 1, 3, 0.0" for n in c["fixed"]] + [f"{n}, 1, 2, 0.0" for n in c.get("fix12", [])]
    L += ["*CLOAD"] + [f"{n}, {c['load_dof']}, {np.sign(c['expect']) * FTOT / 4}" for n in c["load_nodes"]]
    L += ["*NODE PRINT, NSET=NALL", "U", "*NODE FILE", "U", "*EL FILE", "S", "*END STEP"]
    open(path, "w").write("\n".join(L) + "\n")


def write_fistr(c, d):
    L = ["!HEADER", " ana004", "!NODE"] + [f"{n}, {x}, {y}, {z}" for n, x, y, z in c["nodes"]]
    L += ["!ELEMENT, TYPE=361, EGRP=ESOLID"] + [f"{100 + i}, " + ", ".join(map(str, e)) for i, e in enumerate(c["elems"])]
    mats_cnt = []
    if c["kind"] in ("dir", "axial"):
        L.append("!ELEMENT, TYPE=511, EGRP=ESPR")
        for i, (a, b) in enumerate(c["pairs"]):
            L.append(f"{200 + i}, {a}, {b}")
    L += ["!MATERIAL, NAME=STEEL, ITEM=1", "!ITEM=1, SUBITEM=2", f"{E}, {NU}", "!SECTION, TYPE=SOLID, EGRP=ESOLID, MATERIAL=STEEL", "1.0"]
    if c["kind"] in ("dir", "axial"):
        # cnt の !MATERIAL と同じ名前のダミー材料がメッシュ側に必要
        L += ["!MATERIAL, NAME=SPMAT, ITEM=1", "!ITEM=1, SUBITEM=2", "1.0, 0.3", "!SECTION, TYPE=INTERFACE, EGRP=ESPR, MATERIAL=SPMAT", "1.0"]
        if c["kind"] == "dir":
            mats_cnt = ["!MATERIAL, NAME=SPMAT", "!SPRING_D"] + [f" {dof}, {dof}, {K}" for dof in (1, 2, 3)]
        else:
            mats_cnt = ["!MATERIAL, NAME=SPMAT", "!SPRING_A", f" {K}"]
    L += ["!NGROUP, NGRP=LOADN"] + [", ".join(map(str, c["load_nodes"]))]
    if c["fixed"]:
        L += ["!NGROUP, NGRP=FIX"] + [", ".join(map(str, c["fixed"]))]
    if c.get("fix12"):
        L += ["!NGROUP, NGRP=FIX12"] + [", ".join(map(str, c["fix12"]))]
    if c["kind"] == "ground":
        L += ["!NGROUP, NGRP=SPRN"] + [", ".join(map(str, c["spr_nodes"]))]
    L.append("!END")
    open(os.path.join(d, "model.msh"), "w").write("\n".join(L) + "\n")
    C = ["!VERSION", " 3", "!SOLUTION, TYPE=STATIC"] + mats_cnt
    C += ["!BOUNDARY"] + ([" FIX, 1, 3, 0.0"] if c["fixed"] else []) + ([" FIX12, 1, 2, 0.0"] if c.get("fix12") else [])
    if c["kind"] == "ground":
        # !SPRING: 節点（グループ）と地面の間のばね。「節点, 自由度, ばね定数」
        C += ["!SPRING"] + [f" SPRN, {dof}, {K}" for dof in (1, 2, 3)]
    C += ["!CLOAD"] + [f" {n}, {c['load_dof']}, {np.sign(c['expect']) * FTOT / 4}" for n in c["load_nodes"]]
    C += ["!SOLVER, METHOD=DIRECT", "!WRITE,RESULT", "!WRITE,VISUAL", "!VISUAL, method=PSR", "!surface_num=1", "!surface 1",
          "!output_type=VTK", "!END"]
    if C[-1] == "!BOUNDARY":
        pass
    open(os.path.join(d, "model.cnt"), "w").write("\n".join(C) + "\n")
    open(os.path.join(d, "hecmw_ctrl.dat"), "w").write(
        "!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE\n model.msh\n!CONTROL, NAME=fstrCNT\n model.cnt\n"
        "!RESULT, NAME=fstrRES, IO=OUT\n model.res\n!RESULT, NAME=vis_out, IO=OUT\n model_vis\n")


def read_ccx_dat(path):
    out, on = {}, False
    for ln in open(path):
        if "displacements" in ln:
            on = True; continue
        p = ln.split()
        if on and len(p) == 4:
            out[int(p[0])] = np.array(list(map(float, p[1:])))
    return out


env = dict(os.environ)
R = os.path.expanduser("~/opt/ccxroot/usr/lib/x86_64-linux-gnu")
env["LD_LIBRARY_PATH"] = f"{R}/blas:{R}/lapack"
rows = ["| ケース | ばね | 理論値 [mm] | CalculiX [mm] | FrontISTR [mm] |", "|---|---|---|---|---|"]
for name, c in CASES.items():
    dc, df = os.path.join(D, name, "ccx"), os.path.join(D, name, "fistr")
    os.makedirs(dc, exist_ok=True); os.makedirs(df, exist_ok=True)
    write_ccx(c, os.path.join(dc, "model.inp")); write_fistr(c, df)
    subprocess.run([os.path.expanduser("~/opt/ccx.sh"), "-i", "model"], cwd=dc, stdout=open(os.path.join(dc, "ccx.log"), "w"), stderr=subprocess.STDOUT)
    subprocess.run([os.path.expanduser("~/opt/frontistr/bin/fistr1")], cwd=df, env=env, stdout=open(os.path.join(df, "fistr.out"), "w"), stderr=subprocess.STDOUT)
    uc = read_ccx_dat(os.path.join(dc, "model.dat"))[c["mon"]][c["load_dof"] - 1]
    import pyvista as pv
    vt = sorted(glob.glob(os.path.join(df, "model_vis_psf.*.pvtu")))
    m = pv.read(vt[-1]); ids = m.point_data.get("NODE_ID", None)
    if ids is None:
        ids = np.array([n for n, *_ in c["nodes"]])
    uf = m.point_data["DISPLACEMENT"][list(ids).index(c["mon"])][c["load_dof"] - 1]
    lab = {"ground": "地面ばね（SPRING1 / !SPRING）", "dir": "2節点・方向指定（SPRING2 / 511+!SPRING_D）", "axial": "2節点・軸方向（SPRINGA / 511+!SPRING_A）"}[c["kind"]]
    rows.append(f"| {name} | {lab} | {c['expect']:.4f} | {uc:.4f} | {uf:.4f} |")
    print(name, c["expect"], uc, uf)
os.makedirs(os.path.join(D, "results"), exist_ok=True)
open(os.path.join(D, "results", "check.md"), "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
