"""簡易工作機械モデル（ベッド＋コラム＋リニアガイド＋ボールねじ＋テーブル）を
20節点六面体で作り、CalculiX(.inp) と FrontISTR(.msh/.cnt) の入力を書き出す。

単位系: mm, N, t(トン), s  →  応力 MPa, 密度 t/mm^3, ばね N/mm

ばねの考え方（面と面をつなぐばね）:
  レール上面（ベッド側）とキャリッジ下面（テーブル側）は同じ位置・同じ分割で、別々の節点を持つ。
  向かい合う節点の組を、方向ごとの線形ばね（CalculiX: SPRING2 / FrontISTR: 511 + !SPRING_D）でつなぐ。
  面全体のばね定数 K を、各節点が受け持つ面積の割合で配分する（四角形の面の 4 隅に 1/4 ずつ）。
  中間節点（辺の中点）にはばねを付けない。
"""
from __future__ import annotations
import itertools
from dataclasses import dataclass, field
import numpy as np

# ---------------------------------------------------------------- 形状（mm）
# 各部品は直方体の集まり。(x0, x1, y0, y1, z0, z1)
BED = [
    (0, 1000, 0, 900, 0, 200),          # ベッド本体
    (0, 1000, 130, 170, 200, 230),      # レール1（ベッドと一体）
    (0, 1000, 430, 470, 200, 230),      # レール2
    (460, 540, 280, 320, 200, 230),     # ボールねじ支持の代わりの台（ナットの真下）
    (350, 650, 650, 850, 200, 1000),    # コラム（ベッドと一体）
]
TABLE = [
    (250, 750, 100, 500, 270, 350),     # テーブル本体
    (270, 350, 130, 170, 230, 270),     # キャリッジ（レール1）
    (650, 730, 130, 170, 230, 270),
    (270, 350, 430, 470, 230, 270),     # キャリッジ（レール2）
    (650, 730, 430, 470, 230, 270),
    (460, 540, 280, 320, 230, 270),     # ナットブラケット
]
Z_IF = 230.0  # ばねでつなぐ面の高さ
# ばねでつなぐ面: 名前, (x0,x1,y0,y1), {方向: 面全体のばね定数 N/mm}
CARRIAGE_K = {2: 5.6e5, 3: 8.0e5}       # 横 560 N/um、ラジアル 800 N/um（1キャリッジあたり）
SCREW_K = {1: 2.0e5}                    # ボールねじ軸方向 200 N/um
INTERFACES = [
    ("CAR1", (270, 350, 130, 170), CARRIAGE_K),
    ("CAR2", (650, 730, 130, 170), CARRIAGE_K),
    ("CAR3", (270, 350, 430, 470), CARRIAGE_K),
    ("CAR4", (650, 730, 430, 470), CARRIAGE_K),
    ("NUT", (460, 540, 280, 320), SCREW_K),
]
MAT = dict(E=1.2e5, nu=0.27, rho=7.3e-9)   # 鋳鉄（FC300 相当）
H = 50.0                                     # 要素の最大寸法
MONITOR = (500.0, 300.0, 350.0)              # テーブル上面中央（加振点・応答点）


H_IF = 20.0                                  # ばね面の分割寸法


def _lines(boxes, axis, h):
    pts = {b[2 * axis] for b in boxes} | {b[2 * axis + 1] for b in boxes}
    if axis < 2:   # ばね面は細かく分割（x, y 方向）
        for _, r, _ in INTERFACES:
            a, b = r[2 * axis], r[2 * axis + 1]
            pts |= set(np.linspace(a, b, int(round((b - a) / H_IF)) + 1))
    pts = sorted(pts)
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(1, int(np.ceil((b - a) / h - 1e-9)))
        out += list(np.linspace(a, b, n + 1)[1:])
    return np.array(out)


@dataclass
class Part:
    name: str
    nodes: dict = field(default_factory=dict)   # key(coord) -> id
    coords: list = field(default_factory=list)
    elems: list = field(default_factory=list)   # 20 node ids


def _key(p):
    return tuple(int(round(v * 1000)) for v in p)


def build_part(name, boxes, lines, start_id):
    part = Part(name)
    xs, ys, zs = lines
    nid = [start_id]

    def node(p):
        k = _key(p)
        if k not in part.nodes:
            part.nodes[k] = nid[0]; part.coords.append((nid[0], *p)); nid[0] += 1
        return part.nodes[k]

    inside = lambda c: any(b[0] <= c[0] <= b[1] and b[2] <= c[1] <= b[3] and b[4] <= c[2] <= b[5] for b in boxes)
    for i, j, k in itertools.product(range(len(xs) - 1), range(len(ys) - 1), range(len(zs) - 1)):
        x0, x1, y0, y1, z0, z1 = xs[i], xs[i + 1], ys[j], ys[j + 1], zs[k], zs[k + 1]
        if not inside(((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2)):
            continue
        c = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
             (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
        ids = [node(p) for p in c] + [node(tuple((np.array(c[a]) + np.array(c[b])) / 2)) for a, b in edges]
        part.elems.append(ids)
    return part, nid[0]


@dataclass
class Model:
    parts: list
    springs: list          # (name, dof, k, node_bed, node_table)
    fixed: list            # ベッド底面の節点
    monitor: int
    nnode: int


def build(h=H):
    allb = BED + TABLE
    lines = [_lines(allb, a, h) for a in range(3)]
    bed, nxt = build_part("BED", BED, lines, 1)
    table, nxt = build_part("TABLE", TABLE, lines, nxt)
    # ばね: 界面のメッシュの角節点に面積で配分
    springs = []
    xs, ys = lines[0], lines[1]
    for name, (x0, x1, y0, y1), kdict in INTERFACES:
        w = {}
        area = (x1 - x0) * (y1 - y0)
        for i in range(len(xs) - 1):
            for j in range(len(ys) - 1):
                if not (x0 <= (xs[i] + xs[i + 1]) / 2 <= x1 and y0 <= (ys[j] + ys[j + 1]) / 2 <= y1):
                    continue
                a = (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j])
                for p in [(xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]), (xs[i], ys[j + 1])]:
                    w[p] = w.get(p, 0.0) + a / 4 / area
        assert abs(sum(w.values()) - 1) < 1e-9
        for p, wi in w.items():
            nb = bed.nodes[_key((p[0], p[1], Z_IF))]
            nt = table.nodes[_key((p[0], p[1], Z_IF))]
            for dof, K in kdict.items():
                springs.append((name, dof, K * wi, nb, nt))
    fixed = [nid for (nid, x, y, z) in bed.coords if abs(z) < 1e-6]
    mon = table.nodes[_key(MONITOR)]
    return Model([bed, table], springs, fixed, mon, nxt - 1)


def _chunks(seq, n=8):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def _spring_groups(model):
    """同じ (面, 方向, ばね定数) ごとにまとめる"""
    g = {}
    for name, dof, k, nb, nt in model.springs:
        g.setdefault((name, dof, round(k, 6)), []).append((nb, nt))
    return g


# ---------------------------------------------------------------- CalculiX
def write_ccx(model, path, nmodes=20, freq=None, alpha=0.0, beta=0.0, load_dof=1):
    L = ["*HEADING", "simple machine tool: bed+column / linear guides & ball screw as springs / table"]
    L.append("*NODE, NSET=NALL")
    for p in model.parts:
        L += [f"{n}, {x:.6g}, {y:.6g}, {z:.6g}" for n, x, y, z in p.coords]
    L += [f"{n}, {x:.6g}, {y:.6g}, {z:.6g}" for n, x, y, z in getattr(model, "extra_nodes", [])]
    eid = 1
    for p in model.parts:
        L.append(f"*ELEMENT, TYPE=C3D20, ELSET=E{p.name}")
        for e in p.elems:
            L.append(f"{eid}, " + ", ".join(map(str, e[:15])) + ",")
            L.append(", ".join(map(str, e[15:])))
            eid += 1
    for gi, ((name, dof, k), pairs) in enumerate(_spring_groups(model).items()):
        es = f"SP_{name}_{dof}_{gi}"
        L.append(f"*ELEMENT, TYPE=SPRING2, ELSET={es}")
        for nb, nt in pairs:
            L.append(f"{eid}, {nb}, {nt}"); eid += 1
        L += [f"*SPRING, ELSET={es}", f"{dof}, {dof}", f"{k:.6e}"]
    for eq in getattr(model, "equations", []):      # 最初の項が従属自由度
        L += ["*EQUATION", f"{len(eq)}"]
        for c in _chunks(eq, 4):
            L.append(", ".join(f"{n}, {d}, {a:.10g}" for n, d, a in c))
    L.append("*NSET, NSET=FIX")
    L += [", ".join(map(str, c)) + "," for c in _chunks(model.fixed)]
    L += ["*NSET, NSET=MON", f"{model.monitor},"]
    L += ["*MATERIAL, NAME=CASTIRON", "*ELASTIC", f"{MAT['E']}, {MAT['nu']}", "*DENSITY", f"{MAT['rho']}"]
    L += ["*SOLID SECTION, ELSET=EBED, MATERIAL=CASTIRON", "*SOLID SECTION, ELSET=ETABLE, MATERIAL=CASTIRON"]
    L += ["*BOUNDARY", "FIX, 1, 3"]
    L += ["*STEP", "*FREQUENCY, STORAGE=YES", f"{nmodes}", "*NODE FILE", "U", "*END STEP"]
    if freq:
        f0, f1, n = freq
        L += ["*STEP", "*STEADY STATE DYNAMICS", f"{f0}, {f1}, {n}, 1.0",
              "*MODAL DAMPING, RAYLEIGH", f",,{alpha:.6e},{beta:.6e}",
              "*CLOAD", f"MON, {load_dof}, 1.0",
              "*NODE PRINT, NSET=MON", "U", "*END STEP"]
    open(path, "w").write("\n".join(L) + "\n")


# ---------------------------------------------------------------- FrontISTR
def write_fistr_mesh(model, path):
    L = ["!HEADER", " simple machine tool", "!NODE"]
    for p in model.parts:
        L += [f"{n}, {x:.6g}, {y:.6g}, {z:.6g}" for n, x, y, z in p.coords]
    L += [f"{n}, {x:.6g}, {y:.6g}, {z:.6g}" for n, x, y, z in getattr(model, "extra_nodes", [])]
    eid = 1
    for p in model.parts:
        L.append(f"!ELEMENT, TYPE=362, EGRP=E{p.name}")
        for e in p.elems:
            L.append(f"{eid}, " + ", ".join(map(str, e))); eid += 1
    groups = _spring_groups(model)
    for gi, ((name, dof, k), pairs) in enumerate(groups.items()):
        L.append(f"!ELEMENT, TYPE=511, EGRP=SP{gi}")
        for nb, nt in pairs:
            L.append(f"{eid}, {nb}, {nt}"); eid += 1
    L += ["!MATERIAL, NAME=CASTIRON, ITEM=2", "!ITEM=1, SUBITEM=2", f"{MAT['E']}, {MAT['nu']}",
          "!ITEM=2, SUBITEM=1", f"{MAT['rho']}"]
    L += ["!SECTION, TYPE=SOLID, EGRP=EBED, MATERIAL=CASTIRON", "1.0",
          "!SECTION, TYPE=SOLID, EGRP=ETABLE, MATERIAL=CASTIRON", "1.0"]
    for gi in range(len(groups)):
        # cnt 側の !MATERIAL は msh 側に同名の材料が無いとエラーになるので、ダミーを置く
        L += [f"!MATERIAL, NAME=SPMAT{gi}, ITEM=1", "!ITEM=1, SUBITEM=2", "1.0, 0.3"]
        L += [f"!SECTION, TYPE=INTERFACE, EGRP=SP{gi}, MATERIAL=SPMAT{gi}", "1.0"]
    eqs = getattr(model, "equations", [])
    if eqs:
        L.append("!EQUATION")
        for eq in eqs:
            L.append(f"{len(eq)}, 0.0")
            for c in _chunks(eq, 7):
                L.append(", ".join(f"{n}, {d}, {a:.10g}" for n, d, a in c))
    L.append("!NGROUP, NGRP=FIX")
    L += [", ".join(map(str, c)) for c in _chunks(model.fixed)]
    L += ["!NGROUP, NGRP=MON", f"{model.monitor}", "!END"]
    open(path, "w").write("\n".join(L) + "\n")
    # ばねの材料は cnt 側に書く（!SPRING_D は制御ファイルのキーワード）
    mats = []
    for gi, ((name, dof, k), pairs) in enumerate(groups.items()):
        mats += [f"!MATERIAL, NAME=SPMAT{gi}", "!SPRING_D", f"{dof}, {dof}, {k:.6e}"]
    return mats


def write_fistr_eigen_cnt(path, spring_mats, nmodes=20):
    L = ["!VERSION", " 3", "!SOLUTION, TYPE=EIGEN", "!EIGEN", f" {nmodes}, 1.0E-8, 300",
         "!BOUNDARY", " FIX, 1, 3, 0.0"] + spring_mats + [
         "!SOLVER, METHOD=DIRECT", "!WRITE,RESULT", "!WRITE,VISUAL",
         "!VISUAL, method=PSR", "!surface_num=1", "!surface 1", "!output_type=VTK", "!END"]
    open(path, "w").write("\n".join(L) + "\n")


def write_fistr_freq_cnt(path, spring_mats, f0, f1, n, alpha, beta, monitor, nmodes, load_dof=1):
    L = ["!VERSION", " 3", "!SOLUTION, TYPE=DYNAMIC", "!DYNAMIC", " 11, 2",
         f" {f0}, {f1}, {n}, {f0}", " 0.0, 1.0", f" 1, 1, {alpha:.6e}, {beta:.6e}",
         f" 1, 2, {monitor}", " 1, 1, 1, 1, 1, 1",
         "!EIGENREAD", " eigen_0.log", f" 1, {nmodes}",
         "!BOUNDARY", " FIX, 1, 3, 0.0",
         "!FLOAD, LOAD CASE=1", f" MON, {load_dof}, 1.0"] + spring_mats + [
         "!SOLVER, METHOD=DIRECT", "!END"]
    open(path, "w").write("\n".join(L) + "\n")


def write_fistr_ctrl(path, mesh, cnt, res, eig_res=None, dyna=None):
    L = ["!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE", mesh, "!CONTROL, NAME=fstrCNT", cnt,
         "!RESULT, NAME=fstrRES, IO=OUT", res, "!RESULT, NAME=vis_out, IO=OUT", res.replace(".res", "_vis")]
    if dyna:
        L += ["!RESULT, NAME=fstrDYNA, IO=OUT", dyna]
    if eig_res:
        L += ["!RESULT, NAME=result-in, IO=IN", eig_res]
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    m = build()
    ne = sum(len(p.elems) for p in m.parts)
    print(f"nodes={m.nnode}, hex20={ne}, springs={len(m.springs)}, fixed={len(m.fixed)}, monitor={m.monitor}")
    for name, _, kd in INTERFACES:
        for dof in kd:
            print(name, dof, sum(s[2] for s in m.springs if s[0] == name and s[1] == dof))
