"""部品ごとに Hex20 メッシュを作り、接合面を「代表点＋MPC＋ばね」でつなぎ、CalculiX 入力を書く。

接合面 1 つにつき、両側の部品にそれぞれ
  REF 節点（並進）と ROT 節点（並進の自由度 1,2,3 を回転角 θx,θy,θz として使う補助節点）
を面の中心に置き、面上の全節点を  u_i = u_REF + θ × r_i  で REF/ROT に剛体結合する（*EQUATION）。
REF どうしを並進ばね、ROT どうしを回転ばね（面に一様に分布したばねと同じモーメント剛性）でつなぐ。
"""
import itertools
from dataclasses import dataclass, field
import numpy as np
import params as P
from lib.geometry import AXIS, build as build_geometry


def _key(p):
    return tuple(int(round(v * 1000)) for v in p)


def grid_lines(values, h):
    pts = sorted(set(round(v, 6) for v in values))
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        n = max(1, int(np.ceil((b - a) / h - 1e-9)))
        out += list(np.linspace(a, b, n + 1)[1:])
    return np.array(out)


@dataclass
class Model:
    coords: dict = field(default_factory=dict)        # id -> (x,y,z)
    part_nodes: dict = field(default_factory=dict)    # part -> {key: id}
    elems: list = field(default_factory=list)         # (part, material, [20 ids])
    springs: list = field(default_factory=list)       # (name, dof, k, n1, n2)
    equations: list = field(default_factory=list)
    fixed: list = field(default_factory=list)
    info: dict = field(default_factory=dict)
    nnode: int = 0

    def new_node(self, xyz):
        self.nnode += 1; self.coords[self.nnode] = tuple(map(float, xyz)); return self.nnode


def mesh_part(m, name, boxes, extra, h):
    lines = [grid_lines([b[2 * a] for b, _ in boxes] + [b[2 * a + 1] for b, _ in boxes] + extra[a], h) for a in range(3)]
    nodes = {}
    def node(p):
        k = _key(p)
        if k not in nodes:
            nodes[k] = m.new_node(p)
        return nodes[k]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    xs, ys, zs = lines
    for i, j, k in itertools.product(range(len(xs) - 1), range(len(ys) - 1), range(len(zs) - 1)):
        c = ((xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2, (zs[k] + zs[k + 1]) / 2)
        mat = next((mt for b, mt in boxes if b[0] <= c[0] <= b[1] and b[2] <= c[1] <= b[3] and b[4] <= c[2] <= b[5]), None)
        if mat is None:
            continue
        x0, x1, y0, y1, z0, z1 = xs[i], xs[i + 1], ys[j], ys[j + 1], zs[k], zs[k + 1]
        cr = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0), (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
        ids = [node(p) for p in cr] + [node(tuple((np.array(cr[a]) + np.array(cr[b])) / 2)) for a, b in edges]
        m.elems.append((name, mat, ids))
    m.part_nodes[name] = nodes


def face_nodes(m, part, normal, c, rect):
    a = AXIS[normal]; out = []
    for n in m.part_nodes[part].values():
        p = m.coords[n]
        if abs(p[a] - c) > 1e-6:
            continue
        if all(lo - 1e-6 <= p[AXIS[ax]] <= hi + 1e-6 for ax, (lo, hi) in rect.items()):
            out.append(n)
    return out


def rigid_patch(m, nodes, center):
    """面の節点を REF/ROT に剛体結合。u_i - u_ref - θ×r = 0"""
    ref, rot = m.new_node(center), m.new_node(center)
    for n in nodes:
        rx, ry, rz = np.array(m.coords[n]) - np.array(center)
        for eq in ([(n, 1, 1.), (ref, 1, -1.), (rot, 2, -rz), (rot, 3, ry)],
                   [(n, 2, 1.), (ref, 2, -1.), (rot, 3, -rx), (rot, 1, rz)],
                   [(n, 3, 1.), (ref, 3, -1.), (rot, 1, -ry), (rot, 2, rx)]):
            m.equations.append([t for t in eq if abs(t[2]) > 1e-12])
    return ref, rot


def rot_stiffness(normal, rect, k):
    """面に一様分布した並進ばね（合計 k）と同じモーメント剛性。返り値: 全体の回転軸(1,2,3) → N·mm/rad"""
    n = AXIS[normal]; (p, (p0, p1)), (q, (q0, q1)) = rect.items(); ip, iq = AXIS[p], AXIS[q]
    Lp, Lq = p1 - p0, q1 - q0
    kn, kp, kq = k.get(n + 1, 0.), k.get(ip + 1, 0.), k.get(iq + 1, 0.)
    kr = {ip + 1: kn * Lq ** 2 / 12, iq + 1: kn * Lp ** 2 / 12, n + 1: kp * Lq ** 2 / 12 + kq * Lp ** 2 / 12}
    return {d: v for d, v in kr.items() if v > 0}


def build(dy=0.0, dz=0.0, h=P.H):
    parts, joints, feet, contact = build_geometry(dy, dz)
    extra = {name: [[], [], []] for name in parts}
    for j in joints:                      # 接合面の輪郭を両側のメッシュの分割線に入れる
        for side in (j.part_a, j.part_b):
            extra[side][AXIS[j.normal]].append(j.c)
            for ax, (lo, hi) in j.rect.items():
                extra[side][AXIS[ax]] += [lo, hi]
    for x0, x1, y0, y1 in feet:
        extra["BED"][0] += [x0, x1]; extra["BED"][1] += [y0, y1]
    for side in ("HEAD", "TABLE"):
        extra[side][0] += list(contact["rect"]["x"]); extra[side][1] += list(contact["rect"]["y"]) + [contact["xyz"][1]]
        extra[side][0] += [contact["xyz"][0]]
    m = Model()
    for name, boxes in parts.items():
        mesh_part(m, name, boxes, extra[name], h)
    m.info["parts"] = {name: [list(b) + [mt] for b, mt in boxes] for name, boxes in parts.items()}
    m.info["joints"] = []
    for j in joints:
        center = [0., 0., 0.]; center[AXIS[j.normal]] = j.c
        for ax, (lo, hi) in j.rect.items():
            center[AXIS[ax]] = (lo + hi) / 2
        refs = []
        for side in (j.part_a, j.part_b):
            nodes = face_nodes(m, side, j.normal, j.c, j.rect)
            assert len(nodes) >= 4, (j.name, side, len(nodes))
            refs.append(rigid_patch(m, nodes, center) + (len(nodes),))
        kr = rot_stiffness(j.normal, j.rect, j.k)
        for d, k in j.k.items():
            m.springs.append((j.name, d, k, refs[0][0], refs[1][0]))
        for d, k in kr.items():
            m.springs.append((j.name + "_ROT", d, k, refs[0][1], refs[1][1]))
        m.info["joints"].append(dict(name=j.name, a=j.part_a, b=j.part_b, normal=j.normal, c=j.c, rect=j.rect, center=center,
                                     k=j.k, kr=kr, ref_a=refs[0][0], rot_a=refs[0][1], n_a=refs[0][2],
                                     ref_b=refs[1][0], rot_b=refs[1][1], n_b=refs[1][2]))
    m.info["feet"] = []
    for i, (x0, x1, y0, y1) in enumerate(feet, 1):        # 基礎：ベッド底面の 4 か所 → 地面（固定）へばね
        center = ((x0 + x1) / 2, (y0 + y1) / 2, 0.)
        nodes = face_nodes(m, "BED", "z", 0., {"x": (x0, x1), "y": (y0, y1)})
        ref, rot = rigid_patch(m, nodes, center); g, grot = m.new_node(center), m.new_node(center)
        kr = rot_stiffness("z", {"x": (x0, x1), "y": (y0, y1)}, P.K_FOOT)
        for d, k in P.K_FOOT.items():
            m.springs.append((f"FOOT{i}", d, k, ref, g))
        for d, k in kr.items():
            m.springs.append((f"FOOT{i}_ROT", d, k, rot, grot))
        m.fixed += [g, grot]
        m.info["feet"].append(dict(name=f"FOOT{i}", center=center, ref=ref, rot=rot, ground=g, n=len(nodes), kr=kr))
    xyz = contact["xyz"]
    m.info["contact"] = dict(xyz=list(xyz),
                             tool_nodes=face_nodes(m, "HEAD", "z", xyz[2], contact["rect"]),
                             work_nodes=face_nodes(m, "TABLE", "z", xyz[2], contact["rect"]))
    assert len(m.info["contact"]["tool_nodes"]) >= 4 and len(m.info["contact"]["work_nodes"]) >= 4
    mass = {}
    for part, mat, ids in m.elems:
        p = np.array([m.coords[n] for n in ids[:8]]); v = np.prod(np.ptp(p, axis=0))
        mass[part] = mass.get(part, 0.) + v * P.MATERIALS[mat]["rho"] * 1000
    m.info["mass_kg"] = mass; m.info["dy"], m.info["dz"] = dy, dz
    m.info["nnode"], m.info["nelem"] = m.nnode, len(m.elems)
    return m


def write_ccx(m, path, nmodes=P.NMODES):
    L = ["*HEADING", f"ana006 VMC parts joined by guide springs dy={m.info['dy']} dz={m.info['dz']}", "*NODE, NSET=NALL"]
    L += [f"{n}, {x:.6g}, {y:.6g}, {z:.6g}" for n, (x, y, z) in m.coords.items()]
    sets = {}
    L.append("*ELEMENT, TYPE=C3D20, ELSET=EALL")
    for eid, (part, mat, ids) in enumerate(m.elems, 1):
        L.append(f"{eid}, " + ", ".join(map(str, ids[:15])) + ","); L.append(", ".join(map(str, ids[15:])))
        sets.setdefault("M_" + mat, []).append(eid); sets.setdefault("P_" + part, []).append(eid)
    eid = len(m.elems)
    for i, (name, d, k, a, b) in enumerate(m.springs):
        eid += 1; es = f"SP{i}_{name}_{d}"
        L += [f"*ELEMENT, TYPE=SPRING2, ELSET={es}", f"{eid}, {a}, {b}", f"*SPRING, ELSET={es}", f"{d}, {d}", f"{k:.6e}"]
    for name, ids in sets.items():
        L.append(f"*ELSET, ELSET={name}")
        L += [", ".join(map(str, ids[i:i + 16])) + "," for i in range(0, len(ids), 16)]
    for eq in m.equations:
        L += ["*EQUATION", str(len(eq))]
        for i in range(0, len(eq), 4):
            L.append(", ".join(f"{n}, {d}, {a:.10g}" for n, d, a in eq[i:i + 4]))
    L += ["*NSET, NSET=FIX", ", ".join(map(str, m.fixed))]
    for mat, v in P.MATERIALS.items():
        L += [f"*MATERIAL, NAME={mat}", "*ELASTIC", f"{v['E']}, {v['nu']}", "*DENSITY", f"{v['rho']}"]
        if "M_" + mat in sets:
            L.append(f"*SOLID SECTION, ELSET=M_{mat}, MATERIAL={mat}")
    L += ["*BOUNDARY", "FIX, 1, 3",
          "*STEP", "*FREQUENCY, STORAGE=NO", f"{nmodes}", "*NODE FILE", "U", "*END STEP"]
    open(path, "w").write("\n".join(L) + "\n")
