"""面と面のばねを「代表点＋MPC（拘束式）」でつなぐモデルを作る。

kind="rbe2": 面を剛体とみなす。面の全節点 i について  u_i = u_ref + θ × r_i
             （従属＝面の節点、独立＝代表点 REF と回転用の補助節点 ROT）
kind="rbe3": 代表点の動き＝面の節点の重み付き平均（面は剛体にならない）
             u_ref = Σ w_i u_i,  θ = 面の節点の動きから最小二乗で求めた回転
             （従属＝代表点 REF と ROT）
ROT 節点: 並進の自由度 1,2,3 を「回転角 θx,θy,θz」として使う補助節点（要素なし）。
ばね: REF 同士に並進ばね、ROT 同士に回転ばね（分布ばねと同じ合計剛性・同じモーメント剛性）。
"""
import copy
import numpy as np
import machine_model as mm


def _patch_nodes(part, rect, z):
    x0, x1, y0, y1 = rect
    return [(n, x, y) for n, x, y, zz in part.coords
            if abs(zz - z) < 1e-6 and x0 - 1e-6 <= x <= x1 + 1e-6 and y0 - 1e-6 <= y <= y1 + 1e-6]


def build(kind="rbe2", h=mm.H):
    m = mm.build(h)
    bed, tab = m.parts
    coord = {n: (x, y, z) for p in m.parts for n, x, y, z in p.coords}
    nid = m.nnode + 1
    extra, eqs, springs = [], [], []
    for name, rect, kdict in mm.INTERFACES:
        # 分布ばねモデルの重み（角節点）から、合計の並進剛性と回転剛性を求める
        sp = [s for s in m.springs if s[0] == name]
        xc, yc = (rect[0] + rect[1]) / 2, (rect[2] + rect[3]) / 2
        Kt = {d: sum(s[2] for s in sp if s[1] == d) for d in (1, 2, 3)}
        Kr = {1: sum(s[2] * (coord[s[3]][1] - yc) ** 2 for s in sp if s[1] == 3),          # θx: Z ばね × y^2
              2: sum(s[2] * (coord[s[3]][0] - xc) ** 2 for s in sp if s[1] == 3),          # θy: Z ばね × x^2
              3: sum(s[2] * (coord[s[3]][0] - xc) ** 2 for s in sp if s[1] == 2)           # θz: Y ばね × x^2
                 + sum(s[2] * (coord[s[3]][1] - yc) ** 2 for s in sp if s[1] == 1)}        #     + X ばね × y^2
        w_corner = {}
        for s in sp:
            w_corner[s[3]] = s[2] / Kt[s[1]]   # ベッド側の節点の重み（方向によらず同じ）
        refs = {}
        for side, part in (("B", bed), ("T", tab)):
            ref, rot = nid, nid + 1; nid += 2
            extra += [(ref, xc, yc, mm.Z_IF), (rot, xc, yc, mm.Z_IF)]
            refs[side] = (ref, rot)
            nodes = _patch_nodes(part, rect, mm.Z_IF)
            if kind == "rbe2":
                for n, x, y in nodes:           # u_i - u_ref - θ×r = 0（r_z = 0）
                    rx, ry = x - xc, y - yc
                    eqs.append([(n, 1, 1.0), (ref, 1, -1.0), (rot, 3, ry)])            # ux = uref_x - θz*ry
                    eqs.append([(n, 2, 1.0), (ref, 2, -1.0), (rot, 3, -rx)])           # uy = uref_y + θz*rx
                    eqs.append([(n, 3, 1.0), (ref, 3, -1.0), (rot, 1, -ry), (rot, 2, rx)])  # uz = uref_z + θx*ry - θy*rx
            else:
                # 角節点のみ、重み w（ベッド側の節点番号 → テーブル側は同じ座標の節点）
                cn = {}
                for s in sp:
                    nb, nt = s[3], s[4]
                    cn[nb if side == "B" else nt] = w_corner[nb]
                items = [(n, coord[n][0] - xc, coord[n][1] - yc, w) for n, w in cn.items()]
                for d in (1, 2, 3):
                    eqs.append([(ref, d, -1.0)] + [(n, d, w) for n, rx, ry, w in items])
                Sy = sum(w * ry * ry for n, rx, ry, w in items); Sx = sum(w * rx * rx for n, rx, ry, w in items)
                eqs.append([(rot, 1, -1.0)] + [(n, 3, w * ry / Sy) for n, rx, ry, w in items if abs(ry) > 1e-9])
                eqs.append([(rot, 2, -1.0)] + [(n, 3, -w * rx / Sx) for n, rx, ry, w in items if abs(rx) > 1e-9])
                eqs.append([(rot, 3, -1.0)] + [(n, 2, w * rx / (Sx + Sy)) for n, rx, ry, w in items if abs(rx) > 1e-9]
                           + [(n, 1, -w * ry / (Sx + Sy)) for n, rx, ry, w in items if abs(ry) > 1e-9])
        for d, K in Kt.items():
            if K > 0:
                springs.append((name, d, K, refs["B"][0], refs["T"][0]))
        for d, K in Kr.items():
            if K > 0:
                springs.append((name + "_ROT", d, K, refs["B"][1], refs["T"][1]))
    m2 = copy.copy(m)
    m2.springs = springs; m2.extra_nodes = extra; m2.equations = eqs; m2.nnode = nid - 1
    return m2
