"""ana005 2.3：比べる 3 ケース（底面固定・基礎ばね・テーブル結合 10 倍）の境界条件とばねを並べて描く（分解図）
   元：render_springs.py
   ・ベッド底面の 4 か所の面 → MPC（紫）→ 代表点（黄）→ 基礎ばね X・Y・Z ＋回転（赤）→ 地面（灰の四角）
   ・ガイド 4 個とボールねじナットの面 → MPC → 代表点どうしのばね（青）
   ・底面固定モデルとの違い（右）
   出力：results/foundation_springs.png、ParaView 用 results/foundation_springs.vtm"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent / "tools"), str(D.parent / "ana003_ccx_fistr_mpc_spring")]
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from vmc_model import build
import vtu_tools as vt
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
m = build("foundation"); DROP, LIFT = 160.0, 160.0
coord = {n: np.array((x, y, z)) for p in m.parts for n, x, y, z in p.coords}
coord.update({n: np.array((x, y, z)) for n, x, y, z in m.extra_nodes})
tab = {c[0] for c in m.parts[1].coords}
g, ids = vt.make_grid(type("M", (), {"parts": m.parts, "springs": []})())
idx = {n: i for i, n in enumerate(ids)}; pts = g.points.copy(); pts[[idx[n] for n in tab], 2] += LIFT; g.points = pts
bed = g.extract_cells(np.where(g.cell_data["part"] == 0)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
tbl = g.extract_cells(np.where(g.cell_data["part"] == 1)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
foot = {f["ref"]: f for f in m.foundations}; footrot = {f["rot"] for f in m.foundations}; ground = {f["ground"] for f in m.foundations} | {f["ground_rot"] for f in m.foundations}


tab_aux = set()
for eq in m.equations:
    if eq[0][0] in tab: tab_aux |= {n for n, _, _ in eq[1:]}


def pos(n):
    p = coord[n].copy()
    if n in tab or n in tab_aux: p[2] += LIFT
    if n in ground: p[2] -= DROP
    return p


# 拘束式：従属節点（先頭）→ 代表点
refs_ifc = {s_[3] for s_ in m.springs if not s_[0].startswith("FOOT")} | {s_[4] for s_ in m.springs if not s_[0].startswith("FOOT")}
mpc_f, mpc_i = [], []
for eq in m.equations:
    dep = eq[0][0]; refn = [n for n, _, _ in eq[1:] if n in foot or n in refs_ifc]
    if not refn: continue
    seg = [pos(dep), pos(refn[0])]
    (mpc_f if refn[0] in foot else mpc_i).extend(seg)
spr_f, spr_i = [], []
for name, d, k, a, b in m.springs:
    if a in footrot or b in footrot: continue
    (spr_f if name.startswith("FOOT") else spr_i).extend([pos(a), pos(b)])
spr_i = [p + (0, 0, 0) for p in spr_i]

from PIL import ImageChops
bot = [coord[n] for n, x, y, z in m.parts[0].coords if abs(z) < 1e-8]


def render(case):
    p = pv.Plotter(off_screen=True, window_size=(1100, 1000)); p.set_background("white")
    p.add_mesh(bed, color="#b0b7bf", opacity=0.35, show_edges=True, edge_color="#777777", line_width=0.3)
    p.add_mesh(tbl, color="#7fb3d5", opacity=0.45, show_edges=True, edge_color="#2e4053", line_width=0.3)
    if mpc_i: p.add_mesh(pv.line_segments_from_points(np.array(mpc_i)), color="#5dade2", line_width=0.8)
    p.add_mesh(pv.line_segments_from_points(np.array(spr_i)), color="#27ae60" if case == "stiff" else "#1f618d", line_width=11 if case == "stiff" else 5)
    p.add_points(np.array([pos(n) for n in refs_ifc]), color="#f39c12", point_size=12, render_points_as_spheres=True)
    if case == "fixed":
        b = np.array(bot); p.add_points(b, color="#2c3e50", point_size=7, render_points_as_spheres=True)
        p.add_mesh(pv.Plane(center=(500, 450, -1), direction=(0, 0, 1), i_size=1000, j_size=900), color="#2c3e50", opacity=0.55)
        p.add_mesh(pv.Plane(center=(500, 450, -DROP - 21), direction=(0, 0, 1), i_size=1300, j_size=1200), color="#d5d8dc", opacity=0.0)
    else:
        if mpc_f: p.add_mesh(pv.line_segments_from_points(np.array(mpc_f)), color="#8e44ad", line_width=1.2)
        p.add_mesh(pv.line_segments_from_points(np.array(spr_f)), color="#c0392b", line_width=6)
        p.add_points(np.array([pos(n) for n in foot]), color="#f1c40f", point_size=18, render_points_as_spheres=True)
        for f in m.foundations:
            gp = pos(f["ground"]); p.add_mesh(pv.Cube(center=gp - (0, 0, 10), x_length=110, y_length=110, z_length=20), color="#566573")
        p.add_mesh(pv.Plane(center=(500, 450, -DROP - 21), direction=(0, 0, 1), i_size=1300, j_size=1200), color="#d5d8dc", opacity=0.5)
    p.camera_position = [(2400, -1900, 900), (500, 450, 150), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.1)
    im = Image.fromarray(p.screenshot(return_img=True)); p.close(); return im


ims = [render(c) for c in ("fixed", "foundation", "stiff")]
bb = [ImageChops.difference(i, Image.new("RGB", i.size, (255, 255, 255))).getbbox() for i in ims]
box = (min(b[0] for b in bb) - 10, min(b[1] for b in bb) - 10, max(b[2] for b in bb) + 10, max(b[3] for b in bb) + 10)
ims = [i.crop(box) for i in ims]; w, h = ims[0].size
TITLE = [("① 底面固定", "ベッド底面の全節点を固定（黒い面）", "ガイド・ナット：元の剛性（青）"),
         ("② 基礎ばね（標準）", "底面 4 か所の面 → 拘束式（紫）→ 代表点（黄）", "→ 基礎ばね（赤）→ 地面。ガイド・ナット：元の剛性（青）"),
         ("③ テーブル結合 10 倍", "基礎は②と同じ", "ガイド・ナットの並進・回転剛性を 10 倍（太い緑）")]
out = Image.new("RGB", (3 * w, h + 190), "white"); d = ImageDraw.Draw(out)
d.text((14, 8), "比べる 3 ケースの違い（分解図：テーブルを 160 mm 持ち上げ、地面を 160 mm 下げて表示。実際はどちらも長さ 0）", fill=(20, 20, 20), font=font(26))
for k, (im, (t1, t2, t3)) in enumerate(zip(ims, TITLE)):
    out.paste(im, (k * w, 190)); d.text((k * w + 14, 60), t1, fill=(20, 20, 20), font=font(30))
    d.text((k * w + 14, 106), t2, fill=(60, 60, 60), font=font(19)); d.text((k * w + 14, 134), t3, fill=(60, 60, 60), font=font(19))
    if k: d.line((k * w, 60, k * w, h + 180), fill=(200, 200, 200), width=2)
out.save(R / "case_compare_bc.png"); print("ok", out.size)
