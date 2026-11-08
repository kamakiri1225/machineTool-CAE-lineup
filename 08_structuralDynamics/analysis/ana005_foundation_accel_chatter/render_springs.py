"""ana005 2.3：基礎ばねモデルのばねの入れ方を図にする（分解図）
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
p = pv.Plotter(off_screen=True, window_size=(1500, 1100)); p.set_background("white")
p.add_mesh(bed, color="#b0b7bf", opacity=0.35, show_edges=True, edge_color="#777777", line_width=0.3)
p.add_mesh(tbl, color="#7fb3d5", opacity=0.45, show_edges=True, edge_color="#2e4053", line_width=0.3)
if mpc_f: p.add_mesh(pv.line_segments_from_points(np.array(mpc_f)), color="#8e44ad", line_width=1.2)
if mpc_i: p.add_mesh(pv.line_segments_from_points(np.array(mpc_i)), color="#5dade2", line_width=0.8)
p.add_mesh(pv.line_segments_from_points(np.array(spr_f)), color="#c0392b", line_width=6)
p.add_mesh(pv.line_segments_from_points(np.array(spr_i)), color="#1f618d", line_width=5)
p.add_points(np.array([pos(n) for n in foot]), color="#f1c40f", point_size=18, render_points_as_spheres=True)
p.add_points(np.array([pos(n) for n in refs_ifc]), color="#f39c12", point_size=12, render_points_as_spheres=True)
for f in m.foundations:
    gp = pos(f["ground"]); p.add_mesh(pv.Cube(center=gp - (0, 0, 10), x_length=110, y_length=110, z_length=20), color="#566573")
p.add_mesh(pv.Plane(center=(500, 450, -DROP - 21), direction=(0, 0, 1), i_size=1300, j_size=1200), color="#d5d8dc", opacity=0.5)
p.camera_position = [(2400, -1900, 900), (500, 450, 150), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.15)
im = Image.fromarray(p.screenshot(return_img=True)); p.close()
from PIL import ImageChops
bb = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()   # 余白を切り落とす
im = im.crop((max(bb[0] - 20, 0), max(bb[1] - 20, 0), min(bb[2] + 20, im.width), min(bb[3] + 20, im.height)))
W, H = im.size; W = max(W, 1500); out = Image.new("RGB", (W, H + 120), "white"); out.paste(im, ((W - im.width) // 2, 120)); d = ImageDraw.Draw(out)
d.text((14, 8), "基礎ばねモデルのばねの入れ方（分解図：テーブルを 160 mm 持ち上げ、地面を 160 mm 下げて表示。実際はどちらも長さ 0）", fill=(20, 20, 20), font=font(24))
items = [((192, 57, 43), "基礎ばね：脚の代表点 ― 地面（X・Y 10 N/µm、Z 100 N/µm、＋回転）"), ((142, 68, 173), "拘束式：ベッド底面 4 か所の面 → 脚の代表点"),
         ((241, 196, 15), "脚の代表点（4 か所）"), ((31, 97, 141), "ガイド 4 個・ボールねじナットのばね（ana001〜003 と同じ）"), ((86, 101, 115), "地面（固定）")]
for k, (c, t) in enumerate(items):
    x, y = 14 + (k % 2) * 740, 48 + (k // 2) * 26
    d.ellipse((x, y + 4, x + 16, y + 20), fill=c); d.text((x + 24, y), t, fill=(30, 30, 30), font=font(18))
out.save(R / "foundation_springs.png")
pv.MultiBlock({"bed": bed, "table": tbl, "foot_MPC": pv.line_segments_from_points(np.array(mpc_f)),
               "foot_springs": pv.line_segments_from_points(np.array(spr_f)), "guide_springs": pv.line_segments_from_points(np.array(spr_i))}).save(R / "foundation_springs.vtm")
print("ok", len(m.foundations), len(spr_f) // 2, len(spr_i) // 2)
