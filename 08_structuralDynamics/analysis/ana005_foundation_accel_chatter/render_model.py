"""ana005：計算に使った立型マシニングセンタのモデルを、部品ごとに色分け・名前付きで描く（vmc_model.png と ParaView 用 vmc_model.vtm）"""
from pathlib import Path
import sys, numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from vmc_model import build, CONTACT, TOOL, HEAD, SPINDLE, WORK, TOOL_KIND
D = Path(__file__).resolve().parent; R = D / "results"; R.mkdir(exist_ok=True); sys.path.insert(0, str(D.parent / "tools"))
import vtu_tools as vt
import machine_model as mm
m = build("foundation"); g, ids = vt.make_grid(type("M", (), {"parts": m.parts, "springs": []})())
g.point_data["node_id"] = ids; g.cell_data["component"] = m.element_parts
cc = g.cell_centers().points; part = np.array(g.cell_data["part"])
# 分解図：テーブル側（キャリッジ・ナットブラケットを含む）と、z ≥ 260 の部分（コラム上部・主軸頭・主軸・テストバー）を 120 mm 持ち上げる
LIFT = 120.0; gd = g.copy(); P = gd.points.copy(); tabpts = np.unique(np.concatenate([gd.get_cell(c).point_ids for c in np.where(part == 1)[0]]))
lift = np.zeros(len(P), bool); lift[tabpts] = True; lift |= P[:, 2] >= 260; P[lift, 2] += LIFT; gd.points = P
inside = lambda b: (cc[:, 0] > b[0]) & (cc[:, 0] < b[1]) & (cc[:, 1] > b[2]) & (cc[:, 1] < b[3]) & (cc[:, 2] > b[4]) & (cc[:, 2] < b[5])
# 部品の名前・色・ラベル位置（mm）・ラベルのずらし（px）
B, T = mm.BED, mm.TABLE
PARTS = [("ベッド", (0, B[0]), "#9aa3ad", (150, 100, 100), (-230, 20)),
         ("レール（2 本、ベッドと一体）", (0, B[1]), "#5d6d7e", (60, 150, 222), (-330, -90)),
         ("", (0, B[2]), "#5d6d7e", None, None),
         ("ナット受け台（ベッド側）", (0, B[3]), "#d35400", (500, 300, 215), (-420, 200)),
         ("コラム", (0, B[4]), "#aab7b8", (500, 800, 950), (-260, -60)),
         ("主軸頭", (0, HEAD), "#85929e", (420, 500, 800), (-260, -40)),
         ("主軸（鋼）", (0, SPINDLE), "#5dade2", (500, 400, 580), (190, -60)),
         ("テストバー（鋼 60×60×100）" if TOOL_KIND == "bar" else "工具（超硬 20×20×100）", (0, TOOL), "#c0392b", (500, 400, 450), (200, -10)),
         ("テーブル", (1, T[0]), "#2e86c1", (700, 300, 350), (190, 20)),
         ("キャリッジ（4 個、テーブル側）", (1, T[1]), "#f1c40f", (690, 150, 230), (170, 110)),
         ("", (1, T[2]), "#f1c40f", None, None), ("", (1, T[3]), "#f1c40f", None, None), ("", (1, T[4]), "#f1c40f", None, None),
         ("ナットブラケット（テーブル側）", (1, T[5]), "#e67e22", (500, 300, 232), (-470, 0)),
         ("ワーク（鋳鉄 100×100×50）", (1, WORK), "#a04000", (550, 420, 380), (180, 40))]
p = pv.Plotter(off_screen=True, window_size=(1400, 1050)); p.set_background("white")
used = np.zeros(g.n_cells, bool)
for name, (pi, box), col, _, _ in PARTS[::-1]:
    sel = np.where((part == pi) & inside(box) & ~used)[0]; used[sel] = True
    if len(sel): p.add_mesh(gd.extract_cells(sel).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"), color=col, show_edges=False)
p.add_mesh(pv.Plane(center=(500, 450, -70), direction=(0, 0, 1), i_size=1150, j_size=1050), color="#d9e4cd", opacity=0.6)
p.add_points(CONTACT + np.array([0, 0, LIFT]), color="#ff1111", point_size=14, render_points_as_spheres=True)
p.camera_position = [(2100, -1800, 1600), (500, 450, 450), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(0.9); p.add_axes(); p.render()
H = 1050
def screen(q):
    q = np.array(q, float); q[2] += LIFT if q[2] >= 230 else 0
    p.renderer.SetWorldPoint(*q, 1); p.renderer.WorldToDisplay(); x, y, z = p.renderer.GetDisplayPoint(); return x, H - y
pos = [(n, screen(q), o) for n, _, _, q, o in PARTS if q is not None]
im = Image.fromarray(p.screenshot(return_img=True)); p.close(); dr = ImageDraw.Draw(im)
font = lambda n: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), n)
for t, (x, y), (dx, dy) in pos:
    tx, ty = x + dx, y + dy; dr.line((x, y, tx + (0 if dx > 0 else font(18).getlength(t)), ty + 13), fill="#333333", width=2)
    dr.rectangle((tx - 4, ty - 2, tx + font(18).getlength(t) + 5, ty + 26), fill="white", outline="#999999"); dr.text((tx, ty), t, font=font(18), fill="black")
dr.rectangle((0, 0, 1400, 70), fill="white")
dr.text((15, 10), "計算に使った立型マシニングセンタのモデル（分解図：テーブル側と上部を 120 mm 持ち上げて表示。赤点 = 切削点・打撃点）", font=font(24), fill="black")
dr.text((15, 990), "キャリッジとレールのあいだ、ナットブラケットとナット受け台のあいだが、ガイド・ボールねじのばね（長さ 0）。マシニングセンタなので心押し台（テールストック）はない。", font=font(17), fill="black")
im.save(R / "vmc_model.png")
g.save(R / "vmc_model.vtu"); pv.MultiBlock({"VMC_SOLIDS": g}).save(R / "vmc_model.vtm"); print("ok")
