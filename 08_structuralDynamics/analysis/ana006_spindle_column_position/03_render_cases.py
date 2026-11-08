"""ケースのモデル図（部品の色分け・接合面の代表点）を ParaView と同じ VTK（pyvista）で描く → results/cases_model.png"""
import sys
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; sys.path.insert(0, str(D)); R = D / "results"
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import params as P
from lib.fem import build
from lib.view import grid, surface, PART_ID, PART_COLOR, MAT_COLOR
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
JCOL = {"X": "#1f77b4", "Y": "#d62728", "Z": "#9467bd"}
SHOW = [("base", "基準（dy = 0, dz = 0）"), ("z300", "主軸頭 +300 mm 上（dz = 300）"), ("ym100", "コラム 100 mm 手前（dy = −100）"), ("yp050", "コラム 50 mm 奥（dy = +50）")]


def shot(name, out, side=False):
    dy, dz = P.CASES[name]; m = build(dy, dz); g, ids = grid(m)
    pl = pv.Plotter(window_size=(900, 900)); pl.set_background("white")
    for part, pid in PART_ID.items():
        cond = (g.cell_data["part"] == pid) & (g.cell_data["mat"] == 0)
        pl.add_mesh(surface(g, cond), color=PART_COLOR[part], show_edges=True, edge_color="#666666", line_width=0.3, opacity=0.55 if part == "COLUMN" else 1)
    for mi, mat in ((1, "SPINDLE"), (2, "TOOL")):
        pl.add_mesh(surface(g, g.cell_data["mat"] == mi), color=MAT_COLOR[mat], show_edges=True, edge_color="#222222", line_width=0.3)
    for j in m.info["joints"]:
        pl.add_mesh(pv.Sphere(radius=14, center=j["center"]), color=JCOL[j["name"][0]])
    for f in m.info["feet"]:
        pl.add_mesh(pv.Sphere(radius=18, center=f["center"]), color="#ff7f0e")
    pl.add_mesh(pv.Sphere(radius=16, center=m.info["contact"]["xyz"]), color="#ffdd00")
    if side:   # 横（+X 側）から平行投影：Y と Z の位置の違いが見える
        pl.camera_position = [(3000, 560, 620), (500, 560, 620), (0, 0, 1)]; pl.enable_parallel_projection(); pl.camera.parallel_scale = 650
        for yy in (400,): pl.add_mesh(pv.Line((1001, yy, 0), (1001, yy, 1250)), color="#c9a800", line_width=2)
    else:
        pl.camera_position = [(-1700, -2300, 2100), (500, 560, 600), (0, 0, 1)]; pl.camera.zoom(1.05)
    pl.show_axes(); pl.screenshot(out); pl.close()


imgs = []
for side in (False, True):
    for name, lab in SHOW:
        p = R / f"model_{name}.png"; shot(name, str(p), side); im = Image.open(p).convert("RGB"); d = ImageDraw.Draw(im)
        d.text((20, 15), lab if not side else lab + "：横から", fill="black", font=font(34)); imgs.append(im)
W, H = imgs[0].size; canvas = Image.new("RGB", (W * 4, 2 * H + 110), "white")
for i, im in enumerate(imgs):
    canvas.paste(im, ((i % 4) * W, (i // 4) * H))
d = ImageDraw.Draw(canvas); x = 30; y = 2 * H + 30
for lab, col in [("ベッド", PART_COLOR["BED"]), ("テーブル（治具・ワーク）", PART_COLOR["TABLE"]), ("コラム", PART_COLOR["COLUMN"]), ("主軸頭", PART_COLOR["HEAD"]),
                 ("● X 案内（テーブル↔ベッド）", JCOL["X"]), ("● Y 案内（コラム↔ベッド）", JCOL["Y"]), ("● Z 案内（主軸頭↔コラム）", JCOL["Z"]), ("● 基礎の脚", "#ff7f0e"), ("● 刃先―ワーク（横図の黄線 = 基準の Y = 400）", "#c9a800")]:
    if lab.startswith("●"):
        d.text((x, y), lab, fill=col, font=font(30)); x += int(d.textlength(lab, font=font(30))) + 40
    else:
        d.rectangle([x, y + 4, x + 34, y + 34], fill=col, outline="black"); d.text((x + 44, y), lab, fill="black", font=font(30)); x += 44 + int(d.textlength(lab, font=font(30))) + 40
canvas = canvas.resize((canvas.width // 2, canvas.height // 2), Image.LANCZOS); canvas.save(R / "cases_model.png"); print(canvas.size)
for name, _ in SHOW:
    (R / f"model_{name}.png").unlink()
