"""ana001: 固有モードを GIF アニメーションにする（斜め・上から・正面の 3 方向を並べる）
   実行: ~/opt/venv/bin/python make_gifs.py   → results/gif/mode0N.gif, results/gif/modes_all.gif"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm
from frd import read_frd
import vtu_tools as vt
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from gif_tools import save_gif

D = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(D, "results", "gif"); os.makedirs(OUT, exist_ok=True)
FONT = os.path.expanduser("~/.fonts/meiryo.ttc")
NF, SCALE, VIEW = 20, 60.0, (560, 440)

DESC = ["コラムが前後（Y）に曲がる ＝ お辞儀する曲げ",
        "テーブルが送り方向（X）に並進 ＝ ボールねじが伸び縮み",
        "コラムが左右（X）に曲がる ＝ 横に倒れる曲げ",
        "テーブルが横（Y）に並進 ＝ ガイドの横ばねが伸び縮み",
        "テーブルが Z 軸まわりにねじれる（ヨーイング）",
        "テーブルが上下（Z）に動き、曲がる"]
KIND = ["曲げ", "並進", "曲げ", "並進", "ねじれ（回転）", "上下＋曲げ"]
CAMS = {"斜めから": [(2600, -1900, 1900), (500, 400, 400), (0, 0, 1)],
        "上から（XY 面）": [(500, 400, 4200), (500, 400, 300), (0, 1, 0)],
        "正面から（XZ 面）": [(500, -3600, 500), (500, 400, 500), (0, 0, 1)]}

model = mm.build()
g, ids = vt.make_grid(model)
_, blocks = read_frd(os.path.join(D, "ccx", "machine.frd"))
vt.add_modes(g, ids, blocks)
fields = [k for k in g.point_data.keys() if k.startswith("mode")]
font = lambda s: ImageFont.truetype(FONT, s)


def solid_with(field):
    u = g.point_data[field]; mx = np.linalg.norm(u, axis=1).max()
    gg = g.copy(); gg.point_data["U"] = u / mx; gg.point_data["|U|"] = np.linalg.norm(u, axis=1) / mx
    return gg.extract_cells(np.where(gg.cell_data["part"] < 2)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")


def frames_for(field, cams, size=VIEW, scale=SCALE):
    s = solid_with(field); base = s.points.copy(); U = s.point_data["U"]
    per_view = []
    for cpos in cams:
        p = pv.Plotter(off_screen=True, window_size=size)
        p.add_mesh(s.copy(), style="wireframe", color="gray", opacity=0.25, line_width=0.5)
        w = s.copy()
        p.add_mesh(w, scalars="|U|", cmap="turbo", clim=(0, 1), show_edges=True, edge_color="#333333", line_width=0.2,
                   show_scalar_bar=False)
        p.camera_position = cpos
        if cpos[0][2] > 4000 or cpos[0][1] < -3000:
            p.enable_parallel_projection(); p.camera.zoom(1.45)
        p.add_axes(line_width=3, labels_off=False)
        imgs = []
        for i in range(NF):
            w.points = base + U * scale * np.sin(2 * np.pi * i / NF)
            p.render(); imgs.append(Image.fromarray(p.screenshot(return_img=True)))
        p.close(); per_view.append(imgs)
    return per_view


def compose(per_view, labels, header, sub):
    W, H = per_view[0][0].size; top = 92
    out = []
    for i in range(NF):
        im = Image.new("RGB", (W * len(per_view), H + top), "white"); d = ImageDraw.Draw(im)
        d.text((14, 8), header, fill=(20, 20, 20), font=font(26))
        d.text((14, 46), sub, fill=(192, 57, 43), font=font(22))
        for j, v in enumerate(per_view):
            im.paste(v[i], (W * j, top)); d.text((W * j + 10, top + 6), labels[j], fill=(40, 40, 40), font=font(20))
        out.append(im)
    return out


labels = list(CAMS.keys())
for k in range(6):
    f = float(fields[k].split("_")[1][:-2])
    pvw = frames_for(fields[k], list(CAMS.values()))
    fr = compose(pvw, labels, f"モード {k+1}：{f:.0f} Hz 【{KIND[k]}】", DESC[k] + "（変形は大きく拡大、灰色の線は元の形）")
    save_gif(fr, os.path.join(OUT, f"mode{k+1:02d}.gif"))
    print("mode", k + 1)

# 6 モードを 1 枚に（斜めから）
tiles = [frames_for(fields[k], [CAMS["斜めから"]], size=(460, 380))[0] for k in range(6)]
out = []
for i in range(NF):
    im = Image.new("RGB", (460 * 3, (380 + 40) * 2), "white"); d = ImageDraw.Draw(im)
    for k in range(6):
        x, y = 460 * (k % 3), (380 + 40) * (k // 3)
        f = float(fields[k].split("_")[1][:-2])
        d.text((x + 10, y + 6), f"{k+1}次 {f:.0f} Hz【{KIND[k]}】", fill=(20, 20, 20), font=font(20))
        im.paste(tiles[k][i], (x, y + 40))
    out.append(im)
save_gif(out, os.path.join(OUT, "modes_all.gif"))
print("done")
