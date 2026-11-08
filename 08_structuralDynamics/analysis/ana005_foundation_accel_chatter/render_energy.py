"""ana005 5.3 補足：びびりを決めるモードで、どこに変形のエネルギー（ひずみエネルギー密度）がたまっているかを色で示す
   foundation_energy/model.frd（基礎ばね・テストバー、20 モード、*EL FILE ENER）を使う"""
import sys, json
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent / "tools")]
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)


def read_blocks(path, want):
    out = {}; mode = None; cur = None; name = None
    with open(path) as f:
        for ln in f:
            if ln.startswith("  100CL"):
                mode = int(ln[58:63]) if ln[58:63].strip().isdigit() else None; freq = float(ln[12:24])
            elif ln.startswith(" -4"):
                name = ln.split()[1]; cur = {} if name in want else None
                if cur is not None: out.setdefault(name, []).append((freq, cur))
            elif cur is not None and ln.startswith(" -1"):
                cur[int(ln[3:13])] = [float(ln[13 + 12 * i:25 + 12 * i]) for i in range((len(ln.rstrip()) - 13) // 12)]
            elif ln.startswith(" -3"): cur = None
    return out


z = np.load(R / "foundation_modal_basis.npz"); dvec = z["modes"][int(z["work_index"])] - z["modes"][int(z["tool_index"])]
B = read_blocks(D / "foundation_energy" / "model.frd", {"DISP", "ENER", "STRESS"})
g = pv.read(R / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"); nid = surf.point_data["node_id"]
SEL = [(7, "モード 8（182 Hz）：テーブルの送り方向"), (10, "モード 11（431 Hz）：テーブルの横方向"), (17, "モード 18（910 Hz）：主軸・テストバーの曲げ")]
tiles = []
for r, lab in SEL:
    fq, U = B["DISP"][r]; _, E = B["ENER"][r]; _, S = B["STRESS"][r]
    u = np.array([U.get(int(n), [0, 0, 0])[:3] for n in nid]); e = np.array([E.get(int(n), [0])[0] for n in nid])
    sc = 1e-3 / np.linalg.norm(dvec[:2, r])      # 刃先とワークの相対振幅（XY）が 1 µm になるようにモードを拡大縮小
    s = np.array(S_ := [S.get(int(n), [0] * 6) for n in nid]) * sc; e = e * sc**2 * 1e6   # 応力 MPa、エネルギー密度 J/m³
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2) + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    row = []
    for vals, cm, ttl, unit in ((e, "inferno", "ひずみエネルギー密度", "J/m³"), (vm, "jet", "ミーゼス応力", "MPa")):
        top = float(np.percentile(vals, 99.5))
        sf = surf.copy(); sf.points = surf.points + u / np.abs(u).max() * 40; sf.point_data["v"] = vals
        p = pv.Plotter(off_screen=True, window_size=(560, 480)); p.set_background("white")
        p.add_mesh(surf.copy(), style="wireframe", color="#cccccc", opacity=0.15)
        p.add_mesh(sf, scalars="v", cmap=cm, clim=(0, top), scalar_bar_args=dict(title=unit, color="black", n_labels=3, fmt="%.2g", vertical=True,
                   position_x=0.84, position_y=0.15, height=0.6, width=0.06, title_font_size=16, label_font_size=14))
        p.camera_position = [(2300, -1700, 1500), (500, 450, 380), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.15)
        im = Image.fromarray(p.screenshot(return_img=True)); p.close()
        t = Image.new("RGB", (560, 520), "white"); t.paste(im, (0, 40)); ImageDraw.Draw(t).text((10, 6), f"{ttl}（最大 {top:.2g} {unit}）", fill=(40, 40, 40), font=font(18)); row.append(t)
    # 部品ごとの集計（表面節点のエネルギー密度の合計の割合、目安）
    tiles.append((lab, row))
W = 560 * 2 + 360; H = 520 * len(tiles) + 60; out = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(out)
d.text((10, 12), "びびりを決めるモードで、変形のエネルギーと応力がどこに集まるか（刃先とワークの相対振幅 1 µm のとき。変形は拡大）", fill=(20, 20, 20), font=font(22))
for i, (lab, row) in enumerate(tiles):
    y = 60 + 520 * i; d.text((10, y + 220), lab.replace("：", "\n"), fill=(20, 20, 20), font=font(21))
    for j, t in enumerate(row): out.paste(t, (360 + 560 * j, y))
out.save(R / "mode_energy.png"); print("ok")
