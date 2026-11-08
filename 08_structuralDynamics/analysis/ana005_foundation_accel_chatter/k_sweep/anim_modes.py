"""ana005 4.2c：脚の横剛性 1/4・1・4 倍の 1〜3 次モードを揺れるアニメーション（GIF）にする"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; P = D.parent
sys.path[:0] = [str(P.parent / "tools")]
from frd import read_frd
from gif_tools import save_gif
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
CASES = [("k025", "脚の横剛性 1/4", (192, 57, 43)), ("k1", "今", (31, 97, 141)), ("k4", "脚の横剛性 4 倍", (30, 132, 73))]
NF, TW, TH, LW = 20, 420, 380, 230
g = pv.read(P / "results" / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"); nid = surf.point_data["node_id"]
views = {}
for c, lab, col in CASES:
    _, blocks = read_frd(str(D / c / "eig" / "model.frd"))
    for r in range(3):
        name, fq, dd = blocks[r]; U = np.array([dd.get(int(n), (0, 0, 0)) for n in nid]); U = U / np.abs(U).max()
        sf = surf.copy(); sf.point_data["u"] = np.linalg.norm(U, axis=1) / np.linalg.norm(U, axis=1).max()
        p = pv.Plotter(off_screen=True, window_size=(TW, TH)); p.set_background("white")
        p.add_mesh(surf.copy(), style="wireframe", color="#bbbbbb", opacity=0.25)
        p.add_mesh(sf, scalars="u", cmap="turbo", clim=(0, 1), show_scalar_bar=False)
        p.camera_position = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.05)
        imgs = []
        for i in range(NF):
            sf.points = surf.points + U * 90 * np.sin(2 * np.pi * i / NF); p.render(); imgs.append(Image.fromarray(p.screenshot(return_img=True)))
        p.close(); views[(c, r)] = (fq, imgs)
frames = []
for i in range(NF):
    im = Image.new("RGB", (LW + TW * 3, 60 + (TH + 40) * 3), "white"); d = ImageDraw.Draw(im)
    d.text((10, 12), "脚の横剛性を変えたときの 1〜3 次モード（全体 FEM、変形は拡大。色 = 動きの大きさ）", fill=(20, 20, 20), font=font(24))
    for k, (c, lab, col) in enumerate(CASES):
        y = 60 + (TH + 40) * k
        d.text((10, y + TH // 2), lab, fill=col, font=font(22))
        for r in range(3):
            fq, imgs = views[(c, r)]; x = LW + TW * r
            d.text((x + 10, y + 4), f"{r+1} 次  {fq:.1f} Hz", fill=col, font=font(22)); im.paste(imgs[i], (x, y + 40))
    frames.append(im)
save_gif(frames, str(D / "k_sweep_modes.gif"), duration=80); print("ok")
