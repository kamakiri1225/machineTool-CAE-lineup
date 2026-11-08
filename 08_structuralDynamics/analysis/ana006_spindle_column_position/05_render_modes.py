"""最小の切込み限界の点で効いているモード（割合 1 位）の形を描く → results/modes_at_min.png、results/mode_base.gif"""
import sys, json
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; sys.path.insert(0, str(D)); R = D / "results"
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import params as P
from lib.fem import build
from lib.modal import read_frd
from lib.view import grid
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
S = json.loads((R / "summary.json").read_text())
SHOW = [("base", "基準"), ("z300", "主軸頭 +300 mm"), ("ym100", "コラム −100 mm（手前）"), ("yp050", "コラム +50 mm（奥）")]
CAM = [(-1700, -2300, 2100), (500, 560, 600), (0, 0, 1)]


def frame(m, g, ids, U, scale, out, phase=1.0):
    pts0 = np.array([m.coords[n] for n in ids]); u = np.array([U.get(n, (0, 0, 0)) for n in ids])
    gg = g.copy(); gg.points = pts0 + scale * phase * u; gg.point_data["|u|"] = np.linalg.norm(u, axis=1) / np.abs(u).max()
    pl = pv.Plotter(window_size=(900, 900)); pl.set_background("white")
    pl.add_mesh(g.extract_surface(algorithm="dataset_surface"), style="wireframe", color="#bbbbbb", line_width=0.3, opacity=0.25)
    pl.add_mesh(gg.extract_surface(algorithm="dataset_surface"), scalars="|u|", cmap="turbo", clim=(0, 1), show_scalar_bar=False, show_edges=False)
    pl.camera_position = CAM; pl.camera.zoom(1.05); pl.screenshot(out); pl.close()


tiles = []
for name, lab in SHOW:
    if name not in S:
        continue
    dy, dz = P.CASES[name]; m = build(dy, dz); g, ids = grid(m); top = S[name]["top_modes"][0]
    fr, blocks = read_frd(D / "cases" / name / "model.frd", modes={top["mode"] - 1}); U = blocks[top["mode"] - 1]
    scale = 120 / max(np.linalg.norm(v) for v in U.values())
    p = R / f"_mode_{name}.png"; frame(m, g, ids, U, scale, str(p)); im = Image.open(p).convert("RGB"); p.unlink()
    d = ImageDraw.Draw(im); d.text((20, 15), f"{lab}：モード {top['mode']}（{top['f_Hz']:.0f} Hz）", fill="black", font=font(36))
    d.text((20, 65), f"最小 a_lim = {S[name]['amin_mm']:.2f} mm のときの揺れの {top['share']*100:.0f} %", fill="#333333", font=font(30)); tiles.append(im)
    if name == "base":
        frames = []
        for ph in np.sin(np.linspace(0, 2 * np.pi, 24, endpoint=False)):
            q = R / "_f.png"; frame(m, g, ids, U, scale, str(q), ph); fi = Image.open(q).convert("RGB").resize((560, 560), Image.LANCZOS); q.unlink()
            ImageDraw.Draw(fi).text((12, 10), f"基準：モード {top['mode']}（{top['f_Hz']:.0f} Hz）", fill="black", font=font(24)); frames.append(fi)
        frames[0].save(R / "mode_base.gif", save_all=True, append_images=frames[1:], duration=80, loop=0)
W, H = tiles[0].size; canvas = Image.new("RGB", (W * len(tiles), H + 70), "white")
for i, im in enumerate(tiles):
    canvas.paste(im, (i * W, 0))
ImageDraw.Draw(canvas).text((20, H + 10), "色 = 揺れの大きさ（各図の最大を 1 として、青 0 → 赤 1）。灰色の線 = 揺れる前の形。変形は見やすいよう拡大", fill="black", font=font(32))
canvas = canvas.resize((canvas.width // 2, canvas.height // 2), Image.LANCZOS); canvas.save(R / "modes_at_min.png"); print(canvas.size)
