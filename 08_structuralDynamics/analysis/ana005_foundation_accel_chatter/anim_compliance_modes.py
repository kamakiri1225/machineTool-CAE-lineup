"""ana005 3.4a：相対コンプライアンスのピークを作る 8 モードを揺れるアニメーション（GIF）にする"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent / "tools")]
from gif_tools import save_gif
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; idx = {n: i for i, n in enumerate(z["node_ids"])}
rows = [l.split("|")[1:9] for l in (R / "compliance_modes.md").read_text().splitlines()[2:]]
sel = [(int(r[0]), float(r[1]), r[2].strip(), r[4].strip(), r[5].strip()) for r in rows]
g = pv.read(R / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"); pid = np.array([idx[n] for n in surf.point_data["node_id"]])
NF, TW, TH = 20, 470, 420
cols = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189), (140, 86, 75), (227, 119, 194), (127, 127, 127)]
views = []
for no, f, parts, xx, yy in sel:
    r = int(np.argmin(abs(fr - f))); U = phi[pid, :, r]; U = U / np.abs(U).max()
    sf = surf.copy(); sf.point_data["u"] = np.linalg.norm(U, axis=1) / np.linalg.norm(U, axis=1).max()
    p = pv.Plotter(off_screen=True, window_size=(TW, TH)); p.set_background("white")
    p.add_mesh(surf.copy(), style="wireframe", color="#bbbbbb", opacity=0.2); p.add_mesh(sf, scalars="u", cmap="turbo", clim=(0, 1), show_scalar_bar=False)
    p.add_points(np.array([[500.0, 400.0, 400.0]]), color="red", point_size=11, render_points_as_spheres=True)
    p.camera_position = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.15)
    imgs = []
    for i in range(NF):
        sf.points = surf.points + U * 60 * np.sin(2 * np.pi * i / NF); p.render(); imgs.append(Image.fromarray(p.screenshot(return_img=True)))
    p.close(); views.append((no, f, parts, xx, yy, imgs))
frames = []
for i in range(NF):
    im = Image.new("RGB", (TW * 4, 60 + (TH + 90) * 2), "white"); d = ImageDraw.Draw(im)
    d.text((10, 12), "相対コンプライアンスのピークを作るモード（基礎ばねモデル、変形は拡大。色 = 動きの大きさ、赤点 = 切削点）", fill=(20, 20, 20), font=font(24))
    for k, (no, f, parts, xx, yy, imgs) in enumerate(views):
        x, y = TW * (k % 4), 60 + (TH + 90) * (k // 4)
        d.text((x + 10, y), f"モード {no}：{f:.0f} Hz", fill=cols[k], font=font(24)); d.text((x + 10, y + 34), f"大きく動く：{parts}", fill=(40, 40, 40), font=font(17))
        d.text((x + 10, y + 60), f"寄与 X-X {float(xx):.3f}・Y-Y {float(yy):.3f} µm/N", fill=(90, 90, 90), font=font(16)); im.paste(imgs[i], (x, y + 90))
    frames.append(im)
save_gif(frames, str(R / "compliance_mode_shapes.gif"), duration=80); print("ok")
