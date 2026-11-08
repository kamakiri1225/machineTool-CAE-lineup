"""ana005 5.3：線図上の 1 点（約 6,000 min⁻¹）のびびりを、FRF のピークと揺れるモードで見せるアニメーション
   上段：安定限界線図（3 ケース）＋選んだ点、相対コンプライアンス（X-X・Y-Y）にびびり周波数と効いているモードを表示
   下段：びびり振動の形（全モードの合計）と、割合の大きい 2 モードを同じ位相で揺らす"""
import sys
import numpy as np
from pathlib import Path
from io import BytesIO
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
import chatter_by_rpm as C          # ETA, S, rpm, lim, fc, fr, PH, surf などを再利用
from figstyle import plt
from gif_tools import save_gif
import pyvista as pv
from PIL import Image, ImageDraw
k = int(np.nanargmin(C.lim)); o = np.argsort(C.S[k])[::-1][:2]
fig, ax = plt.subplots(1, 2, figsize=(14, 4.2), gridspec_kw=dict(width_ratios=[1, 1.2]))
for name, lab, col in (("fixed", "底面固定", "#1f77b4"), ("foundation", "基礎ばね", "#ff7f0e"), ("foundation_table_stiff", "基礎ばね・テーブル結合10倍", "#2ca02c")):
    a = np.loadtxt(R / f"milling_{name}.csv", delimiter=","); ax[0].plot(a[:, 0], a[:, 1], color=col, lw=1.5, label=lab)
ax[0].plot(C.rpm[k], C.lim[k], "*", color="red", ms=16); ax[0].text(C.rpm[k] + 300, C.lim[k] - 0.6, f"この点\nn = {C.rpm[k]:.0f} min⁻¹\n$a_p$ = {C.lim[k]:.2f} mm", color="red", fontsize=10)
ax[0].set_ylim(0, 30); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8.5, loc="lower right"); ax[0].set_xlabel("主軸回転数 [min⁻¹]"); ax[0].set_ylabel("切込み限界 $a_p$ [mm]")
ax[0].set_title("① 安定限界線図", loc="left", fontsize=12)
ax[1].semilogy(C.fg_, np.abs(C.Gxx) * 1e3, color="#1f77b4", lw=1.8, label="相対 X-X"); ax[1].semilogy(C.fg_, np.abs(C.Gyy) * 1e3, color="#2ca02c", lw=1.8, label="相対 Y-Y")
ax[1].axvline(C.fc[k], color="red", lw=2.2); ax[1].text(C.fc[k] * 1.04, 0.75, f"びびり {C.fc[k]:.0f} Hz", color="red", fontsize=10.5)
for j, (r, c) in enumerate(zip(o, ("#d62728", "#ff7f0e"))):
    yy = np.abs((C.Gxx if j == 0 else C.Gyy))[np.argmin(abs(C.fg_ - C.fr[r]))] * 1e3
    ax[1].plot(C.fr[r], yy, "o", color=c, ms=10); ax[1].annotate(f"モード {r+1}（{C.fr[r]:.0f} Hz）", (C.fr[r], yy), xytext=((C.fr[r] * 1.6, 0.4) if j == 0 else (C.fr[r] * 1.6, 0.02)), color=c, fontsize=10, arrowprops=dict(arrowstyle="->", color=c))
ax[1].set_xscale('log'); ax[1].set_xlim(50, 1500); ax[1].set_ylim(0.005, 1); ax[1].grid(alpha=0.3, which="both"); ax[1].legend(fontsize=9, loc="upper left")
ax[1].set_xlabel("周波数 [Hz]"); ax[1].set_ylabel("µm/N"); ax[1].set_title("② 相対コンプライアンスのどのピークでびびるか", loc="left", fontsize=12)
fig.tight_layout(); buf = BytesIO(); fig.savefig(buf, format="png", dpi=95); plt.close(fig); top = Image.open(buf).convert("RGB")
CAM = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]


def make():
    p = pv.Plotter(off_screen=True, window_size=(440, 420)); p.set_background("white")
    sf = C.surf.copy(); sf.point_data["u"] = np.zeros(sf.n_points)
    p.add_mesh(C.surf.copy(), style="wireframe", color="#999999", opacity=0.25)
    p.add_mesh(sf, scalars="u", cmap="turbo", clim=(0, 1), show_scalar_bar=False, show_edges=False)
    p.add_points(np.array([[500.0, 400.0, 400.0]]), color="red", point_size=12, render_points_as_spheres=True)
    p.camera_position = CAM; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.05); return p, sf


views = [(C.PH @ C.ETA[k], "③ びびり振動の形（全モードの合計）", (200, 30, 30)),
         (C.PH[:, :, o[0]], f"④ モード {o[0]+1}（{C.fr[o[0]]:.0f} Hz）、割合 {C.S[k,o[0]]*100:.0f} %", (214, 39, 40)),
         (C.PH[:, :, o[1]], f"⑤ モード {o[1]+1}（{C.fr[o[1]]:.0f} Hz）、割合 {C.S[k,o[1]]*100:.0f} %", (230, 110, 0))]
plots = [make() for _ in views]; frames = []
for i in range(16):
    ph = 2 * np.pi * i / 16; tiles = []
    for (p, sf), (U, lab, cl) in zip(plots, views):
        U = U / np.abs(U).max(); mag = np.linalg.norm(np.abs(U), axis=1)
        sf.points = C.surf.points + (U * np.exp(1j * ph)).real * 60; sf.point_data["u"] = mag / mag.max(); p.render()
        tiles.append(Image.fromarray(p.screenshot(return_img=True)))
    W = max(top.width, 440 * 3); im = Image.new("RGB", (W, 50 + top.height + 50 + 420), "white"); dr = ImageDraw.Draw(im)
    dr.text((10, 10), "立型MC（基礎ばねモデル）：線図の 1 点のびびりは、どの FRF ピーク・どのモードか（色 = 動きの大きさ、赤点 = 切削点）", fill=(20, 20, 20), font=C.font(19))
    im.paste(top, (0, 50)); y = 50 + top.height
    for t, (tile, (_, lab, cl)) in enumerate(zip(tiles, views)):
        dr.text((440 * t + 8, y + 15), lab, fill=cl, font=C.font(15)); im.paste(tile, (440 * t, y + 50))
    frames.append(im)
[p.close() for p, _ in plots]
save_gif(frames, str(R / "milling_with_modes.gif"), duration=90); print("ok")
