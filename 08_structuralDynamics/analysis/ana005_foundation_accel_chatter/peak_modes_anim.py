"""ana005 5.4：①のコンプライアンスの山ごとに、どのモードか（揺れるアニメーション・ひずみエネルギー密度・ミーゼス応力）を並べる
   基礎ばねモデル・テストバー。foundation_energy/model.frd（20 モード、DISP・ENER・STRESS）を使う。
   応力・エネルギーは、そのモードの最大変位が 1 µm になる大きさにそろえて表示（値の比較用）"""
import sys
import numpy as np
from pathlib import Path
from io import BytesIO
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
from gif_tools import save_gif
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)


def read_blocks(path, want):
    out = {}; cur = None
    with open(path) as f:
        for ln in f:
            if ln.startswith("  100CL"): freq = float(ln[12:24])
            elif ln.startswith(" -4"):
                name = ln.split()[1]; cur = {} if name in want else None
                if cur is not None: out.setdefault(name, []).append((freq, cur))
            elif cur is not None and ln.startswith(" -1"):
                cur[int(ln[3:13])] = [float(ln[13 + 12 * i:25 + 12 * i]) for i in range((len(ln.rstrip()) - 13) // 12)]
            elif ln.startswith(" -3"): cur = None
    return out


B = read_blocks(D / "foundation_energy" / "model.frd", {"DISP", "ENER", "STRESS"})
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2
f = np.geomspace(5, 3000, 3000); w = 2 * np.pi * f; den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
H = {"A：刃先": (phi[it, 0] ** 2 / den).sum(1), "B：ワーク": (phi[iw, 0] ** 2 / den).sum(1), "C：相対": ((phi[iw, 0] - phi[it, 0]) ** 2 / den).sum(1)}
COL = {"A：刃先": "#1f618d", "B：ワーク": "#e67e22", "C：相対": "#c0392b"}
PEAKS = [(1, "モード 2（22 Hz）", "機械全体が脚の上で X に揺れる。A・B の山、C には出ない"),
         (7, "モード 8（182 Hz）", "テーブルが送り方向に揺れる（ボールねじ）。B・C の山"),
         (9, "モード 10（346 Hz）", "主軸頭・コラムが前後に曲がる。A・C の山"),
         (17, "モード 18（910 Hz）", "主軸とテストバーが曲がる。A・C の山")]
g = pv.read(R / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"); nid = surf.point_data["node_id"]
CAM = [(2300, -1700, 1500), (500, 450, 400), (0, 0, 1)]; TW, TH, NF = 400, 340, 16


def view(vals, cm, clim, unit, U=None):
    p = pv.Plotter(off_screen=True, window_size=(TW, TH)); p.set_background("white")
    sf = surf.copy(); sf.point_data["v"] = vals
    p.add_mesh(surf.copy(), style="wireframe", color="#cccccc", opacity=0.15)
    p.add_mesh(sf, scalars="v", cmap=cm, clim=clim, scalar_bar_args=dict(title=unit, color="black", n_labels=3, fmt="%.2g", vertical=True,
               position_x=0.83, position_y=0.12, height=0.65, width=0.06, title_font_size=13, label_font_size=12))
    p.camera_position = CAM; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.15)
    if U is None:
        im = Image.fromarray(p.screenshot(return_img=True)); p.close(); return [im] * NF
    ims = []
    for i in range(NF):
        sf.points = surf.points + U * 60 * np.sin(2 * np.pi * i / NF); p.render(); ims.append(Image.fromarray(p.screenshot(return_img=True)))
    p.close(); return ims


rows = []
for r, title, desc in PEAKS:
    fq, Ud = B["DISP"][r]; _, E = B["ENER"][r]; _, S = B["STRESS"][r]
    U = np.array([Ud.get(int(n), [0, 0, 0])[:3] for n in nid]); sc = 1e-3 / np.linalg.norm(U, axis=1).max()
    e = np.array([E.get(int(n), [0])[0] for n in nid]) * sc**2 * 1e6
    s = np.array([S.get(int(n), [0] * 6) for n in nid]) * sc
    vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2) + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
    mag = np.linalg.norm(U, axis=1)
    anim = view(mag / mag.max(), "turbo", (0, 1), "動き", U / np.abs(U).max())
    en = view(e, "inferno", (0, float(np.percentile(e, 99.5))), "J/m³"); st = view(vm, "jet", (0, float(np.percentile(vm, 99.5))), "MPa")
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for k, h in H.items(): ax.loglog(f, np.abs(h) * 1e3, color=COL[k], lw=1.6, label=k)
    ax.axvline(fq, color="k", lw=2); ax.set_xlim(5, 3000); ax.set_ylim(1e-4, 1); ax.grid(alpha=0.3, which="both")
    ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("X-X [µm/N]"); ax.legend(fontsize=8, loc="lower left"); ax.set_title(f"{fq:.0f} Hz の山", fontsize=11)
    fig.tight_layout(); buf = BytesIO(); fig.savefig(buf, format="png", dpi=85); plt.close(fig); fim = Image.open(buf).convert("RGB")
    rows.append((title, desc, fim, anim, en, st))
W = 400 + TW * 3; RH = TH + 46; frames = []
for i in range(NF):
    im = Image.new("RGB", (W, 60 + RH * len(rows)), "white"); d = ImageDraw.Draw(im)
    d.text((10, 8), "①のコンプライアンスの山はどのモードか（基礎ばねモデル・テストバー）", fill=(20, 20, 20), font=font(24))
    d.text((10, 38), "左から：FRF のどの山か ／ モードの揺れ ／ ひずみエネルギー密度 ／ ミーゼス応力（最大変位 1 µm にそろえた値）", fill=(80, 80, 80), font=font(16))
    for k, (title, desc, fim, anim, en, st) in enumerate(rows):
        y = 60 + RH * k
        d.text((10, y + 4), f"{title}：{desc}", fill=(20, 20, 20), font=font(19))
        im.paste(fim.resize((400, int(400 * fim.height / fim.width))), (0, y + 36))
        im.paste(anim[i], (400, y + 36)); im.paste(en[i], (400 + TW, y + 36)); im.paste(st[i], (400 + 2 * TW, y + 36))
    frames.append(im)
save_gif(frames, str(R / "peak_modes.gif"), duration=90); print("ok")
