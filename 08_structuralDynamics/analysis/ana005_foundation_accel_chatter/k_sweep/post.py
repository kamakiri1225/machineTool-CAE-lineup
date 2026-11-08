"""ana005 4.2c：脚の横剛性 1/4・1・4 倍の全体 FEM の結果をまとめる（固有振動数・コンプライアンス・モード形状）"""
import sys, re, json
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; P = D.parent
sys.path[:0] = [str(P.parent / "tools"), str(P.parent.parent / "study" / "scripts")]
from figstyle import plt, BLUE, RED, GRAY, GREEN
from frd import read_frd
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
CASES = [("k025", "脚の横剛性 1/4（2.5 N/µm/脚）", RED), ("k1", "今（10 N/µm/脚）", BLUE), ("k4", "4 倍（40 N/µm/脚）", GREEN)]


def eig(path):
    out = []; on = False
    for ln in open(path):
        if "E I G E N V A L U E   O U T P U T" in ln: on = True
        elif on and "P A R T I C I P" in ln: break
        elif on:
            p = ln.split()
            if len(p) == 5 and p[0].isdigit(): out.append(float(p[3]))
    return np.array(out)


def frf(path):
    lines = open(path).read().splitlines(); blk = []
    for i, ln in enumerate(lines):
        m = re.match(r"\s*displacements \(vx,vy,vz\) for set (\w+) and time\s+(\S+)", ln)
        if m: blk.append((m.group(1), float(m.group(2)), np.array([float(v) for v in lines[i + 2].split()[1:4]])))
    f, data = [], {}
    for k in range(0, len(blk), 8):
        f.append(blk[k][1])
        for j in range(4):
            data.setdefault(blk[k + j][0], []).append(blk[k + j][2] + 1j * blk[k + 4 + j][2])
    return np.array(f), {k: np.array(v) for k, v in data.items()}


res = {}
for c, lab, col in CASES:
    fr = eig(D / c / "model.dat"); f, d = frf(D / c / "model.dat"); res[c] = (fr, f, d)
    print(c, np.round(fr[:6], 1))
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
for c, lab, col in CASES:
    fr, f, d = res[c]
    axs[0].semilogy(f, np.abs(d["BEDPT"][:, 0]) * 1e3, color=col, lw=2, label=lab)
    axs[1].semilogy(f, np.abs(d["COLTOP"][:, 0]) * 1e3, color=col, lw=2, label=lab)
    axs[2].semilogy(f, np.abs(d["WORK"][:, 0] - d["TOOL"][:, 0]) * 1e3, color=col, lw=2, label=lab)
    for a in axs: a.axvline(fr[0], color=col, ls=":", lw=1)
for a, t in zip(axs, ("① ベッド側面の X（本体全体の揺れ）", "② コラム上端の X", "③ 刃先―ワークの相対 X（加工誤差になる量）")):
    a.set_title(t, loc="left", fontsize=12); a.set_xlabel("周波数 [Hz]"); a.grid(alpha=0.3, which="both"); a.set_xlim(1, 300)
axs[0].set_ylabel("駆動反力 1 N あたりの変位 [µm/N]"); axs[0].legend(fontsize=10)
fig.suptitle("全体 FEM（CalculiX）：ナットの駆動力（テーブル +1 N・ベッド −1 N、X）に対するコンプライアンス。点線 = 1 次固有振動数", fontsize=13)
fig.tight_layout(); fig.savefig(D / "k_sweep_frf.png"); plt.close(fig)
# モード形状
g = pv.read(P / "results" / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface"); nid = surf.point_data["node_id"]
rows = []
for c, lab, col in CASES:
    _, blocks = read_frd(str(D / c / "eig" / "model.frd"))
    tiles = []
    for r in range(3):
        name, fq, dd = blocks[r]; U = np.array([dd.get(int(n), (0, 0, 0)) for n in nid]); s = 80 / np.abs(U).max()
        sf = surf.copy(); sf.points = surf.points + U * s; sf.point_data["u"] = np.linalg.norm(U, axis=1) / np.linalg.norm(U, axis=1).max()
        p = pv.Plotter(off_screen=True, window_size=(460, 400)); p.set_background("white")
        p.add_mesh(surf, style="wireframe", color="#bbbbbb", opacity=0.2); p.add_mesh(sf, scalars="u", cmap="turbo", clim=(0, 1), show_scalar_bar=False)
        p.camera_position = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.1)
        im = Image.fromarray(p.screenshot(return_img=True)); p.close()
        t = Image.new("RGB", (460, 450), "white"); t.paste(im, (0, 50)); ImageDraw.Draw(t).text((10, 8), f"{r+1} 次 {fq:.1f} Hz", fill=tuple(int(255 * v) for v in plt.matplotlib.colors.to_rgb(col)), font=font(24))
        tiles.append(t)
    rows.append((lab, tiles))
W = 460 * 3 + 330; sheet = Image.new("RGB", (W, 450 * 3 + 60), "white"); dr = ImageDraw.Draw(sheet)
dr.text((10, 10), "モード形状（全体 FEM、変形は拡大）：1・2 次 = 機械全体が脚の上で X・Y に揺れる、3 次 = Z 軸まわりにねじれる", fill=(20, 20, 20), font=font(24))
for i, (lab, tiles) in enumerate(rows):
    dr.text((10, 60 + 450 * i + 200), lab.replace("（", "\n（"), fill=(20, 20, 20), font=font(22))
    for j, t in enumerate(tiles): sheet.paste(t, (330 + 460 * j, 60 + 450 * i))
sheet.save(D / "k_sweep_modes.png")
md = ["| ケース | 1 次 Hz | 2 次 Hz | 3 次 Hz | 1 次の比（今 = 1） | 1 自由度モデルの f_n |", "|---|---:|---:|---:|---:|---:|"]
for (c, lab, col), fn1 in zip(CASES, (11.5, 22.9, 45.8)):
    fr = res[c][0]; md.append(f"| {lab} | {fr[0]:.1f} | {fr[1]:.1f} | {fr[2]:.1f} | {fr[0]/res['k1'][0][0]:.2f} | {fn1} |")
(D / "k_sweep.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
