"""ana005 3.6：どのコンプライアンス（相対 FRF のピーク）が、どのモード（揺れ方）から来ているか
   基礎ばねモデルの 80 モード（質量正規化、results/foundation_modal_basis.npz）で、刃先―ワーク相対 FRF をモードの和に分解する
     G_ij(ω) = Σ_r d_ri d_rj / (ω_r² − ω² + 2i ζ_r ω_r ω),  d_r = φ_r(ワーク) − φ_r(刃先)
   出力: results/compliance_modes.png（FRF と寄与）、results/compliance_mode_shapes.png（モード形状の一覧）、results/compliance_modes.md"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
pv.OFF_SCREEN = True
FONT = str(Path.home() / ".fonts/meiryo.ttc"); font = lambda s: ImageFont.truetype(FONT, s)
PART = {0: "コラム", 1: "テーブル", 2: "主軸", 3: "テストバー", 4: "ワーク", 5: "ベッド", 6: "主軸頭"}
RHO = {0: 7.3e-9, 1: 7.3e-9, 2: 7.8e-9, 3: 14.5e-9, 4: 7.3e-9, 5: 7.3e-9, 6: 7.3e-9}

d = np.load(R / "foundation_modal_basis.npz")
fr, phi, iw, it = d["frequency_hz"], d["modes"], int(d["work_index"]), int(d["tool_index"])
wr = 2 * np.pi * fr; zr = d["alpha"] / (2 * wr) + d["beta"] * wr / 2
dvec = phi[iw] - phi[it]                                  # (3, nmodes) 相対モードベクトル
f = np.geomspace(5, 3000, 4000); w = 2 * np.pi * f
den = wr[None, :] ** 2 - w[:, None] ** 2 + 2j * zr[None, :] * wr[None, :] * w[:, None]
contrib = {(i, j): dvec[i][None, :] * dvec[j][None, :] / den for i in (0, 1) for j in (0, 1)}   # (nw, nm)
G = {k: v.sum(1) for k, v in contrib.items()}
# CalculiX の直接計算と照合
chk = []
for j, name in ((0, "X"), (1, "Y")):
    c = np.loadtxt(R / f"foundation_relative_{j+1}.csv", delimiter=",")
    Hc = c[:, 1 + j] + 1j * c[:, 4 + j]
    Hm = np.interp(c[:, 0], f, G[(j, j)].real) + 1j * np.interp(c[:, 0], f, G[(j, j)].imag)
    chk.append(np.max(np.abs(Hm - Hc)) / np.max(np.abs(Hc)))
print("modal sum vs CalculiX max rel diff:", [f"{x:.2e}" for x in chk])

# 寄与の大きいモード（X-X と Y-Y のピーク寄与）
peak = {k: np.abs(dvec[k[0]] * dvec[k[1]]) / (2 * zr * wr**2) for k in ((0, 0), (1, 1))}
sel = sorted(set(list(np.argsort(peak[(0, 0)])[::-1][:4]) + list(np.argsort(peak[(1, 1)])[::-1][:4])))
# どの部品が動くか：モードの運動エネルギーの部品別の割合
g = pv.read(R / "vmc_model.vtu"); nid = g.point_data["node_id"]; comp = np.array(g.cell_data["component"]).copy()
cc = g.cell_centers().points
head = (cc[:, 0] > 380) & (cc[:, 0] < 620) & (cc[:, 1] > 340) & (cc[:, 1] < 750) & (cc[:, 2] > 650) & (cc[:, 2] < 850)
comp[(comp == 0) & head] = 6; comp[(comp == 0) & (cc[:, 2] < 230)] = 5      # 鋳物を主軸頭・ベッド・コラムに分ける
idx = {n: i for i, n in enumerate(d["node_ids"])}; pidx = np.array([idx[n] for n in nid])
vol = g.compute_cell_sizes(length=False, area=False)["Volume"]
cells = [g.get_cell(c).point_ids for c in range(g.n_cells)]


cpts = {k: np.unique(np.concatenate([cells[c] for c in np.where(comp == k)[0]])) for k in PART if (comp == k).any()}


def maxmove(r):
    """部品ごとの最大の動き（モード全体の最大を 1）"""
    u = np.linalg.norm(phi[pidx, :, r], axis=1); return {k: u[v].max() / u.max() for k, v in cpts.items()}


def share(r):
    u2 = (phi[pidx, :, r] ** 2).sum(1)
    e = np.array([u2[cells[c]].mean() * vol[c] * RHO[comp[c]] for c in range(g.n_cells)])
    return {k: e[comp == k].sum() / e.sum() for k in PART}


rows = []
for n, r in enumerate(sel, 1):
    sh = share(r); mv = maxmove(r); main = max(mv, key=mv.get)
    rows.append(dict(no=r + 1, r=r, f=fr[r], xx=peak[(0, 0)][r] * 1e3, yy=peak[(1, 1)][r] * 1e3, share=sh, move=mv, main=main,
                     dw=np.linalg.norm(phi[iw, :2, r]), dt=np.linalg.norm(phi[it, :2, r])))
    print(n, f"{fr[r]:.1f} Hz", PART[main], {PART[k]: f"{v*100:.0f}%" for k, v in sh.items() if v > 0.05}, f"XX {peak[(0,0)][r]*1e3:.4f} YY {peak[(1,1)][r]*1e3:.4f}")

# ---------------------------------------------------------------- 図 1：FRF と寄与
cols = plt.cm.tab10(np.arange(len(sel)))
fig, axs = plt.subplots(2, 2, figsize=(16, 9.5), gridspec_kw=dict(width_ratios=[1.6, 1]))
for row, (k, lab) in enumerate((((0, 0), "X-X"), ((1, 1), "Y-Y"))):
    ax = axs[row, 0]
    ax.loglog(f, np.abs(G[k]) * 1e3, color="k", lw=2.4, label="合計（80 モード）")
    for c, rw in zip(cols, rows):
        ax.loglog(f, np.abs(contrib[k][:, rw["r"]]) * 1e3, "--", color=c, lw=1.4)
        if rw["xx" if k == (0, 0) else "yy"] > 0.02 * max(r_["xx" if k == (0, 0) else "yy"] for r_ in rows):
            y = np.abs(G[k])[np.argmin(abs(f - rw["f"]))] * 1e3
            ax.annotate(f"{rw['no']}", (rw["f"], y), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=12,
                        color=c, weight="bold")
    ax.set_xlim(5, 3000); ax.set_ylim(1e-4, 1); ax.grid(alpha=0.3, which="both")
    ax.set_ylabel(f"相対コンプライアンス {lab} [µm/N]"); ax.set_title(f"刃先―ワーク相対 FRF {lab}（基礎ばねモデル）と各モードの寄与（破線）", loc="left", fontsize=12)
    if row == 1: ax.set_xlabel("周波数 [Hz]")
    ax = axs[row, 1]
    v = [rw["xx" if k == (0, 0) else "yy"] for rw in rows]
    ax.barh(range(len(rows)), v, color=cols); ax.invert_yaxis()
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{rw['no']}: {rw['f']:.0f} Hz {PART[rw['main']]}" for rw in rows], fontsize=10)
    ax.set_xlabel(f"ピークでの寄与 {lab} [µm/N]"); ax.grid(alpha=0.3, axis="x")
    ax.set_title(f"どのモードが {lab} のピークを作るか", loc="left", fontsize=12)
fig.suptitle(f"コンプライアンスをモードに分解（モード重ね合わせ。CalculiX の直接計算との差 {max(chk)*100:.2f} % 以下）", fontsize=13)
fig.tight_layout(); fig.savefig(R / "compliance_modes.png"); plt.close(fig)

# ---------------------------------------------------------------- 図 2：モード形状の一覧
surf_all = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
sid = surf_all.point_data["node_id"]; spidx = np.array([idx[n] for n in sid])
tiles = []
for c, rw in zip(cols, rows):
    U = phi[spidx, :, rw["r"]]; s = 60.0 / np.abs(U).max()
    sf = surf_all.copy(); sf.points = surf_all.points + U * s; sf.point_data["|u|"] = np.linalg.norm(U, axis=1) / np.abs(U).max()
    p = pv.Plotter(off_screen=True, window_size=(520, 460)); p.set_background("white")
    p.add_mesh(surf_all, style="wireframe", color="#bbbbbb", opacity=0.25)
    p.add_mesh(sf, scalars="|u|", cmap="turbo", clim=(0, 1), show_scalar_bar=False)
    p.add_points(np.array([[500.0, 400.0, 400.0]]) + 0, color="red", point_size=12, render_points_as_spheres=True)
    p.camera_position = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.15)
    im = Image.fromarray(p.screenshot(return_img=True)); p.close()
    W, H = im.size; t = Image.new("RGB", (W, H + 110), "white"); t.paste(im, (0, 110)); dr = ImageDraw.Draw(t)
    rgb = tuple(int(255 * v) for v in c[:3])
    dr.text((10, 4), f"モード {rw['no']}：{rw['f']:.1f} Hz", fill=rgb, font=font(26))
    top = [k for k, v in sorted(rw["move"].items(), key=lambda kv: -kv[1]) if v > 0.5][:3]
    dr.text((10, 40), "大きく動く部品：" + "・".join(PART[k] for k in top), fill=(30, 30, 30), font=font(17))
    dr.text((10, 70), f"寄与 X-X {rw['xx']:.3f}・Y-Y {rw['yy']:.3f} µm/N", fill=(80, 80, 80), font=font(17))
    tiles.append(t)
ncol = 4; nrow = int(np.ceil(len(tiles) / ncol)); W, H = tiles[0].size
sheet = Image.new("RGB", (W * ncol, H * nrow + 70), "white"); dr = ImageDraw.Draw(sheet)
dr.text((12, 10), "各ピークのモード形状（基礎ばねモデル、変形は拡大。色＝動きの大きさ、赤点＝切削点）", fill=(20, 20, 20), font=font(26))
for k, t in enumerate(tiles):
    sheet.paste(t, ((k % ncol) * W, 70 + (k // ncol) * H))
sheet.save(R / "compliance_mode_shapes.png")

# ---------------------------------------------------------------- 図 3：3 つの支持条件での比較（CalculiX の直接計算）
fig, axs = plt.subplots(1, 2, figsize=(16, 5))
for j, lab in ((0, "X-X"), (1, "Y-Y")):
    ax = axs[j]
    for case, col, name in (("fixed", GRAY, "底面固定"), ("foundation", BLUE, "基礎ばね"), ("foundation_table_stiff", RED, "基礎ばね＋テーブル結合 10 倍")):
        c = np.loadtxt(R / f"{case}_relative_{j+1}.csv", delimiter=",")
        ax.loglog(c[:, 0], np.hypot(c[:, 1 + j], c[:, 4 + j]) * 1e3, color=col, lw=2 if case == "foundation" else 1.6, label=name)
    for c_, rw in zip(cols, rows):
        ax.axvline(rw["f"], color=c_, lw=0.8, ls=":"); ax.text(rw["f"], 1.3, str(rw["no"]), color=c_, ha="center", fontsize=10, weight="bold")
    ax.set_xlim(5, 3000); ax.set_ylim(1e-4, 2); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=10, loc="lower left")
    ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel(f"相対コンプライアンス {lab} [µm/N]")
    ax.set_title(f"{lab}：支持条件でピークがどう変わるか（番号 = 図のモード、点線は基礎ばねモデルの固有振動数）", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(R / "compliance_cases.png"); plt.close(fig)

md = ["| モード番号 | 固有振動数 Hz | 大きく動く部品（最大の動きの 50 % 以上） | 運動エネルギーの多い部品 | X-X 寄与 µm/N | Y-Y 寄与 µm/N | ワーク側の動き | 刃先側の動き |", "|---|---:|---|---|---:|---:|---:|---:|"]
for rw in rows:
    top = sorted(rw["share"].items(), key=lambda kv: -kv[1])[:2]
    mvk = [PART[k] for k, v in sorted(rw["move"].items(), key=lambda kv: -kv[1]) if v > 0.5][:3]
    md.append(f"| {rw['no']} | {rw['f']:.1f} | " + "・".join(mvk) + " | " + "、".join(f"{PART[k]} {v*100:.0f}%" for k, v in top) +
              f" | {rw['xx']:.4f} | {rw['yy']:.4f} | {rw['dw']:.3g} | {rw['dt']:.3g} |")
(R / "compliance_modes.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
