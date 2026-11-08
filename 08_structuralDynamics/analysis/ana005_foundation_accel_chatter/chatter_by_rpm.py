"""ana005 5.6：回転数ごとのびびりモード（刃先―ワーク相対 FRF、基礎ばねモデル）
   各回転数の限界点 (n, a_lim, f_c) で、特性方程式 det[I + a K_t G(ω_c) B0 (1−e^{−iω_c τ})] = 0 の固有ベクトル v から
   切削力の向き F = B0 v を求め、各モードの応答 η_r = d_r·F / (ω_r² − ω_c² + 2iζ_r ω_r ω_c) を重ねて全体の振動形状を作る。
   出力：chatter_by_rpm.png（回転数ごとのモード寄与）、chatter_by_rpm.gif（回転数を動かしたときの振動形状）"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
from figstyle import plt, GRAY
from gif_tools import save_gif
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from io import BytesIO
pv.OFF_SCREEN = True
font = lambda s: ImageFont.truetype(str(Path.home() / ".fonts/meiryo.ttc"), s)
TEETH, KT, KR = 4, 1000.0, 0.3
B0 = TEETH / 4 * np.array([[KR, 1.0], [-1.0, KR]])
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2; d = phi[iw, :2, :] - phi[it, :2, :]
L = np.load(R / "stability_compare.npz"); rpm, lim, fc = L["rpm"], L["lim_C"], L["fc_C"]


def shape(n, a, f):
    w = 2 * np.pi * f; den = wr**2 - w**2 + 2j * zr * wr * w
    G = np.einsum("ir,jr,r->ij", d, d, 1 / den); tau = 60 / (TEETH * n); h = 1 - np.exp(-1j * w * tau)
    M = np.eye(2) + KT * a * (G @ B0) * h
    ev, V = np.linalg.eig(M); v = V[:, np.argmin(abs(ev))]          # 特性行列をゼロにする向き
    F = B0 @ v; eta = (d.T @ F) / den                                 # 各モードの振幅
    share = np.abs(eta * np.linalg.norm(d, axis=0)); share /= share.sum()
    return eta, share, abs(ev).min()


ok = np.isfinite(lim); S = np.zeros((len(rpm), len(fr))); ETA = np.zeros((len(rpm), len(fr)), complex)
for k in np.where(ok)[0]:
    ETA[k], S[k], res = shape(rpm[k], lim[k], fc[k])
top = np.argsort(S.sum(0))[::-1][:5]
NAME = {r: f"モード {r+1}（{fr[r]:.0f} Hz）" for r in range(len(fr))}
cols = plt.cm.tab10(np.arange(6))
fig, axs = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
axs[0].plot(rpm, lim, color="k", lw=2); axs[0].set_ylabel("切込み限界 $a_p$ [mm]"); axs[0].set_ylim(0, None); axs[0].grid(alpha=0.3)
axs[0].set_title("安定限界線図（刃先―ワーク相対 FRF、基礎ばねモデル）", loc="left", fontsize=12)
axs[1].plot(rpm, fc, color="k", lw=1.5)
for r, c in zip(top, cols): axs[1].axhline(fr[r], color=c, ls=":", lw=1.2); axs[1].text(16100, fr[r], NAME[r], color=c, fontsize=9, va="center")
axs[1].set_ylabel("びびり周波数 [Hz]"); axs[1].grid(alpha=0.3); axs[1].set_title("びびり周波数（点線 = 固有振動数）", loc="left", fontsize=12)
base = np.zeros(len(rpm))
for r, c in zip(top, cols):
    axs[2].fill_between(rpm, base, base + S[:, r], color=c, alpha=0.8, label=NAME[r]); base += S[:, r]
axs[2].fill_between(rpm, base, 1, color=GRAY, alpha=0.3, label="その他"); axs[2].set_ylim(0, 1)
axs[2].set_ylabel("びびり振動に占める割合"); axs[2].set_xlabel("主軸回転数 [min⁻¹]"); axs[2].legend(fontsize=9, loc="lower right", ncol=3); axs[2].grid(alpha=0.3)
axs[2].set_title("回転数ごとに、どのモードがびびり振動をつくっているか", loc="left", fontsize=12)
fig.tight_layout(); fig.savefig(R / "chatter_by_rpm.png"); plt.close(fig)

# アニメーション：各回転数で「どの FRF ピーク」「どのモード」がびびりをつくるかを並べる
g = pv.read(R / "vmc_model.vtu"); surf = g.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
idx = {n: i for i, n in enumerate(z["node_ids"])}; pid = np.array([idx[n] for n in surf.point_data["node_id"]])
PH = phi[pid]
fg_ = np.geomspace(50, 1500, 2000); wg = 2 * np.pi * fg_
dg = wr[None] ** 2 - wg[:, None] ** 2 + 2j * zr[None] * wr[None] * wg[:, None]
Gxx = (d[0] ** 2 / dg).sum(1); Gyy = (d[1] ** 2 / dg).sum(1)
CAM = [(2300, -1700, 1500), (500, 450, 420), (0, 0, 1)]


def plotter():
    p = pv.Plotter(off_screen=True, window_size=(430, 400)); p.set_background("white")
    sf = surf.copy(); sf.point_data["u"] = np.zeros(sf.n_points)
    p.add_mesh(surf.copy(), style="wireframe", color="#bbbbbb", opacity=0.2)
    p.add_mesh(sf, scalars="u", cmap="turbo", clim=(0, 1), show_scalar_bar=False)
    p.add_points(np.array([[500.0, 400.0, 400.0]]), color="red", point_size=10, render_points_as_spheres=True)
    p.camera_position = CAM; p.enable_parallel_projection(); p.reset_camera(); p.camera.zoom(1.05); return p, sf


def draw(p, sf, U, ph):
    U = U / np.abs(U).max(); mag = np.linalg.norm(np.abs(U), axis=1)
    sf.points = surf.points + (U * np.exp(1j * ph)).real * 60; sf.point_data["u"] = mag / mag.max(); p.render()
    return Image.fromarray(p.screenshot(return_img=True))


pc, sc = plotter(); pm = [plotter() for _ in range(2)]
sel = np.where(ok)[0][::12]; frames = []; NPH = 4
for i, k in enumerate(sel):
    o = np.argsort(S[k])[::-1][:2]
    fg, ax = plt.subplots(1, 2, figsize=(12.6, 3.7), gridspec_kw=dict(width_ratios=[1, 1.15]))
    ax[0].plot(rpm, lim, "k", lw=1.5); ax[0].plot(rpm[k], lim[k], "o", color="red", ms=9); ax[0].set_ylim(0, None); ax[0].grid(alpha=0.3)
    ax[0].set_xlabel("主軸回転数 [min⁻¹]"); ax[0].set_ylabel("$a_p$ [mm]"); ax[0].set_title(f"① 安定限界：n = {rpm[k]:.0f} min⁻¹、$a_p$ = {lim[k]:.2f} mm", loc="left", fontsize=11)
    ax[1].semilogy(fg_, np.abs(Gxx) * 1e3, color="#1f77b4", lw=1.6, label="相対 X-X"); ax[1].semilogy(fg_, np.abs(Gyy) * 1e3, color="#2ca02c", lw=1.6, label="相対 Y-Y")
    ax[1].axvline(fc[k], color="red", lw=2); ax[1].text(fc[k] + 10, 0.5, f"びびり {fc[k]:.0f} Hz", color="red", fontsize=10)
    for j, (r, c) in enumerate(zip(o, ("#d62728", "#ff7f0e"))):
        ax[1].axvline(fr[r], color=c, ls=":", lw=1.5); ax[1].text(fr[r] + 8, 0.012 * (3 ** j), f"モード {r+1}\n{fr[r]:.0f} Hz", color=c, fontsize=9)
    ax[1].set_xscale('log'); ax[1].set_ylim(0.005, 1); ax[1].set_xlim(50, 1500); ax[1].grid(alpha=0.3, which="both"); ax[1].legend(fontsize=9, loc="upper left")
    ax[1].set_xlabel("周波数 [Hz]"); ax[1].set_ylabel("µm/N"); ax[1].set_title("② 相対コンプライアンス：どのピークの近くでびびるか", loc="left", fontsize=11)
    fg.tight_layout(); buf = BytesIO(); fg.savefig(buf, format="png", dpi=90); plt.close(fg); top_im = Image.open(buf).convert("RGB")
    for j in range(NPH):
        ph = 2 * np.pi * j / NPH
        tiles = [draw(pc, sc, PH @ ETA[k], ph)] + [draw(pp, ss, PH[:, :, r], ph) for (pp, ss), r in zip(pm, o)]
        W = max(top_im.width, 430 * 3); im = Image.new("RGB", (W, 50 + top_im.height + 60 + 400), "white"); dr = ImageDraw.Draw(im)
        dr.text((10, 10), "回転数ごとのびびり：どの FRF ピーク・どのモードか（基礎ばねモデル、刃先―ワーク相対）", fill=(20, 20, 20), font=font(22))
        im.paste(top_im, (0, 50)); y = 50 + top_im.height
        labs = [("③ びびり振動の形（全モードの合計）", (200, 30, 30)), (f"④ モード {o[0]+1}（{fr[o[0]]:.0f} Hz）割合 {S[k, o[0]]*100:.0f} %", (214, 39, 40)),
                (f"⑤ モード {o[1]+1}（{fr[o[1]]:.0f} Hz）割合 {S[k, o[1]]*100:.0f} %", (255, 127, 14))]
        for t, (tile, (lb, cl)) in enumerate(zip(tiles, labs)):
            dr.text((430 * t + 10, y + 20), lb, fill=cl, font=font(18)); im.paste(tile, (430 * t, y + 60))
        frames.append(im)
pc.close(); [pp.close() for pp, _ in pm]
save_gif(frames, str(R / "chatter_by_rpm.gif"), duration=150)
md = ["| 回転数 [min⁻¹] | 切込み限界 [mm] | びびり周波数 [Hz] | 1 番目のモード | 2 番目のモード |", "|---:|---:|---:|---|---|"]
for n0 in (1000, 2000, 3000, 4000, 6000, 8000, 10000, 12000, 14000, 16000):
    k = int(np.argmin(abs(rpm - n0))); o = np.argsort(S[k])[::-1]
    md.append(f"| {rpm[k]:.0f} | {lim[k]:.2f} | {fc[k]:.0f} | {NAME[o[0]]} {S[k,o[0]]*100:.0f} % | {NAME[o[1]]} {S[k,o[1]]*100:.0f} % |")
(R / "chatter_by_rpm.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
