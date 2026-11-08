"""ana003: MPC（拘束式）で RBE2 を作る仕組みの説明図
  1) mpc_formula_to_input.png : 数式 → 入力ファイルの 1 行への書き換え（CAR1 の面、実際の節点番号・座標）
  2) gif/mpc_rbe2_car1.gif     : CalculiX の RBE2 モデルの 5次モード（ねじれ）で、CAR1 の面・代表点・拘束線が一緒に動く様子（ParaView 風）"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm, mpc_springs
from frd import read_frd
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, "results"); os.makedirs(os.path.join(R, "gif"), exist_ok=True)
FONT = os.path.expanduser("~/.fonts/meiryo.ttc"); font = lambda s: ImageFont.truetype(FONT, s)
m = mpc_springs.build("rbe2")
coord = {n: np.array((x, y, z)) for p in m.parts for n, x, y, z in p.coords}
coord.update({n: np.array((x, y, z)) for n, x, y, z in m.extra_nodes})
REF_B, ROT_B, REF_T, ROT_T = [e[0] for e in m.extra_nodes[:4]]       # CAR1: ベッド側 REF/ROT、テーブル側 REF/ROT
x0, x1, y0, y1 = mm.INTERFACES[0][1]; xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
face_B = [n for n, x, y in mpc_springs._patch_nodes(m.parts[0], mm.INTERFACES[0][1], mm.Z_IF)]
face_T = [n for n, x, y in mpc_springs._patch_nodes(m.parts[1], mm.INTERFACES[0][1], mm.Z_IF)]
N0 = 2912

# ============================================================ 1) 数式 → 入力ファイル
from figstyle_local import plt
fig = plt.figure(figsize=(16, 8.2))
ax = fig.add_axes([0.03, 0.08, 0.36, 0.8])
P = np.array([coord[n][:2] for n in face_B])
ax.scatter(P[:, 0], P[:, 1], s=40, color="#7d3c98", zorder=3, label="面の節点（ベッド側 CAR1 の面、45 個）")
for p in P:
    ax.plot([xc, p[0]], [yc, p[1]], color="#7d3c98", lw=0.5, alpha=0.5)
ax.scatter([xc], [yc], s=260, color="#f1c40f", edgecolor="k", zorder=4, label=f"代表点 REF（節点 {REF_B}）")
p0 = coord[N0][:2]
ax.scatter([p0[0]], [p0[1]], s=180, color="#c0392b", edgecolor="k", zorder=5)
ax.annotate("", xy=p0, xytext=(xc, yc), arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.5))
ax.text(p0[0] + 2, p0[1] - 5, f"節点 {N0}\n$r = (r_x, r_y) = ({p0[0]-xc:.0f},\\ {p0[1]-yc:.0f})$ mm", color="#c0392b", fontsize=12, va="top")
th = np.linspace(0.3, 1.4, 30); ax.plot(xc + 14 * np.cos(th), yc + 14 * np.sin(th), color="#1f618d", lw=2)
ax.annotate("", xy=(xc + 14 * np.cos(1.45), yc + 14 * np.sin(1.45)), xytext=(xc + 14 * np.cos(1.3), yc + 14 * np.sin(1.3)),
            arrowprops=dict(arrowstyle="-|>", color="#1f618d", lw=2))
ax.text(xc + 12, yc + 12, "$\\theta_z$（ROT 節点\n" f"{ROT_B} の 3 番目の自由度）", color="#1f618d", fontsize=11)
ax.set_aspect("equal"); ax.set_xlim(x0 - 15, x1 + 25); ax.set_ylim(y0 - 30, y1 + 22)
ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.grid(alpha=0.3); ax.legend(loc="upper left", fontsize=10)
ax.set_title("上から見た CAR1 の面（実際のモデルの節点）", loc="left", fontsize=13)

ax = fig.add_axes([0.42, 0.02, 0.57, 0.92]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ry = p0[1] - yc
T = [(9.5, "① 剛体の動き（RBE2 の意味）：面の点は、代表点と一緒に「平行移動＋回転」する", "k", 13),
     (8.75, "$u_{x}(i) = u_{x}(\\mathrm{REF}) - \\theta_z\\, r_y$　（$i$ = 面の節点、$r$ = REF から節点 $i$ へのベクトル）", "k", 17),
     (7.85, f"② この節点は $r_y = {ry:.0f}$ mm なので代入する", "k", 13),
     (7.1, f"$u_x({N0}) = u_x({REF_B}) - \\theta_z \\times ({ry:.0f})$", "k", 17),
     (6.2, "③ 全部を左辺に集めて「… = 0」の形にする（ソルバの線形拘束式の形）", "k", 13),
     (5.45, f"$1 \\cdot u_x({N0})\\ -\\ 1 \\cdot u_x({REF_B})\\ {'-' if ry < 0 else '+'}\\ {abs(ry):.0f} \\cdot \\theta_z = 0$", "k", 17),
     (4.5, "④ 項ごとに「節点, 自由度, 係数」を並べると、入力ファイルの 1 行になる", "k", 13)]
for y, t, c, fs in T:
    ax.text(0.2, y, t, fontsize=fs, color=c)
cols = ["#c0392b", "#d4ac0d", "#1f618d"]
terms = [(f"{N0}, 1, 1", f"節点 {N0} の\n自由度 1（$u_x$）\n係数 +1"), (f"{REF_B}, 1, -1", f"代表点 REF の\n自由度 1（$u_x$）\n係数 −1"),
         (f"{ROT_B}, 3, {ry:.0f}" if False else f"{ROT_B}, 3, -20", f"ROT 節点の自由度 3\n＝ $\\theta_z$ として使う\n係数 $-20$")]
for k, ((code, ex), c) in enumerate(zip(terms, cols)):
    x = 0.3 + 3.2 * k
    ax.add_patch(plt.Rectangle((x, 3.15), 2.9, 0.75, fc=c, alpha=0.2, ec=c, lw=2))
    ax.text(x + 1.45, 3.52, code, ha="center", va="center", fontsize=15, family="monospace")
    ax.text(x + 1.45, 2.95, ex, ha="center", va="top", fontsize=11, color=c)
ax.text(0.2, 1.05, "CalculiX:   *EQUATION / 3 / " + f"{N0}, 1, 1, {REF_B}, 1, -1, {ROT_B}, 3, -20", fontsize=12, family="monospace", color="#1f618d")
ax.text(0.2, 0.45, "FrontISTR:  !EQUATION / 3, 0.0 / " + f"{N0}, 1, 1, {REF_B}, 1, -1, {ROT_B}, 3, -20", fontsize=12, family="monospace", color="#c0392b")
ax.text(0.2, 0.0, "（/ は改行。3 = 項の数、0.0 = 右辺。面の 45 節点 × 3 方向 × 10 面 = 1,110 行を Python で自動生成）", fontsize=11, color="#555")
fig.savefig(os.path.join(R, "mpc_formula_to_input.png")); plt.close(fig)
print("fig1")

# ============================================================ 2) ParaView 風アニメーション（5次モード、CAR1 周り）
_, blocks = read_frd(os.path.join(D, "rbe2", "ccx", "machine.frd"))
name, f5, d5 = blocks[4]
LIFT = 60.0
import vtu_tools as vt
g, ids = vt.make_grid(type("M", (), {"parts": m.parts, "springs": []})())
idx = {n: i for i, n in enumerate(ids)}
U = np.array([d5.get(n, (0, 0, 0)) for n in ids]); umax = np.abs(U).max(); SC = 25.0 / umax
tab_ids = {c[0] for c in m.parts[1].coords}
base = g.points.copy(); base[[idx[n] for n in tab_ids], 2] += LIFT
box = (x0 - 70, x1 + 70, y0 - 60, y1 + 60, 150, 330)
pv.OFF_SCREEN = True


# 表示する範囲は「変形前」の形で 1 回だけ決める（毎コマ切り直すと、境界の要素が出たり消えたりしてちらつく）
gb = g.copy(); gb.points = base; gb.point_data["gid"] = np.arange(gb.n_points)
cen = gb.cell_centers().points
inside = np.where((cen[:, 0] > box[0]) & (cen[:, 0] < box[1]) & (cen[:, 1] > box[2]) & (cen[:, 1] < box[3]) &
                  (cen[:, 2] > box[4]) & (cen[:, 2] < box[5]))[0]
sub = gb.extract_cells(inside)
surf0 = sub.extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
orig_s = np.asarray(surf0.point_data["gid"])
surf0.point_data["|U|"] = np.linalg.norm(U[orig_s], axis=1) / umax


def frame(s):
    surf = surf0.copy(); surf.points = base[orig_s] + U[orig_s] * SC * s
    p = pv.Plotter(off_screen=True, window_size=(900, 720)); p.set_background("#52576e", top="#1f2130")
    p.add_mesh(surf, scalars="|U|", cmap="coolwarm", clim=(0, 1), show_edges=True, edge_color="#222222", opacity=0.6,
               scalar_bar_args=dict(title="|U| (norm.)", color="white", n_labels=3, fmt="%.1f"))
    def pos(n, side):
        u = np.array(d5.get(n, (0, 0, 0))) * SC * s
        return coord[n] + (0, 0, LIFT if side == "T" else 0) + u
    for side, ref, face in (("B", REF_B, face_B), ("T", REF_T, face_T)):
        rp = pos(ref, side) + (0, 0, 15 if side == "B" else -15)
        segs = []
        for n in face:
            segs += [rp, pos(n, side)]
        p.add_mesh(pv.line_segments_from_points(np.array(segs)), color="#bb8fce", line_width=2)
        p.add_points(rp[None], color="#f1c40f", point_size=20, render_points_as_spheres=True)
        p.add_points(np.array([pos(n, side) for n in face]), color="#d2b4de", point_size=7, render_points_as_spheres=True)
    p.add_mesh(pv.Line(pos(REF_B, "B") + (0, 0, 15), pos(REF_T, "T") + (0, 0, -15)), color="#ff4d4d", line_width=7)
    p.add_axes(color="white")
    p.camera_position = [(xc + 330, yc - 380, 420), (xc, yc, 255), (0, 0, 1)]
    img = Image.fromarray(p.screenshot(return_img=True)); p.close(); return img


frames = []
for i in range(20):
    im = frame(np.sin(2 * np.pi * i / 20)); W, H = im.size
    c = Image.new("RGB", (W, H + 96), "white"); c.paste(im, (0, 96)); d = ImageDraw.Draw(c)
    d.text((12, 6), f"RBE2（CalculiX）5次 {f5:.0f} Hz：CAR1 のまわりを拡大（テーブル側を {LIFT:.0f} mm 持ち上げた分解図）", fill=(20, 20, 20), font=font(19))
    d.text((12, 36), "紫の線 = 拘束式（*EQUATION）。面の節点は代表点（黄）と一緒に剛体として動く", fill=(125, 60, 152), font=font(18))
    d.text((12, 64), "赤 = 代表点どうしのばね（SPRING2）。変形は拡大表示", fill=(192, 57, 43), font=font(18))
    frames.append(c)
from gif_tools import save_gif
save_gif(frames, os.path.join(R, "gif", "mpc_rbe2_car1.gif"), duration=80)
frames[5].save(os.path.join(R, "mpc_rbe2_car1.png"))
print("gif")
