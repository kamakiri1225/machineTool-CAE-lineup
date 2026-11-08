"""ana004: ParaView 風の図（pyvista = ParaView と同じ VTK）と、ParaView で開ける VTU を作る
   VTU: results/<case>.vtu（立方体 = VTK_HEXAHEDRON、ばね = VTK_LINE、点データ DISPLACEMENT）
   ParaView では  開く → Warp By Vector（DISPLACEMENT）→ Color by DISPLACEMENT  で同じ図になる"""
import os, sys
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
spec = importlib.util.spec_from_file_location("mr", os.path.join(os.path.dirname(__file__), "make_and_run.py"))
src = open(spec.origin).read().split("env = dict(os.environ)")[0]     # モデル定義だけ使う（再計算しない）
ns = {"__file__": spec.origin}; exec(src, ns)
CASES, read_ccx_dat = ns["CASES"], ns["read_ccx_dat"]

D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, "results")
pv.OFF_SCREEN = True
FONT = os.path.expanduser("~/.fonts/meiryo.ttc"); font = lambda s: ImageFont.truetype(FONT, s)
PV_BG, PV_BG2 = "#52576e", "#1f2130"          # ParaView 5 の既定背景に近い色
SCALE = 100.0


def zigzag(p0, p1, n=8, w=6.0):
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float); d = p1 - p0; L = np.linalg.norm(d)
    if L < 1e-9:
        return np.array([p0, p1])
    t = d / L; nrm = np.cross(t, [0, 1, 0] if abs(t[1]) < 0.9 else [1, 0, 0]); nrm /= np.linalg.norm(nrm)
    pts = [p0, p0 + t * L * 0.15]
    for i, s in enumerate(np.linspace(0.15, 0.85, 2 * n + 1)[1:-1]):
        pts.append(p0 + t * L * s + nrm * w * (1 if i % 2 == 0 else -1))
    pts += [p0 + t * L * 0.85, p1]
    return np.array(pts)


def make_vtu(c, U):
    ids = [n for n, *_ in c["nodes"]]; idx = {n: i for i, n in enumerate(ids)}
    pts = np.array([[x, y, z] for _, x, y, z in c["nodes"]])
    cells, types = [], []
    for e in c["elems"]:
        cells += [8] + [idx[n] for n in e]; types.append(pv.CellType.HEXAHEDRON)
    extra = []
    if c["kind"] == "ground":                    # 地面側の点（ばねの反対端、表示用）を追加
        for n in c["spr_nodes"]:
            extra.append(pts[idx[n]] + (0, 0, -60))
        for i, n in enumerate(c["spr_nodes"]):
            cells += [2, idx[n], len(ids) + i]; types.append(pv.CellType.LINE)
    else:
        for a, b in c["pairs"]:
            cells += [2, idx[a], idx[b]]; types.append(pv.CellType.LINE)
    allp = np.vstack([pts] + ([np.array(extra)] if extra else []))
    g = pv.UnstructuredGrid(np.array(cells), np.array(types, dtype=np.uint8), allp)
    disp = np.zeros((len(allp), 3)); disp[:len(ids)] = np.array([U[n] for n in ids])
    g.point_data["DISPLACEMENT"] = disp
    g.cell_data["kind"] = np.array([0 if t == pv.CellType.HEXAHEDRON else 1 for t in types])   # 0=立方体, 1=ばね
    return g, idx, len(ids)


def scene(g, nsolid, c, deformed, lift=0.0):
    """ParaView 風: 背景グラデーション、Cool to Warm、スカラーバー、座標軸"""
    p = pv.Plotter(off_screen=True, window_size=(820, 760))
    p.set_background(PV_BG, top=PV_BG2)
    gg = g.copy(); pts = gg.points.copy()
    if lift:                                     # 長さ 0 のばねを見せるため上の立方体を持ち上げる（表示だけ）
        up = np.arange(8, nsolid); pts[up, 2] += lift
    base = pts.copy()
    if deformed:
        pts += gg.point_data["DISPLACEMENT"] * SCALE
    gg.points = pts
    solid = gg.extract_cells(np.where(gg.cell_data["kind"] == 0)[0])
    mag = np.linalg.norm(solid.point_data["DISPLACEMENT"], axis=1)
    solid.point_data["DISPLACEMENT Magnitude"] = mag
    if deformed:
        g0 = g.copy(); g0.points = base
        ref = g0.extract_cells(np.where(g0.cell_data["kind"] == 0)[0])
        p.add_mesh(ref, style="wireframe", color="white", opacity=0.35, line_width=1)
        p.add_mesh(solid, scalars="DISPLACEMENT Magnitude", cmap="coolwarm", show_edges=True, edge_color="#222222",
                   clim=(0, 0.25), scalar_bar_args=dict(title="|DISPLACEMENT| [mm]", color="white", n_labels=3, fmt="%.2f",
                                                        vertical=True, position_x=0.82, position_y=0.25, height=0.5, width=0.05))
    else:
        p.add_mesh(solid, color="#d9d9d9", show_edges=True, edge_color="#333333", opacity=0.85)
    # ばね（ジグザグで描く）
    for i in np.where(gg.cell_data["kind"] == 1)[0]:
        cid = gg.get_cell(i).point_ids
        a, b = gg.points[cid[0]], gg.points[cid[1]]
        p.add_mesh(pv.lines_from_points(zigzag(b, a) if c["kind"] == "ground" else zigzag(a, b)), color="#ff4d4d", line_width=4)
    if c["kind"] == "ground":
        p.add_mesh(pv.Plane(center=(50, 50, -60), direction=(0, 0, 1), i_size=220, j_size=220), color="#8a8a8a", opacity=0.6)
    # 荷重の矢印
    for n in c["load_nodes"]:
        q = gg.points[ns_idx[n]]
        dv = np.zeros(3); dv[c["load_dof"] - 1] = np.sign(c["expect"])
        p.add_mesh(pv.Arrow(start=q - dv * 60, direction=dv, scale=55), color="#ffd24d")
    # 拘束
    fx = c["fixed"] if c["fixed"] else []
    if fx:
        p.add_points(np.array([gg.points[ns_idx[n]] for n in fx]), color="#4dd2ff", point_size=14, render_points_as_spheres=True)
    p.add_axes(color="white")
    p.camera_position = [(420, -330, 330 + (lift or 0)), (50, 50, 120 if c["kind"] != "ground" else 40), (0, 0, 1)]
    img = Image.fromarray(p.screenshot(return_img=True)); p.close()
    return img


TXT = {
    "A_ground": ("A. 地面ばね", "CalculiX: SPRING1 ＋ *SPRING（自由度 1 つ）", "FrontISTR: !SPRING（cnt に「節点, 自由度, 値」）",
                 "立方体の下面 4 節点を X・Y・Z の地面ばね（1 N/µm）で支え、上面に −Z 1000 N"),
    "B_spring2": ("B. 2 節点・方向を指定（長さ 0）", "CalculiX: SPRING2 ＋ *SPRING（「自由度, 自由度」）", "FrontISTR: 511 ＋ !SPRING_D（「自由度, 自由度, 値」）",
                  "同じ位置の節点 4 組を X・Y・Z のばねでつなぐ（表示では上を 40 mm 持ち上げ）。上の下面に +X 1000 N"),
    "C_springa": ("C. 2 節点・節点を結ぶ線の向き", "CalculiX: SPRINGA ＋ *SPRING（2 行目は空）", "FrontISTR: 511 ＋ !SPRING_A（値だけ）",
                  "50 mm 離れた節点 4 組を縦のばねでつなぐ（X・Y は拘束）。上の下面に −Z 1000 N"),
}
os.makedirs(R, exist_ok=True)
for name, c in CASES.items():
    U = read_ccx_dat(os.path.join(D, name, "ccx", "model.dat"))
    g, ns_idx, nsolid = make_vtu(c, U)
    g.save(os.path.join(R, f"{name}.vtu"))
    lift = 40.0 if c["kind"] == "dir" else 0.0
    a = scene(g, nsolid, c, False, lift); b = scene(g, nsolid, c, True, lift)
    t = TXT[name]; W, H = a.size; top = 150
    im = Image.new("RGB", (2 * W, H + top), "white"); d = ImageDraw.Draw(im)
    d.text((14, 8), t[0], fill=(20, 20, 20), font=font(30))
    d.text((14, 52), t[1], fill=(31, 97, 141), font=font(22)); d.text((14 + W, 52), t[2], fill=(192, 57, 43), font=font(22))
    d.text((14, 90), t[3], fill=(60, 60, 60), font=font(20))
    im.paste(a, (0, top)); im.paste(b, (W, top))
    d.text((12, top + 8), "設定（赤 = ばね、黄 = 荷重、水色 = 固定）", fill="white", font=font(20))
    d.text((W + 12, top + 8), f"結果（CalculiX、変形 ×{SCALE:.0f}、白線 = 元の形）", fill="white", font=font(20))
    im.save(os.path.join(R, f"pv_{name}.png")); print(name)
