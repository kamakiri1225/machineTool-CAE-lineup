"""ばね接続モデルを自分で作る手順（docs/ana003_ばね接続モデルを自分で作る手順.md）の図を作る。
   tutorial/ccx/model.inp（実際に解いた入力）を読み、手順ごとの図（ParaView 風、pyvista）と ParaView 用の VTU を出す。
   実行: ~/opt/venv/bin/python make_tutorial_figs.py   （先に make_tutorial.py で入力を作って解いておく）"""
import os, sys, re
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from frd import read_frd

D = os.path.dirname(os.path.abspath(__file__)); T = os.path.join(D, "tutorial")
OUT = os.path.join(D, "results", "tutorial"); os.makedirs(OUT, exist_ok=True)
VTU = os.path.join(T, "vtu"); os.makedirs(VTU, exist_ok=True)
pv.OFF_SCREEN = True
FONT = os.path.expanduser("~/.fonts/meiryo.ttc"); font = lambda s: ImageFont.truetype(FONT, s)

# ------------------------------------------------------------ 入力ファイルを読む
txt = open(os.path.join(T, "ccx", "model.inp")).read()
blocks = re.split(r"\n(?=\*)", txt)
nodes, elems, nsets, eqs, springs = {}, [], {}, [], []
for b in blocks:
    lines = [l for l in b.strip().split("\n") if l.strip()]
    key = lines[0].upper()
    if key.startswith("*NODE") and "PRINT" not in key and "FILE" not in key:
        for l in lines[1:]:
            p = [float(v) for v in l.split(",")]; nodes[int(p[0])] = np.array(p[1:4])
    elif key.startswith("*ELEMENT") and "C3D8" in key:
        for l in lines[1:]:
            p = [int(v) for v in l.split(",")]; elems.append(p[1:])
    elif key.startswith("*NSET"):
        nsets[re.search(r"NSET=(\w+)", key).group(1)] = [int(v) for v in ",".join(lines[1:]).split(",") if v.strip()]
    elif key.startswith("*EQUATION"):
        p = [v.strip() for v in lines[2].split(",")]
        eqs.append([(int(p[i]), int(p[i + 1]), float(p[i + 2])) for i in range(0, len(p), 3)])
    elif key.startswith("*ELEMENT") and "SPRING2" in key:
        p = [int(v) for v in lines[1].split(",")]; springs.append([p[1], p[2], None, None])
    elif key.startswith("*SPRING"):
        springs[-1][2] = int(lines[1].split(",")[0]); springs[-1][3] = float(lines[2])
REF_B, ROT_B, REF_T, ROT_T = 101, 102, 201, 202
U = {}
try:
    _, bl = read_frd(os.path.join(T, "ccx", "model.frd")); U = {n: np.array(v) for n, v in bl[-1][2].items()}
except Exception as e:
    print("frd なし:", e)

# ------------------------------------------------------------ 表示用の座標（分解図）
LIFT = 90.0
UPPER = set(range(9, 17))


def disp(n):
    """表示位置: 上の直方体を LIFT 持ち上げ、代表点・補助点はずらして見せる"""
    p = nodes[n].copy()
    if n in UPPER: p[2] += LIFT
    if n == REF_B: p[2] += 15
    if n == REF_T: p[2] += LIFT - 15
    if n == ROT_B: p += (0, -38, 15)
    if n == ROT_T: p += (0, -38, LIFT - 15)
    return p


COL = dict(FACE_B="#e74c3c", FACE_T="#3498db", FIX="#48c9b0", REF="#f4d03f", ROT="#58d68d", MPC="#c39bd3", SPR="#ff5252", SPRR="#82e0aa")


def zigzag(a, b, n=7, w=7.0):
    a, b = np.asarray(a, float), np.asarray(b, float); d = b - a; L = np.linalg.norm(d); t = d / L
    nrm = np.cross(t, [1, 0, 0]); nrm /= np.linalg.norm(nrm)
    pts = [a, a + t * L * 0.15]
    for i, s in enumerate(np.linspace(0.15, 0.85, 2 * n + 1)[1:-1]):
        pts.append(a + t * L * s + nrm * w * (1 if i % 2 == 0 else -1))
    return np.array(pts + [a + t * L * 0.85, b])


def blocks_mesh(shift=None, scale=0.0):
    pts = np.array([disp(n) + (U.get(n, 0) * scale if shift else 0) for n in range(1, 17)])
    cells = []
    for e in elems:
        cells += [8] + [n - 1 for n in e]
    g = pv.UnstructuredGrid(np.array(cells), np.array([pv.CellType.HEXAHEDRON] * len(elems), dtype=np.uint8), pts)
    g.cell_data["part"] = np.array([0, 1])
    return g


def base_plotter():
    p = pv.Plotter(off_screen=True, window_size=(1100, 860)); p.set_background("#52576e", top="#1f2130")
    p.camera_position = [(470, -420, 400), (50, 50, 145), (0, 0, 1)]
    p.add_axes(color="white")
    return p


def add_blocks(p, opacity=0.35, deformed=False, scale=0.0, color_parts=True):
    g = blocks_mesh(deformed, scale)
    for k, c in ((0, "#b2babb"), (1, "#aed6f1")):
        p.add_mesh(g.extract_cells([k]), color=c if color_parts else "#cccccc", opacity=opacity, show_edges=True,
                   edge_color="#ffffff", line_width=1.5)


def label_nodes(p, ns, color, size=16, offset=(0, 0, 0)):
    p.add_point_labels(np.array([disp(n) + offset for n in ns]), [str(n) for n in ns], font_size=size, text_color=color,
                       shape_color="#1f2130", shape_opacity=0.6, point_color=color, point_size=14, render_points_as_spheres=True,
                       always_visible=True)


def finish(p, path, title, lines, legend=None):
    img = Image.fromarray(p.screenshot(return_img=True)); p.close()
    lines = lines[:1]                       # 画像の中は短く。詳しい説明はブログ本文に書く
    top = 46 + 34 * len(lines); W, H = img.size
    im = Image.new("RGB", (W, H + top), "white"); im.paste(img, (0, top)); d = ImageDraw.Draw(im)
    d.text((14, 6), title, fill=(20, 20, 20), font=font(26))
    for k, (t, c) in enumerate(lines):
        d.text((14, 46 + 34 * k), t, fill=c, font=font(19))
    if legend:
        x0, y0 = W - 330, top + 14
        d.rectangle((x0, y0, W - 12, y0 + 14 + 32 * len(legend)), fill=(255, 255, 255), outline=(150, 150, 150))
        for k, (c, t) in enumerate(legend):
            yy = y0 + 10 + 32 * k
            d.ellipse((x0 + 12, yy + 3, x0 + 32, yy + 23), fill=c); d.text((x0 + 42, yy), t, fill=(30, 30, 30), font=font(17))
    im.save(path); return im


NAV = (60, 60, 60); RED = (192, 57, 43); BLUE = (31, 97, 141)
# ------------------------------------------------------------ 手順 1：2 つの部品を別々の節点で作る
p = base_plotter(); add_blocks(p, 0.45)
label_nodes(p, [1, 2, 3, 4], "#ffffff"); label_nodes(p, [5, 6, 7, 8], COL["FACE_B"]); label_nodes(p, [9, 10, 11, 12], COL["FACE_T"])
label_nodes(p, [13, 14, 15, 16], "#ffffff")
finish(p, os.path.join(OUT, "step1_parts.png"), "手順 1：2 つの部品を、節点を共有せずに作る",
       [("5〜8 と 9〜12 は同じ座標の「別の節点」（図は上を 90 mm 持ち上げて表示）", NAV),
        (f"図では見やすいように上の直方体を {LIFT:.0f} mm 持ち上げている（実際は z = 100 で接している）", RED),
        ("節点をマージ（共有）すると 2 部品が直接くっつき、ばねを入れる意味がなくなる", BLUE)],
       [(COL["FACE_B"], "下の部品の上面の節点"), (COL["FACE_T"], "上の部品の下面の節点")])

# ------------------------------------------------------------ 手順 2：面上の節点集合
p = base_plotter(); add_blocks(p, 0.25)
for s, c in (("FACE_B", COL["FACE_B"]), ("FACE_T", COL["FACE_T"]), ("FIX", COL["FIX"])):
    pts = np.array([disp(n) for n in nsets[s]])
    p.add_mesh(pv.PolyData(np.vstack([pts, pts[:1]])).delaunay_2d() if False else pv.Polygon(), opacity=0) if False else None
    p.add_mesh(pv.lines_from_points(np.vstack([pts, pts[:1]])), color=c, line_width=5)
    label_nodes(p, nsets[s], c, 18)
finish(p, os.path.join(OUT, "step2_sets.png"), "手順 2：つなぐ面の節点を集める（節点集合）",
       [("FACE_B = 5〜8（下の上面）、FACE_T = 9〜12（上の下面）、FIX = 1〜4（底面）", NAV),
        ("選び方：「その部品の節点」かつ「z = 面の高さ」かつ「x, y が取付面の範囲内」。座標だけで選ぶと両部品が混ざる", RED),
        ("CalculiX: *NSET, NSET=FACE_B　／　FrontISTR: !NGROUP, NGRP=FACE_B（集合を作っただけではまだつながっていない）", BLUE)],
       [(COL["FACE_B"], "FACE_B"), (COL["FACE_T"], "FACE_T"), (COL["FIX"], "FIX（固定）")])

# ------------------------------------------------------------ 手順 3：代表点・回転用補助点
p = base_plotter(); add_blocks(p, 0.25)
for s, c in (("FACE_B", COL["FACE_B"]), ("FACE_T", COL["FACE_T"])):
    p.add_points(np.array([disp(n) for n in nsets[s]]), color=c, point_size=14, render_points_as_spheres=True)
label_nodes(p, [REF_B, REF_T], COL["REF"], 20); label_nodes(p, [ROT_B, ROT_T], COL["ROT"], 18)
for a, b in ((REF_B, ROT_B), (REF_T, ROT_T)):
    p.add_mesh(pv.Line(disp(a), disp(b)), color="#ffffff", line_width=1)
finish(p, os.path.join(OUT, "step3_ref.png"), "手順 3：面の中心に代表点（REF）と回転用の補助点（ROT）を足す",
       [("101・201 = 代表点（並進）、102・202 = 補助点（回転角）。本当は 4 点とも同じ座標", NAV),
        ("4 点とも本当は同じ座標 (50, 50, 100)。図ではずらして描いている。要素は付けない（節点を *NODE / !NODE に足すだけ）", RED),
        ("下の面用（101, 102）と上の面用（201, 202）は必ず別の番号にする", BLUE)],
       [(COL["REF"], "代表点 REF 101 / 201"), (COL["ROT"], "回転用補助点 ROT 102 / 202")])

# ------------------------------------------------------------ 手順 4：MPC（拘束式）で面の節点を代表点に結ぶ
p = base_plotter(); add_blocks(p, 0.2)
segs = []
for eq in eqs:
    dep = eq[0][0]
    ref = [n for n, _, _ in eq if n in (REF_B, REF_T)][0]
    segs += [disp(dep), disp(ref)]
p.add_mesh(pv.line_segments_from_points(np.array(segs)), color=COL["MPC"], line_width=3)
for s, c in (("FACE_B", COL["FACE_B"]), ("FACE_T", COL["FACE_T"])):
    p.add_points(np.array([disp(n) for n in nsets[s]]), color=c, point_size=14, render_points_as_spheres=True)
p.add_points(np.array([disp(REF_B), disp(REF_T)]), color=COL["REF"], point_size=22, render_points_as_spheres=True)
p.add_points(np.array([disp(ROT_B), disp(ROT_T)]), color=COL["ROT"], point_size=16, render_points_as_spheres=True)
p.add_mesh(pv.Arrow(disp(REF_B), disp(5) - disp(REF_B), scale="auto", shaft_radius=0.02, tip_radius=0.05), color="#ff8a65")
label_nodes(p, [5], COL["FACE_B"], 20)
finish(p, os.path.join(OUT, "step4_mpc.png"), f"手順 4：面の節点を代表点に拘束式（MPC）で結ぶ（{len(eqs)} 式）",
       [("紫の線 1 本 = 拘束式 3 本（X・Y・Z）。橙の矢印 = 節点 5 の位置ベクトル r", NAV),
        ("例：節点 5 は代表点から r = (−50, −50) → X の式「5, 1, 1, 101, 1, −1, 102, 3, −50」（橙の矢印が r）", RED),
        ("面の節点 4 個 × 3 方向 × 2 面 = 24 式。CalculiX: *EQUATION　FrontISTR: !EQUATION（.msh に書く）", BLUE)],
       [(COL["MPC"], "拘束式（MPC）"), (COL["REF"], "代表点 REF"), (COL["ROT"], "補助点 ROT")])

# ------------------------------------------------------------ 手順 5：代表点どうしをばねでつなぐ
p = base_plotter(); add_blocks(p, 0.15)
p.add_mesh(pv.line_segments_from_points(np.array(segs)), color=COL["MPC"], line_width=1.5, opacity=0.5)
p.add_mesh(pv.lines_from_points(zigzag(disp(REF_B), disp(REF_T))), color=COL["SPR"], line_width=6)
p.add_mesh(pv.lines_from_points(zigzag(disp(ROT_B), disp(ROT_T), n=5, w=5)), color=COL["SPRR"], line_width=5)
label_nodes(p, [REF_B, REF_T], COL["REF"], 18); label_nodes(p, [ROT_B, ROT_T], COL["ROT"], 16)
finish(p, os.path.join(OUT, "step5_spring.png"), "手順 5：代表点どうしをばねでつなぐ",
       [("赤 = 並進ばね（101―201、X・Y・Z）、緑 = 回転ばね（102―202、θx・θy・θz）", NAV),
        ("緑 = 回転ばね：102 ― 202 を θx・θy・θz の 3 本（各 1,000,000 N·mm/rad）。ガイドのモーメント剛性はここ", RED),
        ("CalculiX: SPRING2 ＋ *SPRING「1, 1 / 1000.0」　FrontISTR: 511 要素（.msh）＋ !SPRING_D「1, 1, 1000.0」（.cnt）", BLUE)],
       [(COL["SPR"], "並進ばね（REF ― REF）"), (COL["SPRR"], "回転ばね（ROT ― ROT）")])

# ------------------------------------------------------------ 手順 6：固定と荷重
p = base_plotter(); add_blocks(p, 0.3)
p.add_mesh(pv.lines_from_points(zigzag(disp(REF_B), disp(REF_T))), color=COL["SPR"], line_width=5)
p.add_points(np.array([disp(n) for n in nsets["FIX"]]), color=COL["FIX"], point_size=20, render_points_as_spheres=True)
for n in nsets["FIX"]:
    p.add_mesh(pv.Cone(center=disp(n) - (0, 0, 12), direction=(0, 0, 1), height=20, radius=9), color=COL["FIX"])
p.add_mesh(pv.Arrow(start=disp(13) - (0, 75, 0), direction=(0, 1, 0), scale=65), color="#f39c12")
label_nodes(p, [13], "#f39c12", 20)
finish(p, os.path.join(OUT, "step6_bc.png"), "手順 6：固定と荷重を入れる",
       [("底面 FIX を固定し、上のブロックの節点 13（上面の角）に Y 方向 1 N", NAV),
        ("上の代表点 201 に Y 方向 1 N（橙の矢印）。CalculiX: *CLOAD「201, 2, 1.0」　FrontISTR: !CLOAD「201, 2, 1.0」", RED),
        ("代表点に力をかけると、拘束式を通して上の面の 4 節点に配られる", BLUE)],
       [(COL["FIX"], "固定（FIX）"), ("#f39c12", "荷重 1 N（節点 13、Y）")])

# ------------------------------------------------------------ 手順 7：解いて確かめる
if U:
    SC = 3000.0
    p = base_plotter()
    g0 = blocks_mesh(); p.add_mesh(g0, style="wireframe", color="#ffffff", opacity=0.4)
    g = blocks_mesh(True, SC)
    mag = np.array([np.linalg.norm(U.get(n, np.zeros(3))) for n in range(1, 17)])
    g.point_data["|U| [µm]"] = mag * 1e3
    p.add_mesh(g, scalars="|U| [µm]", cmap="coolwarm", clim=(0, float(mag.max() * 1e3)), show_edges=True, edge_color="#222222", opacity=0.85,
               scalar_bar_args=dict(title="|U| [um]", color="white", n_labels=3, fmt="%.1f", vertical=True, position_x=0.86, position_y=0.2, height=0.5, width=0.05))
    a = disp(REF_B) + U[REF_B] * SC; b = disp(REF_T) + U[REF_T] * SC
    p.add_mesh(pv.lines_from_points(zigzag(a, b)), color=COL["SPR"], line_width=5)
    du = U[REF_T][1] - U[REF_B][1]
    finish(p, os.path.join(OUT, "step7_result.png"), "手順 7：解いて、ばねの両端の変位の差を確かめる（CalculiX の結果）",
           [(f"ばねの伸び u_Y(201) − u_Y(101) = {du:.6f} mm（手計算 F/k = 0.001）。上のブロックは傾く。変形 {SC:.0f} 倍", NAV),
            (f"u_Y(201) − u_Y(101) = {du:.6f} mm　←　手計算 F / k = 1 / 1000 = 0.001 mm と一致", RED),
            ("FrontISTR でも同じ（0.000999999967 mm）。上の部品の絶対変位には下の部品のわずかな変形も含まれるので、差で比べる", BLUE)])

# ------------------------------------------------------------ 全体図（1 枚目）
p = base_plotter(); add_blocks(p, 0.3)
p.add_mesh(pv.line_segments_from_points(np.array(segs)), color=COL["MPC"], line_width=2.5)
p.add_mesh(pv.lines_from_points(zigzag(disp(REF_B), disp(REF_T))), color=COL["SPR"], line_width=6)
p.add_mesh(pv.lines_from_points(zigzag(disp(ROT_B), disp(ROT_T), n=5, w=5)), color=COL["SPRR"], line_width=4)
for s, c in (("FACE_B", COL["FACE_B"]), ("FACE_T", COL["FACE_T"]), ("FIX", COL["FIX"])):
    p.add_points(np.array([disp(n) for n in nsets[s]]), color=c, point_size=14, render_points_as_spheres=True)
p.add_points(np.array([disp(REF_B), disp(REF_T)]), color=COL["REF"], point_size=22, render_points_as_spheres=True)
p.add_points(np.array([disp(ROT_B), disp(ROT_T)]), color=COL["ROT"], point_size=16, render_points_as_spheres=True)
finish(p, os.path.join(OUT, "step0_goal.png"), "完成形：面 →（MPC）→ 代表点 ―（ばね）― 代表点 ←（MPC）← 面",
       [("手順 1〜7 で 1 つずつ作る（分解図：上の部品を持ち上げて表示）", NAV)],
       [(COL["FACE_B"], "面の節点（下）"), (COL["FACE_T"], "面の節点（上）"), (COL["MPC"], "拘束式 MPC"), (COL["REF"], "代表点 REF"),
        (COL["ROT"], "補助点 ROT"), (COL["SPR"], "並進ばね"), (COL["SPRR"], "回転ばね"), (COL["FIX"], "固定")])

# ------------------------------------------------------------ ParaView 用 VTU（実座標と分解図の 2 つ）
for name, posf in (("tutorial_model.vtu", lambda n: nodes[n]), ("tutorial_model_exploded.vtu", disp)):
    ids = sorted(nodes); idx = {n: i for i, n in enumerate(ids)}
    pts = np.array([posf(n) for n in ids])
    cells, types, kind = [], [], []
    for e in elems:
        cells += [8] + [idx[n] for n in e]; types.append(pv.CellType.HEXAHEDRON); kind.append(0)
    for eq in eqs:
        dep = eq[0][0]; ref = [n for n, _, _ in eq if n in (REF_B, REF_T)][0]
        cells += [2, idx[dep], idx[ref]]; types.append(pv.CellType.LINE); kind.append(1)
    for a, b, dof, k in springs:
        cells += [2, idx[a], idx[b]]; types.append(pv.CellType.LINE); kind.append(2 if a == REF_B else 3)
    g = pv.UnstructuredGrid(np.array(cells), np.array(types, dtype=np.uint8), pts)
    g.cell_data["kind"] = np.array(kind)            # 0=ソリッド, 1=MPC, 2=並進ばね, 3=回転ばね
    setid = np.zeros(len(ids), int)
    for s, v in (("FIX", 3), ("FACE_B", 1), ("FACE_T", 2)):
        for n in nsets[s]: setid[idx[n]] = v
    setid[idx[REF_B]] = setid[idx[REF_T]] = 4; setid[idx[ROT_B]] = setid[idx[ROT_T]] = 5
    g.point_data["set"] = setid                     # 0=その他, 1=FACE_B, 2=FACE_T, 3=FIX, 4=REF, 5=ROT
    g.point_data["node_id"] = np.array(ids)
    du = np.array([U.get(n, np.zeros(3)) for n in ids]); du[[idx[ROT_B], idx[ROT_T]]] = 0   # ROT の「変位」は回転角なので表示では 0
    g.point_data["DISPLACEMENT"] = du
    g.save(os.path.join(VTU, name))
print("ok", len(eqs), "equations", len(springs), "springs")
