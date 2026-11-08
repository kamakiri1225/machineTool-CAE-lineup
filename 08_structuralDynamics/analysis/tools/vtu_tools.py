"""モデルと結果を VTU（ParaView で開ける形式）に書き出し、pyvista で図を作る"""
import numpy as np
import pyvista as pv

pv.OFF_SCREEN = True
PART_ID = {"BED": 0, "TABLE": 1}


def make_grid(model):
    """20節点六面体（VTK_QUADRATIC_HEXAHEDRON）＋ばね（線）"""
    ids = sorted(n for p in model.parts for n in [c[0] for c in p.coords])
    idx = {n: i for i, n in enumerate(ids)}
    pts = np.zeros((len(ids), 3))
    for p in model.parts:
        for n, x, y, z in p.coords:
            pts[idx[n]] = (x, y, z)
    cells, types, part = [], [], []
    for p in model.parts:
        for e in p.elems:   # 節点順は Abaqus/CalculiX/FrontISTR と VTK で同じ
            cells += [20] + [idx[n] for n in e]; types.append(pv.CellType.QUADRATIC_HEXAHEDRON); part.append(PART_ID[p.name])
    seen = set()
    for name, dof, k, nb, nt in model.springs:
        if (nb, nt) in seen:
            continue
        seen.add((nb, nt))
        cells += [2, idx[nb], idx[nt]]; types.append(pv.CellType.LINE); part.append(2)
    g = pv.UnstructuredGrid(np.array(cells), np.array(types, dtype=np.uint8), pts)
    g.cell_data["part"] = np.array(part)   # 0=ベッド・コラム, 1=テーブル, 2=ばね
    g.point_data["node_id"] = np.array(ids)
    return g, ids


def add_modes(g, ids, blocks, freqs=None):
    for k, (name, f, d) in enumerate(blocks):
        g.point_data[f"mode{k+1:02d}_{f:.1f}Hz"] = np.array([d.get(n, (0, 0, 0)) for n in ids])
    return g


def render_mode(g, field, path=None, scale=60.0, title="", plotter=None, cpos=None, show_springs=True):
    u = g.point_data[field]; mx = np.linalg.norm(u, axis=1).max()
    g = g.copy(); g.point_data["U"] = u / mx; g.point_data["|U|"] = np.linalg.norm(u, axis=1) / mx
    solid = g.extract_cells(np.where(g.cell_data["part"] < 2)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
    warped = solid.warp_by_vector("U", factor=scale)
    own = plotter is None
    p = plotter or pv.Plotter(off_screen=True, window_size=(900, 750))
    p.add_mesh(solid, style="wireframe", color="gray", opacity=0.15, line_width=0.5)
    p.add_mesh(warped, scalars="|U|", cmap="turbo", clim=(0, 1), show_edges=True, edge_color="#333333",
               line_width=0.2, scalar_bar_args=dict(title="|U| (normalized)", n_labels=3, fmt="%.1f"))
    if show_springs:
        sp = g.extract_cells(np.where(g.cell_data["part"] == 2)[0]).warp_by_vector("U", factor=scale)
        p.add_mesh(sp, color="black", line_width=3)
    if title:
        p.add_text(title, font_size=12)
    p.camera_position = cpos or [(2600, -1900, 1900), (500, 400, 400), (0, 0, 1)]
    if own:
        p.screenshot(path); p.close()


def render_springs(model, path, lift=180.0):
    """ばねの位置の分解図: テーブル側を lift だけ持ち上げ、ばねを方向別の色の線で描く"""
    g, ids = make_grid(model)
    idx = {n: i for i, n in enumerate(ids)}
    pts = g.points.copy()
    tab_nodes = [idx[c[0]] for c in model.parts[1].coords]
    pts[tab_nodes, 2] += lift
    g.points = pts
    bed = g.extract_cells(np.where(g.cell_data["part"] == 0)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
    tab = g.extract_cells(np.where(g.cell_data["part"] == 1)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
    p = pv.Plotter(off_screen=True, window_size=(1500, 1000))
    p.add_mesh(bed, color="#b0b7bf", show_edges=True, edge_color="#555555", line_width=0.3, opacity=1.0)
    p.add_mesh(tab, color="#7fb3d5", show_edges=True, edge_color="#2e4053", line_width=0.3, opacity=0.55)
    col = {1: "#c0392b", 2: "#27ae60", 3: "#2471a3"}
    lab = {1: "X spring (ball screw axial)", 2: "Y spring (guide lateral)", 3: "Z spring (guide radial)"}
    off = {1: (-4, 0), 2: (0, 0), 3: (4, 0)}
    for dof in (1, 2, 3):
        segs = []
        for name, d, k, nb, nt in model.springs:
            if d != dof:
                continue
            a = pts[idx[nb]].copy(); b = pts[idx[nt]].copy()
            a[0] += off[dof][0]; b[0] += off[dof][0]
            segs += [a, b]
        if segs:
            lines = pv.line_segments_from_points(np.array(segs))
            p.add_mesh(lines, color=col[dof], line_width=4, label=lab[dof])
    p.add_legend(bcolor="white", size=(0.42, 0.16), loc="lower right")
    for name, (x0, x1, y0, y1), kd in __import__("machine_model").INTERFACES:
        p.add_point_labels([((x0 + x1) / 2, (y0 + y1) / 2, 230 + lift + 120)], [name], font_size=16,
                           point_size=1, shape_opacity=0.7, always_visible=True)
    p.add_text(f"Exploded view: table lifted {lift:.0f} mm. Lines = springs (actually zero length at z=230)", font_size=11)
    p.camera_position = [(1900, -1500, 1500), (500, 320, 300), (0, 0, 1)]
    p.screenshot(path); p.close()


def render_mpc(model, path, kind, lift=180.0):
    """代表点＋MPC の分解図: 面の節点 → 代表点 の拘束を細線、代表点どうしのばねを太線で描く"""
    import machine_model as mm
    coord = {n: np.array((x, y, z)) for p in model.parts for n, x, y, z in p.coords}
    coord.update({n: np.array((x, y, z)) for n, x, y, z in model.extra_nodes})
    tab_ids = {c[0] for c in model.parts[1].coords}
    # 代表点は 2 個ずつ（B: ref, rot / T: ref, rot）の順で作っている
    ex = [e[0] for e in model.extra_nodes]
    tab_aux = set(ex[2::4]) | set(ex[3::4])
    lifted = lambda n: coord[n] + (np.array((0, 0, lift)) if (n in tab_ids or n in tab_aux) else 0)
    g, ids = make_grid(type("M", (), {"parts": model.parts, "springs": []})())
    idx = {n: i for i, n in enumerate(ids)}
    pts = g.points.copy(); pts[[idx[n] for n in tab_ids], 2] += lift; g.points = pts
    bed = g.extract_cells(np.where(g.cell_data["part"] == 0)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
    tab = g.extract_cells(np.where(g.cell_data["part"] == 1)[0]).extract_surface(nonlinear_subdivision=0, algorithm="dataset_surface")
    p = pv.Plotter(off_screen=True, window_size=(1500, 1000))
    p.add_mesh(bed, color="#b0b7bf", show_edges=True, edge_color="#555555", line_width=0.3)
    p.add_mesh(tab, color="#7fb3d5", show_edges=True, edge_color="#2e4053", line_width=0.3, opacity=0.45)
    refs = set(ex[0::2]); rots = set(ex[1::2])
    links = []
    for eq in model.equations:
        nodes = [t[0] for t in eq]
        r = [n for n in nodes if n in refs]
        others = [n for n in nodes if n not in refs and n not in rots]
        for a in r:
            for b in others:
                links += [lifted(a) + (0, 0, -20 if a in tab_aux else 20), lifted(b)]
    if links:
        p.add_mesh(pv.line_segments_from_points(np.array(links)), color="#7d3c98", line_width=1.2,
                   label=f"MPC ({kind.upper()})")
    sp = []
    for name, d, k, a, b in model.springs:
        if a in refs:
            sp += [lifted(a) + (0, 0, 20), lifted(b) + (0, 0, -20)]
    p.add_mesh(pv.line_segments_from_points(np.array(sp)), color="#c0392b", line_width=6, label="springs")
    rp = np.array([lifted(n) + (0, 0, -20 if n in tab_aux else 20) for n in refs])
    p.add_points(rp, color="#f1c40f", point_size=18, render_points_as_spheres=True, label="REF points")
    p.camera_position = [(1500, -900, 1000), (500, 300, 330), (0, 0, 1)]
    img = p.screenshot(return_img=True); p.close()
    # 凡例・説明は Pillow で日本語を重ねる（VTK の凡例は文字が切れるため）
    from PIL import Image, ImageDraw, ImageFont
    import os as _os
    font = lambda s: ImageFont.truetype(_os.path.expanduser("~/.fonts/meiryo.ttc"), s)
    im = Image.fromarray(img); d = ImageDraw.Draw(im)
    d.rectangle((0, 0, im.width, 44), fill="white")
    d.text((12, 6), f"{kind.upper()}：分解図（テーブルを {lift:.0f} mm 持ち上げ、代表点は見やすいよう 20 mm ずらして表示）"
           "　モデルデータは CalculiX / FrontISTR の入力と同じ", fill=(20, 20, 20), font=font(20))
    x0, y0 = im.width - 480, 60
    d.rectangle((x0, y0, im.width - 12, y0 + 130), fill="white", outline=(120, 120, 120))
    items = [((125, 60, 152), f"MPC（{kind.upper()}）：面の節点 → 代表点の拘束式"), ((192, 57, 43), "ばね：代表点どうし（並進・回転）"),
             ((241, 196, 15), "代表点（REF）")]
    for k, (c, t) in enumerate(items):
        yy = y0 + 14 + 38 * k
        d.ellipse((x0 + 14, yy + 4, x0 + 34, yy + 24), fill=c); d.text((x0 + 44, yy), t, fill=(30, 30, 30), font=font(19))
    im.save(path)
