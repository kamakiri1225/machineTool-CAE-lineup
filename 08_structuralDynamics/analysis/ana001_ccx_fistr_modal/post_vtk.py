"""ana001: 結果を VTU に書き出し、pyvista（VTK）で図を作る"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm
from frd import read_frd
import vtu_tools as vt
import pyvista as pv

D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, "results"); os.makedirs(R, exist_ok=True)
model = mm.build()
g, ids = vt.make_grid(model)
coords, blocks = read_frd(os.path.join(D, "ccx", "machine.frd"))
vt.add_modes(g, ids, blocks)
g.save(os.path.join(R, "ccx_modes.vtu"))            # ParaView でそのまま開ける

NAMES = ["Column bends in Y", "Table moves in X (ball screw)", "Column bends in X",
         "Table moves in Y (guides)", "Table twists about Z", "Table moves in Z + bends"]
fields = [k for k in g.point_data.keys() if k.startswith("mode")]

# 6モードを 1枚に
p = pv.Plotter(off_screen=True, shape=(2, 3), window_size=(1800, 1200), border=False)
for k in range(6):
    p.subplot(k // 3, k % 3)
    f = float(fields[k].split("_")[1][:-2])
    vt.render_mode(g, fields[k], plotter=p, title=f"Mode {k+1}: {f:.0f} Hz\n{NAMES[k]}")
p.screenshot(os.path.join(R, "pv_mode_shapes.png")); p.close()

# ばね周りの拡大（テーブル X モード、ばねを見せるためにテーブルを横から）
p = pv.Plotter(off_screen=True, window_size=(1200, 700))
vt.render_mode(g, fields[1], plotter=p, scale=40, title="Mode 2 close-up: springs (black lines) between rails and carriages",
               cpos=[(500, -900, 420), (500, 300, 250), (0, 0, 1)])
p.screenshot(os.path.join(R, "pv_mode2_closeup.png")); p.close()

# FrontISTR の結果（同じ 2次モード）を並べて比較
fis = pv.read(os.path.join(D, "fistr", "machine_eigen_vis_psf.0002.pvtu"))
fis.cell_data["part"] = np.where(fis.celltypes == pv.CellType.HEXAHEDRON, 0, 2)
fis.point_data["mode02"] = fis.point_data["DISPLACEMENT"]
p = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1600, 700), border=False)
p.subplot(0, 0); vt.render_mode(g, fields[1], plotter=p, title=f"CalculiX  mode 2: {blocks[1][1]:.1f} Hz")
p.subplot(0, 1); vt.render_mode(fis, "mode02", plotter=p, title="FrontISTR  mode 2: 181.2 Hz")
p.screenshot(os.path.join(R, "pv_ccx_vs_fistr_mode2.png")); p.close()
vt.render_springs(model, os.path.join(R, "pv_spring_locations.png"))
print("saved")
