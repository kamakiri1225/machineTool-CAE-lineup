"""ana003: 代表点＋MPC の分解図（pv_rbe2_links.png / pv_rbe3_links.png）を作る。
   図は CalculiX・FrontISTR に渡したのと同じモデルデータ（mpc_springs.build）から pyvista（ParaView と同じ VTK）で描く"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import mpc_springs, vtu_tools as vt
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
for kind in ("rbe2", "rbe3"):
    vt.render_mpc(mpc_springs.build(kind=kind), os.path.join(R, f"pv_{kind}_links.png"), kind)
print("ok")
