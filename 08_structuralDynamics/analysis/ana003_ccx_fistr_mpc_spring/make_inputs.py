"""ana003: ばねのつなぎ方の比較（分布ばね vs 代表点＋MPC の RBE2 / RBE3）固有値解析"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm
import mpc_springs as ms

D = os.path.dirname(os.path.abspath(__file__))
NMODES = 15
for kind in ("rbe2", "rbe3"):
    m = ms.build(kind)
    for s in ("ccx", "fistr"):
        os.makedirs(os.path.join(D, kind, s), exist_ok=True)
    mm.write_ccx(m, os.path.join(D, kind, "ccx", "machine.inp"), nmodes=NMODES)
    mats = mm.write_fistr_mesh(m, os.path.join(D, kind, "fistr", "machine.msh"))
    mm.write_fistr_eigen_cnt(os.path.join(D, kind, "fistr", "machine_eigen.cnt"), mats, nmodes=NMODES)
    mm.write_fistr_ctrl(os.path.join(D, kind, "fistr", "hecmw_ctrl.dat"), "machine.msh", "machine_eigen.cnt", "machine_eigen.res")
print("written")
