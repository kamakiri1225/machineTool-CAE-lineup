"""ana001: 簡易工作機械モデルの固有値解析（CalculiX と FrontISTR で同じモデル）"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm

D = os.path.dirname(os.path.abspath(__file__))
NMODES = 15
m = mm.build()
mm.write_ccx(m, os.path.join(D, "ccx", "machine.inp"), nmodes=NMODES)
mats = mm.write_fistr_mesh(m, os.path.join(D, "fistr", "machine.msh"))
mm.write_fistr_eigen_cnt(os.path.join(D, "fistr", "machine_eigen.cnt"), mats, nmodes=NMODES)
mm.write_fistr_ctrl(os.path.join(D, "fistr", "hecmw_ctrl.dat"), "machine.msh", "machine_eigen.cnt", "machine_eigen.res")
print("written")
