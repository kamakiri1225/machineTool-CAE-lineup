"""ana002: 簡易工作機械モデルの周波数応答解析（モード重ね合わせ法）
テーブル上面中央（MON）を X/Y/Z に 1 N で加振し、同じ点の変位（コンプライアンス）を求める。"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import machine_model as mm

D = os.path.dirname(os.path.abspath(__file__))
NMODES = 30
F0, F1, NPT = 50.0, 1000.0, 400
# Rayleigh 減衰: zeta(w) = alpha/(2w) + beta*w/2 を 150 Hz と 600 Hz で 2 % にする
ZETA, FA, FB = 0.02, 150.0, 600.0
wa, wb = 2 * np.pi * FA, 2 * np.pi * FB
ALPHA = 2 * ZETA * wa * wb / (wa + wb)
BETA = 2 * ZETA / (wa + wb)

if __name__ == "__main__":
    m = mm.build()
    # CalculiX: 固有値ステップ + 方向ごとの定常応答ステップ（1 ファイル）
    mm.write_ccx(m, os.path.join(D, "ccx", "machine.inp"), nmodes=NMODES)
    txt = open(os.path.join(D, "ccx", "machine.inp")).read()
    for dof in (1, 2, 3):
        txt += "\n".join(["*STEP", "*STEADY STATE DYNAMICS", f"{F0}, {F1}, 20, 1.0",
                          "*MODAL DAMPING, RAYLEIGH", f",,{ALPHA:.6e},{BETA:.6e}",
                          "*CLOAD", f"MON, {dof}, 1.0", "*NODE PRINT, NSET=MON", "U", "*END STEP"]) + "\n"
    open(os.path.join(D, "ccx", "machine.inp"), "w").write(txt)
    # FrontISTR: 固有値解析 → 方向ごとに周波数応答
    mats = mm.write_fistr_mesh(m, os.path.join(D, "fistr", "machine.msh"))
    mm.write_fistr_eigen_cnt(os.path.join(D, "fistr", "eigen.cnt"), mats, nmodes=NMODES)
    mm.write_fistr_ctrl(os.path.join(D, "fistr", "hecmw_ctrl_eigen.dat"), "machine.msh", "eigen.cnt", "eigen.res")
    for dof in (1, 2, 3):
        mm.write_fistr_freq_cnt(os.path.join(D, "fistr", f"freq_{dof}.cnt"), mats, F0, F1, NPT, ALPHA, BETA,
                                m.monitor, NMODES, load_dof=dof)
        mm.write_fistr_ctrl(os.path.join(D, "fistr", f"hecmw_ctrl_freq_{dof}.dat"), "machine.msh", f"freq_{dof}.cnt",
                            f"freq_{dof}.res", eig_res="eigen.res", dyna=f"freq_{dof}_dyna.res")
    print(f"alpha={ALPHA:.4g}, beta={BETA:.4g}, monitor node={m.monitor}")
