"""ana004 追加: ケース B（上下の立方体を同じ位置の節点 4 組でつなぐ）で、
   「ばね（SPRING2）」と「MPC（*EQUATION で u_上 = u_下）」を比べる。ばね定数を変えて解き、MPC は「無限に硬いばね」になることを示す。
   実行: ~/opt/venv/bin/python spring_vs_mpc.py"""
import os, subprocess, shutil
import numpy as np
ns = {"__file__": os.path.join(os.path.dirname(os.path.abspath(__file__)), "make_and_run.py")}
exec(open(ns["__file__"]).read().split("env = dict(os.environ)")[0], ns)
CASES, write_ccx, write_fistr, read_ccx_dat = ns["CASES"], ns["write_ccx"], ns["write_fistr"], ns["read_ccx_dat"]
D = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(D, "spring_vs_mpc"); os.makedirs(W, exist_ok=True)
c = CASES["B_spring2"]
CCX = os.path.expanduser("~/opt/ccx.sh")


def run_ccx(path_dir):
    subprocess.run([CCX, "-i", "model"], cwd=path_dir, stdout=open(os.path.join(path_dir, "ccx.log"), "w"), stderr=subprocess.STDOUT)
    return read_ccx_dat(os.path.join(path_dir, "model.dat"))[c["mon"]][0]


ks = np.logspace(1, 8, 15); us = []
for k in ks:
    d = os.path.join(W, f"k_{k:.0e}"); os.makedirs(d, exist_ok=True)
    ns["K"] = k; write_ccx(c, os.path.join(d, "model.inp")); us.append(run_ccx(d))
ns["K"] = 1000.0
# MPC: ばね要素の代わりに *EQUATION（上の節点の変位 − 下の節点の変位 = 0）
d = os.path.join(W, "mpc"); os.makedirs(d, exist_ok=True)
write_ccx(c, os.path.join(d, "model.inp"))
txt = open(os.path.join(d, "model.inp")).read()
head, rest = txt.split("*ELEMENT, TYPE=SPRING2", 1)
rest = "*NSET" + rest.split("*NSET", 1)[1]
eq = []
for a, b in c["pairs"]:
    for dof in (1, 2, 3):
        eq += ["*EQUATION", "2", f"{b}, {dof}, 1.0, {a}, {dof}, -1.0"]
open(os.path.join(d, "model.inp"), "w").write(head + "\n".join(eq) + "\n" + rest)
u_mpc = run_ccx(d)
# FrontISTR でも MPC（!EQUATION）
df = os.path.join(W, "mpc_fistr"); os.makedirs(df, exist_ok=True)
write_fistr(c, df)
msh = open(os.path.join(df, "model.msh")).read()
a0 = msh.index("!ELEMENT, TYPE=511"); a1 = msh.index("!MATERIAL, NAME=STEEL")
msh = msh[:a0] + msh[a1:]
msh = msh.replace("!MATERIAL, NAME=SPMAT, ITEM=1\n!ITEM=1, SUBITEM=2\n1.0, 0.3\n!SECTION, TYPE=INTERFACE, EGRP=ESPR, MATERIAL=SPMAT\n1.0\n", "")
eqf = []
for a, b in c["pairs"]:
    for dof in (1, 2, 3):
        eqf += ["!EQUATION", "2, 0.0", f"{b}, {dof}, 1.0, {a}, {dof}, -1.0"]
msh = msh.replace("!END", "\n".join(eqf) + "\n!END")
open(os.path.join(df, "model.msh"), "w").write(msh)
cnt = open(os.path.join(df, "model.cnt")).read().replace("!MATERIAL, NAME=SPMAT\n!SPRING_D\n 1, 1, 1000.0\n 2, 2, 1000.0\n 3, 3, 1000.0\n", "")
open(os.path.join(df, "model.cnt"), "w").write(cnt)
R = os.path.expanduser("~/opt/ccxroot/usr/lib/x86_64-linux-gnu")
subprocess.run([os.path.expanduser("~/opt/frontistr/bin/fistr1")], cwd=df, env=dict(os.environ, LD_LIBRARY_PATH=f"{R}/blas:{R}/lapack"),
               stdout=open(os.path.join(df, "fistr.out"), "w"), stderr=subprocess.STDOUT)
import pyvista as pv, glob
m = pv.read(sorted(glob.glob(os.path.join(df, "model_vis_psf.*.pvtu")))[-1])
ids = [n for n, *_ in c["nodes"]]
u_mpc_f = m.point_data["DISPLACEMENT"][ids.index(c["mon"])][0]
print("MPC ccx", u_mpc, "fistr", u_mpc_f)

import sys
sys.path.insert(0, os.path.join(D, "..", "..", "study", "scripts"))
from figstyle import plt, BLUE, RED, GRAY, GREEN, ORANGE
fig, ax = plt.subplots(figsize=(11, 5.6))
ax.loglog(ks, us, "o-", color=BLUE, lw=2, label="ばね（SPRING2）でつないだとき：CalculiX")
ax.loglog(ks, 1000 / (4 * ks), "--", color=GRAY, label="ばねだけの手計算 $u = F/(4k)$")
msg = "MPC（*EQUATION：上の節点の変位 = 下の節点の変位）でつなぐと\n変位は %.0f mm（CalculiX）・%.0e mm（FrontISTR）\n＝ 固定された下の立方体と一体になり、まったく動かない" % (abs(u_mpc), abs(u_mpc_f))
ax.annotate(msg, xy=(1e8, 3e-6), xytext=(12, 4e-5), color=RED, fontsize=11, arrowprops=dict(arrowstyle="->", color=RED))
ax.axvline(1000, color=GREEN, ls=":"); ax.text(1100, 3e-1, "ana004 の値\n$k$ = 1 N/µm", color=GREEN, fontsize=11)
ax.axvspan(1e5, 1e6, color=ORANGE, alpha=0.1); ax.text(1.1e5, 3e-1, "ガイドの\n実際の値の範囲\n（数百〜数千 N/µm）", color=ORANGE, fontsize=10)
ax.set_xlabel("ばね定数 $k$ [N/mm]（1 本あたり）"); ax.set_ylabel("上の立方体の X 変位 [mm]（1000 N）")
ax.set_title("ばねを硬くしていくと MPC の値に近づく ＝ MPC は「無限に硬いばね」", loc="left", fontsize=12)
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=10, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(D, "results", "spring_vs_mpc.png")); plt.close(fig)
open(os.path.join(D, "results", "spring_vs_mpc.md"), "w").write(
    "| つなぎ方 | X 変位 [mm] |\n|---|---|\n" + "".join(f"| ばね k = {k:.0e} N/mm | {u:.6f} |\n" for k, u in zip(ks, us)) +
    f"| MPC（CalculiX *EQUATION） | {u_mpc:.6f} |\n| MPC（FrontISTR !EQUATION） | {u_mpc_f:.6f} |\n")
print("ok")
