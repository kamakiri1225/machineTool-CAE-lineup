"""FrontISTR 5.9 と 6.0 を、ana001〜004 と同じ入力で動かして比べる。
   入力は各 ana の fistr フォルダからコピーし、run/<版>/<ケース>/ で実行する（元の結果は変えない）。
   実行: ~/opt/venv/bin/python compare.py"""
import os, re, shutil, subprocess, time, glob
import numpy as np

D = os.path.dirname(os.path.abspath(__file__)); A = os.path.dirname(D)
BIN = {"5.9": os.path.expanduser("~/opt/frontistr59/bin/fistr1"), "6.0": os.path.expanduser("~/opt/frontistr6/bin/fistr1")}
R = os.path.expanduser("~/opt/ccxroot/usr/lib/x86_64-linux-gnu")
ENV = dict(os.environ, LD_LIBRARY_PATH=f"{R}/blas:{R}/lapack", OMP_NUM_THREADS="8")
CASES = {  # ケース名: (元フォルダ, コピーするファイルのパターン, 実行する全体制御ファイルの順番)
    "ana001_eigen": ("ana001_ccx_fistr_modal/fistr", ["machine.msh", "machine_eigen.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
    "ana002_frf_X": ("ana002_ccx_fistr_frf/fistr", ["machine.msh", "eigen.cnt", "freq_1.cnt", "hecmw_ctrl_eigen.dat", "hecmw_ctrl_freq_1.dat"],
                     ["hecmw_ctrl_eigen.dat", "hecmw_ctrl_freq_1.dat"]),
    "ana003_rbe2": ("ana003_ccx_fistr_mpc_spring/rbe2/fistr", ["machine.msh", "machine_eigen.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
    "ana003_rbe3": ("ana003_ccx_fistr_mpc_spring/rbe3/fistr", ["machine.msh", "machine_eigen.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
    "ana004_A": ("ana004_ccx_fistr_spring_types/A_ground/fistr", ["model.msh", "model.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
    "ana004_B": ("ana004_ccx_fistr_spring_types/B_spring2/fistr", ["model.msh", "model.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
    "ana004_C": ("ana004_ccx_fistr_spring_types/C_springa/fistr", ["model.msh", "model.cnt", "hecmw_ctrl.dat"], ["hecmw_ctrl.dat"]),
}


def eig_freqs(log):
    f = []
    for ln in open(log):
        p = ln.split()
        if len(p) >= 4 and p[0].isdigit() and "E+" in p[1] and len(f) < 30:
            try:
                f.append(float(p[3]))
            except ValueError:
                pass
    return np.array(f)


def frf(log):
    out = [(float(m.group(1)), float(m.group(2))) for m in
           (re.match(r"\s*([\d.E+-]+)\s+\[Hz\]\s*:\s*([\d.E+-]+)", ln) for ln in open(log)) if m]
    return np.array(out)


res = {}
for ver, b in BIN.items():
    for name, (src, files, ctrls) in CASES.items():
        d = os.path.join(D, "run", ver.replace(".", ""), name)
        shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
        for f in files:
            shutil.copy(os.path.join(A, src, f), d)
        t0 = time.time(); ok = True
        for k, c in enumerate(ctrls):
            shutil.copy(os.path.join(d, c), os.path.join(d, "hecmw_ctrl.dat")) if c != "hecmw_ctrl.dat" else None
            r = subprocess.run([b], cwd=d, env=ENV, stdout=open(os.path.join(d, f"out_{k}.txt"), "w"), stderr=subprocess.STDOUT)
            ok &= r.returncode == 0
            if name == "ana002_frf_X" and k == 0:
                shutil.move(os.path.join(d, "0.log"), os.path.join(d, "eigen_0.log"))
        res[(ver, name)] = dict(t=time.time() - t0, ok=ok, d=d)
        print(ver, name, "ok" if ok else "FAIL", f"{time.time() - t0:.1f}s", flush=True)

L = ["# FrontISTR 5.9 と 6.0 の比較（同じ入力）", "", "| ケース | 5.9 時間 [s] | 6.0 時間 [s] | 終了 5.9 / 6.0 | 比較 |", "|---|---|---|---|---|"]
for name in CASES:
    a, b = res[("5.9", name)], res[("6.0", name)]
    cmp = ""
    try:
        if name.startswith(("ana001", "ana003")):
            fa, fb = eig_freqs(os.path.join(a["d"], "0.log")), eig_freqs(os.path.join(b["d"], "0.log"))
            n = min(len(fa), len(fb), 10)
            cmp = "1〜6次 [Hz] 5.9: " + ", ".join(f"{x:.1f}" for x in fa[:6]) + " / 6.0: " + ", ".join(f"{x:.1f}" for x in fb[:6]) + \
                  f"（10次までの最大差 {np.max(np.abs(fa[:n] - fb[:n]) / fa[:n]) * 100:.3f} %）"
        elif name.startswith("ana002"):
            ha, hb = frf(os.path.join(a["d"], "0.log")), frf(os.path.join(b["d"], "0.log"))
            ia, ib = np.argmax(ha[:, 1]), np.argmax(hb[:, 1])
            cmp = (f"ピーク 5.9: {ha[ia,0]:.1f} Hz, {ha[ia,1]*1e3:.4f} µm/N / 6.0: {hb[ib,0]:.1f} Hz, {hb[ib,1]*1e3:.4f} µm/N、"
                   f"50 Hz の値 5.9: {ha[0,1]*1e3:.5f} / 6.0: {hb[0,1]*1e3:.5f} µm/N")
            np.save(os.path.join(D, "frf_59.npy"), ha); np.save(os.path.join(D, "frf_60.npy"), hb)
        else:
            import pyvista as pv
            mon = {"ana004_A": (5, 2), "ana004_B": (11, 0), "ana004_C": (11, 2)}[name]
            v = []
            for x in (a, b):
                m = pv.read(sorted(glob.glob(os.path.join(x["d"], "*_vis_psf.*.pvtu")))[-1])
                if "NODE_ID" in m.point_data:
                    ids = list(m.point_data["NODE_ID"])
                else:   # メッシュファイルの !NODE の順番
                    txt = open(os.path.join(x["d"], "model.msh")).read().split("!NODE")[1].split("!ELEMENT")[0]
                    ids = [int(l.split(",")[0]) for l in txt.strip().splitlines()[0:] if l.strip() and l[0].isdigit()]
                k = ids.index(mon[0])
                v.append(m.point_data["DISPLACEMENT"][k][mon[1]])
            cmp = f"変位 5.9: {v[0]:.5f} mm / 6.0: {v[1]:.5f} mm（理論 ±0.25）"
    except Exception as e:
        cmp = f"比較できず: {e}"
    L.append(f"| {name} | {a['t']:.1f} | {b['t']:.1f} | {'OK' if a['ok'] else 'NG'} / {'OK' if b['ok'] else 'NG'} | {cmp} |")
open(os.path.join(D, "compare.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
