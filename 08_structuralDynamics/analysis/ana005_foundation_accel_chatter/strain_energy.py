"""ana005 5.3 補足：びびりを決めるモードで、ひずみエネルギー（変形のエネルギー）がどこにたまっているか
   質量正規化モードでは全ひずみエネルギー = ω_r²/2。
   ばね：各結合面（キャリッジ 4・ナット・脚 4）の面の平均変位の差 Δ から ½ k Δ²（並進）。回転分は面の平均回転の差から ½ kθ Δθ²
   ソリッド：残り（全体 − ばね）。どの部品が曲がっているかは、要素ごとに求める代わりに部品ごとの変形量の目安で示す"""
import sys, json
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent / "tools"), str(D.parent / "ana003_ccx_fistr_mpc_spring"), str(D.parent.parent / "study" / "scripts")]
import machine_model as mm, mpc_springs as ms
from vmc_model import build
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]
m = build("foundation"); coord = {n: np.array((x, y, zz)) for p in m.parts for n, x, y, zz in p.coords}
ifc = []
for name, rect, kd in mm.INTERFACES:
    fb = [n for n, x, y in ms._patch_nodes(m.parts[0], rect, mm.Z_IF)]; ft = [n for n, x, y in ms._patch_nodes(m.parts[1], rect, mm.Z_IF)]
    kt = {d: k for nm, d, k, a, b in m.springs if nm == name}; kr = {d: k for nm, d, k, a, b in m.springs if nm == name + "_ROT"}
    ifc.append((("ガイド " + name[-1]) if name.startswith("CAR") else "ボールねじ（ナット）", fb, ft, kt, kr))
for f in m.foundations:
    kt = {d: k for nm, d, k, a, b in m.springs if nm == f["name"]}; kr = {d: k for nm, d, k, a, b in m.springs if nm == f["name"] + "_ROT"}
    ifc.append(("脚 " + f["name"][-1], f["surface_nodes"], None, kt, kr))


def rigid(face, U):
    """面の節点の変位から、平均の並進と回転（最小二乗）を求める"""
    X = np.array([coord[n] for n in face]); c = X.mean(0); A = []
    for rx, ry, rz in X - c: A += [[1, 0, 0, 0, rz, -ry], [0, 1, 0, -rz, 0, rx], [0, 0, 1, ry, -rx, 0]]
    sol = np.linalg.lstsq(np.array(A), U[[n - 1 for n in face]].reshape(-1), rcond=None)[0]; return sol[:3], sol[3:]


def energies(r):
    U = phi[:, :, r]; total = 0.5 * (2 * np.pi * fr[r]) ** 2; out = {}
    for name, fb, ft, kt, kr in ifc:
        tb, rb = rigid(fb, U)
        if ft is None: tt, rt = np.zeros(3), np.zeros(3)
        else: tt, rt = rigid(ft, U)
        e = sum(0.5 * k * (tt[d - 1] - tb[d - 1]) ** 2 for d, k in kt.items()) + sum(0.5 * k * (rt[d - 1] - rb[d - 1]) ** 2 for d, k in kr.items())
        out[name.split(" ")[0] if name.startswith("ガイド") or name.startswith("脚") else name] = out.get(name.split(" ")[0] if name.startswith(("ガイド", "脚")) else name, 0) + e
    out["鋳物・主軸などのソリッド（曲げ・ねじれ）"] = total - sum(out.values())
    return {k: v / total for k, v in out.items()}


res = {}
for r, lab in ((7, "モード 8（182 Hz）"), (10, "モード 11（431 Hz）"), (17, "モード 18（910 Hz）"), (9, "モード 10（346 Hz）")):
    res[lab] = energies(r); print(lab, {k: round(v * 100, 1) for k, v in res[lab].items()})
keys = list(next(iter(res.values())).keys()); cols = [ORANGE, RED, GRAY, BLUE]
fig, ax = plt.subplots(figsize=(13, 4.8)); left = np.zeros(len(res))
for k, c in zip(keys, [RED, ORANGE, "#8e44ad", BLUE]):
    v = np.array([res[l][k] for l in res]); ax.barh(range(len(res)), v, left=left, color=c, label=k)
    for i, (vv, ll) in enumerate(zip(v, left)):
        if vv > 0.06: ax.text(ll + vv / 2, i, f"{vv*100:.0f} %", ha="center", va="center", color="white", fontsize=12)
    left += v
ax.set_yticks(range(len(res))); ax.set_yticklabels(list(res), fontsize=12); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_xlabel("そのモードの変形エネルギーに占める割合"); ax.legend(fontsize=11, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4)
fig.tight_layout(); fig.savefig(R / "strain_energy.png"); json.dump(res, open(R / "strain_energy.json", "w"), ensure_ascii=False, indent=1)
