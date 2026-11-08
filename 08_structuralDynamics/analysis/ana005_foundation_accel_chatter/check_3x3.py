"""5.2 補足：ねじれ刃で F_z が出る場合、3×3（Z 列を含む）にすると限界がどれだけ変わるか
   F = K_t a h [(e_t + K_r e_r); K_a/K_t]、h は X・Y の揺れだけで決まる → B0 は 3×2、Z 行の平均は (Z/π)(K_a/K_t)(1, 0)
   Φ = G[xy, xyz] · B0(3×2) を 2×2 の特性方程式に入れる"""
import json, sys
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
from analyze import boundary
from figstyle import plt
TEETH, KT, KR = 4, 1000.0, 0.3
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2; d = phi[iw] - phi[it]          # (3, nm)
grid = np.linspace(5, 3000, 24000); rpm = np.linspace(500, 16000, 600)
w = 2 * np.pi * grid; den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
G = np.einsum("ir,jr,wr->wij", d, d, 1 / den)                                                      # (nw, 3, 3)
print("max |G_xz|/|G_xx| =", float(np.max(np.abs(G[:, 0, 2])) / np.max(np.abs(G[:, 0, 0]))), " max |G_zz|/|G_xx| =", float(np.max(np.abs(G[:, 2, 2])) / np.max(np.abs(G[:, 0, 0]))))


def lobes(ka):
    B = np.zeros((3, 2)); B[:2] = TEETH / 4 * np.array([[KR, 1.0], [-1.0, KR]]); B[2] = TEETH / np.pi * ka * np.array([1.0, 0.0])
    lam = np.linalg.eigvals(G[:, :2, :] @ B)
    for i in range(1, len(lam)):
        if np.sum(abs(lam[i] - lam[i - 1])) > np.sum(abs(lam[i, ::-1] - lam[i - 1])): lam[i] = lam[i, ::-1]
    br = [boundary(interp1d(grid, lam[:, j], bounds_error=True), rpm, grid, Kc=KT, teeth=TEETH) for j in (0, 1)]
    v = np.stack([b[0] for b in br], 1); return np.nanmin(np.where(np.isfinite(v), v, np.nan), 1)


out = {}
fig, ax = plt.subplots(figsize=(12, 4.8))
for ka, c, lab in [(0.0, "k", "2×2（$F_z=0$、今回）"), (0.3, "#c0392b", "3×3、$K_a/K_t=+0.3$"), (-0.3, "#1f618d", "3×3、$K_a/K_t=-0.3$")]:
    lim = lobes(ka); out[str(ka)] = dict(min=float(np.nanmin(lim)), rpm=float(rpm[np.nanargmin(lim)]), maxdiff_pct=None); ax.plot(rpm, lim, color=c, lw=1.8, label=lab); out[str(ka)]["lim"] = lim
for k in ("0.3", "-0.3"):
    out[k]["maxdiff_pct"] = float(np.nanmax(np.abs(out[k]["lim"] / out["0.0"]["lim"] - 1)) * 100)
for k in out: print(k, round(out[k]["min"], 2), "mm @", round(out[k]["rpm"]), "rpm, max diff %", out[k]["maxdiff_pct"]); out[k].pop("lim")
ax.set_ylim(0, 40); ax.grid(alpha=0.3); ax.legend(fontsize=11); ax.set_xlabel("主軸回転数 [min⁻¹]"); ax.set_ylabel("軸方向切込み限界 [mm]")
ax.set_title("ねじれ刃の軸方向力 $F_z$ を入れた 3×3 との比較（基礎ばねモデル・テストバー、相対 FRF）", loc="left", fontsize=13)
fig.tight_layout(); fig.savefig(R / "check_3x3.png", dpi=100); json.dump(out, open(R / "check_3x3.json", "w"), indent=1)
