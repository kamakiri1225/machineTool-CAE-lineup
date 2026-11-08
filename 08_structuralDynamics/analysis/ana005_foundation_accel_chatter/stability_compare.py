"""ana005 5 章：同じ基礎ばねモデルで、どの FRF を使うかによって安定限界線図とびびりモードがどう変わるか
   A：刃先だけを叩いた FRF（ワークは剛と仮定）   G = H_tt
   B：ワーク（テーブル）だけを叩いた FRF（工具は剛と仮定） G = H_ww
   C：刃先―ワークの相対 FRF（正しい）             G = H_ww + H_tt − H_wt − H_tw
   モード基底（results/foundation_modal_basis.npz、80 モード、質量正規化）から作る。切削条件・計算法は milling.py と同じ"""
import json, sys
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
from analyze import boundary
from figstyle import plt, BLUE, RED, GRAY, GREEN, ORANGE
TEETH, KT, KR = 4, 1000.0, 0.3
B0 = TEETH / 4 * np.array([[KR, 1.0], [-1.0, KR]])
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2
pw, pt = phi[iw, :2, :], phi[it, :2, :]                           # (2, nm)
DEF = {"A": ("A：刃先の FRF だけ", pt, pt, BLUE), "B": ("B：テーブル（ワーク）の FRF だけ", pw, pw, ORANGE), "C": ("C：刃先―ワーク相対 FRF", pw - pt, pw - pt, RED)}
grid = np.linspace(5, 3000, 24000); rpm = np.linspace(500, 16000, 600)


def G_of(f, d):
    w = 2 * np.pi * np.atleast_1d(f); den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
    return np.einsum("ir,jr,wr->wij", d, d, 1 / den)


def contrib(f, d):
    w = 2 * np.pi * f; den = wr**2 - w**2 + 2j * zr * wr * w
    c = np.einsum("ir,jr,r->rij", d, d, 1 / den); n = np.linalg.norm(c.reshape(len(fr), -1), axis=1); return n / n.sum()


res = {}
for key, (lab, d, _, col) in DEF.items():
    Gg = G_of(grid, d); lam = np.linalg.eigvals(Gg @ B0)
    for i in range(1, len(lam)):
        if np.sum(abs(lam[i] - lam[i - 1])) > np.sum(abs(lam[i, ::-1] - lam[i - 1])): lam[i] = lam[i, ::-1]
    br = [boundary(interp1d(grid, lam[:, j], bounds_error=True), rpm, grid, Kc=KT, teeth=TEETH) for j in (0, 1)]
    v = np.stack([b[0] for b in br], 1); fq = np.stack([b[1] for b in br], 1); w_ = np.argmin(np.where(np.isfinite(v), v, np.inf), 1)
    lim = v[np.arange(len(rpm)), w_]; fc = fq[np.arange(len(rpm)), w_]
    k = np.nanargmin(lim); sh = contrib(fc[k], d); top = np.argsort(sh)[::-1][:3]
    res[key] = dict(lab=lab, col=col, lim=lim, fc=fc, kmin=k, top=[(int(r), float(fr[r]), float(sh[r])) for r in top], Gg=Gg)
    print(key, f"min {lim[k]:.2f} mm @ {rpm[k]:.0f} rpm, chatter {fc[k]:.0f} Hz, modes", [(r + 1, round(f_, 1), round(s_ * 100)) for r, f_, s_ in res[key]["top"]])
fig, axs = plt.subplots(3, 1, figsize=(13, 13), gridspec_kw=dict(height_ratios=[0.9, 1, 0.8]))
ax = axs[0]
for key, r in res.items():
    ax.loglog(grid, np.abs(r["Gg"][:, 0, 0]) * 1e3, color=r["col"], lw=1.8, label=r["lab"] + "（X-X）")
ax.set_xlim(5, 3000); ax.set_ylim(1e-4, 1); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=10); ax.set_ylabel("コンプライアンス [µm/N]")
ax.set_title("① 使う FRF（X-X）：A・B は床に対する動きなので 22 Hz（脚の上で機械全体が揺れる）のピークを持つ。C ではそれが消える", loc="left", fontsize=12)
ax = axs[1]
for key, r in res.items():
    ax.plot(rpm, r["lim"], color=r["col"], lw=2, label=r["lab"]); k = r["kmin"]; ax.plot(rpm[k], r["lim"][k], "o", color=r["col"], ms=9)
ax.set_ylim(0, 16); ax.grid(alpha=0.3); ax.legend(fontsize=10); ax.set_ylabel("軸方向切込み限界 $a_p$ [mm]")
ax.set_title("② 安定限界線図（4 刃フルスロット、$K_t$ = 1000 N/mm²、$K_r$ = 0.3）。丸 = 最小の点。B は 16 mm 以上を省略（最大 136 mm）", loc="left", fontsize=12)
ax = axs[2]
for key, r in res.items(): ax.plot(rpm, r["fc"], color=r["col"], lw=1.5, label=r["lab"])
ax.grid(alpha=0.3); ax.set_xlabel("主軸回転数 [min⁻¹]"); ax.set_ylabel("境界のびびり周波数 [Hz]"); ax.set_yscale("log")
ax.set_title("③ そのときのびびり周波数", loc="left", fontsize=12)
fig.tight_layout(); fig.savefig(R / "stability_compare.png")
md = ["| 使う FRF | 最小の切込み限界 | その回転数 | びびり周波数 | 効いているモード（寄与） |", "|---|---:|---:|---:|---|"]
for key, r in res.items():
    k = r["kmin"]; md.append(f"| {r['lab']} | {r['lim'][k]:.2f} mm | {rpm[k]:.0f} min⁻¹ | {r['fc'][k]:.0f} Hz | " + "、".join(f"モード {m+1}（{f_:.0f} Hz, {s_*100:.0f} %）" for m, f_, s_ in r["top"]) + " |")
np.savez(R / "stability_compare.npz", rpm=rpm, **{f"lim_{k}": r["lim"] for k, r in res.items()}, **{f"fc_{k}": r["fc"] for k, r in res.items()})
(R / "stability_compare.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
