"""ana005 4.2d：モードが変わっても同じ議論ができるか。全体 FEM のモード（10 個）で加減速の時刻歴を計算し、1 自由度モデルと比べる
   力：テーブルを加速する反力。本体（ベッド）から見ると、テーブルの全質量に慣性力 −m_t a(t) がかかるのと同じ（1 自由度モデルの右辺と同じ）。
       テーブル（＋ワーク）の節点に一様に配る。各モード r は
     q_r'' + 2ζ_r ω_r q_r' + ω_r² q_r = φ_r^T f(t)   （質量正規化モード、ζ_r は Rayleigh）
   観測：ベッド側面 X、コラム上端 X、刃先―ワーク相対 X"""
import sys, json
import numpy as np
from pathlib import Path
from scipy.signal import lsim
D = Path(__file__).resolve().parent; P = D.parent
sys.path[:0] = [str(P), str(P.parent / "tools"), str(P.parent.parent / "study" / "scripts"), str(P.parent / "ana003_ccx_fistr_mpc_spring")]
from frd import read_frd
from figstyle import plt, BLUE, RED, GRAY, GREEN, ORANGE
from frequency_response import ALPHA, BETA
import fn_effect as FE
info = json.loads((D / "info.json").read_text()); meta = json.loads((P / "metadata.json").read_text())["foundation"]
WP, TP = meta["contact_patches"]["WORK"]["nodes"], meta["contact_patches"]["TOOL"]["nodes"]   # 代表点は frd に出ないので面の節点の平均を使う
CASES = [("k025", "1/4", RED, 11.5), ("k1", "今", BLUE, 22.9), ("k4", "4 倍", GREEN, 45.8)]
mt = 125.122e-3                      # t（mm・N・s 系）
from vmc_model import build
_m = build("foundation"); TABLE = [n for n, *_ in _m.parts[1].coords]
out = {}
for tj, tjlab in ((0.0, "ジャーク 0 ms"), (0.04, "ジャーク 40 ms")):
    t, a, tm = FE.profile(tj); a_mm = a * 1000          # mm/s²
    for c, lab, col, fn1 in CASES:
        _, blocks = read_frd(str(D / c / "eig" / "model.frd"))
        nb, nt = info[c]["nut"]; bp, ct = info[c]["nodes"]["BEDPT"], info[c]["nodes"]["COLTOP"]
        resp = {k: np.zeros_like(t) for k in ("bed", "col", "tool", "work", "rel")}
        for name, fq, dd in blocks:
            w = 2 * np.pi * fq; z = ALPHA / (2 * w) + BETA * w / 2
            g = -np.mean([dd.get(n, (0, 0, 0))[0] for n in TABLE])             # φ^T f / (m_t a)（テーブルに一様な −m_t a）
            _, q, _ = lsim(([1.0], [1.0, 2 * z * w, w * w]), g * mt * a_mm, t)
            resp["bed"] += dd.get(bp, (0, 0, 0))[0] * q
            resp["col"] += dd.get(ct, (0, 0, 0))[0] * q
            uw = np.mean([dd.get(n, (0, 0, 0))[0] for n in WP]); ut = np.mean([dd.get(n, (0, 0, 0))[0] for n in TP])
            resp["tool"] += ut * q; resp["work"] += uw * q; resp["rel"] += (uw - ut) * q
        _, x1, _, _, _ = FE.respond(fn1, tj)
        out[(tj, c)] = dict(t=t, tm=tm, a=a, x1=x1 * 1e6, **{k: v * 1000 for k, v in resp.items()}, f=[b[1] for b in blocks])
fig, axs = plt.subplots(5, 2, figsize=(16, 17), sharex=True)
for j, tj in enumerate((0.0, 0.04)):
    for i, (key, ttl) in enumerate((("bed", "ベッド側面の X [µm]"), ("col", "コラム上端の X [µm]"), ("tool", "刃先（テストバー先端）の X [µm]"), ("work", "ワーク（テーブル上）の X [µm]"), ("rel", "刃先―ワークの相対 X [µm]"))):
        ax = axs[i, j]
        for c, lab, col, fn1 in CASES:
            o = out[(tj, c)]
            ax.plot(o["t"], o[key], color=col, lw=1.8, label=f"脚 {lab}（FEM、1 次 {o['f'][0]:.1f} Hz）")
            if key == "bed": ax.plot(o["t"], o["x1"], ":", color=col, lw=1.4, label=f"1 自由度モデル {fn1} Hz")
        ax.axvline(out[(tj, 'k1')]["tm"], color=GRAY, lw=0.8); ax.grid(alpha=0.3); ax.set_ylabel(ttl)
        if i == 0: ax.set_title(f"{'ジャーク 0 ms（段差）' if tj == 0 else 'ジャーク 40 ms'}：同じ 50 mm 移動・2 m/s²", loc="left", fontsize=12)
        if i == 0 and j == 0: ax.legend(fontsize=8.5, ncol=2, loc="lower right")
        if i == 1 and j == 0: ax.legend(fontsize=9, loc="lower right")
    axs[4, j].set_xlabel("時刻 [s]"); axs[4, j].set_xlim(0, 0.8)
fig.suptitle("全体 FEM のモード（10 個）で加減速の揺れを計算：ベッド・コラム・刃先・ワークの床に対する揺れと、刃先―ワークの相対", fontsize=13)
fig.tight_layout(); fig.savefig(D / "transient_modal.png")
rows = ["| ケース | ジャーク | 刃先の最大 µm | ワークの最大 µm | 刃先―ワーク相対の最大 µm | 停止後のベッド揺れ最大 µm | 停止後のコラム上端揺れ最大 µm | 1 自由度モデル µm | ベッドが ±1 µm に収まるまで s | コラム上端が ±1 µm に収まるまで s |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
for tj in (0.0, 0.04):
    for c, lab, col, fn1 in CASES:
        o = out[(tj, c)]; m = o["t"] >= o["tm"]
        def settle(x):
            late = o["t"][m][np.abs(x[m]) > 1]; return 0.0 if len(late) == 0 else late[-1] - o["tm"]
        rows.append(f"| 脚 {lab} | {tj*1000:.0f} ms | {np.abs(o['tool']).max():.2f} | {np.abs(o['work']).max():.2f} | {np.abs(o['rel']).max():.3f} | {np.abs(o['bed'][m]).max():.2f} | {np.abs(o['col'][m]).max():.2f} | {np.abs(o['x1'][m]).max():.2f} | {settle(o['bed']):.3f} | {settle(o['col']):.3f} |")
(D / "transient_modal.md").write_text("\n".join(rows) + "\n"); print("\n".join(rows))
