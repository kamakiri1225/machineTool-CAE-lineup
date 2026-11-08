"""主軸頭を上げた効果を「刃先側（主軸頭・コラム）」と「ワーク側（治具ブロック）」に分ける
   ・静コンプライアンス：刃先側 Σ t t^T/ω²、ワーク側 Σ w w^T/ω²
   ・刃先側だけの FRF（ワークは剛と仮定）で引いた安定限界 → 主軸の高さそのものの効果
   → results/z_split.png, results/z_split.json"""
import sys, json
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; sys.path[:0] = [str(D), str(D.parent.parent / "study" / "scripts")]; R = D / "results"
import params as P
from figstyle import plt
from lib.chatter import lobes, RPM
names = P.SWEEPS["Z"]; dz = [P.CASES[n][1] for n in names]; out = {}
fig, axs = plt.subplots(1, 3, figsize=(17, 5))
cols = plt.cm.viridis(np.linspace(0, 0.85, len(names)))
for n, c in zip(names, cols):
    z = np.load(R / f"case_{n}.npz"); w2 = (2 * np.pi * z["fr"]) ** 2
    st = lambda v: (v[:2] ** 2 / w2).sum(1) * 1e6        # nm/N（xx, yy）
    lim_t, fc_t, _ = lobes(z["fr"], -z["tool"])
    hs = RPM > 3000
    out[n] = dict(dz=P.CASES[n][1], static_tool_nm_N=st(z["tool"]).tolist(), static_work_nm_N=st(z["work"]).tolist(), static_rel_nm_N=st(z["d"]).tolist(),
                  amin_rel=float(np.nanmin(z["lim"])), amin_rel_high=float(np.nanmin(np.where(hs, z["lim"], np.nan))),
                  amin_tool=float(np.nanmin(lim_t)), amin_tool_high=float(np.nanmin(np.where(hs, lim_t, np.nan))))
    axs[1].plot(RPM, lim_t, color=c, lw=1.8, label=f"dz = {P.CASES[n][1]:+d} mm")
k = lambda key, i: [out[n][key][i] for n in names]
axs[0].plot(dz, k("static_tool_nm_N", 0), "o-", color="#1f77b4", label="刃先側 X"); axs[0].plot(dz, k("static_tool_nm_N", 1), "o--", color="#1f77b4", mfc="white", label="刃先側 Y")
axs[0].plot(dz, k("static_work_nm_N", 0), "s-", color="#ff7f0e", label="ワーク側 X"); axs[0].plot(dz, k("static_work_nm_N", 1), "s--", color="#ff7f0e", mfc="white", label="ワーク側 Y")
axs[0].set_title("① 床に対する静コンプライアンス（刃先側・ワーク側）", loc="left"); axs[0].set_xlabel("主軸頭の高さ dz [mm]"); axs[0].set_ylabel("[nm/N]"); axs[0].legend(fontsize=9); axs[0].grid(alpha=.3); axs[0].set_ylim(bottom=0)
axs[1].set_ylim(0, 12); axs[1].grid(alpha=.3); axs[1].legend(fontsize=9); axs[1].set_xlabel("主軸回転数 [min⁻¹]"); axs[1].set_ylabel("$a_{\\lim}$ [mm]")
axs[1].set_title("② 刃先側だけの FRF で引いた線図（ワーク・治具は剛）", loc="left")
axs[2].plot(dz, [out[n]["amin_rel"] for n in names], "o-", color="#d62728", label="相対 FRF（治具込み）：全回転数での最小")
axs[2].plot(dz, [out[n]["amin_rel_high"] for n in names], "o--", color="#d62728", mfc="white", label="相対 FRF：3,000 min⁻¹ 以上での最小")
axs[2].plot(dz, [out[n]["amin_tool"] for n in names], "s-", color="#1f77b4", label="刃先側だけ：全回転数での最小")
axs[2].set_title("③ いちばん小さい切込み限界", loc="left"); axs[2].set_xlabel("主軸頭の高さ dz [mm]"); axs[2].set_ylabel("$a_{\\lim,\\min}$ [mm]"); axs[2].set_ylim(0, 4); axs[2].grid(alpha=.3); axs[2].legend(fontsize=9)
fig.tight_layout(); fig.savefig(R / "z_split.png", dpi=90)
(R / "z_split.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))
