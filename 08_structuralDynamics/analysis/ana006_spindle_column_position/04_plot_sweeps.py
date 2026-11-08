"""Z（主軸頭の高さ）と Y（コラム位置）のスイープ結果をまとめる → results/sweep_Z.png, sweep_Y.png, sweep_trend.png, summary.md"""
import sys, json
from pathlib import Path
import numpy as np
D = Path(__file__).resolve().parent; sys.path[:0] = [str(D), str(D.parent.parent / "study" / "scripts")]; R = D / "results"
import params as P
from figstyle import plt
from lib.chatter import GRID
S = json.loads((R / "summary.json").read_text())
LAB = {"Z": ("主軸頭の高さ", "dz", "dz = {:+d} mm"), "Y": ("コラムの Y 位置", "dy", "dy = {:+d} mm")}
for key, names in P.SWEEPS.items():
    title, par, fmt = LAB[key]; cmap = plt.cm.viridis(np.linspace(0, 0.85, len(names)))
    fig, axs = plt.subplots(3, 1, figsize=(12, 13), gridspec_kw=dict(height_ratios=[1, 1, 0.7]))
    for c, n in zip(cmap, names):
        z = np.load(R / f"case_{n}.npz"); v = P.CASES[n][0 if key == "Y" else 1]; lab = fmt.format(int(v))
        axs[0].loglog(GRID, np.abs(z["Gxx"]) * 1e3, color=c, lw=1.6, label=lab + "（$G_{xx}$）")
        axs[0].loglog(GRID, np.abs(z["Gyy"]) * 1e3, color=c, lw=1.2, ls="--", label=lab + "（$G_{yy}$）")
        axs[1].plot(z["rpm"], z["lim"], color=c, lw=1.8, label=lab); k = np.nanargmin(z["lim"]); axs[1].plot(z["rpm"][k], z["lim"][k], "o", color=c, ms=8)
        axs[2].plot(z["rpm"], z["fc"], color=c, lw=1.3, label=lab)
    axs[0].set_xlim(20, 3000); axs[0].set_ylim(1e-4, 1); axs[0].grid(alpha=.3, which="both"); axs[0].set_ylabel("コンプライアンス [µm/N]"); axs[0].legend(fontsize=9, ncol=2)
    axs[0].set_title(f"① 刃先―ワーク相対コンプライアンス（実線 X-X、破線 Y-Y）：{title}を変えたとき", loc="left")
    axs[1].set_ylim(0, 12); axs[1].grid(alpha=.3); axs[1].set_ylabel("軸方向切込み限界 $a_{\\lim}$ [mm]"); axs[1].legend(fontsize=10)
    axs[1].set_title("② 安定限界線図（4 刃フルスロット、$K_t$ = 1000 N/mm², $K_r$ = 0.3）。丸 = 最小の点", loc="left")
    axs[2].set_yscale("log"); axs[2].grid(alpha=.3, which="both"); axs[2].set_xlabel("主軸回転数 [min⁻¹]"); axs[2].set_ylabel("びびり周波数 [Hz]")
    axs[2].set_title("③ 境界のびびり周波数", loc="left")
    fig.tight_layout(); fig.savefig(R / f"sweep_{key}.png", dpi=90); plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
for key, names, col, mk in (("Z", P.SWEEPS["Z"], "#1f77b4", "o"), ("Y", P.SWEEPS["Y"], "#d62728", "s")):
    x = [P.CASES[n][1 if key == "Z" else 0] for n in names]
    lab = "主軸頭の高さ dz" if key == "Z" else "コラム位置 dy"
    axs[0].plot(x, [S[n]["amin_mm"] for n in names], mk + "-", color=col, label=lab)
    axs[1].plot(x, [S[n]["static_um_per_N"]["xx"] * 1e3 for n in names], mk + "-", color=col, label=lab + "：X")
    axs[1].plot(x, [S[n]["static_um_per_N"]["yy"] * 1e3 for n in names], mk + "--", color=col, mfc="white", label=lab + "：Y")
    axs[2].plot(x, [S[n]["fc_at_min_Hz"] for n in names], mk + "-", color=col, label=lab)
for ax, t, yl in zip(axs, ["① いちばん小さい切込み限界", "② 静剛性の逆数（静コンプライアンス）", "③ そのときのびびり周波数"],
                     ["$a_{\\lim,\\min}$ [mm]", "[nm/N]", "[Hz]"]):
    ax.set_title(t, loc="left"); ax.set_ylabel(yl); ax.set_xlabel("位置の変化量 [mm]"); ax.grid(alpha=.3); ax.legend(fontsize=9)
axs[0].set_ylim(bottom=0); axs[1].set_ylim(bottom=0)
fig.tight_layout(); fig.savefig(R / "sweep_trend.png", dpi=90); plt.close(fig)

md = ["| ケース | dy [mm] | dz [mm] | 最小 $a_{\\lim}$ [mm] | そのときの回転数 [min⁻¹] | びびり周波数 [Hz] | 静コンプライアンス X / Y [nm/N] | 主なモード（割合） |", "|---|---:|---:|---:|---:|---:|---:|---|"]
for n in P.CASES:
    if n not in S:
        continue
    s = S[n]; md.append(f"| {n} | {s['dy']:+d} | {s['dz']:+d} | {s['amin_mm']:.2f} | {s['rpm_at_min']:.0f} | {s['fc_at_min_Hz']:.0f} | {s['static_um_per_N']['xx']*1e3:.1f} / {s['static_um_per_N']['yy']*1e3:.1f} | "
              + "、".join(f"モード {t['mode']}（{t['f_Hz']:.0f} Hz, {t['share']*100:.0f} %）" for t in s["top_modes"][:3]) + " |")
(R / "summary.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
