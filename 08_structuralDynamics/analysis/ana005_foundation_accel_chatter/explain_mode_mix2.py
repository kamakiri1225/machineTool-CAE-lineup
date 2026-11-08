"""5.3 補足：線図の 1 点は「1 つの周波数」を選んだだけで「1 つのモード」を選んだのではない。
   その周波数での揺れは、全モードの応答の足し合わせになる"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D.parent.parent / "study" / "scripts")]
from figstyle import plt
TEETH, KT, KR = 4, 1000.0, 0.3; B0 = TEETH / 4 * np.array([[KR, 1.0], [-1.0, KR]])
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2; d = phi[iw, :2] - phi[it, :2]
fc = 185.0; wc = 2 * np.pi * fc; den_c = wr**2 - wc**2 + 2j * zr * wr * wc
G = np.einsum("ir,jr,r->ij", d, d, 1 / den_c); lam, V = np.linalg.eig(G @ B0); j = np.argmax((-1 / lam).real); v = V[:, j]
F = B0 @ v; eta = (d.T @ F) / den_c; share = np.abs(eta) * np.linalg.norm(d, axis=0); share /= share.sum()
top = np.argsort(share)[::-1][:6]
f = np.geomspace(60, 2000, 3000); w = 2 * np.pi * f; den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
fig = plt.figure(figsize=(17, 9.6)); gsp = fig.add_gridspec(2, 2, height_ratios=[0.42, 1], width_ratios=[1.25, 1])
from matplotlib.patches import FancyBboxPatch
axf = fig.add_subplot(gsp[0, :]); axf.set_xlim(0, 17); axf.set_ylim(0, 2.6); axf.axis("off")
steps = [("線図を作る", "各回転数で最小の $a_{\\lim}$ を選ぶ\n（3,890 min$^{-1}$ で 9.2 mm）", "#7f8c8d"),
         ("周波数が 1 つ決まる", "その $a_{\\lim}$ を出したときの\nびびり周波数 $\\omega_c$ = 185 Hz", "#c0392b"),
         ("その周波数の FRF", "$G(\\omega_c)=\\sum_r$（各モードの分）\n全モードの和（①）", "#1f618d"),
         ("揺れの形 → 割合", "行列式 0 の向き $\\boldsymbol{Q}$ → 力 $\\boldsymbol{F}$\n→ 各モードの $\\eta_r$ → 割合（②）", "#27ae60")]
for k, (h, b, c) in enumerate(steps):
    x = 0.2 + k * 4.2
    axf.add_patch(FancyBboxPatch((x, 0.2), 3.7, 2.1, boxstyle="round,pad=0.05", fc=c, alpha=0.12, ec=c, lw=2))
    axf.text(x + 1.85, 1.85, h, ha="center", fontsize=14, weight="bold", color=c); axf.text(x + 1.85, 0.95, b, ha="center", va="center", fontsize=12)
    if k < 3: axf.annotate("", xy=(x + 4.15, 1.25), xytext=(x + 3.8, 1.25), arrowprops=dict(arrowstyle="-|>", lw=2.5, mutation_scale=22))
axf.set_title("線図で選ぶのは「1 つの周波数」で、「1 つのモード」ではない。その周波数で揺れると、構造は全モードの足し合わせで応答する", fontsize=14, loc="left")
axs = [fig.add_subplot(gsp[1, 0]), fig.add_subplot(gsp[1, 1])]
ax = axs[0]; cols = plt.cm.tab10(np.arange(6))
ax.loglog(f, np.linalg.norm(np.einsum("ir,jr,wr->wij", d, d, 1 / den).reshape(len(f), -1), axis=1) * 1e3, color="k", lw=3, label="相対 FRF 全体（X・Y の 2×2、全モードの合計）")
for c, r in zip(cols, sorted(top[:5], key=lambda r: fr[r])):
    ax.loglog(f, (d[:, r] @ d[:, r]) / np.abs(den[:, r]) * 1e3, color=c, lw=1.6, ls="--", label=f"モード {r + 1}（{fr[r]:.0f} Hz）だけの分")
ax.axvline(fc, color="#c0392b", lw=2.5); ax.text(fc * 1.04, 3e-1, "線図の最小点で選ばれた\nびびり周波数 185 Hz", color="#c0392b", fontsize=12.5, va="top")
ax.set_xlim(60, 2000); ax.set_ylim(1e-5, 1); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=10.5, loc="lower left")
ax.set_xlabel("周波数 [Hz]", fontsize=13); ax.set_ylabel("コンプライアンス [µm/N]", fontsize=13)
ax.set_title("① 185 Hz の縦線の上で、どのモードの線も 0 ではない\n→ 185 Hz で揺らすと、全モードが少しずつ一緒に揺れる", fontsize=13.5, loc="left")
ax = axs[1]
lab = [f"モード {r + 1}\n{fr[r]:.0f} Hz" for r in top]; ax.bar(range(len(top)), share[top] * 100, color=[cols[sorted(top[:5], key=lambda r: fr[r]).index(r)] if r in top[:5] else "gray" for r in top])
for k, r in enumerate(top): ax.text(k, share[r] * 100 + 1.5, f"{share[r] * 100:.0f} %", ha="center", fontsize=12.5)
ax.set_xticks(range(len(top))); ax.set_xticklabels(lab, fontsize=11.5); ax.set_ylabel("びびりの揺れに占める割合 [%]", fontsize=13); ax.set_ylim(0, 100)
ax.set_title("② 185 Hz のびびりの揺れを、モードごとに分けた割合\n（共振に近いモード 8 が大半、残りは離れたモードの裾）", fontsize=13.5, loc="left"); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(R / "mode_mix_why.png", dpi=90); print([(r + 1, round(fr[r]), round(share[r] * 100, 1)) for r in top])
