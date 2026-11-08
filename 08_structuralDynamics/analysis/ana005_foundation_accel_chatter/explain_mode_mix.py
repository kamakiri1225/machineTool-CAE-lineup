"""ana005 5.3 補足：びびり振動の形が、なぜいくつかのモードの重ね合わせになるか（線図の最小点で数値を示す）"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent.parent / "study" / "scripts")]
import chatter_by_rpm as C
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
k = int(np.nanargmin(C.lim)); n, a, f = C.rpm[k], C.lim[k], C.fc[k]
w = 2 * np.pi * f; den = C.wr**2 - w**2 + 2j * C.zr * C.wr * w
G = np.einsum("ir,jr,r->ij", C.d, C.d, 1 / den); tau = 60 / (C.TEETH * n)
M = np.eye(2) + C.KT * a * (G @ C.B0) * (1 - np.exp(-1j * w * tau)); ev, V = np.linalg.eig(M); q = V[:, np.argmin(abs(ev))]
q = q / q[np.argmax(abs(q))]; F = C.B0 @ q
gain = np.abs(1 / den)                       # モードの増幅（共振からの近さ）
drive = np.abs(C.d.T @ F)                    # 力がそのモードをどれだけ押せるか（形と力の向き）
eta = C.d.T @ F / den; share = np.abs(eta) * np.linalg.norm(C.d, axis=0); share /= share.sum()
top = np.argsort(share)[::-1][:6]
fig, axs = plt.subplots(1, 4, figsize=(17, 4.8))
lab = [f"{r+1}\n{C.fr[r]:.0f} Hz" for r in top]
for ax, v, t, c in zip(axs, (gain[top] / gain[top].max(), drive[top] / drive[top].max(), share[top]),
                        ("① 共振からの近さ\n$1/|\\omega_r^2-\\omega_c^2+2i\\zeta_r\\omega_r\\omega_c|$", "② 力の向きとモードの形の合い方\n$|\\boldsymbol{d}_r^{T}\\boldsymbol{F}|$", "③ ①×②：びびり振動に占める割合"), (BLUE, ORANGE, RED)):
    ax.bar(range(len(top)), v, color=c); ax.set_xticks(range(len(top))); ax.set_xticklabels(lab, fontsize=10); ax.set_title(t, loc="left", fontsize=12); ax.grid(alpha=0.3, axis="y")
axs[0].set_ylabel("最大を 1 とした比"); axs[2].set_ylabel("割合")
ax = axs[3]; ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.axis("off")
ax.annotate("", xy=(q.real[0], q.real[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=3))
Fr = F / np.abs(F).max()
ax.annotate("", xy=(Fr.real[0], Fr.real[1]), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=RED, lw=3))
ax.text(q.real[0], q.real[1] + 0.12, "揺れ $q$", color=BLUE, fontsize=13, ha="center"); ax.text(Fr.real[0], Fr.real[1] - 0.2, "切削力 $F=B_0q$", color=RED, fontsize=13, ha="center")
ax.axhline(0, color=GRAY, lw=0.5); ax.axvline(0, color=GRAY, lw=0.5); ax.text(1.25, 0.05, "X", fontsize=12); ax.text(0.05, 1.2, "Y", fontsize=12)
ax.set_title("④ X に揺れると、力は Y にも出る", loc="left", fontsize=12)
fig.suptitle(f"線図の最小点（n = {n:.0f} min⁻¹、$a_p$ = {a:.2f} mm、びびり {f:.0f} Hz）で、各モードがびびり振動にどれだけ混ざるか", fontsize=13)
fig.tight_layout(); fig.savefig(R / "mode_mix.png")
print("q", np.round(q, 3), "F", np.round(F / np.abs(F).max(), 3)); [print(r + 1, round(C.fr[r]), round(gain[r] / gain[top].max(), 3), round(drive[r] / drive[top].max(), 3), round(share[r], 3)) for r in top]
