"""ana002 後処理（追加）: クロス方向を含む 3×3 の周波数応答（加振方向 i → 揺れ方向 j）と、
   各ピークがどの固有モード（ana001）から来ているか（モード寄与）を図にする。
   H_ij(ω) = Σ_r φ_ri φ_rj / (ω_r² − ω² + 2 i ζ_r ω_r ω)   （φ は質量正規化した固有モードの MON 点の成分）"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from plot_modes import plt
from post import read_ccx, read_fistr
from make_inputs import ALPHA, BETA
import machine_model as mm

D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, "results")
MON = mm.build().monitor
DIR = "XYZ"
NAMES = {1: "コラム Y 曲げ", 2: "テーブル X 並進", 3: "コラム X 曲げ", 4: "テーブル Y 並進",
         5: "テーブル Z ねじれ", 6: "テーブル Z＋曲げ"}
COL = ["#1f618d", "#c0392b", "#1e8449", "#e67e22", "#8e44ad", "#7f8c8d", "#17a589", "#b7950b"]


CORNER = 22255      # テーブル上面の角 (750, 500, 350)：中心から X +250, Y +200 ずれた点


def read_modes(path, nodes, nmax=30):
    """frd の先頭の固有モード（DISP ブロック）から指定節点の成分と固有振動数を読む"""
    fr, phi = [], {n: [] for n in nodes}
    with open(path) as f:
        freq = None; disp = False; cnt = 0
        for ln in f:
            if ln.startswith("  100CL"):
                freq = float(ln[12:24]); disp = False
            elif ln.startswith(" -4  DISP"):
                disp = True; fr.append(freq); cnt += 1
            elif disp and ln.startswith(" -1") and int(ln[3:13]) in phi:
                phi[int(ln[3:13])].append([float(ln[13:25]), float(ln[25:37]), float(ln[37:49])])
            elif disp and ln.startswith(" -3"):
                disp = False
                if cnt == nmax:
                    break
    return np.array(fr), {n: np.array(v) for n, v in phi.items()}


fr, PH = read_modes(os.path.join(D, "ccx", "machine.frd"), [MON, CORNER])
phi, phc = PH[MON], PH[CORNER]
wr = 2 * np.pi * fr
zr = ALPHA / (2 * wr) + BETA * wr / 2
ccx = read_ccx(os.path.join(D, "ccx", "machine.dat"))
f = ccx[0][0]; w = 2 * np.pi * f
Hccx = np.zeros((3, 3, len(f)), complex)
for i in range(3):
    Hccx[i] = ccx[i][1].T                       # [加振 i][揺れ j]，単位 mm/N
contrib = phi[:, :, None, None] * phi[:, None, :, None] / (wr[:, None, None, None]**2 - w**2 + 2j * zr[:, None, None, None] * wr[:, None, None, None] * w)
Hmod = contrib.sum(0)
contc = phc[:, :, None, None] * phc[:, None, :, None] / (wr[:, None, None, None]**2 - w**2 + 2j * zr[:, None, None, None] * wr[:, None, None, None] * w)
Hcor = contc.sum(0)                              # 角で加振・角で測る（モード重ね合わせで計算）
err = np.abs(np.abs(Hmod) - np.abs(Hccx)).max() / np.abs(Hccx).max()
print("modes", len(fr), "max rel err modal-sum vs ccx:", err)

# ------------------------------------------------ 1. 3×3 のマトリクス
fig, axs = plt.subplots(3, 3, figsize=(15, 11), sharex=True, sharey=True)
for i in range(3):
    ff, af = read_fistr(os.path.join(D, "fistr", f"out_{i+1}", "0.log"))
    for j in range(3):
        ax = axs[i, j]
        ax.semilogy(f, np.abs(Hccx[i, j]) * 1e3, color="#1f618d", lw=2.2, label="中心 MON（CalculiX）")
        ax.semilogy(f, np.abs(Hcor[i, j]) * 1e3, color="#e67e22", lw=1.8, label="角 CORNER（モード重ね合わせ）")
        if i == j:
            ax.semilogy(ff, af * 1e3, "--", color="k", lw=1, label="中心 MON（FrontISTR）")
        for r in range(6):
            ax.axvline(fr[r], color="gray", lw=0.6, ls=":")
        ax.set_title(f"{DIR[i]}-{DIR[j]}：{DIR[i]} 方向に加振 → {DIR[j]} 方向の揺れ" + ("（直接）" if i == j else "（クロス）"),
                     fontsize=11, loc="left", color="k" if i == j else "#c0392b")
        ax.grid(alpha=0.3, which="both")
        if i == 0 and j == 0:
            ax.legend(fontsize=9, loc="lower left")
for j in range(3):
    axs[2, j].set_xlabel("周波数 [Hz]")
    axs[i, 0].set_ylabel("コンプライアンス [µm/N]")
for i in range(3):
    axs[i, 0].set_ylabel("コンプライアンス [µm/N]")
for r in range(6):
    axs[0, 0].text(fr[r], 3e1, f"{r+1}", ha="center", fontsize=9, color="gray")
axs[0, 0].set_ylim(1e-5, 1e2); axs[0, 0].set_xlim(50, 1000)
fig.suptitle("3×3 周波数応答（加振方向 − 揺れ方向、加振点 = 測定点）。青 = テーブル上面中央、橙 = テーブル上面の角。灰色の点線と番号 = 固有モード",
             fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(R, "frf_matrix_3x3.png")); plt.close(fig)

# ------------------------------------------------ 2. 固有モードが MON をどの向きに動かすか
n = 8
fig, axs2 = plt.subplots(1, 2, figsize=(14, 6))
for ax, P, ttl in [(axs2[0], phi, "テーブル中央 MON"), (axs2[1], phc, "テーブルの角 CORNER")]:
    A = np.abs(P[:n]); A = A / max(np.abs(phi[:n]).max(), np.abs(phc[:n]).max())
    im = ax.imshow(A, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for r in range(n):
        for c in range(3):
            ax.text(c, r, f"{A[r, c]:.2f}", ha="center", va="center", fontsize=11, color="white" if A[r, c] > 0.5 else "k")
    ax.set_xticks(range(3)); ax.set_xticklabels(["X 成分", "Y 成分", "Z 成分"], fontsize=11)
    ax.set_yticks(range(n)); ax.set_yticklabels([f"{r+1}次 {fr[r]:.0f} Hz " + NAMES.get(r + 1, "") for r in range(n)], fontsize=10)
    ax.set_title(f"{ttl}：各モードで点がどの向きに動くか", loc="left", fontsize=12)
fig.suptitle("$|\\phi_{r,j}|$（質量正規化モードの成分、2 点の最大を 1）。$H_{ij}$ のピーク高さは $\\phi_{r,i}\\,\\phi_{r,j}$ に比例", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(R, "mode_participation_mon.png")); plt.close(fig)

# ------------------------------------------------ 3. ピークのモード寄与（直接 3 つ＋代表的なクロス 2 つ）
pairs = [(0, 0, "MON"), (1, 1, "MON"), (2, 2, "MON"), (0, 1, "CORNER"), (1, 2, "CORNER")]
fig, axs = plt.subplots(len(pairs), 1, figsize=(12, 17), sharex=True)
for ax, (i, j, pt) in zip(axs, pairs):
    H_, C_ = (Hccx, contrib) if pt == "MON" else (Hcor, contc)
    ax.semilogy(f, np.abs(H_[i, j]) * 1e3, color="k", lw=2.5, label="合計")
    peak = np.abs(C_[:, i, j, :]).max(1)
    top = np.argsort(peak)[::-1][:4]
    for c, r in zip(COL, sorted(top)):
        ax.semilogy(f, np.abs(C_[r, i, j]) * 1e3, "--", color=c, lw=1.6,
                    label=f"{r+1}次 {fr[r]:.0f} Hz " + NAMES.get(r + 1, "") + " の寄与")
    ax.set_ylim(np.abs(H_[i, j]).max() * 1e3 * 1e-4, np.abs(H_[i, j]).max() * 1e3 * 3)
    ax.set_title(f"{DIR[i]}-{DIR[j]}（{DIR[i]} に加振 → {DIR[j]} の揺れ）＠" + ("テーブル中央" if pt == "MON" else "テーブルの角"), loc="left", fontsize=12)
    ax.legend(fontsize=9, loc="upper right", ncol=2); ax.grid(alpha=0.3, which="both"); ax.set_ylabel("µm/N")
axs[-1].set_xlabel("周波数 [Hz]"); axs[-1].set_xlim(50, 1000)
fig.suptitle("周波数応答はモードの足し算：$H_{ij}(\\omega) = \\sum_r \\dfrac{\\phi_{r,i}\\,\\phi_{r,j}}{\\omega_r^2 - \\omega^2 + 2i\\zeta_r\\omega_r\\omega}$",
             fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(R, "frf_modal_contribution.png")); plt.close(fig)

# ------------------------------------------------ 4. 表（ピーク）
L = []
for title, H_, C_ in [("テーブル上面中央 MON (500, 300, 350)", Hccx, contrib), ("テーブル上面の角 CORNER (750, 500, 350)", Hcor, contc)]:
    L += [f"### {title}", "", "| 加振 − 揺れ | 50 Hz の値 [µm/N] | 最大ピーク [Hz] | ピーク値 [µm/N] | 主に効いているモード |", "|---|---|---|---|---|"]
    for i in range(3):
        for j in range(3):
            a = np.abs(H_[i, j]) * 1e3; k = np.argmax(a)
            r = np.argmax(np.abs(C_[:, i, j, k]))
            mode = f"{r+1}次 {fr[r]:.0f} Hz {NAMES.get(r+1, '')}" if a[k] > 1e-6 else "―（対称のため 0）"
            L.append(f"| {DIR[i]}-{DIR[j]} | {a[0]:.2e} | {f[k]:.0f} | {a[k]:.2e} | {mode} |")
    L.append("")
open(os.path.join(R, "frf_matrix_summary.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
