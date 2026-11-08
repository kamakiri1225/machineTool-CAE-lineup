"""5.9 と 6.0 の比較図（compare.py の後に実行）"""
import os, sys, re
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(D, "..", "..", "study", "scripts"))
from figstyle import plt, BLUE, RED, GRAY
a, b = np.load(os.path.join(D, "frf_59.npy")), np.load(os.path.join(D, "frf_60.npy"))
rows = [l for l in open(os.path.join(D, "compare.md")) if l.startswith("| ana")]
names = [r.split("|")[1].strip() for r in rows]; t59 = [float(r.split("|")[2]) for r in rows]; t60 = [float(r.split("|")[3]) for r in rows]
fig, axs = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw=dict(width_ratios=[1.3, 1]))
ax = axs[0]
ax.semilogy(a[:, 0], a[:, 1] * 1e3, color=BLUE, lw=3, label="FrontISTR 5.9")
ax.semilogy(b[:, 0], b[:, 1] * 1e3, "--", color=RED, lw=1.8, label="FrontISTR 6.0")
ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("コンプライアンス [µm/N]"); ax.grid(alpha=0.3, which="both"); ax.legend()
ax.set_title(f"ana002 X-X の周波数応答：曲線は完全に重なる（最大差 {np.max(np.abs(a[:,1]-b[:,1])/a[:,1])*100:.2e} %）", loc="left", fontsize=12)
ax = axs[1]; x = np.arange(4); w = 0.38
ax.bar(x - w / 2, t59[:4], w, color=BLUE, label="5.9"); ax.bar(x + w / 2, t60[:4], w, color=RED, label="6.0")
for i in range(4):
    ax.text(i + w / 2, t60[i] + 0.5, f"{(t60[i]/t59[i]-1)*100:+.0f} %", ha="center", fontsize=10, color=RED)
ax.set_xticks(x); ax.set_xticklabels(names[:4], fontsize=10); ax.set_ylabel("計算時間 [s]（OpenMP 8 スレッド）")
ax.set_title("計算時間（同じ入力・同じ PC）", loc="left", fontsize=12); ax.legend(); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(os.path.join(D, "v59_v60_compare.png")); print("ok")
