"""ana002 後処理: 加振点のコンプライアンス（変位/力）を CalculiX と FrontISTR で比較"""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from plot_modes import plt
from make_inputs import ZETA, FA, FB, ALPHA, BETA

D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, "results"); os.makedirs(R, exist_ok=True)


def read_ccx(path):
    """*NODE PRINT の出力: 周波数ごとに実部・虚部の 2 ブロック。ステップが変わると周波数が戻る"""
    lines = open(path).read().splitlines()
    rec = []
    for i, ln in enumerate(lines):
        m = re.match(r"\s*displacements \(vx,vy,vz\) for set MON and time\s+([\d.E+-]+)", ln)
        if m:
            v = [float(x) for x in lines[i + 2].split()[1:4]]
            rec.append((float(m.group(1)), v))
    steps, cur = [], []
    for k in range(0, len(rec), 2):
        f, re_ = rec[k]; _, im = rec[k + 1]
        if cur and f < cur[-1][0]:
            steps.append(cur); cur = []
        cur.append((f, np.array(re_) + 1j * np.array(im)))
    steps.append(cur)
    return [(np.array([c[0] for c in s]), np.array([c[1] for c in s])) for s in steps]


def read_fistr(path):
    f, a = [], []
    for ln in open(path):
        m = re.match(r"\s*([\d.E+-]+)\s+\[Hz\]\s*:\s*([\d.E+-]+)", ln)
        if m:
            f.append(float(m.group(1))); a.append(float(m.group(2)))
    return np.array(f), np.array(a)


def eig_ccx(path):
    out = []
    txt = open(path).read()
    blk = txt.split("E I G E N V A L U E   O U T P U T")[1].split("P A R T I C I P A T I O N")[0]
    for ln in blk.splitlines():
        p = ln.split()
        if len(p) == 5 and p[0].isdigit():
            out.append(float(p[3]))
    return np.array(out)


ccx = read_ccx(os.path.join(D, "ccx", "machine.dat"))
fe = eig_ccx(os.path.join(D, "ccx", "machine.dat"))
names = {1: "X（送り方向）", 2: "Y（横方向）", 3: "Z（上下方向）"}
fig, axs = plt.subplots(3, 1, figsize=(11, 11), sharex=True)
summary = ["| 加振方向 | 静コンプライアンス [µm/N]（50 Hz） | 静剛性 [N/µm] | 最大ピーク [Hz] | ピーク値 [µm/N] | 動剛性の最小 [N/µm] |", "|---|---|---|---|---|---|"]
for k, dof in enumerate((1, 2, 3)):
    ax = axs[k]
    f, U = ccx[k]
    amp_dir = np.abs(U[:, dof - 1]) * 1000          # µm/N（加振方向成分）
    amp_all = np.linalg.norm(np.abs(U), axis=1) * 1000
    ff, af = read_fistr(os.path.join(D, "fistr", f"out_{dof}", "0.log"))
    ax.semilogy(f, amp_dir, "-", color="#1f618d", lw=2, label="CalculiX（加振方向の成分）")
    ax.semilogy(ff, af * 1000, "--", color="#c0392b", lw=1.6, label="FrontISTR（モニタ節点の変位振幅）")
    for e in fe[fe < 1000]:
        ax.axvline(e, color="gray", lw=0.6, ls=":")
    ax.set_ylabel("コンプライアンス [µm/N]")
    ax.set_title(f"テーブル上面中央を {names[dof]} に 1 N で加振 → 同じ点の {names[dof][0]} 方向変位", loc="left", fontsize=12)
    ax.grid(alpha=0.3, which="both"); ax.legend(loc="upper right", fontsize=10)
    i = np.argmax(amp_dir)
    summary.append(f"| {names[dof]} | {amp_dir[0]:.4f} | {1/amp_dir[0]:.0f} | {f[i]:.0f} | {amp_dir[i]:.4f} | {1/amp_dir[i]:.0f} |")
    ax.annotate(f"{f[i]:.0f} Hz", xy=(f[i], amp_dir[i]), xytext=(f[i] + 40, amp_dir[i] * 0.55), fontsize=10,
                arrowprops=dict(arrowstyle="->"))
axs[-1].set_xlabel("周波数 [Hz]（灰色の点線 = 固有振動数）")
axs[-1].set_xlim(50, 1000)
fig.suptitle(f"周波数応答（Rayleigh 減衰：{FA:.0f} Hz と {FB:.0f} Hz で ζ = {ZETA*100:.0f} %）", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(R, "frf_compare.png")); plt.close(fig)

# 減衰比の周波数特性
w = 2 * np.pi * np.linspace(20, 1200, 400)
fig, ax = plt.subplots(figsize=(8, 3.6))
ax.plot(w / 2 / np.pi, 100 * (ALPHA / (2 * w) + BETA * w / 2), color="#1f618d", lw=2)
ax.axhline(ZETA * 100, color="gray", ls=":"); ax.plot([FA, FB], [ZETA * 100] * 2, "o", color="#c0392b")
ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("減衰比 ζ [%]"); ax.set_ylim(0, 8); ax.grid(alpha=0.3)
ax.set_title(f"Rayleigh 減衰 $\\zeta = \\alpha/(2\\omega) + \\beta\\omega/2$　（$\\alpha$={ALPHA:.1f}, $\\beta$={BETA:.2e}）", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(R, "rayleigh_damping.png")); plt.close(fig)
open(os.path.join(R, "summary.md"), "w").write("\n".join(summary) + "\n")
print("\n".join(summary))
