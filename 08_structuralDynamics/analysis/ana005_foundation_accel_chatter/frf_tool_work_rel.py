"""ana005 3.3：刃先・ワーク・相対の 3 つのコンプライアンスを並べる（テストバー、3 つの支持条件）
   各ケースのモード基底（results/{case}_modal_basis.npz）から
     刃先 H_tt（刃先だけを叩いて刃先で測る）、ワーク H_ww（ワークだけを叩いてワークで測る）、相対 G = H_ww + H_tt − H_wt − H_tw
   出力：frf_tool_work_rel.png（XX・YY・ZZ × 3 ケース）、relative_matrix.png（基礎ばねの 9 成分、3 種類を重ねる）"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent; R = D / "results"
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt
CASES = [("fixed", "底面固定"), ("foundation", "基礎ばね"), ("foundation_table_stiff", "基礎ばね・テーブル結合 10 倍")]
f = np.geomspace(5, 3000, 5000); w = 2 * np.pi * f
COL = {"tool": "#1f77b4", "work": "#ff7f0e", "rel": "#d62728"}
LAB = {"tool": "刃先（テストバー先端）", "work": "ワーク（テーブル上）", "rel": "刃先―ワーク相対"}


def frfs(case):
    z = np.load(R / f"{case}_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
    wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2
    den = wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]
    H = lambda a, b: np.einsum("ir,jr,wr->wij", a, b, 1 / den)
    pt, pw = phi[it], phi[iw]
    return {"tool": H(pt, pt), "work": H(pw, pw), "rel": H(pw - pt, pw - pt)}


data = {c: frfs(c) for c, _ in CASES}
fig, axs = plt.subplots(3, 3, figsize=(17, 12), sharex=True)
for j, (c, lab) in enumerate(CASES):
    for i, d in enumerate("XYZ"):
        ax = axs[i, j]
        for k in ("tool", "work", "rel"):
            ax.loglog(f, np.abs(data[c][k][:, i, i]) * 1e3, color=COL[k], lw=2.2 if k == "rel" else 1.5, ls="-" if k == "rel" else "--", label=LAB[k])
        ax.set_xlim(5, 3000); ax.set_ylim(1e-4, 1); ax.grid(alpha=0.3, which="both")
        ax.set_title(f"{lab}：{d}-{d}", loc="left", fontsize=12)
        if j == 0: ax.set_ylabel("コンプライアンス [µm/N]")
        if i == 2: ax.set_xlabel("周波数 [Hz]")
axs[0, 0].legend(fontsize=10, loc="lower left")
fig.suptitle("刃先・ワーク・相対のコンプライアンス（テストバー 60×60×100 mm、点線 = 片側だけを叩いたとき、実線 = 刃先とワークの相対）", fontsize=13)
fig.tight_layout(); fig.savefig(R / "frf_tool_work_rel.png"); plt.close(fig)
d = data["foundation"]
fig, axs = plt.subplots(3, 3, figsize=(17, 12), sharex=True)
for i in range(3):
    for j in range(3):
        ax = axs[i, j]
        for k in ("tool", "work", "rel"):
            ax.loglog(f, np.abs(d[k][:, i, j]) * 1e3, color=COL[k], lw=2.2 if k == "rel" else 1.5, ls="-" if k == "rel" else "--", label=LAB[k])
        ax.set_xlim(5, 3000); ax.set_ylim(1e-5, 1); ax.grid(alpha=0.3, which="both")
        ax.set_title(f"{'XYZ'[i]}{'XYZ'[j]}：{'XYZ'[j]} 方向に叩く → {'XYZ'[i]} 方向の揺れ", loc="left", fontsize=11)
        if j == 0: ax.set_ylabel("µm/N")
        if i == 2: ax.set_xlabel("周波数 [Hz]")
axs[0, 0].legend(fontsize=9, loc="lower left")
fig.suptitle("基礎ばねモデルの 9 成分：刃先・ワーク・相対（テストバー）", fontsize=13)
fig.tight_layout(); fig.savefig(R / "relative_matrix.png"); plt.close(fig)
for c, lab in CASES:
    for k in ("tool", "work", "rel"):
        a = np.abs(data[c][k][:, 0, 0]); m = f > 100; i = np.argmax(a[m])
        print(lab, k, f"X-X 100 Hz 以上の最大 {a[m][i]*1e3:.3f} µm/N @ {f[m][i]:.0f} Hz, 静的 {a[0]*1e3:.4f}")
