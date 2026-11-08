"""5.2 補足：主軸の正転・逆転、ねじれの向きは安定限界に効くか（図解と計算）"""
import json, sys
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
D = Path(__file__).resolve().parent; R = D / "results"
sys.path[:0] = [str(D), str(D.parent / "tools"), str(D.parent.parent / "study" / "scripts")]
from analyze import boundary
from figstyle import plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Wedge, Polygon
TEETH, KT, KR = 4, 1000.0, 0.3
z = np.load(R / "foundation_modal_basis.npz"); fr = z["frequency_hz"]; phi = z["modes"]; iw, it = int(z["work_index"]), int(z["tool_index"])
wr = 2 * np.pi * fr; zr = z["alpha"] / (2 * wr) + z["beta"] * wr / 2; d = phi[iw] - phi[it]
grid = np.linspace(5, 3000, 24000); rpm = np.linspace(500, 16000, 600)
w = 2 * np.pi * grid; G = np.einsum("ir,jr,wr->wij", d, d, 1 / (wr[None] ** 2 - w[:, None] ** 2 + 2j * zr[None] * wr[None] * w[:, None]))


def lobes(s, ka):
    B = np.zeros((3, 2)); B[:2] = TEETH / 4 * np.array([[KR, s], [-s, KR]]); B[2] = TEETH / np.pi * ka * np.array([1.0, 0.0])
    lam = np.linalg.eigvals(G[:, :2, :] @ B)
    for i in range(1, len(lam)):
        if np.sum(abs(lam[i] - lam[i - 1])) > np.sum(abs(lam[i, ::-1] - lam[i - 1])): lam[i] = lam[i, ::-1]
    br = [boundary(interp1d(grid, lam[:, j], bounds_error=True), rpm, grid, Kc=KT, teeth=TEETH) for j in (0, 1)]
    v = np.stack([b[0] for b in br], 1); return np.nanmin(np.where(np.isfinite(v), v, np.nan), 1)


fig = plt.figure(figsize=(17, 10.5)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.9])
BL, RD, GR, OR = "#1f618d", "#c0392b", "#7f8c8d", "#f5b041"


def top(ax, s, title):
    ax.add_patch(Rectangle((-2.6, -1.6), 5.2, 3.2, fc="#d5d8dc", ec="none"))                 # ワーク（全体）
    ax.add_patch(Rectangle((-2.6, -1.0), 2.6, 2.0, fc="white", ec="none"))                    # 削り終わった溝（空間）
    ax.plot([-2.6, 0], [1, 1], "k", lw=1); ax.plot([-2.6, 0], [-1, -1], "k", lw=1)
    ax.text(-1.6, 0.0, "削り終わった溝\n（空間）", ha="center", va="center", fontsize=11, color="#555")
    ax.text(-1.3, 1.25, "溝の壁（ワーク）", ha="center", fontsize=10.5, color="#555"); ax.text(-1.3, -1.4, "溝の壁（ワーク）", ha="center", fontsize=10.5, color="#555")
    ax.add_patch(Wedge((0, 0), 1.0, -90, 90, fc="#f5b041", alpha=0.7, ec="none"))           # 刃が削っている所
    ax.text(1.75, -0.95, "これから削る所", ha="center", fontsize=11, color="#555")
    ax.add_patch(Circle((0, 0), 0.93, fc="#5d6d7e", ec="k"))
    ax.add_patch(FancyArrowPatch((0.55 * np.sin(-0.6), 0.55 * np.cos(-0.6) * s), (0.55 * np.sin(0.9), 0.55 * np.cos(0.9) * s), connectionstyle=f"arc3,rad={-0.5 * s}",
                                 arrowstyle="-|>", mutation_scale=22, color="white", lw=2.5))
    ax.text(0, -0.15, "時計回り" if s > 0 else "反時計回り", color="white", ha="center", fontsize=12)
    p = np.array([1.0, 0.0]); f = np.array([KR, -s]) * 1.1
    ax.annotate("", xy=p + f, xytext=p, arrowprops=dict(arrowstyle="-|>", color=RD, lw=3, mutation_scale=22))
    ax.text(*(p + f + ([-0.2, -0.3] if s > 0 else [0.15, -0.35])), "ワークが受ける力", color=RD, fontsize=12.5, va="center")
    ax.annotate("", xy=(2.3, 1.35), xytext=(1.3, 1.35), arrowprops=dict(arrowstyle="-|>", color="green", lw=2.5, mutation_scale=18)); ax.text(1.3, 1.45, "送り（X）", color="green", fontsize=12)
    ax.text(-2.5, -1.75, "刃が X 側で削っているとき\n力は " + ("−Y（下）" if s > 0 else "+Y（上）") + " 向き", fontsize=12.5, va="top")
    ax.set_xlim(-2.6, 3.4); ax.set_ylim(-2.6, 1.9); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=14.5, loc="left")
    ax.text(-2.5, 1.7, "上から見た図（灰：ワーク、白：削った溝、橙：刃が削っている所）", fontsize=10.5, color="#555")


top(fig.add_subplot(gs[0, 0]), 1, "(a) 正転：X に揺れると −Y に力")
top(fig.add_subplot(gs[0, 1]), -1, "(b) 逆転：X に揺れると +Y に力（鏡に映した形）")

ax = fig.add_subplot(gs[0, 2])
for x0, up, lab in [(-1.3, True, "ねじれと回転の\n組み合わせ ①"), (1.3, False, "組み合わせ ②\n（どちらかを逆に）")]:
    ax.add_patch(Rectangle((x0 - 0.45, 0.0), 0.9, 2.6, fc="#5d6d7e", ec="k"))
    for yy in np.arange(0.1, 2.5, 0.45): ax.plot([x0 - 0.45, x0 + 0.45], [yy, yy + 0.35] if up else [yy + 0.35, yy], color="white", lw=2)
    ax.annotate("", xy=(x0, 3.4 if up else -0.8), xytext=(x0, 2.7 if up else -0.1), arrowprops=dict(arrowstyle="-|>", color=RD, lw=3, mutation_scale=22))
    ax.text(x0 + 0.15, 3.1 if up else -0.5, "$F_z$", color=RD, fontsize=15)
    ax.text(x0, -1.3, lab, ha="center", va="top", fontsize=12)
ax.add_patch(Rectangle((-2.4, -0.1), 4.8, 0.1, fc="#d5d8dc", ec="none"))
ax.set_xlim(-2.6, 2.6); ax.set_ylim(-2.5, 3.8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(c) ねじれ刃：$F_z$ が上か下か（横から見た図）", fontsize=14.5, loc="left")

ax = fig.add_subplot(gs[1, :]); out = {}
for s, ka, c, ls, lab in [(1, 0, "k", "-", "2×2 正転"), (-1, 0, "#f39c12", "--", "2×2 逆転（正転と完全に重なる）"),
                          (1, 0.3, RD, "-", "3×3 ねじれ刃 $F_z$ 向き ①"), (1, -0.3, BL, "-", "3×3 ねじれ刃 $F_z$ 向き ②")]:
    l = lobes(s, ka); out[lab] = dict(min=float(np.nanmin(l)), rpm=float(rpm[np.nanargmin(l)])); ax.plot(rpm, l, color=c, ls=ls, lw=2.2 if ls == "--" else 1.8, label=f"{lab}：最小 {np.nanmin(l):.2f} mm")
ax.set_ylim(0, 30); ax.grid(alpha=0.3); ax.legend(fontsize=11.5, loc="upper left"); ax.set_xlabel("主軸回転数 [min⁻¹]", fontsize=12); ax.set_ylabel("軸方向切込み限界 [mm]", fontsize=12)
ax.set_title("(d) 計算結果：このモデルでは正転・逆転は同じ（X と Y が連成しない $G_{xy}\\approx0$ ため）。$F_z$ の向きは ±1 割効く", fontsize=14, loc="left")
fig.tight_layout(); fig.savefig(R / "rotation_direction.png", dpi=95); json.dump(out, open(R / "rotation_direction.json", "w"), ensure_ascii=False, indent=1); print(out)
