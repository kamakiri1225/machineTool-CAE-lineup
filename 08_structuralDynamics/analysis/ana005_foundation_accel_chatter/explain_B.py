"""ana005 5.2：方向係数 B(φ) と B0 がどう出てくるかを 4 段の図で示す
   ① 刃の位置と向き（e_r, e_t）  ② 食い違いが削る厚さを変える Δh = e_r·Δq  ③ 力の向き ΔF = K_t a_p Δh (e_t + K_r e_r)
   ④ φ を 0〜π で動かしたときの B の 4 成分と、その平均（= B0）"""
import sys
import numpy as np
from pathlib import Path
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import Circle, Wedge, FancyArrowPatch, Arc
Kr, Z = 0.3, 4
arr = lambda ax, a, b, c, lw=3: ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=22, color=c, lw=lw))
phi = np.radians(55); er = np.array([np.sin(phi), np.cos(phi)]); et = np.array([np.cos(phi), -np.sin(phi)])
fig = plt.figure(figsize=(20, 12))
# ① 
ax = fig.add_axes([0.0, 0.5, 0.25, 0.46]); ax.set_xlim(-1.6, 1.8); ax.set_ylim(-1.6, 1.8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("① 刃の位置と 2 つの向き", loc="left", fontsize=19)
ax.add_patch(Circle((0, 0), 1, fill=False, lw=2)); ax.add_patch(Wedge((0, 0), 1, -90, 90, fc=ORANGE, alpha=0.25))
ax.plot([0, er[0]], [0, er[1]], "k", lw=2); ax.plot(*er, "o", color="k", ms=9)
ax.plot([0, 0], [0, 1.3], color=GRAY, ls="--"); ax.add_patch(Arc((0, 0), 0.7, 0.7, theta1=90 - 55, theta2=90, color="k", lw=1.5)); ax.text(0.12, 0.45, r"$\phi$", fontsize=20)
arr(ax, er, er + 0.6 * er, BLUE); ax.text(*(er + 0.62 * er + [0.05, 0]), r"$\boldsymbol{e}_r=\left(\sin\phi,\ \cos\phi\right)$" + "\n半径方向", color=BLUE, fontsize=15)
arr(ax, er, er + 0.6 * et, RED); ax.text(*(er + 0.6 * et + [-0.9, -0.75]), r"$\boldsymbol{e}_t=\left(\cos\phi,\ -\sin\phi\right)$" + "\n接線方向（刃が進む向き）", color=RED, fontsize=15)
ax.text(-1.55, -1.5, "刃が削るのは φ = 0〜π（オレンジ）", fontsize=14, color=ORANGE)
ax.text(1.6, 0.05, "X", fontsize=15); ax.plot([0, 1.55], [0, 0], color=GRAY, ls="--"); ax.text(-0.1, 1.55, "Y", fontsize=15)
# ② 
ax = fig.add_axes([0.25, 0.5, 0.25, 0.46]); ax.set_xlim(-1.6, 1.8); ax.set_ylim(-1.6, 1.8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("② 食い違い Δq は、半径方向の成分だけ\n　 削る厚さを変える", loc="left", fontsize=19)
ax.add_patch(Circle((0, 0), 1, fill=False, lw=2, ls="--", ec=GRAY)); dq = np.array([0.55, 0.15])
ax.add_patch(Circle(dq, 1, fill=False, lw=2)); arr(ax, (0, 0), dq, GREEN); ax.text(0.15, -0.3, r"$\Delta\boldsymbol{q}=(\Delta x,\Delta y)$", color=GREEN, fontsize=16)
p = er; dh = er @ dq; arr(ax, p, p + dh * er, BLUE); ax.text(*(p + dh * er + [0.05, 0.05]), r"$\Delta h=\boldsymbol{e}_r\cdot\Delta\boldsymbol{q}$" + "\n" + r"$=\sin\phi\,\Delta x+\cos\phi\,\Delta y$", color=BLUE, fontsize=16)
ax.text(-1.55, -1.5, "点線：揺れる前の刃先の円　実線：揺れた後", fontsize=14, color=GRAY)
# ③
ax = fig.add_axes([0.5, 0.5, 0.25, 0.46]); ax.set_xlim(-1.6, 1.8); ax.set_ylim(-1.6, 1.8); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("③ 厚さが変わると、2 方向の力が変わる", loc="left", fontsize=19)
ax.add_patch(Circle((0, 0), 1, fill=False, lw=2)); ax.plot(*er, "o", color="k", ms=9)
arr(ax, er, er + 0.9 * et, RED); ax.text(*(er + 0.9 * et + [-0.6, -0.35]), r"$\Delta F_t=K_t\,a_p\,\Delta h$", color=RED, fontsize=16)
arr(ax, er, er + 0.9 * Kr * er, BLUE); ax.text(*(er + 0.3 * er + [0.1, 0.1]), r"$\Delta F_r=K_r\,\Delta F_t$", color=BLUE, fontsize=16)
F = 0.9 * (et + Kr * er); arr(ax, er, er + F, "k", 2); ax.text(*(er + F + [0.05, -0.25]), "合力", fontsize=15)
ax.text(-1.55, -1.25, r"$\Delta\boldsymbol{F}=K_t a_p\left(\boldsymbol{e}_t+K_r\boldsymbol{e}_r\right)\,\boldsymbol{e}_r^{T}\Delta\boldsymbol{q}$", fontsize=17)
ax.text(-1.55, -1.6, r"$=K_t a_p\,B(\phi)\,\Delta\boldsymbol{q}$", fontsize=17)
# ③' B(φ) の形
ax = fig.add_axes([0.75, 0.5, 0.25, 0.46]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title("④ 書き出すと B(φ)", loc="left", fontsize=19)
ax.text(0.2, 7.6, r"$B(\phi)=\left(\boldsymbol{e}_t+K_r\boldsymbol{e}_r\right)\boldsymbol{e}_r^{T}$", fontsize=19)
rows = [r"$B_{xx}=\sin\phi\cos\phi+K_r\sin^2\phi$", r"$B_{xy}=\cos^2\phi+K_r\sin\phi\cos\phi$",
        r"$B_{yx}=-\sin^2\phi+K_r\sin\phi\cos\phi$", r"$B_{yy}=-\sin\phi\cos\phi+K_r\cos^2\phi$"]
for i, r in enumerate(rows): ax.text(0.4, 5.9 - i * 1.25, r, fontsize=17)
ax.text(0.2, 0.4, "刃が回ると φ が変わるので、B も時々刻々変わる", fontsize=14, color=GRAY)
# ⑤ 平均
f = np.linspace(0, np.pi, 400); s, c = np.sin(f), np.cos(f)
comp = {"$B_{xx}$": s * c + Kr * s**2, "$B_{xy}$": c**2 + Kr * s * c, "$B_{yx}$": -s**2 + Kr * s * c, "$B_{yy}$": -s * c + Kr * c**2}
cols = [BLUE, ORANGE, GREEN, RED]
for i, ((k, v), col) in enumerate(zip(comp.items(), cols)):
    ax = fig.add_axes([0.04 + i * 0.235, 0.07, 0.2, 0.33])
    ax.plot(np.degrees(f), v, color=col, lw=2.5); ax.axhline(v.mean(), color=col, ls="--", lw=2)
    ax.fill_between(np.degrees(f), v, v.mean(), color=col, alpha=0.12); ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_xlim(0, 180); ax.set_ylim(-1.15, 1.15); ax.set_xticks([0, 90, 180]); ax.set_xlabel(r"刃の角度 $\phi$ [°]", fontsize=13); ax.grid(alpha=0.3)
    ax.set_title(f"{k}：半周の平均 = {v.mean():+.2f}", loc="left", fontsize=16, color=col)
fig.text(0.04, 0.44, r"⑤ φ = 0〜π で平均（$\sin\phi\cos\phi$ は平均 0、$\sin^2\phi$・$\cos^2\phi$ は平均 1/2）→ 刃が削るのは 1 回転の半分、刃は Z 枚なので × Z/2", fontsize=17)
fig.text(0.04, 0.005, "⑥ " + r"$B_0=\frac{Z}{2}\times$" + "（上の破線の値）" + r"$=\frac{Z}{4}\,[[K_r,\ 1],\ [-1,\ K_r]]$" + "　4 刃なら " + r"$B_{xx}=B_{yy}=0.3$、$B_{xy}=1$、$B_{yx}=-1$", fontsize=18)
fig.savefig(D / "results" / "B_derivation.png"); print("ok")
