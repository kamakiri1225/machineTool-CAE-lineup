"""ref004 の解説図（模式図・減衰の図）: python make_ref004_figs2.py"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon, Circle

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref004")
os.makedirs(OUT, exist_ok=True)


def spring(ax, p0, p1, n=6, w=0.12, color="k", lw=2.0, lead=0.18):
    """p0→p1 のジグザグばね"""
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0; L = np.linalg.norm(d); t = d / L; nrm = np.array([-t[1], t[0]])
    s = [0, lead]
    zz = np.linspace(lead, 1 - lead, 2 * n + 1)
    pts = [p0, p0 + t * L * lead]
    for i, z in enumerate(zz[1:-1]):
        pts.append(p0 + t * L * z + nrm * w * (1 if i % 2 == 0 else -1))
    pts += [p0 + t * L * (1 - lead), p1]
    pts = np.array(pts)
    ax.plot(pts[:, 0], pts[:, 1], color=color, lw=lw, solid_joinstyle="round")


def box(ax, x, y, w, h, fc, ec="k", text=None, fs=11, tc="k", **kw):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=1.5, **kw))
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)


def arrow(ax, p0, p1, color="k", lw=2, ms=18, style="-|>"):
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=ms))


def ground(ax, x0, x1, y):
    ax.plot([x0, x1], [y, y], "k", lw=2)
    for x in np.arange(x0, x1, 0.15):
        ax.plot([x, x - 0.1], [y, y - 0.1], "k", lw=0.8)


# ===================================================================== 1. ばねの配置
fig, axs = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw=dict(width_ratios=[1.35, 1]))
# --- 側面図（X-Z）
ax = axs[0]; ax.set_xlim(-0.3, 7.3); ax.set_ylim(-0.6, 4.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("側面から見る（送り方向 X が左右）", loc="left", fontsize=13)
box(ax, 0, 0, 7, 0.9, "#d5d8dc", text="ベッド")
box(ax, 0.3, 0.9, 6.4, 0.25, "#aab7b8")          # レール
ax.text(6.75, 1.02, "レール", fontsize=10, va="center")
box(ax, 1.2, 2.3, 4.0, 0.8, "#aed6f1", text="テーブル", fs=13)
for cx in (1.5, 4.2):                               # キャリッジ
    spring(ax, (cx + 0.35, 1.15), (cx + 0.35, 2.3), color=BLUE, w=0.13)
    ax.add_patch(Rectangle((cx, 1.15), 0.7, 0.2, fc=BLUE, alpha=0.35, ec=BLUE))
ax.text(4.75, 2.0, "ガイド 上下 $k_z$", color=BLUE, fontsize=11, va="center", ha="left",
        bbox=dict(fc="white", ec="none", alpha=0.85))
# ボールねじ
ax.add_patch(Rectangle((0.2, 1.55), 6.5, 0.18, fc="#f5b7b1", ec=RED))
box(ax, -0.25, 1.35, 0.45, 0.6, "#555", text="M", tc="white")
ax.text(-0.02, 1.15, "モータ", ha="center", va="top", fontsize=9)
spring(ax, (3.0, 1.64), (3.0, 2.3), color=RED, w=0.0, lead=0.5)   # ナットブラケット（剛）
box(ax, 2.65, 1.45, 0.7, 0.38, RED, alpha=0.5)
spring(ax, (5.6, 3.6), (6.9, 3.6), color=RED, w=0.13)
ax.text(5.5, 4.0, "ボールねじ $k_x$（送り方向だけ）", color=RED, fontsize=11)
ax.text(3.42, 1.3, "ナット", color=RED, fontsize=10, ha="left")
arrow(ax, (3.3, 3.6), (5.4, 3.6), color="k")
ax.text(3.3, 3.85, "テーブルを X に押すと…", fontsize=11)
ax.text(1.2, -0.45, "ガイドは X 方向に転がるだけ → X は 0。X を支えるのはボールねじだけ", fontsize=12,
        color=RED, weight="bold")
# --- 正面図（Y-Z）
ax = axs[1]; ax.set_xlim(-0.3, 5.3); ax.set_ylim(-0.6, 4.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("正面から見る（横方向 Y が左右）", loc="left", fontsize=13)
box(ax, 0, 0, 5, 0.9, "#d5d8dc", text="ベッド")
box(ax, 0.6, 2.3, 3.8, 0.8, "#aed6f1", text="テーブル", fs=13)
for rx in (0.9, 3.5):
    box(ax, rx, 0.9, 0.6, 0.35, "#aab7b8")
    spring(ax, (rx + 0.3, 1.25), (rx + 0.3, 2.3), color=BLUE, w=0.13)
    spring(ax, (rx - 0.55 if rx < 2 else rx + 1.15, 1.75), (rx + 0.3, 1.75), color=GREEN, w=0.1, n=4)
ax.text(0.15, 2.0, "横 $k_y$", color=GREEN, fontsize=11)
ax.text(1.45, 2.05, "上下 $k_z$", color=BLUE, fontsize=11)
ax.add_patch(Circle((2.5, 1.65), 0.2, fc="#f5b7b1", ec=RED))
ax.text(2.5, 1.2, "ねじ軸\n（紙面奥へ $k_x$）", color=RED, fontsize=10, ha="center", va="top")
ax.text(0.0, -0.45, "1 キャリッジに $k_z$（上下）と $k_y$（横）の 2 本", fontsize=12, color=BLUE, weight="bold")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "spring_map.png")); plt.close(fig)

# ===================================================================== 2. ガイドの荷重-変位と取付部
fig, axs = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw=dict(width_ratios=[1.1, 1]))
ax = axs[0]
Fp = 9.28; K = 2.0                                 # kN, kN/µm (RA45 中予圧)
F = np.linspace(0, 45, 400); Fk = 2.8 * Fp
d = np.where(F < Fk, F / K, Fk / K + (F - Fk) / (0.55 * K))
ax.plot(d, F, color=BLUE, lw=2.8)
ax.axvspan(0, Fk / K, color=BLUE, alpha=0.08)
ax.text(1.2, 38, "予圧が効いている範囲\n（ここの傾き = 剛性）", color=BLUE, fontsize=12)
ax.plot([0, 5], [0, 10], "--", color=RED, lw=1.5)
ax.plot([5, 5], [0, 10], ":", color=RED); ax.plot([0, 5], [10, 10], ":", color=RED)
ax.annotate("10,000 N で 5 µm\n→ $k = 10000/5 = 2000$ N/µm", xy=(5, 10), xytext=(7, 4), fontsize=12, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.annotate("予圧の約 2.8 倍で予圧が抜ける\n→ 傾きが寝る（柔らかくなる）", xy=(Fk / K, Fk), xytext=(15.5, 22), fontsize=11,
            arrowprops=dict(arrowstyle="->"))
ax.set_xlabel("変位 [µm]"); ax.set_ylabel("荷重 [kN]"); ax.set_xlim(0, 30); ax.set_ylim(0, 45); ax.grid(alpha=0.3)
ax.set_title("カタログの荷重−変位線図の読み方（RA45・中予圧のイメージ）", loc="left", fontsize=12)

ax = axs[1]; ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
ax.set_title("実機では取付部のばねが直列に入る → 実効剛性は下がる", loc="left", fontsize=12)
ground(ax, 0.3, 2.3, 0.6); ax.text(1.3, 0.25, "ベッド", ha="center", fontsize=10)
spring(ax, (1.3, 0.6), (1.3, 2.0), color=GRAY, w=0.18); ax.text(1.6, 1.3, "取付面・ボルト\nベッド局部 $k_{\\mathrm{mount}}$", fontsize=10, va="center")
spring(ax, (1.3, 2.0), (1.3, 3.6), color=BLUE, w=0.18); ax.text(1.6, 2.8, "ガイド（カタログ）$k_{\\mathrm{cat}}$", fontsize=10, color=BLUE, va="center")
spring(ax, (1.3, 3.6), (1.3, 4.6), color=GRAY, w=0.18); ax.text(1.6, 4.1, "テーブル局部", fontsize=10, va="center")
box(ax, 0.5, 4.6, 1.6, 0.6, "#aed6f1", text="テーブル", fs=10)
ax.text(5.2, 5.3, r"$\dfrac{1}{k_{\mathrm{eff}}} = \dfrac{1}{k_{\mathrm{cat}}} + \dfrac{1}{k_{\mathrm{mount}}} + \cdots$", fontsize=16)
vals = [("カタログ値", 2000, BLUE), ("取付部", 2000, GRAY), ("実効値", 1000, RED)]
for i, (nm, v, c) in enumerate(vals):
    y = 3.6 - i * 1.0
    ax.add_patch(Rectangle((5.2, y), v / 2000 * 4.2, 0.55, fc=c, alpha=0.75))
    ax.text(5.1, y + 0.27, nm, ha="right", va="center", fontsize=11)
    ax.text(5.3 + v / 2000 * 4.2, y + 0.27, f"{v:,} N/µm", va="center", fontsize=11)
ax.text(5.2, 0.6, "同じくらいの硬さが直列だと半分になる\n→ まずカタログ値 × 0.5 から", fontsize=12, color=RED)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "guide_catalog_mount.png")); plt.close(fig)

# ===================================================================== 3. ボールねじの直列ばねと支持方式
fig, axs = plt.subplots(2, 1, figsize=(14, 7.6))
for ax, kind in zip(axs, ["fs", "ff"]):
    ax.set_xlim(-0.5, 13.5); ax.set_ylim(-1.3, 2.4); ax.set_aspect("equal"); ax.axis("off")
    box(ax, -0.4, 0.2, 0.8, 0.8, "#555", text="M", tc="white")
    # 固定側軸受
    box(ax, 0.8, 0.0, 0.5, 1.2, ORANGE, alpha=0.6); ax.text(1.05, -0.35, "固定側\n軸受", ha="center", va="top", fontsize=9)
    ground(ax, 0.6, 1.5, 0.0)
    # 軸
    ax.add_patch(Rectangle((1.3, 0.5), 11.0, 0.2, fc="#f5b7b1", ec=RED))
    xn = 6.0
    box(ax, xn, 0.25, 0.9, 0.7, RED, alpha=0.55, text="ナット", fs=9)
    box(ax, xn - 0.6, 1.4, 2.1, 0.7, "#aed6f1", text="テーブル", fs=10)
    ax.plot([xn + 0.45, xn + 0.45], [0.95, 1.4], "k", lw=2)
    ax.annotate("", xy=(xn, -0.25), xytext=(1.3, -0.25), arrowprops=dict(arrowstyle="<->"))
    ax.text((xn + 1.3) / 2, -0.55, "$a$（ここだけ伸び縮みする）", ha="center", va="top", fontsize=11, color=RED)
    if kind == "fs":
        box(ax, 12.3, 0.1, 0.5, 1.0, "#d5d8dc"); ax.text(12.55, -0.35, "支持側\n（軸方向フリー）", ha="center", va="top", fontsize=9)
        ax.set_title("固定 − 支持：力は固定側だけで受ける  →  $K_{\\mathrm{shaft}} = AE/a$（モータから遠いほど柔らかい）",
                     loc="left", fontsize=12)
        # 直列の模式
        y0 = 1.9
        ax.text(8.6, 2.15, "直列：", fontsize=11)
        spring(ax, (9.4, y0 + 0.25), (10.5, y0 + 0.25), color=ORANGE, w=0.1, n=4)
        spring(ax, (10.5, y0 + 0.25), (11.8, y0 + 0.25), color=RED, w=0.1, n=4)
        spring(ax, (11.8, y0 + 0.25), (12.9, y0 + 0.25), color=BLUE, w=0.1, n=4)
        ax.text(9.95, y0 - 0.05, "軸受", color=ORANGE, ha="center", fontsize=9)
        ax.text(11.15, y0 - 0.05, "ねじ軸", color=RED, ha="center", fontsize=9)
        ax.text(12.35, y0 - 0.05, "ナット", color=BLUE, ha="center", fontsize=9)
    else:
        box(ax, 12.3, 0.0, 0.5, 1.2, ORANGE, alpha=0.6); ground(ax, 12.1, 13.0, 0.0)
        ax.text(12.55, -0.35, "固定側\n軸受", ha="center", va="top", fontsize=9)
        ax.annotate("", xy=(12.3, -0.25), xytext=(xn + 0.9, -0.25), arrowprops=dict(arrowstyle="<->"))
        ax.text((xn + 0.9 + 12.3) / 2, -0.55, "$L-a$", ha="center", va="top", fontsize=11, color=BLUE)
        ax.set_title("固定 − 固定：両側から支える（並列） →  $K_{\\mathrm{shaft}} = AE/a + AE/(L-a)$（中央が最小）",
                     loc="left", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "ballscrew_support.png")); plt.close(fig)

# ===================================================================== 4. テーブル位置ごとに解く
E = 206e3; A = np.pi / 4 * 34.0**2


def kfs(a):
    return 1 / (a / (A * E) + 2 / 1000e3) / 1e3


fig, axs = plt.subplots(3, 1, figsize=(12, 7.2))
for ax, a, lab in zip(axs, [150, 500, 850], ["モータ側", "中央", "反対側"]):
    ax.set_xlim(-0.6, 11); ax.set_ylim(-0.3, 1.6); ax.set_aspect("equal"); ax.axis("off")
    box(ax, -0.5, 0.15, 0.5, 0.5, "#555", text="M", tc="white", fs=9)
    ax.add_patch(Rectangle((0, 0.32), 10, 0.14, fc="#f5b7b1", ec=RED))
    xt = a / 100
    box(ax, xt - 1.2, 0.65, 2.4, 0.6, "#aed6f1", text="テーブル", fs=10)
    box(ax, xt - 0.25, 0.2, 0.5, 0.38, RED, alpha=0.55)
    ax.plot([0, xt], [0.08, 0.08], color=RED, lw=4, alpha=0.7)
    k = kfs(a)
    ax.text(10.3, 0.75, f"{lab}\n$a$ = {a} mm\n$K$ ≈ {k:.0f} N/µm", fontsize=11, va="center",
            color=RED if a == 850 else "k")
axs[0].set_title("テーブル位置ごとにばね値（と部品位置）を変えて解く（固定−支持、φ40 の例）", loc="left", fontsize=12)
axs[2].text(0, -0.25, "赤い太線 = 伸び縮みするねじ軸の長さ $a$。長いほど柔らかい → 反対側が最悪値", fontsize=11, color=RED)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "position_cases.png")); plt.close(fig)

# ===================================================================== 5. 減衰比とピーク
f = np.linspace(50, 350, 3000); fn = 181.0; r = f / fn
fig, axs = plt.subplots(1, 2, figsize=(15, 5.2))
ax = axs[0]
for z, c in [(0.01, RED), (0.02, BLUE), (0.05, GREEN)]:
    H = 1 / np.sqrt((1 - r**2) ** 2 + (2 * z * r) ** 2)
    ax.semilogy(f, H, color=c, lw=2.3, label=f"$\\zeta$ = {z*100:.0f} %  → ピーク ≈ $1/(2\\zeta)$ = {1/(2*z):.0f} 倍")
ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("コンプライアンス（静的変位を 1 とした倍率）")
ax.set_title("減衰比 $\\zeta$ はピークの高さを決める（固有振動数は変わらない）", loc="left", fontsize=12)
ax.legend(fontsize=11); ax.grid(alpha=0.3, which="both")
ax = axs[1]
z = 0.03; H = 1 / np.sqrt((1 - r**2) ** 2 + (2 * z * r) ** 2); Hm = H.max()
ax.plot(f, H, color=BLUE, lw=2.5)
lvl = Hm / np.sqrt(2); idx = np.where(H >= lvl)[0]; f1, f2 = f[idx[0]], f[idx[-1]]
ax.axhline(lvl, color=GRAY, ls="--"); ax.text(122, lvl * 1.03, "ピーク × $1/\\sqrt{2}$", color=GRAY, fontsize=11)
ax.annotate("", xy=(f2, lvl), xytext=(f1, lvl), arrowprops=dict(arrowstyle="<->", color=RED, lw=2))
ax.text(fn, lvl * 0.88, f"$f_1$ = {f1:.0f}, $f_2$ = {f2:.0f} Hz", ha="center", color=RED, fontsize=11, va="top")
ax.text(205, Hm * 0.62, r"$\zeta \approx \dfrac{f_2 - f_1}{2 f_n}$" + f"\n= ({f2:.0f} − {f1:.0f}) / (2×{fn:.0f})\n≈ {(f2-f1)/2/fn*100:.1f} %",
        fontsize=14, color=RED)
ax.set_xlim(120, 260); ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("コンプライアンス倍率")
ax.set_title("打撃試験の FRF から $\\zeta$ を読む（半値幅法）", loc="left", fontsize=12); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "damping_peak.png")); plt.close(fig)

# ===================================================================== 6. 減衰の出どころ
src = [("鋳物・鋼（材料）", 0.0003, 0.002, GRAY), ("ボールねじ・軸受", 0.01, 0.04, ORANGE),
       ("リニアガイド（ころがり）", 0.02, 0.06, BLUE), ("ボルト結合部・接合面", 0.02, 0.08, GREEN),
       ("すべり案内（角ガイド）", 0.05, 0.15, RED)]
fig, ax = plt.subplots(figsize=(11, 4.2))
for i, (nm, lo, hi, c) in enumerate(src):
    ax.barh(i, hi - lo, left=lo, color=c, alpha=0.75, height=0.55)
    ax.text(hi * 1.1, i, f"{lo*100:g} 〜 {hi*100:g} %", va="center", fontsize=10)
ax.set_yticks(range(len(src))); ax.set_yticklabels([s[0] for s in src], fontsize=11)
ax.set_xscale("log"); ax.set_xlim(1e-4, 0.5); ax.set_xlabel("その部分が主役のときの減衰比 $\\zeta$ のおおよその幅（対数）")
ax.axvspan(0.01, 0.05, color=BLUE, alpha=0.07); ax.text(0.011, -0.42, "機械全体の構造モード 1〜5 %", color=BLUE, fontsize=10)
ax.set_title("減衰はほとんど「材料」ではなく「つなぎ目」から生まれる（目安・イメージ）", loc="left", fontsize=12)
ax.grid(alpha=0.3, axis="x", which="both"); ax.set_ylim(4.5, -0.7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "damping_sources.png")); plt.close(fig)

# ===================================================================== 7. 値の決め方の流れ
steps = [("① カタログ値を読む", "ガイド: 上下剛性\nねじ: ナット・軸受", BLUE),
         ("② 取付部を見込む", "ガイド × 0.5\nねじ軸は $AE/a$ で直列", ORANGE),
         ("③ FEM で解く", "位置 3 か所\n減衰 $\\zeta$ = 2 %", GREEN),
         ("④ 打撃試験", "固有振動数\nピーク高さ", GRAY),
         ("⑤ 合わせ込み", "振動数のずれ → $k$\nピークのずれ → $\\zeta$", RED)]
fig, ax = plt.subplots(figsize=(15, 3.6)); ax.set_xlim(0, 15); ax.set_ylim(-0.6, 3); ax.axis("off")
for i, (t, s, c) in enumerate(steps):
    x = 0.2 + i * 3.0
    ax.add_patch(FancyBboxPatch((x, 0.3), 2.5, 2.2, boxstyle="round,pad=0.05", fc=c, alpha=0.15, ec=c, lw=2))
    ax.text(x + 1.25, 2.1, t, ha="center", fontsize=12, weight="bold", color=c)
    ax.text(x + 1.25, 1.15, s, ha="center", va="center", fontsize=11)
    if i < 4:
        arrow(ax, (x + 2.58, 1.4), (x + 2.95, 1.4))
ax.annotate("", xy=(7.45, 0.25), xytext=(13.45, 0.25),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.8, connectionstyle="arc3,rad=-0.12", mutation_scale=16))
ax.text(10.4, -0.55, "合うまで繰り返す（固有振動数は $\\sqrt{k}$ に比例）", ha="center", color=RED, fontsize=11)
fig.savefig(os.path.join(OUT, "workflow.png")); plt.close(fig)
print("done")
