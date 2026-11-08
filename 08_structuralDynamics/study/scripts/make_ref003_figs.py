"""ref003 の図を生成する: python make_ref003_figs.py
2慣性系（モータ ― ボールねじ(ばね) ― テーブル）の位置・速度ループで、
機械の1次固有振動数とサーボゲインの上限の関係を示す。
"""
import os
import numpy as np
from figstyle import plt, BLUE, RED, ORANGE, GRAY, GREEN
from matplotlib.patches import Rectangle, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "..", "images", "ref003")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- モデル
Jm, JL = 2e-3, 4e-3            # モータ慣性, テーブル換算慣性 [kg m^2]
ZETA = 0.02                    # 機械の減衰比
TD = 0.5e-3                    # 制御の遅れ（電流ループ・サンプリング等をまとめたむだ時間）[s]


def stiffness_for(f_res):
    """共振周波数 f_res [Hz] になるねじり剛性 K"""
    jeq = Jm * JL / (Jm + JL)
    return (2 * np.pi * f_res) ** 2 * jeq


def damping_for(K):
    jeq = Jm * JL / (Jm + JL)
    return 2 * ZETA * np.sqrt(K * jeq)


def plant(w, K):
    """トルク → モータ速度 の周波数応答"""
    s = 1j * w
    c = damping_for(K)
    return (JL * s**2 + c * s + K) / (s * (Jm * JL * s**2 + c * (Jm + JL) * s + K * (Jm + JL)))


KV = 1.0                       # 速度ループ比例ゲイン [N m s/rad]（モータ側で検出）


def _inner(w, K):
    """速度指令 → テーブル速度 の周波数応答（速度ループはモータ側で閉じる）"""
    s = 1j * w
    c = damping_for(K)
    wi = KV / (Jm + JL) / 5
    C = KV * (1 + wi / s) * np.exp(-s * TD)
    Gm = plant(w, K)
    Tv = C * Gm / (1 + C * Gm)
    return Tv * (c * s + K) / (JL * s**2 + c * s + K)


def open_loop(w, K, Kp):
    """位置ループ（テーブル位置で閉じる）の一巡伝達関数"""
    return Kp / (1j * w) * _inner(w, K)


def gain_margin(K, Kp):
    w = 2 * np.pi * np.logspace(-0.5, 3.5, 200000)
    L = open_loop(w, K, Kp)
    ph = np.unwrap(np.angle(L)); mag = np.abs(L)
    gms = []
    for n in range(4):
        idx = np.where(np.diff(np.sign(ph + np.pi + 2 * np.pi * n)))[0]
        gms += [1 / mag[i] for i in idx]
    return min(gms)


def simulate(K, Kp, T=0.3, dt=5e-6, xref=1.0):
    c = damping_for(K)
    wi = KV / (Jm + JL) / 5
    n = int(T / dt); nd = int(TD / dt)
    buf = np.zeros(nd); k = 0
    wm = wl = tw = xl = ie = 0.0
    t = np.arange(n) * dt
    out = np.zeros(n)
    for i in range(n):
        wref = Kp * (xref - xl)               # 位置ループ（P制御）
        e = wref - wm                          # 速度ループ（PI制御、モータ側）
        ie += e * dt
        tau = buf[k]; buf[k] = KV * (e + wi * ie); k = (k + 1) % nd
        tau = min(max(tau, -50.0), 50.0)
        f_spring = K * tw + c * (wm - wl)
        wm += (tau - f_spring) / Jm * dt
        wl += f_spring / JL * dt
        tw += (wm - wl) * dt
        xl += wl * dt
        out[i] = xl
    return t, out


def arrow(ax, p, q, color=BLUE, style="<|-|>", lw=1.5, ms=12):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, color=color, lw=lw, mutation_scale=ms))


# ---------------------------------------------------------------- 図1 固有振動数とは
def natural_frequency():
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw=dict(width_ratios=[0.8, 1, 1.2]))

    # (a) ばね-質量
    ax = axs[0]; ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(0, 6); ax.set_ylim(0, 10)
    ax.set_title("(a) いちばん単純なモデル", loc="left", fontsize=13)
    ax.add_patch(Rectangle((1, 9), 3, 0.3, fc="#566573"))
    ys = np.linspace(9, 5.5, 13); xs = 2.5 + 0.4 * np.array([0, 1, -1] * 4 + [0])
    ax.plot(xs, ys, color=BLUE, lw=2)
    ax.text(3.1, 7.2, "ばね $k$\n（剛性）", color=BLUE, fontsize=12)
    ax.add_patch(Rectangle((1.5, 3.8), 2, 1.7, fc="#aab7b8", ec="#566573"))
    ax.text(2.5, 4.65, "質量 $m$", ha="center", va="center", fontsize=12)
    arrow(ax, (4.2, 3.6), (4.2, 5.7), color=RED)
    ax.text(4.4, 4.5, "ゆれる", color=RED, fontsize=11)
    ax.text(0.2, 2.6, "$f_n = \\dfrac{1}{2\\pi}\\sqrt{\\dfrac{k}{m}}$  [Hz]", fontsize=17)
    ax.text(0.2, 0.6, "かたい（$k$ 大）・軽い（$m$ 小）ほど\n速くゆれる ＝ 固有振動数が高い", fontsize=11)

    # (b) 片持ちはりのモード
    ax = axs[1]
    ax.set_title("(b) 構造には複数の固有振動数がある（片持ちはりの例）", loc="left", fontsize=13)
    bl = [1.8751, 4.6941, 7.8548]
    x = np.linspace(0, 1, 300)
    cols = [RED, ORANGE, GREEN]
    for i, b in enumerate(bl):
        sg = (np.cosh(b) + np.cos(b)) / (np.sinh(b) + np.sin(b))
        phi = np.cosh(b * x) - np.cos(b * x) - sg * (np.sinh(b * x) - np.sin(b * x))
        phi = phi / np.max(np.abs(phi))
        off = -2.6 * i
        ax.plot(x, off + 0 * x, color=GRAY, lw=0.8, ls="--")
        ax.plot(x, off + 0.9 * phi, color=cols[i], lw=2.5)
        ax.text(1.04, off, f"{i+1}次モード\n$f_{i+1} = {(b/bl[0])**2:.2f}\\,f_1$", color=cols[i], va="center", fontsize=11)
        ax.add_patch(Rectangle((-0.06, off - 0.6), 0.06, 1.2, fc="#566573"))
    ax.set_xlim(-0.1, 1.45); ax.set_ylim(-6.3, 1.9); ax.axis("off")
    ax.text(0, 1.45, "根元を固定。いちばん低い＝ゆれやすいのが「1次」", fontsize=11)

    # (c) 周波数応答
    ax = axs[2]
    ax.set_title("(c) 揺する周波数を変えたときのゆれ（周波数応答）", loc="left", fontsize=13)
    r = np.linspace(0.01, 3, 2000)
    for z, ls in [(0.02, "-"), (0.1, "--")]:
        H = 1 / np.sqrt((1 - r**2) ** 2 + (2 * z * r) ** 2)
        ax.semilogy(r, H, ls, color=BLUE, lw=2, label=f"減衰比 $\\zeta = {z}$")
    ax.axvline(1, color=RED, lw=1, ls=":")
    ax.text(1.12, 8, "固有振動数で\nゆれが何十倍にも\nなる（共振）", color=RED, fontsize=11)
    ax.set_xlabel("揺する周波数 / 固有振動数　$f / f_n$")
    ax.set_ylabel("ゆれの大きさ（静的たわみに対する倍率）")
    ax.grid(alpha=0.3, which="both"); ax.legend(loc="upper right"); ax.set_ylim(0.05, 50)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "natural_frequency.png")); plt.close(fig)


# ---------------------------------------------------------------- 図2 ゲインの上限
def servo_gain_limit():
    cases = [
        ("A: 1次 30 Hz・位置ループゲイン 20 /s", 30, 20, BLUE, "-"),
        ("B: 1次 30 Hz・位置ループゲイン 100 /s", 30, 100, RED, "-"),
        ("C: 1次 120 Hz・位置ループゲイン 100 /s", 120, 100, GREEN, "--"),
    ]
    for lab, fr, Kp, *_ in cases:
        print(f"{lab}: gain margin = {20*np.log10(gain_margin(stiffness_for(fr), Kp)):+.1f} dB")

    fig = plt.figure(figsize=(16, 7.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], height_ratios=[1.3, 1])
    axm = fig.add_subplot(gs[0, 0]); axp = fig.add_subplot(gs[1, 0], sharex=axm)
    axs = fig.add_subplot(gs[:, 1])

    w = 2 * np.pi * np.logspace(0, 2.7, 4000)
    for lab, fr, Kp, col, ls in cases:
        L = open_loop(w, stiffness_for(fr), Kp)
        axm.semilogx(w / 2 / np.pi, 20 * np.log10(abs(L)), ls, color=col, lw=2, label=lab)
        axp.semilogx(w / 2 / np.pi, np.degrees(np.unwrap(np.angle(L))), ls, color=col, lw=2)
    axm.axhline(0, color="k", lw=1)
    axm.set_ylabel("ゲイン [dB]"); axm.set_ylim(-50, 30); axm.grid(alpha=0.3, which="both")
    axm.legend(loc="lower left", fontsize=11)
    axm.set_title("位置ループの一巡伝達関数（ボード線図）", loc="left", fontsize=13)
    axm.annotate("B は 0 dB を横切る所（約 22 Hz）で\n位相がすでに −180° より遅れている\n→ 発振する",
                 xy=(22, 0), xytext=(45, 12), color=RED, fontsize=11,
                 arrowprops=dict(arrowstyle="->", color=RED))
    axp.axhline(-180, color="k", lw=1)
    axp.set_ylabel("位相 [deg]"); axp.set_xlabel("周波数 [Hz]"); axp.set_ylim(-540, 0)
    axp.grid(alpha=0.3, which="both")
    axp.text(1.1, -170, "−180°", fontsize=10)

    for lab, fr, Kp, col, ls in cases:
        t, x = simulate(stiffness_for(fr), Kp)
        axs.plot(t * 1000, x, ls, color=col, lw=2, label=lab)
    axs.axhline(1, color=GRAY, lw=1, ls=":")
    axs.set_ylim(-0.2, 2.0); axs.set_xlim(0, 300)
    axs.set_xlabel("時間 [ms]"); axs.set_ylabel("テーブルの位置（指令 = 1）")
    axs.set_title("位置指令をステップで与えたときのテーブルの動き", loc="left", fontsize=13)
    axs.grid(alpha=0.3); axs.legend(loc="lower right", fontsize=11)
    axs.text(90, 1.9, "A: 安定だが目標に着くのが遅く、揺れも残る\nB: ゲインを上げたら発振（振動が成長）\nC: 1次固有振動数を上げれば\n　 同じ高ゲインでも速く安定",
             fontsize=11, va="top", bbox=dict(fc="white", ec=GRAY, alpha=0.95))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "servo_gain_limit.png")); plt.close(fig)


# ---------------------------------------------------------------- 図3 3重ループ
def _box(ax, x0, y0, w, h, text, fc="#eaf2f8", ec=BLUE, fs=11):
    ax.add_patch(Rectangle((x0, y0), w, h, fc=fc, ec=ec, lw=1.5, zorder=3))
    ax.text(x0 + w / 2, y0 + h / 2, text, ha="center", va="center", fontsize=fs, zorder=4)


def _sum(ax, x, y, r=0.22):
    from matplotlib.patches import Circle
    ax.add_patch(Circle((x, y), r, fc="white", ec="k", lw=1.2, zorder=3))
    ax.text(x - 0.42, y + 0.18, "+", fontsize=10); ax.text(x + 0.05, y - 0.55, "−", fontsize=12)


def _arr(ax, pts, color="k", lw=1.4):
    for a, b in zip(pts[:-2], pts[1:-1]):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=lw, zorder=2)
    arrow(ax, pts[-2], pts[-1], color=color, style="-|>", lw=lw, ms=12)


def cascade_loop():
    fig, ax = plt.subplots(figsize=(16, 9.5))
    ax.set_xlim(0, 17); ax.set_ylim(-0.3, 10); ax.axis("off"); ax.set_aspect("equal")
    y = 7.2
    # ループの枠
    for (x0, y0, x1, y1, col, lab) in [
        (1.5, 3.6, 16.7, 9.3, RED, "位置ループ（外側）：位置のずれ → 速度指令"),
        (5.2, 4.4, 13.4, 8.7, GREEN, "速度ループ：速度のずれ → 電流（トルク）指令"),
        (8.7, 5.2, 11.75, 8.1, ORANGE, "電流ループ（内側）"),
    ]:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=col, ls="--", lw=1.8))
        ax.text(x0 + 0.1, y1 + 0.08, lab, color=col, fontsize=11, va="bottom")
    # ブロック
    ax.text(0.0, y + 0.35, "位置指令\n$x_{\\mathrm{ref}}$", fontsize=11, va="bottom")
    _sum(ax, 2.2, y)
    _box(ax, 2.8, y - 0.55, 1.7, 1.1, "位置制御\n$K_p$", fc="#fdedec", ec=RED)
    _sum(ax, 5.9, y)
    _box(ax, 6.5, y - 0.55, 1.7, 1.1, "速度制御\nPI", fc="#e9f7ef", ec=GREEN)
    _sum(ax, 9.3, y)
    _box(ax, 9.8, y - 0.55, 1.7, 1.1, "電流制御\n＋アンプ", fc="#fef5e7", ec=ORANGE)
    _box(ax, 12.0, y - 0.55, 1.2, 1.1, "モータ", fc="#f2f3f4", ec="#566573")
    _box(ax, 13.8, y - 0.55, 2.0, 1.1, "ボールねじ\n＋テーブル", fc="#f2f3f4", ec="#566573")
    # 前向きの矢印
    _arr(ax, [(0.9, y), (1.98, y)])
    _arr(ax, [(2.42, y), (2.8, y)])
    _arr(ax, [(4.5, y), (5.68, y)]); ax.text(4.65, y - 0.45, "速度指令", fontsize=10)
    _arr(ax, [(6.12, y), (6.5, y)])
    _arr(ax, [(8.2, y), (9.08, y)]); ax.text(8.25, y - 0.45, "電流指令", fontsize=10)
    _arr(ax, [(9.52, y), (9.8, y)])
    _arr(ax, [(11.5, y), (12.0, y)]); ax.text(11.45, y + 0.6, "電流", fontsize=10)
    _arr(ax, [(13.2, y), (13.8, y)]); ax.text(13.15, y + 0.6, "トルク", fontsize=10)
    _arr(ax, [(15.8, y), (16.9, y)]); ax.text(15.95, y + 0.15, "位置 $x$", fontsize=10)
    # フィードバック
    _arr(ax, [(11.6, y), (11.6, 5.6), (9.3, 5.6), (9.3, y - 0.22)], color=ORANGE)
    ax.text(9.5, 5.3, "電流センサ", color=ORANGE, fontsize=10)
    _arr(ax, [(12.6, y - 0.55), (12.6, 4.8), (5.9, 4.8), (5.9, y - 0.22)], color=GREEN)
    ax.text(6.1, 4.5, "モータの速度（モータ後ろのエンコーダ）", color=GREEN, fontsize=10)
    _arr(ax, [(16.3, y), (16.3, 4.0), (2.2, 4.0), (2.2, y - 0.22)], color=RED)
    ax.text(2.4, 3.7, "テーブルの位置（リニアスケール、またはモータのエンコーダから換算）", color=RED, fontsize=10)

    # 下段: 機械の絵
    yb = 1.3
    ax.add_patch(Rectangle((0.6, yb - 0.6), 0.5, 1.2, fc="#d6eaf8", ec=BLUE))          # エンコーダ
    ax.text(0.85, yb - 1.0, "エンコーダ", ha="center", fontsize=10, color=GREEN)
    ax.add_patch(Rectangle((1.1, yb - 0.8), 2.0, 1.6, fc="#aab7b8", ec="#566573"))      # モータ
    ax.text(2.1, yb, "サーボ\nモータ", ha="center", va="center", fontsize=11)
    ax.add_patch(Rectangle((3.1, yb - 0.3), 0.7, 0.6, fc="#f5cba7", ec="#a04000"))      # カップリング
    ax.text(3.45, yb + 0.5, "カップリング", ha="center", fontsize=9)
    ax.add_patch(Rectangle((3.8, yb - 0.5), 0.5, 1.0, fc="#566573"))                     # 軸受
    ax.add_patch(Rectangle((15.2, yb - 0.5), 0.5, 1.0, fc="#566573"))
    ax.text(4.05, yb - 0.95, "軸受", ha="center", fontsize=9); ax.text(15.45, yb - 0.95, "軸受", ha="center", fontsize=9)
    ax.add_patch(Rectangle((4.3, yb - 0.15), 10.9, 0.3, fc="#d5d8dc", ec="#566573"))     # ねじ軸
    for xx in np.arange(4.4, 15.1, 0.25):
        ax.plot([xx, xx + 0.15], [yb - 0.15, yb + 0.15], color="#566573", lw=0.7)
    ax.text(4.6, yb - 1.0, "ボールねじ（回転 → 直進に変える。ここが「ばね」になる）", fontsize=10)
    ax.add_patch(Rectangle((8.6, yb - 0.45), 1.2, 0.9, fc="#f5b041", ec="#b9770e"))      # ナット
    ax.add_patch(Rectangle((7.4, yb + 0.45), 3.6, 0.9, fc="#aed6f1", ec=BLUE))           # テーブル
    ax.text(9.2, yb + 0.9, "テーブル（工作物を載せる）", ha="center", va="center", fontsize=10)
    ax.text(10.0, yb + 0.15, "ナット", fontsize=9)
    ax.plot([5, 15], [yb + 1.55, yb + 1.55], color="#7d3c98", lw=3)                    # スケール
    ax.text(11.3, yb + 1.65, "リニアスケール（テーブル位置を直接測る）", color="#7d3c98", fontsize=10)
    arrow(ax, (11.4, yb + 0.9), (12.6, yb + 0.9), color=RED, style="-|>", lw=2.5, ms=16)
    ax.text(12.7, yb + 0.8, "送り", color=RED, fontsize=11)
    ax.text(0.0, 2.95, "実際の機械（横から見た図）", fontsize=12, fontweight="bold")
    fig.savefig(os.path.join(OUT, "cascade_loop.png")); plt.close(fig)


# ---------------------------------------------------------------- 図4 ゲインを上げるとは
def gain_meaning():
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1.05, 1, 1]))

    ax = axs[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_title("(a) 位置制御は「見えないばね」", loc="left", fontsize=13)
    for i, (k, lab, col) in enumerate([(3, "ゲイン低 ＝ やわらかいばね", BLUE), (9, "ゲイン高 ＝ かたいばね", RED)]):
        y0 = 7.2 - 4.0 * i
        ax.plot([1.0, 1.0], [y0 - 0.8, y0 + 0.8], color="k", lw=2)
        ax.text(0.2, y0 + 1.1, "目標位置", fontsize=10)
        n = k
        xs = np.linspace(1.0, 6.0, 2 * n + 1); ys = y0 + 0.35 * np.array([0] + [1, -1] * n)[:len(xs)]
        ys[-1] = y0
        ax.plot(xs, ys, color=col, lw=2)
        ax.add_patch(Rectangle((6.0, y0 - 0.6), 1.8, 1.2, fc="#aed6f1", ec=BLUE))
        ax.text(6.9, y0, "テーブル", ha="center", va="center", fontsize=10)
        ax.text(1.0, y0 - 1.5, lab, color=col, fontsize=11)
    ax.text(0.2, 0.2, "ずれ $e$ に比例して引き戻す：$v_{\\mathrm{cmd}} = K_p\\,e$\n$K_p$ を大きくする ＝ 同じずれでも強く・速く戻す", fontsize=11)

    ax = axs[1]
    ax.set_title("(b) 目標へ動かしたとき（理想的な機械）", loc="left", fontsize=13)
    t = np.linspace(0, 0.3, 600)
    for Kp, col in [(10, BLUE), (30, ORANGE), (100, RED)]:
        ax.plot(t * 1000, 1 - np.exp(-Kp * t), color=col, lw=2, label=f"$K_p = {Kp}\\ \\mathrm{{s^{{-1}}}}$（{1000/Kp:.0f} ms で 63 %）")
    ax.axhline(1, color=GRAY, ls=":")
    ax.set_xlabel("時間 [ms]"); ax.set_ylabel("位置（目標 = 1）"); ax.grid(alpha=0.3); ax.legend(loc="lower right")
    ax.text(110, 0.35, "時定数 $T = 1/K_p$\nゲインが高いほど速く着く", fontsize=11)

    ax = axs[2]
    ax.set_title("(c) 一定速度で動かしたとき（追従遅れ）", loc="left", fontsize=13)
    F = 1000 / 60  # mm/s
    t = np.linspace(0, 0.3, 600)
    ax.plot(t * 1000, F * t, color="k", lw=1.5, ls="--", label="指令（1000 mm/min）")
    for Kp, col in [(10, BLUE), (30, ORANGE), (100, RED)]:
        x = F * (t - (1 - np.exp(-Kp * t)) / Kp)
        ax.plot(t * 1000, x, color=col, lw=2, label=f"$K_p = {Kp}$：遅れ {F/Kp:.2f} mm")
    ax.set_xlabel("時間 [ms]"); ax.set_ylabel("位置 [mm]"); ax.grid(alpha=0.3); ax.legend(loc="upper left")
    ax.text(150, 0.6, "遅れ $e = F / K_p$\nゲインが高いほど\n指令にぴったり付いてくる", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "gain_meaning.png")); plt.close(fig)



# ---------------------------------------------------------------- 図5 発振の仕組み（4コマ）
def oscillation_mechanism():
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.22)

    # (1) ばねの先のテーブルは遅れてついてくる
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0, 0], hspace=0.55)
    fn, z = 30.0, 0.05
    for r, (f, lab) in enumerate([(10, "ゆっくり（10 Hz）：ほぼ同時に動く"), (30, "固有振動数（30 Hz）：1/4 周期遅れ、大きくゆれる"),
                                   (60, "速い（60 Hz）：半周期遅れ ＝ 逆向きに動く")]):
        ax = fig.add_subplot(sub[r])
        w, wn = 2 * np.pi * f, 2 * np.pi * fn
        H = (wn**2 + 2j * z * wn * w) / (wn**2 - w**2 + 2j * z * wn * w)
        t = np.linspace(0, 3 / f, 400)
        ax.plot(t * f, np.sin(w * t), color="k", lw=1.5, label="モータ（ねじの根元）")
        ax.plot(t * f, abs(H) * np.sin(w * t + np.angle(H)), color=RED, lw=2, label="テーブル（ばねの先）")
        ax.set_title(lab + f"（遅れ {-np.degrees(np.angle(H)):.0f}°）", loc="left", fontsize=11)
        ax.set_yticks([]); ax.set_xlim(0, 3); ax.grid(alpha=0.3)
        if r == 0:
            ax.legend(loc="upper right", fontsize=9, ncol=2)
        if r == 2:
            ax.set_xlabel("時間（周期の何倍か）")
    fig.text(0.07, 0.94, "① テーブルはボールねじ（ばね）の先にぶら下がっているので、モータより遅れてついてくる（固有振動数 30 Hz の例）", fontsize=12)

    # (2) ゲインを上げると反応する範囲が広がる
    ax = fig.add_subplot(gs[0, 1])
    w = 2 * np.pi * np.logspace(0, 2.5, 3000)
    K = stiffness_for(30)
    for Kp, col in [(20, BLUE), (100, RED)]:
        L = open_loop(w, K, Kp); mag = 20 * np.log10(abs(L)); fc = w[np.argmin(abs(mag))] / 2 / np.pi
        ax.semilogx(w / 2 / np.pi, mag, color=col, lw=2.2, label=f"位置ループゲイン {Kp} /s")
        ax.plot(fc, 0, "o", color=col, ms=9)
        ax.annotate(f"0 dB を横切る所 {fc:.1f} Hz", xy=(fc, 0), xytext=(fc * 1.3, 12 if Kp == 100 else -18),
                    color=col, fontsize=11, arrowprops=dict(arrowstyle="->", color=col))
    ax.axhline(0, color="k", lw=1)
    ax.axhspan(0, 40, color="#d4efdf", alpha=0.5); ax.text(1.1, 25, "0 dB より上 ＝ 制御が強く効く周波数", fontsize=11, color=GREEN)
    ax.set_ylim(-40, 35); ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("ゲイン [dB]"); ax.grid(alpha=0.3, which="both"); ax.legend(loc="lower left")
    ax.set_title("② ゲインを上げると、制御が効く範囲が高い周波数まで広がる", loc="left", fontsize=12)

    # (3) 半周期遅れると逆向きに効く
    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 0], hspace=0.5)
    t = np.linspace(0, 3, 600)
    x = np.sin(2 * np.pi * t)
    for r, (delay, lab) in enumerate([(0.0, "遅れなし：ずれと反対向きに押す → ゆれが収まる"),
                                       (0.5, "半周期（180°）遅れ：ずれと同じ向きに押してしまう → ゆれが大きくなる")]):
        ax = fig.add_subplot(sub[r])
        F = -0.75 * np.sin(2 * np.pi * (t - delay))
        ax.plot(t, x, color="k", lw=1.8, label="テーブルのずれ")
        ax.plot(t, F, color=RED if delay else GREEN, lw=2, ls="--", label="制御が押す力（戻そうとする力）")
        same = np.sign(F) == np.sign(x)
        ax.fill_between(t, -1.2, 1.2, where=same, color=RED, alpha=0.12, step="mid")
        ax.set_ylim(-1.3, 1.3); ax.set_yticks([]); ax.set_xlim(0, 3); ax.grid(alpha=0.3)
        ax.set_title(lab, loc="left", fontsize=11)
        if r == 0:
            ax.legend(loc="upper right", fontsize=9, ncol=2)
        else:
            ax.text(0.05, -1.15, "赤い帯 ＝ 押し戻すつもりが押し出している時間", color=RED, fontsize=10)
            ax.set_xlabel("時間（周期の何倍か）")
    fig.text(0.07, 0.475, "③ 制御が効く周波数で遅れが 180° を超えると、直そうとした操作が逆向きに効いて発振する", fontsize=12)

    # (4) 固有振動数を上げると遅れが高い周波数へ逃げる
    ax = fig.add_subplot(gs[1, 1])
    for fr, col, ls in [(30, RED, "-"), (120, GREEN, "--")]:
        L = open_loop(w, stiffness_for(fr), 100)
        ax.semilogx(w / 2 / np.pi, np.degrees(np.unwrap(np.angle(L))), ls, color=col, lw=2.2, label=f"1次固有振動数 {fr} Hz の機械")
        mag = 20 * np.log10(abs(L)); i = np.argmin(abs(mag)); fc = w[i] / 2 / np.pi
        ph = np.degrees(np.unwrap(np.angle(L)))[i]
        ax.plot(fc, ph, "o", color=col, ms=10)
        ax.annotate(f"0 dB の所 {fc:.0f} Hz で {ph:.0f}°", xy=(fc, ph), xytext=(fc * 1.4, ph + (-60 if fr == 30 else 50)),
                    color=col, fontsize=11, arrowprops=dict(arrowstyle="->", color=col))
    ax.axhline(-180, color="k", lw=1.2); ax.axhspan(-600, -180, color="#fadbd8", alpha=0.5)
    ax.text(1.1, -220, "−180° より下 ＝ 逆向きに効く領域", color=RED, fontsize=11)
    ax.set_ylim(-450, 0); ax.set_xlabel("周波数 [Hz]"); ax.set_ylabel("位相（遅れ）[deg]"); ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower left")
    ax.set_title("④ 1次固有振動数を上げると、大きな遅れが高い周波数へ逃げる（ゲイン 100 /s）", loc="left", fontsize=12)
    fig.savefig(os.path.join(OUT, "oscillation_mechanism.png")); plt.close(fig)

if __name__ == "__main__":
    natural_frequency()
    servo_gain_limit()
    cascade_loop()
    gain_meaning()
    oscillation_mechanism()
    print("saved to", os.path.abspath(OUT))
