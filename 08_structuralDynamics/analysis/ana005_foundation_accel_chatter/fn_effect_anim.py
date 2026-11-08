"""ana005 4.2a のアニメーション（GIF）
   1) fn_effect_anim.gif      : 固有振動数 11.5 / 22.9 / 45.8 Hz の 3 台を並べ、同じ加減速（ジャーク 0）で本体がどう揺れるか
   2) jerk_effect_anim.gif    : 22.9 Hz の機械で、ジャーク 0 と 44 ms（周期と同じ）の比較
   本体の揺れ x_b は µm なので、図では 1 µm = 0.012（テーブル送り 50 mm = 1.6）として大きく拡大して描く"""
import sys
import numpy as np
from pathlib import Path
from io import BytesIO
from PIL import Image
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts")); sys.path.insert(0, str(D.parent / "tools"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
from matplotlib.patches import Rectangle
from gif_tools import save_gif
import fn_effect as FE                     # 同じモデル（profile / respond）を使う

UM = 0.012; RS = 1.6 / 0.05               # 表示倍率：1 µm → 0.012、テーブル 1 m → 32
T_END, DT_FRAME = 0.8, 0.01


def cases_sim(cases):
    out = []
    for fn, tj, col, lab in cases:
        t, x, tm, peak, settle = FE.respond(fn, tj)
        _, a, _ = FE.profile(tj)
        v = np.cumsum(a) * FE.dt; r = np.cumsum(v) * FE.dt  # 送り速度・位置 r(t)
        out.append(dict(t=t, x=x * 1e6, r=r, a=a, v=v, tm=tm, settle=tm + settle, col=col, lab=lab, peak=peak * 1e6))
    return out


def draw_machine(ax, xb, r, col):
    ax.set_xlim(-0.6, 3.6); ax.set_ylim(-0.45, 2.3); ax.set_aspect("equal"); ax.axis("off")
    ax.plot([-0.6, 3.6], [0, 0], "k", lw=2)
    for xx in np.arange(-0.6, 3.6, 0.2): ax.plot([xx, xx - 0.12], [0, -0.12], "k", lw=0.6)
    s = xb * UM
    for xf in (0.2, 2.6):                                  # 脚と横ばね
        ax.add_patch(Rectangle((xf - 0.15, 0), 0.3, 0.1, fc="#7f8c8d"))
        xs = np.linspace(xf - 0.15, xf + 0.15 + s, 9); yy = 0.2 + np.array([0, .05, -.05, .05, -.05, .05, -.05, .05, 0])
        ax.plot(xs, yy, color=BLUE, lw=1.8)
    ax.add_patch(Rectangle((0 + s, 0.3), 3.0, 0.45, fc="#d5d8dc", ec="k"))          # ベッド
    ax.add_patch(Rectangle((2.2 + s, 0.75), 0.6, 1.3, fc="#d5d8dc", ec="k"))        # コラム
    ax.add_patch(Rectangle((1.5 + s, 1.75), 0.7, 0.3, fc="#d5d8dc", ec="k"))        # 主軸頭
    ax.add_patch(Rectangle((0.1 + s + r * RS, 0.75), 0.7, 0.25, fc="#f5cba7", ec="k"))  # テーブル
    ax.plot([0, 0], [0.25, 2.2], color=GRAY, ls=":", lw=1)                            # 元の位置
    ax.annotate("", xy=(s, 2.15), xytext=(0, 2.15), arrowprops=dict(arrowstyle="-|>", color=col, lw=2.5, mutation_scale=14)) if abs(s) > 0.02 else None


def make_gif(cases, title, path, sub):
    sims = cases_sim(cases); n = len(sims)
    frames = []
    for tf in np.arange(0, T_END + 1e-9, DT_FRAME):
        fig = plt.figure(figsize=(max(5.2 * n, 11), 9.6))
        gs = fig.add_gridspec(3, n, height_ratios=[1.05, 0.62, 1])
        axc = fig.add_subplot(gs[1, :]); axv = axc.twinx()
        axb = fig.add_subplot(gs[2, :], sharex=axc)
        for k, sm in enumerate(sims):
            i = min(np.searchsorted(sm["t"], tf), len(sm["t"]) - 1)
            ax = fig.add_subplot(gs[0, k]); draw_machine(ax, sm["x"][i], sm["r"][i], sm["col"])
            moving = tf < sm["tm"]; ok = tf >= sm["settle"]
            state = "送り中" if moving else ("±1 µm 以内 → 加工できる" if ok else "停止したが揺れている")
            ax.set_title(f"{sm['lab']}\n{state}", fontsize=12, color=sm["col"] if not ok or moving else GREEN)
            ax.text(-0.55, -0.42, f"$x_b$ = {sm['x'][i]:+6.1f} µm", fontsize=11, family="monospace", color=sm["col"])
            m = sm["t"] <= tf
            ls = "-" if k == 0 else "--"
            axc.plot(sm["t"], sm["a"], color=sm["col"], lw=1, alpha=0.25)
            axc.plot(sm["t"][m], sm["a"][m], color=sm["col"], lw=2, ls=ls, label=f"加速度指令 {sm['lab']}")
            axv.plot(sm["t"][m], sm["v"][m] * 1000, color=sm["col"], lw=1.2, ls=":", alpha=0.9)
            axb.plot(sm["t"][m], sm["x"][m], color=sm["col"], lw=2, label=sm["lab"])
            axb.plot([tf], [sm["x"][i]], "o", color=sm["col"], ms=6)
            axb.axvline(sm["settle"], color=sm["col"], ls="--", lw=1, alpha=0.6 if tf >= sm["settle"] else 0)
        axc.set_ylim(-2.6, 2.6); axc.set_ylabel("加速度指令 $\\ddot r$ [m/s²]"); axc.grid(alpha=0.3); axc.axvline(tf, color="k", lw=0.8)
        axv.set_ylim(-260, 260); axv.set_ylabel("速度 [mm/s]（点線）", color=GRAY)
        axc.legend(loc="upper right", fontsize=9); axc.tick_params(labelbottom=False)
        axc.set_title("送り指令（実線 = 加速度、点線 = 速度。薄い線 = これからの指令）", loc="left", fontsize=11)
        axb.axhspan(-1, 1, color=GREEN, alpha=0.12); axb.text(T_END - 0.01, 1.5, "±1 µm", color=GREEN, ha="right")
        lim = max(abs(s["x"]).max() for s in sims) * 1.1
        axb.set_xlim(0, T_END); axb.set_ylim(-lim, lim); axb.grid(alpha=0.3); axb.legend(loc="upper right", fontsize=10)
        axb.set_xlabel("時刻 [s]"); axb.set_ylabel("本体の X 変位 $x_b$ [µm]"); axb.axvline(tf, color="k", lw=0.8)
        axb.set_title("縦の破線 = 揺れが ±1 µm に収まった時刻（次の加工を始められる）", loc="left", fontsize=11)
        fig.suptitle(f"{title}　t = {tf:.2f} s\n{sub}", fontsize=13)
        fig.tight_layout(); buf = BytesIO(); fig.savefig(buf, format="png", dpi=80); plt.close(fig)
        frames.append(Image.open(buf).convert("RGB"))
    frames += [frames[-1]] * 15
    save_gif(frames, str(path), duration=80); print("saved", path, len(frames))


make_gif([(11.5, 0.0, RED, "$f_n$ = 11.5 Hz"), (22.9, 0.0, BLUE, "$f_n$ = 22.9 Hz（今）"), (45.8, 0.0, GREEN, "$f_n$ = 45.8 Hz")],
         "固有振動数で比べる（同じ 50 mm 移動・2 m/s²・ジャーク 0）", D / "results" / "fn_effect_anim.gif",
         "本体の揺れは µm を大きく拡大して表示（テーブルの送りは縮小）")
make_gif([(22.9, 0.0, RED, "ジャーク 0 ms（段差）"), (22.9, 0.044, GREEN, "ジャーク 44 ms（周期と同じ）")],
         "ジャークで比べる（$f_n$ = 22.9 Hz の機械）", D / "results" / "jerk_effect_anim.gif",
         "ジャークを周期と同じ長さにすると、揺れがほぼ打ち消される（移動は 0.04 s 長くなる）")
