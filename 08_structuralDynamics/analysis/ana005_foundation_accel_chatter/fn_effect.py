"""ana005 4.4：基礎の固有振動数を上げると、加減速後の揺れと待ち時間はどう変わるか（4.1 と同じ 1 自由度モデル）
   M・m_t・送り条件（50 mm、0.2 m/s、2 m/s^2）は同じ。脚の横剛性 k だけを変えて f_n を 11.5 / 23 / 46 Hz にする。減衰比は 10.5 % で固定"""
import sys, json
import numpy as np
from pathlib import Path
from scipy.signal import lsim
D = Path(__file__).resolve().parent
sys.path.insert(0, str(D.parent.parent / "study" / "scripts"))
from figstyle import plt, BLUE, RED, GRAY, ORANGE, GREEN
acc = json.loads((D / "results" / "acceleration.json").read_text())
M, mt, z = acc["M_kg"], acc["moving_mass_kg"], acc["zeta"]
L, V, A = 0.05, 0.2, 2.0
dt = 2e-5


def profile(tj):
    ta = V / A + tj; tc = L / V - ta; tm = 2 * ta + tc
    t = np.arange(0, tm + 1.5, dt)
    def a1(q):
        if tj == 0: return A
        return np.where(q < tj, A * q / tj, np.where(q < ta - tj, A, A * (ta - q) / tj))
    a = np.zeros_like(t)
    m1 = t < ta; a[m1] = a1(t[m1])
    m3 = (t >= ta + tc) & (t < tm); a[m3] = -a1(t[m3] - ta - tc)
    return t, a, tm


def respond(fn, tj):
    w = 2 * np.pi * fn; k = M * w**2; c = 2 * z * np.sqrt(k * M)
    t, a, tm = profile(tj)
    _, x, _ = lsim(([-mt], [M, c, k]), a, t)
    after = t >= tm
    peak = np.abs(x[after]).max()
    late = t[after][np.abs(x[after]) > 1e-6]
    settle = 0.0 if len(late) == 0 else late[-1] - tm
    return t, x, tm, peak, settle


fns = [11.5, 22.9, 45.8]; cols = [RED, BLUE, GREEN]
tjs = np.linspace(0, 0.1, 51)
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
rows = []
for fn, col in zip(fns, cols):
    pk, st, tot = [], [], []
    for tj in tjs:
        _, _, tm, p, s = respond(fn, tj); pk.append(p * 1e6); st.append(s); tot.append(tm + s)
    lab = f"$f_n$ = {fn:.1f} Hz（周期 {1000/fn:.0f} ms）"
    axs[0].plot(tjs * 1000, pk, color=col, lw=2.2, label=lab)
    axs[1].plot(tjs * 1000, tot, color=col, lw=2.2, label=lab)
    for kk in (1, 2):
        axs[0].axvline(1000 / fn * kk, color=col, ls=":", lw=1)
    i = int(np.argmin(tot)); rows.append((fn, pk[0], st[0], tot[0], tjs[i] * 1000, tot[i]))
axs[0].set_yscale("log"); axs[0].set_ylim(1e-2, 50); axs[0].set_xlim(0, 100); axs[0].axhline(1, color=GRAY, ls="--"); axs[0].text(98, 1.1, "1 µm", ha="right", color=GRAY)
axs[0].set_xlabel("ジャーク立上げ時間 [ms]"); axs[0].set_ylabel("停止後の本体の揺れの最大 [µm]"); axs[0].grid(alpha=0.3, which="both"); axs[0].legend(fontsize=10)
axs[0].set_title("① 停止後に残る揺れ（点線 = 周期の 1・2 倍）", loc="left", fontsize=12)
axs[1].set_xlabel("ジャーク立上げ時間 [ms]"); axs[1].set_ylabel("移動 ＋ 揺れが ±1 µm に収まるまで [s]"); axs[1].grid(alpha=0.3); axs[1].legend(fontsize=10)
axs[1].set_title("② 次の加工を始められるまでの時間", loc="left", fontsize=12)
ax = axs[2]
for fn, col in zip(fns, cols):
    t, x, tm, p, s = respond(fn, 0.0)
    ax.plot(t, x * 1e6, color=col, lw=1.6, label=f"$f_n$ = {fn:.1f} Hz")
ax.axhline(1, color=GRAY, ls=":"); ax.axhline(-1, color=GRAY, ls=":"); ax.set_xlim(0, 1); ax.grid(alpha=0.3); ax.legend(fontsize=10)
ax.set_xlabel("時刻 [s]"); ax.set_ylabel("本体の X 変位 $x_b$ [µm]")
ax.set_title("③ 加速度段差（ジャーク 0）のときの揺れ", loc="left", fontsize=12)
fig.suptitle("基礎の固有振動数を上げると：同じ 50 mm 移動（2 m/s²）で、揺れが小さく・早く収まり・短いジャークで済む", fontsize=13)
fig.tight_layout(); fig.savefig(D / "results" / "fn_effect.png")
md = ["| 固有振動数 | ジャーク 0 の残留揺れ [µm] | ジャーク 0 の整定 [s] | ジャーク 0 の移動＋整定 [s] | 最短になるジャーク [ms] | そのときの移動＋整定 [s] |", "|---|---|---|---|---|---|"]
md += [f"| {r[0]:.1f} Hz | {r[1]:.2f} | {r[2]:.3f} | {r[3]:.3f} | {r[4]:.0f} | {r[5]:.3f} |" for r in rows]
(D / "results" / "fn_effect.md").write_text("\n".join(md) + "\n"); print("\n".join(md))
