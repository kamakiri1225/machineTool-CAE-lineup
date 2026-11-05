#!/usr/bin/env python3
"""0-D inflator tank model with temperature-dependent specific heat (NASA / JANAF 7-coefficient polynomials):
    cp(T)/R_s = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4          (4th-order polynomial in T)
    h(T)/R_s  = a1 T + a2 T^2/2 + a3 T^3/3 + a4 T^4/4 + a5 T^5/5 + a6
Low/high coefficient sets switch at Tcommon = 1000 K (same data as OpenFOAM's janaf thermo, so the CFD can use them too).
Combustion chamber: gas enters at Tc (h_gas(Tc)); chamber temperature is computed from its energy balance (adiabatic).
Tank: adiabatic, air + combustion gas, T solved from the internal energy with Newton iterations.
usage: zero_d_inflator_janaf.py <outdir>
"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17, "axes.titlesize": 19, "axes.labelsize": 17})
OUT = sys.argv[1] if len(sys.argv) > 1 else "zeroD_janaf"; os.makedirs(OUT, exist_ok=True)
JANAF = {  # W [g/mol], low (200/250-1000 K), high (1000-5000 K)
    "N2": (28.0134, [3.298677, 0.0014082404, -3.963222e-06, 5.641515e-09, -2.444854e-12, -1020.8999, 3.950372],
           [2.92664, 0.0014879768, -5.68476e-07, 1.0097038e-10, -6.753351e-15, -922.7977, 5.980528]),
    "O2": (31.9988, [3.78245636, -0.00299673416, 9.84730201e-06, -9.68129509e-09, 3.24372837e-12, -1063.94356, 3.65767573],
           [3.28253784, 0.00148308754, -7.57966669e-07, 2.09470555e-10, -2.16717794e-14, -1088.45772, 5.45323129]),
    "CO2": (44.0095, [2.27572, 0.00992207, -1.04091e-05, 6.86669e-09, -2.11728e-12, -48373.1, 10.1885],
            [4.45362, 0.00314017, -1.27841e-06, 2.394e-10, -1.66903e-14, -48967, -0.955396]),
    "H2O": (18.0153, [3.38684, 0.00347498, -6.3547e-06, 6.96858e-09, -2.50659e-12, -30208.1, 2.59023],
            [2.67215, 0.00305629, -8.73026e-07, 1.201e-10, -6.39162e-15, -29899.2, 6.86282]),
}
Ru = 8.314462618
P = dict(M_gas=0.025, tau=0.006, Vc=49.48e-6, Tc=1000.0, Do=11.0e-3, Cd=0.85, Vt=0.060, pt0=101325.0, Tt0=293.15, tEnd=0.1)
Ygas = {"N2": 0.70, "CO2": 0.20, "H2O": 0.10}; Yair = {"N2": 0.767, "O2": 0.233}


def Rs(sp): return Ru / (JANAF[sp][0] * 1e-3)


def cp_s(sp, T):
    a = JANAF[sp][1] if T < 1000 else JANAF[sp][2]
    return Rs(sp) * (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3 + a[4] * T**4)


def hs_s(sp, T):          # sensible enthalpy relative to 298.15 K
    def h(T):
        a = JANAF[sp][1] if T < 1000 else JANAF[sp][2]
        return Rs(sp) * (a[0] * T + a[1] * T**2 / 2 + a[2] * T**3 / 3 + a[3] * T**4 / 4 + a[4] * T**5 / 5 + a[5])
    return h(T) - h(298.15)


def mix(Y, f, T): return sum(Y[k] * f(k, T) for k in Y)
def Rmix(Y): return sum(Y[k] * Rs(k) for k in Y)


def T_from_e(Y, e, Tg):   # e = hs - R T  (sensible internal energy rel. 298.15 K, consistent with hs)
    R = Rmix(Y); T = Tg
    for _ in range(30):
        f = mix(Y, hs_s, T) - R * T - e; T -= f / (mix(Y, cp_s, T) - R)
        if abs(f) < 1e-6: break
    return T


Ao = np.pi * P["Do"] ** 2 / 4
dt = 2e-6; n = int(P["tEnd"] / dt) + 1; t = np.arange(n) * dt
# chamber: mass, energy, composition (initially air at 293 K, 1 atm)
mc_a = P["pt0"] * P["Vc"] / (Rmix(Yair) * P["Tt0"]); mc_g = 0.0
Ec = mc_a * (mix(Yair, hs_s, P["Tt0"]) - Rmix(Yair) * P["Tt0"])
ma = P["pt0"] * P["Vt"] / (Rmix(Yair) * P["Tt0"]); mg = 0.0
Et = ma * (mix(Yair, hs_s, P["Tt0"]) - Rmix(Yair) * P["Tt0"])
hgen = mix(Ygas, hs_s, P["Tc"])
pc = np.zeros(n); Tch = np.zeros(n); mo_ = np.zeros(n); pt = np.zeros(n); Tt = np.zeros(n)
Tcg, Ttg = P["Tt0"], P["Tt0"]
for i, ti in enumerate(t):
    gen = P["M_gas"] * ti / P["tau"] ** 2 * np.exp(-ti / P["tau"])
    mc = mc_a + mc_g; Yc = {k: (mc_a * Yair.get(k, 0) + mc_g * Ygas.get(k, 0)) / mc for k in ("N2", "O2", "CO2", "H2O")}
    Tcg = T_from_e(Yc, Ec / mc, Tcg); Rc = Rmix(Yc); p_c = mc * Rc * Tcg / P["Vc"]
    mt = ma + mg; Yt = {k: (ma * Yair.get(k, 0) + mg * Ygas.get(k, 0)) / mt for k in ("N2", "O2", "CO2", "H2O")}
    Ttg = T_from_e(Yt, Et / mt, Ttg); p_t = mt * Rmix(Yt) * Ttg / P["Vt"]
    cpc = mix(Yc, cp_s, Tcg); g = cpc / (cpc - Rc)
    pr = p_t / p_c; crit = (2 / (g + 1)) ** (g / (g - 1))
    if p_c <= p_t: mo = 0.0
    elif pr <= crit: mo = P["Cd"] * Ao * p_c * np.sqrt(g / (Rc * Tcg)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    else: mo = P["Cd"] * Ao * p_c * np.sqrt(2 * g / (Rc * Tcg * (g - 1)) * (pr ** (2 / g) - pr ** ((g + 1) / g)))
    pc[i], Tch[i], mo_[i], pt[i], Tt[i] = p_c, Tcg, mo, p_t, Ttg
    hc = mix(Yc, hs_s, Tcg)                                  # outflow carries the chamber enthalpy
    fa = mc_a / mc
    mc_a += -mo * fa * dt; mc_g += (gen - mo * (1 - fa)) * dt
    Ec += (gen * hgen - mo * hc) * dt
    ma += mo * fa * dt; mg += mo * (1 - fa) * dt; Et += mo * hc * dt
np.savetxt(os.path.join(OUT, "zero_d_result.csv"), np.c_[t, pc, mo_, pt, Tt, Tch], delimiter=",",
           header="t[s],p_chamber[Pa],mdot_orifice[kg/s],p_tank[Pa],T_tank[K],T_chamber[K]", comments="")
json.dump(dict(P, peak_chamber_MPa=pc.max() / 1e6, peak_tank_kPa=(pt.max() - P["pt0"]) / 1e3, T_tank_end=Tt[-1]),
          open(os.path.join(OUT, "zero_d_summary.json"), "w"), indent=1)
# cp(T) figure for the doc
TT = np.linspace(250, 1500, 300)
fig, ax = plt.subplots(figsize=(10, 6.2))
for sp, c in (("N2", "#1565C0"), ("O2", "#2E7D32"), ("CO2", "#C62828"), ("H2O", "#EF6C00")):
    ax.plot(TT, [cp_s(sp, x) for x in TT], lw=3, color=c, label=sp)
ax.plot(TT, [mix(Yair, cp_s, x) for x in TT], "k--", lw=2.5, label="空気（N₂+O₂）")
ax.set_xlabel("温度 [K]"); ax.set_ylabel("定圧比熱 cp [J/kgK]"); ax.set_title("JANAF の 4 次式で表した比熱の温度依存性"); ax.grid(alpha=0.3)
ax.legend(frameon=False, ncol=3, loc="upper left"); fig.tight_layout(); fig.savefig(os.path.join(OUT, "cp_janaf.png"), dpi=90, facecolor="white")
print(f"chamber peak {pc.max()/1e6:.2f} MPa, tank peak {(pt.max()-P['pt0'])/1e3:.1f} kPa(G), T_tank end {Tt[-1]:.0f} K, "
      f"air cp(293)={mix(Yair, cp_s, 293):.0f}, gas cp(1000)={mix(Ygas, cp_s, 1000):.0f}")
