#!/usr/bin/env python3
"""0-D model of a pyrotechnic airbag inflator fired into a 60 L tank (tank test).

Combustion chamber (volume Vc, gas temperature Tc after the filter, fixed):
    d(m_c)/dt = mdot_gen(t) - mdot_out
    p_c = m_c R Tc / Vc
    mdot_gen(t) = M_gas * (t/tau^2) * exp(-t/tau)            (gas generation by the burning propellant)
Orifice (area Ao, Cd): isentropic nozzle, choked while p_c/p_t > critical ratio.
Tank (volume Vt, initially air at 101.3 kPa, 293 K, adiabatic walls):
    mass and internal energy balance with the hot inflator gas (enthalpy cp*Tc), ideal-gas mixture.
Writes: inflator_chamber_p0.csv (t, p_c) for the CFD inlet, and zero_d_result.csv / .png.
"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 17, "axes.titlesize": 19, "axes.labelsize": 17})

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)
P = dict(M_gas=0.025, tau=0.006, Vc=50e-6, Tc=1000.0, Do=11.0e-3, Cd=0.85, Vt=0.060, pt0=101325.0, Tt0=293.15, tEnd=0.1,
         Y={"N2": 0.70, "CO2": 0.20, "H2O": 0.10})
Ru = 8.314462618
W = {"N2": 28.0134, "CO2": 44.0095, "H2O": 18.0153, "O2": 31.9988}
cp = {"N2": 1170.0, "CO2": 1230.0, "H2O": 2290.0, "O2": 1090.0}        # J/kgK, representative 300-1400 K
Ri = sum(P["Y"][k] / (W[k] * 1e-3) for k in P["Y"]) * Ru
cpi = sum(P["Y"][k] * cp[k] for k in P["Y"])
gi = cpi / (cpi - Ri)
air = {"N2": 0.767, "O2": 0.233}
Ra = sum(air[k] / (W[k] * 1e-3) for k in air) * Ru
cpa = sum(air[k] * (1005.0 if k else 0) for k in air) / 1.0
cpa = sum(air[k] * cp[k] for k in air)          # same species cp as the CFD (N2 1170, O2 1090) -> 1151 J/kgK
Ao = np.pi * P["Do"] ** 2 / 4

dt = 2e-6
n = int(P["tEnd"] / dt) + 1
t = np.arange(n) * dt
mc = 0.0
pc = np.zeros(n); mout = np.zeros(n); pt = np.zeros(n); Tt = np.zeros(n)
ma = P["pt0"] * P["Vt"] / (Ra * P["Tt0"]); mg = 0.0
Ua = ma * (cpa - Ra) * P["Tt0"]; Ug = 0.0
for i, ti in enumerate(t):
    gen = P["M_gas"] * ti / P["tau"] ** 2 * np.exp(-ti / P["tau"])
    p_c = mc * Ri * P["Tc"] / P["Vc"]
    # tank state: mixture
    Rmix = (ma * Ra + mg * Ri) / (ma + mg); cvmix = (ma * (cpa - Ra) + mg * (cpi - Ri)) / (ma + mg)
    T_t = (Ua + Ug) / ((ma + mg) * cvmix); p_t = (ma + mg) * Rmix * T_t / P["Vt"]
    pr = p_t / max(p_c, 1.0); crit = (2 / (gi + 1)) ** (gi / (gi - 1))
    if p_c <= p_t:
        mo = 0.0
    elif pr <= crit:
        mo = P["Cd"] * Ao * p_c * np.sqrt(gi / (Ri * P["Tc"])) * (2 / (gi + 1)) ** ((gi + 1) / (2 * (gi - 1)))
    else:
        mo = P["Cd"] * Ao * p_c * np.sqrt(2 * gi / (Ri * P["Tc"] * (gi - 1)) * (pr ** (2 / gi) - pr ** ((gi + 1) / gi)))
    pc[i], mout[i], pt[i], Tt[i] = p_c, mo, p_t, T_t
    mc = max(mc + (gen - mo) * dt, 0.0)
    mg += mo * dt; Ug += mo * cpi * P["Tc"] * dt          # enthalpy inflow
np.savetxt(os.path.join(OUT, "inflator_chamber_p0.csv"), np.c_[t[::50], pc[::50] + 0.0], delimiter=",", header="t[s],p_chamber[Pa]", comments="")
np.savetxt(os.path.join(OUT, "zero_d_result.csv"), np.c_[t, pc, mout, pt, Tt], delimiter=",",
           header="t[s],p_chamber[Pa],mdot_orifice[kg/s],p_tank[Pa],T_tank[K]", comments="")
fig, ax = plt.subplots(1, 3, figsize=(18, 6))
ax[0].plot(t * 1e3, pc / 1e6, lw=3, color="#C62828"); ax[0].set_ylabel("燃焼室の圧力 [MPa]"); ax[0].set_title("インフレーターの燃焼室")
ax[1].plot(t * 1e3, mout * 1e3, lw=3, color="#EF6C00"); ax[1].set_ylabel("オリフィスの流量 [g/s]"); ax[1].set_title("オリフィスを通る流量")
ax[2].plot(t * 1e3, (pt - P["pt0"]) / 1e3, lw=3, color="#1565C0"); ax[2].set_ylabel("タンクのゲージ圧 [kPa]"); ax[2].set_title("60 L タンクの圧力（タンク試験）")
for a in ax:
    a.set_xlabel("時間 [ms]"); a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "zero_d_result.png"), dpi=90, facecolor="white")
json.dump(dict(P, Ri=Ri, cp=cpi, gamma=gi, peak_chamber_MPa=pc.max() / 1e6, peak_tank_kPa=(pt.max() - P["pt0"]) / 1e3,
               gas_out_g=mg * 1e3, T_tank_peak=Tt.max()), open(os.path.join(OUT, "zero_d_summary.json"), "w"), indent=1)
print(f"chamber peak {pc.max()/1e6:.2f} MPa at {t[pc.argmax()]*1e3:.1f} ms, tank peak {(pt.max()-P['pt0'])/1e3:.1f} kPa(G), "
      f"T_tank max {Tt.max():.0f} K, gas out {mg*1e3:.1f} g, gamma {gi:.3f}")
