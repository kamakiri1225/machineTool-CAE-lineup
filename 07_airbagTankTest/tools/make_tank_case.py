#!/usr/bin/env python3
"""Pyrotechnic inflator tank test, axisymmetric (5 deg wedge), reactingFoam (no chemistry, transient).

Domain (closed, adiabatic walls) [mm], axis = y:
  combustion chamber  r 0-15,  y -75..-5   (50 cc) : gas generation source (mass, energy, CO2/H2O)
  orifice             r 0-5.5, y -5..0     (d 11 mm, the total port area of the inflator)
  tank                r 0-170, y 0..660    (60 L, initially air 101.3 kPa, 293 K)
Gas generation: mdot(t) = M_gas t/tau^2 exp(-t/tau), gas at Tc (after the filter), mass fractions N2/CO2/H2O,
the same as the 0-D comparison model (tools/zero_d_inflator.py).
usage: python3 make_tank_case.py <caseDir> [json]
"""
import json, math, os, sys
import numpy as np

case = os.path.abspath(sys.argv[1])
P = dict(M_gas=0.025, tau=0.006, Tc=1000.0, YCO2=0.20, YH2O=0.10, tEnd=0.05, thermo="hConst")
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    object      {obj};\n}}\n\n'


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path); os.makedirs(os.path.dirname(full), exist_ok=True)
    open(full, "w").write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n")


R = [0, 5.5, 15, 170]
Y = [-75, -5, 0, 60, 660]
nR = [10, 10, 60]; gR = [1, 1, 10]
nY = [50, 10, 40, 80]; gY = [1, 1, 1, 6]
present = {(0, 0): "chamber", (1, 0): "chamber", (0, 1): "orifice",
           (0, 2): "tank", (1, 2): "tank", (2, 2): "tank", (0, 3): "tank", (1, 3): "tank", (2, 3): "tank"}
A = math.radians(2.5)
verts, vidx = [], {}


def V(i, j, s):
    r = R[i]; key = (i, j, 0 if r == 0 else s)
    if key not in vidx:
        vidx[key] = len(verts); verts.append((r * math.cos(A), Y[j], s * r * math.sin(A)))
    return vidx[key]


H = {k: [V(k[0], k[1], -1), V(k[0] + 1, k[1], -1), V(k[0] + 1, k[1] + 1, -1), V(k[0], k[1] + 1, -1),
         V(k[0], k[1], 1), V(k[0] + 1, k[1], 1), V(k[0] + 1, k[1] + 1, 1), V(k[0], k[1] + 1, 1)] for k in present}
FACE = {"ymin": (0, 4, 5, 1), "ymax": (3, 2, 6, 7), "rmin": (0, 3, 7, 4), "rmax": (1, 5, 6, 2), "front": (0, 1, 2, 3), "back": (4, 7, 6, 5)}
nb = {"ymin": (0, -1), "ymax": (0, 1), "rmax": (1, 0)}
walls = []
for (i, j), h in H.items():
    for side, (di, dj) in nb.items():
        if (i + di, j + dj) not in present:
            walls.append("(" + " ".join(str(h[k]) for k in FACE[side]) + ")")
blk = "\n".join(f"    hex ({' '.join(map(str, h))}) {present[k]} ({nR[k[0]]} {nY[k[1]]} 1) simpleGrading ({gR[k[0]] if k[1] >= 2 else 1} {gY[k[1]]} 1)" for k, h in H.items())
pat = lambda n, t, fs: f"    {n}\n    {{\n        type {t};\n        faces\n        (\n" + "\n".join("            " + x for x in fs) + "\n        );\n    }"
w("system/blockMeshDict", "scale 1e-3;\n\nvertices\n(\n" + "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in verts)
  + f"\n);\n\nblocks\n(\n{blk}\n);\n\nedges ();\n\nboundary\n(\n" + "\n".join([
      pat("walls", "wall", walls),
      pat("front", "wedge", ["(" + " ".join(str(h[k]) for k in FACE["front"]) + ")" for h in H.values()]),
      pat("back", "wedge", ["(" + " ".join(str(h[k]) for k in FACE["back"]) + ")" for h in H.values()])]) + "\n);\n")

# ---------------------------------------------------------------- thermo: air + combustion gas
JAN = {
    "N2":  (28.0134, "250", "5000", "( 2.92664 0.0014879768 -5.68476e-07 1.0097038e-10 -6.753351e-15 -922.7977 5.980528 )", "( 3.298677 0.0014082404 -3.963222e-06 5.641515e-09 -2.444854e-12 -1020.8999 3.950372 )", 3.2e-05),
    "O2":  (31.9988, "200", "3500", "( 3.28253784 0.00148308754 -7.57966669e-07 2.09470555e-10 -2.16717794e-14 -1088.45772 5.45323129 )", "( 3.78245636 -0.00299673416 9.84730201e-06 -9.68129509e-09 3.24372837e-12 -1063.94356 3.65767573 )", 3.6e-05),
    "CO2": (44.0095, "200", "5000", "( 4.45362 0.00314017 -1.27841e-06 2.394e-10 -1.66903e-14 -48967 -0.955396 )", "( 2.27572 0.00992207 -1.04091e-05 6.86669e-09 -2.11728e-12 -48373.1 10.1885 )", 3.0e-05),
    "H2O": (18.0153, "200", "5000", "( 2.67215 0.00305629 -8.73026e-07 1.201e-10 -6.39162e-15 -29899.2 6.86282 )", "( 3.38684 0.00347498 -6.3547e-06 6.96858e-09 -2.50659e-12 -30208.1 2.59023 )", 2.5e-05),
}
janaf = P["thermo"] == "janaf"
w("constant/thermophysicalProperties", f'''thermoType
{{
    type            hePsiThermo;
    mixture         reactingMixture;
    transport       {"sutherland" if janaf else "const"};
    thermo          {"janaf" if janaf else "hConst"};
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}}
inertSpecie     N2;
chemistryReader foamChemistryReader;
foamChemistryFile "<constant>/reactions";
foamChemistryThermoFile "<constant>/thermo.species";
''')
if janaf:   # cp(T)/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
    w("constant/thermo.species", "\n".join(
        f"{k}\n{{\n    specie\n    {{\n        molWeight       {v[0]};\n    }}\n"
        f"    thermodynamics\n    {{\n        Tlow            {v[1]};\n        Thigh           {v[2]};\n        Tcommon         1000;\n"
        f"        highCpCoeffs    {v[3]};\n        lowCpCoeffs     {v[4]};\n    }}\n"
        f"    transport\n    {{\n        As              1.512e-06;\n        Ts              120;\n    }}\n}}\n"
        for k, v in JAN.items()))
else:
    CONST = {"N2": (28.0134, 1170, 3.2e-05, 0.72), "O2": (31.9988, 1090, 3.6e-05, 0.72),
             "CO2": (44.0095, 1230, 3.0e-05, 0.72), "H2O": (18.0153, 2290, 2.5e-05, 0.9)}
    w("constant/thermo.species", "\n".join(
        f"{k}\n{{\n    specie\n    {{\n        molWeight       {v[0]};\n    }}\n"
        f"    thermodynamics\n    {{\n        Cp              {v[1]};\n        Hf              0;\n    }}\n"
        f"    transport\n    {{\n        mu              {v[2]};\n        Pr              {v[3]};\n    }}\n}}\n"
        for k, v in CONST.items()))
w("constant/reactions", "species\n(\n    N2\n    O2\n    CO2\n    H2O\n);\n\nreactions\n{\n}\n")
w("constant/chemistryProperties", "chemistryType\n{\n    solver          none;\n}\n\nchemistry       off;\n"
  "initialChemicalTimeStep 1e-07;\n")
w("constant/combustionProperties", "combustionModel none;\n")
w("constant/turbulenceProperties", "simulationType  RAS;\n\nRAS\n{\n    RASModel        kOmegaSST;\n    turbulence      on;\n"
  "    printCoeffs     on;\n}\n")

# gas generation source (per unit volume of the chamber): table in time
Vc = math.pi * (15e-3) ** 2 * 70e-3
tt = np.linspace(0, P["tEnd"], 401)
md = P["M_gas"] * tt / P["tau"] ** 2 * np.exp(-tt / P["tau"]) / Vc                     # kg/m3/s
hs = sum(Yk * cpk for Yk, cpk in ((1 - P["YCO2"] - P["YH2O"], 1170), (P["YCO2"], 1230), (P["YH2O"], 2290))) * (P["Tc"] - 298.15)
if P.get("thermo") == "janaf":
    import importlib.util
    spec = importlib.util.spec_from_file_location("zj", os.path.join(os.path.dirname(os.path.abspath(__file__)), "janaf_props.py"))
    zj = importlib.util.module_from_spec(spec); spec.loader.exec_module(zj)
    hs = sum(Yk * zj.hs_s(k, P["Tc"]) for k, Yk in (("N2", 1 - P["YCO2"] - P["YH2O"]), ("CO2", P["YCO2"]), ("H2O", P["YH2O"])))
tab = lambda f: "table\n                (\n" + "\n".join(f"                    ({a:.6g} {b:.6g})" for a, b in zip(tt, f)) + "\n                )"
w("constant/fvOptions", f'''gasGeneration
{{
    type            scalarSemiImplicitSource;
    selectionMode   cellZone;
    cellZone        chamber;
    volumeMode      specific;           // per m3 of the chamber
    sources
    {{
        rho
        {{
            explicit        {tab(md)};
            implicit        0;
        }}
        h
        {{
            explicit        {tab(md * hs)};
            implicit        0;
        }}
        CO2
        {{
            explicit        {tab(md * P["YCO2"])};
            implicit        0;
        }}
        H2O
        {{
            explicit        {tab(md * P["YH2O"])};
            implicit        0;
        }}
    }}
}}
limitT
{{
    type            limitTemperature;
    selectionMode   all;
    min             150;
    max             3000;
}}
''')
w("system/controlDict", f'''application     reactingFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {P["tEnd"]};
deltaT          1e-8;
writeControl    adjustable;
writeInterval   0.0005;
purgeWrite      0;
writeFormat     binary;
writePrecision  8;
timeFormat      general;
timePrecision   6;
adjustTimeStep  yes;
maxCo           0.5;
maxDeltaT       1e-5;
runTimeModifiable true;
functions
{{
    probes
    {{
        type            probes;
        libs            (sampling);
        fields          (p T);
        probeLocations  ((0.005 -0.040 0) (0.150 0.330 0) (0.010 0.600 0));
        writeControl    timeStep;
        writeInterval   20;
    }}
    tankAvg
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      cellZone;
        name            tank;
        operation       volAverage;
        fields          (p T CO2);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
    chamberAvg
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      cellZone;
        name            chamber;
        operation       volAverage;
        fields          (p T);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
    MachNo
    {{
        type            MachNo;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeControl    writeTime;
    }}
}}
''')
w("system/fvSchemes", '''ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwindV grad(U);
    div(phi,Yi_h)   Gauss upwind;
    div(phi,K)      Gauss upwind;
    div(phid,p)     Gauss upwind;
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}
''')
w("system/fvSolution", '''solvers
{
    "rho.*"
    {
        solver          diagonal;
    }
    "p.*"
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-10;      // strict: a loose tolerance leaves p unsolved in this case
        relTol          0.01;
        maxIter         200;
    }
    "pFinal"
    {
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       1e-10;
        relTol          0.001;
        maxIter         200;
    }
    "(U|h|k|omega|Yi).*"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-9;
        relTol          0;
    }
}
PIMPLE
{
    momentumPredictor   yes;
    transonic           yes;
    nOuterCorrectors    2;
    nCorrectors         2;
    nNonOrthogonalCorrectors 0;
}

relaxationFactors
{
    equations
    {
        ".*"            1;
    }
}
''')


def field(name, cls, dim, internal, wall):
    body = (f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n    walls\n    {{\n"
            + "".join(f"        {k:<16} {v};\n" for k, v in wall.items()) + "    }\n    \"(front|back)\"\n    {\n        type            wedge;\n    }\n}\n")
    w(f"0.orig/{name}", body, cls=cls)


S = "volScalarField"
field("p", S, "[1 -1 -2 0 0 0 0]", "uniform 101325", {"type": "zeroGradient"})
field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {"type": "noSlip"})
field("T", S, "[0 0 0 1 0 0 0]", "uniform 293.15", {"type": "zeroGradient"})
for sp, v in (("N2", 0.767), ("O2", 0.233), ("CO2", 0), ("H2O", 0), ("Ydefault", 0)):
    field(sp, S, "[0 0 0 0 0 0 0]", f"uniform {v}", {"type": "zeroGradient"})
field("k", S, "[0 2 -2 0 0 0 0]", "uniform 1e-4", {"type": "kqRWallFunction", "value": "uniform 1e-4"})
field("omega", S, "[0 0 -1 0 0 0 0]", "uniform 10", {"type": "omegaWallFunction", "value": "uniform 10"})
field("nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {"type": "nutkWallFunction", "value": "uniform 0"})
field("alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {"type": "compressible::alphatWallFunction", "value": "uniform 0"})
json.dump(dict(P, Vc=Vc), open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
open(os.path.join(case, "Allrun"), "w").write("#!/bin/sh\ncd \"${0%/*}\" || exit\nrm -rf 0 && cp -r 0.orig 0\nblockMesh > log.blockMesh 2>&1\ncheckMesh > log.checkMesh 2>&1\nreactingFoam > log.reactingFoam 2>&1\n")
os.chmod(os.path.join(case, "Allrun"), 0o755)
print("ok", case, "cells", sum(nR[i] * nY[j] for i, j in present), "Vc[cc]", Vc * 1e6)
