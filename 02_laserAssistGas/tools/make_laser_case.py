#!/usr/bin/env python3
"""Generate an axisymmetric laser-cutting assist-gas case (reactingFoam, LTS pseudo-steady).

Geometry (axis = y, radial = x, 2.5 deg wedge in z), all in mm:
    nozzle bore: conical, inlet radius RIN at y = H+LN  ->  exit radius RN at y = H
    nozzle tip face: r = RN..RNOZ at y = H, outer side r = RNOZ (y = H..H+LN)
    stand-off gap: y = 0..H
    steel plate: y = -T..0, kerf modelled as a round hole r < RK (axisymmetric approximation of the slot)
    region below the plate: y = -T-LB..-T (melt/gas exit side)
Gas: O2 / N2 mixture (constant-cp perfect gas, valid at the low temperatures of supersonic expansion).

usage: python3 make_laser_case.py <caseDir> <json-params>
"""
import json
import math
import os
import sys

case = os.path.abspath(sys.argv[1])
P = dict(RK=0.25, RN=0.75, RIN=1.5, RNOZ=3.0, ROUT=8.0, H=1.0, LN=3.0, T=6.0, LB=4.0,
         p0_gauge=0.1e6, YO2_in=1.0, endTime=8000, title="")
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
PAMB, T0 = 101325.0, 300.0
P0 = PAMB + P["p0_gauge"]
YO2_AIR, YN2_AIR = 0.233, 0.767

HDR = '''FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

'''


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w") as f:
        f.write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n// ************************************************************************* //\n")


# ------------------------------------------------------------------ mesh
RK, RN, RIN, RNOZ, ROUT = P["RK"], P["RN"], P["RIN"], P["RNOZ"], P["ROUT"]
H, LN, T, LB = P["H"], P["LN"], P["T"], P["LB"]
Y0, Y1, Y2, Y3, Y4 = -T - LB, -T, 0.0, H, H + LN
RK4 = RK * RIN / RN                     # radial split line inside the cone at the nozzle inlet
A = math.radians(2.5)
verts, vidx = [], {}


def V(r, y, s):
    key = (round(r, 9), round(y, 9), 0 if r == 0 else s)
    if key not in vidx:
        vidx[key] = len(verts)
        verts.append((r * math.cos(A), y, s * r * math.sin(A)))
    return vidx[key]


def hexq(r00, r10, r11, r01, ya, yb):
    """quad in (r,y): bottom (r00 -> r10 at ya), top (r01 -> r11 at yb)"""
    return [V(r00, ya, -1), V(r10, ya, -1), V(r11, yb, -1), V(r01, yb, -1),
            V(r00, ya, 1), V(r10, ya, 1), V(r11, yb, 1), V(r01, yb, 1)]


def n_(L, d):
    return max(4, int(round(L / d)))


dr_k = RK / 10
cells = {"rk": 10, "rn": n_(RN - RK, dr_k * 1.2), "rnoz": n_(RNOZ - RN, 0.06), "rout": n_(ROUT - RNOZ, 0.25),
         "below": n_(LB, 0.12), "kerf": n_(T, 0.04), "gap": n_(H, 0.025), "noz": n_(LN, 0.05)}
B = []   # (name, hex, (nr, ny), grading)
radial = [(0, RK, "rk", "1"), (RK, RN, "rn", "1"), (RN, RNOZ, "rnoz", "3"), (RNOZ, ROUT, "rout", "4")]
for r0, r1, key, g in radial:
    B.append((f"below_{key}", hexq(r0, r1, r1, r0, Y0, Y1), (cells[key], cells["below"]), f"{g} 0.25 1"))
B.append(("kerf", hexq(0, RK, RK, 0, Y1, Y2), (cells["rk"], cells["kerf"]), "1 1 1"))
for r0, r1, key, g in radial:
    B.append((f"gap_{key}", hexq(r0, r1, r1, r0, Y2, Y3), (cells[key], cells["gap"]), f"{g} 1 1"))
B.append(("noz_rk", hexq(0, RK, RK4, 0, Y3, Y4), (cells["rk"], cells["noz"]), "1 1 1"))
B.append(("noz_rn", hexq(RK, RN, RIN, RK4, Y3, Y4), (cells["rn"], cells["noz"]), "1 1 1"))
B.append(("above", hexq(RNOZ, ROUT, ROUT, RNOZ, Y3, Y4), (cells["rout"], cells["noz"]), "4 2 1"))
H_ = {n: h for n, h, _, _ in B}

FACE = {"ymin": (0, 4, 5, 1), "ymax": (3, 2, 6, 7), "rmin": (0, 3, 7, 4), "rmax": (1, 5, 6, 2),
        "front": (0, 1, 2, 3), "back": (4, 7, 6, 5)}


def F(name, side):
    h = H_[name]
    return "(" + " ".join(str(h[i]) for i in FACE[side]) + ")"


patches = {
    "inlet": ("patch", [F("noz_rk", "ymax"), F("noz_rn", "ymax")]),
    "nozzle": ("wall", [F("noz_rn", "rmax"), F("gap_rnoz", "ymax"), F("above", "rmin")]),
    "plate": ("wall", [F("gap_rn", "ymin"), F("gap_rnoz", "ymin"), F("gap_rout", "ymin"), F("kerf", "rmax"),
                       F("below_rn", "ymax"), F("below_rnoz", "ymax"), F("below_rout", "ymax")]),
    "ambient": ("patch", [F("gap_rout", "rmax"), F("above", "rmax"), F("above", "ymax"),
                          F("below_rout", "rmax")] + [F(f"below_{k}", "ymin") for k in ("rk", "rn", "rnoz", "rout")]),
    "front": ("wedge", [F(n, "front") for n in H_]),
    "back": ("wedge", [F(n, "back") for n in H_]),
}
blk = "\n".join(f"    hex ({' '.join(map(str, h))}) ({nr} {ny} 1) simpleGrading ({g})" for _, h, (nr, ny), g in B)
pat = "\n".join(f"    {n}\n    {{\n        type {t};\n        faces\n        (\n" + "\n".join("            " + x for x in fs) + "\n        );\n    }"
                for n, (t, fs) in patches.items())
w("system/blockMeshDict", "scale   1e-3;\n\nvertices\n(\n" + "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in verts)
  + f"\n);\n\nblocks\n(\n{blk}\n);\n\nedges\n(\n);\n\nboundary\n(\n{pat}\n);\n")

# ------------------------------------------------------------------ constant
w("constant/thermophysicalProperties", '''thermoType
{
    type            hePsiThermo;
    mixture         reactingMixture;
    transport       const;
    thermo          hConst;
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}

inertSpecie     N2;
chemistryReader foamChemistryReader;
foamChemistryFile "<constant>/reactions";
foamChemistryThermoFile "<constant>/thermo.hConst";
''')
w("constant/thermo.hConst", '''O2
{
    specie          { molWeight 31.9988; }
    thermodynamics  { Cp 918; Hf 0; }
    transport       { mu 2.06e-05; Pr 0.72; }
}
N2
{
    specie          { molWeight 28.0134; }
    thermodynamics  { Cp 1040; Hf 0; }
    transport       { mu 1.79e-05; Pr 0.71; }
}
''')
w("constant/reactions", "species\n(\n    O2\n    N2\n);\n\nreactions\n{\n}\n")
w("constant/chemistryProperties", "chemistryType\n{\n    solver          none;\n}\n\nchemistry       off;\n\ninitialChemicalTimeStep 1e-07;\n")
w("constant/combustionProperties", "combustionModel none;\n")
w("constant/turbulenceProperties", "simulationType  RAS;\n\nRAS\n{\n    RASModel        kOmegaSST;\n    turbulence      on;\n    printCoeffs     on;\n}\n")

# ------------------------------------------------------------------ system
w("system/controlDict", f'''application     reactingFoam;

startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {P["endTime"]};      // LTS pseudo-time iterations
deltaT          1;
writeControl    timeStep;
writeInterval   1000;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    MachNo
    {{
        type            MachNo;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeControl    writeTime;
    }}
    massFlow
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            inlet;
        operation       sum;
        fields          (phi);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
    platePressure
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            plate;
        operation       max;
        fields          (p);
        writeFields     false;
        writeControl    timeStep;
        writeInterval   20;
    }}
}}
''')
w("system/fvSchemes", '''ddtSchemes      { default localEuler; }
gradSchemes     { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
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
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.5; }
wallDist        { method meshWave; }
''')
w("system/fvSolution", '''solvers
{
    "rho.*"         { solver diagonal; }
    "p.*"           { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.05; }
    "(U|h|k|omega|Yi).*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.1; }
}

PIMPLE
{
    momentumPredictor   yes;
    transonic           yes;
    nOuterCorrectors    1;
    nCorrectors         2;
    nNonOrthogonalCorrectors 1;

    maxCo               0.2;
    maxDeltaT           1e-6;
    alphaTemp           0.05;
    alphaY              0.05;
    Yref { O2 0.1; }
    rDeltaTSmoothingCoeff   0.02;
    rDeltaTDampingCoeff     1;
}

relaxationFactors { equations { ".*" 1; } }
''')
w("system/decomposeParDict", "numberOfSubdomains 4;\nmethod          scotch;\n")

# ------------------------------------------------------------------ 0.orig
def field(name, cls, dim, internal, bc):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if "(" in p else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16}{val};\n" for k, val in d.items()) + "    }\n"
    body += '    "(front|back)"\n    {\n        type            wedge;\n    }\n}\n'
    w(f"0.orig/{name}", body, cls=cls)


S, WALL = "volScalarField", "(nozzle|plate)"
ramp = f"table ((0 {PAMB}) (1500 {P0}) (1e6 {P0}))"
field("p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {
    "inlet": {"type": "uniformTotalPressure", "p0": ramp, "value": f"uniform {PAMB}"},
    "ambient": {"type": "totalPressure", "p0": f"uniform {PAMB}", "value": f"uniform {PAMB}"},
    WALL: {"type": "zeroGradient"}})
field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
    "(inlet|ambient)": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    WALL: {"type": "noSlip"}})
field("T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {
    "inlet": {"type": "totalTemperature", "gamma": "1.4", "T0": f"uniform {T0}", "value": f"uniform {T0}"},
    "ambient": {"type": "inletOutlet", "inletValue": f"uniform {T0}", "value": f"uniform {T0}"},
    WALL: {"type": "fixedValue", "value": f"uniform {T0}"}})
yin = P["YO2_in"]
for sp, a, b in (("O2", yin, YO2_AIR), ("N2", 1 - yin, YN2_AIR)):
    field(sp, S, "[0 0 0 0 0 0 0]", f"uniform {b}", {
        "inlet": {"type": "fixedValue", "value": f"uniform {a}"},
        "ambient": {"type": "inletOutlet", "inletValue": f"uniform {b}", "value": f"uniform {b}"},
        WALL: {"type": "zeroGradient"}})
field("Ydefault", S, "[0 0 0 0 0 0 0]", "uniform 0", {
    "(inlet|ambient)": {"type": "inletOutlet", "inletValue": "uniform 0", "value": "uniform 0"},
    WALL: {"type": "zeroGradient"}})
field("k", S, "[0 2 -2 0 0 0 0]", "uniform 1", {
    "inlet": {"type": "turbulentIntensityKineticEnergyInlet", "intensity": "0.05", "value": "uniform 1"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 1e-3", "value": "uniform 1"},
    WALL: {"type": "kqRWallFunction", "value": "uniform 1"}})
field("omega", S, "[0 0 -1 0 0 0 0]", "uniform 1e4", {
    "inlet": {"type": "turbulentMixingLengthFrequencyInlet", "mixingLength": "1e-4", "value": "uniform 1e4"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 100", "value": "uniform 1e4"},
    WALL: {"type": "omegaWallFunction", "value": "uniform 1e4"}})
field("nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {
    "(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "nutkWallFunction", "value": "uniform 0"}})
field("alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {
    "(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "compressible::alphatWallFunction", "value": "uniform 0"}})

with open(os.path.join(case, "Allrun"), "w") as f:
    f.write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
rm -rf 0 && cp -r 0.orig 0
runApplication blockMesh
runApplication checkMesh
runApplication $(getApplication)
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
json.dump(P, open(os.path.join(case, "caseParams.json"), "w"), indent=1, ensure_ascii=False)
open(os.path.join(case, "post.foam"), "w").close()
print(case, "cells:", sum(nr * ny for _, _, (nr, ny), _ in B), "p0 abs =", P0)
