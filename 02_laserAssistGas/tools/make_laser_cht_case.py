#!/usr/bin/env python3
"""Axisymmetric laser-cutting assist-gas case WITH the steel plate solved as a solid (conjugate heat transfer),
steady, chtMultiRegionSimpleFoam.  Same geometry / conditions as tools/make_laser_case.py (run003: high-pressure air).

Regions
  air   : fluid, one gas (air: assist gas = surroundings, so no species are needed), perfect gas, k-omega SST
  plate : solid steel (rho 7850, cp 460, kappa 45). The cutting front (cells next to the kerf wall, RK..RK+DFRONT)
          is held at TFRONT (about the melting point of steel) with fvOptions fixedTemperatureConstraint,
          standing in for the laser that keeps the front molten.
Notes
  chtMultiRegionSimpleFoam has no 'transonic' pressure formulation and no species; this case checks how far it
  gets for the supersonic air jet, compared with reactingFoam (run003).

usage: python3 make_laser_cht_case.py <caseDir> [json-params]
"""
import json, math, os, sys

case = os.path.abspath(sys.argv[1])
P = dict(RK=0.15, RN=1.0, RIN=2.0, RNOZ=3.0, ROUT=8.0, H=0.8, LN=3.0, T=2.0, LB=4.0, p0_gauge=0.8e6,
         TFRONT=1800.0, DFRONT=0.04, endTime=8000, ramp=2000)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
PAMB, T0 = 101325.0, 300.0
P0 = PAMB + P["p0_gauge"]
AIR = dict(W=28.85, Cp=1011.6, mu=1.84e-5, Pr=0.71)          # same mixture as run003 (O2 0.233 / N2 0.767)
STEEL = dict(rho=7850, Cp=460, kappa=45)
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    location    "{loc}";\n    object      {obj};\n}}\n\n'


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w") as f:
        f.write(HDR.format(cls=cls, loc=os.path.dirname(path), obj=os.path.basename(path)) + body)


def foam(d, ind=0):
    """python dict -> OpenFOAM dictionary text (one entry per line, indented)."""
    pad, out = " " * ind, ""
    for k, v in d.items():
        if isinstance(v, dict):
            out += f"{pad}{k}\n{pad}{{\n{foam(v, ind + 4)}{pad}}}\n"
        elif v is None:
            out += "\n"
        else:
            out += f"{pad}{k:<{max(16, len(k) + 1)}}{v};\n"
    return out


# ------------------------------------------------------------------ mesh (blocks carry the zone name air / plate)
RK, RN, RIN, RNOZ, ROUT = P["RK"], P["RN"], P["RIN"], P["RNOZ"], P["ROUT"]
H, LN, T, LB = P["H"], P["LN"], P["T"], P["LB"]
Y0, Y1, Y2, Y3, Y4 = -T - LB, -T, 0.0, H, H + LN
RK4 = RK * RIN / RN
A = math.radians(2.5)
verts, vidx = [], {}


def V(r, y, s):
    key = (round(r, 9), round(y, 9), 0 if r == 0 else s)
    if key not in vidx:
        vidx[key] = len(verts)
        verts.append((r * math.cos(A), y, s * r * math.sin(A)))
    return vidx[key]


def hexq(r00, r10, r11, r01, ya, yb):
    return [V(r00, ya, -1), V(r10, ya, -1), V(r11, yb, -1), V(r01, yb, -1),
            V(r00, ya, 1), V(r10, ya, 1), V(r11, yb, 1), V(r01, yb, 1)]


def n_(L, d):
    return max(4, int(round(L / d)))


dr_k = RK / 10
cells = {"rk": 10, "rn": n_(RN - RK, dr_k * 1.2), "rnoz": n_(RNOZ - RN, 0.06), "rout": n_(ROUT - RNOZ, 0.25),
         "below": n_(LB, 0.12), "kerf": n_(T, 0.04), "gap": n_(H, 0.025), "noz": n_(LN, 0.05)}
B = []
radial = [(0, RK, "rk", "1"), (RK, RN, "rn", "1"), (RN, RNOZ, "rnoz", "3"), (RNOZ, ROUT, "rout", "4")]
for r0, r1, key, g in radial:
    B.append((f"below_{key}", "air", hexq(r0, r1, r1, r0, Y0, Y1), (cells[key], cells["below"]), f"{g} 0.25 1"))
B.append(("kerf", "air", hexq(0, RK, RK, 0, Y1, Y2), (cells["rk"], cells["kerf"]), "1 1 1"))
for r0, r1, key, g in radial[1:]:                         # the plate itself (solid), r = RK..ROUT, y = -T..0
    B.append((f"plate_{key}", "plate", hexq(r0, r1, r1, r0, Y1, Y2), (cells[key], cells["kerf"]), f"{g} 1 1"))
for r0, r1, key, g in radial:
    B.append((f"gap_{key}", "air", hexq(r0, r1, r1, r0, Y2, Y3), (cells[key], cells["gap"]), f"{g} 1 1"))
B.append(("noz_rk", "air", hexq(0, RK, RK4, 0, Y3, Y4), (cells["rk"], cells["noz"]), "1 1 1"))
B.append(("noz_rn", "air", hexq(RK, RN, RIN, RK4, Y3, Y4), (cells["rn"], cells["noz"]), "1 1 1"))
B.append(("above", "air", hexq(RNOZ, ROUT, ROUT, RNOZ, Y3, Y4), (cells["rout"], cells["noz"]), "4 2 1"))
H_ = {n: h for n, _, h, _, _ in B}
FACE = {"ymin": (0, 4, 5, 1), "ymax": (3, 2, 6, 7), "rmin": (0, 3, 7, 4), "rmax": (1, 5, 6, 2),
        "front": (0, 1, 2, 3), "back": (4, 7, 6, 5)}


def F(name, side):
    return "(" + " ".join(str(H_[name][i]) for i in FACE[side]) + ")"


patches = {
    "inlet": ("patch", [F("noz_rk", "ymax"), F("noz_rn", "ymax")]),
    "nozzle": ("wall", [F("noz_rn", "rmax"), F("gap_rnoz", "ymax"), F("above", "rmin")]),
    "ambient": ("patch", [F("gap_rout", "rmax"), F("above", "rmax"), F("above", "ymax"),
                          F("below_rout", "rmax")] + [F(f"below_{k}", "ymin") for k in ("rk", "rn", "rnoz", "rout")]),
    "plateEdge": ("wall", [F("plate_rout", "rmax")]),      # far edge of the plate (r = ROUT), held at 300 K
    "front": ("wedge", [F(n, "front") for n in H_]),
    "back": ("wedge", [F(n, "back") for n in H_]),
}
blk = "\n".join(f"    hex ({' '.join(map(str, h))}) {z} ({nr} {ny} 1) simpleGrading ({g})" for _, z, h, (nr, ny), g in B)
pat = "\n".join(f"    {n}\n    {{\n        type            {t};\n        faces\n        (\n"
                + "\n".join("            " + x for x in fs) + "\n        );\n    }" for n, (t, fs) in patches.items())
w("system/blockMeshDict", "scale           1e-3;\n\nvertices\n(\n" + "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in verts)
  + f"\n);\n\nblocks\n(\n{blk}\n);\n\nedges\n(\n);\n\nboundary\n(\n{pat}\n);\n")

# ------------------------------------------------------------------ constant
w("constant/regionProperties", "regions\n(\n    fluid       (air)\n    solid       (plate)\n);\n")
w("constant/g", "dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 0 0);    // gravity is negligible for this jet\n",
  cls="uniformDimensionedVectorField")
w("constant/air/thermophysicalProperties", foam({
    "thermoType": {"type": "heRhoThermo", "mixture": "pureMixture", "transport": "const", "thermo": "hConst",
                   "equationOfState": "perfectGas", "specie": "specie", "energy": "sensibleEnthalpy"},
    "_": None,
    "mixture": {"specie": {"molWeight": AIR["W"]},
                "thermodynamics": {"Cp": AIR["Cp"], "Hf": 0},
                "transport": {"mu": AIR["mu"], "Pr": AIR["Pr"]}}}))
w("constant/air/turbulenceProperties", foam({"simulationType": "RAS",
                                             "RAS": {"RASModel": "kOmegaSST", "turbulence": "on", "printCoeffs": "on"}}))
w("constant/air/radiationProperties", foam({"radiation": "off", "radiationModel": "none"}))
w("constant/air/fvOptions", foam({"limitT": {"type": "limitTemperature", "active": "yes", "selectionMode": "all",
                                            "min": 80, "max": 2500}}))   # steady start-up guard (no transonic formulation)
w("constant/plate/thermophysicalProperties", foam({
    "thermoType": {"type": "heSolidThermo", "mixture": "pureMixture", "transport": "constIso", "thermo": "hConst",
                   "equationOfState": "rhoConst", "specie": "specie", "energy": "sensibleEnthalpy"},
    "_": None,
    "mixture": {"specie": {"molWeight": 55.85},
                "transport": {"kappa": STEEL["kappa"]},
                "thermodynamics": {"Hf": 0, "Cp": STEEL["Cp"]},
                "equationOfState": {"rho": STEEL["rho"]}}}))
w("constant/plate/radiationProperties", foam({"radiation": "off", "radiationModel": "none"}))
w("constant/plate/fvOptions", foam({"laserFront": {
    "type": "fixedTemperatureConstraint", "active": "yes", "selectionMode": "cellZone", "cellZone": "laserFront",
    "mode": "uniform", "temperature": f"constant {P['TFRONT']}"}}))

# ------------------------------------------------------------------ system
w("system/controlDict", foam({
    "application": "chtMultiRegionSimpleFoam", "startFrom": "latestTime", "startTime": 0, "stopAt": "endTime",
    "endTime": P["endTime"], "deltaT": 1, "writeControl": "timeStep", "writeInterval": 1000, "purgeWrite": 0,
    "writeFormat": "ascii", "writePrecision": 8, "timeFormat": "general", "timePrecision": 6,
    "runTimeModifiable": "true", "_": None,
    "functions": {
        "massFlow": {"type": "surfaceFieldValue", "libs": "(fieldFunctionObjects)", "region": "air",
                     "regionType": "patch", "name": "inlet", "operation": "sum", "fields": "(phi)",
                     "writeFields": "false", "writeControl": "timeStep", "writeInterval": 20},
        "heatToGas": {"type": "wallHeatFlux", "libs": "(fieldFunctionObjects)", "region": "air",
                      "patches": "(air_to_plate)", "writeControl": "writeTime"},
        "MachNo": {"type": "MachNo", "libs": "(fieldFunctionObjects)", "region": "air",
                   "executeControl": "writeTime", "writeControl": "writeTime"}}}))
w("system/fvSchemes", "".join(f"{k}\n{{\n}}\n\n" for k in ("ddtSchemes", "gradSchemes", "divSchemes", "laplacianSchemes", "interpolationSchemes", "snGradSchemes")))
w("system/fvSolution", "")
w("system/decomposeParDict", foam({"numberOfSubdomains": 4, "method": "scotch"}))
w("system/plate/topoSetDict", f"""actions
(
    {{
        name    laserFrontSet;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     (0 {-T * 1e-3 - 1e-6} -1) ({(RK + P['DFRONT']) * 1e-3} 1e-6 1);   // plate cells next to the kerf wall
    }}
    {{
        name    laserFront;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     laserFrontSet;
    }}
);
""")
w("system/air/fvSchemes", foam({
    "ddtSchemes": {"default": "steadyState"},
    "gradSchemes": {"default": "Gauss linear", "grad(U)": "cellLimited Gauss linear 1"},
    "divSchemes": {"default": "none", "div(phi,U)": "bounded Gauss upwind", "div(phi,K)": "bounded Gauss upwind",
                   "div(phi,h)": "bounded Gauss upwind", "div(phi,k)": "bounded Gauss upwind",
                   "div(phi,omega)": "bounded Gauss upwind",
                   "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear"},
    "laplacianSchemes": {"default": "Gauss linear limited corrected 0.5"},
    "interpolationSchemes": {"default": "linear"},
    "snGradSchemes": {"default": "limited corrected 0.5"},
    "wallDist": {"method": "meshWave"}}))
w("system/air/fvSolution", foam({
    "solvers": {"rho": {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-7", "relTol": 0},
                "p_rgh": {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-8", "relTol": 0.01},
                '"(U|h|k|omega)"': {"solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": "1e-8", "relTol": 0.1}},
    "_1": None,
    "SIMPLE": {"momentumPredictor": "yes", "nNonOrthogonalCorrectors": 1, "rhoMin": 0.05, "rhoMax": 30},
    "_2": None,
    "relaxationFactors": {"fields": {"rho": 0.5, "p_rgh": 0.1},
                          "equations": {"U": 0.3, "h": 0.3, '"(k|omega)"': 0.3}}}))
w("system/plate/fvSchemes", foam({
    "ddtSchemes": {"default": "steadyState"}, "gradSchemes": {"default": "Gauss linear"},
    "divSchemes": {"default": "none"}, "laplacianSchemes": {"default": "none", "laplacian(alpha,h)": "Gauss linear corrected"},
    "interpolationSchemes": {"default": "linear"}, "snGradSchemes": {"default": "corrected"}}))
w("system/plate/fvSolution", foam({
    "solvers": {"h": {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-8", "relTol": 0.1}},
    "_1": None, "SIMPLE": {"nNonOrthogonalCorrectors": 0}, "_2": None,
    "relaxationFactors": {"equations": {"h": 0.9}}}))
for r in ("air", "plate"):
    w(f"system/{r}/decomposeParDict", foam({"numberOfSubdomains": 4, "method": "scotch"}))


# ------------------------------------------------------------------ 0.orig (patch names after splitMeshRegions)
def field(region, name, cls, dim, internal, bc):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if "(" in p or "*" in p else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16}{v};\n" for k, v in d.items()) + "    }\n"
    body += '    "(front|back)"\n    {\n        type            wedge;\n    }\n}\n'
    w(f"0.orig/{region}/{name}", body, cls)


S = "volScalarField"
ramp = f"table ((0 {PAMB}) ({P['ramp']} {P0}) (1e6 {P0}))"
WALL = "(nozzle|air_to_plate)"
field("air", "p_rgh", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {
    "inlet": {"type": "uniformTotalPressure", "p0": ramp, "value": f"uniform {PAMB}"},
    "ambient": {"type": "totalPressure", "p0": f"uniform {PAMB}", "value": f"uniform {PAMB}"},
    WALL: {"type": "fixedFluxPressure", "value": f"uniform {PAMB}"}})
field("air", "p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {".*": {"type": "calculated", "value": f"uniform {PAMB}"}})
field("air", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
    "(inlet|ambient)": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    WALL: {"type": "noSlip"}})
field("air", "T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {
    "inlet": {"type": "totalTemperature", "gamma": 1.4, "T0": f"uniform {T0}", "value": f"uniform {T0}"},
    "ambient": {"type": "inletOutlet", "inletValue": f"uniform {T0}", "value": f"uniform {T0}"},
    "nozzle": {"type": "zeroGradient"},
    "air_to_plate": {"type": "compressible::turbulentTemperatureCoupledBaffleMixed", "Tnbr": "T",
                     "kappaMethod": "fluidThermo", "value": f"uniform {T0}"}})
field("air", "k", S, "[0 2 -2 0 0 0 0]", "uniform 1", {
    "inlet": {"type": "turbulentIntensityKineticEnergyInlet", "intensity": 0.05, "value": "uniform 1"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 1e-3", "value": "uniform 1"},
    WALL: {"type": "kqRWallFunction", "value": "uniform 1"}})
field("air", "omega", S, "[0 0 -1 0 0 0 0]", "uniform 1e4", {
    "inlet": {"type": "turbulentMixingLengthFrequencyInlet", "mixingLength": "1e-4", "value": "uniform 1e4"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 100", "value": "uniform 1e4"},
    WALL: {"type": "omegaWallFunction", "value": "uniform 1e4"}})
field("air", "nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {
    "(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "nutkWallFunction", "value": "uniform 0"}})
field("air", "alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {
    "(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "compressible::alphatJayatillekeWallFunction", "Prt": 0.85, "value": "uniform 0"}})
field("plate", "T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {
    "plateEdge": {"type": "fixedValue", "value": f"uniform {T0}"},
    "plate_to_air": {"type": "compressible::turbulentTemperatureCoupledBaffleMixed", "Tnbr": "T",
                     "kappaMethod": "solidThermo", "value": f"uniform {T0}"}})
field("plate", "p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {".*": {"type": "calculated", "value": f"uniform {PAMB}"}})

with open(os.path.join(case, "Allrun"), "w") as f:
    f.write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
export HWLOC_COMPONENTS=-gl        # WSL2: keeps mpirun from hanging
rm -rf 0 constant/polyMesh constant/air/polyMesh constant/plate/polyMesh
runApplication blockMesh
runApplication splitMeshRegions -cellZones -overwrite
runApplication topoSet -region plate
cp -r 0.orig/* 0/ 2>/dev/null || cp -r 0.orig 0
runApplication checkMesh -region air
runApplication decomposePar -allRegions -force
runParallel $(getApplication)
runApplication reconstructPar -allRegions -latestTime
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
json.dump(P, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
print(case, "cells:", sum(nr * ny for _, _, _, (nr, ny), _ in B), "p0 abs =", P0)
