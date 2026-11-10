#!/usr/bin/env python3
"""Pyrotechnic inflator tank test WITH the steel tank wall and the outside air (conjugate heat transfer),
transient, chtMultiRegionFoam.  Same inflator, gas generation and tank as tools/make_tank_case.py (run001).

Regions (axisymmetric 5 deg wedge, axis = y) [mm]
  gas     : fluid  chamber r 0-15 y -75..-5, orifice r 0-5.5 y -5..0, tank r 0-170 y 0..660
            mixture N2/O2/CO2/H2O (constant cp), k-omega SST, propellant gas enters through the chamber bottom (gasInlet, mass flow table)
  wall    : solid  steel tank wall, WT thick: bottom plate (y -WT..0, r 5.5..170+WT), side (r 170..170+WT),
            lid (y 660..660+WT)
  outside : fluid  air around the tank, r 170+WT..170+WT+DO and above the lid (laminar, gravity on)
Heat transfer between the gas, the wall and the outside air is solved, so no heat transfer coefficient is assumed.
The inflator housing (chamber wall) is adiabatic as in run001.

usage: python3 make_tank_cht_case.py <caseDir> [json]
"""
import json, math, os, sys
import numpy as np

case = os.path.abspath(sys.argv[1])
P = dict(M_gas=0.025, tau=0.006, Tc=1000.0, YCO2=0.20, YH2O=0.10, tEnd=0.05, WT=5.0, DO=60.0, writeInterval=0.0005)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
WT, DO = P["WT"], P["DO"]
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    location    "{loc}";\n    object      {obj};\n}}\n\n'
CONST = {"N2": (28.0134, 1170, 3.2e-05, 0.72), "O2": (31.9988, 1090, 3.6e-05, 0.72),
         "CO2": (44.0095, 1230, 3.0e-05, 0.72), "H2O": (18.0153, 2290, 2.5e-05, 0.9)}
STEEL = dict(rho=7850, Cp=460, kappa=45)
T0, PAMB = 293.15, 101325.0


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path); os.makedirs(os.path.dirname(full), exist_ok=True)
    open(full, "w").write(HDR.format(cls=cls, loc=os.path.dirname(path), obj=os.path.basename(path)) + body)


def foam(d, ind=0):
    pad, out = " " * ind, ""
    for k, v in d.items():
        if isinstance(v, dict):
            out += f"{pad}{k}\n{pad}{{\n{foam(v, ind + 4)}{pad}}}\n"
        elif v is None:
            out += "\n"
        else:
            out += f"{pad}{k:<{max(16, len(k) + 1)}}{v};\n"
    return out


# ---------------------------------------------------------------- mesh
R = [0, 5.5, 15, 170, 170 + WT, 170 + WT + DO]
Y = [-75, -5 - 0.0, -WT if WT < 5 else -5, 0, 60, 660, 660 + WT, 660 + WT + DO]
Y = [-75, -WT, 0, 60, 660, 660 + WT, 660 + WT + DO]          # chamber top = bottom plate underside (y = -WT)
nR = [10, 10, 60, 6, 20]
gR = [1, 1, 1, 1, 4]
nY = [50, 10, 40, 80, 6, 15]
gY = [1, 1, 1, 6, 1, 4]
Z = {}   # (i, j) -> region
for i in (0, 1):
    Z[(i, 0)] = "gas"                                              # chamber
Z[(0, 1)] = "gas"                                                  # orifice
for i in (1, 2, 3):
    Z[(i, 1)] = "wall"                                             # bottom plate
for i in (0, 1, 2):
    Z[(i, 2)] = Z[(i, 3)] = "gas"                                  # tank
Z[(3, 2)] = Z[(3, 3)] = "wall"                                     # side wall
for i in (0, 1, 2, 3):
    Z[(i, 4)] = "wall"                                             # lid
    Z[(i, 5)] = "outside"                                          # air above the lid
for j in (1, 2, 3, 4, 5):
    Z[(4, j)] = "outside"                                          # air around the side
A = math.radians(2.5)
verts, vidx = [], {}


def V(i, j, s):
    r = R[i]; key = (i, j, 0 if r == 0 else s)
    if key not in vidx:
        vidx[key] = len(verts); verts.append((r * math.cos(A), Y[j], s * r * math.sin(A)))
    return vidx[key]


H = {k: [V(k[0], k[1], -1), V(k[0] + 1, k[1], -1), V(k[0] + 1, k[1] + 1, -1), V(k[0], k[1] + 1, -1),
         V(k[0], k[1], 1), V(k[0] + 1, k[1], 1), V(k[0] + 1, k[1] + 1, 1), V(k[0], k[1] + 1, 1)] for k in Z}
FACE = {"ymin": (0, 4, 5, 1), "ymax": (3, 2, 6, 7), "rmin": (0, 3, 7, 4), "rmax": (1, 5, 6, 2),
        "front": (0, 1, 2, 3), "back": (4, 7, 6, 5)}


def F(k, side):
    return "(" + " ".join(str(H[k][i]) for i in FACE[side]) + ")"


patches = {
    "gasInlet": ("patch", [F((0, 0), "ymin"), F((1, 0), "ymin")]),                       # gas from the propellant enters here
    "inflator": ("wall", [F((1, 0), "rmax")]),                                            # inflator housing, adiabatic
    "wallBottom": ("wall", [F((2, 1), "ymin"), F((3, 1), "ymin")]),                       # tank underside (on a fixture), adiabatic
    "farField": ("wall", [F((4, j), "rmax") for j in range(1, 6)] + [F((i, 5), "ymax") for i in range(5)]
                 + [F((4, 1), "ymin")]),                                                   # outer edge of the air, 293 K
    "front": ("wedge", [F(k, "front") for k in Z]),
    "back": ("wedge", [F(k, "back") for k in Z]),
}
blk = "\n".join(f"    hex ({' '.join(map(str, H[k]))}) {Z[k]} ({nR[k[0]]} {nY[k[1]]} 1) "
                f"simpleGrading ({gR[k[0]]} {gY[k[1]]} 1)" for k in Z)
pat = "\n".join(f"    {n}\n    {{\n        type            {t};\n        faces\n        (\n"
                + "\n".join("            " + x for x in fs) + "\n        );\n    }" for n, (t, fs) in patches.items())
w("system/blockMeshDict", "scale           1e-3;\n\nvertices\n(\n" + "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in verts)
  + f"\n);\n\nblocks\n(\n{blk}\n);\n\nedges\n(\n);\n\nboundary\n(\n{pat}\n);\n")

# ---------------------------------------------------------------- constant
w("constant/regionProperties", "regions\n(\n    fluid       (gas outside)\n    solid       (wall)\n);\n")
w("constant/g", "dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 -9.81 0);\n", cls="uniformDimensionedVectorField")
species = foam({"species": None}).strip()
for r in ("gas", "outside"):
    w(f"constant/{r}/thermophysicalProperties", foam({
        "thermoType": {"type": "heRhoThermo", "mixture": "reactingMixture", "transport": "const", "thermo": "hConst",
                       "energy": "sensibleEnthalpy", "equationOfState": "perfectGas", "specie": "specie"},
        "_": None, "inertSpecie": "N2", "chemistryReader": "foamChemistryReader",
        "foamChemistryFile": f'"<constant>/{r}/reactions"', "foamChemistryThermoFile": f'"<constant>/{r}/thermo.species"'}))
    w(f"constant/{r}/reactions", "species\n(\n    N2\n    O2\n    CO2\n    H2O\n);\n\nreactions\n{\n}\n")
    w(f"constant/{r}/thermo.species", "\n".join(foam({k: {
        "specie": {"molWeight": v[0]}, "thermodynamics": {"Cp": v[1], "Hf": 0},
        "transport": {"mu": v[2], "Pr": v[3]}}}) for k, v in CONST.items()))
    w(f"constant/{r}/combustionProperties", foam({"combustionModel": "none"}))
    w(f"constant/{r}/chemistryProperties", foam({"chemistryType": {"solver": "none"}, "chemistry": "off",
                                                 "initialChemicalTimeStep": "1e-07"}))
    w(f"constant/{r}/radiationProperties", foam({"radiation": "off", "radiationModel": "none"}))
w("constant/gas/turbulenceProperties", foam({"simulationType": "RAS",
                                             "RAS": {"RASModel": "kOmegaSST", "turbulence": "on", "printCoeffs": "on"}}))
w("constant/outside/turbulenceProperties", foam({"simulationType": "laminar"}))
w("constant/wall/thermophysicalProperties", foam({
    "thermoType": {"type": "heSolidThermo", "mixture": "pureMixture", "transport": "constIso", "thermo": "hConst",
                   "equationOfState": "rhoConst", "specie": "specie", "energy": "sensibleEnthalpy"},
    "_": None,
    "mixture": {"specie": {"molWeight": 55.85}, "transport": {"kappa": STEEL["kappa"]},
                "thermodynamics": {"Hf": 0, "Cp": STEEL["Cp"]}, "equationOfState": {"rho": STEEL["rho"]}}}))
w("constant/wall/radiationProperties", foam({"radiation": "off", "radiationModel": "none"}))

# gas generation in the chamber (same as run001)
Vc = math.pi * (15e-3) ** 2 * 70e-3
tt = np.linspace(0, P["tEnd"], 401)
md = P["M_gas"] * tt / P["tau"] ** 2 * np.exp(-tt / P["tau"]) / Vc
hs = sum(Yk * CONST[k][1] for k, Yk in (("N2", 1 - P["YCO2"] - P["YH2O"]), ("CO2", P["YCO2"]), ("H2O", P["YH2O"]))) * (P["Tc"] - 298.15)
tab = lambda f: "table\n                (\n" + "\n".join(f"                    ({a:.6g} {b:.6g})" for a, b in zip(tt, f)) + "\n                )"
MDOT_PATCH = md * Vc * 2.5 * 2 / 360.0      # kg/s through the 5-degree wedge patch
mtab = "table\n        (\n" + "\n".join(f"            ({a:.6g} {b:.6g})" for a, b in zip(tt, MDOT_PATCH)) + "\n        )"
w("constant/gas/fvOptions", """limitT
{
    type            limitTemperature;
    selectionMode   all;
    min             150;
    max             3000;
}
""")


# ---------------------------------------------------------------- system
w("system/controlDict", foam({
    "application": "chtMultiRegionFoam", "startFrom": "latestTime", "startTime": 0, "stopAt": "endTime",
    "endTime": P["tEnd"], "deltaT": "1e-08", "writeControl": "adjustable", "writeInterval": P["writeInterval"],
    "purgeWrite": 0, "writeFormat": "binary", "writePrecision": 8, "timeFormat": "general", "timePrecision": 6,
    "runTimeModifiable": "true", "adjustTimeStep": "yes", "maxCo": 0.5, "maxDi": 10.0, "maxDeltaT": "1e-05",
    "_": None,
    "functions": {
        "tankAvg": {"type": "volFieldValue", "libs": "(fieldFunctionObjects)", "region": "gas", "regionType": "cellZone",
                    "name": "tank", "operation": "volAverage", "fields": "(p T CO2)", "writeFields": "false",
                    "writeControl": "timeStep", "writeInterval": 20},
        "chamberAvg": {"type": "volFieldValue", "libs": "(fieldFunctionObjects)", "region": "gas", "regionType": "cellZone",
                       "name": "chamber", "operation": "volAverage", "fields": "(p T)", "writeFields": "false",
                       "writeControl": "timeStep", "writeInterval": 20},
        "wallAvg": {"type": "volFieldValue", "libs": "(fieldFunctionObjects)", "region": "wall", "regionType": "all",
                    "operation": "volAverage", "fields": "(T)", "writeFields": "false",
                    "writeControl": "timeStep", "writeInterval": 20},
        "heatToWall": {"type": "wallHeatFlux", "libs": "(fieldFunctionObjects)", "region": "gas",
                       "patches": "(gas_to_wall)", "writeControl": "writeTime"},
        "innerWallT": {"type": "surfaceFieldValue", "libs": "(fieldFunctionObjects)", "region": "gas",
                       "regionType": "patch", "name": "gas_to_wall", "operation": "areaAverage", "fields": "(T)",
                       "writeFields": "false", "writeControl": "timeStep", "writeInterval": 20},
        "MachNo": {"type": "MachNo", "libs": "(fieldFunctionObjects)", "region": "gas",
                   "executeControl": "writeTime", "writeControl": "writeTime"}}}))
w("system/fvSchemes", "".join(f"{k}\n{{\n}}\n\n" for k in ("ddtSchemes", "gradSchemes", "divSchemes", "laplacianSchemes",
                                                           "interpolationSchemes", "snGradSchemes")))
w("system/fvSolution", foam({"PIMPLE": {"nOuterCorrectors": 2}}))
w("system/decomposeParDict", foam({"numberOfSubdomains": 8, "method": "scotch"}))
fluidSchemes = {
    "ddtSchemes": {"default": "Euler"},
    "gradSchemes": {"default": "Gauss linear", "grad(U)": "cellLimited Gauss linear 1"},
    "divSchemes": {"default": "none", "div(phi,U)": "Gauss linearUpwindV grad(U)", "div(phi,Yi_h)": "Gauss upwind",
                   "div(phi,K)": "Gauss upwind", "div(phi,h)": "Gauss upwind", "div(phi,k)": "Gauss upwind",
                   "div(phi,omega)": "Gauss upwind", "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear"},
    "laplacianSchemes": {"default": "Gauss linear corrected"},
    "interpolationSchemes": {"default": "linear"},
    "snGradSchemes": {"default": "corrected"},
    "wallDist": {"method": "meshWave"}}
fluidSolution = {
    "solvers": {'"rho.*"': {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-8", "relTol": 0},
                "p_rgh": {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-10", "relTol": 0.01, "maxIter": 200},
                "p_rghFinal": {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-10", "relTol": 0, "maxIter": 200},
                '"(U|h|k|omega|Yi).*"': {"solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": "1e-9", "relTol": 0}},
    "_1": None,
    "PIMPLE": {"momentumPredictor": "yes", "nCorrectors": 2, "nNonOrthogonalCorrectors": 0,
               "pRefPoint": "REFPOINT", "pRefValue": PAMB},
    "_2": None,
    "relaxationFactors": {"equations": {'".*"': 1}}}
for r in ("gas", "outside"):
    w(f"system/{r}/fvSchemes", foam(fluidSchemes))
    ref = "(0.01 0.6 0)" if r == "gas" else f"({(170 + WT + DO / 2) * 1e-3} 0.3 0)"   # a point inside each closed region
    w(f"system/{r}/fvSolution", foam(fluidSolution).replace("REFPOINT", ref))
w("system/wall/fvSchemes", foam({
    "ddtSchemes": {"default": "Euler"}, "gradSchemes": {"default": "Gauss linear"}, "divSchemes": {"default": "none"},
    "laplacianSchemes": {"default": "none", "laplacian(alpha,h)": "Gauss linear corrected"},
    "interpolationSchemes": {"default": "linear"}, "snGradSchemes": {"default": "corrected"}}))
w("system/wall/fvSolution", foam({
    "solvers": {'"h.*"': {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-9", "relTol": 0}},
    "_": None, "PIMPLE": {"nNonOrthogonalCorrectors": 0}}))
w("system/gas/topoSetDict", """actions
(
    {
        name    chamberSet;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     (-1 -1 -1) (1 -0.0049 1);      // y < -4.9 mm: the combustion chamber
    }
    {
        name    chamber;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     chamberSet;
    }
    {
        name    tankSet;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     (-1 0 -1) (1 1 1);              // y > 0: the tank
    }
    {
        name    tank;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     tankSet;
    }
);
""")
for r in ("gas", "outside", "wall"):
    w(f"system/{r}/decomposeParDict", foam({"numberOfSubdomains": 8, "method": "scotch"}))


# ---------------------------------------------------------------- 0.orig
def field(region, name, cls, dim, internal, bc):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if any(c in p for c in "(|*.") else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16}{v};\n" for k, v in d.items()) + "    }\n"
    body += '    "(front|back)"\n    {\n        type            wedge;\n    }\n}\n'
    w(f"0.orig/{region}/{name}", body, cls)


S = "volScalarField"
coupled = lambda km: {"type": "compressible::turbulentTemperatureCoupledBaffleMixed", "Tnbr": "T",
                      "kappaMethod": km, "value": f"uniform {T0}"}
YG = {"N2": 1 - P["YCO2"] - P["YH2O"], "O2": 0, "CO2": P["YCO2"], "H2O": P["YH2O"], "Ydefault": 0}
for r in ("gas", "outside"):
    inl = {"gasInlet": None} if r == "gas" else {}
    field(r, "p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {".*": {"type": "calculated", "value": f"uniform {PAMB}"}})
    field(r, "p_rgh", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {".*": {"type": "fixedFluxPressure", "value": f"uniform {PAMB}"}})
    field(r, "U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", dict(
        ({"gasInlet": {"type": "flowRateInletVelocity", "massFlowRate": mtab, "rho": "rho", "rhoInlet": 1.0,
                       "value": "uniform (0 0 0)"}} if r == "gas" else {}), **{".*": {"type": "noSlip"}}))
    for sp, v in (("N2", 0.767), ("O2", 0.233), ("CO2", 0), ("H2O", 0), ("Ydefault", 0)):
        field(r, sp, S, "[0 0 0 0 0 0 0]", f"uniform {v}", dict(
            ({"gasInlet": {"type": "fixedValue", "value": f"uniform {YG[sp]}"}} if r == "gas" else {}),
            **{".*": {"type": "zeroGradient"}}))
field("gas", "T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {"gasInlet": {"type": "fixedValue", "value": f"uniform {P['Tc']}"},
                                                          "inflator": {"type": "zeroGradient"},
                                                          "gas_to_wall": coupled("fluidThermo")})
field("outside", "T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {"farField": {"type": "fixedValue", "value": f"uniform {T0}"},
                                                              "outside_to_wall": coupled("fluidThermo")})
field("gas", "k", S, "[0 2 -2 0 0 0 0]", "uniform 1e-4", {"gasInlet": {"type": "fixedValue", "value": "uniform 1"},
                                                           ".*": {"type": "kqRWallFunction", "value": "uniform 1e-4"}})
field("gas", "omega", S, "[0 0 -1 0 0 0 0]", "uniform 10", {"gasInlet": {"type": "fixedValue", "value": "uniform 1e4"},
                                                            ".*": {"type": "omegaWallFunction", "value": "uniform 10"}})
field("gas", "nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {"gasInlet": {"type": "calculated", "value": "uniform 0"},
                                                         ".*": {"type": "nutkWallFunction", "value": "uniform 0"}})
field("gas", "alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {
    "gasInlet": {"type": "calculated", "value": "uniform 0"},
    "inflator": {"type": "compressible::alphatWallFunction", "value": "uniform 0"},
    "gas_to_wall": {"type": "compressible::alphatJayatillekeWallFunction", "Prt": 0.85, "value": "uniform 0"}})
field("outside", "alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {".*": {"type": "calculated", "value": "uniform 0"}})
field("wall", "T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {"wallBottom": {"type": "zeroGradient"},
                                                           "wall_to_gas": coupled("solidThermo"),
                                                           "wall_to_outside": coupled("solidThermo")})
field("wall", "p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {".*": {"type": "calculated", "value": f"uniform {PAMB}"}})

open(os.path.join(case, "Allrun"), "w").write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
export HWLOC_COMPONENTS=-gl        # WSL2: keeps mpirun from hanging
rm -rf 0 constant/polyMesh constant/*/polyMesh
runApplication blockMesh
runApplication splitMeshRegions -cellZones -overwrite
runApplication topoSet -region gas
rm -rf 0 && cp -r 0.orig 0
runApplication checkMesh -region gas
runApplication decomposePar -allRegions -force
runParallel $(getApplication)
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
json.dump(dict(P, Vc=Vc), open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
print("ok", case, "cells", sum(nR[i] * nY[j] for i, j in Z))
