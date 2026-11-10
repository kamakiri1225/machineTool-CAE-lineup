#!/usr/bin/env python3
"""Laser assist gas (O2/N2 mixture, supersonic) + steel plate, conjugate heat transfer, STEADY,
for the Foundation release OpenFOAM 13 (openfoam.org): foamMultiRun with
    air   : multicomponentFluid  (species O2/N2, transonic pressure, LTS 'localEuler' -> pseudo-steady)
    plate : solid                (ddtSchemes steadyState: the solid module does not take LTS, it is solved steady)
The mesh (blockMeshDict with zones air / plate) is made by tools/make_laser_cht_case.py, so the geometry and
the cutting-front heating (cells next to the kerf wall held at TFRONT) are the same as the ESI v2512 cases.

usage (OpenFOAM 13 environment):
    source /opt/openfoam13/etc/bashrc
    python3 make_laser_cht_foundation.py <caseDir> '{"YO2_in": 0.0}'     # 0 = N2 cutting, 1 = O2, 0.233 = air
    cd <caseDir> && ./Allrun
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
case = os.path.abspath(sys.argv[1])
P = dict(YO2_in=0.0, endTime=12000, ramp=1500, maxCo=0.2, nProc=4)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
geo = {k: v for k, v in P.items() if k not in ("YO2_in", "endTime", "ramp", "maxCo", "nProc")}
subprocess.run([sys.executable, os.path.join(HERE, "make_laser_cht_case.py"), case, json.dumps(geo)], check=True,
               stdout=subprocess.DEVNULL)
C = json.load(open(os.path.join(case, "caseParams.json")))
# keep only the mesh and the topoSetDict from the ESI generator; everything else is rewritten for OpenFOAM 13
for d in ("constant", "system", "0.orig"):
    for root, _, files in os.walk(os.path.join(case, d)):
        for f in files:
            if f not in ("blockMeshDict", "topoSetDict"):
                os.remove(os.path.join(root, f))
PAMB, T0 = 101325.0, 300.0
P0 = PAMB + C["p0_gauge"]
YO2_AIR = 0.233
HDR = 'FoamFile\n{{\n    format      ascii;\n    class       {cls};\n    location    "{loc}";\n    object      {obj};\n}}\n\n'


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


# ---------------------------------------------------------------- controlDict: one LTS iteration = one 'time' step
w("system/controlDict", foam({
    "application": "foamMultiRun",
    "regionSolvers": {"air": "multicomponentFluid", "plate": "solid"},
    "_": None,
    "startFrom": "latestTime", "startTime": 0, "stopAt": "endTime", "endTime": P["endTime"], "deltaT": 1,
    "writeControl": "timeStep", "writeInterval": 1000, "purgeWrite": 0, "writeFormat": "ascii", "writePrecision": 8,
    "timeFormat": "general", "timePrecision": 6, "runTimeModifiable": "true",
    "_2": None,
    "functions": {
        "massFlow": {"type": "surfaceFieldValue", "libs": '("libfieldFunctionObjects.so")', "region": "air",
                     "fields": "(phi)", "operation": "sum", "select": "patch", "patch": "inlet",
                     "writeFields": "false", "writeControl": "timeStep", "writeInterval": 20},
        "plateT": {"type": "volFieldValue", "libs": '("libfieldFunctionObjects.so")', "region": "plate",
                   "fields": "(T)", "operation": "volAverage", "cellZone": "all", "writeFields": "false",
                   "writeControl": "timeStep", "writeInterval": 20},
        "heatToGas": {"type": "wallHeatFlux", "libs": '("libfieldFunctionObjects.so")', "region": "air",
                      "patches": "(air_to_plate)", "writeControl": "writeTime"}}}))
w("system/fvSolution", foam({"PIMPLE": {"nOuterCorrectors": 1}}))
w("system/decomposeParDict", foam({"numberOfSubdomains": P["nProc"], "method": "scotch"}))

# ---------------------------------------------------------------- air: O2/N2 mixture, LTS + transonic (same as reactingFoam run002/003)
w("constant/air/physicalProperties", foam({
    "thermoType": {"type": "heRhoThermo", "mixture": "multicomponentMixture", "transport": "const", "thermo": "hConst",
                   "energy": "sensibleEnthalpy", "equationOfState": "perfectGas", "specie": "specie"},
    "_": None, "defaultSpecie": "N2"}) + "\nspecies\n(\n    O2\n    N2\n);\n\n"
  + foam({"O2": {"specie": {"molWeight": 31.9988}, "thermodynamics": {"Cp": 918, "hf": 0},
                 "transport": {"mu": 2.06e-05, "Pr": 0.72}}}) + "\n"
  + foam({"N2": {"specie": {"molWeight": 28.0134}, "thermodynamics": {"Cp": 1040, "hf": 0},
                 "transport": {"mu": 1.79e-05, "Pr": 0.71}}}))
w("constant/air/momentumTransport", foam({"simulationType": "RAS", "RAS": {"model": "kOmegaSST", "turbulence": "on",
                                                                          "printCoeffs": "on"}}))
w("constant/air/thermophysicalTransport", foam({"RAS": {"model": "unityLewisEddyDiffusivity", "Prt": 0.85}}))
w("system/air/fvSchemes", foam({
    "ddtSchemes": {"default": "localEuler"},
    "gradSchemes": {"default": "Gauss linear", "grad(U)": "cellLimited Gauss linear 1"},
    "divSchemes": {"default": "none", "div(phi,U)": "Gauss linearUpwindV grad(U)", "div(phi,Yi_h)": "Gauss upwind",
                   "div(phi,K)": "Gauss upwind", "div(phid,p)": "Gauss upwind", "div(phi,k)": "Gauss upwind",
                   "div(phi,omega)": "Gauss upwind", "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear"},
    "laplacianSchemes": {"default": "Gauss linear limited corrected 0.5"},
    "interpolationSchemes": {"default": "linear"},
    "snGradSchemes": {"default": "limited corrected 0.5"},
    "wallDist": {"method": "meshWave"}}))
w("system/air/fvSolution", foam({
    "solvers": {'"rho.*"': {"solver": "diagonal"},
                '"p.*"': {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-8", "relTol": 0.05},
                '"(U|h|k|omega|Yi).*"': {"solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": "1e-8", "relTol": 0.1}},
    "_1": None,
    "PIMPLE": {"momentumPredictor": "yes", "transonic": "yes", "nCorrectors": 2, "nNonOrthogonalCorrectors": 1,
               "maxCo": P["maxCo"], "maxDeltaT": "1e-6", "alphaTemp": 1, "rDeltaTSmoothingCoeff": 0.02,
               "rDeltaTDampingCoeff": 1},
    "_2": None,
    "relaxationFactors": {"equations": {'".*"': 1}}}))
w("system/air/fvConstraints", foam({"limitT": {"type": "limitTemperature", "select": "all", "min": 50, "max": 2500}}))

# ---------------------------------------------------------------- plate: steady conduction, cutting front held at TFRONT
w("constant/plate/physicalProperties", foam({"thermoType": "constSolidThermo",
                                             "rho": {"type": "uniform", "value": 7850},
                                             "Cv": {"type": "uniform", "value": 460},
                                             "kappa": {"type": "uniform", "value": 45}}))
w("system/plate/fvConstraints", foam({"laserFront": {"type": "fixedTemperature", "cellZone": "laserFront",
                                                     "mode": "uniform", "temperature": f"constant {C['TFRONT']}"}}))
w("system/plate/fvSchemes", foam({
    "ddtSchemes": {"default": "steadyState"}, "gradSchemes": {"default": "Gauss linear"}, "divSchemes": {"default": "none"},
    "laplacianSchemes": {"default": "Gauss linear corrected"}, "interpolationSchemes": {"default": "linear"},
    "snGradSchemes": {"default": "corrected"}}))
w("system/plate/fvSolution", foam({
    "solvers": {'"e.*"': {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-9", "relTol": 0.01}},
    "_": None, "PIMPLE": {"nNonOrthogonalCorrectors": 0},
    "_2": None, "relaxationFactors": {"equations": {'".*"': 0.9}}}))
for r in ("air", "plate"):
    w(f"system/{r}/decomposeParDict", foam({"numberOfSubdomains": P["nProc"], "method": "scotch"}))


# ---------------------------------------------------------------- 0.orig
def field(region, name, internal, bc, cls="volScalarField", dim="[0 0 0 0 0 0 0]"):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if any(c in p for c in "(|*") else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16}{v};\n" for k, v in d.items()) + "    }\n"
    body += '    "(front|back)"\n    {\n        type            wedge;\n    }\n}\n'
    w(f"0.orig/{region}/{name}", body, cls)


WALL, PR = "(nozzle|air_to_plate)", "[1 -1 -2 0 0 0 0]"
ramp = f"table ((0 {PAMB}) ({P['ramp']} {P0}) (1e9 {P0}))"
field("air", "p", f"uniform {PAMB}", {
    "inlet": {"type": "uniformTotalPressure", "p0": ramp, "value": f"uniform {PAMB}"},
    "ambient": {"type": "totalPressure", "p0": f"uniform {PAMB}", "value": f"uniform {PAMB}"},
    WALL: {"type": "zeroGradient"}}, dim=PR)
field("air", "U", "uniform (0 0 0)", {"(inlet|ambient)": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
                                     WALL: {"type": "noSlip"}}, cls="volVectorField", dim="[0 1 -1 0 0 0 0]")
field("air", "T", f"uniform {T0}", {
    "inlet": {"type": "totalTemperature", "gamma": 1.4, "T0": f"uniform {T0}", "value": f"uniform {T0}"},
    "ambient": {"type": "inletOutlet", "inletValue": f"uniform {T0}", "value": f"uniform {T0}"},
    "nozzle": {"type": "zeroGradient"},
    "air_to_plate": {"type": "coupledTemperature", "value": f"uniform {T0}"}}, dim="[0 0 0 1 0 0 0]")
for sp, a, b in (("O2", P["YO2_in"], YO2_AIR), ("N2", 1 - P["YO2_in"], 1 - YO2_AIR)):
    field("air", sp, f"uniform {b}", {"inlet": {"type": "fixedValue", "value": f"uniform {a}"},
                                      "ambient": {"type": "inletOutlet", "inletValue": f"uniform {b}", "value": f"uniform {b}"},
                                      WALL: {"type": "zeroGradient"}})
field("air", "Ydefault", "uniform 0", {"(inlet|ambient)": {"type": "inletOutlet", "inletValue": "uniform 0", "value": "uniform 0"},
                                       WALL: {"type": "zeroGradient"}})
field("air", "k", "uniform 1", {
    "inlet": {"type": "turbulentIntensityKineticEnergyInlet", "intensity": 0.05, "value": "uniform 1"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 1e-3", "value": "uniform 1"},
    WALL: {"type": "kqRWallFunction", "value": "uniform 1"}}, dim="[0 2 -2 0 0 0 0]")
field("air", "omega", "uniform 1e4", {
    "inlet": {"type": "turbulentMixingLengthFrequencyInlet", "mixingLength": "1e-4", "value": "uniform 1e4"},
    "ambient": {"type": "inletOutlet", "inletValue": "uniform 100", "value": "uniform 1e4"},
    WALL: {"type": "omegaWallFunction", "value": "uniform 1e4"}}, dim="[0 0 -1 0 0 0 0]")
field("air", "nut", "uniform 0", {"(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
                                  WALL: {"type": "nutkWallFunction", "value": "uniform 0"}}, dim="[0 2 -1 0 0 0 0]")
field("air", "alphat", "uniform 0", {"(inlet|ambient)": {"type": "calculated", "value": "uniform 0"},
                                     "nozzle": {"type": "compressible::alphatWallFunction", "value": "uniform 0"},
                                     "air_to_plate": {"type": "compressible::alphatJayatillekeWallFunction", "Prt": 0.85,
                                                      "value": "uniform 0"}}, dim="[1 -1 -1 0 0 0 0]")
field("plate", "T", f"uniform {T0}", {"plateEdge": {"type": "fixedValue", "value": f"uniform {T0}"},
                                      "plate_to_air": {"type": "coupledTemperature", "value": f"uniform {T0}"}},
      dim="[0 0 0 1 0 0 0]")

open(os.path.join(case, "Allrun"), "w").write(f'''#!/bin/sh
# OpenFOAM 13 (Foundation): source /opt/openfoam13/etc/bashrc ParaView_TYPE=none   (ParaView check hangs on WSL)
cd "${{0%/*}}" || exit
. $WM_PROJECT_DIR/bin/tools/RunFunctions
export HWLOC_COMPONENTS=-gl        # WSL2: keeps mpirun from hanging
rm -rf 0 constant/polyMesh constant/air/polyMesh constant/plate/polyMesh processor*
runApplication blockMesh
runApplication splitMeshRegions -cellZones -overwrite
runApplication topoSet -region plate
rm -rf 0 && cp -r 0.orig 0
runApplication decomposePar -allRegions -force
runParallel foamMultiRun
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
C.update(dict(solver="foamMultiRun (OpenFOAM 13): multicomponentFluid LTS + solid steadyState", YO2_in=P["YO2_in"],
              endTime=P["endTime"], ramp=P["ramp"], maxCo=P["maxCo"]))
json.dump(C, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
print("OpenFOAM 13 case ready:", case)
