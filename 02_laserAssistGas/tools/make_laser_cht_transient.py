#!/usr/bin/env python3
"""Laser assist gas + steel plate (conjugate heat transfer) with chtMultiRegionFoam, run in time until the flow and
the plate temperature are steady.  Builds on tools/make_laser_cht_case.py (same mesh, regions, boundary names) and
replaces what differs:

  * chtMultiRegionSimpleFoam (steady) failed when the nozzle started to choke (no transonic pressure formulation),
    so the transient chtMultiRegionFoam is used: its pressure equation keeps the compressibility term psi*dp/dt.
  * gas = O2/N2 mixture (rhoReactionThermo, species transported), the same as the reactingFoam case run003.
  * plate density is divided by RHOSCALE so that the plate reaches its steady temperature in about a millisecond
    (the steady temperature field does not depend on rho*cp; only the time to get there does).

usage: python3 make_laser_cht_transient.py <caseDir> [json-params]
       params as make_laser_cht_case.py, plus endTime [s] (default 1.5e-3), ramp [s] (2e-4), RHOSCALE (1e4)
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
case = os.path.abspath(sys.argv[1])
P = dict(endTime=1.8e-3, ramp=5e-4, RHOSCALE=1e4, YO2_in=0.233, maxCo=0.3, nOuter=2)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
base = {k: v for k, v in P.items() if k not in ("endTime", "ramp", "RHOSCALE", "YO2_in", "maxCo", "nOuter")}
subprocess.run([sys.executable, os.path.join(HERE, "make_laser_cht_case.py"), case, json.dumps(base)], check=True)
sys.argv = [sys.argv[0], case]
ns = {}
exec(open(os.path.join(HERE, "make_laser_cht_case.py")).read().split("# ------------------------------------------------------------------ mesh")[0], ns)
w, foam, field, PAMB, T0 = ns["w"], ns["foam"], None, ns["PAMB"], ns["T0"]
C = json.load(open(os.path.join(case, "caseParams.json")))
P0 = PAMB + C["p0_gauge"]
YO2_AIR = 0.233

# ---------------------------------------------------------------- gas: O2/N2 mixture, same as run003
w("constant/air/thermophysicalProperties", foam({
    "thermoType": {"type": "heRhoThermo", "mixture": "reactingMixture", "transport": "const", "thermo": "hConst",
                   "energy": "sensibleEnthalpy", "equationOfState": "perfectGas", "specie": "specie"},
    "_": None, "inertSpecie": "N2", "chemistryReader": "foamChemistryReader",
    "foamChemistryFile": '"<constant>/air/reactions"', "foamChemistryThermoFile": '"<constant>/air/thermo.species"'}))
w("constant/air/reactions", "species\n(\n    O2\n    N2\n);\n\nreactions\n{\n}\n")
w("constant/air/thermo.species", foam({"O2": {"specie": {"molWeight": 31.9988}, "thermodynamics": {"Cp": 918, "Hf": 0},
                                              "transport": {"mu": 2.06e-05, "Pr": 0.72}}}) + "\n"
  + foam({"N2": {"specie": {"molWeight": 28.0134}, "thermodynamics": {"Cp": 1040, "Hf": 0},
                 "transport": {"mu": 1.79e-05, "Pr": 0.71}}}))
w("constant/air/combustionProperties", foam({"combustionModel": "none"}))
w("constant/air/chemistryProperties", foam({"chemistryType": {"solver": "none"}, "chemistry": "off",
                                            "initialChemicalTimeStep": "1e-07"}))
# ---------------------------------------------------------------- plate: fast thermal response (steady result unchanged)
STEEL = dict(rho=7850 / P["RHOSCALE"], Cp=460, kappa=45)
w("constant/plate/thermophysicalProperties", foam({
    "thermoType": {"type": "heSolidThermo", "mixture": "pureMixture", "transport": "constIso", "thermo": "hConst",
                   "equationOfState": "rhoConst", "specie": "specie", "energy": "sensibleEnthalpy"},
    "_": None,
    "mixture": {"specie": {"molWeight": 55.85}, "transport": {"kappa": STEEL["kappa"]},
                "thermodynamics": {"Hf": 0, "Cp": STEEL["Cp"]},
                "equationOfState": {"rho": f"{STEEL['rho']:.6g}"}}}))
# ---------------------------------------------------------------- system
cd = open(os.path.join(case, "system/controlDict")).read()
w("system/controlDict", foam({
    "application": "chtMultiRegionFoam", "startFrom": "latestTime", "startTime": 0, "stopAt": "endTime",
    "endTime": P["endTime"], "deltaT": "1e-10", "writeControl": "adjustable", "writeInterval": "1e-4",
    "purgeWrite": 0, "writeFormat": "binary", "writePrecision": 8, "timeFormat": "general", "timePrecision": 8,
    "runTimeModifiable": "true", "adjustTimeStep": "yes", "maxCo": P["maxCo"], "maxDi": 1000, "maxDeltaT": "1e-6",
    "_": None,
    "functions": {
        "massFlow": {"type": "surfaceFieldValue", "libs": "(fieldFunctionObjects)", "region": "air",
                     "regionType": "patch", "name": "inlet", "operation": "sum", "fields": "(phi)",
                     "writeFields": "false", "writeControl": "timeStep", "writeInterval": 20},
        "heatToGas": {"type": "wallHeatFlux", "libs": "(fieldFunctionObjects)", "region": "air",
                      "patches": "(air_to_plate)", "writeControl": "writeTime"},
        "plateT": {"type": "volFieldValue", "libs": "(fieldFunctionObjects)", "region": "plate", "regionType": "all",
                   "operation": "volAverage", "fields": "(T)", "writeFields": "false",
                   "writeControl": "timeStep", "writeInterval": 50},
        "MachNo": {"type": "MachNo", "libs": "(fieldFunctionObjects)", "region": "air",
                   "executeControl": "writeTime", "writeControl": "writeTime"}}}))
w("system/fvSolution", foam({"PIMPLE": {"nOuterCorrectors": P["nOuter"]}}))
w("constant/air/fvOptions", foam({"limitT": {"type": "limitTemperature", "active": "yes", "selectionMode": "all", "min": 50, "max": 2500}}))
w("system/air/fvSchemes", foam({
    "ddtSchemes": {"default": "Euler"},
    "gradSchemes": {"default": "Gauss linear", "grad(U)": "cellLimited Gauss linear 1"},
    "divSchemes": {"default": "none", "div(phi,U)": "Gauss linearUpwindV grad(U)", "div(phi,Yi_h)": "Gauss upwind",
                   "div(phi,K)": "Gauss upwind", "div(phi,h)": "Gauss upwind", "div(phi,k)": "Gauss upwind",
                   "div(phi,omega)": "Gauss upwind", "div(((rho*nuEff)*dev2(T(grad(U)))))": "Gauss linear"},
    "laplacianSchemes": {"default": "Gauss linear limited corrected 0.5"},
    "interpolationSchemes": {"default": "linear"},
    "snGradSchemes": {"default": "limited corrected 0.5"},
    "wallDist": {"method": "meshWave"}}))
w("system/air/fvSolution", foam({
    "solvers": {'"rho.*"': {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-8", "relTol": 0},
                "p_rgh": {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-9", "relTol": 0.01, "maxIter": 200},
                "p_rghFinal": {"solver": "GAMG", "smoother": "GaussSeidel", "tolerance": "1e-9", "relTol": 0, "maxIter": 200},
                '"(U|h|k|omega|Yi).*"': {"solver": "PBiCGStab", "preconditioner": "DILU", "tolerance": "1e-9", "relTol": 0}},
    "_1": None,
    "PIMPLE": {"momentumPredictor": "yes", "nCorrectors": 2, "nNonOrthogonalCorrectors": 1},
    "_2": None,
    "relaxationFactors": {"equations": {'".*"': 1}}}))
w("system/plate/fvSchemes", foam({
    "ddtSchemes": {"default": "Euler"}, "gradSchemes": {"default": "Gauss linear"}, "divSchemes": {"default": "none"},
    "laplacianSchemes": {"default": "none", "laplacian(alpha,h)": "Gauss linear corrected"},
    "interpolationSchemes": {"default": "linear"}, "snGradSchemes": {"default": "corrected"}}))
w("system/plate/fvSolution", foam({
    "solvers": {'"h.*"': {"solver": "PCG", "preconditioner": "DIC", "tolerance": "1e-9", "relTol": 0}},
    "_": None, "PIMPLE": {"nNonOrthogonalCorrectors": 0}}))
# ---------------------------------------------------------------- inlet pressure ramp in time; species fields


def fld(region, name, internal, bc, cls="volScalarField", dim="[0 0 0 0 0 0 0]"):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if "(" in p or "*" in p else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16}{v};\n" for k, v in d.items()) + "    }\n"
    body += '    "(front|back)"\n    {\n        type            wedge;\n    }\n}\n'
    w(f"0.orig/{region}/{name}", body, cls)


ramp = f"table ((0 {PAMB}) ({P['ramp']} {P0}) (1 {P0}))"
fld("air", "p_rgh", f"uniform {PAMB}", {
    "inlet": {"type": "uniformTotalPressure", "p0": ramp, "value": f"uniform {PAMB}"},
    "ambient": {"type": "totalPressure", "p0": f"uniform {PAMB}", "value": f"uniform {PAMB}"},
    "(nozzle|air_to_plate)": {"type": "fixedFluxPressure", "value": f"uniform {PAMB}"}}, dim="[1 -1 -2 0 0 0 0]")
for sp, a, b in (("O2", P["YO2_in"], YO2_AIR), ("N2", 1 - P["YO2_in"], 1 - YO2_AIR)):
    fld("air", sp, f"uniform {b}", {
        "inlet": {"type": "fixedValue", "value": f"uniform {a}"},
        "ambient": {"type": "inletOutlet", "inletValue": f"uniform {b}", "value": f"uniform {b}"},
        "(nozzle|air_to_plate)": {"type": "zeroGradient"}})
fld("air", "Ydefault", "uniform 0", {"(inlet|ambient)": {"type": "inletOutlet", "inletValue": "uniform 0", "value": "uniform 0"},
                                     "(nozzle|air_to_plate)": {"type": "zeroGradient"}})
C.update(dict(solver="chtMultiRegionFoam", endTime=P["endTime"], ramp=P["ramp"], RHOSCALE=P["RHOSCALE"], YO2_in=P["YO2_in"]))
json.dump(C, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
print("transient CHT case ready:", case)
