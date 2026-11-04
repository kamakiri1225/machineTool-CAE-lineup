#!/usr/bin/env python3
"""spindle (shaft + bearings + housing) heat-up, axisymmetric (5 deg wedge), solidFoam.

Geometry [mm] (axis = y, radial = x):
  shaft    r 0-25,  y 0-300
  bearings r 25-45, front y 40-60 (Q_front), rear y 230-250 (Q_rear)   -> heat sources
  housing  r 45-90, y 20-270, outer surface cooled by a water jacket (h_jacket, T_cool)
  other surfaces: air convection (h_air, T_air)
All parts steel (one solid region). Transient 0..tEnd s.
usage: python3 make_spindle_case.py <caseDir> [json]
"""
import json, math, os, sys

case = os.path.abspath(sys.argv[1])
P = dict(Qf=150.0, Qr=80.0, hJ=800.0, hA=15.0, Tc=293.15, Ta=293.15, T0=293.15, tEnd=3600, dt=5)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    object      {obj};\n}}\n\n'


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path); os.makedirs(os.path.dirname(full), exist_ok=True)
    open(full, "w").write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n")


R = [0, 25, 45, 90]
Y = [0, 20, 40, 60, 230, 250, 270, 300]
nR = [12, 8, 16]
nY = [4, 4, 4, 34, 4, 4, 6]
A = math.radians(2.5)
verts, vidx = [], {}


def V(ri, yi, s):
    r = R[ri]; key = (ri, yi, 0 if r == 0 else s)
    if key not in vidx:
        vidx[key] = len(verts); verts.append((r * math.cos(A), Y[yi], s * r * math.sin(A)))
    return vidx[key]


def hexb(i, j):
    return [V(i, j, -1), V(i + 1, j, -1), V(i + 1, j + 1, -1), V(i, j + 1, -1),
            V(i, j, 1), V(i + 1, j, 1), V(i + 1, j + 1, 1), V(i, j + 1, 1)]


present = set()
for j in range(7):
    present.add((0, j))                                  # shaft
for j in (2, 4):
    present.add((1, j))                                  # bearings
for j in range(1, 6):
    present.add((2, j))                                  # housing
zone = {(1, 2): "bearingFront", (1, 4): "bearingRear"}
H = {k: hexb(*k) for k in sorted(present)}
FACE = {"ymin": (0, 4, 5, 1), "ymax": (3, 2, 6, 7), "rmin": (0, 3, 7, 4), "rmax": (1, 5, 6, 2),
        "front": (0, 1, 2, 3), "back": (4, 7, 6, 5)}
nb = {"ymin": (0, -1), "ymax": (0, 1), "rmin": (-1, 0), "rmax": (1, 0)}
jacket, air = [], []
for (i, j), h in H.items():
    for side, (di, dj) in nb.items():
        ni, nj = i + di, j + dj
        if (ni, nj) in present or (side == "rmin" and i == 0):
            continue
        f = "(" + " ".join(str(h[k]) for k in FACE[side]) + ")"
        (jacket if (side == "rmax" and i == 2) else air).append(f)
blk = "\n".join(f"    hex ({' '.join(map(str, h))}) {zone.get(k, 'steel')} ({nR[k[0]]} {nY[k[1]]} 1) simpleGrading (1 1 1)" for k, h in H.items())
pat = lambda n, t, fs: f"    {n}\n    {{\n        type {t};\n        faces\n        (\n" + "\n".join("            " + x for x in fs) + "\n        );\n    }"
w("system/blockMeshDict", "scale 1e-3;\n\nvertices\n(\n" + "\n".join(f"    ({x:.9g} {y:.9g} {z:.9g})" for x, y, z in verts)
  + f"\n);\n\nblocks\n(\n{blk}\n);\n\nedges ();\n\nboundary\n(\n" + "\n".join([
      pat("jacket", "wall", jacket), pat("air", "wall", air),
      pat("front", "wedge", ["(" + " ".join(str(h[k]) for k in FACE["front"]) + ")" for h in H.values()]),
      pat("back", "wedge", ["(" + " ".join(str(h[k]) for k in FACE["back"]) + ")" for h in H.values()])]) + "\n);\n")
w("constant/thermophysicalProperties", '''thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie          { molWeight 55.8; }
    transport       { kappa 45; }
    thermodynamics  { Hf 0; Cp 460; }
    equationOfState { rho 7850; }
}
''')
w("constant/radiationProperties", "radiation off;\nradiationModel none;\n")
vol = math.pi * ((45e-3) ** 2 - (25e-3) ** 2) * 20e-3
w("constant/fvOptions", f'''bearingFrontHeat
{{
    type            scalarSemiImplicitSource;
    selectionMode   cellZone;
    cellZone        bearingFront;
    volumeMode      specific;           // [W/m3] (Q_front / bearing volume)
    sources {{ h {{ explicit {P["Qf"]/vol:.6g}; implicit 0; }} }}
}}
bearingRearHeat
{{
    type            scalarSemiImplicitSource;
    selectionMode   cellZone;
    cellZone        bearingRear;
    volumeMode      specific;
    sources {{ h {{ explicit {P["Qr"]/vol:.6g}; implicit 0; }} }}
}}
''')
w("system/controlDict", f'''application     solidFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {P["tEnd"]};
deltaT          {P["dt"]};
writeControl    runTime;
writeInterval   120;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
functions
{{
    Tprobe
    {{
        type            probes;
        libs            (sampling);
        fields          (T);
        probeLocations  ((0.001 0.0 0) (0.035 0.05 0) (0.035 0.24 0) (0.089 0.15 0));
        writeControl    timeStep;
        writeInterval   2;
    }}
}}
''')
w("system/fvSchemes", '''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''')
w("system/fvSolution", '''solvers
{
    "h.*" { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0; }
}
PIMPLE { nOuterCorrectors 1; nNonOrthogonalCorrectors 0; }
''')
ext = lambda h, Ta: {"type": "externalWallHeatFluxTemperature", "mode": "coefficient", "kappaMethod": "solidThermo",
                     "h": f"uniform {h}", "Ta": f"constant {Ta}", "value": f"uniform {P['T0']}"}
bc = lambda d: "".join(f"        {k:<16} {v};\n" for k, v in d.items())
w("0/T", f"dimensions [0 0 0 1 0 0 0];\n\ninternalField uniform {P['T0']};\n\nboundaryField\n{{\n"
  f"    jacket\n    {{\n{bc(ext(P['hJ'], P['Tc']))}    }}\n    air\n    {{\n{bc(ext(P['hA'], P['Ta']))}    }}\n"
  "    \"(front|back)\" { type wedge; }\n}\n", cls="volScalarField")
w("0/p", "dimensions [1 -1 -2 0 0 0 0];\n\ninternalField uniform 1e5;\n\nboundaryField\n{\n    \"(jacket|air)\" { type calculated; value uniform 1e5; }\n    \"(front|back)\" { type wedge; }\n}\n", cls="volScalarField")
json.dump(P, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
open(os.path.join(case, "Allrun"), "w").write("#!/bin/sh\ncd \"${0%/*}\" || exit\nblockMesh > log.blockMesh 2>&1 && checkMesh > log.checkMesh 2>&1 && solidFoam > log.solidFoam 2>&1\n")
os.chmod(os.path.join(case, "Allrun"), 0o755)
print("ok", case)
