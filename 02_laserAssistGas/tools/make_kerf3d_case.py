#!/usr/bin/env python3
"""3-D half model of the assist-gas flow into a laser-cutting slot (kerf), after
Matsumura, Kiwata, Yamaguchi (JSME 2017 annual meeting, S0520603) and the Tokai Univ. CFD page.

Coordinates [mm]: x = cutting direction (the slot is already cut for x < 0, the plate is still
uncut for x > 0, i.e. the cutting front is the plane x = 0), y = vertical (plate top y = 0),
z = across the slot. Only z >= 0 is meshed (symmetry plane z = 0 = slot centre plane).

  convergent nozzle: exit diameter DE at y = H, bore inlet diameter DIN at y = H + LN,
                     outer radius RNOZ, nozzle axis at (x, z) = (XOFF, 0)
  plate: thickness T, slot width W (x < 0)
Solver: rhoCentralFoam (density based, compressible, air, realizable k-epsilon), run to a steady
state in time. The assist-gas fraction (N2 or O2 coming from the nozzle) is a passive scalar
"Yassist" (scalarTransport function object); air properties are used for the mixture.
usage: python3 make_kerf3d_case.py <caseDir> '<json>'
"""
import json
import math
import os
import sys

case = os.path.abspath(sys.argv[1])
P = dict(p0=2.0e6, DE=1.5, DIN=3.0, LN=3.0, RNOZ=3.0, H=1.0, T=1.2, W=0.3, XOFF=0.0,
         gas="N2", endTime=4000, dx0=0.25, title="")
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
PAMB, T0 = 101325.0, 300.0
DE, DIN, LN, RNOZ, H, T, W, XOFF = (P[k] for k in ("DE", "DIN", "LN", "RNOZ", "H", "T", "W", "XOFF"))
XMIN, XMAX, ZMAX = -4.0, 4.0, 4.0
YMIN, YMAX = -T - 2.0, H + LN
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    object      {obj};\n}}\n\n'


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w") as f:
        f.write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n")


# ------------------------------------------------------------------ STL geometry (mm)
def stl(name, tris):
    full = os.path.join(case, "constant/triSurface", name + ".stl")
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w") as f:
        f.write(f"solid {name}\n")
        for a, b, c in tris:
            u = [b[i] - a[i] for i in range(3)]; v = [c[i] - a[i] for i in range(3)]
            n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
            L = math.sqrt(sum(x * x for x in n)) or 1.0
            f.write(f"  facet normal {n[0]/L:g} {n[1]/L:g} {n[2]/L:g}\n    outer loop\n")
            for p in (a, b, c):
                f.write(f"      vertex {p[0]:.6g} {p[1]:.6g} {p[2]:.6g}\n")
            f.write("    endloop\n  endfacet\n")
        f.write(f"endsolid {name}\n")


def revolve(profile, n=96):
    """closed polygon in (r, y) revolved about the axis (XOFF, *, 0)"""
    tris = []
    for k in range(n):
        a0, a1 = 2 * math.pi * k / n, 2 * math.pi * (k + 1) / n
        for i in range(len(profile)):
            (r1, y1), (r2, y2) = profile[i], profile[(i + 1) % len(profile)]
            p = lambda r, y, a: (XOFF + r * math.cos(a), y, r * math.sin(a))
            tris.append((p(r1, y1, a0), p(r2, y2, a0), p(r2, y2, a1)))
            tris.append((p(r1, y1, a0), p(r2, y2, a1), p(r1, y1, a1)))
    return [t for t in tris if len({t[0], t[1], t[2]}) == 3]


def box(x0, x1, y0, y1, z0, z1):
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0), (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    q = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (0, 4, 7, 3), (1, 2, 6, 5)]
    return [t for a, b, c, d in q for t in ((v[a], v[b], v[c]), (v[a], v[c], v[d]))]


ytop = YMAX + 0.5
stl("nozzle", revolve([(DE / 2, H), (RNOZ, H), (RNOZ, ytop), (DIN / 2 + (DIN - DE) / 2 * 0.5 / LN, ytop)]))
stl("plate", box(XMIN - 1, XMAX + 1, -T, 0.0, W / 2, ZMAX + 1) + box(0.0, XMAX + 1, -T, 0.0, -1.0, W / 2 + 0.01))

# ------------------------------------------------------------------ background mesh + snappy
dx = P["dx0"]
nx, ny, nz = round((XMAX - XMIN) / dx), round((YMAX - YMIN) / dx), round(ZMAX / dx)
w("system/blockMeshDict", f'''scale 1e-3;
vertices
(
    ({XMIN} {YMIN} 0) ({XMAX} {YMIN} 0) ({XMAX} {YMAX} 0) ({XMIN} {YMAX} 0)
    ({XMIN} {YMIN} {ZMAX}) ({XMAX} {YMIN} {ZMAX}) ({XMAX} {YMAX} {ZMAX}) ({XMIN} {YMAX} {ZMAX})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
edges ();
boundary
(
    symmetryPlane {{ type symmetryPlane; faces ((0 3 2 1)); }}
    top           {{ type patch; faces ((3 7 6 2)); }}
    bottom        {{ type patch; faces ((1 5 4 0)); }}
    ambient       {{ type patch; faces ((2 6 5 1) (0 4 7 3) (4 5 6 7)); }}
);
''')
s = 1e-3
w("system/snappyHexMeshDict", f'''castellatedMesh true;
snap            true;
addLayers       false;
geometry
{{
    nozzle.stl {{ type triSurfaceMesh; name nozzle; scale {s}; }}
    plate.stl  {{ type triSurfaceMesh; name plate;  scale {s}; }}
    jet   {{ type searchableCylinder; point1 ({XOFF*s} {YMIN*s} 0); point2 ({XOFF*s} {YMAX*s} 0); radius {DE/2*1.7*s}; }}
    slot  {{ type searchableBox; min ({-2.0*s} {(-T-0.3)*s} 0); max ({0.3*s} {0.2*s} {(W/2+0.1)*s}); }}
}}
castellatedMeshControls
{{
    maxLocalCells 2000000; maxGlobalCells 4000000; minRefinementCells 0; nCellsBetweenLevels 2;
    maxLoadUnbalance 0.1; allowFreeStandingZoneFaces true; resolveFeatureAngle 30;
    features ();
    refinementSurfaces
    {{
        nozzle {{ level (2 2); patchInfo {{ type wall; }} }}
        plate  {{ level (2 3); patchInfo {{ type wall; }} }}
    }}
    refinementRegions
    {{
        jet  {{ mode inside; levels ((1e15 2)); }}
        slot {{ mode inside; levels ((1e15 3)); }}
    }}
    locationInMesh ({(XMAX-0.37)*s} {(YMAX-1.03)*s} {(ZMAX-0.37)*s});
}}
snapControls
{{
    nSmoothPatch 3; tolerance 2.0; nSolveIter 50; nRelaxIter 5;
    nFeatureSnapIter 10; implicitFeatureSnap true; explicitFeatureSnap false; multiRegionFeatureSnap false;
}}
addLayersControls {{ relativeSizes true; layers {{}}; expansionRatio 1.2; finalLayerThickness 0.3; minThickness 0.1;
    nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10;
    maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedialAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50; }}
meshQualityControls {{ #includeEtc "caseDicts/meshQualityDict" }}
mergeTolerance 1e-6;
''')
w("system/topoSetDict", f'''actions
(
    {{ name inletFaces; type faceSet; action new; source patchToFace; patch top; }}
    {{ name inletFaces; type faceSet; action subset; source cylinderToFace;
       point1 ({XOFF*s} {(YMAX-1)*s} 0); point2 ({XOFF*s} {(YMAX+1)*s} 0); radius {(DIN/2+0.3)*s}; }}
);
''')
w("system/createPatchDict", '''pointSync false;
patches ( { name inlet; patchInfo { type patch; } constructFrom set; set inletFaces; } );
''')
w("system/meshQualityDict", '#includeEtc "caseDicts/meshQualityDict"\n')

# ------------------------------------------------------------------ physics
w("constant/thermophysicalProperties", '''thermoType
{
    type            hePsiThermo;
    mixture         pureMixture;
    transport       sutherland;
    thermo          hConst;
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleInternalEnergy;
}
mixture
{
    specie          { molWeight 28.96; }
    thermodynamics  { Cp 1004.5; Hf 0; }
    transport       { As 1.458e-06; Ts 110.4; Pr 0.71; }
}
''')
w("constant/turbulenceProperties", "simulationType RAS;\nRAS { RASModel realizableKE; turbulence on; printCoeffs on; }\n")
w("system/controlDict", f'''application     rhoCentralFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {P["endTime"]};
deltaT          1;          // LTS: 1 iteration per 'time step'
writeControl    timeStep;
writeInterval   500;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   8;
adjustTimeStep  no;
maxCo           0.3;          // local Courant number used for the local time step
maxDeltaT       1e-5;         // upper limit of the local (pseudo) time step [s]
rDeltaTSmoothingCoeff 0.02;
functions
{{
    Yassist
    {{
        type            scalarTransport;
        libs            (solverFunctionObjects);
        field           Yassist;
        phi             phi;
        rho             rho;
        alphaD          1;
        alphaDt         1.43;
        schemesField    Yassist;
        resetOnStartUp  false;
        writeControl    writeTime;
    }}
    MachNo      {{ type MachNo; libs (fieldFunctionObjects); executeControl writeTime; writeControl writeTime; }}
    wallShearStress {{ type wallShearStress; libs (fieldFunctionObjects); patches (plate); writeControl writeTime; }}
    slotFlow    {{ type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name bottom;
                   operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 20; }}
    inletFlow   {{ type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name inlet;
                   operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 20; }}
}}
''')
w("system/fvSchemes", '''fluxScheme      Kurganov;
ddtSchemes      { default localEuler; }   // local time stepping -> pseudo-steady
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(tauMC)      Gauss linear;
    div(phi,k)      Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,Yassist) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes
{
    default         linear;
    reconstruct(rho) vanLeer;
    reconstruct(U)  vanLeerV;
    reconstruct(T)  vanLeer;
}
snGradSchemes   { default limited corrected 0.5; }
wallDist        { method meshWave; }
''')
w("system/fvSolution", '''solvers
{
    "(rho|rhoU|rhoE)" { solver diagonal; }
    U  { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-9; relTol 0.01; }
    "(e|k|epsilon|Yassist)" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-9; relTol 0.01; }
}
''')


def field(name, cls, dim, internal, bc):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if "(" in p else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16} {val};\n" for k, val in d.items()) + "    }\n"
    body += "    symmetryPlane\n    {\n        type             symmetryPlane;\n    }\n}\n"
    w(f"0.orig/{name}", body, cls=cls)


S, WALL, AMB = "volScalarField", "(nozzle|plate)", "(top|bottom|ambient)"
P0 = P["p0"] + PAMB          # p0 is the gauge nozzle inlet pressure (as in the JSME paper)
field("p", S, "[1 -1 -2 0 0 0 0]", f"uniform {PAMB}", {
    "inlet": ({"type": "uniformTotalPressure", "p0": f"table ((0 {PAMB}) (500 {P0}) (1e9 {P0}))", "value": f"uniform {PAMB}"}
              if P.get("ramp", True) else {"type": "totalPressure", "p0": f"uniform {P0}", "value": f"uniform {P0}"}),
    AMB: {"type": "totalPressure", "p0": f"uniform {PAMB}", "value": f"uniform {PAMB}"},
    WALL: {"type": "zeroGradient"}})
field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
    "inlet": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    AMB: {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    WALL: {"type": "noSlip"}})
field("T", S, "[0 0 0 1 0 0 0]", f"uniform {T0}", {
    "inlet": {"type": "totalTemperature", "gamma": "1.4", "T0": f"uniform {T0}", "value": f"uniform {T0}"},
    AMB: {"type": "inletOutlet", "inletValue": f"uniform {T0}", "value": f"uniform {T0}"},
    WALL: {"type": "zeroGradient"}})
field("Yassist", S, "[0 0 0 0 0 0 0]", "uniform 0", {
    "inlet": {"type": "fixedValue", "value": "uniform 1"},
    AMB: {"type": "inletOutlet", "inletValue": "uniform 0", "value": "uniform 0"},
    WALL: {"type": "zeroGradient"}})
field("k", S, "[0 2 -2 0 0 0 0]", "uniform 1", {
    "inlet": {"type": "turbulentIntensityKineticEnergyInlet", "intensity": "0.05", "value": "uniform 1"},
    AMB: {"type": "inletOutlet", "inletValue": "uniform 1", "value": "uniform 1"},
    WALL: {"type": "kqRWallFunction", "value": "uniform 1"}})
field("epsilon", S, "[0 2 -3 0 0 0 0]", "uniform 1e4", {
    "inlet": {"type": "turbulentMixingLengthDissipationRateInlet", "mixingLength": "1e-4", "value": "uniform 1e4"},
    AMB: {"type": "inletOutlet", "inletValue": "uniform 1e4", "value": "uniform 1e4"},
    WALL: {"type": "epsilonWallFunction", "value": "uniform 1e4"}})
field("nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {
    "(inlet|top|bottom|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "nutkWallFunction", "value": "uniform 0"}})
field("alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", {
    "(inlet|top|bottom|ambient)": {"type": "calculated", "value": "uniform 0"},
    WALL: {"type": "compressible::alphatWallFunction", "value": "uniform 0"}})

with open(os.path.join(case, "Allrun"), "w") as f:
    f.write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
rm -rf 0 constant/polyMesh
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication topoSet
runApplication createPatch -overwrite
cp -r 0.orig 0
runApplication checkMesh
runApplication rhoCentralFoam
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
json.dump(P, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
print(case, "background", nx, ny, nz, "p0 abs", P0)
