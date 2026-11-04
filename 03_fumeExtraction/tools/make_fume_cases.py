#!/usr/bin/env python3
"""Laser-cutting fume extraction (downdraft / under-table suction), two-step model.

  flow/      simpleFoam (steady, incompressible, k-omega SST): air drawn in through a roof opening,
             around the sheet edges and through the cut, and sucked out of the dust chamber
             below the sheet by two side suction openings.
  particles/ icoUncoupledKinematicParcelFoam (transient): fume / dross particles injected at the
             cut point ride on the frozen steady velocity field (one-way coupling).

Geometry [m]: cabin x 0..1.2, y -0.3..0.6 (sheet at y = 0), z 0..0.8
  sheet (baffle wall) x 0.1..1.1, z 0.1..0.7 with the cut opening 50 x 50 mm at (0.6, 0, 0.4)
  dust chamber y -0.3..0, suction openings on x = 0 and x = 1.2 walls (y -0.3..-0.1, z 0.3..0.5)
usage: python3 make_fume_cases.py <runDir> [json]
"""
import json
import os
import sys

run = os.path.abspath(sys.argv[1])
P = dict(Q=0.30, cell=0.025, endTime=3.0)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})

HDR = '''FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

'''


def w(case, path, body, cls="dictionary"):
    full = os.path.join(run, case, path)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w") as f:
        f.write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n// ************************************************************************* //\n")


c = P["cell"]
nx, ny, nz = round(1.2 / c), round(0.9 / c), round(0.8 / c)
# ------------------------------------------------------------------ flow case
w("flow", "system/blockMeshDict", f'''vertices
(
    (0 -0.3 0) (1.2 -0.3 0) (1.2 0.6 0) (0 0.6 0)
    (0 -0.3 0.8) (1.2 -0.3 0.8) (1.2 0.6 0.8) (0 0.6 0.8)
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
edges ();
boundary
(
    roofOpening {{ type patch; faces ((3 7 6 2)); }}
    walls       {{ type wall;  faces ((2 6 5 1) (0 4 7 3) (0 3 2 1) (4 5 6 7) (1 5 4 0)); }}
);
''')
e = 1e-3
w("flow", "system/topoSetDict", f'''actions
(
    {{ name sheetFaces; type faceSet; action new;      source boxToFace; box (0.1 {-e} 0.1) (1.1 {e} 0.7); }}
    {{ name sheetFaces; type faceSet; action subtract; source boxToFace; box (0.575 {-e} 0.375) (0.625 {e} 0.425); }}
    {{ name sheet;      type faceZoneSet; action new;  source setToFaceZone; faceSet sheetFaces; }}
    {{ name suctionFaces; type faceSet; action new;    source boxToFace; box ({-e} -0.3 0.3) ({e} -0.1 0.5); }}
    {{ name suctionFaces; type faceSet; action add;    source boxToFace; box ({1.2-e} -0.3 0.3) ({1.2+e} -0.1 0.5); }}
);
''')
w("flow", "system/createBafflesDict", '''internalFacesOnly true;
baffles
{
    sheet
    {
        type        faceZone;
        zoneName    sheet;
        patchPairs
        {
            type        wall;
        }
    }
}
''')
w("flow", "system/createPatchDict", '''pointSync false;
patches
(
    { name suction; patchInfo { type patch; } constructFrom set; set suctionFaces; }
);
''')
w("flow", "system/controlDict", '''application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         1500;
deltaT          1;
writeControl    timeStep;
writeInterval   500;
purgeWrite      1;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
runTimeModifiable true;
functions
{
    suctionFlow { type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name suction;
                  operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 10; }
    roofFlow    { type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name roofOpening;
                  operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 10; }
}
''')
w("flow", "system/fvSchemes", '''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwind grad(U);
    div(phi,k) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
''')
w("flow", "system/fvSolution", '''solvers
{
    p { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.05; }
    "(U|k|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
    residualControl { p 1e-4; U 1e-5; "(k|omega)" 1e-5; }
}
relaxationFactors { equations { U 0.9; ".*" 0.9; } }
''')
w("flow", "constant/transportProperties", "transportModel  Newtonian;\nnu              1.5e-05;\n")
w("flow", "constant/turbulenceProperties", "simulationType RAS;\nRAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }\n")


def field(case, name, cls, dim, internal, bc):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for p, d in bc.items():
        pn = f'"{p}"' if "(" in p else p
        body += f"    {pn}\n    {{\n" + "".join(f"        {k:<16} {val};\n" for k, val in d.items()) + "    }\n"
    w(case, f"0.orig/{name}", body + "}\n", cls=cls)


S, WL = "volScalarField", "(walls|sheet.*)"
field("flow", "U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
    "roofOpening": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    "suction": {"type": "flowRateOutletVelocity", "volumetricFlowRate": f"constant {P['Q']}", "value": "uniform (0 0 0)"},
    WL: {"type": "noSlip"}})
field("flow", "p", S, "[0 2 -2 0 0 0 0]", "uniform 0", {
    "roofOpening": {"type": "totalPressure", "p0": "uniform 0", "value": "uniform 0"},
    "suction": {"type": "zeroGradient"}, WL: {"type": "zeroGradient"}})
field("flow", "k", S, "[0 2 -2 0 0 0 0]", "uniform 1e-3", {
    "roofOpening": {"type": "inletOutlet", "inletValue": "uniform 1e-3", "value": "uniform 1e-3"},
    "suction": {"type": "zeroGradient"}, WL: {"type": "kqRWallFunction", "value": "uniform 1e-3"}})
field("flow", "omega", S, "[0 0 -1 0 0 0 0]", "uniform 1", {
    "roofOpening": {"type": "inletOutlet", "inletValue": "uniform 1", "value": "uniform 1"},
    "suction": {"type": "zeroGradient"}, WL: {"type": "omegaWallFunction", "value": "uniform 1"}})
field("flow", "nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", {
    "(roofOpening|suction)": {"type": "calculated", "value": "uniform 0"},
    WL: {"type": "nutkWallFunction", "value": "uniform 0"}})

# ------------------------------------------------------------------ particle case (mesh + U copied by Allrun)
w("particles", "system/controlDict", f'''application     icoUncoupledKinematicParcelFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {P["endTime"]};
deltaT          2e-3;
writeControl    adjustable;
writeInterval   0.05;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
''')
w("particles", "system/fvSchemes", '''ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
''')
w("particles", "system/fvSolution", "solvers {}\nPIMPLE { nOuterCorrectors 1; }\n")
w("particles", "constant/transportProperties", "transportModel  Newtonian;\nnu              1.5e-05;\nrhoInf          1.2;\n")
w("particles", "constant/turbulenceProperties", "simulationType laminar;\n")
w("particles", "constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 -9.81 0);\n", cls="uniformDimensionedVectorField")
size = '''sizeDistribution
            {
                type        RosinRammler;
                RosinRammlerDistribution
                {
                    minValue    2e-06;
                    maxValue    5e-05;
                    lambda      1.5e-05;
                    n           2;
                }
            }'''
w("particles", "constant/kinematicCloudProperties", f'''solution
{{
    active          true;
    coupled         false;
    transient       yes;
    cellValueSourceCorrection off;
    maxCo           0.3;
    sourceTerms {{ schemes {{}} }}
    interpolationSchemes {{ rho cell; U cellPoint; mu cell; }}
    integrationSchemes {{ U Euler; }}
}}

constantProperties
{{
    rho0            5000;      // iron oxide fume / dross
}}

subModels
{{
    particleForces
    {{
        sphereDrag;
        gravity;
    }}

    injectionModels
    {{
        plumeUp        // fume rising from the cut (thermal plume)
        {{
            type            coneInjection;
            SOI             0;
            duration        2.0;
            positionAxis    ( ((0.6 0.01 0.4) (0 1 0)) );
            parcelsPerInjector 6000;
            flowRateProfile constant 1;
            massTotal       1e-5;
            parcelBasisType fixed;
            nParticle       1;
            Umag            constant 1.5;
            thetaInner      constant 0;
            thetaOuter      constant 35;
            {size}
        }}
        kerfDown       // fume / dross blown down through the kerf by the assist gas
        {{
            type            coneInjection;
            SOI             0;
            duration        2.0;
            positionAxis    ( ((0.6 -0.01 0.4) (0 -1 0)) );
            parcelsPerInjector 6000;
            flowRateProfile constant 1;
            massTotal       1e-5;
            parcelBasisType fixed;
            nParticle       1;
            Umag            constant 20;
            thetaInner      constant 0;
            thetaOuter      constant 10;
            {size}
        }}
    }}

    dispersionModel none;

    patchInteractionModel localInteraction;
    localInteractionCoeffs
    {{
        patches
        (
            "(roofOpening|suction)" {{ type escape; }}
            "(walls|sheet.*)"       {{ type stick; }}
        );
    }}

    surfaceFilmModel none;
    stochasticCollisionModel none;
    collisionModel none;
}}

cloudFunctions
{{
    particlePostProcessing1
    {{
        type            particlePostProcessing;
        maxStoredParcels 100000;
        patches         (suction roofOpening);
    }}
}}
''')

with open(os.path.join(run, "Allrun"), "w") as f:
    f.write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
# ---- 1) steady flow
( cd flow && rm -rf 0 \\
  && runApplication blockMesh && runApplication topoSet \\
  && runApplication createPatch -overwrite && runApplication createBaffles -overwrite \\
  && cp -r 0.orig 0 && runApplication checkMesh && runApplication simpleFoam && touch post.foam )
# ---- 2) particles on the frozen flow field
L=$(cd flow && foamListTimes -latestTime)
( cd particles && rm -rf 0 constant/polyMesh && cp -r 0.orig 0 \\
  && cp -r ../flow/constant/polyMesh constant/ && cp ../flow/$L/U 0/U \\
  && runApplication icoUncoupledKinematicParcelFoam && touch post.foam )
''')
os.chmod(os.path.join(run, "Allrun"), 0o755)
# rho, mu for the particle solver
for name, val, dim in (("rho", 1.2, "[1 -3 0 0 0 0 0]"), ("mu", 1.8e-5, "[1 -1 -1 0 0 0 0]")):
    w("particles", f"0.orig/{name}", f"dimensions      {dim};\n\ninternalField   uniform {val};\n\nboundaryField\n{{\n    \".*\"\n    {{\n        type calculated;\n        value uniform {val};\n    }}\n}}\n", cls="volScalarField")
json.dump(P, open(os.path.join(run, "caseParams.json"), "w"), indent=1)
print("cells", nx * ny * nz)
