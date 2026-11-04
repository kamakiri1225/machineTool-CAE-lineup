#!/usr/bin/env python3
"""run001: two machine tools in a warm factory (after Klimant et al. 2026 / Franke: 20 degC machine in 30 degC air,
1 m/s forced convection from an open door, floor kept at 20 degC). Machine A: bare surface. Machine B: covered by
30 mm insulation (PUR, k = 0.03 W/mK, modelled as a thin thermal resistance on its surface).
chtMultiRegionFoam: stage 1 (0-300 s) transient buoyant flow; stage 2 (300 s - tEnd) frozenFlow, large time steps.
usage: python3 make_factory_case.py <caseDir> [json]
"""
import json, os, re, shutil, sys

case = os.path.abspath(sys.argv[1])
P = dict(Tair=303.15, Tm=293.15, Tfloor=293.15, Udoor=1.0, tEnd=21600, tFlow=300, insThick=0.03, insK=0.03)
P.update(json.loads(sys.argv[2]) if len(sys.argv) > 2 else {})
TUT = "/usr/lib/openfoam/openfoam2512/tutorials/heatTransfer/chtMultiRegionFoam/multiRegionHeater"
if os.path.exists(case):
    for d in ("0.orig", "constant", "system"):
        shutil.rmtree(os.path.join(case, d), ignore_errors=True)
os.makedirs(case, exist_ok=True)
for d in ("constant", "system"):
    shutil.copytree(os.path.join(TUT, d), os.path.join(case, d), symlinks=False)
shutil.copytree(os.path.join(TUT, "0.orig"), os.path.join(case, "0.orig"))
for r in ("bottomWater", "heater", "leftSolid", "rightSolid", "topAir"):
    for d in ("constant", "system"):
        p = os.path.join(case, d, r)
        if os.path.isdir(p) and r not in ("topAir", "heater"):
            shutil.rmtree(p)
for d in ("constant", "system"):
    os.rename(os.path.join(case, d, "topAir"), os.path.join(case, d, "air"))
    shutil.copytree(os.path.join(case, d, "heater"), os.path.join(case, d, "machineA"))
    os.rename(os.path.join(case, d, "heater"), os.path.join(case, d, "machineB"))
HDR = 'FoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       {cls};\n    object      {obj};\n}}\n\n'


def w(path, body, cls="dictionary"):
    full = os.path.join(case, path); os.makedirs(os.path.dirname(full), exist_ok=True)
    open(full, "w").write(HDR.format(cls=cls, obj=os.path.basename(path)) + body + "\n")


w("constant/regionProperties", "regions\n(\n    fluid (air)\n    solid (machineA machineB)\n);\n")
w("constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 -9.81 0);\n", cls="uniformDimensionedVectorField")
for m in ("machineA", "machineB"):     # cast iron
    p = os.path.join(case, "constant", m, "thermophysicalProperties"); s = open(p).read()
    s = re.sub(r"kappa\s+[\d.]+;", "kappa   50;", s); s = re.sub(r"Cp\s+[\d.]+;", "Cp      500;", s); s = re.sub(r"rho\s+[\d.]+;", "rho     7200;", s)
    open(p, "w").write(s)
p = os.path.join(case, "constant/air/thermophysicalProperties"); s = open(p).read()
s = re.sub(r"Cp\s+[\d.]+;", "Cp              1005;", s); open(p, "w").write(s)

w("system/blockMeshDict", '''vertices
(
    (0 0 0) (6 0 0) (6 3.6 0) (0 3.6 0)
    (0 0 4.2) (6 0 4.2) (6 3.6 4.2) (0 3.6 4.2)
);
blocks ( hex (0 1 2 3 4 5 6 7) (40 24 28) simpleGrading (1 1 1) );
edges ();
boundary
(
    floor   { type wall; faces ((1 5 4 0)); }
    ceiling { type wall; faces ((3 7 6 2)); }
    walls   { type wall; faces ((0 4 7 3) (2 6 5 1) (0 3 2 1) (4 5 6 7)); }
);
''')
w("system/topoSetDict", '''actions
(
    { name machineA; type cellSet; action new; source boxToCell; box (1.66 -1 1.36) (2.84 1.79 2.84); }
    { name machineA; type cellZoneSet; action new; source setToCellZone; set machineA; }
    { name machineB; type cellSet; action new; source boxToCell; box (3.46 -1 1.36) (4.64 1.79 2.84); }
    { name machineB; type cellZoneSet; action new; source setToCellZone; set machineB; }
    { name air; type cellSet; action new; source cellToCell; set machineA; }
    { name air; type cellSet; action add; source cellToCell; set machineB; }
    { name air; type cellSet; action invert; }
    { name air; type cellZoneSet; action new; source setToCellZone; set air; }
    { name doorFaces; type faceSet; action new; source boxToFace; box (-0.01 0 1.5) (0.01 2.1 2.7); }
    { name ventFaces; type faceSet; action new; source boxToFace; box (5.99 2.7 1.5) (6.01 3.3 2.7); }
);
''')
w("system/createPatchDict", '''pointSync false;
patches
(
    { name door; patchInfo { type patch; } constructFrom set; set doorFaces; }
    { name vent; patchInfo { type patch; } constructFrom set; set ventFaces; }
);
''')
Ta, Tm, Tf, Ud = P["Tair"], P["Tm"], P["Tfloor"], P["Udoor"]
w("system/air/changeDictionaryDict", f'''U
{{
    internalField   uniform (0 0 0);
    boundaryField
    {{
        ".*"  {{ type noSlip; }}
        door  {{ type fixedValue; value uniform ({Ud} 0 0); }}
        vent  {{ type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }}
    }}
}}
T
{{
    internalField   uniform {Ta};
    boundaryField
    {{
        ".*"    {{ type zeroGradient; value uniform {Ta}; }}
        floor   {{ type fixedValue; value uniform {Tf}; }}
        door    {{ type fixedValue; value uniform {Ta}; }}
        vent    {{ type inletOutlet; inletValue uniform {Ta}; value uniform {Ta}; }}
        "air_to_.*" {{ type compressible::turbulentTemperatureRadCoupledMixed; Tnbr T; kappaMethod fluidThermo; value uniform {Ta}; }}
    }}
}}
epsilon
{{
    internalField   uniform 0.01;
    boundaryField
    {{
        ".*"  {{ type epsilonWallFunction; value uniform 0.01; }}
        door  {{ type fixedValue; value uniform 0.01; }}
        vent  {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    }}
}}
k
{{
    internalField   uniform 0.01;
    boundaryField
    {{
        ".*"  {{ type kqRWallFunction; value uniform 0.01; }}
        door  {{ type fixedValue; value uniform 0.01; }}
        vent  {{ type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }}
    }}
}}
p_rgh
{{
    internalField   uniform 1e5;
    boundaryField
    {{
        ".*"  {{ type fixedFluxPressure; value uniform 1e5; }}
        vent  {{ type fixedValue; value uniform 1e5; }}
    }}
}}
p
{{
    internalField   uniform 1e5;
    boundaryField
    {{
        ".*"  {{ type calculated; value uniform 1e5; }}
    }}
}}
''')
for m, extra in (("machineA", ""), ("machineB", f"thicknessLayers ({P['insThick']}); kappaLayers ({P['insK']});")):
    w(f"system/{m}/changeDictionaryDict", f'''T
{{
    internalField   uniform {Tm};
    boundaryField
    {{
        ".*"   {{ type zeroGradient; value uniform {Tm}; }}
        floor  {{ type fixedValue; value uniform {Tf}; }}
        "{m}_to_.*" {{ type compressible::turbulentTemperatureRadCoupledMixed; Tnbr T; kappaMethod solidThermo; {extra} value uniform {Tm}; }}
    }}
}}
''')
cd = os.path.join(case, "system/controlDict"); s = open(cd).read()
for k, v in (("startTime", "0"), ("endTime", str(P["tFlow"])), ("deltaT", "0.01"), ("writeInterval", "60"), ("maxCo", "0.8"), ("maxDi", "100.0")):
    s = re.sub(rf"\n{k}\s+[^;]+;", f"\n{k}            {v};", s)
s = s.replace("maxDi           100.0;", "maxDi           100.0;\n\nmaxDeltaT       0.5;")
s = re.sub(r"\nfunctions\s*\{.*", "", s, flags=re.S)
s += '''
functions
{
    TCP
    {
        type            probes;
        libs            (sampling);
        region          machineA;
        fields          (T);
        probeLocations  ((2.25 0.9 2.1) (2.25 1.75 2.1));
        writeControl    timeStep;
        writeInterval   1;
    }
    TCPB
    {
        type            probes;
        libs            (sampling);
        region          machineB;
        fields          (T);
        probeLocations  ((4.05 0.9 2.1) (4.05 1.75 2.1));
        writeControl    timeStep;
        writeInterval   1;
    }
}
'''
open(cd, "w").write(s)
open(os.path.join(case, "Allrun"), "w").write(f'''#!/bin/sh
cd "${{0%/*}}" || exit
. ${{WM_PROJECT_DIR:?}}/bin/tools/RunFunctions
runApplication blockMesh
runApplication topoSet
runApplication createPatch -overwrite
restore0Dir
runApplication splitMeshRegions -cellZones -overwrite
for region in machineA machineB; do remove0DirFields -region "$region" -- nut alphat epsilon k U p_rgh; done
for region in air machineA machineB; do runApplication -s $region changeDictionary -region $region; done
# stage 1: buoyant flow develops ({P["tFlow"]} s)
runApplication -s stage1 chtMultiRegionFoam
# stage 2: freeze the flow, march the temperatures with large time steps
foamDictionary -entry PIMPLE.frozenFlow -add yes system/air/fvSolution
foamDictionary -entry endTime -set {P["tEnd"]} system/controlDict
foamDictionary -entry writeInterval -set 600 system/controlDict
foamDictionary -entry maxCo -set 1e5 system/controlDict
foamDictionary -entry maxDeltaT -set 30 system/controlDict
runApplication -s stage2 chtMultiRegionFoam
touch post.foam
''')
os.chmod(os.path.join(case, "Allrun"), 0o755)
json.dump(P, open(os.path.join(case, "caseParams.json"), "w"), indent=1)
open(os.path.join(case, "post.foam"), "w").close()
print("ok", case)
