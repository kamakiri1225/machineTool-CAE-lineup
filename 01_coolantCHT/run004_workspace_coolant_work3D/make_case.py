#!/usr/bin/env python3
"""run004: machine-tool workspace (air) + warm coolant jet hitting a workpiece on a steel table (solid), 3-D.

Built from the tutorial heatTransfer/chtMultiRegionTwoPhaseEulerFoam/solidQuenching2D
(Euler-Euler gas/liquid two-phase flow coupled with heat conduction in a solid region).

Simplified model
  - enclosure 600 x 600(height, y) x 400 mm, 20 mm cells (30 x 30 x 20)
  - steel table x 140-460, y 240-300 (mid-height), z 100-300 mm + workpiece (steel) x 240-360, y 300-360, z 140-260 mm
    (table and workpiece form one solid region "solid")
  - coolant jet centre (x 300, z 200) is directly above the workpiece centre
  - coolant inlet 40 x 40 mm in the ceiling above the table, liquid at 1.0 m/s downward,
    308 K (35 degC, warmed by the pump); about 96 L/min
  - drains: floor strips x < 100 and x > 500 mm; ceiling = air vent (open)
  - air / table / enclosure walls 293 K
  - gas phase = air (perfect gas), liquid = water-based coolant (water properties)
Run: python3 make_case.py && ./Allrun
"""
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
TUT = "/usr/lib/openfoam/openfoam2512/tutorials/heatTransfer/chtMultiRegionTwoPhaseEulerFoam/solidQuenching2D"
T_AIR, T_COOL, U_IN = 293, 308, 1.0

for d in ("0.orig", "constant", "system"):
    if os.path.exists(os.path.join(HERE, d)):
        shutil.rmtree(os.path.join(HERE, d))
    shutil.copytree(os.path.join(TUT, d), os.path.join(HERE, d))
shutil.copy(os.path.join(TUT, "Tsat_water_1_2bar.csv"), HERE)


def rw(path, fn):
    p = os.path.join(HERE, path)
    s = fn(open(p).read())
    open(p, "w").write(s)


# ---------------------------------------------------------------- mesh
def mesh(s):
    i = s.index("vertices")
    j = s.index("// ****", i)
    return s[:i] + '''vertices
(
    (0   0   0) (600 0   0) (600 600 0) (0   600 0)
    (0   0   400) (600 0   400) (600 600 400) (0   600 400)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (30 30 20) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    maxY    // ceiling: air vent (open)
    {
        type patch;
        faces ((3 7 6 2));
    }
    maxX    // enclosure side walls
    {
        type wall;
        faces ((2 6 5 1) (0 4 7 3) (0 3 2 1) (4 5 6 7));
    }
    minY    // floor
    {
        type wall;
        faces ((1 5 4 0));
    }
);

mergePatchPairs
(
);

''' + s[j:]


rw("system/blockMeshDict", mesh)

rw("system/topoSetDict", lambda s: s.replace(
    "box     (-0.01 29e-3 -1) (4.77e-3 70e-3 1);", "box     (0.141 0.241 0.101) (0.459 0.299 0.299);\n    }\n    {\n        name    heaterCellSet;\n        type    cellSet;\n        action  add;\n        source  boxToCell;\n        box     (0.241 0.299 0.141) (0.359 0.359 0.259);").replace(
    "actions\n(", '''actions
(
    {
        name    inletFaces;
        type    faceSet;
        action  new;
        source  boxToFace;
        box     (0.28 0.599 0.18) (0.32 0.601 0.22);
    }
    {
        name    drainFaces;
        type    faceSet;
        action  new;
        source  boxToFace;
        box     (-0.01 -0.001 -0.01) (0.1 0.001 0.41);
    }
    {
        name    drainFaces;
        type    faceSet;
        action  add;
        source  boxToFace;
        box     (0.5 -0.001 -0.01) (0.61 0.001 0.41);
    }
''', 1))

with open(os.path.join(HERE, "system/createPatchDict"), "w") as f:
    f.write(open(os.path.join(HERE, "system/topoSetDict")).read().split("FoamFile")[0] + '''FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      createPatchDict;
}

pointSync false;

patches
(
    {
        name inlet;
        patchInfo { type patch; }
        constructFrom set;
        set inletFaces;
    }
    {
        name drain;
        patchInfo { type patch; }
        constructFrom set;
        set drainFaces;
    }
);
''')

# no initial liquid: whole enclosure filled with air
rw("system/water/setFieldsDict", lambda s: re.sub(r"regions\s*\(.*?\n\);", "regions\n(\n);", s, flags=re.S))

rw("system/controlDict", lambda s: re.sub(r"\nwriteInterval\s+[^;]+;", "\nwriteInterval   0.25;",
                                          re.sub(r"\nendTime\s+[^;]+;", "\nendTime         3;",
                                                 re.sub(r"\nmaxDeltaT\s+[^;]+;", "\nmaxDeltaT       2e-3;", s))))

# ---------------------------------------------------------------- gas phase = air
def air(s):
    s = re.sub(r"molWeight\s+[\d.]+;", "molWeight       28.96;", s)
    s = re.sub(r"Cp\s+[^;]+;[^\n]*", "Cp          1007;", s)
    s = re.sub(r"mu\s+[^;]+;", "mu          1.81e-5;", s)
    s = re.sub(r"Pr\s+[^;]+;", "Pr          0.71;", s)
    return s


rw("constant/water/thermophysicalProperties.gas", air)

# ---------------------------------------------------------------- boundary conditions
def setT(s, val):
    return re.sub(r"uniform 3\d\d(?:\.\d+)?", f"uniform {val}", s)


extra = {   # field: (inlet, drain)
    "T.gas": (f"type fixedValue; value uniform {T_COOL};",
              f"type inletOutlet; inletValue uniform {T_AIR}; value uniform {T_AIR};"),
    "T.liquid": (f"type fixedValue; value uniform {T_COOL};",
                 f"type inletOutlet; inletValue uniform {T_COOL}; value uniform {T_COOL};"),
    "U.gas": (f"type fixedValue; value uniform (0 {-U_IN} 0);",
              "type pressureInletOutletVelocity; phi phi.gas; inletValue uniform (0 0 0); value uniform (0 0 0);"),
    "U.liquid": (f"type fixedValue; value uniform (0 {-U_IN} 0);",
                 "type pressureInletOutletVelocity; phi phi.liquid; inletValue uniform (0 0 0); value uniform (0 0 0);"),
    "alpha.gas": ("type fixedValue; value uniform 0;",
                  "type inletOutlet; phi phi.gas; inletValue uniform 1; value uniform 1;"),
    "alpha.liquid": ("type fixedValue; value uniform 1;", "type calculated; value uniform 0;"),
    "alphat.gas": ("type calculated; value uniform 1e-8;", "type calculated; value uniform 1e-8;"),
    "alphat.liquid": ("type calculated; value uniform 1e-8;", "type calculated; value uniform 1e-8;"),
    "nut.gas": ("type calculated; value uniform 0;", "type calculated; value uniform 0;"),
    "nut.liquid": ("type calculated; value uniform 0;", "type calculated; value uniform 0;"),
    "k.gas": ("type fixedValue; value uniform 3.75e-3;", "type inletOutlet; inletValue $internalField; value $internalField;"),
    "k.liquid": ("type fixedValue; value uniform 3.75e-3;", "type inletOutlet; inletValue $internalField; value $internalField;"),
    "epsilon.gas": ("type fixedValue; value uniform 1e-2;", "type inletOutlet; inletValue $internalField; value $internalField;"),
    "epsilon.liquid": ("type fixedValue; value uniform 1e-2;", "type inletOutlet; inletValue $internalField; value $internalField;"),
    "f.gas.bubbles": ("type fixedValue; value uniform 1;", "type inletOutlet; inletValue uniform 1; value uniform 1;"),
    "p_rgh": ("type fixedFluxPressure; value $internalField;",
              "type prghTotalPressure; p $internalField; p0 $internalField; U U.gas; value $internalField;"),
}
for fld, (inl, drn) in extra.items():
    def add(s, inl=inl, drn=drn, fld=fld):
        if fld.startswith("T."):
            s = setT(s, T_AIR)
            s = s.replace("internalField   uniform 333;", f"internalField   uniform {T_AIR};")
        if fld == "alpha.gas":
            s = re.sub(r"internalField\s+uniform\s+\S+;", "internalField   uniform 1;", s)
        if fld == "alpha.liquid":
            s = re.sub(r"internalField\s+uniform\s+\S+;", "internalField   uniform 0;", s)
        k = s.rindex("}")
        return s[:k] + f"    inlet\n    {{\n        {inl}\n    }}\n    drain\n    {{\n        {drn}\n    }}\n" + s[k:]
    rw(f"0.orig/water/{fld}", add)
rw("0.orig/water/T.liquid", lambda s: s.replace("uniform 1123", f"uniform {T_AIR}"))
rw("0.orig/water/T.gas", lambda s: s.replace("uniform 1123", f"uniform {T_AIR}"))
rw("0.orig/solid/T", lambda s: s.replace("internalField   uniform 1123;", f"internalField   uniform {T_AIR};")
   .replace("value           uniform 1123;", f"value           uniform {T_AIR};")
   .replace("    maxY\n    {\n        type            zeroGradient;\n    }",
            f"    minY\n    {{\n        type            fixedValue;\n        value           uniform {T_AIR};\n    }}"))

with open(os.path.join(HERE, "Allrun"), "w") as f:
    f.write('''#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
runApplication blockMesh
runApplication topoSet
runApplication createPatch -overwrite
restore0Dir
runApplication splitMeshRegions -cellZonesOnly -overwrite
runApplication setFields -region water
runApplication checkMesh -region water
runApplication $(getApplication)
touch post.foam
''')
os.chmod(os.path.join(HERE, "Allrun"), 0o755)
print("ok")
