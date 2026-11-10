#!/usr/bin/env python3
"""Make the STEADY 3-region CHT case (gas + plate + nozzle, OpenFOAM 13 foamMultiRun) from a transient case made by
tools/make_nozzle_temperature_case.py (run011 type).  This is exactly how run013 was made from run011:
  * gas (air)          : ddtSchemes localEuler (LTS) -> pseudo-steady
  * plate, nozzle      : ddtSchemes steadyState
  * controlDict        : deltaT 1 (one LTS iteration), endTime 12000, write only the final iteration
  * inlet total pressure ramp: over 1500 iterations instead of 0.2 ms
usage: python3 make_nozzle_steady_case.py <transientCase> <newSteadyCase> [iterations]
       (the transient case must already have its mesh: constant/<region>/polyMesh and 0/<region>)"""
import os, re, shutil, sys

src, dst = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
N = int(sys.argv[3]) if len(sys.argv) > 3 else 12000
if os.path.exists(dst):
    sys.exit(f"{dst} exists; choose a new name (this script never overwrites)")
os.makedirs(dst)
for d in ("constant", "system", "0.orig"):
    shutil.copytree(os.path.join(src, d), os.path.join(dst, d))
shutil.copytree(os.path.join(src, "0"), os.path.join(dst, "0"))
for f in ("caseParams.json", "probe_locations.csv"):
    if os.path.exists(os.path.join(src, f)):
        shutil.copy(os.path.join(src, f), dst)
for junk in ("0/nozzle/wallHeatFlux", "0/uniform"):
    p = os.path.join(dst, junk)
    shutil.rmtree(p) if os.path.isdir(p) else (os.remove(p) if os.path.exists(p) else None)


def edit(rel, fn):
    p = os.path.join(dst, rel); s = open(p).read(); open(p, "w").write(fn(s))


def control(s):
    s = re.sub(r"\nstopAt\s+\w+;", "\nstopAt          endTime;", s)
    s = re.sub(r"\nendTime\s+[^;]+;", f"\nendTime         {N};      // number of LTS iterations", s)
    s = re.sub(r"\ndeltaT\s+[^;]+;", "\ndeltaT          1;", s)
    s = re.sub(r"\nwriteControl\s+\w+;", "\nwriteControl    timeStep;", s, count=1)
    s = re.sub(r"\nwriteInterval\s+[^;]+;", f"\nwriteInterval   {N};      // steady: keep only the final iteration", s, count=1)
    s = re.sub(r"\nadjustTimeStep\s+\w+;", "\nadjustTimeStep  no;", s)
    s = re.sub(r"\ntimePrecision\s+\d+;", "\ntimePrecision   6;", s)
    s = re.sub(r"\nwriteFormat\s+\w+;", "\nwriteFormat     ascii;", s)
    # monitors: every 20 iterations; the wall heat flux FIELD only at the final write (else one folder per 20 iterations)
    s = s.replace("writeControl    runTime;\n        writeInterval   1e-5;", "writeControl    timeStep;\n        writeInterval   20;")
    i = s.find("heatNozzle")
    if i >= 0:
        j = s.index("}", i); b = s[i:j]
        b = re.sub(r"writeControl\s+\w+;", "writeControl    writeTime;", b); b = re.sub(r"\n\s*writeInterval\s+[^;]+;", "", b)
        s = s[:i] + b + s[j:]
    return s


edit("system/controlDict", control)
edit("system/air/fvSchemes", lambda s: re.sub(r"(ddtSchemes\s*\{\s*default\s+)\w+;", r"\1localEuler;", s))
for r in ("plate", "nozzle"):
    edit(f"system/{r}/fvSchemes", lambda s: re.sub(r"(ddtSchemes\s*\{\s*default\s+)\w+;", r"\1steadyState;", s))
for f in ("0/air/p", "0.orig/air/p"):
    if os.path.exists(os.path.join(dst, f)):
        edit(f, lambda s: re.sub(r"values\s*\(\(0 101325\) \([0-9.e-]+ ([0-9.]+)\) \(1 [0-9.]+\)\);",
                                 r"values ((0 101325) (1500 \1) (1e9 \1));   // ramp over 1500 LTS iterations", s))
open(os.path.join(dst, "post.foam"), "w").close()
print("steady case ready:", dst)
print("next: source /opt/openfoam13/etc/bashrc ParaView_TYPE=none; cd", dst,
      "; decomposePar -allRegions -force; mpirun -np 4 foamMultiRun -parallel > log.foamMultiRun 2>&1")
