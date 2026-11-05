#!/usr/bin/env python3
"""Build (and optionally run) a compressible multi-species gas-jet case for OpenFOAM v2512 from a YAML file.

  Solver : reactingFoam with chemistry off (mixture of inert gases, e.g. O2/N2/Ar/CO2/He),
           hePsiThermo + constant cp, LTS (localEuler) pseudo-steady, transonic pressure equation.
  YAML   : mesh source, gas list, ambient state, and one entry per boundary patch:
           inlet (total pressure or mass flow, composition in mole% or mass%), outlet, symmetry, wall.
           Patches not listed in the YAML become walls ("default: wall").

usage:
  python3 make_case_from_yaml.py case.yaml                 # mesh + boundary + 0/ constant/ system/
  python3 make_case_from_yaml.py case.yaml --run           # ... and run (serial or mpirun -np N)
  python3 make_case_from_yaml.py case.yaml --check         # only print the resolved settings

needs: OpenFOAM v2512 sourced in the shell (blockMesh etc. on PATH), python3, PyYAML (pip install pyyaml)
"""
import argparse, os, re, shutil, subprocess, sys
import yaml

# --------------------------------------------------------------------------- gas data (constant cp, ~300 K)
GAS = {  #        W [g/mol], cp [J/kgK], mu [Pa s], Pr
    "N2":  (28.0134, 1040.0, 1.79e-5, 0.71),
    "O2":  (31.9988,  918.0, 2.06e-5, 0.72),
    "Ar":  (39.948,   520.3, 2.27e-5, 0.67),
    "CO2": (44.0095,  846.0, 1.48e-5, 0.77),
    "He":  (4.0026,  5193.0, 1.99e-5, 0.68),
    "H2O": (18.0153, 1864.0, 1.00e-5, 0.95),   # vapour, only as a minor species
}
PATM = 101325.0
HDR = 'FoamFile\n{{\n    version 2.0;\n    format ascii;\n    class {cls};\n    location "{loc}";\n    object {obj};\n}}\n\n'


def die(msg):
    sys.exit(f"[error] {msg}")


def write(case, rel, body, cls="dictionary"):
    p = os.path.join(case, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(HDR.format(cls=cls, loc=os.path.dirname(rel), obj=os.path.basename(rel)) + body)


def sh(case, cmd, log):
    print(f"  $ {cmd}   (log.{log})")
    with open(os.path.join(case, f"log.{log}"), "w") as f:
        r = subprocess.run(cmd, shell=True, cwd=case, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        die(f"'{cmd}' failed, see {case}/log.{log}")


# --------------------------------------------------------------------------- composition
def mass_fractions(comp, basis, species):
    """comp: {gas: percent or fraction}, basis 'mole' or 'mass'  ->  {gas: mass fraction} over all species."""
    for k in comp:
        if k not in species:
            die(f"gas '{k}' in a composition is not in 'gases: {species}'")
    tot = sum(comp.values())
    if tot <= 0:
        die(f"composition {comp} sums to zero")
    x = {k: v / tot for k, v in comp.items()}
    if basis == "mole":
        m = {k: x[k] * GAS[k][0] for k in x}
        s = sum(m.values()); x = {k: m[k] / s for k in m}
    elif basis != "mass":
        die(f"basis must be 'mole' or 'mass', got '{basis}'")
    return {k: x.get(k, 0.0) for k in species}


# --------------------------------------------------------------------------- mesh
def make_mesh(case, M, base):
    t = M.get("type", "blockMesh")
    src = lambda k: os.path.abspath(os.path.join(base, M[k]))
    os.makedirs(os.path.join(case, "constant"), exist_ok=True)
    if t == "blockMesh":
        shutil.copy(src("blockMeshDict"), os.path.join(case, "system/blockMeshDict"))
        sh(case, "blockMesh", "blockMesh")
    elif t == "snappyHexMesh":                       # a directory with blockMeshDict, snappyHexMeshDict, (meshQualityDict), triSurface/
        d = src("dir")
        for f in ("blockMeshDict", "snappyHexMeshDict", "meshQualityDict", "surfaceFeatureExtractDict"):
            if os.path.exists(os.path.join(d, f)):
                shutil.copy(os.path.join(d, f), os.path.join(case, "system", f))
        if os.path.isdir(os.path.join(d, "triSurface")):
            shutil.copytree(os.path.join(d, "triSurface"), os.path.join(case, "constant/triSurface"), dirs_exist_ok=True)
        sh(case, "blockMesh", "blockMesh")
        if os.path.exists(os.path.join(case, "system/surfaceFeatureExtractDict")):
            sh(case, "surfaceFeatureExtract", "surfaceFeatureExtract")
        sh(case, "snappyHexMesh -overwrite", "snappyHexMesh")
    elif t == "polyMesh":                            # an existing constant/polyMesh directory
        shutil.copytree(src("dir"), os.path.join(case, "constant/polyMesh"), dirs_exist_ok=True)
    elif t in ("fluent", "gmsh", "unv"):
        conv = {"fluent": "fluentMeshToFoam", "gmsh": "gmshToFoam", "unv": "ideasUnvToFoam"}[t]
        sh(case, f'{conv} "{src("file")}"', conv)
    else:
        die(f"mesh.type '{t}' unknown (blockMesh, snappyHexMesh, polyMesh, fluent, gmsh, unv)")
    s = float(M.get("scale", 1.0))
    if s != 1.0:
        sh(case, f"transformPoints -scale '({s} {s} {s})'", "transformPoints")


def read_patches(case):
    txt = open(os.path.join(case, "constant/polyMesh/boundary")).read()
    body = txt[txt.index("(", txt.index("}")) + 1:txt.rindex(")")]
    return [(m.group(1), m.group(2)) for m in re.finditer(r"(\S+)\s*\{[^}]*?type\s+(\w+);", body)]


def set_patch_types(case, types):
    """rewrite 'type' (and drop inGroups) of each patch in constant/polyMesh/boundary."""
    p = os.path.join(case, "constant/polyMesh/boundary")
    txt = open(p).read()
    for name, typ in types.items():
        pat = re.compile(r"(\n\s*" + re.escape(name) + r"\s*\{)([^}]*)\}")
        m = pat.search(txt)
        blk = re.sub(r"\n\s*inGroups[^;]*;", "", m.group(2))
        blk = re.sub(r"type\s+\w+;", f"type            {typ};", blk)
        txt = txt[:m.start()] + m.group(1) + blk + "}" + txt[m.end():]
    open(p, "w").write(txt)


# --------------------------------------------------------------------------- settings resolution
def resolve(C, patches):
    species = C["gases"]
    for g in species:
        if g not in GAS:
            die(f"gas '{g}' not in the built-in table {list(GAS)}")
    A = C.get("ambient", {})
    amb = dict(p=float(A.get("p", PATM)), T=float(A.get("T", 300.0)),
               Y=mass_fractions(A.get("composition", {"N2": 79, "O2": 21}), A.get("basis", "mole"), species))
    default = C.get("default", "wall")
    names = [n for n, _ in patches]
    B = C.get("boundaries", {})
    out = {}
    for pat_name, spec in B.items():                       # keys may be regular expressions, e.g. "inlet_.*"
        hit = [n for n in names if re.fullmatch(pat_name, n)]
        if not hit:
            die(f"boundary '{pat_name}' in YAML matches no mesh patch. Mesh patches: {names}")
        for n in hit:
            out[n] = dict(spec)
    for n, t in patches:
        if n not in out:
            if t.startswith("cyclic") or t == "processor":
                die(f"patch '{n}' is {t}; cyclic/processor patches are not supported")
            out[n] = {"type": default}                      # everything not listed (incl. blockMesh 'defaultFaces') -> default
    for n, s in out.items():
        t = s["type"]
        if t == "inlet":
            s["Y"] = mass_fractions(s["composition"], s.get("basis", "mole"), species)
            if "p_total_gauge" in s:
                s["p0"] = amb["p"] + float(s["p_total_gauge"])
            elif "p_total" in s:
                s["p0"] = float(s["p_total"])
            elif "mass_flow" not in s:
                die(f"inlet '{n}' needs p_total_gauge, p_total (absolute) or mass_flow [kg/s]")
            s["T0"] = float(s.get("T_total", s.get("T", amb["T"])))
        elif t == "outlet":
            s["p"] = float(s.get("p", amb["p"]))
            s["T_back"] = float(s.get("T_backflow", amb["T"]))
            s["Y_back"] = (mass_fractions(s["composition_backflow"], s.get("basis", "mole"), species)
                           if "composition_backflow" in s else amb["Y"])
        elif t not in ("wall", "symmetryPlane", "symmetry", "wedge", "empty"):
            die(f"patch '{n}': type '{t}' unknown (inlet, outlet, wall, symmetryPlane, symmetry, wedge, empty)")
    return species, amb, out


# --------------------------------------------------------------------------- field files
def bc_block(name, d):
    inner = "".join(f"        {k:<16} {v};\n" for k, v in d.items())
    return f"    {name}\n    {{\n{inner}    }}\n"


def field(case, name, cls, dim, internal, bcs):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    body += "".join(bc_block(n, d) for n, d in bcs.items()) + "}\n"
    write(case, f"0/{name}", body, cls)


def write_fields(case, species, amb, P, C):
    turb = C.get("turbulence", {}).get("model", "realizableKE")
    om = turb in ("kOmegaSST", "kOmega")
    k0, eps0 = 1.0, 1e4
    S, V = "volScalarField", "volVectorField"
    ramp = int(C.get("solver", {}).get("pressure_ramp_iterations", 500))
    f = {k: {} for k in ["p", "U", "T", "k", "eps", "nut", "alphat"] + species}
    for n, s in P.items():
        t = s["type"]
        if t in ("symmetryPlane", "symmetry", "wedge", "empty"):
            for k in f:
                f[k][n] = {"type": t}
            continue
        if t == "wall":
            f["p"][n] = {"type": "zeroGradient"}
            f["U"][n] = {"type": "noSlip"}
            f["T"][n] = ({"type": "fixedValue", "value": f"uniform {float(s['T'])}"} if "T" in s else {"type": "zeroGradient"})
            for g in species:
                f[g][n] = {"type": "zeroGradient"}
            f["k"][n] = {"type": "kqRWallFunction", "value": f"uniform {k0}"}
            f["eps"][n] = ({"type": "omegaWallFunction", "value": f"uniform {eps0/(0.09*k0)}"} if om
                           else {"type": "epsilonWallFunction", "value": f"uniform {eps0}"})
            f["nut"][n] = {"type": "nutkWallFunction", "value": "uniform 0"}
            f["alphat"][n] = {"type": "compressible::alphatWallFunction", "value": "uniform 0"}
        elif t == "inlet":
            if "p0" in s:
                f["p"][n] = {"type": "uniformTotalPressure",
                             "p0": f"table ((0 {amb['p']}) ({ramp} {s['p0']}) (1e9 {s['p0']}))" if ramp > 0 else f"constant {s['p0']}",
                             "psi": "thermo:psi", "gamma": 1.4,      # compressible (isentropic) total pressure
                             "value": f"uniform {amb['p']}"}
                f["U"][n] = {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"}
                f["T"][n] = {"type": "totalTemperature", "gamma": 1.4, "T0": f"uniform {s['T0']}", "value": f"uniform {s['T0']}"}
            else:
                f["p"][n] = {"type": "zeroGradient"}
                f["U"][n] = {"type": "flowRateInletVelocity", "massFlowRate": f"constant {float(s['mass_flow'])}",
                             "rho": "rho", "rhoInlet": 1.0, "value": "uniform (0 0 0)"}
                f["T"][n] = {"type": "fixedValue", "value": f"uniform {s['T0']}"}
            for g in species:
                f[g][n] = {"type": "fixedValue", "value": f"uniform {s['Y'][g]:.6g}"}
            I, L = float(s.get("turbulent_intensity", 0.05)), float(s.get("mixing_length", 1e-4))
            f["k"][n] = {"type": "turbulentIntensityKineticEnergyInlet", "intensity": I, "value": f"uniform {k0}"}
            f["eps"][n] = ({"type": "turbulentMixingLengthFrequencyInlet", "mixingLength": L, "value": f"uniform {eps0/(0.09*k0)}"} if om
                           else {"type": "turbulentMixingLengthDissipationRateInlet", "mixingLength": L, "value": f"uniform {eps0}"})
            f["nut"][n] = {"type": "calculated", "value": "uniform 0"}
            f["alphat"][n] = {"type": "calculated", "value": "uniform 0"}
        elif t == "outlet":
            f["p"][n] = ({"type": "fixedValue", "value": f"uniform {s['p']}"} if s.get("pressure_type", "total") == "static"
                         else {"type": "totalPressure", "p0": f"uniform {s['p']}", "value": f"uniform {s['p']}"})
            f["U"][n] = {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"}
            f["T"][n] = {"type": "inletOutlet", "inletValue": f"uniform {s['T_back']}", "value": f"uniform {s['T_back']}"}
            for g in species:
                f[g][n] = {"type": "inletOutlet", "inletValue": f"uniform {s['Y_back'][g]:.6g}", "value": f"uniform {s['Y_back'][g]:.6g}"}
            f["k"][n] = {"type": "inletOutlet", "inletValue": f"uniform {k0}", "value": f"uniform {k0}"}
            e = eps0 / (0.09 * k0) if om else eps0
            f["eps"][n] = {"type": "inletOutlet", "inletValue": f"uniform {e}", "value": f"uniform {e}"}
            f["nut"][n] = {"type": "calculated", "value": "uniform 0"}
            f["alphat"][n] = {"type": "calculated", "value": "uniform 0"}
    field(case, "p", S, "[1 -1 -2 0 0 0 0]", f"uniform {amb['p']}", f["p"])
    field(case, "U", V, "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", f["U"])
    field(case, "T", S, "[0 0 0 1 0 0 0]", f"uniform {amb['T']}", f["T"])
    field(case, "k", S, "[0 2 -2 0 0 0 0]", f"uniform {k0}", f["k"])
    if om:
        field(case, "omega", S, "[0 0 -1 0 0 0 0]", f"uniform {eps0/(0.09*k0)}", f["eps"])
    else:
        field(case, "epsilon", S, "[0 2 -3 0 0 0 0]", f"uniform {eps0}", f["eps"])
    field(case, "nut", S, "[0 2 -1 0 0 0 0]", "uniform 0", f["nut"])
    field(case, "alphat", S, "[1 -1 -1 0 0 0 0]", "uniform 0", f["alphat"])
    for g in species:
        field(case, g, S, "[0 0 0 0 0 0 0]", f"uniform {amb['Y'][g]:.6g}", f[g])
    ydef = {n: ({"type": d["type"]} if d["type"] in ("symmetryPlane", "symmetry", "wedge", "empty")
                else {"type": "zeroGradient"}) for n, d in P.items()}
    field(case, "Ydefault", S, "[0 0 0 0 0 0 0]", "uniform 0", ydef)


# --------------------------------------------------------------------------- constant / system
def write_dicts(case, species, P, C):
    SV = C.get("solver", {})
    turb = C.get("turbulence", {}).get("model", "realizableKE")
    inert = "N2" if "N2" in species else species[0]
    write(case, "constant/thermophysicalProperties",
          "thermoType\n{\n    type hePsiThermo;\n    mixture reactingMixture;\n    transport const;\n    thermo hConst;\n"
          "    energy sensibleEnthalpy;\n    equationOfState perfectGas;\n    specie specie;\n}\n\n"
          f"inertSpecie {inert};\nchemistryReader foamChemistryReader;\n"
          'foamChemistryFile "<constant>/reactions";\nfoamChemistryThermoFile "<constant>/thermo.hConst";\n')
    write(case, "constant/reactions", "species\n(\n" + "".join(f"    {g}\n" for g in species) + ");\n\nreactions\n{\n}\n")
    write(case, "constant/thermo.hConst", "".join(
        f"{g}\n{{\n    specie {{ molWeight {GAS[g][0]}; }}\n    thermodynamics {{ Cp {GAS[g][1]}; Hf 0; }}\n"
        f"    transport {{ mu {GAS[g][2]}; Pr {GAS[g][3]}; }}\n}}\n" for g in species))
    write(case, "constant/combustionProperties", "combustionModel none;\n")
    write(case, "constant/chemistryProperties", "chemistryType { solver none; }\nchemistry off;\ninitialChemicalTimeStep 1e-07;\n")
    write(case, "constant/turbulenceProperties",
          "simulationType laminar;\n" if turb == "laminar" else
          f"simulationType RAS;\nRAS {{ RASModel {turb}; turbulence on; printCoeffs on; }}\n")
    it, wi = int(SV.get("iterations", 4000)), int(SV.get("write_interval", 500))
    flows = "".join(f"    flow_{n} {{ type surfaceFieldValue; libs (fieldFunctionObjects); regionType patch; name {n}; "
                    f"operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 20; }}\n"
                    for n, s in P.items() if s["type"] in ("inlet", "outlet"))
    write(case, "system/controlDict",
          f"application reactingFoam;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\n"
          f"endTime {it};   // LTS: number of pseudo-time iterations\ndeltaT 1;\nwriteControl timeStep;\n"
          f"writeInterval {wi};\npurgeWrite 0;\nwriteFormat ascii;\nwritePrecision 8;\ntimeFormat general;\n"
          f"timePrecision 6;\nrunTimeModifiable true;\n\nfunctions\n{{\n"
          "    MachNo { type MachNo; libs (fieldFunctionObjects); executeControl writeTime; writeControl writeTime; }\n"
          + flows + "}\n")
    k = "omega" if turb in ("kOmegaSST", "kOmega") else "epsilon"
    write(case, "system/fvSchemes",
          "ddtSchemes { default localEuler; }\n"
          "gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }\n"
          "divSchemes\n{\n    default none;\n    div(phi,U) Gauss linearUpwindV grad(U);\n    div(phi,Yi_h) Gauss upwind;\n"
          f"    div(phi,K) Gauss upwind;\n    div(phid,p) Gauss upwind;\n    div(phi,k) Gauss upwind;\n    div(phi,{k}) Gauss upwind;\n"
          "    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;\n}\n"
          "laplacianSchemes { default Gauss linear limited corrected 0.5; }\ninterpolationSchemes { default linear; }\n"
          "snGradSchemes { default limited corrected 0.5; }\nwallDist { method meshWave; }\n")
    write(case, "system/fvSolution",
          'solvers\n{\n    "rho.*" { solver diagonal; }\n'
          '    "p.*" { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.05; }\n'
          f'    "(U|h|k|{k}|Yi).*" {{ solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.1; }}\n}}\n\n'
          "PIMPLE\n{\n    momentumPredictor yes;\n    transonic yes;\n    nOuterCorrectors 1;\n    nCorrectors 2;\n"
          f"    nNonOrthogonalCorrectors {int(SV.get('non_orthogonal_correctors', 1))};\n"
          f"    maxCo {float(SV.get('maxCo', 0.2))};\n    maxDeltaT 1e-6;\n    alphaTemp 0.05;\n    alphaY 0.05;\n"
          f"    Yref {{ {species[0]} 0.1; }}\n    rDeltaTSmoothingCoeff 0.02;\n    rDeltaTDampingCoeff 1;\n}}\n\n"
          'relaxationFactors { equations { ".*" 1; } }\n')
    lo, hi = SV.get("T_limits", [50, 3000])
    write(case, "constant/fvOptions",
          f"limitT\n{{\n    type limitTemperature;\n    active yes;\n    selectionMode all;\n    min {lo};\n    max {hi};\n}}\n")
    write(case, "system/decomposeParDict", f"numberOfSubdomains {int(C.get('parallel', 1))};\nmethod scotch;\n")


def write_allrun(case, C):
    n = int(C.get("parallel", 1))
    run = (f"runApplication decomposePar -force\nrunParallel reactingFoam\nrunApplication reconstructPar -latestTime\n"
           if n > 1 else "runApplication reactingFoam\n")
    with open(os.path.join(case, "Allrun"), "w") as f:
        f.write('#!/bin/sh\ncd "${0%/*}" || exit\n. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n'
                "export HWLOC_COMPONENTS=-gl   # WSL2: keeps mpirun from hanging\nrm -f log.reactingFoam log.decomposePar\n"
                + run + "touch post.foam\n")
    os.chmod(os.path.join(case, "Allrun"), 0o755)


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("yaml"); ap.add_argument("--run", action="store_true", help="run the solver after set-up")
    ap.add_argument("--check", action="store_true", help="print the resolved settings only (needs an existing mesh)")
    ap.add_argument("--force", action="store_true", help="overwrite an existing case directory")
    a = ap.parse_args()
    base = os.path.dirname(os.path.abspath(a.yaml))
    C = yaml.safe_load(open(a.yaml, encoding="utf-8"))
    case = os.path.abspath(os.path.join(base, C["case"]))
    if shutil.which("blockMesh") is None:
        die("OpenFOAM is not sourced. Run e.g.  source /usr/lib/openfoam/openfoam2512/etc/bashrc")
    if not a.check:
        if os.path.exists(case):
            if not a.force:
                die(f"{case} exists. Use --force to overwrite it.")
            shutil.rmtree(case)
        os.makedirs(os.path.join(case, "system"))
        write(case, "system/controlDict", "application reactingFoam;\nstartTime 0;\nendTime 1;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n")
        write(case, "system/fvSchemes", "")
        write(case, "system/fvSolution", "")
        print("[1/4] mesh"); make_mesh(case, C["mesh"], base)
    patches = read_patches(case)
    species, amb, P = resolve(C, patches)
    print(f"gases {species}, ambient p {amb['p']:.0f} Pa, T {amb['T']} K, Y {{{', '.join(f'{k}: {v:.4f}' for k, v in amb['Y'].items())}}}")
    print(f"{'patch':<20}{'type':<15}settings")
    for n, s in P.items():
        extra = ""
        if s["type"] == "inlet":
            extra = (f"p0 {s['p0']/1e6:.4f} MPa(abs)" if "p0" in s else f"mass flow {s['mass_flow']} kg/s") + \
                    f", T0 {s['T0']} K, Y {{{', '.join(f'{k}: {v:.4f}' for k, v in s['Y'].items())}}}"
        elif s["type"] == "outlet":
            extra = f"p {s['p']:.0f} Pa, backflow T {s['T_back']} K"
        elif s["type"] == "wall":
            extra = f"T {s['T']} K" if "T" in s else "adiabatic"
        print(f"{n:<20}{s['type']:<15}{extra}")
    if a.check:
        return
    print("[2/4] boundary types")
    ptype = {"inlet": "patch", "outlet": "patch", "wall": "wall"}
    set_patch_types(case, {n: ptype.get(s["type"], s["type"]) for n, s in P.items()})
    print("[3/4] 0/ constant/ system/")
    write_fields(case, species, amb, P, C)
    write_dicts(case, species, P, C)
    write_allrun(case, C)
    sh(case, "checkMesh", "checkMesh")
    open(os.path.join(case, "post.foam"), "w").close()
    shutil.copy(a.yaml, os.path.join(case, "case.yaml"))
    print(f"[4/4] ready: {case}")
    if a.run:
        print("running ./Allrun ..."); subprocess.run("./Allrun", shell=True, cwd=case, check=True)
        print("done. see log.reactingFoam, postProcessing/flow_*")
    else:
        print(f"next:  cd {case} && ./Allrun")


if __name__ == "__main__":
    main()
