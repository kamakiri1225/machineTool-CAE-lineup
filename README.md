# machineTool-CAE-lineup

A lineup of simple CFD / heat-transfer models around machine tools, built with **OpenFOAM v2512 (ESI)**.
Each case contains the analysis setup (`system/`, `constant/`, `0.orig/`, or a Python case generator in `tools/`),
documents (Japanese) and animations. Heavy results (time directories, meshes, logs) are not included; run the case to regenerate them.

**Browse the documents with figures:** https://kamakiri1225.github.io/machineTool-CAE-lineup/

| Folder | Topic | Solver |
|---|---|---|
| `01_coolantCHT` | Warm coolant on a workpiece / table (gas-liquid Euler-Euler + solid conduction) | chtMultiRegionTwoPhaseEulerFoam |
| `02_laserAssistGas` | Laser-cutting assist gas (compressible, supersonic, O2/N2 mixture; axisymmetric and 3-D kerf) | reactingFoam (LTS, transonic) |
| `03_fumeExtraction` | Cutting fume extraction under the table (steady flow + particle tracking) | simpleFoam + icoUncoupledKinematicParcelFoam |
| `04_factoryAmbient` | Factory ambient temperature and machine temperature with / without a cover | chtMultiRegionFoam (frozenFlow) |
| `05_spindleBearingHeat` | Spindle bearing heat-up and cooling (axisymmetric solid) | solidFoam |
| `06_mistChips` | Coolant mist and chips in the machining enclosure | simpleFoam + icoUncoupledKinematicParcelFoam |

### Compute time (Xeon Gold 6212U, WSL2)

| Case | Cells | Processes | Wall time |
|---|---|---|---|
| 01 coolant run004 (3-D, 3 s) | 18,000 | 1 | 12 min |
| 02 assist gas run002 / run003 (axisymmetric) | 10,525 / 12,270 | 1 | 24 / 42 min |
| 02 assist gas run004 (3-D kerf, 0.5 / 2.0 MPa) | 588,724 | 4 | ~7.5 h each |
| 03 fume extraction (flow + particles) | 55,296 | 1 | 6–11 min + 3 min |
| 04 factory ambient (flow 120 s + heat 24 h) | 26,880 | 1 | 23 + 3 min |
| 05 spindle bearing heat (1 h) | 1,584 | 1 | 18 s |
| 06 mist and chips (flow + 2 particle runs) | 73,728 | 1 | ~25 + 4 + 3 min |

All models use generic dimensions or conditions taken from published papers; they are not the design data of any specific machine.

### Notes
- Parallel runs on WSL2: `export HWLOC_COMPONENTS=-gl` before `mpirun` (otherwise Open MPI hangs at start-up).
- Figures are rendered with python VTK scripts in each `tools/` folder (ParaView `pvbatch` hangs in the same environment).
- Papers cited in the documents are listed with links only (PDFs are not redistributed).
