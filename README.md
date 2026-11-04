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
| `04_factory_spindle_mist` | Factory ambient temperature, spindle bearing heat, coolant mist and chips | chtMultiRegionFoam, solidFoam, simpleFoam + particles |

All models use generic dimensions or conditions taken from published papers; they are not the design data of any specific machine.

### Notes
- Parallel runs on WSL2: `export HWLOC_COMPONENTS=-gl` before `mpirun` (otherwise Open MPI hangs at start-up).
- Figures are rendered with python VTK scripts in each `tools/` folder (ParaView `pvbatch` hangs in the same environment).
- Papers cited in the documents are listed with links only (PDFs are not redistributed).
