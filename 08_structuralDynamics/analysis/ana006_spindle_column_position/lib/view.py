"""pyvista 表示用：Model → UnstructuredGrid（表示は 8 節点の六面体で十分）"""
import numpy as np
import pyvista as pv
pv.OFF_SCREEN = True
PART_ID = {"BED": 0, "TABLE": 1, "COLUMN": 2, "HEAD": 3}
PART_COLOR = {"BED": "#b9c0c8", "TABLE": "#8fb8de", "COLUMN": "#e3c08d", "HEAD": "#a7d3a1"}
MAT_COLOR = {"SPINDLE": "#6c7a89", "TOOL": "#444444"}


def grid(m):
    ids = sorted(m.coords); idx = {n: i for i, n in enumerate(ids)}
    pts = np.array([m.coords[n] for n in ids])
    cells = np.hstack([[8] + [idx[n] for n in e[:8]] for _, _, e in m.elems])
    g = pv.UnstructuredGrid(cells, np.full(len(m.elems), pv.CellType.HEXAHEDRON), pts)
    g.cell_data["part"] = [PART_ID[p] for p, _, _ in m.elems]
    g.cell_data["mat"] = [0 if mt == "FRAME" else (1 if mt == "SPINDLE" else 2) for _, mt, _ in m.elems]
    return g, ids


def surface(g, cond):
    return g.extract_cells(np.where(cond)[0]).extract_surface(algorithm="dataset_surface")
