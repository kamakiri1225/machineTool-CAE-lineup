"""CalculiX の .frd から固有振動数とモード形状（質量正規化）を読む。刃先・ワークの相対変位ベクトル d_r を作る。"""
import json
import numpy as np


def read_frd(path, nodes=None, modes=None):
    """返り値: freq (nm,), U {node: (3, nm)}。nodes=None なら全節点、modes=None なら全モード"""
    want = None if nodes is None else set(int(n) for n in nodes)
    freq, blocks = [], []
    with open(path) as s:
        for line in s:
            if not line.startswith("  100CL"):
                continue
            f = float(line[12:24]); name = next(s)[5:13].strip()
            if name != "DISP":
                continue
            j = len(freq); freq.append(f); take = modes is None or j in modes; d = {}
            for row in s:
                if row.startswith(" -4") or row.startswith(" -5"):
                    continue
                if row.startswith(" -3"):
                    break
                if take and row.startswith(" -1"):
                    n = int(row[3:13])
                    if want is None or n in want:
                        d[n] = (float(row[13:25]), float(row[25:37]), float(row[37:49]))
            blocks.append(d)
    return np.array(freq), blocks


def contact_vectors(case_dir):
    """刃先・ワークの接触面の節点を剛体フィット → 接触点の並進 (3, nm)。d = ワーク − 刃先"""
    info = json.loads((case_dir / "model_info.json").read_text())
    c = info["contact"]; nodes = c["tool_nodes"] + c["work_nodes"]
    freq, blocks = read_frd(case_dir / "model.frd", nodes)
    coords = {}
    for line in open(case_dir / "model.inp"):
        if line.startswith("*ELEMENT"):
            break
        p = line.split(",")
        if len(p) == 4 and p[0].strip().isdigit() and int(p[0]) in set(nodes):
            coords[int(p[0])] = np.array(list(map(float, p[1:4])))
    out = {}
    for side, ns in (("tool", c["tool_nodes"]), ("work", c["work_nodes"])):
        A = []
        for n in ns:
            rx, ry, rz = coords[n] - np.array(c["xyz"])
            A += [[1, 0, 0, 0, rz, -ry], [0, 1, 0, -rz, 0, rx], [0, 0, 1, ry, -rx, 0]]
        A = np.array(A); U = np.array([[blocks[j][n] for j in range(len(freq))] for n in ns])   # (nn, nm, 3)
        U = U.transpose(0, 2, 1).reshape(-1, len(freq)); sol = np.linalg.lstsq(A, U, rcond=None)[0]
        # 接触面は剛体にしていない（実際の面のまま）。剛体フィットの残差（モードごと）を記録する
        out[side] = sol[:3]; out[side + "_fit_residual"] = np.linalg.norm(A @ sol - U, axis=0) / np.linalg.norm(U, axis=0)
    return freq, out["work"] - out["tool"], out
