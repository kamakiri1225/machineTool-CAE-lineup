"""CalculiX の .frd（ASCII）から節点座標と変位（モードごと）を読む簡易リーダ"""
import numpy as np


def read_frd(path):
    coords, blocks = {}, []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("    2C"):
            i += 1
            while not lines[i].startswith(" -3"):
                t = lines[i]
                coords[int(t[3:13])] = (float(t[13:25]), float(t[25:37]), float(t[37:49])); i += 1
        elif ln.startswith("    1PSTEP"):
            pass
        elif ln.startswith("  100CL"):
            hdr = ln
            freq = float(hdr[12:24])
            i += 1
            name = lines[i][5:13].strip()
            data = {}
            while not lines[i].startswith(" -1"):
                i += 1
            while lines[i].startswith(" -1"):
                t = lines[i]
                data[int(t[3:13])] = (float(t[13:25]), float(t[25:37]), float(t[37:49])); i += 1
            blocks.append((name, freq, data))
            continue
        i += 1
    return coords, blocks
