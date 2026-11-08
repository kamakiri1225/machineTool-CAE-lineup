"""モード形状の図（外表面を変形させて色付け）"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "study", "scripts"))
from figstyle import plt  # 日本語フォント・TeX 数式の共通スタイル
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm, colors


def exterior_faces(parts):
    faces = []
    loc = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    for p in parts:
        cnt = {}
        for e in p.elems:
            for f in loc:
                q = tuple(e[i] for i in f)
                k = tuple(sorted(q))
                cnt[k] = (cnt.get(k, (0, q))[0] + 1, q)
        faces += [q for n, q in cnt.values() if n == 1]
    return faces


def plot_modes(model, coords, disp_list, titles, path, scale=120.0, ncol=3):
    faces = exterior_faces(model.parts)
    n = len(disp_list); nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(5.4 * ncol, 4.8 * nrow))
    for k, (d, ttl) in enumerate(zip(disp_list, titles)):
        ax = fig.add_subplot(nrow, ncol, k + 1, projection="3d")
        mx = max(np.linalg.norm(v) for v in d.values())
        polys, vals = [], []
        for q in faces:
            P = np.array([np.array(coords[i]) + scale * np.array(d[i]) / mx for i in q])
            polys.append(P); vals.append(np.mean([np.linalg.norm(d[i]) / mx for i in q]))
        pc = Poly3DCollection(polys, facecolors=cm.turbo(np.array(vals)), edgecolors="k", linewidths=0.05)
        ax.add_collection3d(pc)
        ax.set_xlim(-100, 1100); ax.set_ylim(-100, 1000); ax.set_zlim(0, 1100)
        ax.set_box_aspect((1.2, 1.1, 1.1)); ax.view_init(elev=22, azim=-60)
        ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
        ax.tick_params(labelsize=7)
        ax.set_title(ttl, fontsize=12)
    sm = cm.ScalarMappable(norm=colors.Normalize(0, 1), cmap="turbo")
    fig.colorbar(sm, ax=fig.axes, shrink=0.5, label="変位の大きさ（最大を 1 に正規化）")
    fig.savefig(path); plt.close(fig)
