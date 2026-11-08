"""代表点どうしを「6 成分のばね（並進 kx, ky, kz ＋ 回転 krx, kry, krz）」でつなぐ道具。
   使う人は「面 A の節点・面 B の節点・代表点の座標・6 つのばね定数」を渡すだけでよい。
   中では、ソリッドの節点に回転の自由度がないため、次を自動で作る（CalculiX・FrontISTR 共通）:
     代表点 REF（並進）と回転用補助点 ROT（自由度 1〜3 を θx, θy, θz として使う）を面ごとに 1 組
     面の各節点を REF・ROT に剛体で結ぶ拘束式（RBE2 相当）
     REF どうしに並進ばね 3 本、ROT どうしに回転ばね 3 本
   例:
     from rep_spring import RepSpring
     rs = RepSpring(first_node=101, first_elem=101)
     rs.add("CAR1", faceA, faceB, center=(50, 50, 100), k=(1000, 1000, 1000, 1e6, 1e6, 1e6), coords=coords)
     ccx_lines = rs.ccx()                 # *NODE の追加分・*EQUATION・SPRING2 を .inp に貼る
     msh_lines, cnt_lines = rs.fistr()    # !NODE 追加分・!EQUATION・511 要素（.msh）と !SPRING_D（.cnt）
"""
import numpy as np


class RepSpring:
    def __init__(self, first_node, first_elem):
        self.n, self.e = first_node, first_elem
        self.nodes, self.eqs, self.springs, self.info = [], [], [], []

    def _rigid(self, face, ref, rot, c, coords):
        for n in face:
            rx, ry, rz = np.asarray(coords[n], float) - np.asarray(c, float)
            for eq in ([(n, 1, 1.), (ref, 1, -1.), (rot, 2, -rz), (rot, 3, ry)],
                       [(n, 2, 1.), (ref, 2, -1.), (rot, 3, -rx), (rot, 1, rz)],
                       [(n, 3, 1.), (ref, 3, -1.), (rot, 1, -ry), (rot, 2, rx)]):
                self.eqs.append([t for t in eq if abs(t[2]) > 1e-12])

    def add(self, name, faceA, faceB, center, k, coords):
        """faceA, faceB: 節点番号のリスト。center: 代表点の座標。k: (kx, ky, kz, krx, kry, krz) [N/mm, N·mm/rad]"""
        refA, rotA, refB, rotB = self.n, self.n + 1, self.n + 2, self.n + 3; self.n += 4
        for q in (refA, rotA, refB, rotB): self.nodes.append((q, *center))
        self._rigid(faceA, refA, rotA, center, coords); self._rigid(faceB, refB, rotB, center, coords)
        for d in range(3):
            if k[d] > 0: self.springs.append((f"{name}_T{d+1}", self.e, refA, refB, d + 1, k[d])); self.e += 1
            if k[d + 3] > 0: self.springs.append((f"{name}_R{d+1}", self.e, rotA, rotB, d + 1, k[d + 3])); self.e += 1
        self.info.append(dict(name=name, refA=refA, rotA=rotA, refB=refB, rotB=rotB))
        return self.info[-1]

    def ccx(self):
        L = ["** --- 代表点ばね（rep_spring.py で生成） ---", "*NODE"] + [f"{q}, {x:g}, {y:g}, {z:g}" for q, x, y, z in self.nodes]
        for eq in self.eqs:
            L += ["*EQUATION", str(len(eq)), ", ".join(f"{n}, {d}, {c:g}" for n, d, c in eq)]
        for name, e, a, b, d, k in self.springs:
            L += [f"*ELEMENT, TYPE=SPRING2, ELSET={name}", f"{e}, {a}, {b}", f"*SPRING, ELSET={name}", f"{d}, {d}", f"{float(k):.6E}"]
        return L

    def fistr(self):
        M = ["!NODE"] + [f"{q}, {x:g}, {y:g}, {z:g}" for q, x, y, z in self.nodes] + ["!EQUATION"]
        for eq in self.eqs:
            M += [f"{len(eq)}, 0.0", ", ".join(f"{n}, {d}, {c:g}" for n, d, c in eq)]
        C = []
        for name, e, a, b, d, k in self.springs:
            M += [f"!ELEMENT, TYPE=511, EGRP={name}", f"{e}, {a}, {b}", f"!MATERIAL, NAME={name}, ITEM=1", "!ITEM=1, SUBITEM=2", "1.0, 0.3",
                  f"!SECTION, TYPE=INTERFACE, EGRP={name}, MATERIAL={name}", "1.0"]
            C += [f"!MATERIAL, NAME={name}", "!SPRING_D", f"{d}, {d}, {float(k):.6E}"]
        return M, C
