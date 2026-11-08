"""パラメータ (dy, dz) から、部品ごとの直方体と、部品どうしをつなぐ接合面を作る。

部品は 4 つ。部品どうしはメッシュを共有しない（節点は別々）。接合面ごとに
「代表点＋MPC（剛体結合）」を両側に作り、代表点どうしをばねでつなぐ。
  BED    ベッド（X レール・X ねじ台・Y レール・Y ねじ台を含む）
  TABLE  テーブル（X キャリッジ・ナットブラケット・治具・ワーク）
  COLUMN コラム（Y キャリッジ・Y ナット・Z レール・Z ねじ台）
  HEAD   主軸頭（Z キャリッジ・Z ナット・主軸・テストバー）
"""
from dataclasses import dataclass
import params as P

AXIS = {"x": 0, "y": 1, "z": 2}


@dataclass
class Joint:
    name: str
    part_a: str          # 固定側
    part_b: str          # 動く側
    normal: str          # 接合面の法線方向 'x' 'y' 'z'
    c: float             # 接合面の座標（法線方向）
    rect: dict           # 面内 2 方向の範囲 {'x':(x0,x1), 'y':(y0,y1)} など
    k: dict              # 全体座標の方向(1,2,3) → 並進ばね N/mm


def sh(b, dx=0, dy=0, dz=0):
    return (b[0] + dx, b[1] + dx, b[2] + dy, b[3] + dy, b[4] + dz, b[5] + dz)


def carriage_k(normal, rail_axis):
    lat = ({"x", "y", "z"} - {normal, rail_axis}).pop()
    return {AXIS[normal] + 1: P.K_CAR_NORMAL, AXIS[lat] + 1: P.K_CAR_LATERAL}


def build(dy=0.0, dz=0.0):
    parts = {}
    # ---------------- BED（動かない）
    parts["BED"] = [((0, 1000, 0, 1100, 0, 200), "FRAME"),
                    ((0, 1000, 130, 170, 200, 230), "FRAME"),     # X レール 1
                    ((0, 1000, 430, 470, 200, 230), "FRAME"),     # X レール 2
                    ((460, 540, 280, 320, 200, 230), "FRAME"),    # X ねじ（ナットの真下の台）
                    ((370, 410, 520, 1100, 200, 230), "FRAME"),   # Y レール 1
                    ((590, 630, 520, 1100, 200, 230), "FRAME"),   # Y レール 2
                    ((480, 520, 520, 1100, 200, 230), "FRAME")]   # Y ねじ（帯）
    # ---------------- TABLE（X だけ動く。今回は中央に固定）
    tab = [((250, 750, 100, 500, 270, 350), "FRAME")]
    xcar = [(270, 350, 130, 170), (650, 730, 130, 170), (270, 350, 430, 470), (650, 730, 430, 470)]
    tab += [((x0, x1, y0, y1, 230, 270), "FRAME") for x0, x1, y0, y1 in xcar]
    tab += [((460, 540, 280, 320, 230, 270), "FRAME")]            # ナットブラケット
    yw = 400 + dy                                                   # ワーク（刃先）の Y
    if dz > 0:                                                      # 治具ブロック（刃先が上がった分）
        tab += [((350, 650, max(100, yw - 100), min(500, yw + 100), 350, 350 + dz), "FRAME")]
    tab += [((450, 550, yw - 50, yw + 50, 350 + dz, 400 + dz), "FRAME")]   # ワーク 100×100×50
    parts["TABLE"] = tab
    # ---------------- COLUMN（Y に dy 動く）
    ycar = [(370, 410, 660, 740), (370, 410, 760, 840), (590, 630, 660, 740), (590, 630, 760, 840)]
    col = [(sh((350, 650, 650, 850, 270, 1200), dy=dy), "FRAME")]                       # コラム本体
    col += [(sh((x0, x1, y0, y1, 230, 270), dy=dy), "FRAME") for x0, x1, y0, y1 in ycar]  # Y キャリッジ
    col += [(sh((480, 520, 730, 770, 230, 270), dy=dy), "FRAME")]                        # Y ナット
    col += [(sh((390, 430, 630, 650, 300, 1200), dy=dy), "FRAME"),                      # Z レール 1
            (sh((570, 610, 630, 650, 300, 1200), dy=dy), "FRAME"),                      # Z レール 2
            (sh((480, 520, 630, 650, 300, 1200), dy=dy), "FRAME")]                      # Z ねじ（帯）
    parts["COLUMN"] = col
    # ---------------- HEAD（Y に dy、Z に dz 動く）
    zcar = [(390, 430, 660, 740), (390, 430, 760, 840), (570, 610, 660, 740), (570, 610, 760, 840)]
    hd = [(sh((470, 530, 370, 430, 400, 500), dy=dy, dz=dz), "TOOL"),
          (sh((450, 550, 350, 450, 500, 650), dy=dy, dz=dz), "SPINDLE"),
          (sh((380, 620, 340, 610, 650, 850), dy=dy, dz=dz), "FRAME")]
    hd += [(sh((x0, x1, 610, 630, z0, z1), dy=dy, dz=dz), "FRAME") for x0, x1, z0, z1 in zcar]
    hd += [(sh((480, 520, 610, 630, 720, 780), dy=dy, dz=dz), "FRAME")]
    parts["HEAD"] = hd

    joints = []
    for i, (x0, x1, y0, y1) in enumerate(xcar, 1):
        joints.append(Joint(f"XCAR{i}", "BED", "TABLE", "z", 230, {"x": (x0, x1), "y": (y0, y1)}, carriage_k("z", "x")))
    joints.append(Joint("XNUT", "BED", "TABLE", "z", 230, {"x": (460, 540), "y": (280, 320)}, {1: P.K_NUT}))
    for i, (x0, x1, y0, y1) in enumerate(ycar, 1):
        joints.append(Joint(f"YCAR{i}", "BED", "COLUMN", "z", 230, {"x": (x0, x1), "y": (y0 + dy, y1 + dy)}, carriage_k("z", "y")))
    joints.append(Joint("YNUT", "BED", "COLUMN", "z", 230, {"x": (480, 520), "y": (730 + dy, 770 + dy)}, {2: P.K_NUT}))
    for i, (x0, x1, z0, z1) in enumerate(zcar, 1):
        joints.append(Joint(f"ZCAR{i}", "COLUMN", "HEAD", "y", 630 + dy, {"x": (x0, x1), "z": (z0 + dz, z1 + dz)}, carriage_k("y", "z")))
    joints.append(Joint("ZNUT", "COLUMN", "HEAD", "y", 630 + dy, {"x": (480, 520), "z": (720 + dz, 780 + dz)}, {3: P.K_NUT}))

    feet = [(50, 150, 50, 150), (850, 950, 50, 150), (50, 150, 950, 1050), (850, 950, 950, 1050)]
    contact = dict(xyz=(500.0, 400.0 + dy, 400.0 + dz), rect={"x": (470, 530), "y": (370 + dy, 430 + dy)})
    return parts, joints, feet, contact
