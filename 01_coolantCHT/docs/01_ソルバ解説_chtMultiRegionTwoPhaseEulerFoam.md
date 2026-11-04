# 01 ソルバ解説：chtMultiRegionTwoPhaseEulerFoam

作成日：2026-11-03　OpenFOAM v2512（ESI）
対象：このフォルダの run001・run004（`system/controlDict` の `application chtMultiRegionTwoPhaseEulerFoam;`）

---

## 1. ひとことで言うと

**「気液二相流（オイラー–オイラー法）」と「固体の熱伝導」を、領域（region）ごとに分けて解き、
境界で熱をやり取りさせる非定常ソルバ**。

ソースの説明文（`applications/solvers/heatTransfer/chtMultiRegionFoam/chtMultiRegionTwoPhaseEulerFoam/`）：

> Transient solver for buoyant, turbulent fluid flow and solid heat conduction with conjugate heat
> transfer between solid and fluid regions. It solves a two-phase Euler approach on the fluid region.

| 名前の部分 | 意味 |
|---|---|
| cht | Conjugate Heat Transfer（共役熱伝達：流体と固体の熱連成） |
| MultiRegion | 流体領域と固体領域を別々のメッシュ（region）として持つ |
| TwoPhaseEuler | 流体領域を「気相＋液相」の 2 つの連続体として解く（オイラー–オイラー法） |
| Foam | OpenFOAM のソルバ |

---

## 2. 何を解いているか

### 2.1 流体領域（例：`water` 領域）— オイラー–オイラー二相流

気相（gas）と液相（liquid）を、それぞれ**体積分率 α と速度・温度をもつ連続体**として重ね合わせる。

- 体積分率：α<sub>gas</sub> + α<sub>liquid</sub> = 1
- 各相ごとに：質量保存、運動量保存、エネルギー（エンタルピー）保存
- 圧力 p は両相で共通（`p_rgh` = p − ρgh を解く）
- 相の間のやり取り（`constant/<流体領域>/phaseProperties`）
  - 抗力 `drag`（IshiiZuber）、仮想質量 `virtualMass`、揚力 `lift`、壁潤滑 `wallLubrication`、乱流分散
  - 相間熱伝達 `heatTransfer`（RanzMarshall など）
  - どちらが「連続相」でどちらが「分散相」かは α で自動的に切り替わる（`blending`）
    - `(gas in liquid)`：液の中の気泡
    - `(liquid in gas)`：空気の中の液滴（クーラントの噴流・飛沫はこちら）
- 乱流：相ごとに k-ε（`turbulenceProperties.gas` / `.liquid`）

**VOF（interFoam）との違い**：VOF は 1 つの速度場で界面を鋭く追う。
オイラー–オイラー法は界面を追わず、セルの中に気と液が「混ざって」存在してよい。
噴霧・飛沫・気泡流に向き、薄い液膜や界面の形状そのものは苦手。

### 2.2 固体領域（例：`solid` 領域）— 熱伝導

$$
\rho c_p \frac{\partial T}{\partial t} = \nabla\cdot(\kappa\nabla T) + Q
$$

物性は `constant/solid/thermophysicalProperties`（本ケースは鋼）。

### 2.3 流体と固体の境界 — 熱の受け渡し

`splitMeshRegions` が作る境界 `water_to_solid` / `solid_to_water` に、
`compressible::turbulentTemperatureTwoPhaseRadCoupledMixed` という境界条件を置く。

- 両側の温度と熱流束が一致するように、毎反復で互いの値を渡し合う
- 流体側では、気相・液相それぞれの熱伝達（壁関数 `alphat`）を α で重み付けして合計する
  （`kappaMethod phaseSystem;`）
- 液が濡れている所は熱伝達が大きく、空気だけの所は小さくなる
  → **クーラントがかかる場所の熱伝達率が自然に大きくなる**

---

### 2.4 基礎式（まとめ）

記号：相 k（g = 気相＝空気、l = 液相＝クーラント）、α<sub>k</sub> 体積分率、ρ<sub>k</sub> 密度、
**u**<sub>k</sub> 速度、h<sub>k</sub> 比エンタルピー、p 圧力（両相で共通）、**g** 重力加速度。

**(1) 体積分率の拘束**

$$
\alpha_g + \alpha_l = 1
$$

**(2) 各相の質量保存**（相変化なしなので右辺 0）

$$
\frac{\partial (\alpha_k \rho_k)}{\partial t} + \nabla\cdot(\alpha_k \rho_k \mathbf{u}_k) = 0
$$

**(3) 各相の運動量保存**

$$
\frac{\partial (\alpha_k \rho_k \mathbf{u}_k)}{\partial t} + \nabla\cdot(\alpha_k \rho_k \mathbf{u}_k \mathbf{u}_k)
= -\alpha_k \nabla p + \nabla\cdot\left(\alpha_k \boldsymbol{\tau}_k^{\mathrm{eff}}\right) + \alpha_k \rho_k \mathbf{g} + \mathbf{M}_k ,
\qquad \mathbf{M}_g = -\mathbf{M}_l
$$

- **τ**<sup>eff</sup> は分子粘性＋乱流粘性（k-ε）の応力
- **M**<sub>k</sub> は相間で受け渡す力（`phaseProperties` で選ぶ）

$$
\mathbf{M}_l = \underbrace{K(\mathbf{u}_g - \mathbf{u}_l)}_{\text{抗力}}
+ \underbrace{C_{vm}\alpha_d\rho_c\left(\frac{D\mathbf{u}_c}{Dt}-\frac{D\mathbf{u}_d}{Dt}\right)}_{\text{仮想質量}}
+ \text{揚力} + \text{壁潤滑} + \text{乱流分散},
\qquad
K = \frac{3}{4}\,C_D\,\frac{\alpha_d \rho_c\,|\mathbf{u}_g-\mathbf{u}_l|}{d}
$$

（c = 連続相、d = 分散相。空気の中の液滴なら c = 空気・d = 液。C<sub>D</sub> は Ishii–Zuber の式。
どちらが連続相かは α で自動的に切り替わる＝ `blending`）

**(4) 各相のエネルギー保存**

$$
\frac{\partial (\alpha_k \rho_k h_k)}{\partial t} + \nabla\cdot(\alpha_k \rho_k \mathbf{u}_k h_k)
= \alpha_k \frac{\partial p}{\partial t} + \nabla\cdot\left(\alpha_k \kappa_k^{\mathrm{eff}} \nabla T_k\right) + Q_k ,
\qquad Q_l = h_i\, a_i\, (T_g - T_l) = -Q_g
$$

- 界面面積濃度 a<sub>i</sub> = 6α<sub>d</sub>/d、界面の熱伝達率 h<sub>i</sub> は Ranz–Marshall の式
  Nu = 2 + 0.6 Re<sup>1/2</sup> Pr<sup>1/3</sup>

**(5) 固体（テーブル・ワーク）の熱伝導**

$$
\rho_s c_{p,s} \frac{\partial T_s}{\partial t} = \nabla\cdot(\kappa_s \nabla T_s)
$$

（流れも変形も解かない。温度だけ）

**(6) 流体と固体の境界の条件**（`turbulentTemperatureTwoPhaseRadCoupledMixed`）

$$
T_{\text{流体側}} = T_s ,\qquad
\sum_{k=g,l} \alpha_k\,\kappa_k^{\mathrm{eff}} \frac{\partial T_k}{\partial n} = \kappa_s \frac{\partial T_s}{\partial n}
$$

境界の温度が等しく、熱流束がつり合う。流体側の熱流束は相ごとの壁関数の熱伝達を α で重み付けして足すので、
**液が濡れている所ほど熱がよく伝わる**。

## 3. 1 タイムステップの流れ

`chtMultiRegionTwoPhaseEulerFoam.C` の本体は次の順で回っている。

1. 流体の Courant 数・固体の拡散数から Δt を決める（`maxCo`, `maxDi`, `maxDeltaT`）
2. **流体領域を解く**（`fluid/solveFluid.H`）
   - α（体積分率）→ 運動量 → 圧力補正 → エネルギー → 乱流
3. **固体領域を解く**（`solid/solveSolid.H`）：熱伝導方程式
4. 必要なら 2〜3 を外側反復（`nOuterCorrectors`）で繰り返し、境界の温度・熱流束を合わせる
5. 書き出し

---

## 4. ケースの作り方（run001・run004 共通）

```
blockMesh                         # 全体を 1 つのメッシュで作る
topoSet                           # 固体にするセル(cellZone)、入口・排出口の面(faceSet)を選ぶ
createPatch -overwrite            # （run004）面集合から inlet / drain パッチを作る
splitMeshRegions -cellZonesOnly   # cellZone ごとに water / solid の 2 領域に分ける
setFields -region water           # 初期の液の分布
chtMultiRegionTwoPhaseEulerFoam   # 本計算
```

主な設定ファイル

| ファイル | 内容 |
|---|---|
| `constant/regionProperties` | `fluid (water); solid (solid);` どの領域が流体/固体か |
| `constant/water/phaseProperties` | 相の組み合わせ、相間の力・熱伝達、表面張力、飽和温度 |
| `constant/water/thermophysicalProperties.gas/.liquid` | 気相・液相の物性 |
| `constant/solid/thermophysicalProperties` | 固体の物性 |
| `system/controlDict` | 終了時刻、`maxCo`、`maxDi`（固体の拡散数）、`maxDeltaT` |
| `system/<領域>/fvSchemes, fvSolution` | 領域ごとの離散化・解法 |
| `0.orig/<領域>/*` | 領域ごとの初期値・境界条件 |

---

## 5. 各ケースでの使い方

| ケース | 内容 | 気相 | 液相 | 固体 |
|---|---|---|---|---|
| run001 | 公式例題 solidQuenching2D（2 次元・焼入れ、沸騰あり）。動作確認用 | 水蒸気 | 水 | 加熱体 1123 K |
| run004 | 加工エリア：空気の中で、温かいクーラントをテーブル中央のワークにかける | **空気** | クーラント（水物性） | 鋼テーブル＋ワーク |

（途中で作った run002〈run001 の 3 次元版〉と run003〈テーブルが床にある旧版〉は削除済み）

run004 では、気相の物性を空気（分子量 28.96、c<sub>p</sub> 1007、μ 1.81e-5、Pr 0.71）に変え、
クーラントを 35 ℃（ポンプで温まった想定）、空気・テーブル・囲いを 20 ℃ にしている。
液は沸点（`Tsat_water_1_2bar.csv`）よりずっと低いので、例題に入っていた沸騰のモデルは働かない。

---

## 6. 注意点・限界

- 液の噴流や液膜を「界面として」は表さない（α が 0〜1 の間でぼやけた分布になる）。
  液膜の厚さや跳ね返りを詳しく見たいときは VOF 系（interFoam）が必要だが、
  VOF と固体の熱連成は標準ソルバにはない。
- 囲いの壁は固体として解かず、温度一定（20 ℃）の壁にしている。
  壁の中の熱伝導も見たい場合は、壁を固体領域として追加する。
- 並列計算（mpirun）がこの環境では起動したまま止まるため、現状は 1 並列で計算している。
