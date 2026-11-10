# 05 手順書：OpenFOAM 13 で「混合ガス＋固体」の固体連成（CHT）を定常で解く

作成日：2026-11-11　対象ケース：`run013_nozzleTemperature_OF13_steady`（結果は 02 の 5.9）

この手順書は、レーザー切断のアシストガス（窒素）と、加工する板・ノズル（固体）を同時に解いた計算の**設定ファイルを 1 つずつ読めるようにする**ためのものです。
ファイルの中身は、実際に計算したケースからそのまま引用しています。

---

## 0. 全体像

### 0.1 何を解いているか

| 領域（region） | 種類 | 中身 | OpenFOAM 13 の部品（モジュール） | 時間の扱い |
|---|---|---|---|---|
| `air` | 流体 | 窒素の噴流と周りの空気（N₂/O₂ の混合ガス）。超音速 | `multicomponentFluid` | LTS（局所時間刻み）で定常へ |
| `plate` | 固体 | 加工する板（鋼）。切断面は 1,800 K に固定 | `solid` | 定常 |
| `nozzle` | 固体 | ノズル（仮に銅）。取り付け面は 300 K | `solid` | 定常 |

- 計算は 5° のくさび形（軸対称）。中心軸は y 軸、単位は m。
- 3 つの領域は、境目の面で温度と熱流束をやりとりする（`coupledTemperature`）。

### 0.2 なぜ OpenFOAM 13 なのか

ESI 版（v2512・v2606）の CHT ソルバでは、「定常」「混合ガス」「超音速（`transonic`）」の 3 つを同時に満たせなかった（02 の 5.2）。
Foundation 版の OpenFOAM 13 では、領域ごとに部品を選ぶ `foamMultiRun` になっていて、次の組み合わせができる。
- 流体：`multicomponentFluid`（混合ガス・`transonic`・LTS に対応）
- 固体：`solid`（定常の熱伝導）

### 0.3 フォルダの中身（ケースの構成）

```
run013_nozzleTemperature_OF13_steady/
├── system/
│   ├── controlDict            ← 使う部品（regionSolvers）、反復回数、途中経過の出力
│   ├── fvSolution             ← 全体の PIMPLE（外側の反復回数）
│   ├── blockMeshDict          ← メッシュ（全領域をまとめて作る）
│   ├── decomposeParDict       ← 並列の分割数
│   ├── air/    fvSchemes, fvSolution, fvConstraints
│   ├── plate/  fvSchemes, fvSolution, fvConstraints, topoSetDict
│   └── nozzle/ fvSchemes, fvSolution
├── constant/
│   ├── air/    physicalProperties, momentumTransport, thermophysicalTransport
│   ├── plate/  physicalProperties
│   └── nozzle/ physicalProperties
└── 0/
    ├── air/    p, U, T, N2, O2, Ydefault, k, omega, nut, alphat
    ├── plate/  T
    └── nozzle/ T
```

**ポイント**：領域ごとの設定は、`system/<領域名>/`・`constant/<領域名>/`・`0/<領域名>/` に分けて置く。

---

## 1. 環境の準備

```bash
source /opt/openfoam13/etc/bashrc ParaView_TYPE=none
```

- `ParaView_TYPE=none` を付ける理由：この WSL では、環境の読み込み中に ParaView の確認（`pvserver --version`）が止まってしまう。
- ESI 版（`/usr/lib/openfoam/openfoam2512`）と同じシェルで混ぜないこと。

---

## 2. メッシュを作る

### 2.1 手順

```bash
cd 02_laserAssistGas
python3 tools/make_nozzle_temperature_case.py run011_nozzleTemperature_OF13_transient   # 非定常の元ケースを作る
cd run011_nozzleTemperature_OF13_transient
blockMesh                                   # 全領域を 1 つのメッシュとして作る
splitMeshRegions -cellZones -overwrite      # セルの名前（cellZone）ごとに領域へ分ける
topoSet -region plate                       # 板の中の「切断面」のセルに名前を付ける
cp -r 0.orig 0                              # 初期値と境界条件
```

### 2.2 `system/blockMeshDict`：ブロックに領域名を付ける

```
blocks
(
    hex (0 1 2 3 0 4 5 3)     air    (10 33 1) simpleGrading (1 0.25 1)
    hex (1 6 7 2 4 8 9 5)     air    (47 33 1) simpleGrading (1 0.25 1)
    ...
    hex (...)                 plate  (...)     simpleGrading (...)
    hex (...)                 nozzle (...)     simpleGrading (...)
);
```

- `hex (頂点8個)` のあとの `air`・`plate`・`nozzle` が **cellZone の名前**。
- `splitMeshRegions -cellZones` は、この名前ごとにメッシュを分けて、`constant/air/polyMesh` などを作る。
- 分けたときにできる境目の面には、`air_to_plate`・`plate_to_air` のように **「自分の領域_to_相手の領域」** という名前が自動で付く。

### 2.3 できる境界（patch）の一覧

| 領域 | 境界の名前 | 意味 |
|---|---|---|
| `air` | `inlet` | ノズルの入口（窒素が入る） |
| | `ambient` | 周りの空気との境（出入り自由） |
| | `air_to_plate` | 板との境目（熱をやりとり） |
| | `air_to_nozzle` | ノズルとの境目（熱をやりとり） |
| | `front`・`back` | くさびの両面（`wedge`） |
| `plate` | `plate_to_air` | ガスとの境目 |
| | `plateEdge` | 板の外周（半径 8 mm）。300 K 固定 |
| `nozzle` | `nozzle_to_air` | ガスとの境目 |
| | `nozzleMount` | ノズルの取り付け面。300 K 固定 |

### 2.4 `system/plate/topoSetDict`：切断面のセルに名前を付ける

```
actions
(
    {
        name    laserFrontSet;
        type    cellSet;
        action  new;
        source  boxToCell;
        box     (0 -0.002001 -1) (0.00019 1e-6 1);    // 切断溝に接する板のセル（半径 0.15〜0.19 mm）
    }
    {
        name    laserFront;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     laserFrontSet;
    }
);
```

- 箱（`box`）の中にある板のセルを選んで、`laserFront` という名前の cellZone にする。
- この名前は、あとで「温度を 1,800 K に固定する」設定（4.3）で使う。

---

## 3. 非定常のケースから定常のケースを作る

```bash
cd 02_laserAssistGas
python3 tools/make_nozzle_steady_case.py run011_nozzleTemperature_OF13_transient run013_nozzleTemperature_OF13_steady
```

このスクリプトがする書き換えは、次の 4 つだけ。

| ファイル | 非定常（run011） | 定常（run013） |
|---|---|---|
| `system/air/fvSchemes` | `ddtSchemes { default Euler; }` | `ddtSchemes { default localEuler; }` |
| `system/plate/fvSchemes`・`system/nozzle/fvSchemes` | `default Euler;` | `default steadyState;` |
| `system/controlDict` | `deltaT 1e-9;`、`endTime 0.05;`（秒） | `deltaT 1;`、`endTime 12000;`（反復回数） |
| `0/air/p` の入口の圧力の上げ方 | 0.2 ms かけて上げる | 1,500 回かけて上げる |

LTS では「時刻」が「反復回数」の意味になる。そのため、`endTime` と圧力を上げる表の横軸も、反復回数で書く。

---

## 4. 設定ファイルの読み方

### 4.1 `system/controlDict`：領域ごとに部品を選ぶ

```
application     foamMultiRun;

regionSolvers
{
    air             multicomponentFluid;    // 混合ガスの流体
    plate           solid;                  // 固体の熱伝導
    nozzle          solid;
}

startFrom       latestTime;
stopAt          endTime;
endTime         12000;          // LTS の反復回数
deltaT          1;              // 1 回 = 1 反復
writeControl    timeStep;
writeInterval   12000;          // 定常なので最後の 1 回だけ保存
adjustTimeStep  no;
```

- **`regionSolvers`** が、OpenFOAM 13 の CHT の要。左が領域名、右が使う部品。
- 途中経過（`functions`）は、20 回ごとに 4 つを記録している。
  - `nozzleProbes`：ノズルの中の 4 点の温度
  - `nozzleAverage`：ノズルの平均温度
  - `inletFlow`：入口の流量
  - `fieldBounds`：ガスの最低の圧力・温度

```
nozzleProbes
{
    type            probes;
    libs            ("libsampling.so");
    region          nozzle;                 // どの領域の値か
    fields          (T);
    interpolationScheme cellPoint;
    writeControl    timeStep;
    writeInterval   20;
    probeLocations
    (
        (0.00115 0.0009 0)      // P1 先端
        (0.00165 0.0023 0)      // P2 内壁側
        (0.00285 0.0023 0)      // P3 外壁側
        (0.0022  0.0035 0)      // P4 取り付け側
    );
}
```

- `region` を書かないと、最初の領域（`air`）の値になる。固体の値を記録するときは必ず書く。
- **注意**：壁の熱流束の場（`wallHeatFlux`）を 20 回ごとに書き出すと、そのたびに時刻のフォルダができる（最初の計算で約 600 個できた）。壁の熱流束は `writeControl writeTime;`（保存のときだけ）にする。

### 4.2 流体（`air`）

**物性 `constant/air/physicalProperties`**

```
thermoType
{
    type            heRhoThermo;            // 密度を状態方程式から求める
    mixture         multicomponentMixture;  // 複数の成分（混合ガス）
    transport       const;                  // 粘性一定
    thermo          hConst;                 // 比熱一定
    energy          sensibleEnthalpy;
    equationOfState perfectGas;             // 理想気体
    specie          specie;
}

defaultSpecie   N2;     // 残りの成分（N2 = 1 − ほかの成分）

species
(
    O2
    N2
);

N2
{
    specie
    {
        molWeight       28.0134;
    }
    thermodynamics
    {
        Cp              1040;       // [J/(kg·K)]
        hf              0;
    }
    transport
    {
        mu              1.79e-05;   // [Pa·s]
        Pr              0.71;
    }
}
// O2 も同じ形（molWeight 31.9988、Cp 918、mu 2.06e-05、Pr 0.72）
```

- ESI 版の `thermophysicalProperties`＋`thermo.*` が、OpenFOAM 13 では 1 つのファイル `physicalProperties` にまとまっている。
- 比熱一定にした理由：超音速で膨張すると 150 K 近くまで冷え、JANAF の多項式の範囲（200 K〜）を外れるため。

**乱流 `constant/air/momentumTransport`・`thermophysicalTransport`**

```
simulationType  RAS;

RAS
{
    model           kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
```

```
RAS
{
    model           unityLewisEddyDiffusivity;  // 熱と成分の乱流拡散を同じ速さにする
    Prt             0.85;
}
```

**解き方 `system/air/fvSchemes`・`fvSolution`**

```
ddtSchemes
{
    default         localEuler;     // LTS：セルごとに許される一番大きな時間刻みで進める
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss upwind;
    div(phi,Yi_h)   Gauss upwind;   // 成分とエンタルピー
    div(phi,K)      Gauss upwind;   // 運動エネルギー
    div(phid,p)     Gauss upwind;   // 圧力の対流項（transonic のとき使う）
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
```

```
PIMPLE
{
    momentumPredictor   yes;
    transonic           yes;        // 超音速・衝撃波を扱う圧力の式（最重要）
    nCorrectors         2;
    nNonOrthogonalCorrectors 1;

    maxCo               0.2;        // セルごとの Courant 数の上限（LTS の時間刻みを決める）
    maxDeltaT           1e-6;
    rDeltaTSmoothingCoeff 0.02;     // 隣のセルとの時間刻みの差をならす
    rDeltaTDampingCoeff 1;
    alphaTemp           1;          // 温度による時間刻みの制限は使わない
}
```

- `transonic yes` がないと、ノズルのチョーク（音速で詰まる）と超音速の膨張で計算が止まる。
- `div(phid,p)` は、`transonic yes` のときに使われる圧力の対流項。

**温度の上下限 `system/air/fvConstraints`**

```
limitT
{
    type            limitTemperature;
    select          all;
    min             50;
    max             2500;
}
```

- 起動直後に温度が一瞬暴れても、計算が止まらないようにする安全装置。収束したあとは効いていないことをログで確かめる。

### 4.3 固体（`plate`・`nozzle`）

**物性 `constant/plate/physicalProperties`**

```
thermoType      constSolidThermo;   // 物性一定の固体

rho
{
    type            uniform;
    value           7850;           // 密度 [kg/m³]
}
Cv
{
    type            uniform;
    value           460;            // 比熱 [J/(kg·K)]
}
kappa
{
    type            uniform;
    value           45;             // 熱伝導率 [W/(m·K)]
}
```

- ノズル（`constant/nozzle/physicalProperties`）は同じ形で、銅の値（8,960、385、400）。
- 定常では、密度と比熱は結果に効かない（熱伝導率だけが効く）。

**解き方 `system/plate/fvSchemes`**

```
ddtSchemes
{
    default         steadyState;    // 定常の熱伝導
}

laplacianSchemes
{
    default         Gauss linear corrected;
}
```

- `solid` は LTS に対応していない。ガスを LTS で回すときは、固体を `steadyState` にする。
- 固体の温度の時刻歴がほしいときは、ガスも固体も `Euler`（本当の時間）で解く必要がある（run011）。

**切断面の温度を固定 `system/plate/fvConstraints`**

```
laserFront
{
    type            fixedTemperature;
    cellZone        laserFront;         // 2.4 で作ったセルの名前
    mode            uniform;
    temperature     constant 1800.0;    // 鋼の融点に近い温度 [K]
}
```

- レーザーで溶けている切断面の代わり。レーザーの熱の入り方そのものは解いていない。

### 4.4 境界条件（`0/<領域>/`）

**ガス `0/air/`**

| 境界 | `p` | `U` | `T` | `N2`・`O2` |
|---|---|---|---|---|
| `inlet` | `uniformTotalPressure`（全圧を上げる） | `pressureInletOutletVelocity` | `totalTemperature`（300 K） | N₂ = 1、O₂ = 0 |
| `ambient` | `totalPressure`（101,325 Pa） | `pressureInletOutletVelocity` | `inletOutlet`（300 K） | `inletOutlet`（空気：O₂ 0.233） |
| `air_to_plate`・`air_to_nozzle` | `zeroGradient` | `noSlip` | **`coupledTemperature`** | `zeroGradient` |

入口の全圧（`0/air/p`）

```
inlet
{
    type            uniformTotalPressure;
    p0
    {
        type    table;
        values  ((0 101325) (1500 901325.0) (1e9 901325.0));   // 1,500 回かけて 0.8 MPa（ゲージ）へ
    }
    value           uniform 101325.0;
}
```

- いきなり 0.9 MPa にすると、最初の数回で計算が壊れやすい。少しずつ上げる。

ガス側の熱の受け渡し（`0/air/T`）

```
air_to_plate
{
    type            coupledTemperature;     // 相手（板）と温度・熱流束をつり合わせる
    value           uniform 300.0;
}
```

壁の近くの乱流の熱拡散（`0/air/alphat`）

```
air_to_plate
{
    type            compressible::alphatJayatillekeWallFunction;
    Prt             0.85;
    value           uniform 0;
}
```

- 壁の近くの温度分布の形から、乱流による熱の伝わりやすさを決める（熱伝達率を入力しなくてよい理由）。

**固体 `0/plate/T`・`0/nozzle/T`**

```
plate_to_air
{
    type            coupledTemperature;     // ガスとの境目
    value           uniform 300.0;
}
plateEdge
{
    type            fixedValue;             // 板の外周は 300 K
    value           uniform 300.0;
}
```

```
nozzleMount
{
    type            fixedValue;             // ノズルの取り付け面は 300 K（仮定）
    value           uniform 300.0;
}
```

ESI 版との名前の違い

| 内容 | ESI v2512 | OpenFOAM 13 |
|---|---|---|
| 熱の受け渡しの境界条件 | `compressible::turbulentTemperatureCoupledBaffleMixed`（`Tnbr`・`kappaMethod` も書く） | `coupledTemperature`（それだけでよい） |
| 物性のファイル | `thermophysicalProperties`、固体は `heSolidThermo` | `physicalProperties`、固体は `constSolidThermo` |
| 乱流のファイル | `turbulenceProperties`（`RASModel`） | `momentumTransport`（`model`） |
| 温度の固定・上下限 | `fvOptions` | `fvConstraints` |

---

## 5. 計算する

```bash
source /opt/openfoam13/etc/bashrc ParaView_TYPE=none
cd run013_nozzleTemperature_OF13_steady
decomposePar -allRegions -force                       # 全領域を 4 つに分ける
export HWLOC_COMPONENTS=-gl                           # WSL2 で mpirun が止まらないようにする
mpirun -np 4 foamMultiRun -parallel > log.foamMultiRun 2>&1
```

- 4 並列で約 28 分（12,000 回）。
- 起動したら、ログの最初のほうに次の 2 行が出ることを確かめる（LTS と定常の組み合わせが受け付けられた印）。

```
Using LTS
PIMPLE: Operating solver in mixed steady-state/transient mode
```

---

## 6. 収束を確かめる

| 見るもの | ファイル | 目安 |
|---|---|---|
| 入口の流量 | `postProcessing/air/inletFlow/0/surfaceFieldValue.dat` | 最後の 1,000 回で変化が 0.1 % 以下（run013：0.01 % 以下） |
| ノズルの温度 | `postProcessing/nozzle/nozzleProbes/0/T` | 4 点とも変化が止まっている |
| ガスの最低温度 | `postProcessing/air/fieldBounds/0/volFieldValue.dat` | 温度の下限 50 K に張り付いていない |

- OpenFOAM 13 では、途中経過のファイルが `postProcessing/<領域名>/<名前>/` に出る（ESI 版とは場所が違う）。

---

## 7. 結果を見る・図を作る

```bash
# 結果を 1 つにまとめる（最後の 1 回だけ）
reconstructPar -allRegions -latestTime

# 斜めのアングルの図（分割したままでも読める）
python3 ../tools/render_cht_regions3d.py . images/regions.png probe_locations.csv   # 領域分け
python3 ../tools/render_cht_fields3d.py  . images/fields.png                        # 流速・温度・マッハ数・O₂
```

- ParaView で見るときは、ケースの `post.foam` を開き、Mesh Regions で `air`・`plate`・`nozzle` の `internalMesh` を選ぶ。
- 軸対称なので、`Rotational Extrusion`（回転体にする）か `Reflect`（折り返す）を使うと見やすい。

---

## 8. よくあるつまずき

| 症状 | 原因 | 対策 |
|---|---|---|
| 環境の読み込みで止まる | ParaView の確認が止まる | `source /opt/openfoam13/etc/bashrc ParaView_TYPE=none` |
| 固体で「LTS に対応していない」と出る | `solid` は LTS 不可 | 固体を `steadyState`（または全体を `Euler`）にする |
| 時刻のフォルダが大量にできる | 壁の熱流束の場を毎回書き出している | 熱流束は `writeControl writeTime;` |
| 固体の値が記録されない | `functions` に `region` がない | `region nozzle;` などを書く |
| 起動直後に温度が負になる | 入口の圧力をいきなり上げた | 入口の全圧を 1,500 回かけて上げる、`limitTemperature` を入れる |
| `volFieldValue` で「cellZone not specified」 | OpenFOAM 13 では対象の指定が必要 | `cellZone all;` を書く |
