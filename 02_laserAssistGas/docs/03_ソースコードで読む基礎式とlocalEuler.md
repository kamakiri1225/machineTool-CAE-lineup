# 03 ソースコードで読む基礎式と localEuler（局所時間刻み）

作成日：2026-11-04　OpenFOAM v2512（ESI）
ソースの場所：`/usr/lib/openfoam/openfoam2512/applications/solvers/combustion/reactingFoam/`
（run002・run003・run004 で使っている `reactingFoam`）

---

## 1. reactingFoam はどんなソルバか（特徴）

| 項目 | 内容 |
|---|---|
| 対象 | **多成分ガス**（化学種ごとに質量分率 Y<sub>i</sub> を解く）の圧縮性流れ。燃焼にも使えるが、本ケースは反応なし（`chemistry off`）で O₂・N₂ を混ぜるだけ |
| 時間の扱い | 非定常（PIMPLE 法）。`localEuler` にすると局所時間刻みで疑似定常（→ 4 章） |
| 解き方 | **圧力ベース**（圧力方程式を立てて速度・密度を直す）。`transonic yes` で圧縮性の強い流れにも対応 |
| 密度 | ψ 型の熱物性（ρ = ψp、`psiReactionThermo`）。混合ガスの ψ は成分ごとの分子量から決まる |
| 混合物の物性 | 比熱・粘性などを成分ごとに持ち、質量分率で混ぜる（`reactingMixture`） |
| 向いている流れ | 低速〜超音速、成分が混ざる流れ（アシストガスと空気の混合など） |
| 苦手な点 | 強い衝撃波は密度ベースの `rhoCentralFoam` よりなまる。成分ごとに方程式が増えるので重い |

**rhoCentralFoam との違い**：`rhoCentralFoam` は密度ベース（Kurganov–Tadmor の中心風上スキームで
質量・運動量・全エネルギーの流束を直接作る）で衝撃波に強いが、**1 成分のガスしか扱えない**。
混合ガスを正しく解くため、run004 は `reactingFoam` を使っている。

---

## 2. 基礎式とソースの対応

1 回の反復で、次の順にファイルが実行される（`reactingFoam.C` の本体）。

```
rhoEqn.H → UEqn.H → YEqn.H → EEqn.H → pEqn.H（圧力を直し、速度・密度・流束を更新）→ 乱流
```

### 2.1 質量保存（`rhoEqn.H`）

$$
\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho\mathbf{u}) = 0
$$

```cpp
fvScalarMatrix rhoEqn( fvm::ddt(rho) + fvc::div(phi) == fvOptions(rho) );
```

`phi` は面を通る質量流束 ρ**u**·**S**。

### 2.2 運動量保存（`UEqn.H`）

$$
\frac{\partial (\rho\mathbf{u})}{\partial t} + \nabla\cdot(\rho\mathbf{u}\mathbf{u})
- \nabla\cdot\left[\mu_{\mathrm{eff}}\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \tfrac{2}{3}(\nabla\cdot\mathbf{u})\mathbf{I}\right)\right]
= -\nabla p
$$

```cpp
fvm::ddt(rho, U) + fvm::div(phi, U) + turbulence->divDevRhoReff(U)  ==  fvOptions(rho, U)
solve(UEqn == -fvc::grad(p));        // momentumPredictor yes のとき
```

| ソース | 式の項 |
|---|---|
| `fvm::ddt(rho, U)` | ∂(ρ**u**)/∂t（`localEuler` のときはセルごとの Δt で離散化 → 4 章） |
| `fvm::div(phi, U)` | 対流項 ∇·(ρ**uu**) |
| `turbulence->divDevRhoReff(U)` | 粘性＋乱流の応力（μ<sub>eff</sub> = μ + μ<sub>t</sub>） |
| `-fvc::grad(p)` | 圧力勾配 |

### 2.3 成分（化学種）の輸送（`YEqn.H`）

$$
\frac{\partial (\rho Y_i)}{\partial t} + \nabla\cdot(\rho\mathbf{u}Y_i) - \nabla\cdot(\mu_{\mathrm{eff}}\nabla Y_i) = \dot{\omega}_i
\qquad (i = \mathrm{O_2}),\qquad Y_{\mathrm{N_2}} = 1 - \sum_{i \ne \mathrm{N_2}} Y_i
$$

```cpp
fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi) - fvm::laplacian(turbulence->muEff(), Yi)
 == reaction->R(Yi) + fvOptions(rho, Yi)
...
Y[inertIndex] = scalar(1) - Yt;      // 不活性種（inertSpecie N2）は残り
```

- 拡散係数は μ<sub>eff</sub>（シュミット数 1 相当）
- `reaction->R(Yi)` は反応による生成（本ケースは 0）
- **ここが「混合ガス」の本体**：O₂ を解き、N₂ は 1 − Y<sub>O₂</sub>。
  密度・比熱・粘性は各セルの組成から決まり、流れ（運動量・圧力）に反映される

### 2.4 エネルギー保存（`EEqn.H`、本ケースは顕エンタルピー h）

$$
\frac{\partial (\rho h)}{\partial t} + \nabla\cdot(\rho\mathbf{u}h)
+ \frac{\partial (\rho K)}{\partial t} + \nabla\cdot(\rho\mathbf{u}K)
- \frac{\partial p}{\partial t}
- \nabla\cdot(\alpha_{\mathrm{eff}}\nabla h) = \dot{Q},
\qquad K = \tfrac12|\mathbf{u}|^2
$$

```cpp
fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
+ fvc::ddt(rho, K) + fvc::div(phi, K)
+ ( he.name() == "e" ? fvc::div(phiv, p) : -dpdt )
- fvm::laplacian(turbulence->alphaEff(), he)
 == Qdot + fvOptions(rho, he)
```

運動エネルギー K の項が入っているので、超音速の膨張で温度が下がる（全温一定で静温が下がる）ことが表せる。

### 2.5 圧力方程式（`pEqn.H`、`transonic yes`）

詳しくは `01_…transonic.md`。ρ = ψp を使い、

$$
\frac{\partial (\psi p)}{\partial t} + \nabla\cdot(\psi_f\,\boldsymbol{\phi}_d\, p) - \nabla\cdot\left(\frac{\rho}{A}\nabla p\right) = 0
$$

```cpp
fvm::ddt(psi, p) + fvm::div(phid, p) - fvm::laplacian(rhorAUf, p)
```

圧力の対流項 `fvm::div(phid, p)` によって、チョーク・衝撃波を含む流れを扱える。

### 2.6 状態方程式と混合物の物性（`constant/thermophysicalProperties`）

$$
p = \rho R_{\mathrm{mix}} T,\qquad R_{\mathrm{mix}} = R_u\sum_i \frac{Y_i}{W_i},\qquad
c_{p,\mathrm{mix}} = \sum_i Y_i\,c_{p,i}
$$

（W<sub>i</sub>：分子量 O₂ 32.0、N₂ 28.0、c<sub>p</sub>：O₂ 918、N₂ 1040 J/kgK、`hConst` = 比熱一定）

---

## 3. 乱流（realizable k-ε、run004）

$$
\frac{\partial(\rho k)}{\partial t} + \nabla\cdot(\rho\mathbf{u}k) = \nabla\cdot\left[\left(\mu+\frac{\mu_t}{\sigma_k}\right)\nabla k\right] + G_k - \rho\varepsilon,
\qquad \mu_t = \rho C_\mu \frac{k^2}{\varepsilon}
$$

realizable 型は C<sub>μ</sub> を流れの変形の速さに応じて変え、噴流の広がりを標準 k-ε より正しく出しやすい
（機械学会の論文と同じモデル）。ソース：`src/TurbulenceModels/turbulenceModels/RAS/realizableKE/`。

---

## 4. localEuler（局所時間刻み）は式の上で何をしているか

### 4.1 ふつうの時間進行（Euler）

全セル共通の Δt で

$$
\frac{(\rho\phi)^{n+1}_P - (\rho\phi)^{n}_P}{\Delta t}\,V_P + \sum_f F_f(\phi^{n+1}) = S_P
$$

Δt は「いちばん小さいセル・いちばん速い流れ」の Courant 数で決まるので、
ノズルの近く（細かいセル・超音速）に合わせた小さな Δt で、遠くの遅い領域まで進めることになる。

### 4.2 localEuler

時間微分の離散化はそのままで、**Δt をセル P ごとに変える**：

$$
\frac{(\rho\phi)^{n+1}_P - (\rho\phi)^{n}_P}{\Delta t_P}\,V_P + \sum_f F_f(\phi^{n+1}) = S_P
$$

各セルの Δt<sub>P</sub>（実際は逆数 1/Δt<sub>P</sub> = `rDeltaT`）は `setRDeltaT.H` で決める：

```cpp
rDeltaT = fvc::surfaceSum(mag(phi)) / ((2*maxCo)*mesh.V()*rho);   // そのセルの Courant 数 = maxCo
rDeltaT.clamp_min(1/maxDeltaT);                                       // Δt の上限
rDeltaT = max(rDeltaT, mag(Qdot)/(alphaTemp*rho*Cp*T));                // 温度の変化が大きすぎないように
（alphaY < 1 なら 成分の変化でも制限）
fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);                          // 隣のセルと急に変わらないように
```

式で書くと

$$
\frac{1}{\Delta t_P} = \max\left(\frac{\sum_f |\phi_f|}{2\,Co_{\max}\,\rho_P V_P},\ \frac{1}{\Delta t_{\max}},\ \frac{|\dot Q_P|}{\alpha_T\,\rho_P c_p T_P},\ \cdots\right)
$$

- 細かいセル・速い流れのセルは小さい Δt<sub>P</sub>、粗いセル・遅い流れのセルは大きい Δt<sub>P</sub>。
- **どのセルも「自分の Courant 数 ≒ maxCo」で進む**ので、全体が少ない反復で落ち着く。

### 4.3 なぜ定常解になるのか

収束すると、どのセルでも (ρφ)<sup>n+1</sup> = (ρφ)<sup>n</sup>、つまり時間微分の項が 0 になり、

$$
\sum_f F_f(\phi) = S_P
$$

と、**定常の方程式そのもの**が残る。Δt<sub>P</sub> の値は最終の解に影響しない（途中の進み方だけを変える）。
一方、途中経過はセルごとに時刻がずれているので、**物理的な時間変化としては意味がない**。
そのため `controlDict` の `endTime`・`deltaT 1` は「反復回数」として使っている。

| | ふつうの時間進行 | localEuler |
|---|---|---|
| Δt | 全セル共通（最も厳しいセルで決まる） | セルごと |
| 途中経過 | 物理的な時間変化 | 意味なし |
| 最終結果 | 定常に達すれば定常解 | 定常解 |
| 速さ | 遅い（遠くの領域がなかなか落ち着かない） | 速い |
| 使えないもの | — | 本当に非定常な現象（渦の放出・振動）は表せない |

`rhoCentralFoam` も同じ仕組み（`rhoCentralFoam/setRDeltaT.H`。こちらは Kurganov の最大波速 `amaxSf` から Δt<sub>P</sub> を決める）。

---

## 5. 設定との対応（run004）

| 設定（`system/fvSolution` の PIMPLE） | 値 | 4 章の式での意味 |
|---|---|---|
| `maxCo` | 0.2 | Co<sub>max</sub> |
| `maxDeltaT` | 1e-6 s | Δt<sub>max</sub> |
| `alphaTemp` / `alphaY` | 0.05 | 温度・成分の 1 回の変化を 5 % 以内に |
| `rDeltaTSmoothingCoeff` | 0.02 | 隣のセルとの Δt の比の制限 |
| `fvSchemes` の `ddtSchemes default localEuler;` | — | 局所時間刻みを使う宣言 |
