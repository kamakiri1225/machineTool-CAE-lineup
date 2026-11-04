# 01 圧縮性ソルバの選び方と transonic オプション

作成日：2026-11-03　対象：run002・run003（OpenFOAM v2512）。※ 最初の試し計算 run001（平板のみ・JANAF 物性）は削除済み

レーザー切断のアシストガス（O₂ と空気の混合ガス）は、ノズル出口で音速に達し、
条件によっては超音速になる。このメモは、そのような流れを OpenFOAM で解くときの
ソルバと設定の選び方をまとめたもの。

---

## 1. 定常か非定常か：LTS による疑似定常

run002・run003 は `reactingFoam`（非定常ソルバ）を使っているが、`system/fvSchemes` で

```
ddtSchemes
{
    default         localEuler;   // 局所時間刻み（LTS）
}
```

としているので、**実質的には定常解析**である。

- セルごとに、そのセルの流速・格子幅から決まる「自分に合った大きな疑似時間刻み」で進める。
- 時間の途中経過には物理的な意味がない。繰り返すうちに時間微分項が 0 に近づき、
  収束した状態が定常解になる。
- `controlDict` の `endTime 6000; deltaT 1;` は「6000 回繰り返す」という意味で、秒ではない。
- 完全な定常ソルバ `rhoSimpleFoam`（SIMPLE 法）は 1 種類のガスしか扱えないため、
  O₂/N₂ の混合ガスを解くためにこの方法を選んだ。

LTS の調整は `system/fvSolution` の `PIMPLE` にある。

| 項目 | 値 | 意味 |
|---|---|---|
| `maxCo` | 0.2 | セルごとのクーラン数の上限 |
| `maxDeltaT` | 1e-6 | 疑似時間刻みの上限。最初は流れが止まっていて刻みが上限で決まるため、小さくしないと発散した |
| `alphaTemp` / `alphaY` | 0.05 | 温度・濃度の変化で刻みを制限する |
| `rDeltaTSmoothingCoeff` | 0.02 | 隣のセルと刻みが急に変わらないようにならす |

また、止まっている気体に急に 2 気圧をかけると最初の数ステップで発散したため、
入口の全圧を `uniformTotalPressure` で最初の 1000 回の間に 1 気圧から 2 気圧まで上げている。

---

## 2. reactingFoam と rhoReactingFoam の違い

どちらも「多成分ガス・圧縮性・非定常（PIMPLE）」のソルバで、違いは**密度の扱い**だけ。

| | reactingFoam | rhoReactingFoam |
|---|---|---|
| 熱物性クラス | `psiReactionThermo` | `rhoReactionThermo` |
| 密度 | 圧縮率 ψ を使い ρ = ψp | 状態方程式から ρ を直接計算 |
| 向いている流れ | 高速・高マッハ数（遷音速〜超音速、衝撃波） | 低マッハ数（燃焼炉など）。液体の状態方程式も使える |
| `transonic` オプション | 使える | 基本的に使わない |

アシストガスは音速を超えるので、**psi 型の `reactingFoam` に `transonic yes` を組み合わせる**のが合う。
（化学反応は `constant/chemistryProperties` の `chemistry off;` で切っているので、O₂ と N₂ が混ざるだけ。）

---

## 3. `transonic yes` を入れると何が変わるか

`system/fvSolution` の `PIMPLE` に `transonic yes;` を書くと、
`reactingFoam` の圧力方程式（`applications/solvers/combustion/reactingFoam/pEqn.H`）が別の形に切り替わる。

### transonic なし（低マッハ数向け）

面の質量流束を「密度 × 速度」で作り、密度は**今わかっている値で固定**して圧力を解く。

$$
\frac{\partial (\psi p)}{\partial t} + \nabla\cdot\big(\underbrace{\rho_f\,\boldsymbol{H}/A}_{\text{既知}}\big) - \nabla\cdot\left(\frac{\rho}{A}\nabla p\right) = 0
$$

圧力は**ラプラシアン（拡散のような形）だけ**で決まる。
圧力の変化が一瞬で全体に伝わる、という低速流れの性質に合った形。

### transonic あり（遷音速・超音速向け）

密度を ρ = ψp と**圧力で書き直し**、質量流束を「ψ × 速度 × 圧力」とする。

$$
\frac{\partial (\psi p)}{\partial t} + \nabla\cdot\big(\underbrace{\psi_f\,\boldsymbol{H}/A}_{\phi_d}\;p\big) - \nabla\cdot\left(\frac{\rho}{A}\nabla p\right) = 0
$$

コードでは次の 2 行が増える。

```cpp
surfaceScalarField phid("phid", fvc::interpolate(psi)*(fvc::flux(HbyA) + ...));
fvm::ddt(psi, p) + fvm::div(phid, p) - fvm::laplacian(rhorAUf, p)
```

`fvm::div(phid, p)` という**圧力の対流項**が加わるのがポイント。

| | transonic なし | transonic あり |
|---|---|---|
| 圧力方程式の形 | 拡散型（楕円型） | 対流＋拡散型（双曲型の性質をもつ） |
| 圧力・密度の情報の伝わり方 | 一瞬で全体へ | 流れに乗って運ばれる（超音速では下流にしか伝わらない） |
| 衝撃波・チョーク（音速での流量頭打ち） | うまく表せない。収束しにくい | 表せる |
| 必要な設定 | なし | `fvSchemes` に `div(phid,p)`（本ケースは `Gauss upwind`） |
| 注意 | 超音速の領域があると不安定・不正確 | 低速だけの流れでは不要。upwind の影響で衝撃波は数セル分なまる |

**超音速では、下流の圧力変化は上流に伝わらない**（音より速く流れているため）。
transonic なしの形ではこれを表せないため、チョークや衝撃波を含む流れでは `transonic yes` が必要になる。

---

## 4. 最初の試し計算（run001、削除済み）で起きたこと・改善

- 入口全圧 2.0 bar（絶対圧）と大気圧の比は 1.97 で、音速になる限界の比 1.89 を超えるため、
  ノズル出口で音速、出た直後に**最大マッハ数 1.13**（1000 回目）になった。
- ただし、物性データ（JANAF）の適用範囲は 250 K 以上で、超音速の膨張で温度がそれより下がるため、
  温度が 250 K に張り付いた（理論的にはマッハ 1.13 で約 239 K）。
- 次の `run002` では、
  1. 物性を比熱一定（`hConst`）に変え、低温でも使えるようにする
  2. 入口圧力を 0.5 MPa（ゲージ圧）程度に上げ、マッハディスクのある超音速噴流にする
  3. ノズル出口〜鋼板の間のメッシュを細かくする
- 衝撃波をもっと鋭く捉えたい場合は、密度ベースの `rhoCentralFoam`（非定常・1 成分のみ）と比較する。
