[← docs 目次](README.md) ｜ [⌂ トップ](../README.md) ｜ [次: ref005 ばね・固有値・周波数応答の設定比較 →](ref005_CalculiX-FrontISTR_ばね・固有値・周波数応答の設定比較.md)

# ref000 オープンソース マルチフィジックス総覧（工作機械向け）

作成: 2026-11-07
目的: 工作機械の主要テーマ（振動・機構・サーボモータ・ボールねじ・リニアガイド・熱変位）を
オープンソースだけで解析・学習するために、使えるツールを整理する。

> バージョンや開発状況は変わるので、実際に試す時は各公式サイトで確認すること。

---

## 1. そもそもマルチフィジックスとは

複数の物理現象を**連成**させて解くこと。工作機械で典型的なのは次の組み合わせ。

| 連成 | 工作機械での例 |
|---|---|
| 構造 ⇔ 熱 | 主軸・モータ発熱 → 熱膨張 → 加工精度の低下（**熱変位**） |
| 電磁場 ⇔ 構造 ⇔ 熱 | サーボモータ／ビルトインモータの銅損・鉄損 → 発熱、電磁力 → 振動 |
| 機構 ⇔ 制御 ⇔ 構造 | サーボ制御ループ ＋ ボールねじ ＋ 弾性テーブル → 位置決め誤差・共振 |
| 構造 ⇔ 切削力 | **びびり振動（チャタ）**：工具・ワークの振動が切削厚みを変え、それが力を変える |
| 流体 ⇔ 熱 | クーラント・主軸冷却油・空気による冷却 |
| 接触 ⇔ 構造 | リニアガイドの転動体、軸受、ボルト結合部の剛性・減衰 |

「何でも一つのソフトで連成」より、**1D（制御・機構）と 3D（FEM）を使い分けて、必要な所だけ繋ぐ**のが現実的。

---

## 2. 分野別ツール一覧

★＝このプロジェクトで最初に触る候補

### 2.1 構造・振動 FEM（3D）

| ツール | 特徴 | 工作機械での使い道 |
|---|---|---|
| ★**CalculiX** | Abaqus類似の入力形式。軽量、Windows/Linux両対応。静解析・モーダル・周波数応答・熱・接触 | フレーム（ベッド・コラム）の固有振動数、静剛性 |
| ★**code_aster**（Salome-Meca同梱） | EDF（仏電力）製。構造動力学が非常に強い。ローターダイナミクス、非線形接触、減衰、ハーモニック、過渡 | 本格的な振動・接触・熱応力。習得コストは高い |
| ★**FrontISTR** | 日本発（東京大学 奥田研究室が中心、FrontISTR Commons が開発、MITライセンス）。**日本語マニュアルあり**。静解析・非線形・接触・固有値・周波数応答・動解析（陰解法・陽解法）・熱伝導。並列計算に強い | 日本語で学べる構造FEMの本命。ばね・ダッシュポット要素（511）でガイド・ボールねじの剛性を入れられる |
| Elmer FEM | 下記マルチフィジックスの方で詳述 | 構造＋熱＋電磁の連成 |
| FEniCSx / deal.II / MFEM | 方程式を自分でコードで書くライブラリ型 | 研究的な独自モデル（切削力モデル連成など） |
| Kratos Multiphysics | 構造・流体・DEM・連成をPythonから組める | 連成解析の発展編 |

### 2.2 マルチフィジックス統合型

| ツール | 特徴 |
|---|---|
| ★**Elmer FEM** | フィンランドCSC製。構造・熱・電磁場・流体・音響を1ソルバで連成。GUI（ElmerGUI）あり。**モータの電磁場→発熱→熱変位**の流れを1つで追える |
| Kratos Multiphysics | アプリケーション（モジュール）の組合せで連成 |
| MOOSE | 米INL製。大規模連成、C++。原子力系が主だが汎用 |
| preCICE | **連成ライブラリ**。OpenFOAM＋CalculiX など別ソフト同士を結合（FSI、共役熱伝達） |

### 2.3 流体・熱（CFD）

| ツール | 特徴 | 使い道 |
|---|---|---|
| ★OpenFOAM | 定番CFD。共役熱伝達（chtMultiRegionFoam） | 主軸冷却ジャケット、クーラント、機内の空気の流れ |
| SU2 | 圧縮性・随伴最適化に強い | 工作機械ではあまり使わない |

### 2.4 電磁場（サーボモータ・リニアモータ・ビルトインモータ）

| ツール | 特徴 |
|---|---|
| ★**FEMM** | 2D磁場の定番（Windows、無償）。トルク・コギング・インダクタンス計算。Pythonから操作可（pyfemm） |
| **pyleecan** | 電動機設計用のPythonフレームワーク。FEMMやGMSHと連携。PMSM（サーボモータ）の解析フローが揃う |
| SyR-e | 同期リラクタンス／PMモータ設計（MATLAB/Octave＋FEMM） |
| GetDP + Gmsh（ONELAB） | 2D/3D電磁場。回転機の例題あり |
| Elmer FEM | 3D電磁場＋熱連成 |

### 2.5 1D システム・制御（サーボループ、駆動系）

ここが**サーボモータ・ボールねじ・テーブル送り系**の主戦場。

| ツール | 特徴 | 使い道 |
|---|---|---|
| ★**OpenModelica** | Modelica言語のOSS環境。標準ライブラリ（MSL）に Electrical.Machines（モータ）, Mechanics.Rotational / Translational（慣性・ばね・摩擦・ギア）, Blocks（PID）が揃う | **モータ→カップリング→ボールねじ→テーブル**の送り系モデル、位置／速度／電流の3重ループ、摩擦・バックラッシ |
| ★**Python（python-control, SciPy）** | 伝達関数・状態方程式・ボード線図・根軌跡 | 制御ゲイン設計、周波数応答、安定判別 |
| Scilab / Xcos | Simulink的なブロック線図 | 同上（GUIで組みたい場合） |
| GNU Octave | MATLAB互換 | 既存MATLABコードの流用 |
| FMPy / OMSimulator | FMI（モデル交換規格）で別ツールのモデルを連成 | Modelicaモデルを Python から回す |

### 2.6 機構・マルチボディ動力学

| ツール | 特徴 | 使い道 |
|---|---|---|
| ★**Project Chrono** | 剛体・弾性体マルチボディ、接触、FEA、Python（PyChrono） | 機構全体の動き、リニアガイド・転動体の接触 |
| MBDyn | ミラノ工科大。弾性マルチボディ、制御連成 | 弾性体を含む機構 |
| MuJoCo | 高速なマルチボディ＋接触（ロボット向け） | 5軸機構の運動学・動力学の学習 |
| Drake | 制御・最適化と一体のマルチボディ | 軌道生成・制御 |
| FreeCAD Assembly | CADの拘束による機構確認 | 干渉チェック、ストローク確認 |

### 2.7 回転体・主軸（ローターダイナミクス）

| ツール | 特徴 |
|---|---|
| ★**ROSS** | Python製ローターダイナミクス。軸・円板・軸受要素で主軸モデルを組み、危険速度・キャンベル線図・不釣合い応答 |
| OpenTorsion | ねじり振動（駆動系）専用のPythonライブラリ |
| code_aster | 3Dローターダイナミクス（ジャイロ効果込み） |

### 2.8 びびり振動（切削の安定性）

専用OSSは少ない。**理論（Altintas / Tlusty の安定限界線図＝SLD）を Python で自作する**のが一番学びになる。

- 入力：工具先端の周波数応答関数（FRF）→ FEM（CalculiX）やROSSで計算、または打撃試験
- 計算：ゼロ次近似（ZOA）法、半離散化法（Semi-Discretization）
- 出力：安定限界線図（主軸回転数 vs 限界切込み）

### 2.9 実験モーダル解析・信号処理（実機と比較するため）

| ツール | 特徴 |
|---|---|
| sdypy（旧 pyEMA ほか） | Python。実験モーダル解析（FRFからの固有値同定） |
| pyFBS | 周波数応答ベースのサブストラクチャリング（工具＋ホルダ＋主軸の結合 ＝ RCSA的な使い方） |
| SciPy.signal | FFT・フィルタ・FRF推定 |

### 2.10 プリ・ポスト・CAD・形状

| ツール | 用途 |
|---|---|
| ★**FreeCAD** | CAD。FEMワークベンチで CalculiX / Elmer を GUI から実行できる → **入門に最適** |
| ★**Gmsh** | メッシュ生成（Pythonで自動化可） |
| Salome | CAD＋メッシュ＋code_aster連携（Salome-Meca） |
| Netgen | メッシュ（FreeCAD内蔵） |
| ★**ParaView** | 結果の可視化 |

### 2.11 工作機械そのもの（NC・CAM）

| ツール | 用途 |
|---|---|
| LinuxCNC | 実際に機械を動かせるOSSのCNC。補間・加減速・サーボループ・運動学（5軸）の実装が読める |
| FreeCAD CAM（旧Path） | 工具経路生成 |
| grbl | 小型マシン用ファームウェア。加減速計画が読みやすい |

### 2.12 最適化・パラメトリック

| ツール | 用途 |
|---|---|
| OpenMDAO | 多分野設計最適化（NASA製） |
| Dakota | パラメータスタディ・不確かさ評価 |
| Optuna | Python製。ソルバを外から回して最適化 |

### 2.13 ばね要素（ボールねじ・リニアガイドの剛性）の入れ方

工作機械の構造解析では、ボールねじやリニアガイドを**ばね剛性**に置き換えるのが定番。各ソフトの対応は次のとおり。

| ソフト | 地面と結ぶばね | 2点間のばね | $6\times6$ 剛性行列 | 備考 |
|---|---|---|---|---|
| CalculiX | `SPRING1` | `SPRINGA`（軸方向）, `SPRING2`（自由度指定）, `DASHPOTA` | ×（ばねを組み合わせて表す） | FreeCAD の GUI からも「Spring」拘束で面に剛性を与えられる |
| FrontISTR | `!SPRING`（静解析・陰解法動解析で使える） | 要素タイプ **511**（2節点のばね・ダッシュポットのコネクタ要素、`!SPRING_A` / `!SPRING_D`） | × | `!SPRING` は固有値解析の制御データには載っていない。固有値解析では511要素を使う想定だが、**実際に試して確認する** |
| code_aster | `DISCRET`（1節点） | `DISCRET`（2節点）: `K_T_D_L`, `K_TR_D_L` | ○ `K_TR_L`（方向同士が影響し合う剛性も可）。質量・減衰行列、すき間付き非線形（`DIS_CHOC`）も | 最も表現力が高い |
| Elmer | 境界のばね係数（弾性床） | △ | × | 離散的なばねは不得意 |
| OpenModelica | `Fixed`＋`Spring` | `SpringDamper`（並進・回転） | × | 1Dの送り系モデル向け |
| Project Chrono | — | `ChLinkTSDA` / `ChLinkRSDA` | ○ ブッシュ要素 | マルチボディ向け |

剛性の値の作り方:
- **ボールねじの軸方向**: 直列ばねとして合成する

$$
\frac{1}{K} = \frac{1}{K_{\mathrm{shaft}}} + \frac{1}{K_{\mathrm{nut}}} + \frac{1}{K_{\mathrm{bearing}}}, \qquad K_{\mathrm{shaft}} = \frac{AE}{L}
$$

  $K_{\mathrm{shaft}}$：ねじ軸、$K_{\mathrm{nut}}$：ナット、$K_{\mathrm{bearing}}$：サポート軸受。$A$ はねじ軸の断面積、$E$ はヤング率、$L$ は固定側軸受からナットまでの距離（ナットの位置で変わる）
- **リニアガイド**: カタログのラジアル・逆ラジアル・横方向の剛性（予圧の等級で変わる）。モーメント剛性は、キャリッジを複数のばねで支えて表す
- キャリッジやナットの取付面は、1つの節点に剛につなぎ（CalculiX `*RIGID BODY`、FrontISTR `!EQUATION` / `!MPC`、code_aster `LIAISON_SOLIDE`）、その節点とレール側の間にばねを入れる

---

## 3. 学習テーマ × ツール対応表

| 学習テーマ | まず試す | 次に試す | 理論で押さえること |
|---|---|---|---|
| フレームの振動（ベッド・コラム） | FreeCAD＋CalculiX | FrontISTR / code_aster | 固有値解析、モード形状、質量・剛性行列 |
| 主軸の振動・危険速度 | ROSS | code_aster | はり要素、軸受剛性、ジャイロ効果、キャンベル線図 |
| びびり振動 | Python自作SLD | CalculiXのFRFと連携 | 再生びびり、安定限界線図 |
| サーボモータ（電磁） | FEMM / pyleecan | Elmer（3D・熱連成） | dq変換、トルク定数、コギング |
| サーボ制御 | python-control | OpenModelica | 電流・速度・位置ループ、帯域、ゲイン余裕 |
| ボールねじ送り系 | OpenModelica | Chrono＋FEM | リード、慣性換算、軸方向剛性、ねじり剛性、危険速度 |
| リニアガイド | 手計算（Hertz接触）＋CalculiX / FrontISTR（ばね要素） | code_aster（6×6）、Chrono（転動体接触） | 予圧、剛性、寿命計算 |
| 熱変位 | Elmer / CalculiX（熱→構造） | OpenFOAM（冷却） | 熱伝導、熱伝達率、熱膨張 |
| 機構・5軸運動学 | Python（同次変換行列） | MuJoCo / LinuxCNC | 順・逆運動学、幾何誤差 |

---

## 4. 学習ロードマップ（案）

段階的に「手計算 → 1D → 3D → 連成」の順で、**毎回理論値と照合**する。

1. **ana001** FreeCAD＋CalculiX：片持ちはりのモーダル解析（理論解と一致確認＝ツールの使い方を覚える）
2. **ana002** CalculiX：簡易コラム形状の固有振動数・静剛性
3. **ana003** Python：ボールねじ送り系の慣性換算・軸方向剛性・危険速度の手計算ツール
4. **ana004** python-control：位置・速度ループのボード線図、2慣性系の共振
5. **ana005** OpenModelica：モータ＋ボールねじ＋テーブルの送り系シミュレーション（摩擦込み）
6. **ana006** ROSS：主軸（軸受2点支持）の危険速度・キャンベル線図
7. **ana007** Python：1自由度びびり振動の安定限界線図
8. **ana008** FEMM / pyleecan：小型PMSMのトルク・コギング
9. **ana009** Elmer：モータ発熱 → 主軸の熱変位（熱-構造連成）
10. **ana010** 連成：送り系（1D）＋フレームFEMの縮約モデル → 加減速時の振動

---

## 5. 環境メモ（WSL2 / Linux 前提）

| ツール | 入手方法の目安 |
|---|---|
| CalculiX | `apt install calculix-ccx` |
| FrontISTR | 公式サイトのバイナリ、公式Dockerイメージ、またはソースビルド（[frontistr.com](https://www.frontistr.com/)） |
| FreeCAD | 公式AppImage または Windows版（GUIはWindows側が楽） |
| Elmer | PPA（`ppa:elmer-csc-ubuntu/elmer-csc-ppa`）またはソースビルド |
| code_aster | Salome-Meca のコンテナ版（Singularity/Apptainer） |
| OpenFOAM | openfoam.com または openfoam.org の apt リポジトリ |
| OpenModelica | 公式aptリポジトリ（Windows版インストーラも可） |
| Gmsh / ParaView | apt または公式バイナリ |
| Python系（ROSS, python-control, pyleecan, pychrono, sdypy, scipy） | `pip` / `conda`（pychrono は conda-forge） |
| FEMM | Windows専用（WSLからは Windows 側で実行） |

---

## 6. 参考になる理論書・キーワード

- Y. Altintas『Manufacturing Automation』— 切削力・びびり・送り駆動系制御の定番
- 工作機械の熱変位、送り駆動系の剛性 — 日本機械学会・精密工学会の論文
- ボールねじ・リニアガイドの選定 — 各メーカー技術資料（THK、NSK、HIWIN 等の総合カタログに剛性・寿命・危険速度の計算式が載っている）
- 主軸テーパ規格 → [study/ref001](../study/ref001_工作機械_主軸テーパ30-40-50番.md)

---

[← docs 目次](README.md) ｜ [⌂ トップ](../README.md) ｜ [次: ref005 ばね・固有値・周波数応答の設定比較 →](ref005_CalculiX-FrontISTR_ばね・固有値・周波数応答の設定比較.md)
