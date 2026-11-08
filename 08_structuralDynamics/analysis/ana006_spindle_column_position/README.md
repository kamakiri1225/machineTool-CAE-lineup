# ana006 主軸頭の高さ（Z）とコラム位置（Y）を変えたときのびびり

説明は [docs/ana006_CalculiX_主軸頭の高さとコラム位置のびびりパラスタ.md](../../docs/ana006_CalculiX_主軸頭の高さとコラム位置のびびりパラスタ.md) にあります。

## フォルダ構成

```
ana006_spindle_column_position/
├── params.py              寸法・ばね定数・切削条件・ケース表（位置を変えるときはここだけ書き換える）
├── lib/
│   ├── geometry.py        (dy, dz) → 部品ごとの直方体と、部品どうしの接合面
│   ├── fem.py             部品ごとに Hex20 メッシュ → 接合面を「代表点＋MPC＋ばね」でつなぐ → CalculiX 入力
│   ├── modal.py           .frd からモードを読み、刃先―ワークの相対変位 d_r を作る
│   ├── chatter.py         相対 FRF（モード和）・安定限界線図・モードの割合（ana005 と同じ式）
│   └── view.py            pyvista 表示用の格子
├── 01_build_and_run.py    全ケースのモデルを作って CalculiX で固有値解析
├── 02_modal_lobes.py      相対 FRF と安定限界 → results/case_*.npz, summary.json
├── 03_render_cases.py     ケースのモデル図 → results/cases_model.png
├── 04_plot_sweeps.py      スイープのグラフ・表 → results/sweep_*.png, summary.md
├── 05_render_modes.py     最小点で効くモードの形 → results/modes_at_min.png, mode_base.gif
├── run_all.sh             上の 01〜05 を順に実行
├── cases/<ケース名>/      model.inp（入力）、model.frd/dat（結果）、model_info.json（接合面・節点番号）
└── results/               図・表・数値
```

## 部品とつなぎ方

| 部品 | 動かし方 | 下の部品とのつなぎ |
|---|---|---|
| ベッド | 動かない | 基礎の脚 4 か所 → 地面（ばね） |
| テーブル | （X は中央で固定） | X 案内：キャリッジ 4 個＋ボールねじナット |
| コラム | Y に `dy` | Y 案内：キャリッジ 4 個＋ボールねじナット |
| 主軸頭 | Y に `dy`、Z に `dz` | Z 案内（コラム前面）：キャリッジ 4 個＋ボールねじナット |

部品どうしは節点を共有しません。接合面ごとに両側の面の節点を代表点へ剛体結合（`*EQUATION`）し、
代表点どうしを並進ばね・回転ばね（`SPRING2`）でつなぎます。部品を動かすと、接合面も一緒に動きます。

## 実行

```bash
bash run_all.sh                      # 全部（7 ケース × 約 110 s）
~/opt/venv/bin/python 01_build_and_run.py z100   # 1 ケースだけ
```

ケースを足すには `params.py` の `CASES`（と `SWEEPS`）に `"名前": (dy, dz)` を追加します。
