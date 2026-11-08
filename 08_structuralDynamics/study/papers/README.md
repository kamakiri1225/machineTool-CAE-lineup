[← study 目次](../README.md) ｜ [⌂ トップ](../../README.md)

# 論文解説（papers）

工作機械・構造最適化・振動に関する**論文を 1 本ずつ深掘りしたノート**を置く。
自分でまとめた知識のノート（`study/refNNN_…`）とは分けている。

## 命名ルール

- ファイル名：`paper連番_第1著者年_内容.md`（連番は 3 桁、通し番号）
  - 例：`paper001_Law2013_位置依存のびびり安定性.md`
- 図：`images/paper連番/`。生成スクリプトは `study/scripts/make_paper連番_*.py`
- 論文の PDF：`pdf/` に置く。**著作物なので公開しない**（ブログ・社外資料に入れない）。ノートには DOI と出版社・著者公開版のリンクを書く
- 各ノートの構成：書誌情報 → ひとことで → 手法の流れ（図）→ 数値 → 学べること・どこまで使えるか → 自分のモデルで試す

## 一覧

| ファイル | 論文 | 内容 |
|---|---|---|
| [paper001_Law2013_位置依存のびびり安定性.md](paper001_Law2013_位置依存のびびり安定性.md) | Law, Altintas, Phani, *IJMTM* 68 (2013) | 部品ごとの縮約＋代表点で速い全機モデル（1 位置 1 分）、位置・送り方向ごとのびびり安定限界、補強 → トポロジー最適化で質量を戻す。自分のモデルで安定限界を再現 |

## 予定

| 予定 | 論文 | 状態 |
|---|---|---|
| paper002 | A comprehensive framework for computationally efficient system-level design optimization of machine tools（*J. Manuf. Syst.*, 2026、[ScienceDirect](https://www.sciencedirect.com/science/article/pii/S027861252600035X)） | PDF を `pdf/` に置いてもらえれば作成 |
