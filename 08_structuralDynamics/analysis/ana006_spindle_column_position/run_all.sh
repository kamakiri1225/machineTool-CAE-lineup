#!/bin/bash
# ana006 の全工程：モデル作成・固有値解析 → 安定限界 → 図
set -e; cd "$(dirname "$0")"; PY=~/opt/venv/bin/python
$PY 01_build_and_run.py
$PY 02_modal_lobes.py
$PY 03_render_cases.py
$PY 04_plot_sweeps.py
$PY 05_render_modes.py
$PY 06_tool_work_split.py
echo DONE
