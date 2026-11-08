#!/bin/bash
# テストバー（60×60×100 mm 鋼）で ana005 の解析を全部やり直す。旧（細い工具）の結果は results_slender_tool/ に保存済み
set -e
FIRST=${FIRST-run_relative_frf.py}
cd "$(dirname "$0")"; export VMC_TOOL=bar PY=~/opt/venv/bin/python
for s in $FIRST analyze.py check_static.py export_modal_basis.py render_model.py plot_relative_matrix.py milling.py \
         compliance_modes.py stability_compare.py chatter_modes.py chatter_by_rpm.py chatter_point_anim.py anim_compliance_modes.py peak_modes_anim.py animate_cutting.py; do
  echo "=== $s $(date +%T)"; $PY $s > results/log_$s.txt 2>&1 || { echo "FAILED $s"; tail -20 results/log_$s.txt; exit 1; }
done
echo DONE
