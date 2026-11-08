#!/bin/bash
# ana002 実行スクリプト（CalculiX と FrontISTR）
set -e
cd "$(dirname "$0")"
~/opt/venv/bin/python make_inputs.py
( cd ccx && ~/opt/ccx.sh -i machine > ccx.log 2>&1 )
R=$HOME/opt/ccxroot/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH=$R/blas:$R/lapack OMP_NUM_THREADS=8
cd fistr
cp hecmw_ctrl_eigen.dat hecmw_ctrl.dat && ~/opt/frontistr/bin/fistr1 > eigen.out 2>&1 && mv 0.log eigen_0.log
for d in 1 2 3; do
  cp hecmw_ctrl_freq_$d.dat hecmw_ctrl.dat && ~/opt/frontistr/bin/fistr1 > freq_$d.out 2>&1
  mkdir -p out_$d && mv -f 0.log dyna_* *.txt out_$d/ 2>/dev/null || true
done
