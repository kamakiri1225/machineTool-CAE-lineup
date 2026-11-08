#!/bin/bash
set -e
cd "$(dirname "$0")"
~/opt/venv/bin/python make_inputs.py
R=$HOME/opt/ccxroot/usr/lib/x86_64-linux-gnu
for k in rbe2 rbe3; do
  ( cd $k/ccx && ~/opt/ccx.sh -i machine > ccx.log 2>&1 )
  ( cd $k/fistr && LD_LIBRARY_PATH=$R/blas:$R/lapack OMP_NUM_THREADS=8 ~/opt/frontistr/bin/fistr1 > fistr.log 2>&1 )
done
