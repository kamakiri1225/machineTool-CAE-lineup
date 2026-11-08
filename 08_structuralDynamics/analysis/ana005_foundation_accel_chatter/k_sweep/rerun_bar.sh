#!/bin/bash
# テストバーで 4.5〜4.6 節の脚剛性スイープをやり直す
set -e; cd "$(dirname "$0")"; export VMC_TOOL=bar; PY=~/opt/venv/bin/python
$PY run_k_sweep.py
for c in k025 k1 k4; do
  mkdir -p $c/eig
  $PY - $c <<'PYEOF'
import sys,re
c=sys.argv[1]; t=open(f"{c}/model.inp").read()
i=t.index("*STEADY STATE DYNAMICS"); t=t[:t.rindex("*STEP",0,i)]
t=re.sub(r"(\*FREQUENCY[^\n]*\n)[^\n]*", r"\g<1>10, 1.0", t)
open(f"{c}/eig/model.inp","w").write(t)
PYEOF
  (cd $c/eig && ~/opt/ccx.sh -i model > ccx.log 2>&1)
  rm -f $c/model.frd
done
$PY post.py; $PY anim_modes.py; $PY transient_modal.py
echo DONE
