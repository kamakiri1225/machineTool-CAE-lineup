"""Insert numerical results into the human-readable report."""
from pathlib import Path
import json,numpy as np
D=Path(__file__).resolve().parent;R=D/'results';doc=D.parent.parent/'docs/ana005_基礎ばねモデルの加減速と再生びびり分析.md'
s=doc.read_text();a=json.loads((R/'acceleration.json').read_text());checks=json.loads((R/'static_checks.json').read_text());mill=json.loads((R/'milling_checks.json').read_text());shapes=json.loads((R/'cutting_shapes.json').read_text())
s=s.replace('<!-- BASE_PARAMETERS -->',f"全質量M={a['M_kg']:.2f} kg、テーブル＋ワーク側質量mt={a['moving_mass_kg']:.3f} kg、k={a['k_N_m']/1e6:.1f} MN/m、c={a['c_Ns_m']:.1f} N·s/m。近似並進固有振動数{a['fn_Hz']:.2f} Hz、減衰比{a['zeta']*100:.2f}%となる。")
s=s.replace('<!-- ACCEL_RESULTS -->',(R/'acceleration.md').read_text()+"\n最大加速度2 m/s²の例では、ジャーク立上げ0→40 msで指令移動時間は0.350→0.390 sとなる一方、停止後の本体変位最大は6.566→0.612 µmとなり、±1 µm基準の移動＋本体整定は0.465→0.390 sとなった。この基礎モデル・指令条件での比較値である。")
labels={'fixed':'底面固定','foundation':'基礎ばね','foundation_table_stiff':'テーブル結合10倍'}
rows=['| ケース | 方向 | 静解析 µm/N | 5 Hz FRF実部 µm/N | 比較差 % |','|---|---|---:|---:|---:|']
for name,c in checks.items():
 for j in range(3):
  v=c['static_diagonal_mm_N'][j];w=c['dynamic_5Hz_diagonal_mm_N'][j];rows.append(f'| {labels[name]} | {"XYZ"[j]} | {v*1000:.6f} | {w*1000:.6f} | {abs(w/v-1)*100:.2f} |')
s=s.replace('<!-- STATIC_CHECKS -->','\n'.join(rows)+"\n\nXYの比較差は約2.7〜3.0%。Zは約14〜16%の差が残っており、今回のフライス安定限界はXY行列に限定する。5 Hzと静解析の違いも含む比較なので、この表だけで厳密なモード数収束を証明したとは扱わない。FEMメッシュ・モード数・周波数点の収束確認は今後の検証項目。")
rows=['| ケース | 掃引範囲内の最小切込み限界 mm | 約6,000 min⁻¹での限界 mm | そのときのびびり周波数 Hz |','|---|---:|---:|---:|']
for name,label in labels.items():
 c=mill[name];p=c['at_6000_rpm'];rows.append(f"| {label} | {c['minimum_axial_depth_mm']:.3f} | {p['axial_depth_mm']:.3f} | {p['chatter_frequency_Hz']:.1f} |")
s=s.replace('<!-- MILLING_RESULTS -->','\n'.join(rows)+"\n\n約6,000 min⁻¹の数値は実際には6,011.69 min⁻¹の計算点。今回の境界のびびり周波数は主に工具側ピークの1.7〜1.8 kHz付近にあり、テーブル結合10倍による限界の変化は小さい。**この条件では、テーブル起因のびびりが加工限界を決めているとは判断しない。** 工具側が先に限界を決める条件になっている。テーブル側に焦点を当てる次の比較では、工具突出しを短くする等により工具曲げの影響を変え、同じ手順で境界とモード寄与を調べる。")
needle='### 3.4 静解析との照合'
addition='''### 3.4 切削力対に対する全体振動

X方向切削力対のピークを、100〜600 Hzの範囲と800〜3,000 Hzの範囲からそれぞれ選び、保存した固有モードを使って全節点の複素変位を改めて出力した。モデルの切り出し・分解は行わず、変位だけを拡大した定常応答を表示する。下段は相対XYZの3曲線で、縦線と丸印が上段の応答周波数。

![181.8 Hzの全体振動と相対FRF](../analysis/ana005_foundation_accel_chatter/results/cutting_low.gif)

![1831.4 Hzの全体振動と相対FRF](../analysis/ana005_foundation_accel_chatter/results/cutting_high.gif)

| 応答周波数 Hz | ワーク側X振幅 µm/N | 刃先側X振幅 µm/N | 相対X振幅 µm/N |
|---:|---:|---:|---:|
'''
for sh in shapes:addition+=f"| {sh['frequency_Hz']:.2f} | {sh['work_amplitude_um_N'][0]:.6f} | {sh['tool_amplitude_um_N'][0]:.6f} | {sh['relative_amplitude_um_N'][0]:.6f} |\n"
addition+='\n181.8 Hzではワーク側の動きが大きく、テーブル結合を変えた感度と併せてテーブル側の寄与を調べられる。1,831.4 Hzではワーク側は約0.00024 µm/N、刃先側は約0.339 µm/Nで、工具側の動きが支配的。振幅を単純加算・減算せず、複素差で相対応答を求めている。\n\n### 3.5 静解析との照合'
s=s.replace(needle,addition)
s=s.replace('テーブル結合だけを強くした場合の変化を、','![相対FRFの9成分](../analysis/ana005_foundation_accel_chatter/results/relative_matrix.png)\n\n行は相対変位の観測方向、列は切削力対の加振方向。各図の縦軸は線形で、上限は成分ごとに異なる。\n\nテーブル結合だけを強くした場合の変化を、')
s=s.replace('~/opt/venv/bin/python milling.py','~/opt/venv/bin/python milling.py\n~/opt/venv/bin/python plot_relative_matrix.py\n~/opt/venv/bin/python animate_cutting.py')
doc.write_text(s)
print('Updated',doc)
