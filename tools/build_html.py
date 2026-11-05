#!/usr/bin/env python3
"""Build browsable HTML from every docs/*.md in the repository (pandoc + MathJax), and a top index.html.
Links to .md are rewritten to .html. Images are referenced in place (relative paths), so the HTML works both
locally (open index.html) and on GitHub Pages.   usage: python3 tools/build_html.py"""
import os, re, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSS = """<style>
body{font-family:'Noto Sans JP','Hiragino Sans','Meiryo',sans-serif;max-width:1100px;margin:2em auto;padding:0 1.2em;line-height:1.75;color:#222}
h1{border-bottom:3px solid #1F3A5F;padding-bottom:.3em;color:#1F3A5F} h2{border-left:6px solid #2F6FB0;padding-left:.5em;margin-top:2em;color:#1F3A5F}
h3{color:#2F6FB0} img{max-width:100%;height:auto;border:1px solid #ddd;margin:.5em 0}
table{border-collapse:collapse;margin:1em 0;font-size:.95em} th,td{border:1px solid #bbb;padding:.35em .7em} th{background:#E3ECF7}
code{background:#f4f4f4;padding:.1em .3em;border-radius:3px} pre{background:#f4f4f4;padding:1em;overflow-x:auto}
nav.top{font-size:.9em;margin-bottom:1em} nav.top a{color:#2F6FB0}
</style>"""


def build(md):
    rel = os.path.relpath(ROOT, os.path.dirname(md))
    html = md[:-3] + ".html"
    src = open(md, encoding="utf-8").read()
    src = re.sub(r"\]\(([^)]+?)\.md\)", r"](\1.html)", src)
    # join multi-line display math so that lines starting with "+" or "-" are not read as list items
    src = re.sub(r"\$\$(.+?)\$\$", lambda m: "$$" + " ".join(l.strip() for l in m.group(1).strip().splitlines()) + "$$", src, flags=re.S)
    title = (re.search(r"^# (.+)$", src, re.M) or [None, os.path.basename(md)[:-3]])[1]
    tmp = md + ".tmp"
    open(tmp, "w", encoding="utf-8").write(src)
    header = f'{CSS}\n<nav class="top"><a href="{rel}/index.html">← 一覧に戻る</a></nav>'
    hp = md + ".hdr"
    open(hp, "w", encoding="utf-8").write(header)
    subprocess.run(["pandoc", tmp, "-f", "gfm+tex_math_dollars", "-t", "html5", "-s", "--mathjax=https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml-full.js",
                    "--metadata", f"title={title}", "--include-in-header", hp, "-o", html], check=True)
    os.remove(tmp); os.remove(hp)
    h = open(html, encoding="utf-8").read()
    h = re.sub(r'<script src="https://polyfill\.io[^"]*"></script>\s*', "", h)   # polyfill.io is no longer safe/available
    open(html, "w", encoding="utf-8").write(h)
    return html, title


items = []
for d, _, fs in sorted(os.walk(ROOT)):
    if "/.git" in d or "/processor" in d:
        continue
    for f in sorted(fs):
        if f.endswith(".md") and ("/docs" in d or d == ROOT):
            if d == ROOT and f == "README.md":
                continue
            h, t = build(os.path.join(d, f))
            items.append((os.path.relpath(h, ROOT), t))
groups = {}
for h, t in items:
    groups.setdefault(h.split("/")[0], []).append((h, t))
NAMES = {"01_coolantCHT": "01 加工エリアのクーラントと熱（気液二相流＋固体の熱連成）",
         "02_laserAssistGas": "02 レーザー切断のアシストガス（圧縮性・超音速、混合ガス）",
         "03_fumeExtraction": "03 切断ヒュームの吸い込み（流れ＋粒子）",
         "04_factoryAmbient": "04 工場の気温と機械の温度（流れ＋固体の熱連成）",
         "05_spindleBearingHeat": "05 主軸と軸受の発熱・冷却（固体の熱伝導）",
         "06_mistChips": "06 加工室のミストと切りくず（流れ＋粒子）",
         "07_airbagTankTest": "07 【参考・テスト用】エアバッグ用インフレーターのタンク試験（燃焼ガスの湧き出し・チョーク・混合ガス）"}
TIME = """## 計算時間の一覧

計算環境：Intel Xeon Gold 6212U 2.40 GHz（WSL2）、OpenFOAM v2512。並列は `export HWLOC_COMPONENTS=-gl` を付けた `mpirun`。

| フォルダ / ケース | ソルバ | セル数 | 並列数 | 実時間 |
|---|---|---|---|---|
| 01 クーラント run001（公式例題 2D、0.2 s） | chtMultiRegionTwoPhaseEulerFoam | 例題の 2D メッシュ | 1 | 12 分 |
| 01 クーラント run004（加工エリア 3D、3 s） | chtMultiRegionTwoPhaseEulerFoam | 流体 17,412 ＋ 固体 588 | 1 | 12 分 |
| 02 アシストガス run002（酸素切断、軸対称） | reactingFoam（LTS） | 10,525 | 1 | 24 分 |
| 02 アシストガス run003（高圧空気切断、軸対称） | reactingFoam（LTS） | 12,270 | 1 | 42 分 |
| 02 アシストガス run004（3D 切断溝、N₂、0.5 / 2.0 MPa） | reactingFoam（LTS） | 588,724 | 4（× 2 本同時） | 各 約 7.5 時間（4,000 回） |
| 03 ヒューム run001 流れ ＋ 粒子（Q = 0.30） | simpleFoam ＋ icoUncoupledKinematicParcelFoam | 55,296 | 1 | 6 分 ＋ 3 分 |
| 03 ヒューム run002 / run003 流れ（Q = 0.15 / 0.60） | simpleFoam | 55,296 | 1 | 8.5 分 / 11 分 |
| 04 工場の気温 run001（流れ 120 s ＋ 温度 24 h） | chtMultiRegionFoam | 26,880 | 1 | 23 分 ＋ 3 分 |
| 05 主軸の発熱 run001（1 h） | solidFoam | 1,584 | 1 | 18 秒 |
| 06 ミスト・切りくず run001（流れ ＋ 粒子 2 種） | simpleFoam ＋ icoUncoupledKinematicParcelFoam | 73,728 | 1 | 約 25 分 ＋ 4 分 ＋ 3 分 |
| 07 【参考】インフレーターのタンク試験 run001（50 ms、非定常） | reactingFoam（軸対称） | 10,700 | 4 | 約 5.6 時間 |

- 2 万セル程度までは 1 並列でも数十分。3D で数十万セルになると、4 並列でも数時間かかる。
- 02 の run004 の 4.0 / 8.0 MPa は、計算が途中で進まなくなったため条件を見直し中。
"""
body = ["# 工作機械まわりの CFD（OpenFOAM v2512）", "",
        "工作機械まわりの流体・熱の現象を、一般的な寸法や論文の条件をなぞった簡易モデルで計算した記録です。",
        "各ケースの解析設定（`system/`・`constant/`・`0.orig/` とケースを作る Python スクリプト）、解説、図とアニメーションを置いています。", "", TIME]
for g in sorted(groups):
    body.append(f"## {NAMES.get(g, g)}"); body.append("")
    for h, t in groups[g]:
        body.append(f"- [{t}]({h})")
    body.append("")
open(os.path.join(ROOT, "index.md"), "w", encoding="utf-8").write("\n".join(body))
hp = os.path.join(ROOT, "index.hdr"); open(hp, "w", encoding="utf-8").write(CSS)
subprocess.run(["pandoc", os.path.join(ROOT, "index.md"), "-f", "gfm", "-t", "html5", "-s", "--metadata", "title=工作機械まわりの CFD",
                "--include-in-header", hp, "-o", os.path.join(ROOT, "index.html")], check=True)
os.remove(hp); os.remove(os.path.join(ROOT, "index.md"))
print(len(items), "pages ->", os.path.join(ROOT, "index.html"))
