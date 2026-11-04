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
         "04_factory_spindle_mist": "04 工場の気温・主軸の発熱・ミストと切りくず"}
body = ["# 工作機械まわりの CFD（OpenFOAM v2512）", "",
        "工作機械まわりの流体・熱の現象を、一般的な寸法や論文の条件をなぞった簡易モデルで計算した記録です。",
        "各ケースの解析設定（`system/`・`constant/`・`0.orig/` とケースを作る Python スクリプト）、解説、図とアニメーションを置いています。", ""]
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
