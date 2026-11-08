"""docs / study の各ページに「目次へ戻る・前・次」のリンクを付ける（何度実行しても同じ結果）"""
import os, glob, re

ROOT = os.path.dirname(os.path.abspath(__file__))
SETS = {
    "docs": ("docs 目次", sorted(glob.glob(os.path.join(ROOT, "docs", "ref*.md"))) + sorted(glob.glob(os.path.join(ROOT, "docs", "ana*.md")))),
    "study": ("study 目次", sorted(glob.glob(os.path.join(ROOT, "study", "ref*.md")))),
    "study/papers": ("論文解説の目次", sorted(glob.glob(os.path.join(ROOT, "study", "papers", "paper*.md")))),
}


def short(path):
    b = os.path.basename(path)[:-3]
    return b.split("_")[0] + " " + "_".join(b.split("_")[2:]) if b.count("_") >= 2 else b


for folder, (toc, files) in SETS.items():
    for i, f in enumerate(files):
        top = os.path.relpath(os.path.join(ROOT, "README.md"), os.path.dirname(f)).replace(os.sep, "/")
        nav = [f"[← {toc}](README.md)", f"[⌂ トップ]({top})"]
        if i > 0:
            nav.append(f"[← 前: {short(files[i-1])}]({os.path.basename(files[i-1])})")
        if i < len(files) - 1:
            nav.append(f"[次: {short(files[i+1])} →]({os.path.basename(files[i+1])})")
        line = " ｜ ".join(nav)
        lines = open(f, encoding="utf-8").read().splitlines()
        # 既存のナビ行（先頭・末尾）を取り除く
        while lines and (lines[0].startswith("[← ") or lines[0].startswith("[⌂") or lines[0].strip() == ""):
            lines.pop(0)
        while lines and (lines[-1].startswith("[← ") or lines[-1].startswith("[⌂") or lines[-1].strip() in ("", "---")):
            lines.pop()
        out = [line, ""] + lines + ["", "---", "", line, ""]
        open(f, "w", encoding="utf-8").write("\n".join(out))
        print(os.path.relpath(f, ROOT))
