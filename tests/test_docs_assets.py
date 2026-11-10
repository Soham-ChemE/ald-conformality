"""The live page (docs/index.html) ships with every local asset it references. Standard library only, so the CI
pages-check job can run it with a bare interpreter:  python -m pytest -q tests/test_docs_assets.py  or  python
tests/test_docs_assets.py."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
INDEX = DOCS / "index.html"


def local_refs(html):
    """Local (non-http, non-data) src= and href= targets, plus the two fallback PNGs of every data-viz-fig."""
    refs = set()
    for m in re.finditer(r'\b(?:src|href)="([^"]+)"', html):
        u = m.group(1)
        if u.startswith(("http://", "https://", "data:", "#", "mailto:")):
            continue
        refs.add(u.split("#")[0].split("?")[0])
    for m in re.finditer(r'data-viz-fig="([^"]+)"', html):
        refs.add(f"assets/{m.group(1)}_light.png")
        refs.add(f"assets/{m.group(1)}_dark.png")
    return sorted(r for r in refs if r)


def missing_assets():
    html = INDEX.read_text(encoding="utf-8")
    return [r for r in local_refs(html) if not (DOCS / r).is_file()]


def test_index_exists():
    assert INDEX.is_file(), INDEX


def test_index_references_its_script_style_and_data():
    html = INDEX.read_text(encoding="utf-8")
    refs = local_refs(html)
    assert "assets/site.js" in refs and "assets/site.css" in refs and "assets/data.js" in refs, refs


def test_every_local_asset_is_present():
    assert missing_assets() == []


if __name__ == "__main__":
    miss = missing_assets()
    print(f"{INDEX.relative_to(ROOT)}: {len(local_refs(INDEX.read_text(encoding='utf-8')))} local references, "
          f"{len(miss)} missing")
    for m in miss:
        print("  missing:", m)
    sys.exit(1 if miss else 0)
