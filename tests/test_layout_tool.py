"""The layout-to-recipe tool on the live page (docs/assets/layout-tool.js) against the same arithmetic in Python and
against the chamber runs logged in paper/assumptions.md (Night 3, part 1).

The Python half always runs. The JavaScript half evaluates docs/assets/layout-tool.js in headless Chromium through
Playwright for Python and is skipped when Playwright or its browser is not installed."""
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
JS = ROOT / "docs" / "assets" / "layout-tool.js"
KB = 1.380649e-23
CHAMBER = dict(x=1.0, Q=300.0, p=300.0, T=300.0, N0=4e18)          # baseline: 1 % TMA, 300 sccm, 300 Pa, 300 C
PRESETS = {
    "dense": [(0, 75, 100, 100, 1.78e13), (75, 100, 100, 100, 8.48e12), (100, 150, 100, 0, 0)],
    "uniform": [(0, 150, 100, 100, 1.78e13)],
    "ramp": [(30 * i, 30 * (i + 1), 100, ar, 1.78e13) for i, ar in enumerate((100, 90, 80, 70, 60))],
}


def supply(ch):
    x, p, T = ch["x"] / 100, ch["p"], ch["T"] + 273.15
    c0 = x * p / (KB * T)
    U = ch["Q"] * 1.6667e-8 * (101325 / p) * (T / 273.15) / (math.pi * 0.15 ** 2)
    return c0, U, ch["N0"] / (U * c0)


def compute(zones, ch, wafer_r=150.0):
    """Same arithmetic as LayoutTool.compute: per-zone wall ratio and local rule, area-weighted average, margin band."""
    c0, U, ts0 = supply(ch)
    rows, wall_area = [], 0.0
    for r_in, r_out, d, ar, den in zones:
        bare = not (den > 0 and ar > 0 and d > 0)
        wall = 0.0 if bare else den * math.pi * (d * 1e-9) * ar * (d * 1e-9)
        area = math.pi * (min(r_out, wafer_r) ** 2 - min(r_in, wafer_r) ** 2)
        wall_area += wall * area
        rows.append(dict(bare=bare, wall=wall, local=1.2 * (1 + wall) * ts0))
    holes = [r for r in rows if not r["bare"]]
    gov = max(holes, key=lambda r: r["local"])
    avg_wall = wall_area / (math.pi * wafer_r ** 2)
    avg_rule = 1.2 * (1 + avg_wall) * ts0
    covered = sum(math.pi * (min(b, wafer_r) ** 2 - min(a, wafer_r) ** 2) for a, b, *_ in zones)
    walls = [r["wall"] for r in holes]
    uniform = len(holes) == len(rows) and covered >= math.pi * wafer_r ** 2 * (1 - 1e-3) and (max(walls) - min(walls)) <= 1e-2 * max(walls)
    lo, hi = (0.98, 0.98) if uniform else (0.82, 0.87)
    return dict(recipe=gov["local"], lo=gov["local"] * lo, hi=gov["local"] * hi, avg_wall=avg_wall, avg_rule=avg_rule,
                under=gov["local"] / avg_rule, uniform=uniform, rows=rows, U=U, c0=c0, ts0=ts0)


def test_python_reproduces_the_logged_runs():
    """Night 3 part 1: uniform 55.9 / 14.39 s (14.05 s measured); dense centre avg 19.2 / 5.10 s against 11.85 s needed,
    local 14.39 s; AR ramp avg 41.0 / 10.62 s against 12.50 s, local 14.36 s. The log rounds U to 0.0501 and c0 to
    3.79e20, which puts its rule 0.01 s above the exact arithmetic."""
    d = compute(PRESETS["dense"], CHAMBER)
    assert d["rows"][0]["wall"] == pytest.approx(55.92, abs=0.01)
    assert d["recipe"] == pytest.approx(14.39, abs=0.02)
    assert d["avg_wall"] == pytest.approx(19.2, abs=0.1)
    assert d["avg_rule"] == pytest.approx(5.10, abs=0.05)
    assert (d["lo"], d["hi"]) == (pytest.approx(11.8, abs=0.05), pytest.approx(12.5, abs=0.05))
    assert d["under"] == pytest.approx(2.82, abs=0.02)
    assert d["lo"] <= 11.85 <= d["hi"]                                # the measured dense-centre time sits inside the band
    u = compute(PRESETS["uniform"], CHAMBER)
    assert u["uniform"] and u["recipe"] == pytest.approx(14.39, abs=0.02) and u["lo"] == pytest.approx(14.09, abs=0.02)
    assert abs(14.05 - u["lo"]) < 0.05
    r = compute(PRESETS["ramp"], CHAMBER)
    assert r["recipe"] == pytest.approx(14.36, abs=0.03)
    assert r["avg_rule"] == pytest.approx(10.62, abs=0.25)          # 5-zone staircase of the linear ramp: 10.42 s
    assert r["rows"][-1]["local"] == pytest.approx(8.76, abs=0.05)  # AR 60 edge: 8.73 s against the log's 8.76 s at AR 60.2
    assert r["lo"] <= 12.50 <= r["hi"]


def test_margin_band_is_two_runs_wide():
    """A 40 % depth variation across the wafer gives a 39 % spread in the rule, as in run (b)."""
    r = compute(PRESETS["ramp"], CHAMBER)
    assert 1 - r["rows"][-1]["local"] / r["rows"][0]["local"] == pytest.approx(0.39, abs=0.01)


@pytest.fixture(scope="module")
def js_tool():
    pw = pytest.importorskip("playwright.sync_api")
    with pw.sync_playwright() as p:
        try:
            browser = p.chromium.launch()
        except Exception as e:  # browser not installed
            pytest.skip(f"Chromium not available: {e}")
        page = browser.new_page()
        page.set_content("<!doctype html><html><body></body></html>")
        page.add_script_tag(path=str(JS))
        page.wait_for_function("() => window.LayoutTool && window.LayoutTool.compute")
        yield page
        browser.close()


def js_zones(zs):
    return [dict(rIn=a, rOut=b, d=d, AR=ar, density=den) for a, b, d, ar, den in zs]


@pytest.mark.parametrize("name", list(PRESETS))
def test_javascript_matches_python(js_tool, name):
    r = js_tool.evaluate("([z, c]) => { const r = window.LayoutTool.compute(z, c); return {recipe: r.recipe.local, lo: r.recipe.lo, hi: r.recipe.hi, avgWall: r.average.wall, avgRule: r.average.rule, under: r.underDose, uniform: r.uniform, walls: r.zones.map(z => z.wall), locals: r.zones.map(z => z.local), U: r.supply.U, c0: r.supply.c0}; }",
                        [js_zones(PRESETS[name]), CHAMBER])
    py = compute(PRESETS[name], CHAMBER)
    assert r["U"] == pytest.approx(py["U"], rel=1e-9) and r["c0"] == pytest.approx(py["c0"], rel=1e-9)
    for k, pk in (("recipe", "recipe"), ("lo", "lo"), ("hi", "hi"), ("avgWall", "avg_wall"), ("avgRule", "avg_rule"), ("under", "under")):
        assert r[k] == pytest.approx(py[pk], rel=1e-9), k
    assert r["uniform"] == py["uniform"]
    assert r["walls"] == pytest.approx([z["wall"] for z in py["rows"]], rel=1e-9)
    assert r["locals"] == pytest.approx([z["local"] for z in py["rows"]], rel=1e-9)


def test_javascript_presets_and_csv_round_trip(js_tool):
    presets = js_tool.evaluate("() => Object.fromEntries(Object.entries(window.LayoutTool.presets).map(([k, p]) => [k, p.zones.map(z => [z.rIn, z.rOut, z.d, z.AR, z.density])]))")
    assert {k: [tuple(z) for z in v] for k, v in presets.items()} == {k: [tuple(map(float, z)) for z in v] for k, v in PRESETS.items()}
    back = js_tool.evaluate("() => { const L = window.LayoutTool, z = L.presets.ramp.zones, p = L.parseCsv(L.toCsv(z)); return [p.errors, p.zones.map(q => [q.rIn, q.rOut, q.d, q.AR, q.density])]; }")
    assert back[0] == [] and [tuple(z) for z in back[1]] == [tuple(map(float, z)) for z in PRESETS["ramp"]]
    parsed = js_tool.evaluate("() => window.LayoutTool.parseCsv('r_in_mm, r_out_mm, d_nm, AR, density_per_m2\\n# comment\\n0, 75, 100, 100, 1.78e13\\n75, 60, 100, 100, 1e13\\nabc\\n')")
    assert len(parsed["zones"]) == 1 and len(parsed["errors"]) == 2
    warn = js_tool.evaluate("() => window.LayoutTool.compute([{rIn:0, rOut:150, d:100, AR:150, density:1.78e13}], {x:0.2, Q:300, p:300, T:300, N0:4e18}).warnings")
    assert any("AR 150" in w for w in warn) and any("TMA 0.2" in w for w in warn)
