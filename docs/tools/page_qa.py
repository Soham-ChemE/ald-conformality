#!/usr/bin/env python3
"""Page QA and fallback snapshots for the live page (Playwright for Python; driven by docs/build_assets.py).

  python docs/tools/page_qa.py --docs docs [--snap] [--qa] [--lighthouse] [--hero-png] [--hero-json out.json] [--pages index.html]

  --snap        render every [data-viz] scene at full progress in both themes to assets/viz_<name>_<theme>.png (static fallback)
  --qa          screenshots in both themes: hero, the first two body scenes, a 375 px hero, plus a console-error check, to qa/<theme>/
  --hero-png    assets/hero.png, the hero canvas alone at the first keyframe (dark theme)
  --hero-json   projected particle positions from window.__hero (used by the profile page's README hero)
  --lighthouse  desktop Lighthouse run against the served page (lighthouse CLI on PATH), report to qa/lighthouse.report.{html,json}

Serves docs/ on a local port, opens it in headless Chromium with software WebGL (so frame timing is not representative) and
window.__forceGL set so the live scenes render rather than the static path. Mirrors docs/tools/page-qa.mjs of the sibling repos."""
from __future__ import annotations
import argparse, http.server, json, os, shutil, socketserver, subprocess, sys, threading, time
from functools import partial
from pathlib import Path

ap = argparse.ArgumentParser(); ap.add_argument("--docs", default="docs"); ap.add_argument("--pages", default="index.html")
for f in ("--snap", "--qa", "--lighthouse", "--hero-png"): ap.add_argument(f, action="store_true")
ap.add_argument("--hero-json"); A = ap.parse_args()
DOCS = Path(A.docs).resolve(); PAGES = A.pages.split(",")

class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *a): pass
    def end_headers(self): self.send_header("cache-control", "no-cache"); super().end_headers()
srv = socketserver.TCPServer(("127.0.0.1", 0), partial(Quiet, directory=str(DOCS))); PORT = srv.server_address[1]
threading.Thread(target=srv.serve_forever, daemon=True).start(); BASE = f"http://127.0.0.1:{PORT}/"

from playwright.sync_api import sync_playwright
qa_json = DOCS / "qa" / "qa.json"; report = json.loads(qa_json.read_text()) if qa_json.exists() else {"pages": {}}; report.setdefault("pages", {})

def open_page(browser, page_name, theme, viewport):
    ctx = browser.new_context(viewport=viewport, device_scale_factor=2 if viewport["width"] < 500 else 1, reduced_motion="no-preference", color_scheme=theme)
    ctx.add_init_script(f"window.__forceGL = true; try {{ localStorage.setItem('theme', '{theme}'); }} catch (e) {{}}")
    page = ctx.new_page(); errors, failed = [], []
    page.on("console", lambda m: errors.append(m.text) if m.type == "error" else None); page.on("pageerror", lambda e: errors.append(str(e))); page.on("requestfailed", lambda r: failed.append(r.url))
    page.goto(BASE + page_name, wait_until="networkidle"); page.wait_for_function("() => window.SITE && window.SITE.ready"); time.sleep(1.2)
    return ctx, page, errors, failed

def scroll_to(page, sel, frac=.5):
    page.evaluate("""([s, fr]) => { const el = document.querySelector(s); const r = el.getBoundingClientRect(); const y = r.top + scrollY - (innerHeight - r.height) * fr;
      if (window.SITE && window.SITE.lenis) window.SITE.lenis.scrollTo(y, { immediate: true, force: true }); window.scrollTo(0, y); if (window.ScrollTrigger) ScrollTrigger.update(); }""", [sel, frac])
    page.evaluate("() => { if (window.ScrollTrigger) { ScrollTrigger.update(); ScrollTrigger.refresh(); } if (window.__viz) window.__viz.force(1); document.querySelectorAll('.reveal').forEach(e => e.classList.add('in')); }")

with sync_playwright() as pw:
    browser = pw.chromium.launch(args=["--use-gl=angle", "--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"])
    for page_name in PAGES:
        stem = page_name.replace(".html", "").replace("/", "_"); pr = report["pages"].setdefault(page_name, {})
        for theme in ("dark", "light"):
            if A.snap:
                ctx, page, errors, _ = open_page(browser, page_name, theme, {"width": 1280, "height": 860})
                names = page.eval_on_selector_all("[data-viz]", "els => els.map(e => e.dataset.viz)")
                for name in names:
                    sel = f'[data-viz="{name}"]'; scroll_to(page, sel); page.evaluate("() => window.__viz && window.__viz.force(1)"); time.sleep(0.9)
                    out = DOCS / "assets" / f"viz_{name}_{theme}.png"; page.query_selector(sel).screenshot(path=str(out), timeout=180000); print("snap", out.name)
                if errors: print("console errors during snap", theme, errors)
                ctx.close()
            if A.qa:
                qd = DOCS / "qa" / theme; qd.mkdir(parents=True, exist_ok=True)
                ctx, page, errors, failed = open_page(browser, page_name, theme, {"width": 1440, "height": 900})
                page.screenshot(path=str(qd / f"{stem}-hero.png"), timeout=180000)
                names = page.eval_on_selector_all("[data-viz]", "els => els.map(e => e.dataset.viz)")[:2]
                for i, name in enumerate(names, 1): scroll_to(page, f'[data-viz="{name}"]'); time.sleep(2.6); page.screenshot(path=str(qd / f"{stem}-scene{i}-{name}.png"), timeout=180000)
                imgs = page.eval_on_selector_all("img", "els => els.filter(i => i.src && !i.src.startsWith('data:')).map(i => ({ ok: i.complete && i.naturalWidth > 0, src: i.src.split('/').pop() }))")
                nbytes = page.evaluate("() => performance.getEntriesByType('resource').reduce((a, r) => a + (r.transferSize || r.encodedBodySize || 0), 0) + (performance.getEntriesByType('navigation')[0]?.transferSize || 0)")
                pr[theme] = {"errors": errors, "failed": failed, "imagesLoaded": f"{sum(1 for i in imgs if i['ok'])}/{len(imgs)}", "transferBytes": nbytes, "scenes": names}; ctx.close()
                mctx, mpage, merr, _ = open_page(browser, page_name, theme, {"width": 375, "height": 812}); mpage.screenshot(path=str(qd / f"{stem}-mobile-375-hero.png"), timeout=180000); pr[theme]["mobileErrors"] = merr; mctx.close()
                print("qa", theme, json.dumps(pr[theme]))
        if A.hero_png:
            ctx, page, _, _ = open_page(browser, page_name, "dark", {"width": 1600, "height": 900}); time.sleep(2.5)
            page.evaluate("() => document.querySelectorAll('nav, .hero .copy, .hero .readout, .hero .shade, .scrollhint').forEach(e => e.style.visibility = 'hidden')"); time.sleep(.5)
            page.query_selector("#heroCanvas").screenshot(path=str(DOCS / "assets" / "hero.png"), timeout=180000); print("hero.png"); ctx.close()
        if A.hero_json and page_name == PAGES[0]:
            ctx, page, _, _ = open_page(browser, page_name, "dark", {"width": 1280, "height": 800}); page.wait_for_function("() => window.__hero")
            data = page.evaluate("() => window.__hero.PROJECTS.map((p, i) => ({ key: p.key, name: p.name, points: window.__hero.project(i, 300) }))")
            Path(A.hero_json).write_text(json.dumps(data)); print("hero json", A.hero_json, [len(d["points"]) for d in data]); ctx.close()
        if A.lighthouse:
            out = DOCS / "qa" / ("lighthouse" if stem == "index" else f"lighthouse.{stem}")
            # Playwright's full Chromium (not the headless shell); WebGL is disabled for the audit so it measures the page's static path,
            # as the sibling pages' reports do: on macOS the software-GL GPU process alone delays first paint by about two seconds
            chrome = next(iter([p for p in sorted((Path.home() / "Library/Caches/ms-playwright").glob("chromium-*/chrome-mac*/*.app/Contents/MacOS/*")) if "Helper" not in str(p) and "Framework" not in str(p)] + sorted(Path("/opt/pw-browsers").glob("chromium-*/chrome-linux*/chrome"))), None)
            env = {**os.environ, "CHROME_PATH": os.environ.get("CHROME_PATH", str(chrome) if chrome else "")}
            cmd = ["lighthouse", BASE + page_name, "--preset=desktop", "--quiet", "--chrome-flags=--headless=new --no-sandbox --disable-gpu", "--output=html", "--output=json", f"--output-path={out}", "--only-categories=performance,accessibility,best-practices,seo", "--max-wait-for-load=30000"]
            subprocess.run(cmd, check=True, env=env, timeout=600)
            j = json.loads(Path(str(out) + ".report.json").read_text()); c, a = j["categories"], j["audits"]
            pr["lighthouse"] = {"performance": round(c["performance"]["score"] * 100), "accessibility": round(c["accessibility"]["score"] * 100), "bestPractices": round(c["best-practices"]["score"] * 100), "seo": round(c["seo"]["score"] * 100),
                                "fcp": a["first-contentful-paint"]["displayValue"], "lcp": a["largest-contentful-paint"]["displayValue"], "tbt": a["total-blocking-time"]["displayValue"], "cls": a["cumulative-layout-shift"]["displayValue"], "si": a["speed-index"]["displayValue"], "bytes": a["total-byte-weight"]["displayValue"]}
            print("lighthouse", json.dumps(pr["lighthouse"]))
    browser.close()
(DOCS / "qa").mkdir(exist_ok=True); qa_json.write_text(json.dumps(report, indent=1)); srv.shutdown()
