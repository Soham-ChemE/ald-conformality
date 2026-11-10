# QA for the live page

Served from docs/ by docs/tools/page_qa.py (Playwright for Python) and loaded in headless Chromium. The screenshots use software WebGL with `window.__forceGL` so the live scenes render (frame timing is therefore not representative); the Lighthouse run has WebGL disabled, so it measures the page's static path, as the sibling pages' reports do. Screenshots in both themes: the hero, the first two body scenes, the 375 px hero.

## index.html

Lighthouse (desktop preset): performance 100, accessibility 100, best practices 100, SEO 100

| Metric | Value |
|---|---|
| First Contentful Paint | 0.2 s |
| Largest Contentful Paint | 0.2 s |
| Speed Index | 0.4 s |
| Total Blocking Time | 0 ms |
| Cumulative Layout Shift | 0 |
| Total byte weight | Total size was 533 KiB |

| Theme | Console errors | Failed requests | Transfer (bytes) | Body scenes captured |
|---|---|---|---|---|
| dark | 0 | 0 | 736,561 | stack, map |
| light | 0 | 0 | 736,561 | stack, map |

Files: `light/` and `dark/` hold `index-hero.png`, `index-scene1-*.png`, `index-scene2-*.png` and `index-mobile-375-hero.png`; `lighthouse.report.html` and `.json` are the Lighthouse run; `qa.json` is the machine-readable summary.
