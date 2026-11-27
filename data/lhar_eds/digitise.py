"""Digitise raster SEM-EDS line scans (Yim et al., PCCP 22, 23107 (2020)).
Run: .venv/bin/python -I data/lhar_eds/digitise.py   (from repo root)
No smoothing, normalisation or shifting. One point per pixel column.
"""
import glob, json, os
import numpy as np, pymupdf
from PIL import Image

OUT = os.path.dirname(os.path.abspath(__file__))
INBOX = os.path.join(OUT, "..", "papers_inbox")
ESI = os.path.join(INBOX, "d0cp03358h1_suppl.pdf")
MAIN = glob.glob(os.path.join(INBOX, "Saturation profile based conformality*.pdf"))[0]


def load(pdf, xref):
    d = pymupdf.open(pdf)
    pix = pymupdf.Pixmap(d, xref)
    if pix.n - pix.alpha > 3:
        pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    a = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
    return a.astype(float)


def line_centres(profile, thr):
    """Weighted centres of runs where profile > thr."""
    out, i, n = [], 0, len(profile)
    while i < n:
        if profile[i] > thr:
            j = i
            while j < n and profile[j] > thr:
                j += 1
            w = profile[i:j]
            out.append(float((np.arange(i, j) * w).sum() / w.sum()))
            i = j
        else:
            i += 1
    return out


def calib(px, vals):
    px, vals = np.asarray(px), np.asarray(vals, float)
    A = np.vstack([px, np.ones_like(px)]).T
    (m, c), *_ = np.linalg.lstsq(A, vals, rcond=None)
    res = vals - (m * px + c)
    return m, c, res


def classify(a, palette, maxdist):
    names = list(palette)
    P = np.array([palette[k] for k in names], float)
    d = np.sqrt(((a[:, :, None, :] - P[None, None]) ** 2).sum(-1))
    idx = d.argmin(-1)
    ok = d.min(-1) < maxdist
    return {k: (idx == i) & ok for i, k in enumerate(names)}


def trace(mask, others, x0, x1, y0, y1):
    """Per column: centre of the vertical run of mask pixels (largest run; if several, the one continuous with the previous column).
    Returns cols, ycentre, flag (1 if another trace touches the run), runlen."""
    rows = []
    prev = None
    for x in range(x0, x1 + 1):
        col = mask[y0:y1 + 1, x]
        ys = np.flatnonzero(col)
        if ys.size == 0:
            continue
        # split into runs (gap tolerance 1 px)
        brk = np.flatnonzero(np.diff(ys) > 2)
        runs = np.split(ys, brk + 1)
        r = max(runs, key=len)
        if len(runs) > 1 and prev is not None:
            # continuity: of the runs at least 2 px long, take the one nearest the previous column's centre
            cand = [q for q in runs if len(q) >= 2] or runs
            r = min(cand, key=lambda q: abs((q[0] + q[-1]) / 2.0 + y0 - prev))
        prev = (r[0] + r[-1]) / 2.0 + y0
        top, bot = r[0] + y0, r[-1] + y0
        o = others[:, x]
        # a real occluding stroke is several px thick; isolated edge-blend pixels are ignored
        above = o[max(top - 4, 0):top].sum() >= 3
        below = o[bot + 1:bot + 5].sum() >= 3
        rows.append((x, top, bot, int(above), int(below), len(runs)))
    return np.array(rows)


def finish(rows, mx, cx, my, cy, name, header):
    x, top, bot, ab, be = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4]
    run = bot - top + 1
    clean = (ab == 0) & (be == 0)
    t = float(np.median(run[clean & (run <= np.percentile(run[clean], 60))])) if clean.any() else 5.0
    yc = (top + bot) / 2.0
    # partially hidden flat stroke: use the visible edge and half the stroke thickness
    part = (~clean) & (run < t - 0.5)
    yc = np.where(part & (be == 1) & (ab == 0), top + (t - 1) / 2.0, yc)
    yc = np.where(part & (ab == 1) & (be == 0), bot - (t - 1) / 2.0, yc)
    flag = (~clean).astype(int)
    xu = mx * x + cx
    yv = my * yc + cy
    p = os.path.join(OUT, name + ".csv")
    with open(p, "w") as f:
        f.write("x_um,intensity,overlap_flag\n")
        for a_, b_, c_ in zip(xu, yv, flag):
            f.write(f"{a_:.2f},{b_:.1f},{c_}\n")
    return dict(id=name, n=int(len(xu)), xmin=float(xu.min()), xmax=float(xu.max()),
                n_flag=int(flag.sum()), stroke_px=t, x=xu, y=yv, flag=flag)


def overlay(a, results, mx, cx, my, cy, path, scale=4):
    im = Image.fromarray(a.astype(np.uint8)).resize((a.shape[1] * scale, a.shape[0] * scale), Image.NEAREST)
    px = im.load()
    for r in results:
        for xu, yv in zip(r["x"], r["y"]):
            X = int(round(((xu - cx) / mx + 0.5) * scale)); Y = int(round(((yv - cy) / my + 0.5) * scale))
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if 0 <= X + dx < im.size[0] and 0 <= Y + dy < im.size[1]:
                        px[X + dx, Y + dy] = (0, 0, 0)
    im.save(path)


def stats(r, lo, hi, clean_only=False):
    m = (r["x"] >= lo) & (r["x"] <= hi)
    if clean_only:
        m &= r["flag"] == 0
    return (float(np.mean(r["y"][m])), float(np.std(r["y"][m])), int(m.sum())) if m.any() else None


summary = {}

# ---------------- ESI Fig. S8 (page 8, xref 105, 726x438) ----------------
a = load(ESI, 105)
neutral = (abs(a[:, :, 0] - a[:, :, 1]) < 12) & (abs(a[:, :, 1] - a[:, :, 2]) < 12) & (a[:, :, 0] < 235) & (a[:, :, 0] > 140)
dark = np.where(neutral, 255 - a[:, :, 0], 0)
rc = line_centres(dark[:, 136:680].sum(1), 0.15 * dark[:, 136:680].sum(1).max())
cc = line_centres(dark[20:270].sum(0), 0.25 * dark[20:270].sum(0).max())
rc = [r for r in rc if 15 < r < 290]
cc = [c for c in cc if 120 < c < 700]
print("S8 grid rows", np.round(rc, 2)); print("S8 grid cols", np.round(cc, 2))
assert len(rc) == 6 and len(cc) == 11
my, cy, ry = calib(rc, [2500, 2000, 1500, 1000, 500, 0])
mx, cx, rx = calib(cc, np.arange(0, 501, 50))
print("S8 x: um/px", mx, "resid", np.round(rx, 2), " y: a.u./px", my, "resid", np.round(ry, 1))
pal = dict(white=(255, 255, 255), grid=(217, 217, 217), axis=(170, 170, 170), black=(60, 60, 60),
           yellow=(255, 192, 0), gold=(191, 143, 0), orange=(237, 125, 49))
m = classify(a, pal, 70)
x0, x1 = int(np.ceil(cc[0])), int(np.floor(cc[-1])); y0, y1 = int(rc[0]) - 2, int(rc[-1]) + 6
res = []
# draw order (legend order, confirmed at crossings): yellow, gold, orange (top).
# Only later-drawn series can hide a trace.
for col, sid, occ in (("yellow", "figS8_Al_V0057_250cyc", ("gold", "orange")), ("gold", "figS8_Al_V0069_500cyc_as_labelled", ("orange",)), ("orange", "figS8_Al_V0031_1000cyc", ())):
    oth = np.zeros_like(m[col])
    for k in occ:
        oth |= m[k]
    rows = trace(m[col], oth, x0, x1, y0, y1)
    res.append(finish(rows, mx, cx, my, cy, sid, ""))
overlay(a, res, mx, cx, my, cy, os.path.join(OUT, "check_figS8_overlay.png"))
Image.fromarray(a.astype(np.uint8)).save(os.path.join(OUT, "source_figS8_xref105.png"))
for r in res:
    r["bg_left"] = stats(r, 0, 40); r["bg_right"] = stats(r, 460, 500)
    r["bg_left_clean"] = stats(r, 0, 40, True); r["bg_right_clean"] = stats(r, 460, 500, True)
    r["open_area_225_275"] = stats(r, 225, 275)
    r["plateau_L_150_200"] = stats(r, 150, 200); r["plateau_R_305_355"] = stats(r, 305, 355)
    summary[r["id"]] = {k: v for k, v in r.items() if k not in ("x", "y", "flag")}
summary["figS8_calib"] = dict(um_per_px=mx, au_per_px=my, x_resid_um=rx.tolist(), y_resid_au=ry.tolist(), grid_cols_px=cc, grid_rows_px=rc)

# ---------------- Main Fig. 5 (page 6, xref 250, 807x1181), panel (b) ----------------
a = load(MAIN, 250)
OFF = 800
b = a[OFF:]
neutral = (abs(b[:, :, 0] - b[:, :, 1]) < 12) & (abs(b[:, :, 1] - b[:, :, 2]) < 12) & (b[:, :, 0] < 240) & (b[:, :, 0] > 185)
dark = np.where(neutral, 255 - b[:, :, 0], 0)
rp = dark[:, 170:788].sum(1); cp = dark[12:285].sum(0)
rc = [r + OFF for r in line_centres(rp, 0.4 * rp.max()) if r < 300]
cc = [c for c in line_centres(cp, 0.4 * cp.max()) if 150 < c < 800]
print("F5 grid rows", np.round(rc, 2)); print("F5 grid cols", np.round(cc, 2))
assert len(rc) == 6 and len(cc) == 11
my, cy, ry = calib(rc, [5000, 4000, 3000, 2000, 1000, 0])
mx, cx, rx = calib(cc, np.arange(-250, 251, 50))
print("F5 x: um/px", mx, "resid", np.round(rx, 2), " y: a.u./px", my, "resid", np.round(ry, 1))
pal = dict(white=(255, 255, 255), grid=(217, 217, 217), black=(60, 60, 60),
           blue=(91, 155, 213), orange=(237, 125, 49), grey=(165, 165, 165))
m = classify(a, pal, 60)
x0, x1 = int(np.ceil(cc[0])), int(np.floor(cc[-1])); y0, y1 = int(rc[0]) - 2, int(rc[-1]) + 6
# the x-axis line is mid-grey and would be mistaken for the Si trace; Si never drops below ~1300 a.u.
m["grey"][int(rc[4]):, :] = False
res = []
# draw order: O (blue), Al (orange), Si (grey, top)
for col, sid, occ in (("orange", "fig5b_Al_sample11_1000cyc", ("grey",)), ("blue", "fig5b_O_sample11_1000cyc", ("orange", "grey")), ("grey", "fig5b_Si_sample11_1000cyc", ())):
    oth = np.zeros_like(m[col])
    for k in occ:
        oth |= m[k]
    rows = trace(m[col], oth, x0, x1, y0, y1)
    res.append(finish(rows, mx, cx, my, cy, sid, ""))
overlay(a[OFF - 40:], [dict(x=r["x"], y=r["y"]) for r in res], mx, cx, my, cy + my * (OFF - 40), os.path.join(OUT, "check_fig5b_overlay.png"))
Image.fromarray(a.astype(np.uint8)).save(os.path.join(OUT, "source_fig5_xref250.png"))
for r in res:
    r["bg_left_-250_-200"] = stats(r, -250, -200); r["bg_right_200_250"] = stats(r, 200, 250)
    r["open_area_-30_30"] = stats(r, -30, 30)
    r["plateau_L_-110_-50"] = stats(r, -110, -50); r["plateau_R_50_110"] = stats(r, 50, 110)
    summary[r["id"]] = {k: v for k, v in r.items() if k not in ("x", "y", "flag")}
summary["fig5b_calib"] = dict(um_per_px=mx, au_per_px=my, x_resid_um=rx.tolist(), y_resid_au=ry.tolist(), grid_cols_px=cc, grid_rows_px=rc)
json.dump(summary, open(os.path.join(OUT, "digitise_summary.json"), "w"), indent=1, default=float)
for k, v in summary.items():
    if "calib" not in k:
        print(k, {a_: (np.round(b_, 1).tolist() if isinstance(b_, tuple) else b_) for a_, b_ in v.items() if a_ != "id"})
