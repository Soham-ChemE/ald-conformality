"""Full 3D chamber with a one-sided pumping port (reactor/chamber3d): per-face saturation times of the hole bottoms
and the flat surface on the wafer patch, azimuthal statistics at fixed radii, the worst face on the wafer against the
axisymmetric wedge, and the gas-side asymmetry (steady radial velocity and surface TMA density, port side against
far side). Writes figures/chamber3d_t90.csv, figures/chamber3d_t90_map.png, figures/chamber3d_azimuthal.png and
figures/chamber3d_summary.txt.
  .venv/bin/python coupling/postprocess_chamber3d.py [--case reactor/chamber3d/pulse] [--no-cache]

Reader: the wafer patch is one patch of 6144 faces whose order is the same in every field file. Face centres come
from the boundaryField of C written by `postProcess -func writeCellCentres -time 0` (boundary values of C are the
patch face centres); the per-face thetaBottom, thetaFlat and value lists are read from the wafer block of each
time directory's T file (the thetaHole list, 200 entries per face, is not read). Radius r = sqrt(y^2 + z^2),
azimuth phi = atan2(z, y) with phi = 0 the port direction (+y), as in the BC and the mesh.
t90 is the 0.9 crossing of the coverage, linearly interpolated between write times (0.1 s apart); the first write
time at or above 0.9 is kept as well (the wedge numbers in earlier parts are first-write-time values at 0.05 s).
"""
import re, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "coupling")]
from patterned_lib import load_patterned

OUT = ROOT / "figures"
CASE = ROOT / "reactor/chamber3d/pulse"
WEDGE = ROOT / "reactor/showerhead_pulse_patterned_AR100_long"
RADII_MM = (50.0, 100.0, 140.0)
BAND_MM = 2.5                       # faces within +- band of the nominal radius (cells are 2-5 mm here)
NBINS = 24                          # azimuthal bins of 15 degrees
PORT_HALF_DEG = 30.0
T_SNAP = (0.3, 8.0, 14.0)           # times at which the surface TMA density / coverage are compared by azimuth
LEVEL = 0.9                         # coverage level defining the saturation time (--level x for testing on short runs)
if "--case" in sys.argv: CASE = Path(sys.argv[sys.argv.index("--case") + 1]).resolve()
if "--level" in sys.argv: LEVEL = float(sys.argv[sys.argv.index("--level") + 1])


# ------------------------------------------------------------------ readers
def _patch_block(txt, name):
    m = re.search(r"\n\s+%s\s*\{(.*?)\n\s+\}" % name, txt, re.S)
    return m.group(1) if m else ""


def _scalar_list(blk, key, n):
    """Per-face list of `key` in a patch block; a uniform entry is broadcast, a missing one (time 0) is zero."""
    m = re.search(r"\b%s\s+nonuniform\s+List<scalar>\s*\n(\d+)\s*\n\(\n(.*?)\n\)" % key, blk, re.S)
    if m: return np.array(m.group(2).split(), float)
    m = re.search(r"\b%s\s+uniform\s+([0-9.eE+-]+)" % key, blk)
    return np.full(n, float(m.group(1)) if m else 0.0)


def _vector_list(txt_or_blk, key="value"):
    m = re.search(r"\b%s\s+nonuniform\s+List<vector>\s*\n(\d+)\s*\n\(\n(.*?)\n\)" % key, txt_or_blk, re.S)
    n = int(m.group(1))
    a = np.array(m.group(2).replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)
    return a


def wafer_face_centres(case):
    txt = (case / "0/C").read_text()
    return _vector_list(_patch_block(txt, "wafer"))


def cell_centres(case):
    return _vector_list((case / "0/C").read_text(), "internalField")


def load_wafer_fields(case, n_faces, use_cache=True):
    """t (ntime), fields dict thetaBottom / thetaFlat / value (ntime x nFaces) from the wafer block of each T file."""
    cache = case / "wafer_fields.npz"
    times = sorted((float(d.name), d) for d in case.iterdir()
                   if re.fullmatch(r"[0-9.]+(e-?[0-9]+)?", d.name) and (d / "T").exists())
    t = np.array([tt for tt, _ in times])
    if use_cache and cache.exists():
        z = np.load(cache)
        if len(z["t"]) == len(t) and np.allclose(z["t"], t):
            return z["t"], {k: z[k] for k in ("thetaBottom", "thetaFlat", "value")}
    keys = ("thetaBottom", "thetaFlat", "value")
    out = {k: [] for k in keys}
    t0 = time.time()
    for i, (_, d) in enumerate(times):
        blk = _patch_block((d / "T").read_text(), "wafer")
        for k in keys: out[k].append(_scalar_list(blk, k, n_faces))
        if i % 20 == 0: print(f"  read {i + 1}/{len(times)} time directories ({time.time() - t0:.0f} s)", flush=True)
    out = {k: np.array(v) for k, v in out.items()}
    np.savez_compressed(cache, t=t, **out)
    return t, out


def t90(theta, t, level=None):
    """LEVEL (0.9) crossing per face: (interpolated time, first write time at or above it); nan where never reached."""
    level = LEVEL if level is None else level
    n = theta.shape[1]
    ti, tw = np.full(n, np.nan), np.full(n, np.nan)
    hit = theta.max(axis=0) >= level
    k = np.argmax(theta >= level, axis=0)
    for j in np.where(hit)[0]:
        kk = k[j]; tw[j] = t[kk]
        if kk == 0: ti[j] = t[0]; continue
        a, b = theta[kk - 1, j], theta[kk, j]
        ti[j] = t[kk - 1] + (level - a) / (b - a) * (t[kk] - t[kk - 1])
    return ti, tw


def band_stats(vals, r_mm, phi_deg, r0):
    sel = np.abs(r_mm - r0) <= BAND_MM
    v, p = vals[sel], phi_deg[sel]
    edges = np.linspace(-180, 180, NBINS + 1)
    idx = np.clip(np.digitize(p, edges) - 1, 0, NBINS - 1)
    binmean = np.array([np.nanmean(v[idx == b]) if np.any(idx == b) else np.nan for b in range(NBINS)])
    # value at exactly r0 per bin: linear fit in r of the bin's faces (removes the radial gradient inside the band)
    rr = r_mm[sel]; at_r0 = np.full(NBINS, np.nan)
    for b in range(NBINS):
        m = (idx == b) & np.isfinite(v)
        if m.sum() >= 3 and np.ptp(rr[m]) > 0.1: at_r0[b] = np.polyval(np.polyfit(rr[m] - r0, v[m], 1), 0.0)
        elif m.any(): at_r0[b] = np.nanmean(v[m])
    port = np.abs(p) <= PORT_HALF_DEG; far = np.abs(p) >= 180 - PORT_HALF_DEG
    mean = np.nanmean(v)
    return dict(n=int(sel.sum()), mean=mean, min=np.nanmin(v), max=np.nanmax(v),
                p2p_pct=100 * (np.nanmax(v) - np.nanmin(v)) / mean,
                bin_centres=0.5 * (edges[1:] + edges[:-1]), bin_mean=binmean,
                bin_p2p_pct=100 * (np.nanmax(binmean) - np.nanmin(binmean)) / mean,
                at_r0=at_r0, at_r0_p2p_pct=100 * (np.nanmax(at_r0) - np.nanmin(at_r0)) / np.nanmean(at_r0),
                at_r0_port=np.nanmean(at_r0[np.abs(0.5 * (edges[1:] + edges[:-1])) <= PORT_HALF_DEG]),
                at_r0_far=np.nanmean(at_r0[np.abs(0.5 * (edges[1:] + edges[:-1])) >= 180 - PORT_HALF_DEG]),
                port=np.nanmean(v[port]), far=np.nanmean(v[far]),
                port_minus_far_pct=100 * (np.nanmean(v[port]) - np.nanmean(v[far])) / mean,
                sel=sel, idx=idx)


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    lines = []
    def say(s=""):
        print(s, flush=True); lines.append(s)

    Cf = wafer_face_centres(CASE)
    y_mm, z_mm = 1e3 * Cf[:, 1], 1e3 * Cf[:, 2]
    r_mm = np.hypot(y_mm, z_mm); phi = np.degrees(np.arctan2(z_mm, y_mm))
    nF = Cf.shape[0]
    t, f = load_wafer_fields(CASE, nF, use_cache="--no-cache" not in sys.argv)
    assert f["thetaBottom"].shape[1] == nF, (f["thetaBottom"].shape, nF)
    pulse = float(re.search(r"\(\(0 1\) \(([0-9.]+) 1\)", (CASE / "0/T").read_text()).group(1))
    say(f"case {CASE.relative_to(ROOT)}: {nF} wafer faces, {len(t)} write times 0..{t[-1]:g} s (every {np.diff(t).mean():.3g} s), "
        f"pulse {pulse:g} s; face radii {r_mm.min():.2f}..{r_mm.max():.2f} mm")

    tb, tbw = t90(f["thetaBottom"], t)
    tf_, tfw = t90(f["thetaFlat"], t)

    # wedge reference (same reader logic on the 150 radial faces of the axisymmetric case)
    wt, wf, wpar, wpulse = load_patterned(WEDGE)
    w_r = np.arange(150) + 0.5
    wtb, wtbw = t90(wf["thetaBottom"], wt); wtf, wtfw = t90(wf["thetaFlat"], wt)
    def wedge_at(arr, r0): return float(np.interp(r0, w_r, arr))
    say(f"wedge reference ({WEDGE.name}, writes every {np.diff(wt).mean():.3g} s, pulse {wpulse:g} s): hole-bottom t90 centre "
        f"{wtb[0]:.3f} s interpolated / {wtbw[0]:.2f} s first write; worst {np.nanmax(wtb):.3f} s at r = {w_r[np.nanargmax(wtb)]:.1f} mm "
        f"(first write {np.nanmax(wtbw):.2f} s); flat centre {wtf[0]:.3f} / {wtfw[0]:.2f} s, flat edge {wtf[149]:.3f} s")

    # ---- whole-wafer numbers
    jc = np.argmin(r_mm); jw = np.nanargmax(tb)
    say(f"3D: hole-bottom t90 at the innermost face (r = {r_mm[jc]:.1f} mm, phi = {phi[jc]:.0f} deg): {tb[jc]:.3f} s interpolated, "
        f"{tbw[jc]:.1f} s first write; wedge centre {wtb[0]:.3f} s ({100*(tb[jc]/wtb[0]-1):+.2f} %)")
    say(f"3D: worst hole-bottom t90 on the wafer {tb[jw]:.3f} s (first write {tbw[jw]:.1f} s) at r = {r_mm[jw]:.1f} mm, "
        f"phi = {phi[jw]:.0f} deg; wedge worst {np.nanmax(wtb):.3f} s ({100*(tb[jw]/np.nanmax(wtb)-1):+.2f} %); "
        f"wafer min {np.nanmin(tb):.3f} s at r = {r_mm[np.nanargmin(tb)]:.1f} mm; "
        f"wafer-wide spread {100*(np.nanmax(tb)-np.nanmin(tb))/np.nanmean(tb):.2f} % of the mean {np.nanmean(tb):.3f} s")
    say(f"3D: flat-surface t90 innermost face {tf_[jc]:.3f} s (wedge {wtf[0]:.3f}), wafer range {np.nanmin(tf_):.3f}..{np.nanmax(tf_):.3f} s; "
        f"faces never reaching 0.9: bottoms {int(np.isnan(tb).sum())}, flat {int(np.isnan(tf_).sum())}")
    # deviation from the wedge radial profile, face by face
    dev_pct = 100 * (tb / np.interp(r_mm, w_r, wtb) - 1)
    say(f"3D minus wedge at the same radius, hole-bottom t90: mean {np.nanmean(dev_pct):+.2f} %, range {np.nanmin(dev_pct):+.2f}..{np.nanmax(dev_pct):+.2f} %")

    # ---- azimuthal statistics at the three radii
    say("\nAzimuthal statistics of the hole-bottom t90 (faces within +-%.1f mm of the radius; azimuth bins of %d deg; port at phi = 0):"
        % (BAND_MM, 360 // NBINS))
    say("| r (mm) | faces | mean (s) | min (s) | max (s) | peak-to-peak (% of mean) | 15-deg bin means p2p (%) | detrended to r exactly: p2p (%) | port side - far side (%) | wedge t90 at r (s) | 3D mean / wedge |")
    say("|---|---|---|---|---|---|---|---|---|---|---|")
    stats = {}
    for r0 in RADII_MM:
        S = band_stats(tb, r_mm, phi, r0); Sf = band_stats(tf_, r_mm, phi, r0); stats[r0] = (S, Sf)
        wv = wedge_at(wtb, r0)
        say(f"| {r0:.0f} | {S['n']} | {S['mean']:.3f} | {S['min']:.3f} | {S['max']:.3f} | {S['p2p_pct']:.2f} | {S['bin_p2p_pct']:.2f} | {S['at_r0_p2p_pct']:.2f} | "
            f"{S['port_minus_far_pct']:+.2f} | {wv:.3f} | {S['mean']/wv:.4f} |")
        say(f"    detrended profile at r = {r0:.0f} mm: port bins {S['at_r0_port']:.3f} s, far bins {S['at_r0_far']:.3f} s, "
            f"min {np.nanmin(S['at_r0']):.3f} at {S['bin_centres'][np.nanargmin(S['at_r0'])]:.0f} deg, max {np.nanmax(S['at_r0']):.3f} at {S['bin_centres'][np.nanargmax(S['at_r0'])]:.0f} deg")
    say("Flat-surface t90, same bands:")
    say("| r (mm) | mean (s) | min (s) | max (s) | peak-to-peak (%) | port - far (%) | wedge flat t90 at r (s) |")
    say("|---|---|---|---|---|---|---|")
    for r0 in RADII_MM:
        S, Sf = stats[r0]
        say(f"| {r0:.0f} | {Sf['mean']:.3f} | {Sf['min']:.3f} | {Sf['max']:.3f} | {Sf['p2p_pct']:.2f} | {Sf['port_minus_far_pct']:+.2f} | {wedge_at(wtf, r0):.3f} |")
    r_out = r_mm.max()
    S = band_stats(tb, r_mm, phi, r_out - 0.5)
    say(f"Outermost row of wafer faces (r = {r_out:.1f} mm, {S['n']} faces): hole-bottom t90 mean {S['mean']:.3f}, min {S['min']:.3f}, max {S['max']:.3f} s, "
        f"p2p {S['p2p_pct']:.2f} %, port side {S['port']:.3f} s, far side {S['far']:.3f} s ({S['port_minus_far_pct']:+.2f} %); "
        f"wedge outermost face (r = 149.5 mm) {wtb[149]:.3f} s interpolated / {wtbw[149]:.2f} s first write")
    # radial contribution to the band spread, from the wedge profile over the same band
    for r0 in RADII_MM:
        wb = wtb[(w_r >= r0 - BAND_MM) & (w_r <= r0 + BAND_MM)]
        say(f"  wedge radial variation of the bottom t90 across the r = {r0:.0f} +- {BAND_MM} mm band: {100*(wb.max()-wb.min())/wb.mean():.2f} % "
            f"(the part of the 3D band spread that is radial, not azimuthal)")

    # ---- gas-side asymmetry: surface TMA density and coverage by azimuth at fixed times
    say("\nSurface TMA density value/c0 on the wafer faces and coverages at fixed times, r = 140 mm band, port side (|phi| <= 30 deg) vs far side (|phi| >= 150 deg):")
    for ts in T_SNAP:
        k = np.argmin(abs(t - ts))
        for key, lab in (("value", "TMA at the surface / c0"), ("thetaFlat", "flat coverage"), ("thetaBottom", "hole-bottom coverage")):
            S = band_stats(f[key][k], r_mm, phi, 140.0)
            Sc = band_stats(f[key][k], r_mm, phi, 50.0)
            say(f"  t = {t[k]:g} s, {lab}: r 140 mm port {S['port']:.4g}, far {S['far']:.4g}, p2p {S['p2p_pct']:.2f} % of mean {S['mean']:.4g}; "
                f"r 50 mm p2p {Sc['p2p_pct']:.2f} % of mean {Sc['mean']:.4g}")

    # ---- steady carrier flow: radial velocity by azimuth from the frozen U field
    try:
        Cc = cell_centres(CASE); Uc = _vector_list((CASE / "0/U").read_text(), "internalField")
        cy, cz = 1e3 * Cc[:, 1], 1e3 * Cc[:, 2]; cr = np.hypot(cy, cz); cphi = np.degrees(np.arctan2(cz, cy))
        ur = (Uc[:, 1] * cy + Uc[:, 2] * cz) / np.maximum(cr, 1e-9)
        say("\nSteady carrier flow (gap-averaged radial velocity u_r, positive outward, from the frozen 3D U field; axisymmetric value Q/(2 pi r g)):")
        for r0 in (100.0, 140.0, 148.0, 175.0, 199.0):
            sel = np.abs(cr - r0) <= BAND_MM
            if not sel.any(): continue
            p = cphi[sel]; v = ur[sel]
            port = np.abs(p) <= PORT_HALF_DEG; far = np.abs(p) >= 180 - PORT_HALF_DEG
            side = (np.abs(np.abs(p) - 90) <= PORT_HALF_DEG)
            uax = 3.5337e-3 / (2 * np.pi * r0 * 1e-3 * 0.015)
            say(f"  r = {r0:.0f} mm: port side {v[port].mean():+.4f} m/s, 90 deg {v[side].mean():+.4f} m/s, far side {v[far].mean():+.4f} m/s; "
                f"axisymmetric {uax:.4f} m/s")
        sel = cr > 185
        umag = np.linalg.norm(Uc[sel], axis=1); p = cphi[sel]
        say(f"  |U| in the outer 15 mm: port sector mean {umag[np.abs(p) <= PORT_HALF_DEG].mean():.3f} m/s, far sector {umag[np.abs(p) >= 150].mean():.4f} m/s; "
            f"port face mean speed from the flow rate 3.5337e-3 m3/s / 3.1415e-3 m2 = {3.5337e-3/3.1415e-3:.3f} m/s")
        say(f"  gap Peclet number U g / D = {0.05*0.015/0.0084:.3f}; radial Peclet at the wafer edge u_r g / D (axisymmetric) = "
            f"{3.5337e-3/(2*np.pi*0.15*0.015)*0.015/0.0084:.2f}; lateral diffusion time over the wafer radius R^2/D = {0.15**2/0.0084:.2f} s")
    except Exception as e:  # noqa
        say(f"steady-flow analysis skipped: {e}")

    # ---- csv
    k03, k8 = np.argmin(abs(t - 0.3)), np.argmin(abs(t - 8.0))
    df = pd.DataFrame(dict(y_mm=y_mm, z_mm=z_mm, r_mm=r_mm, azimuth_deg=phi,
                           t90_bottom_s=tb, t90_bottom_first_write_s=tbw, t90_flat_s=tf_, t90_flat_first_write_s=tfw,
                           dev_from_wedge_pct=dev_pct,
                           **{f"value_c0_{t[k03]:g}s": f["value"][k03], f"value_c0_{t[k8]:g}s": f["value"][k8]}))
    df.to_csv(OUT / "chamber3d_t90.csv", index=False, float_format="%.6g")

    # ---- figures
    COL = {50.0: "#2a78d6", 100.0: "#eb6834", 140.0: "#3a9d5d"}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.6), subplot_kw=dict(projection="polar"))
    th = np.radians(phi)
    for ax, vals, lab, cmap in ((axes[0], tb, "hole-bottom t90 (s)", "viridis"),
                                (axes[1], dev_pct, "hole-bottom t90, deviation from the wedge at the same radius (%)", "viridis")):
        vmin, vmax = np.nanpercentile(vals, 0.5), np.nanpercentile(vals, 99.5)
        sc = ax.scatter(th, r_mm, c=vals, s=6, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True)
        cb = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.07, shrink=0.85); cb.set_label(lab)
        ax.set_rmax(205); ax.set_rticks([50, 100, 150, 200]); ax.set_rlabel_position(135)
        ax.plot(np.linspace(-np.pi/6, np.pi/6, 50), np.full(50, 200), color="#c0392b", lw=5, solid_capstyle="butt")
        ax.text(0.0, 214, "pump port\n(60 deg)", color="#c0392b", ha="left", va="center", fontsize=9, clip_on=False)
        ax.plot(np.linspace(0, 2*np.pi, 200), np.full(200, 150), color="0.4", lw=0.8, ls="--")
        ax.grid(alpha=0.35)
    axes[0].set_title(f"Hole-bottom t90 over the wafer, 3D chamber with a one-sided pump\nmean {np.nanmean(tb):.2f} s, "
                      f"wafer-wide spread {100*(np.nanmax(tb)-np.nanmin(tb))/np.nanmean(tb):.2f} % (wafer edge dashed, port red)", fontsize=10, pad=16)
    axes[1].set_title(f"Deviation from the axisymmetric wedge at the same radius\nmean {np.nanmean(dev_pct):+.2f} %, "
                      f"range {np.nanmin(dev_pct):+.2f} to {np.nanmax(dev_pct):+.2f} %", fontsize=10, pad=16)
    fig.tight_layout(); fig.savefig(OUT / "chamber3d_t90_map.png", dpi=170); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for ax, key, arr, warr, ttl in ((axes[0], "bottom", tb, wtb, "hole-bottom t90"), (axes[1], "flat", tf_, wtf, "flat-surface t90")):
        ax.axvspan(-PORT_HALF_DEG, PORT_HALF_DEG, color="#c0392b", alpha=0.08, lw=0)
        ax.text(0, 0.02, "pump port", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color="#c0392b", fontsize=8)
        for r0 in RADII_MM:
            S = band_stats(arr, r_mm, phi, r0)
            ax.plot(phi[S["sel"]], arr[S["sel"]], ".", color=COL[r0], ms=3, alpha=0.35)
            ax.plot(S["bin_centres"], S["at_r0"], "-", color=COL[r0], lw=2, label=f"r = {r0:.0f} mm")
            ax.axhline(wedge_at(warr, r0), color=COL[r0], lw=1.2, ls="--")
        ax.set_xlabel("azimuth from the pump port (deg)"); ax.set_ylabel(f"{ttl} (s)"); ax.set_xlim(-180, 180)
        ax.set_xticks(np.arange(-180, 181, 60)); ax.grid(alpha=0.3)
        ax.set_title(f"{ttl} against azimuth, 3D chamber with a one-sided pump\n"
                     "dots: faces within +-2.5 mm; lines: 15 deg bin values at the radius; dashed: wedge", fontsize=9)
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
    fig.tight_layout(); fig.savefig(OUT / "chamber3d_azimuthal.png", dpi=170); plt.close(fig)

    (OUT / "chamber3d_summary.txt").write_text("\n".join(lines) + "\n")
    print(f"\nwrote {OUT / 'chamber3d_t90.csv'}, {OUT / 'chamber3d_t90_map.png'}, {OUT / 'chamber3d_azimuthal.png'}, {OUT / 'chamber3d_summary.txt'}")
