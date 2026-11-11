"""Tests of coupling/device_geometries.py: the GAA nanosheet and 3D-NAND word-line cell builders produce the intended
volumes, periodicity and face regions, and a tiny-budget voxel run completes (skipped when libaldmc is not built)."""
import numpy as np
import pytest

import device_geometries as dg
from aldconf import voxel


def test_gaa_stack_volumes_and_periodicity():
    h = 1.25e-9
    g = dg.gaa_stack(n_sheets=3, t_s=5e-9, g_s=10e-9, W=40e-9, Lg=12e-9, gap_x=20e-9, H_top=40e-9, g_bot=10e-9, h=h)
    m = g.meta
    assert g.periodic == (True, False)
    assert g.shape == (48, 12, 69)                       # (40 + 20) / 1.25, 12.5 / 1.25 + 2 pad, 32 + 3*4 + 2*8 + 8 + 1
    occ = g.occupancy
    # trench: gas over the full x range, inside the spacer walls, down to the floor; the spacers and floor are solid
    assert occ[:, 0, :].all() and occ[:, -1, :].all() and occ[:, :, -1].all()
    assert (occ[:, 1:-1, 0] == 0).all()                  # the mouth plane is open across the whole period
    # three sheets of 4 voxels, 32 voxels wide, each spanning the trench length
    x0, x1 = m["sheet_x"]
    for (ka, kb) in m["sheet_k"]:
        assert kb - ka == 4
        assert occ[x0:x1, 1:-1, ka:kb].all()
        assert (occ[:x0, 1:-1, ka:kb] == 0).all() and (occ[x1:, 1:-1, ka:kb] == 0).all()
    # slots of 8 voxels between the sheets and under the lowest one
    for (ka, kb) in m["gap_k"]:
        assert kb - ka == 8 and (occ[x0:x1, 1:-1, ka:kb] == 0).all()
    # gas volume: trench minus the sheets
    gas = int((occ == 0).sum())
    assert gas == 48 * 10 * 68 - 3 * 32 * 10 * 4
    # the smooth area formula counts every exposed voxel face exactly once (axis-aligned geometry)
    faces = 2 * 48 * 68 - 2 * 3 * 32 * 4 + 48 * 10 + 3 * (2 * 32 * 10 + 2 * 4 * 10)
    assert abs(dg.gaa_smooth_area(m) / (faces * h * h) - 1) < 1e-9
    assert abs(m["pitch_x"] - 60e-9) < 1e-15 and abs(m["Lg_actual"] - 12.5e-9) < 1e-15


def test_nand_cell_recesses_and_pillars():
    h = 5e-9
    g = dg.nand_cell(w_s=100e-9, L_s=1e-6, g_r=29e-9, L_r=365e-9, p_r=51e-9, d_h=100e-9, pitch_h=200e-9, rows=2, h=h)
    m = g.meta
    assert g.periodic == (True, False)
    assert g.shape == (40, 20 + 2 * 73 + 2, 201)
    assert np.isclose(m["g_r_actual"], 30e-9) and np.isclose(m["p_r_actual"], 50e-9)
    assert m["n_recesses"] == 20                         # oxide first: 20 whole recesses in 1 um
    occ = g.occupancy
    ya, yb = m["y_slit"]
    # the slit is open from the mouth to the floor except where nothing else: the slit column holds no pillar
    assert (occ[:, ya:yb, :200] == 0).all() and occ[:, ya:yb, 200].all()
    # pillars: four discs of round(pi r^2) voxels per slice, solid through the whole depth, inside the recess bands
    disc = int(round(np.pi * 10 ** 2))
    assert m["pillar_voxels_per_slice"] == 4 * disc
    ka, kb = m["recess_k"][0]
    plan = occ[:, :, (ka + kb) // 2]
    solid_in_recess = plan[:, 1:ya].sum() + plan[:, yb:-1].sum()
    assert solid_in_recess == 4 * disc                   # inside a recess the only solid is the pillars
    kox = m["recess_k"][0][0] - 1                        # an oxide slice is fully solid outside the slit
    assert occ[:, :ya, kox].all() and occ[:, yb:, kox].all()
    # the pillar at x = 0 wraps around the periodic boundary: both edge columns carry solid in the recess plane
    assert plan[0, 1:ya].sum() > 0 and plan[-1, 1:ya].sum() > 0
    # smooth area: the pillars add their walls inside every recess and remove two disc footprints each; the recess
    # is lower than the hole radius, so the cell carries less area than the plain recessed slit
    m0 = dict(m); m0["pillar_xy"] = []; m0["rows"] = 0
    a1, a0 = dg.nand_smooth_area(m), dg.nand_smooth_area(m0)
    per_pillar = m["g_r_actual"] * np.pi * m["d_h"] - 2 * np.pi * m["d_h"] ** 2 / 4
    assert a0 > 0 and np.isclose(a1 - a0, m["n_recesses"] * 4 * per_pillar) and a1 < a0


def test_rule_and_effusion_helpers():
    assert abs(dg.rule_pulse("TMA", 55.9) - 14.38) < 0.05                 # wp2 / wp4 value for the AR-100 wafer
    assert dg.rule_pulse("TiCl4", 55.9) < dg.rule_pulse("TDMAH", 55.9) < dg.rule_pulse("TMA", 55.9)
    assert dg.effusion_limit("TMA", 100.0) > 0


@pytest.mark.skipif(not voxel.is_available(), reason="libaldmc not built: " + voxel.BUILD_COMMAND)
def test_tiny_gaa_run_completes():
    g = dg.build_gaa("gaa_ibm2017")
    dg.THREADS = 2
    rec = dg.run_cell(g, "TMA", "constant", dt=0.02, M=2000, t_end=0.2, regions_fn=dg.gaa_regions, label="tiny",
                      area_correction=False, verbose_every=0, stop_level=2.0)
    assert len(rec["t_hist"]) == 10
    assert set(rec["hist"]) >= {"sheet0_top", "gap0_mid", "gap_bottom_mid", "floor_under_stack"}
    # the plain-face inventory of the axis-aligned cell equals the smooth area
    assert abs(rec["site_area"] / dg.gaa_smooth_area(g.meta) - 1) < 1e-9
    # the top of the top sheet leads the slot centres; every particle is bookkept
    assert rec["hist"]["sheet0_top"][-1] > rec["hist"]["gap_bottom_mid"][-1] > 0
    c = rec["counts"]
    assert c["escaped"] + c["stuck"] == rec["particles"] and c["left"] == 0
    cell, wafer = dg.gaa_wall_ratios(g, rec["site_area"])
    assert cell > wafer > 1


@pytest.mark.skipif(not voxel.is_available(), reason="libaldmc not built: " + voxel.BUILD_COMMAND)
def test_tiny_nand_run_completes():
    g = dg.nand_cell(w_s=100e-9, L_s=0.5e-6, g_r=29e-9, L_r=200e-9, p_r=51e-9, d_h=100e-9, pitch_h=200e-9, rows=1, h=10e-9)
    dg.THREADS = 2
    rec = dg.run_cell(g, "TMA", "constant", dt=0.05, M=2000, t_end=0.15, regions_fn=dg.nand_regions, label="tiny",
                      normal_radius=2, verbose_every=0, stop_level=2.0)
    assert len(rec["t_hist"]) == 3
    assert "recess_deep_deephalf" in rec["hist"] and len(rec["hist_min_list"]["recess_all_deephalf"]) == 3
    c = rec["counts"]
    assert c["escaped"] + c["stuck"] == rec["particles"] and c["left"] == 0
