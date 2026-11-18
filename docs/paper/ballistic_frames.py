"""Ballistic (Machball, Yanguas-Gil and Elam) coverage profiles of a long closed parallel-plate trench, saved frame by
frame so that measured LHAR profiles can be fitted without rerunning the solver. The loop is Machball's solve_ideal
with one change: it stops when the half-coverage front passes Z_STOP instead of when the end wall saturates, so the
end wall (AR 800) stays far from the fitted fronts (z50 = 215-270) and its 'trunk' does not enter the fit.
  PYTHONPATH=src .venv/bin/python docs/paper/ballistic_frames.py 0.02 [Nz]
Writes figures/ballistic_frames/s_<s>_Nz<Nz>.npz with exposure (Pa s), z (gap units), coverage frames.
"""
import sys, time
from pathlib import Path
import numpy as np
from machball import ALDIdeal
from machball.ballistic import Trench
from machball.ballistic.ballistic import ballistic_markov
import os
ROOT = Path(__file__).resolve().parents[2]
s = float(sys.argv[1]); NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
AR, N0, T, MASS, Z_STOP, EP, WT = int(os.environ.get('BF_AR', 800)), 4e18, 573.15, 72.09, float(os.environ.get('BF_ZSTOP', 330.0)), 0.05, 0.002
ald = ALDIdeal(beta0=s, MM=MASS, T=T, p0=1.0, s0=1.0 / N0)          # p0 = 1 Pa: time in s equals exposure in Pa s
st = Trench(AR, NZ); av = np.ones(st.N - 1); t = 0.0; total = np.sum(st.areas[1:]); sav = 1.0
z = (np.arange(NZ) + 0.5) / NZ * AR                      # wall sections; the last coverage entry is the end wall
times, frames = [], []; t0 = time.time(); step = 0
while True:
    sp = np.ones(st.N); sp[1:] = s * av
    flux = ballistic_markov(sp, st, [0]); persite = st.areas[0] * s * flux[1:] / st.areas[1:]
    dt = EP / np.amax(persite); av = av / (1 + persite * dt); t += dt; step += 1
    nsav = np.sum(av * st.areas[1:]) / total
    if sav - nsav > WT:
        cov = 1 - av; times.append(t * ald.t0); frames.append(cov.astype(np.float32)); sav = nsav
        i = np.argmax(cov[:NZ] < 0.5); z50 = z[i]
        if len(times) % 20 == 0: print(f"s {s:g}: frame {len(times)}, step {step}, z50 {z50:.0f}, exposure {times[-1]:.3f} Pa s, {time.time()-t0:.0f} s", flush=True)
        if z50 > Z_STOP: break
out = ROOT / "figures" / "ballistic_frames"; out.mkdir(exist_ok=True)
np.savez_compressed(out / f"s_{s:g}_Nz{NZ}.npz", s=s, AR=AR, Nz=NZ, exposure_Pa_s=np.array(times), z=z, cov=np.array(frames)[:, :NZ], end_wall=np.array(frames)[:, NZ])
print(f"done s {s:g}: {len(times)} frames, {step} steps, {time.time()-t0:.0f} s", flush=True)
