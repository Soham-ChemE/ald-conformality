#!/bin/sh
# Finite sticking on the wafer: kinetic flux j = s (vbar/4) c_face balances diffusion D (c_cell - c_face)/delta.
# As an OpenFOAM mixed BC (refValue 0, refGradient 0): valueFraction f = k / (k + D/delta), k = s vbar/4.
# TMA at 573 K: vbar = sqrt(8 R T / (pi M)) = 410 m/s -> k = 102.6 s m/s. D/delta = 8.4e-3 / 2.5e-4 = 33.6 m/s.
# Usage: ~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c 'sh reactor/sweep_sticking.sh'
set -e
cd "$(dirname "$0")"
for s in 1 0.1 0.01 0.001; do
  f=$(python3 -c "k=$s*410.0/4; print(f'{k/(k+33.6):.6g}')")
  C=showerhead_pulse_s$s
  rm -rf "$C"; mkdir -p "$C"; cp -r showerhead_pulse/0 showerhead_pulse/constant showerhead_pulse/system "$C"/
  python3 - "$C/0/T" "$f" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1]); s = p.read_text()
old = "    wafer      { type fixedValue; value uniform 0; }     // perfectly absorbing wafer (sticking -> 1 limit)"
new = f"    wafer      {{ type mixed; refValue uniform 0; refGradient uniform 0; valueFraction uniform {sys.argv[2]}; value uniform 0; }}  // finite sticking, see sweep_sticking.sh"
assert old in s; p.write_text(s.replace(old, new))
PY
  ( cd "$C" && scalarTransportFoam > log.scalarTransportFoam 2>&1 && echo "s=$s f=$f done: $(grep '^ExecutionTime' log.scalarTransportFoam | tail -1)" )
done
