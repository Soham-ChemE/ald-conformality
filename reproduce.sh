#!/usr/bin/env bash
# One entry point for reproducing ald-conformality, used by the Dockerfile and by .github/workflows/ci.yml.
#
#   ./reproduce.sh --tests     pytest over tests/ (default; the OpenFOAM regression test skips unless --regress ran)
#   ./reproduce.sh --build     compile the boundary-condition library reactor/lib/aldWallKinetics with wmake
#   ./reproduce.sh --flow      steady N2 flow: blockMesh + simpleFoam in reactor/showerhead_steady (2000 iterations,
#                              1-3 min; ALD_STEADY_END_TIME=<n> shortens it for exploration only, the regression
#                              reference was made with 2000)
#   ./reproduce.sh --regress   --build, --flow, then copy U, phi and the mesh into
#                              reactor/showerhead_pulse_layout_regress and run its 2 s TMA pulse (about 30 s), then
#                              the regression test against tests/reference/
#   ./reproduce.sh --quick     --regress, the Python analyses that need no OpenFOAM run (analysis/doe_gp.py,
#                              coupling/validate_supply_rule.py, coupling/temporal_vs_spatial.py,
#                              coupling/peald_ceiling.py) and the full pytest. 10-20 minutes.
#   ./reproduce.sh --full      everything OpenFOAM: the blanket chain (reactor/Allrun), the 16-run pulse-time map,
#                              the operating sweep, the 25-run DOE and the two layout cases, then their
#                              post-processing and --quick. HOURS on one core (the map alone is about 25 min on an
#                              M-series Mac; the DOE several hours).
#
# OpenFOAM: the script sources the first /usr/lib/openfoam/openfoam*/etc/bashrc it finds (the OpenCFD packages and
# Docker image) unless WM_PROJECT_DIR is already set (running inside the macOS launcher:
#   ~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c './reproduce.sh --regress'
# ) or FOAM_BASHRC names another bashrc. The user library goes to ALDCONF_FOAM_USER_DIR (default
# <repo>/.openfoam-user) so no home directory needs to be writable.
# Python: ALDCONF_PYTHON (default <repo>/.venv/bin/python, falling back to python3 on PATH).
set -euo pipefail

# The bash bundled with the macOS openfoam-app (<volume>/env/bin/bash) segfaults part-way through this script
# (reactor/make_layout_cases.py notes the same); re-run under the system bash, which the launcher environment keeps.
case "${BASH:-}" in
    */env/bin/bash) [ -x /bin/bash ] && [ -z "${ALDCONF_REEXEC:-}" ] && ALDCONF_REEXEC=1 exec /bin/bash "$0" "$@" ;;
esac

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODE=${1:---tests}
PY=${ALDCONF_PYTHON:-$ROOT/.venv/bin/python}
[ -x "$PY" ] || PY=$(command -v python3)
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg

say()  { printf '\n==> %s\n' "$*"; }
die()  { printf 'reproduce.sh: %s\n' "$*" >&2; exit 1; }

foam_env() {
    if [ -z "${WM_PROJECT_DIR:-}" ]; then
        local rc=${FOAM_BASHRC:-}
        if [ -z "$rc" ]; then
            rc=$(ls -d /usr/lib/openfoam/openfoam*/etc/bashrc 2>/dev/null | sort | tail -1 || true)
        fi
        [ -n "$rc" ] && [ -f "$rc" ] || die "no OpenFOAM found: set FOAM_BASHRC or run inside the OpenFOAM launcher"
        set +u
        # shellcheck disable=SC1090
        source "$rc" || true          # the bashrc returns non-zero when it cannot detect the shell; harmless
        set -u
    fi
    [ -n "${WM_OPTIONS:-}" ] || die "OpenFOAM environment incomplete (WM_OPTIONS unset)"
    export WM_PROJECT_USER_DIR=${ALDCONF_FOAM_USER_DIR:-$ROOT/.openfoam-user}
    export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
    export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
    mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN"
    export PATH="$FOAM_USER_APPBIN:$PATH"
    export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
    export DYLD_LIBRARY_PATH="$FOAM_USER_LIBBIN${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
    say "OpenFOAM $WM_PROJECT_VERSION ($WM_OPTIONS), user libraries in $FOAM_USER_LIBBIN"
}

# The driver scripts written on the Mac call the openfoam-app launcher by its absolute path; on Linux provide a
# shim with the same interface (launcher -c 'command') so they run unchanged.
launcher_shim() {
    local L="$HOME/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam"
    [ -e "$L" ] && return 0
    mkdir -p "$(dirname "$L")"
    cat > "$L" <<SHIM
#!/usr/bin/env bash
# shim for the macOS openfoam-app launcher, written by reproduce.sh: runs 'launcher -c COMMAND' in this environment
[ "\$1" = "-c" ] || { echo "shim: expected -c COMMAND" >&2; exit 2; }
export PATH="$PATH" LD_LIBRARY_PATH="$LD_LIBRARY_PATH" WM_PROJECT_USER_DIR="$WM_PROJECT_USER_DIR" \\
       FOAM_USER_LIBBIN="$FOAM_USER_LIBBIN" FOAM_USER_APPBIN="$FOAM_USER_APPBIN" PYTHONPATH="$PYTHONPATH"
exec bash -c "\$2"
SHIM
    chmod +x "$L"
    say "launcher shim installed at $L"
}

build_lib() {
    foam_env
    say "wmake reactor/lib/aldWallKinetics"
    ( cd "$ROOT/reactor/lib/aldWallKinetics" && wmake libso . 2>&1 | tee log.wmake | tail -3 )
    # The cases name the library as libaldWallKinetics.dylib (written on macOS). OpenFOAM's dlOpen remaps the
    # extension to the platform's, but a same-named link keeps the load path short and unambiguous on Linux.
    if [ -f "$FOAM_USER_LIBBIN/libaldWallKinetics.so" ] && [ ! -e "$FOAM_USER_LIBBIN/libaldWallKinetics.dylib" ]; then
        ln -s libaldWallKinetics.so "$FOAM_USER_LIBBIN/libaldWallKinetics.dylib"
    fi
    ls -l "$FOAM_USER_LIBBIN"/libaldWallKinetics.* || die "library did not build"
}

clean_case() {   # remove every result of an earlier run, keep 0/ and the inputs
    ( cd "$1" && for d in [0-9]*; do [ "$d" = "0" ] || rm -rf "$d"; done; rm -rf postProcessing log.* dynamicCode )
}

run_flow() {
    foam_env
    local C="$ROOT/reactor/showerhead_steady"
    say "steady N2 flow in reactor/showerhead_steady (endTime ${ALD_STEADY_END_TIME:-2000})"
    clean_case "$C"
    ( cd "$C" && blockMesh > log.blockMesh 2>&1 && checkMesh > log.checkMesh 2>&1 && simpleFoam > log.simpleFoam 2>&1 ) \
        || { tail -20 "$C"/log.simpleFoam 2>/dev/null; die "steady flow failed (see reactor/showerhead_steady/log.*)"; }
    grep -E "^ExecutionTime" "$C/log.simpleFoam" | tail -1
}

last_time() {   # largest numeric time directory of a case
    ( cd "$1" && ls -d [0-9]* 2>/dev/null | grep -E '^[0-9.]+$' | sort -g | tail -1 )
}

run_regress() {
    build_lib
    run_flow
    local S="$ROOT/reactor/showerhead_steady" C="$ROOT/reactor/showerhead_pulse_layout_regress"
    local LAST; LAST=$(last_time "$S")
    [ -n "$LAST" ] && [ -f "$S/$LAST/U" ] && [ -f "$S/$LAST/phi" ] || die "no U/phi written by the steady case"
    say "regression case: frozen flow from showerhead_steady/$LAST, 1 s TMA pulse, 2 s run"
    clean_case "$C"
    rm -rf "$C/constant/polyMesh" && cp -r "$S/constant/polyMesh" "$C/constant/"
    cp "$S/$LAST/U" "$S/$LAST/phi" "$C/0/"
    ( cd "$C" && scalarTransportFoam > log.scalarTransportFoam 2>&1 ) \
        || { tail -20 "$C/log.scalarTransportFoam"; die "scalarTransportFoam failed in the regression case"; }
    [ -f "$C/2/T" ] || die "regression case did not reach t = 2 s"
    say "regression test"
    "$PY" -m pytest -q "$ROOT/tests/test_regression_reference.py"
}

run_python_analyses() {
    say "Python analyses that need no OpenFOAM run"
    ( cd "$ROOT" && "$PY" analysis/doe_gp.py | tail -5 )
    ( cd "$ROOT" && "$PY" coupling/validate_supply_rule.py | tail -3 )
    ( cd "$ROOT" && "$PY" coupling/temporal_vs_spatial.py | tail -3 )
    ( cd "$ROOT" && "$PY" coupling/peald_ceiling.py | tail -3 )
}

run_tests() {
    say "pytest ($("$PY" --version 2>&1))"
    ( cd "$ROOT" && "$PY" -m pytest -q tests )
}

run_full() {
    echo "WARNING: --full reruns every chamber study: the blanket chain, the 16-run pulse-time map (~25 min),"
    echo "         the operating sweep, the 25-run DOE and the layout cases. Expect HOURS on one core."
    build_lib
    launcher_shim
    say "blanket-wafer chain (reactor/Allrun)";          ( cd "$ROOT" && sh reactor/Allrun )
    say "pulse-time map, 16 runs";                        ( cd "$ROOT" && "$PY" reactor/sweep_pattern_map.py )
    say "operating sweep";                                ( cd "$ROOT" && "$PY" reactor/sweep_operating.py )
    say "25-run DOE";                                     ( cd "$ROOT" && "$PY" reactor/doe_patterned.py )
    say "spatially varying layouts";                      ( cd "$ROOT" && "$PY" reactor/make_layout_cases.py --run )
    say "post-processing"
    ( cd "$ROOT" && "$PY" reactor/postprocess_pulse.py | tail -3 )
    ( cd "$ROOT" && "$PY" reactor/plot_pulse_time_map.py | tail -3 )
    ( cd "$ROOT" && "$PY" reactor/plot_sweep_operating.py | tail -3 )
    ( cd "$ROOT" && "$PY" coupling/postprocess_layout.py | tail -3 )
    run_regress
    run_python_analyses
    run_tests
}

case "$MODE" in
    --tests)   run_tests ;;
    --build)   build_lib ;;
    --flow)    run_flow ;;
    --regress) run_regress ;;
    --quick)   run_regress; run_python_analyses; run_tests ;;
    --full)    run_full ;;
    -h|--help) sed -n '2,30p' "$0" ;;
    *) die "unknown mode $MODE (use --tests, --build, --flow, --regress, --quick, --full)" ;;
esac
say "done ($MODE)"
