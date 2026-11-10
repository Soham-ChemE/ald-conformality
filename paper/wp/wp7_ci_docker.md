# Work package 7: continuous integration, Docker image, tests

Repository state: no commits made; everything below is in the working tree for the integrator.

## What was added

| File | Purpose |
|---|---|
| `.github/workflows/ci.yml` | Three jobs on push to `main`, pull requests and manual dispatch: `python` (ubuntu-latest, Python 3.12, `pip install -e '.[test,crosscheck]'`, the console script, `pytest -q tests`), `openfoam` (container `opencfd/openfoam-default:2512`: Python 3.12 from Ubuntu 24.04, `reproduce.sh --build` (wmake), `reproduce.sh --regress` (steady flow, copy of U/phi/mesh, 2 s regression pulse, regression test), the full suite with the regression test active, solver logs uploaded as the artifact `openfoam-logs`), `pages-check` (`python tests/test_docs_assets.py`: `docs/index.html` exists and every local `src=`/`href=` target and every `data-viz-fig` fallback PNG is present; no browser). |
| `Dockerfile` | `FROM opencfd/openfoam-default:2512`; installs python3/venv/pip; copies the repository (filtered by `.dockerignore`); creates `/work/.venv` with `pip install -e '.[test,crosscheck]'`; compiles the boundary-condition library into `/work/.openfoam-user`; `ENTRYPOINT bash /work/reproduce.sh`, `CMD --tests`. |
| `.dockerignore` | Excludes `.git`, `.venv*`, caches, `private/`, `figures/*.vtp`, `docs/qa`, and under `reactor/` at any depth: time directories other than `0/` (`[1-9]*`, `0.[0-9]*`), `log.*`, `polyMesh`, `postProcessing`, `dynamicCode`, the copied `0/U` and `0/phi` of derived cases (the two steady cases' `0/U` inputs are re-included), compiled `Make/darwin*`, `Make/linux*`, `lnInclude`. The 20 GB of local results stay out of the build context. |
| `reproduce.sh` | Modes `--tests` (default), `--build`, `--flow`, `--regress`, `--quick`, `--full` (see the header comment). Sources the first `/usr/lib/openfoam/openfoam*/etc/bashrc` unless `WM_PROJECT_DIR` is already set (so it also runs inside the macOS launcher: `$OF -c './reproduce.sh --regress'`), puts the user library in `ALDCONF_FOAM_USER_DIR` (default `<repo>/.openfoam-user`, so no home directory needs to be writable), and for `--full` installs a shim at `~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam` with the launcher's `-c COMMAND` interface so the Mac driver scripts (`make_layout_cases.py`, `sweep_sticking_sensitivity.py`) run unchanged on Linux. The bash bundled with the openfoam-app (`<volume>/env/bin/bash`, 5.2) segfaults part-way through the script (the quirk `reactor/make_layout_cases.py` already notes; the launcher then reports exit 0), so the script re-executes itself under `/bin/bash` when it finds itself in that shell; on Linux nothing changes. Validated locally end to end through the real launcher (`--regress`: wmake into `.openfoam-user`, blockMesh, simpleFoam 2000 iterations, copy of U/phi/polyMesh, scalarTransportFoam, regression test passed; a second run reproduced the reference, see below). |
| `pyproject.toml` | Existed already (version 0.1.0). Changed: version `0.9.1` (matches the paper draft), `scikit-learn` added to the dependencies (the surrogate refit `analysis/doe_gp.py` and `tests/test_surrogate.py` need it), `[project.scripts] aldconf = "aldconf.cli:main"`, extras `test = ["pytest"]` and `crosscheck = ["machball"]`. `pip install -e .` from the repository root verified in a throwaway venv: `aldconf rule ...` prints 14.385 s, `aldconf surrogate ...` 11.859 s, and the suite imports the package from `src/`. `PYTHONPATH=src` still works (the tests were run that way). |
| `reactor/showerhead_steady/system/controlDict` | `endTime ${ALD_STEADY_END_TIME:-2000};` OpenFOAM expands the environment variable with a default (verified with `foamDictionary`: 2000 unset, 50 when set). CI and the regression reference use 2000; the override is for exploration only, since U and phi at the final iteration feed the regression case. Locally the 2000 iterations take 19 s, so CI does not need it. |
| `coupling/postprocess_layout.py` | Four pure functions extracted (`local_wall_ratio`, `area_weighted_mean`, `supply_rule`, `open_area_fraction`) and the script body moved under `main()` so the module can be imported by the tests without running the cases. Output verified unchanged: `figures/pattern_layout.csv` and `figures/pattern_layout.png` byte-identical before and after. |
| `tests/conftest.py` | Puts `src/` and `coupling/` on `sys.path` so the suite runs both with `PYTHONPATH=src` and after `pip install -e .`. |
| `tests/reference/regress_t2.csv` | Wafer patch of `reactor/showerhead_pulse_layout_regress` at t = 2 s: face, thetaFlat, thetaBottom, value (150 faces, 8 significant digits). Made from the chain CI runs: `reactor/showerhead_steady` to 2000 iterations, its U, phi and mesh copied into the regress case, 1 s pulse, 2 s run (see "Frozen flow" below for why the case's earlier local results were not used). |
| `tests/reference/regress_t1_thetaBottom_first20.csv` | thetaBottom of the first 20 faces at t = 1.0 s, mid-pulse (0.99580 to 0.99586: nearly, not fully, saturated; the discriminating check, since at 2 s both coverages are 1 to 8 digits). |
| `tests/reference/make_regress_reference.py` | Regenerates both CSVs from a local run with `coupling/patterned_lib.load_patterned` (`PYTHONPATH=src .venv/bin/python tests/reference/make_regress_reference.py`). |
| `tests/fixtures/patterned_case/{0,0.5,1}/T` | Static three-time `patternedSurface` case: the author's one-line dictionary with a `holeDensityProfile` table and a 1.5 s pulse (t = 0), every wafer field in the `uniform` scalar form with `ARField`/`holeDensityField` absent (t = 0.5), `nonuniform List<scalar>` fields with the per-face pattern present and known linear ramps in the face index (t = 1). Other patches carry their own `value` entries so the block isolation is tested. |
| `tests/test_patterned_lib.py` | 12 tests of the parser on the fixture: times, pulse, shapes, scalar parameters, zero fill at t = 0, uniform broadcast, fallback of the pattern fields to the scalar AR/holeDensity, face-by-face lists, the profile ramp, block isolation, the `_block`/`_field`/`_scalar` helpers. |
| `tests/test_cli.py` | 9 tests running `python -m aldconf` by subprocess with `PYTHONPATH=src`: `rule` prints 14.38 s (within 0.05) with wall ratio 55.92 and U 0.0501; `--factor`; `surrogate` at the README point x 0.012, p 316 gives 11.86 s within 2 %, at the exact coded centre 11.62 s (MATLAB) within 2 %; `pareto --max-lag 5` header, 41 points, sorted by t90, precursor/min non-increasing and equal to 1/util; `--n`; usage error and `--help`. |
| `tests/test_layout_postprocess.py` | 8 tests of the extracted functions: uniform 1.78e13 /m2, AR 100 gives 55.9; scalar inputs; area weighting (holes inside r < 75 mm give a quarter of 55.92, not half); weights sum to the wafer area; the rule gives 14.39 s for the reference wafer and 0.25 s blanket; elementwise rule on a profile; 14 % open area; `t90_profile` first crossing and mask. |
| `tests/test_regression_reference.py` | 3 tests; skip with the reason "reactor/showerhead_pulse_layout_regress/2/T not found: run the steady case and the regression case first" when the case has not been run. Checks the setup (1 s pulse, AR 50, holeDensity 1, 2 s), thetaFlat/thetaBottom/value at t = 2 against the reference to `ALDCONF_REGRESS_ATOL` (default 1e-6, absolute), and the mid-pulse thetaBottom. |
| `tests/test_rule_consistency.py` | 24 tests: `rule_pulse` against a line-for-line port of the live page's `rule()` (5 parameter sets, all five outputs to 1e-12), U = 0.0501 m/s, the 16 rows of `figures/pulse_time_map.csv` (1.2 x t_supply equals t_pulse within 0.5 % and the wall ratio to 1e-9, one test per row), and the band t90_worst/t_supply in (1.0, 1.3) for every patterned run. |
| `tests/test_docs_assets.py` | 3 tests, standard library only, also runnable as a script (the pages-check job): `docs/index.html` exists, references `assets/site.js`, `site.css`, `data.js`, and all 32 local references (including the light/dark fallback PNG pair of each of the 13 `data-viz-fig` panels) are present. |

## Docker Hub tag decision

Queried `https://hub.docker.com/v2/repositories/opencfd/openfoam-default/tags?page_size=100` on 2026-11-10: 12 tags, `latest`, `2512` (2026-01-23, linux/amd64 + arm64), `2506`, `2412`, `2406`, `2312`, `2306`, `2212`, `2206`, `2112`, `2106`, `2012`. There is no `2606` image. Chosen: **`opencfd/openfoam-default:2512`**, the newest release tag, pinned rather than `latest` so the CI stays reproducible. The image config (read from the registry manifest) shows Ubuntu 24.04, packages `openfoam2512`, `openfoam2512-dev` (wmake and headers) and `openfoam2512-default`, entrypoint `/openfoam/run`, and **no `USER` instruction**: the container runs as root unless `docker run -u` is used (the OpenCFD `openfoam-docker` script does that; GitHub Actions does not). The workflow therefore sets `options: --user root` explicitly (needed for `apt-get`) and still keeps the user library inside the workspace (`ALDCONF_FOAM_USER_DIR`), so nothing depends on a writable home directory. Python is not in the image; the Ubuntu 24.04 `python3` package is 3.12, the project's version.

Library naming: the cases say `libs ("libaldWallKinetics.dylib")`. OpenFOAM's `dlOpen` (src/OSspecific/POSIX/POSIX.C, v2606) replaces a foreign extension with the platform's (`.dylib` -> `.so` on Linux), so the cases need no edit; `reproduce.sh --build` also links `libaldWallKinetics.dylib -> libaldWallKinetics.so` in the user library directory as a second line of defence.

## Tests

Local run, `PYTHONPATH=src .venv/bin/python -m pytest -q tests` (regression case present, machball installed):

```
77 passed in 15.39s
```

Per file: test_scaling 5, test_surrogate 5, test_machball_crosscheck 2 (the 12 that existed); this work package adds
test_patterned_lib 12, test_cli 9, test_layout_postprocess 8, test_regression_reference 3, test_rule_consistency 24,
test_docs_assets 3 (59); `tests/test_layout_tool.py` (6, Python half; its Playwright half skips without a browser) was
added by another work package during the same session and is counted in the 77. `./reproduce.sh --tests` gives the same
77. In the `python` CI job the 3 regression tests skip with the reason above (74 passed, 3 skipped expected); in the
`openfoam` job all 77 run.

## Frozen flow: a finding the integrator should know

The patterned pulse cases (`showerhead_pulse_patterned_check`, `_AR50`, `_AR100_long`, and therefore the regress case
copied from `_check`) carry a frozen flow in `0/U`, `0/phi` that is NOT the 2000-iteration solution of
`reactor/showerhead_steady`: U differs by up to 9.7e-4 m/s (0.26 % of the 0.37 m/s maximum), phi by 6e-9 m3/s. The
blanket case `showerhead_pulse/0/U` IS byte-identical to the fresh steady solution (the `reactor/Allrun` chain). The
patterned cases' flow is shared by every case the map driver built (`reactor/map/*`, the `showerhead_h2o_*` cases, 13,927
identical copies found by a scan of every U file under `reactor/`), but no steady case in the tree holds it as its final
iteration: the steady run it came from was overwritten or never kept. So it cannot be re-derived, which settles the
choice of reference. Consequence: a reference made from the case's existing local results
could never be met by CI, which builds the flow from `showerhead_steady` (the first local `--regress` run failed by
1.7e-4 in `value` and 1.2e-5 in mid-pulse thetaBottom for exactly this reason). The reference CSVs were therefore
regenerated from the CI chain itself, and a second independent `--regress` run reproduced them (3 passed, differences 0). The 0.26 % flow
difference changes the saturation times by well under one write interval and none of the paper's numbers; it is noted
here only so that nobody chases it later.

## README text to insert

After the "Try the surrogate in Python" subsection (or at the end of "Reproduce"):

```markdown
### Continuous integration and Docker

![CI](https://github.com/Soham-ChemE/ald-conformality/actions/workflows/ci.yml/badge.svg)

Every push runs three jobs: the Python feature solver and its 77 tests on Python 3.12; the OpenFOAM chain in
OpenCFD's `opencfd/openfoam-default:2512` container (the two boundary conditions compiled with wmake, the steady
chamber flow, a 2 s TMA pulse over a patterned wafer, and a regression test of the wafer coverage against a stored
reference to 1e-6); and a check that the live page ships every asset it references.

The same chain runs anywhere Docker runs, no OpenFOAM or Python install needed:

```bash
docker build -t ald-conformality . && docker run --rm ald-conformality        # the test suite, ~1 min after the build
docker run --rm ald-conformality --quick                                      # + steady flow, regression pulse and
                                                                              #   the Python analyses, 10-20 min
```

`docker run --rm ald-conformality --full` reruns every chamber study (hours). `reproduce.sh` in the repository root is
the entry point the image uses; it also runs natively inside the OpenFOAM launcher.
```

Also update the test count in the existing Reproduce block: `PYTHONPATH=src .venv/bin/python -m pytest -q      # 7 tests` -> `# 77 tests`.

## For the integrator to check after the first push

1. The CI run itself (not observable from here): https://github.com/Soham-ChemE/ald-conformality/actions. The `openfoam` job is the one that can fail for reasons outside the repository: `apt-get` inside the container (network), the deb package layout (`/usr/lib/openfoam/openfoam2512/etc/bashrc` is what `reproduce.sh` looks for), wmake against v2512 headers (the two BCs use only `mixedFvPatchScalarField`, `Function1`, `interpolateXY`-style code that has been stable for years), and the regression tolerance across compilers (reference made with clang on Apple Silicon, CI compares gcc/x86-64 output to 1e-6 absolute on quantities of order 1 written with 8 significant digits). If only the tolerance fails, the artifact `openfoam-logs` and the test message give the max difference and face; `ALDCONF_REGRESS_ATOL` can be set in the workflow env (1e-5 would still be a meaningful regression test) rather than changing the test.
2. Add `.openfoam-user/` to `.gitignore` (not owned by this work package): the local validation runs created `<repo>/.openfoam-user/` (the compiled library, 240 kB) and it currently shows as untracked. It must not be committed.
3. `pip install -e .` now pulls scikit-learn; the README's venv line already installs everything through `-e .`.
4. `docker build` was not run here (Docker is not installed on this machine); the Dockerfile was checked by reading: every `COPY` source exists (`pyproject.toml`, `README.md`, `LICENSE`, `src`), the base tag exists on Docker Hub for amd64 and arm64, the apt package names are Ubuntu 24.04's, and the build step is the same `reproduce.sh --build` that ran locally. First `docker build` on a Mac or Linux box is the remaining check; expect 2-4 minutes (pip wheels for numpy/scipy/pandas/matplotlib/scikit-learn, then wmake).
5. The uploaded artifact is named `openfoam-logs`; `if-no-files-found: warn` keeps a failed build step from hiding the logs that do exist.
6. `docs/index.html` shows as modified in the working tree; that change is not from this work package (docs is not owned here) and the pages check passes against it (32 local references, 0 missing).
7. Not covered by `--full`: the water half-cycle, temperature, edge-ring and hole-size studies, whose drivers (`make_h2o_cases.py`, `make_h2o_desorption_cases.py`, `make_temperature_cases.py` + `run_temperature_cases.sh`, the `edgering_*` and `holesize_*` cases) document no single command; they remain as documented in CLAIMS.md.
