# aldconf: feature-scale ALD / PEALD conformality model

What dose and purge give 95% step coverage in a 20:1 to 100:1 aspect-ratio feature, how does that dose scale with aspect ratio, and where does plasma ALD stop being fixable by dose because radicals recombine on the walls first.

Model: Knudsen (plus Bosanquet transition) diffusion along a slit or hole, irreversible Langmuir adsorption with sticking coefficient `s` on open sites, optional radical loss with recombination probability `r` on every wall collision, channel narrowing as the film grows. Quasi-steady gas transport, time-dependent coverage, solved cycle by cycle by finite volumes. Validation against public PillarHall LHAR thickness profiles (Zenodo 3974438, CC BY 4.0).

Status: week 1 of 8. Thermal solver, scaling-law tests and raw-data plotting are in; PEALD recombination, Monte Carlo cross-check, DOE and surrogate follow.

Run:
```
pip install -e .
pytest
python notebooks/plot_zenodo_raw.py
```

References: Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018); Arts et al., J. Vac. Sci. Technol. A 37, 030908 (2019); Arts et al., J. Phys. Chem. C (2019, 2021); Yanguas-Gil and Elam; Gordon et al. (2003); Cremers, Puurunen, Knoops, Appl. Phys. Rev. 6, 021302 (2019); Yim et al., PCCP (2020), doi 10.1039/D0CP03358H.

## Mac setup (Apple Silicon)
```
~/.local/bin/python3.12 -m venv .venv
.venv/bin/pip install -e . viennaps machball pytest
sh scripts/fix_viennaps_libomp.sh .venv     # macOS only: one OpenMP runtime, see the script
.venv/bin/python notebooks/viennaps_smoke.py
```
The `objc[...]: Class vtkCocoa... is implemented in both` lines at import are the two wheels' bundled VTK copies; harmless so far.

## OpenFOAM (reactor scale, macOS)
Native Apple Silicon build from https://github.com/gerlero/openfoam-app (v2.2.3, OpenFOAM v2606, GPL-3.0, ad-hoc signed, not Apple-notarized). Installed as `~/Applications/OpenFOAM-v2606.app`; the app mounts its own read-only volume at `/Volumes/OpenFOAM-v2606`. Run a command inside the OpenFOAM environment:
```
~/Applications/OpenFOAM-v2606.app/Contents/Resources/etc/openfoam -c 'blockMesh && simpleFoam'
```
Smoke case: `reactor/smoke_pitzDaily` (tutorial copy, 50 iterations). ParaView is not bundled; figures are made in Python from the written fields instead.
