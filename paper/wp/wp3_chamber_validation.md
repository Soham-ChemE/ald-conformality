# Work package 3: validation at chamber scale (2026-11-10)

## What was done

A. The blanket-wafer saturating TMA pulse of the chamber model (0.345 s at the centre, 0.385 s worst, at 300 Pa,
573 K, 300 sccm N2, 1 % TMA = 3 Pa, s = 5e-3, N0 = 4e18 /m2) was converted to an exposure and set beside published
Al2O3 ALD saturation data from seven independent sources (Ylilammi 2018; Arts 2019; Vandalon and Kessels 2017/2019;
Gakis 2018/2019; Soto and Tysoe 1991 and Dillon 1995 as summarised in Puurunen's 2005 review; the PillarHall record
already used in this repository). Two exposures are kept apart throughout: the feed exposure p_feed x t90 that a recipe
sets (1.03 Pa s = 7,760 L for our run) and the surface exposure, the time integral of the TMA partial pressure at the
wafer up to 90 % coverage (0.146 Pa s = 1,090 L), which is what the Langmuir kinetics see. For every source the
exposure the experiment applied is tabulated as stated and the planar 90 %-coverage exposure that the source's own
sticking coefficient or reaction cross-section implies is computed, E90 = ln(10) N0/(s F1) with F1 the impingement
rate per pascal, or ln(10)/(sigma F1) when a cross-section is given (no site density needed). No new chamber run was
made for this part; the 1/s scaling of the kinetics-limited exposure gives the sticking coefficient the model would
need to match each source. Script analysis/blanket_saturation_literature.py; outputs figures/blanket_saturation_literature.csv and .png.

B. The reactor-scale CFD-plus-measurement case of Gakis, Vergnes, Scheid, Vahlas, Caussat and Boudouvis was located
and read in full (author manuscripts: Chem. Eng. Res. Des. 132, 795 (2018), hal-01730196 / OATAO 19671; Chem. Eng.
Sci. 195, 399 (2019), hal-02329135 / OATAO 21237; the HAL and core.ac.uk servers sit behind a bot check, the copies were
fetched through the Internet Archive). An axisymmetric approximation of their Ultratech Fiji F200 chamber was built
in this repository's framework under reactor/validation_gakis/ (blockMesh wedge, simpleFoam steady Ar flow,
scalarTransportFoam TMA pulses with the langmuirSticking wafer), with the chamber dimensions measured on their
published cross-section because the papers give none, their operating conditions, their computed TMA pulse (amount
and shape), and s and N0 converted from their kinetic model. The predicted coverage profile across the 200 mm wafer,
the non-uniformity, the wafer exposure and the saturation state are compared with their measured thickness profile
(digitised from their Fig. 10) and with their own 3D CFD. Scripts reactor/validation_gakis/make_gakis_cases.py (build and
run) and postprocess_gakis.py; outputs figures/validation_gakis.png and .csv.

## Results

### A. Blanket saturation against published data (figures/blanket_saturation_literature.csv)

Our run (reactor/showerhead_pulse_langmuir_p0.5, 0.5 s pulse): p_feed = 3.00 Pa, t90 = 0.345 s at the centre
(0.385 s at r = 149 mm). Feed exposure p_feed t90 = 1.03 Pa s = 7,760 L (worst point 1.16 Pa s). Surface exposure
(integral of the wafer-face TMA partial pressure to theta = 0.9) = 0.146 Pa s = 1,090 L at the centre and at the
edge alike; the surface sees on average 0.14 of the feed pressure because the wafer draws the gas down (Night 1
part 2: 0.24 c0 at the surface on the blanket wafer). The surface exposure equals the kinetic Langmuir value
ln(10) N0/(s F1) = 0.142 Pa s at s = 5e-3, as it must: coverage depends only on the integrated impingement.

| Source (DOI) | T (C) | What the number is | Exposure (Pa s) | Langmuir | s this model would need |
|---|---|---|---|---|---|
| This work, chamber model, centre | 300 | feed p x t90 | 1.03 | 7,760 | 5e-3 (as run) |
| This work, chamber model, centre | 300 | surface integral p_s dt to theta = 0.9 | 0.146 | 1,090 | 5e-3 (as run) |
| Ylilammi, Ylivaara, Puurunen 2018 (10.1063/1.5028178), Picosun R-150, 300 Pa, 0.1 s pulse | 300 | applied: fitted entrance pressure 147 Pa x 0.1 s (Table I); flat region saturated | 14.7 | 110,000 | upper bound only |
| same, implied by their fitted gamma = 5.72e-3 (q = 5e18 /m2, their worked example) | 300 | kinetic E90 | 0.155 | 1,160 | 4.7e-3 |
| Arts et al. 2019 (10.1116/1.5093620), Oxford OpAL, set-point 275 C | ~220 | applied: 190 mTorr s TMA dose, TMA-limited LHAR profile; planar surface saturated | 25.3 | 190,000 | upper bound only |
| same, 1200 mTorr s (H2O-limited profile) | ~220 | applied | 160 | 1,200,000 | upper bound only |
| same, implied by s0 = 2e-3 / 5e-4 (Table I range; N0 = 4e18 assumed) | ~220 | kinetic E90 | 0.33 / 1.32 | 2,470 / 9,890 | 2.2e-3 / 5.5e-4 |
| Vandalon and Kessels 2019 (10.1021/acs.langmuir.9b01600), SFG, sigma_TMA = 6.5 +/- 0.6 e-18 cm2 (steady growth); temperature-independent 100-300 C per their 2017 paper (10.1116/1.4993597) | 300 (100) | kinetic E90 from sigma, no site density needed | 0.273 (0.221) | 2,050 (1,650) | 2.7e-3 (3.3e-3) |
| Gakis et al. 2018 (10.1016/j.cherd.2018.02.031), Fiji F200, 9.6 Pa, 25 ms valve: wafer impingement exposure, read from the colour scale of their Fig. 12, 0.0118-0.0130 mol/m2 (+/- 0.0002) x 46.45 Pa s per mol/m2 | 300 | CFD exposure that gave GPC 0.975-0.998 A/cycle, i.e. >= 97 % saturated (2019 Fig. 10, digitised, +/- 0.003) | 0.55-0.60 | 4,100-4,500 | upper bound only |
| Gakis et al. 2019 (10.1016/j.ces.2018.09.037), implied by fitted s0 = 0.004 x 86 % reaction probability at 300 C, N0 = 4.1e18 from their GPC 1.0 A/cycle at 3500 kg/m3 | 300 | kinetic E90 | 0.212 | 1,590 | 3.4e-3 |
| PillarHall record used here (Zenodo 3974438), fit coupling/fit_pillarhall.py | 300 | applied: fitted mouth exposure 75 Pa x 0.1 s | 7.5 | 56,000 | upper bound only |
| same, implied by the fitted s = 8.5e-3 | 300 | kinetic E90 | 0.084 | 630 | 8.7e-3 |
| Soto and Tysoe 1991, UHV, as summarised by Puurunen 2005 (10.1063/1.1940727): reaction probability 0.01 | 25 | kinetic E90 (N0 = 4e18 assumed) | 0.051 | 380 | 1.4e-2 |
| Dillon et al. 1995, porous alumina FTIR, as summarised by Puurunen 2005: k_a = 0.0076 /s at 0.01 Torr | 227 | first-order rate on a porous membrane: transport into the membrane, not planar kinetics | 404 | 3,000,000 | (1.8e-6, not meaningful) |

Reading. On the like-for-like basis, surface exposure, our 0.146 Pa s (1,090 L) sits inside the literature spread
of kinetic saturation exposures, 0.05 to 1.3 Pa s (Soto 0.05; PillarHall fit 0.08; Ylilammi 0.16; Gakis 0.21;
Vandalon 0.27; Arts 0.33 to 1.3), in the middle of the cluster at 300 C (0.08 to 0.27 Pa s). The feed exposure
7,760 L also falls inside the range of exposures that saturating recipes actually apply (0.56 to 160 Pa s = 4,200 to
1.2 million L), but that comparison says little: applied exposures are upper bounds chosen to over-dose, and feed
exposures overstate what the surface sees by the drawdown factor (7x in our chamber at 1 % TMA). To match each
source's kinetic value the model would need s = 2.7e-3 (Vandalon), 3.4e-3 (Gakis), 4.7e-3 (Ylilammi), 2.2e-3 to 5.5e-4
(Arts), 8.7e-3 (PillarHall fit): our 5e-3 is a factor 1.8 above the Vandalon/Gakis pair and a factor 1.7 below the
PillarHall fit. The two classic George-group curves (Ott 1997, Groner 2004) could not be read: their full texts are
not open and the secondary sources give pulse times without partial pressures; Dillon 1995 is in the table through
Puurunen's summary but is a porous-sample rate and is excluded from the verdict. Rahtu et al. 2001 (F-120, saturation
at about 3 s, Puurunen Fig. 26) has no partial pressure and is not plotted.

### B. Gakis et al. Fiji F200 reproduction (figures/validation_gakis.csv)

What the papers give. Reactor: Ultratech (Cambridge NanoTech) Fiji F200, 200 mm Si(100) wafer on a heated chuck,
three Ar inlets (top 100 sccm, loading door 50 sccm, side 30 sccm carrying the reactants from a heated feeding line;
MFC reference 24 C and 1 bar), turbo pump through an automatic pressure controller, outlet below the chuck; pressure
72 mTorr = 9.6 Pa (Pirani 50 mTorr corrected by their Eq. 5); substrate 300 C, walls 270 C, loading door 36 C; TMA
valve 25 ms (experiment 11) or 60 ms (experiment 12), H2O valve 60-100 ms, purges 5 s, 550 cycles, 50-60 nm films.
Their feeding-system CFD gives the TMA pulse: 0.63 mg = 8.74e-6 mol per 25 ms opening, flow peaking at about 450
sccm at 15 ms, mole fraction 0.9 falling to 0.1 by 50 ms, tail to 0.3 s (their Fig. 8, read off the figure, about
20 % on each ordinate). Kinetic model (2019, Table 2): reversible TMA adsorption on OH, s0,TMA = 0.004 (fitted),
desorption barrier 0.61 eV, reaction barrier 0.52 eV, pre-exponentials kT/h, reaction probability of adsorbed TMA
86 % at 300 C (their Section 4.1.2); OH site density 1.04e-5 mol/m2 = 6.3e18 /m2 at 300 C (their linear fit to
Haukka and Root 1994, two TMA per three OH); GPC 1.0 A/cycle at 300 C with alumina density 3500 kg/m3, i.e. 4.1e18 Al
(and TMA) per m2 per cycle. Measurements: GPC along the diameter for experiment 11 (2019 Fig. 10, squares): 0.975 at
-95.5 mm rising to 0.998 at +37 mm and 0.983 at +94.7 mm (loading-door side), digitised here by connected-component
centroids of the markers with +/- 0.003 A/cycle; non-uniformity (max - min)/min = 2.37 % digitised against 2.40 %
stated; 60 ms: 0.58 % stated, map 0.98 to 1.01 A/cycle, GPC 1.0 A/cycle. Their own CFD plus kinetics: 2.45 % (25 ms)
and 0.83 % (60 ms); their wafer TMA impingement exposure for the 25 ms pulse 0.0118 to 0.0130 mol/m2 (2018 Fig. 12,
read from the colour scale, +/- 0.0002) = 0.55 to 0.60 Pa s, 10.7 % non-uniform; surface TMA partial pressure at
25 ms 0.005 to 0.028 Torr = 0.7 to 3.7 Pa (2019 Fig. 7b). Not given: any chamber dimension, the 60 ms TMA amount,
the measured 60 ms diameter profile as numbers.

What was built (every approximation). Dimensions measured on the 2018 Fig. 5a cross-section with the 200 mm chuck
top as the scale (1.65 px/mm, +/- 3 mm): inlet tube r = 20 mm, 82 mm long; conical ceiling from r = 20 mm at 199 mm
above the wafer to r = 128 mm at 35 mm above it; wide chamber r = 128 mm from +35 to -50 mm (the real taper at -33 to
-69 mm is a step); chuck a cylinder r = 100 mm, 50 mm deep (the real chuck is a dome 124 mm deep, r = 95 mm at -47 mm
and 59 mm at -96 mm; its lower half is omitted); pumping duct r = 106 mm from -50 to -205 mm, outlet on its floor.
The three inlets become one axial inlet carrying 180 sccm and the TMA; the loading-door box, the off-axis side-inlet
tube and the pulse-driven recirculation are not represented; the carrier flow is frozen at its steady value although
the real pulse is a slug of nearly pure TMA at up to 450 sccm for 40 ms (their outlet pressure rises 18 %); the TMA is
a passive scalar whose peak partial pressure at the inlet (11.5 Pa for the 25 ms amount) exceeds the 9.6 Pa of Ar;
isothermal 573 K (walls 543 K in reality); D_TMA = 0.263 m2/s is the N2 value scaled 1/p and used for Ar; mean free
path 1.4 mm against a 2 mm first cell. Flow: Q = 0.0603 m3/s at chamber conditions (their pump speed 57.1 L/s at the
543 K outlet is the same number), 48 m/s in the inlet tube (Re 4.3; their Re 2.05 at the side inlet during purge),
residence time 0.23 s, 13,020 cells, simpleFoam 2000 iterations, inlet and outlet flux 0.0602 m3/s (balance to 1e-6).
Wafer: langmuirSticking, s = 0.004 x 0.86 = 3.44e-3 (their adsorption probability times the fraction of adsorbed TMA
that reacts rather than desorbs, the irreversible-Langmuir equivalent of their reversible step), N0 = 4.1e18 /m2,
vbar 410 m/s; one case at this repository's s = 5e-3. Pulse: their 25 ms amount with their Fig. 8 shape (table in
the script, integral renormalised to 0.63 mg); 60 ms at 1.5x and 2.4x that amount (headspace depletion against linear
in valve time; the paper does not quantify it). Runs to 1.0 s (three residence times after the tail), dt 0.4 ms.

| Case | theta after the cycle, centre / min / max (|r| <= 95 mm) | non-uniformity (max-min)/min | t90 after valve opening, centre / worst | wafer impingement exposure (Pa s), min-max | exposure non-uniformity | TMA captured by the wafer | peak surface TMA (Pa) |
|---|---|---|---|---|---|---|---|
| 25 ms valve, 0.63 mg, s = 3.44e-3 | 0.9994 / 0.9994 / 0.9994 | 0.0007 % | 0.07 / 0.08 s | 0.671-0.673 | 0.28 % | 2.4 % of 5.26e18 (41 monolayer-equivalents delivered) | 5.6 |
| 60 ms, 1.5 x amount | 1.0000 / 1.0000 / 1.0000 | 0.00001 % | 0.06 / 0.06 s | 1.017-1.019 | 0.18 % | 1.6 % | 8.5 |
| 60 ms, 2.4 x amount | 1.0000 / 1.0000 / 1.0000 | 0.0000 % | 0.05 / 0.05 s | 1.639-1.640 | 0.11 % | 1.0 % | 13.7 |
| 25 ms, s = 5e-3 | 1.0000 / 1.0000 / 1.0000 | 0.00003 % | 0.06 / 0.07 s | 0.671-0.673 | 0.28 % | 2.4 % | 5.6 |
| Measured, 25 ms (2019 Fig. 10, digitised, / 1.00 A per cycle at 60 ms) | 0.996 / 0.975 / 0.998 | 2.37 % (2.40 % stated) | | | | | |
| Measured, 60 ms (stated) | 1.00 | 0.58 % | | | | | |
| Their CFD + kinetics, 25 ms / 60 ms | | 2.45 % / 0.83 % | | 0.55-0.60 (2018 Fig. 12) | 10.7 % | | 0.7-3.7 at 25 ms (2019 Fig. 7b) |

Agreement, honestly. (i) Exposure magnitude: the axisymmetric model puts 0.67 Pa s of TMA impingement on the wafer
for the 25 ms pulse against 0.55 to 0.60 Pa s in their 3D CFD, 12 to 22 % high, with the surface partial pressure
rising from 1.1 Pa (20 ms) through 3.4 Pa (30 ms) to a 5.6 Pa peak at 50 ms against their 0.7 to 3.7 Pa snapshot at
25 ms: the transport part of the chamber, how much of the 41 monolayer-equivalents injected reaches the wafer before
the pump takes it, is reproduced to the accuracy of the dimensions read off a figure. (ii) Saturation state: the
model says the 25 ms pulse saturates the wafer (theta 0.9994, 90 % at 70 ms) and the 60 ms pulse saturates it
completely; the measurement says 0.975 to 0.998 of the 60 ms value, i.e. saturated to within 2.5 %. Both agree that
at 300 C the 25 ms pulse is at the knee of the saturation curve and that 60 ms is past it (their Section 4.1.1:
"at 300 C the OH surface sites were almost totally covered even with the 25 ms TMA pulse"). (iii) Non-uniformity:
not reproduced. The model gives 0.0007 % (25 ms) and 0 % (60 ms) where the measurement gives 2.40 % and 0.58 % and
their 3D CFD 2.45 % and 0.83 %. Two reasons, both outside what an axisymmetric frozen-flow model can carry. First,
their exposure field is 10.7 % non-uniform because the pulse momentum through the off-axis side inlet drives a
recirculation above the wafer (2018 Fig. 10-12) and the loading-door purge sweeps one side; ours is 0.28 % non-uniform
because the single axial inlet and the diffusive gap (Peclet 0.1 over the 35 mm edge gap, 0.8 over the cone) smooth
everything. Second, with irreversible Langmuir kinetics a 10.7 % exposure variation at 3.2 E90 would move coverage by
only 0.1 to 0.5 %; their 2.4 % needs softer kinetics (their reversible adsorption with 14 % desorption and a
three-OH-per-two-TMA site balance), or a contribution from the 21 C temperature gradient across the wafer and the
incompletely purged H2O at the cold loading door that they document (2019 Section 4.3: the 60 ms non-uniformity of
0.58 % persists in a saturated TMA step, so at least that much is not a TMA-dose effect). By the 1/s arithmetic, an
effective sticking of 1.7e-3 would put the 25 ms pulse at the measured minimum (theta 0.975) and 3.0e-3 at the
measured maximum (0.998); their 3.4e-3 sits at the top of that bracket. The residual non-uniformity therefore
reflects the three-dimensional injector and door geometry and the kinetic softness, not the chamber-scale transport,
which the axisymmetric model gets right.

Verdict. The chamber model reproduces the published Fiji F200 case where the published case can be represented:
the TMA exposure delivered to the wafer per pulse agrees with their 3D CFD to 12-22 %, and the saturation state of
the 25 ms and 60 ms pulses agrees with the measured growth per cycle (saturated to within 2.5 % and fully saturated).
The 2.4 % measured non-uniformity and its drop to 0.58 % at 60 ms are not reproduced: they come from an off-axis
pulse recirculation, a cold loading door and softer-than-Langmuir kinetics, none of which an axisymmetric frozen-flow
model with irreversible Langmuir uptake contains. The blanket saturation exposure of our showerhead model, 0.146 Pa s
at the surface, sits inside the 0.05 to 1.3 Pa s spread of kinetic saturation exposures that seven independent
sources imply, and within a factor 1.8 of the three 300 C sources that give an explicit sticking coefficient or
cross-section (Vandalon, Gakis, Ylilammi).

## Figures

- figures/blanket_saturation_literature.png: TMA exposure for saturation against temperature. Open squares: exposures
  the experiments applied (feed or reactor pressure times time); filled circles: planar 90 % exposure implied by each
  source's sticking coefficient or cross-section; diamond: the Gakis 25 ms impingement exposure that gave >= 97 %
  saturation; crosses: the PillarHall fit in this repository; stars: this chamber model, feed exposure (1.03 Pa s,
  7,760 L) and surface exposure (0.146 Pa s, 1,090 L); grey lines: Langmuir E90 at N0 = 4e18 /m2 for s = 1e-3, 5e-3,
  1e-2. Right axis in Langmuir. Data figures/blanket_saturation_literature.csv (source, DOI, temperature, kind of
  exposure, Pa s, L, s, N0, s this model would need, note on what was read from a figure).
- figures/validation_gakis.png: axisymmetric Fiji F200 approximation against Gakis et al. Left: coverage after the
  cycle across the 200 mm wafer for the four cases, with the measured 25 ms GPC profile (2019 Fig. 10, digitised,
  divided by the 1.00 A/cycle of the 60 ms pulse; filled squares the side away from the loading door, open circles
  the loading-door side, +/- 0.003). Middle: coverage against time at the centre (solid) and edge (dashed). Right:
  TMA partial pressure at the wafer centre; the shaded band is the 9.6 Pa of Ar. Data figures/validation_gakis.csv.

## Proposed paper text (new subsection in paper/draft.md)

### 2.5 Chamber-scale checks
The blanket-wafer pulse of Section 3.1 can be compared with published saturation data once the exposure is defined
with care. The recipe sets a feed exposure, 3 Pa of TMA for 0.345 s, 1.03 Pa s or 7,800 Langmuir; the surface sees
the time integral of its own TMA partial pressure, 0.146 Pa s or 1,100 Langmuir, because the wafer draws the gas
down to 0.14 of the feed value on average while it saturates, and that surface exposure equals the Langmuir value
ln(10) N0/(s F1) at s = 5e-3, as it must. Seven independent sources give either the exposure they applied or a
sticking coefficient or reaction cross-section from which the planar saturating exposure follows (Figure,
blanket_saturation_literature.png). The kinetic values span 0.05 to 1.3 Pa s: 0.27 Pa s from the sum-frequency
cross-section of Vandalon and Kessels, 0.21 from the reactor-fitted sticking coefficient of Gakis et al., 0.16 from
the PillarHall fit of Ylilammi et al., 0.33 to 1.3 from Arts et al. at 220 C, 0.05 at room temperature in vacuum
(Soto and Tysoe), 0.08 from our own PillarHall fit. Our 0.146 Pa s sits in the middle of the 300 C cluster, within a
factor 1.8 of each of the three sources that state a sticking coefficient at that temperature; the applied
exposures (0.56 to 160 Pa s) are upper bounds that recipes choose to over-dose and say nothing sharper. Exposure scales
as 1/s in this regime, so the model would reproduce the individual sources at s between 2.7e-3 and 8.7e-3.

A second check uses the one published case that combines reactor-scale CFD with a measured thickness map in a
commercial tool: the Ultratech Fiji F200 of Gakis et al. (200 mm wafer, 9.6 Pa Ar, 180 sccm, 300 C, 0.63 mg of TMA
per 25 ms valve opening, measured growth per cycle 0.975 to 0.998 A/cycle along the diameter against 1.00 at 60 ms,
non-uniformity 2.4 % falling to 0.58 %). We rebuilt the chamber as an axisymmetric wedge with dimensions read off their
cross-section (the papers give none), a single axial inlet standing in for their three, the carrier flow frozen,
their computed pulse amount and shape, and their kinetic model reduced to an irreversible sticking coefficient
s = 0.004 x 0.86 = 3.4e-3 and N0 = 4.1e18 /m2 from their growth per cycle. The model delivers 0.67 Pa s of TMA
impingement to the wafer against 0.55 to 0.60 Pa s in their three-dimensional CFD and saturates the wafer with the
25 ms pulse (coverage 0.9994, 90 % at 70 ms) and completely with the 60 ms pulse, as the measurement does to within
2.5 %. It does not reproduce the 2.4 % non-uniformity: their exposure field is 10.7 % non-uniform because the pulse
momentum through an off-axis inlet drives a recirculation over the wafer and the loading-door purge sweeps one side,
ours is 0.3 % because an axisymmetric inlet into a diffusive gap smooths everything; and with irreversible Langmuir
kinetics even their exposure variation would move coverage by under 0.5 %, so the remainder is their reversible,
softer kinetics and the 21 C temperature gradient and cold loading door they document. The chamber-scale transport
that this paper's conclusions rest on, how much of the injected precursor the wafer takes before the pump does, is
reproduced to the accuracy of the dimensions; the within-wafer uniformity of a particular injector is not a claim of
this paper.

## CLAIMS rows (CLAIMS.md table format)

| Blanket saturating exposure, feed and surface | feed 1.03 Pa s (7,760 L) at the centre, 1.16 Pa s worst; surface integral 0.146 Pa s (1,090 L), equal to the Langmuir ln(10) N0/(s F1) = 0.142 Pa s at s = 5e-3 | figures/blanket_saturation_literature.csv | PYTHONPATH=src .venv/bin/python analysis/blanket_saturation_literature.py (reads reactor/showerhead_pulse_langmuir_p0.5) |
| Blanket surface exposure vs literature | inside the 0.05-1.3 Pa s spread of kinetic saturation exposures from seven sources; within a factor 1.8 of Vandalon 2019 (0.27), Gakis 2019 (0.21), Ylilammi 2018 (0.16); s needed to match each: 2.7e-3 to 8.7e-3 | figures/blanket_saturation_literature.csv, .png | same; literature: DOIs in the csv (values read from figures flagged in the note column) |
| Fiji F200 (Gakis et al.) axisymmetric reproduction: wafer TMA exposure per 25 ms pulse | 0.671-0.673 Pa s vs 0.55-0.60 Pa s in their 3D CFD (2018 Fig. 12, colour scale +/- 0.0002 mol/m2); exposure non-uniformity 0.28 % vs 10.7 % | figures/validation_gakis.csv | .venv/bin/python reactor/validation_gakis/make_gakis_cases.py --run; reactor/validation_gakis/postprocess_gakis.py; literature: 10.1016/j.cherd.2018.02.031 |
| Fiji F200 reproduction: saturation state | 25 ms: coverage 0.9994 (90 % at 0.07 s), 60 ms: 1.0000, vs measured GPC 0.975-0.998 of the 60 ms value (2019 Fig. 10, digitised +/- 0.003) and 1.00 at 60 ms | figures/validation_gakis.csv, .png | same; literature: 10.1016/j.ces.2018.09.037 |
| Fiji F200 reproduction: non-uniformity | 0.0007 % (25 ms) / 0 % (60 ms) vs measured 2.40 % / 0.58 % and their CFD 2.45 % / 0.83 %: not reproduced; off-axis recirculation, cold loading door and reversible kinetics are absent from the axisymmetric model | figures/validation_gakis.csv | same |
| Fiji F200 reproduction: chamber flow | Q = 0.0603 m3/s at 9.6 Pa, 573 K (their pump 57.1 L/s at 543 K); inlet 48 m/s, Re 4.3; inlet/outlet flux 0.0602 m3/s; 13,020 cells, max non-orthogonality 56 deg | reactor/validation_gakis/steady/log.simpleFoam, log.checkMesh | reactor/validation_gakis/make_gakis_cases.py --run |

## Run log (paper/assumptions.md style)

## Night 3, part 4: chamber-scale validation (reactor/validation_gakis, analysis/blanket_saturation_literature.py; figures/blanket_saturation_literature.png, .csv, figures/validation_gakis.png, .csv)
(a) Blanket saturation against published data. No new chamber run. From reactor/showerhead_pulse_langmuir_p0.5:
feed exposure p_feed t90 = 3.00 Pa x 0.345 s = 1.03 Pa s (7,760 L) at the centre, 1.16 Pa s at r = 149 mm; surface
exposure (integral of the face TMA partial pressure to theta = 0.9) 0.146 Pa s (1,090 L) at the centre and 0.146 at
the edge, equal to ln(10) N0/(s F1) = 0.142 Pa s (F1 = 1.30e22 /m2/s/Pa at 573 K). Literature exposures compiled with
their DOIs in the csv: applied exposures 0.56 (Gakis, impingement from their CFD), 7.5 (our PillarHall fit), 14.7
(Ylilammi, fitted entrance pressure), 25 and 160 Pa s (Arts, stated doses); kinetic E90 implied by stated s or sigma:
0.05 (Soto and Tysoe, 25 C, via Puurunen 2005), 0.08 (PillarHall fit, s = 8.5e-3), 0.16 (Ylilammi, gamma = 5.72e-3,
q = 5e18), 0.21 (Gakis, s0 = 0.004 x 0.86), 0.27 (Vandalon and Kessels, sigma = 6.5e-18 cm2, 300 C; 0.22 at 100 C),
0.33 to 1.3 (Arts, s0 = 2e-3 to 5e-4, 220 C); Dillon 1995 via Puurunen (k = 0.0076 /s at 0.01 Torr, 227 C) gives
404 Pa s on a porous membrane and is transport-limited, listed but excluded. Values read from figures: Gakis 2018
Fig. 12 colour scale (0.0118-0.0130 mol/m2, +/- 0.0002), Gakis 2019 Fig. 10 squares (digitised by marker centroids,
+/- 0.003 A/cycle; digitised non-uniformity 2.37 % vs 2.40 % stated), Gakis 2018 Fig. 8 pulse shape (+/- 20 % per
ordinate, integral renormalised). Not obtained: Ott 1997, Groner 2004, Elam 2002 full texts (not open; secondary
sources give pulse times without partial pressures). Reading: our surface exposure sits in the middle of the 300 C
cluster (0.08-0.27 Pa s), within a factor 1.8 of each explicit-sticking source; the feed exposure is 7x larger because
the blanket wafer draws the surface density down to 0.14 c0 while it saturates (the supply-limited chamber), and
literature feed exposures overstate surface exposures for the same reason. s needed to match each source: 2.7e-3
(Vandalon), 3.4e-3 (Gakis), 4.7e-3 (Ylilammi), 2.2e-3 to 5.5e-4 (Arts), 8.7e-3 (PillarHall fit); s = 5e-3 stays.
(b) Gakis et al. Fiji F200. Geometry measured on the 2018 Fig. 5a cross-section (1.65 px/mm on the 200 mm chuck top,
+/- 3 mm): inlet tube r = 20 mm, cone from r = 20 mm at +199 mm to r = 128 mm at +35 mm, wide chamber r = 128 mm to
-50 mm, chuck cylinder r = 100 mm to -50 mm (real dome 124 mm deep), duct r = 106 mm to -205 mm, outlet on its floor.
One axial inlet for their three (top 100, door 50, side 30 sccm), no loading-door box, carrier flow frozen, TMA a
passive scalar (its inlet partial pressure peaks at 11.5 Pa against 9.6 Pa of Ar: the real pulse is a slug of
nearly pure TMA; their outlet pressure rises 18 %), isothermal 573 K, D = 0.263 m2/s (N2 value scaled 1/p),
lambda = 1.4 mm with 2 mm first cells. Q = 0.0603 m3/s (180 sccm at 24 C, 1 bar -> 9.6 Pa, 573 K), U_in = 48 m/s,
Re 4.3, nu = 0.447 m2/s, residence 0.23 s. Mesh 13,020 cells, max non-orthogonality 56 deg, max skewness 2.5
(checkMesh OK). simpleFoam: the showerhead template's consistent SIMPLE at relaxation 0.9 diverged at iteration 98
(11 s); standard SIMPLE with one non-orthogonal corrector, p 0.3 / U 0.7, converged in 2000 iterations (Ux residual
3.8e-7, p 2.5e-8, ClockTime 77 s); inlet and outlet flux 0.0602 m3/s. Pulses: s = 3.44e-3 (their s0 = 0.004 times
the 86 % reaction probability of adsorbed TMA at 300 C), N0 = 4.1e18 /m2 (their 1.0 A/cycle at 3500 kg/m3), vbar
410 m/s; 25 ms amount 8.74e-6 mol (their 0.63 mg) with their Fig. 8 shape (peak at 12 ms, 5 % by 50 ms, tail to
0.3 s), 60 ms at 1.5x and 2.4x (not quantified in the paper), one case at s = 5e-3; Danckwerts inlet fraction
f = 48/(48 + 128) = 0.27; dt 0.4 ms to 1.0 s, writes every 10 ms; 104 to 117 s wall time each (1.7 to 1.9 min,
sharing the CPU with another solver; a first attempt at dt 0.2 ms to 1.5 s was stopped at 0.25 s after 3 min).
| case | theta after the cycle (centre / min / max over |r| <= 95 mm) | non-uniformity | t90 centre / worst | wafer impingement exposure | exposure non-uniformity | captured | peak surface TMA |
| 25 ms, s = 3.44e-3 | 0.9994 / 0.9994 / 0.9994 | 0.0007 % | 0.07 / 0.08 s | 0.671-0.673 Pa s | 0.28 % | 2.4 % | 5.6 Pa at 50 ms |
| 60 ms x1.5 | 1.0000 | 0.00001 % | 0.06 / 0.06 s | 1.017-1.019 | 0.18 % | 1.6 % | 8.5 Pa |
| 60 ms x2.4 | 1.0000 | 0 % | 0.05 / 0.05 s | 1.639-1.640 | 0.11 % | 1.0 % | 13.7 Pa |
| 25 ms, s = 5e-3 | 1.0000 | 0.00003 % | 0.06 / 0.07 s | 0.671-0.673 | 0.28 % | 2.4 % | 5.6 Pa |
| measured (2019 Fig. 10 / stated) | 25 ms: 0.996 / 0.975 / 0.998 of the 60 ms value; 60 ms: 1.00 | 2.40 % / 0.58 % | | their CFD 0.55-0.60 Pa s | 10.7 % | | 0.7-3.7 Pa at 25 ms (2019 Fig. 7b) |
Reading: exposure to the wafer per pulse agrees with their 3D CFD to 12-22 %; saturation state agrees (25 ms at the
knee, 60 ms past it); non-uniformity not reproduced (0.0007 % vs 2.4 %): their recirculating off-axis injection and
cold loading door give a 10.7 % exposure spread we cannot carry axisymmetrically, and irreversible Langmuir turns even
that into < 0.5 % at 3.2 E90, so their reversible kinetics, the 21 C wafer gradient and the H2O purge residue (the
0.58 % at 60 ms is a saturated-TMA value) carry the rest. By 1/s, s_eff = 1.7e-3 to 3.0e-3 would put the 25 ms pulse
at the measured 0.975 to 0.998. Caveats: dimensions from a figure; frozen flow during a pulse that is 50 % of the
gas inventory; passive-scalar TMA above the carrier pressure; isothermal; chuck dome truncated; 60 ms amount assumed.

## Sources

- Gakis, Vergnes, Scheid, Vahlas, Caussat, Boudouvis, Chem. Eng. Res. Des. 132, 795-811 (2018), doi
  10.1016/j.cherd.2018.02.031; author manuscript hal-01730196 / OATAO 19671 (fetched via web.archive.org). Values
  read from figures: Fig. 5a chamber cross-section (dimensions, +/- 3 mm on a 1.65 px/mm scale), Fig. 8 TMA pulse
  flow-rate and mole-fraction history (+/- 20 % per ordinate), Fig. 12 wafer TMA exposure colour scale 0.0118-0.0130
  mol/m2 (+/- 0.0002). Stated in text: Table 1 conditions, 0.63 mg TMA per 25 ms, 57.1 L/s pump, 142 mTorr side
  inlet, 10.7 % exposure non-uniformity, purge times 2.7 / 3.1 s.
- Gakis, Vergnes, Scheid, Vahlas, Boudouvis, Caussat, Chem. Eng. Sci. 195, 399-412 (2019), doi
  10.1016/j.ces.2018.09.037; author manuscript hal-02329135 / OATAO 21237. Read from figures: Fig. 10 squares
  (digitised, +/- 0.003 A/cycle), Fig. 3 (GPC 1.0 at 300 C, +/- 0.01), Fig. 7b surface pressure scale (0.005-0.028
  Torr), Fig. 9 map ranges. Stated in text: Table 1 experiments, Table 2 kinetic parameters, s0,TMA = 0.004, 86 %
  reaction probability at 300 C, non-uniformities 2.40 / 0.58 % measured and 2.45 / 0.83 % modelled, 21 C wafer
  gradient, OH density fit.
- Ylilammi, Ylivaara, Puurunen, J. Appl. Phys. 123, 205301 (2018), doi 10.1063/1.5028178 (Aalto open copy): Table I
  pA0 = 147 Pa, gamma = 0.00572, K = 219 /Pa; text: Picosun R-150, 300 Pa, 0.1 s pulses, 150 sccm per line; q = 5e18
  from their worked example (Fig. 5 caption).
- Arts, Vandalon, Puurunen, Utriainen, Gao, Kessels, Knoops, J. Vac. Sci. Technol. A 37, 030908 (2019), doi
  10.1116/1.5093620 (TU/e open copy): text 190 and 1200 mTorr s TMA doses, Table I s0 ranges and the Vandalon
  (4 +/- 1)e-3 at 300 C.
- Vandalon and Kessels, Langmuir 35, 10374 (2019), doi 10.1021/acs.langmuir.9b01600 (PMC6694409): sigma_TMA =
  6.5 +/- 0.6 e-18 cm2 steady growth, 2.0 e-18 initial, H2O 3.7-4.0 e-20 cm2; site density assumption 6e14 /cm2.
  Vandalon and Kessels, J. Vac. Sci. Technol. A 35, 05C313 (2017), doi 10.1116/1.4993597: TMA cross-section
  independent of temperature 100-300 C (as cited by Gakis 2019 and Arts 2019; abstract only).
- Puurunen, J. Appl. Phys. 97, 121301 (2005), doi 10.1063/1.1940727 (VTT open copy): Section on kinetics quoting
  Soto and Tysoe, J. Vac. Sci. Technol. A 9, 2686 (1991), reaction probability 0.01, and Dillon, Ott, Way, George,
  Surf. Sci. 322, 230 (1995), k_a = 0.0076 /s at 227 C and 0.01 Torr; Fig. 26 (Rahtu, Alaranta, Ritala, Langmuir 17,
  6506 (2001)) saturation at about 3 s pulse, no partial pressure.
- PillarHall record: Zenodo 3974438 (Ylivaara, Puurunen; CC BY 4.0), fit in this repository (coupling/fit_pillarhall.py).
- Not obtained (full text not open): Ott, Klaus, Johnson, George, Thin Solid Films 292, 135 (1997), doi
  10.1016/S0040-6090(96)08934-1; Groner, Fabreguette, Elam, George, Chem. Mater. 16, 639 (2004), doi 10.1021/cm0304546;
  Elam, Groner, George, Rev. Sci. Instrum. 73, 2981 (2002), doi 10.1063/1.1490410; Gakis PhD thesis (NTUA 2019,
  tel-04169507, 293 pages, behind the HAL bot check and not archived).

## Files created
analysis/blanket_saturation_literature.py; figures/blanket_saturation_literature.csv, .png;
reactor/validation_gakis/make_gakis_cases.py, postprocess_gakis.py, steady/ (mesh, 2000/, logs, runtime.txt),
pulse_25ms/, pulse_60ms_x1.5/, pulse_60ms_x2.4/, pulse_25ms_s5e-3/ (each with params.txt, runtime.txt, logs, 101
write times); figures/validation_gakis.csv, .png; this file. Nothing under docs/, CLAIMS.md, assumptions.md,
draft.md, main.tex or README.md was edited.
