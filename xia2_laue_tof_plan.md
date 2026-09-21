# xia2.laue_tof — neutron Laue time-of-flight processing pipeline

## Status (2026-09-21)

**Written and working** (all lint/format/mypy clean, verified on real data):

| file | what it is |
|---|---|
| `src/xia2/Modules/Laue_TOF/laue_tof.py` | phil scope, the `*Params` dataclasses with `from_phil`, `classify_input`, `setup_from_phil`, `run_xia2_laue_tof` |
| `src/xia2/Modules/Laue_TOF/data_integration.py` | the step functions: `bin_run`, `run_import`, `find_spots`, `index`, `refine`, `integrate`, `combine`, `export_shelx`, `process_exposure`, `run_data_integration` |
| `src/xia2/Modules/Laue_TOF/util.py` | `report_timing` (copied, not imported from SSX — nothing here may depend on `xia2.Modules.SSX`) |
| `src/xia2/cli/laue_tof.py` | the `xia2.laue_tof` command: `ArgumentParser` / `setup_logging` / `cleanup` / `write_citations`, mirroring `cli/ssx.py`. Registered in `setup.py` as `xia2.laue_tof=xia2.cli.laue_tof:run` |

Verified end to end **from the shell** on
`ppase_allconfigs/config1/config1_scipp_output_50bins.h5`: 2296 strong spots →
`fft3d` 1927/2296 indexed (83.9 %), RMSD 1.52 px → refined cell
`(105.9, 95.69, 113.94, 90, 98.06, 90) C 1 2 1` → 1359 integrated reflections
carrying `wavelength_cal` over 2.23–3.54 Å, with `Adding Lorentz correction` in
the integration log; 35 s for the run. Seeded indexing gives `1925/2296 (83.8 %)`
and returns the correct centred cell.

Also verified at **N=2** (the same file linked under two names, so two
orientations of identical data): orientation_1 indexes de novo, orientation_2 is
seeded (`known_orientation`), and the combine/export step writes
`scale/combined.{expt,refl}` with 2 experiments, 2 imagesets, unique identifiers
and 2703 reflections, then `scale/dials.hkl` with **one SHELX batch per
orientation** (1359 + 1344 reflections in batches 0 and 1).

**Next, in order:**

1. pointless and lawless (§ *Step 8/9*). The unmerged MTZ converter is written.
2. Integrate everything the geometry allows, not only the observed spots
   (§ *Resolution: observed, then calculated*).
3. Finish `tests/regression/test_laue_tof.py` (SXD, the converter unit test),
   then the docs.

**`laue_tof_mtz.py` is written** (2026-09-21) and wired in as the `unmerged_mtz`
workflow step, writing `scale/unmerged.mtz` from the combined data with gemmi:
H, K, L, `M/ISYM` = 1, `BATCH` in dataset 0; `I`/`SIGI`, `LAMBDA`, `ROT`,
`XDET`/`YDET` in dataset 1; one batch record per orientation with `LDTYPE = 3`,
the cell, `PHIRANGE = 1`, `BSCALE = 1`, the refined `UMAT` (column-major),
`SOURCE`/`S0` anti-parallel to the beam and `ALAMBD`/`DELAMB` from that batch's
wavelength range. `validate_unmerged_mtz` re-reads the file and fails loudly on
each of the silent failures in the Step 7 table; `check_wavelength_column` is
there for checking pointless output. Verified on the two-orientation NMX data:
2703 observations, 2 batches, `ALAMBD` 2.84 Å, `DELAMB` 0.65 Å.

**Both intensities go in one MTZ** (user, 2026-09-21): summation as `I`/`SIGI`
and profile-fitted as `IPR`/`SIGIPR`, so the file serves either choice and the
two can be compared after scaling. Consequences, all implemented:

- A reflection whose profile fit failed gets the MTZ missing-number flag in
  `IPR`/`SIGIPR`, never a zero, so a failed fit cannot be read as a measurement.
- A poor profile-fit fraction no longer stops later orientations from attempting
  the fit - it only decides which intensities the SHELX export uses - since the
  MTZ wants both. Only an outright abort still falls back to summation.
- The combine step drops `intensity.prf.*` from every table when any orientation
  lacks it, because tables with different columns cannot be concatenated and an
  MTZ column cannot cover only some observations. It says which orientation was
  missing.
- SHELX HKLF2 holds one intensity by construction, so `output.intensity` still
  chooses there.

**Integration is a two-axis decision, taken once per run** (2026-09-21):

- `dials.tof_integrate` always writes `intensity.sum.*`; a profile method adds
  `intensity.prf.*` **alongside** it. So a profile run that works gives both, and
  "fall back to summation" is only needed for an outright failure.
- Outright failure is real: `method=profile_1d_ibix` **aborts** on the NMX test
  data (both masks, exit 134) with
  `DIALS_ASSERT(A >= min_bounds[0] && A <= max_bounds[0])` at
  `tof_profile_1d_ibix.h:246`. `A` is hard-coded to 1.0 and `min_bounds[0]` is
  1.0, so the assert can only fail through `max_bounds[0] = 1e4 * intensity_max`
  being NaN — i.e. a shoebox with non-finite intensities. It kills the whole
  run, not one reflection, and the message names no reflection. **Second DIALS
  bug report candidate.** `integration.method` therefore defaults to
  `profile_1d_ibix` with `fallback_method=summation`, and the fallback is
  remembered for the remaining orientations instead of failing on each.
- A profile run finishing is **not** the method working: `profile_3d_gutmann`
  completed on the same data but fitted **13 of 1359** reflections (~4 min,
  against 7 s for summation), leaving an `intensity.prf` column covering 1 % of
  the data. So `integration.min_profile_fraction` (default 0.5) is checked
  after every successful profile run: below it, the summation intensities are
  used, the remaining orientations skip the profile fit, and
  `output.intensity=auto` exports `sum`. (Gutmann is not a priority - user,
  2026-09-21 - it is kept only as a phil choice.)
- The two foreground masks give **different data, not different quality
  metrics of the same data** (summation, Lorentz on, `nproc=8`, one NMX
  orientation):

  | mask | wall | CPU | integrated | mean I/σ | median I/σ | mean foreground px | mean I |
  |---|---|---|---|---|---|---|---|
  | `ellipse` | 7.0 s | 4.5 s | 1359 | 36.5 | 29.8 | 71.6 | 1.49e-3 |
  | `seed_skewness` | 26.5 s | 144 s | 1359 | 59.6 | 46.2 | 322.2 | 3.85e-3 |

  CC between the two sets of intensities is 0.95; `seed_skewness` takes in ~4.5×
  the foreground and 2.6× the intensity. Higher I/σ is **not** proof it is
  better — only merging statistics after lawless can settle that, so
  `mask=both` integrates each way on the **first** orientation, logs the
  comparison, keeps the winner by `mask_metric`, and **reuses that mask for
  every other orientation** so all batches are integrated alike.
- `dials.export format=shelx` refuses its own `intensity=auto` when both
  columns are present (*"Only 1 intensity option can be exported in this
  format"*), so the export step always passes an explicit `intensity=`, chosen
  from `output.intensity` (auto → profile when `intensity.prf.value` exists).

**Notes carried out of the combine step:**

- `dials.export` takes the Laue batch number straight from `imageset_id`, which
  is 0 in every per-orientation file, so `combine` renumbers it to the position
  of the experiment in the combined list (which is also where its imageset now
  sits). Without that every orientation exports as batch 0 — one of the silent
  failures listed in Risks. Batches are therefore **0-based** in the SHELX file;
  the unmerged-MTZ converter must emit **1-based** batch numbers for MTZ/pointless.
- `dials.tof_integrate` leaves **no record of the Lorentz correction** in its
  output: with `corrections.lorentz=True` it multiplies the intensities and
  variances by `L` in the C++ integrator and writes no column and no flag, so
  only the line *"Adding Lorentz correction"* in its log says it happened. The
  integrate step therefore stamps a `lorentz_applied` boolean column on
  `integrated.refl`; it survives the combine (concat) and `dials.export`, the
  combine step logs what it says, and orientations that disagree are rejected
  rather than scaled together. The MTZ converter should carry it into the file
  as well, and it is what drives the lawless `LORENTZ` keyword.
- Identical `image=` values are **deduplicated by phil**, and `FileInput.resolve_paths`
  resolves symlinks, so two links to one file end up as one input path (but still
  two entries, hence two orientations). Worth knowing when constructing tests.

**MANDI (ORNL) test case, added 2026-09-21** — `tests/regression/test_laue_tof.py`:

- Data: `MANDI_13378.nxs.h5` (CuZnSOD, IPTS-37773), 3.2 GB of raw events,
  copied off `/Volumes/Finke_NMX/Mandi/CuZnSOD/` into **`work/mandi/`**, which
  is gitignored (`/work` in `.gitignore`). Too big for `dials_data`, so the
  tests skip unless a copy is there or `XIA2_LAUE_TOF_MANDI` points at one.
  The known cell, from the user's own lawless work on this dataset, is
  **`P 65 2 2`, `65.73 65.73 151.88 90 90 120`**.
- **MANDI is binned in place, not by essnmx** (user, 2026-09-21): dxtbx writes
  the histogram into the event file itself - bin edges in
  `entry/time_of_flight`, counts in `entry/<bank>_events/spectra` - which is
  what `FormatMANDI` reads, and **the event data is kept**
  (`remove_event_data=False`). `FormatMANDI` has no reader for a separately
  reduced file yet; when it gets one, this goes back through `essmandi-reduce`.
  `binning.essmandi_reduce` was removed and a generic `binning.tof_padding`
  (100 us, the dxtbx default) added. **No instrument geometry is repeated in
  xia2**: the panel size comes from `FormatMANDI._get_image_size()`, which can
  be asked before the histogram exists, since the format instantiates on a raw
  event file. `histogram_mandi_run` skips the work if the file already carries a
  histogram, and one binning per file means `binning.integrate_bins` does not
  apply (warned about, `index_bins` wins).
- **Verified end to end on `MANDI_13378.nxs.h5`**: 50 bins of 337.4 us over
  14639-31507 us written into the event file (which grew 3.220 -> 3.274 GB, the
  events kept), then `dials.import` read it with **FormatMANDI, 40 panels, 50
  TOF bins**, and spotfinding found **227 strong spots**; 10 m 35 s, almost all
  of it the histogram. It did not index - not pursued yet; 227 spots off 40
  panels is few, so `spotfinding.min_spot_size` (15, tuned for NMX) and the bin
  count are the first things to try. With the histogram cached both tests run
  in 10 s.
- **Spotfinding and indexing parameters on MANDI** (all with the known cell,
  `indexing.max_cell` from it, `d_min_start=4`):

  | spot list | strong | fft3d | real_space_grid_search | fft1d |
  |---|---|---|---|---|
  | `min_spot_size=15` (NMX default) | 227 | fails | 147 (65 %) | 159 (70 %) |
  | `min_spot_size=3` | 370 | 205 (55 %) | 227 (61 %) | 238 (64 %) |
  | `min_spot_size=3` + `radial_profile` | 603 | 236 (39 %) | 282 (47 %) | 284 (47 %) |

  Every solution refines to the right cell (65.7-66.2, 151.9-152.9 against the
  known 65.73, 151.88). `min_spot_size=15` is **too big for this data**: of 558
  spots found, 331 are smaller than 15 pixels, so it throws away 59 % of them,
  and MANDI spots span few TOF frames at 337 us bins. It still indexes, just off
  a third of the spots. `sigma_strong` makes almost no difference (366 spots at
  2.0 against 370 at 3.0); `radial_profile` finds the most.
  The phil default stays 15, which is what the NMX study chose - say if it
  should change - but the step now **says so**: it reports how much of the spot
  list the size filter took and, when that is over half, there is room to lower
  it, and **fewer than 500 spots survive**, suggests `min_spot_size` and
  `radial_profile`. That last condition matters: thresholding turns up tens of
  thousands of one and two pixel candidates, so on NMX the filter discards 99%
  of 187300 and leaves 2296 spots, which is not a problem, where MANDI loses 59%
  of 558 and is left with 227.
- **Two bugs in the indexing acceptance test, both found here**:
  1. `min_indexed` was an absolute 250, so every good MANDI solution (147-284
     indexed) was rejected and the orientation counted as "not indexed". It is
     now a floor of 50 plus `min_indexed_fraction` (0.25 of the strong spots),
     which travels between datasets: 65 % of 227 passes, while 300 of 2300 NMX
     spots no longer would.
  2. There was **no test on the residuals**, so the first method to index enough
     reflections won. On MANDI that took fft3d at 3.31 px over
     real_space_grid_search at 3.25 px, and with `radial_profile` it accepted
     fft3d at **6.10 px**. Now `max_rmsd_px` (5.0) rejects a bad fit outright,
     `target_rmsd_px` (2.0) stops the ladder when a solution is good enough,
     and otherwise every method is tried and the **best by RMSD is kept** -
     each attempt writes `indexed_<method>.{expt,refl}` and the chosen one is
     copied to `indexed.{expt,refl}`.
     Note this RMSD is over every indexed reflection, so it reads higher than
     the `RMSD_X`/`RMSD_Y` in the dials.index log, which are post-outlier
     rejection: MANDI is 3.25 px here where dials.index reports 1.5/1.6 px.
  Verified: MANDI `min_spot_size=3` now walks fft3d (3.31) ->
  real_space_grid_search (3.25) -> fft1d (5.42, rejected) and keeps
  real_space_grid_search, 51 s. NMX is unchanged - fft3d 1927/2296 (83.9 %),
  RMSD 1.52 px, accepted at once, 26 s.
- **Two things in `FormatMANDI` had to be worked around**, both worth reporting:
  1. `get_time_range_for_panel` returns `event_time_zero + event_time_offset`,
     i.e. **seconds since the start of the run plus microseconds since the
     pulse**. For this 23 h run it gives 14700-113000 where every event
     actually arrives between 14739 and 31407 us, so histogramming over that
     range would silently pile every event into the first few bins. It is also
     a Python loop over all 4.8 M pulses of all 44 banks, minutes of work.
     `_mandi_tof_range` takes min/max of `event_time_offset` instead - the
     quantity `generate_histogram_data_for_panel` actually histograms - in
     seconds, or uses `binning.min_time_bin`/`max_time_bin` when given.
  2. `add_histogram_data_to_nxs_file` derives its bins from that range with a
     `delta_tof` bin *width*. The step calls `write_histogram_data` with
     explicit edges instead, so the bin *count* is the one asked for and the
     event data is not scanned twice.

  A third thing, not worked around: `generate_histogram_data_for_panel` builds a
  pool task per pixel (65536 of them) and hands each one the panel's whole event
  arrays, so a panel takes about a minute and a 40 panel run about 40 minutes.
  The same histogram is a few seconds of `np.bincount` on
  `pixel_index * nbins + bin_index`. Worth replacing, in dxtbx or here.
- `run_import` names the format class in its errors and logs panels/TOF bins on
  success: a reduced NMX file and a reduced MANDI file are both NXlauetof, and
  an instrument-specific class that claims the wrong one fails deep in dxtbx.
  (Seen with the essnmx route: `ess/nmx/nexus.py::_set_default_instrument`
  writes `entry/instrument/name = "NMX"` unconditionally, so a MANDI reduction
  is claimed by `FormatESSNMX`, whose `get_detector` does
  `i = int(panel_name[-1])` into a three-entry NMX dict and raises
  `KeyError: 2`. Fix that in essnmx before returning to that route.)

**Not started**: everything in *Wavelength normalisation and scaling*, which was
deliberately deferred — the module deals only with binning and the DIALS steps.

**Environment**: `.venv/` in this repo is the interpreter for all of it — see
[Environment](#environment-and-external-programs) below. `dials_known_orientation_centred_cell_issue.md`
in the repo root is a ready-to-file DIALS bug report, with
`dials_known_orientation_fix_check.py` as its verification script.

## Context

xia2 has no support for neutron **Laue time-of-flight (TOF)** data. All existing
pipelines assume monochromatic X-rays (a single scalar wavelength per sweep, e.g.
`XWavelength.set_wavelength`), and `load_imagesets` filters to rotation
`ImageSequence` only. Neutron TOF Laue is polychromatic: each reflection's
wavelength is derived from its arrival time, the beam is a `PolychromaticBeam`,
and intensities must be normalised by an incident spectrum (vanadium) and
corrected for background (empty run) and absorption.

DIALS 3.30 (the tested dev environment at `/Users/aaronfinke/dials/dials-v3-30-0`,
into which this xia2 fork is built as `modules/xia2`) already provides the complete
TOF toolchain — this task wires those tools into a single xia2 command so a user
can go from raw images to combined/exported intensities in one call. The final
**wavelength normalisation, scaling and merging** stage is done by two CCP4-side
programs built locally: **POINTLESS 1.16.1** (`/Users/aaronfinke/pointless`) and
**LAWLESS 0.0.1** (`/Users/aaronfinke/aimless`, the Laue branch of aimless).

**Experiment type**: a **single crystal measured as multiple stationary exposures
at different orientations** ("settings"). Each exposure is itself a full Laue-TOF
measurement — a stack of TOF bins imported as one `ImageSequence` — recorded with
the crystal held **stationary during the exposure**. This is the key contrast with
standard rotation crystallography, where the crystal rotates during each exposure.
The pipeline therefore processes each orientation independently (import →
find_spots → index → refine → tof_integrate) and then **combines the integrated
results across orientations**, writes one unmerged Laue MTZ, and passes it through
POINTLESS → LAWLESS to produce scaled/normalised merged and unmerged data. Primary
instrument **ISIS SXD** (the only one with public `dials_data`), but format-agnostic
so MaNDi / ESS NMX work too.

Where the input is **not already TOF-binned**, the pipeline first reduces it with
**essnmx** (`/Users/aaronfinke/ess/packages/essnmx`) — see
[Step 0](#step-0-new-tof-binning-with-essnmx).

The proven DIALS recipe (from the `isis_sxd_nacl_processed` dataset provenance)
is what the pipeline automates per exposure:
```
dials.import  sxd_nacl_run.nxs
dials.find_spots  imported.expt
dials.index  imported.expt strong.refl unit_cell=5.64,5.64,5.64,90,90,90 space_group=225
dials.refine  indexed.expt indexed.refl
dials.tof_integrate  refined.expt refined.refl [corrections.incident_run=vanadium] [corrections.empty_run=empty]
dials.export format=shelx  integrated.expt integrated.refl
```
with `dials.refine_bravais_settings` added after indexing (user request) for
lattice-symmetry determination, a combine step across all orientations, and then
the new scaling tail:
```
<xia2 converter>            combined.{expt,refl} -> unmerged.mtz   (LAMBDA, M/ISYM, batch headers)
pointless HKLIN unmerged.mtz HKLOUT sorted.mtz   (LAUEGROUP / CHOOSE SPACEGROUP given explicitly)
lawless   HKLIN sorted.mtz   HKLOUT scaled.mtz   (PROBE NEUTRON, LAUE NORMGPR ...)
```

## Decisions (confirmed with user)

- Command name: **`xia2.laue_tof`**
- Data: **one crystal, multiple stationary exposures at varying orientations**
  (crystal stationary during each exposure).
- Scope: **full pipeline**, per exposure import → find_spots → index →
  refine_bravais_settings → refine → tof_integrate, then **combine all orientations**,
  write an **unmerged Laue MTZ**, and finish with **POINTLESS (sort/symmetry) →
  LAWLESS (wavelength normalisation, scaling, merging)**. SHELX HKLF export of the
  unscaled combined data is kept as a checkpoint; the scaled/merged MTZ from lawless
  is the primary deliverable. No `dials.scale`/`dials.merge` — those are
  monochromatic-only; lawless replaces them (and LSCALE) for Laue.
- Vanadium (`incident_run`) + empty (`empty_run`) runs: **optional, recommended** —
  passed to `dials.tof_integrate` when given, corrections skipped when absent.
- Regression testing against **ISIS SXD** (`dials_data("isis_sxd_example_data")`).
- Include **`dials.refine_bravais_settings`** as a lattice-determination step.
- **Indexing defaults come from the NMX parameter study**, not from the DIALS/ToF
  auto-defaults, and indexing is a **retry ladder with an acceptance test**, not a
  single call — see [Indexing defaults](#indexing-defaults-from-the-nmx-study).
- **Formats**: `FormatISISSXD`, `FormatMANDI`, `FormatESSNMX`. The `*Laue`
  variants (`FormatMANDILaue`, `FormatESSNMXLaue`) are for TOF-flattened, true
  Laue data and are **out of scope** here.
- **Un-binned input is reduced in-pipeline by essnmx**, chosen on the input file's
  NeXus type (NXlauetof = ready; TOFRAW / NXsnsevent = bin first) — see
  [Step 0](#step-0-new-tof-binning-with-essnmx).
- **Local program builds, not the CCP4 9 ones**: pointless must be **> 1.13.6**
  (the CCP4-9 build drops the wavelength column, which makes `LAUE` fatal in
  lawless); lawless is not in CCP4 at all yet. Both are PHIL-configurable paths.

## DIALS building blocks (all present in 3.30, verified)

- **Import / formats**: `FormatISISSXD`, **`FormatMANDI`**, **`FormatESSNMX`**
  produce an `ImageSequence` (images = TOF bins) with `PolychromaticBeam` +
  `Spectrum` + `Scan`, `Probe = neutron`. `ExperimentList.all_tof()` /
  `all_laue()` identify such data.
  `FormatESSNMX_McStas` (in the same file, keyed on
  `entry/metadata/mcstas_weight2count_scale_factor`) is for the superseded McStas
  format and is **not** used — current McStas output is identical to real NMX data.
  **Not `FormatMANDILaue` / `FormatESSNMXLaue` either.** Those subclass the TOF
  formats but read the `laue_data` group — the data **flattened in the TOF
  dimension** — which is true (wavelength-integrated) Laue, a different processing
  problem. This pipeline is Laue-TOF: the TOF axis is the signal, so it must stay.
  The `understand()` of each `*Laue` class keys on `"laue_data" in entry/instrument/detector_panel_0`,
  so a file carrying both groups can be claimed by the wrong class — pass
  `format.name`/an explicit format class at import if that is ever seen, and assert
  `ExperimentList.all_tof()` on `imported.expt` before going further.
- **Spot finding**: `dials.find_spots` auto-selects `TOFSpotFinder` for TOF
  experiments (adds initial per-reflection wavelengths; TOF reciprocal-space
  proximity filter).
- **Indexing**: `dials.index` — works when unit cell + space group are supplied
  (TOF gives wavelength per spot), and also indexes many cases with no known cell.
  **The defaults matter a great deal and the DIALS/ToF auto-defaults are not the
  right ones** — see [Indexing defaults](#indexing-defaults-from-the-nmx-study)
  below. Expose `indexing.method` across `fft3d | fft1d | real_space_grid_search`
  (plus `pink_indexer | low_res_spot_match` as polychromatic fallbacks) and make it
  cheap to retry a failed orientation with a different method.
- **Bravais settings**: `dials.refine_bravais_settings` **works with polychromatic /
  TOF-Laue experiments in DIALS 3.30** (confirmed by user). xia2 already wraps it at
  `src/xia2/Wrappers/Dials/RefineBravaisSettings.py` (reads `bravais_summary.json`).
- **Integration**: `dials.tof_integrate` (`dials/algorithms/integration/tof`) —
  `method = summation | profile_1d_ibix | profile_3d_gutmann`,
  `integration_type = observed | calculated`,
  `background_model = constant2d | constant3d | linear2d | *linear3d`,
  `mask = ellipse | seed_skewness`, `bbox_tof_padding` / `bbox_xy_padding`,
  `wavelength_range`, `mp.nproc`, and phil **`corrections.incident_run` /
  `corrections.empty_run` / `corrections.lorentz` / `corrections.absorption.{incident,target}_spectrum{…}`**
  (incident-spectrum defaults tuned for SXD vanadium). Note the exact scope names —
  they are under `corrections{}`, not top level. Outputs carry `wavelength_cal`
  per reflection.
  **`corrections.lorentz` defaults to `False`** (`tof_integrate.py:104`): the TOF
  Lorentz factor is *not* applied unless asked for. This couples directly to the
  lawless `LORENTZ` keyword — see [Lorentz: exactly once](#lorentz-exactly-once).
- **Export**: `dials.export format=shelx` handles `all_laue()/all_tof()` — writes
  HKLF2 with per-reflection `wavelength_cal` and Laue batch numbers (batch =
  `imageset_id`/`id`, i.e. one per orientation; `dials/util/export_shelx.py`).
  **`format=mtz` is refused for TOF data** — `dials/command_line/export.py:362`
  raises *"Time-of-flight Laue data cannot be exported as .mtz. Run with
  format=shelx to export as .hkl"*. Even if it were not, `export_mtz.py` writes no
  per-reflection wavelength column, sets `M/ISYM = 0` and `LDTYPE = 2`. **xia2 must
  therefore write the unmerged Laue MTZ itself** (see the converter below); this is
  the "in-between step" the pointless/lawless tail needs.

## Step 0 (new): TOF binning with essnmx

`dials.import` needs a **TOF-binned** file: raw event data is not read by
`FormatESSNMX` or `FormatMANDI`, so it must be reduced first, by **essnmx**
(`/Users/aaronfinke/ess/packages/essnmx`, the `ess.nmx` / `ess.mandi` packages of
the scipp ESS monorepo). McStas simulations now export the same format as real NMX
data, so they take exactly this path — **`FormatESSNMX_McStas` and
`essnmx_reduce_mcstas` are for the old McStas format and are out of scope**.

### Deciding whether to bin

Classify the input file by its NeXus entry before anything else — this is the check
that decides whether Step 0 runs:

| file | marker | action |
|---|---|---|
| **NXlauetof** (already reduced/binned) | `entry.attrs["NX_class"] == "NXlauetof"`, `entry/definitions == "NXlauetof"`; panels carry `data` + `time_of_flight` | none — hand straight to `dials.import` |
| **NXTOFRAW** (raw NMX) | `entry/definition == "TOFRAW"`, or an `NXentry` whose `instrument` holds choppers / `instrument_xml` and event-style panels | `essnmx-reduce` |
| **NXsnsevent** (raw MANDI) | `entry` NX_class / `definition` of `NXsnsevent`; `instrument/name == "MANDI"` | `essmandi-reduce` |

Verified against real files: `ppase_allconfigs/cco/config1_binned{50,200}.h5` classify
as `nxlauetof`, `config1_sampling_output.h5` as `tofraw`. Implement this as one
small `classify_input(path)` helper in `laue_tof_programs.py` (h5py only, no scipp
import needed for the decision), and log the verdict per input file.

### Running the reduction

All three entry points are console scripts of `essnmx`, driven exactly like the
DIALS steps (`subprocess.run`):

```
essnmx-reduce   --input-file raw.h5             --nbins 50 --output-file binned50.h5
essmandi-reduce --input-file mandi_event.nxs.h5 --nbins 50 --output-file binned50.h5
```

Relevant options, all of which should be reachable from xia2 PHIL:
`--nbins` **or** `--time-bin-width` (mutually exclusive), `--min-time-bin` /
`--max-time-bin` / `--time-bin-unit`, `--detector-ids`, `--compression`,
`--output-dir`, `--overwrite`.

This is where the
[binning recommendation](#tof-binning-is-a-pipeline-parameter) is actually
implemented: when the input is raw, xia2 can produce **both** reductions from one
raw file — a coarse one to index against and a fine one to integrate — instead of
making the user pre-bin twice. PHIL:

```
binning {
  enabled = Auto            # Auto = bin only if the input is not NXlauetof
  index_bins = 50
  integrate_bins = 200
  time_bin_unit = *us ms ns
  min_time_bin = None  max_time_bin = None
  detector_ids = None
  essnmx_reduce = None      # executable; default: PATH lookup
  keywords = None           # extra raw CLI args appended
}
```

with `binning.enabled=False` to force "use the file as given", and the reduced files
written into the per-exposure subdirectory so a re-run does not redo them.

## Indexing defaults (from the NMX study)

Source: `/Users/aaronfinke/ppase_allconfigs/NMX_indexing_recommendations.html` and
`NMX_indexing_study.html` (PPase, C 1 2 1, 106.1 95.5 113.7 90 98.1 90; 4,172
indexing trials over eleven simulated ESS NMX detector layouts at two TOF
binnings, DIALS 3.30 + optuna), plus the second study on
`/Users/aaronfinke/ppase_allconfigs/cco` (CcO, P 21 21 21, 182 204 178). Per-run
recipes live in `configN/opt_NNbins/best/` and `best_by_balanced/` as
`find_spots.phil`, `index.phil`, `run.sh`, with every trial in `results.jsonl`.

The search took NMX datasets that autoindex from 4 of 22 to 15 of 22, and best
positional RMSD from 2.3 px to 0.59 px. The changes worth baking in as
`xia2.laue_tof` defaults:

| parameter | DIALS/ToF default | xia2.laue_tof default | why |
|---|---|---|---|
| `spotfinder.filter.min_spot_size` | 6 (ToF auto) | **15** (search 8–25) | the largest single improvement. 6 admits small split peaks; raising it drops positional RMSD from 2.3 px to 0.6–0.8 px on the same data, at some cost in reflection count. Balanced picks land at 14–24 in five of six near layouts |
| `spotfinder.threshold.algorithm` | `radial_profile` for ToF | **`dispersion_extended`**, keep `radial_profile` as first fallback | wins on four of six near PPase layouts. On the CcO data `radial_profile` wins more often — so this is a *default*, not a rule, and both must stay easy to select |
| `indexing.method` | `fft3d` | **`fft3d`** near, **`fft1d`** for distant/sparse layouts | `fft1d` is the finding of the study: the only method that indexed config8 and config10 (1 m panels). Its 1-D projection does not need a well-populated 3-D reciprocal-space grid, which a narrow recorded wedge cannot provide. `real_space_grid_search` is the other near-layout winner |
| `indexing.max_cell` | Auto | **fixed**, ≈ 1.2–1.4 × longest cell edge (120–150 Å for PPase, ≈ 230 Å for CcO) | nearest-neighbour estimation is unreliable on sparse Laue-TOF spot lists; *every* difficult-layout success used a fixed value |
| `refinement.reflections.outlier.algorithm` | `auto` | **`mcd`**, `tukey.iqr_multiplier` 0.6–1.0, `separate_panels=True` for multi-panel | the distant layouts carry a heavy residual tail (median 0.55 px, p99 26 px, p99/p50 = 47.6 against 8.9 for a good layout). A tight cut is what let config8 index at all |
| `indexing.refinement_protocol.d_min_start` | Auto | **≈ 4.0** | starting low-resolution helps the difficult layouts; immaterial for the near ones |

Design consequences for the pipeline:

- **Ship a retry ladder, not one recipe.** Near the indexing boundary success is a
  *rare event*, not a property of the data: config8 landed 2 acceptable solutions in
  212 trials, config11 4 in 190, config10 74 in 330. A single default recipe will
  fail on data that is perfectly indexable. `indexing` should therefore accept an
  ordered list of strategies and try them in turn — vary `method` (`fft3d` →
  `real_space_grid_search` → `fft1d`), then threshold algorithm, then
  `min_spot_size` — logging each attempt and stopping at the first solution that
  passes the acceptance test below. Expose the whole ladder as PHIL so a site can
  replace it, and let `indexing.strategy=<file.phil>` drop in a tuned
  `index.phil` from `configN/opt_NNbins/best/` verbatim.
- **Accept a solution on more than "it indexed".** The study's gate — refined cell
  matches within tolerance, orientation within tolerance of a reference where one
  exists, and a minimum indexed count — exists because maximising indexed *fraction*
  alone selects a small clean spot list that indexes almost completely (88.7 % of
  909 spots beating 72.1 % of 2,837), and because an a↔b swap can index 1,988
  reflections on the wrong cell. Rank on the geometric mean of indexed count,
  indexed fraction and RMSD quality, and log all three per orientation.
- **Use the other orientations as the reference.** In this experiment every exposure
  is the same crystal, so once one orientation indexes, its crystal model is the
  natural seed for the rest (`known_symmetry.A_matrix` / `dials.index` with a
  supplied orientation), and cross-agreement is a real check: six NMX
  configurations indexed independently agreed on the orientation to a median of
  0.04–0.09°. **Seeding also rescues layouts that never autoindex** — config7 and
  config9 produced nothing in 480 and 330 de-novo trials but process normally when
  given an orientation. This folds directly into the existing "consistent indexing
  across orientations" item in Risks: index the easiest exposure first, seed the
  rest, and treat a failure to autoindex as a reason to seed, not to stop.
- **Watch for the two DIALS defects the study surfaced** (both reproduce in 3.30,
  both worth reporting upstream):
  1. `known_symmetry.A_matrix` **corrupts a centred cell**. Supplying an orientation
     routes into `IndexerKnownOrientation`, which expects a primitive-setting
     crystal; the conversion loop in `known_orientation.py` assigns to the loop
     variable and discards the result, so it never happens, and
     `_apply_symmetry_post_indexing` then converts an already-centred cell again
     (C 1 2 1 `106.1 96.1 114.4 · 98.1` → `184.1 158.6 102.1 · 90`). **Convert in
     the caller**, which makes the missing conversion a no-op:
     ```python
     cb_op = crystal.get_space_group().info().change_of_basis_op_to_primitive_setting()
     A = crystal.change_basis(cb_op).get_A()
     ```
     xia2 must do this wherever it seeds an orientation.
  2. **ToF reflection predictor assert on large spot lists** — around 40,000
     reflections trips `DIALS_ASSERT(table.nrows() == h.size())` in
     `reflection_predictor.h` during refinement. Cap the spot list (or catch and
     retry with a stricter spot-finding recipe) rather than letting a long run die
     in refinement.

### TOF binning is a pipeline parameter

Where the input is a rebinnable TOF reduction (NMX/scipp, MaNDi), the number of TOF
bins is a *choice*, and the two studies disagree about which way it points — so it
must be a parameter, with the default stated and overridable:

- **PPase (small cell, 106–114 Å)**: every configuration has a higher indexable
  ceiling at 200 bins (config8 by 174 %, config11 by 497 %), and **eight of the nine
  comparable configurations autoindex *less* often at 200 bins** — config8 has the
  study's largest ceiling, 4,701 reflections, and zero de-novo solutions in 190
  trials. Recommended: **index at 50 bins, integrate at 200**, carrying the
  orientation across. config11 is the lone exception, indexing only at 200.
- **CcO (large cell, 182–204 Å)**: the balanced picks go the other way — 200 bins
  gives both the better indexed fraction (0.70–0.84 against 0.24–0.58) and far
  better orientation agreement (0.05–0.09° against ~2.5°). A larger cell packs
  reflections closer in TOF, so it needs the finer axis to separate them.

So: `input.tof_bins` (or a pre-binned file per stage) with
`indexing.tof_bins` / `integration.tof_bins` overridable independently, default
"index coarse, integrate fine", and a documented note that a layout failing at one
binning is worth retrying at the other before being written off. The mechanism —
finer binning resolves more genuine reflections while quartering counts per frame
admits proportionally more marginal ones, and the basis-vector search needs a peak
above the noise floor rather than more reflections — is a hypothesis consistent with
all 22 datasets, not a measurement.

## Wavelength normalisation and scaling (POINTLESS → LAWLESS)

Reference documentation: `/Users/aaronfinke/aimless/lawless.md` (program write-up)
and `/Users/aaronfinke/aimless/CLAUDE.md` (implementation notes, MTZ-input
pitfalls). Test invocation pattern: `/Users/aaronfinke/aimless/test_data/run_aimless.sh`.

### Why these two programs

`dials.tof_integrate` normalises each run by the measured incident spectrum, but
that leaves per-orientation scale, relative B, absorption/secondary-beam effects
and any **residual wavelength dependence** of the effective spectrum and detector
efficiency. Lawless fits exactly those: it is aimless plus a **wavelength
normalisation term `w(λ)`** in the scale model, fitted either as a Chebyshev
polynomial in `ln f(λ)` or as a Gaussian process, plus a `PROBE NEUTRON` switch
that sets neutron-appropriate defaults. Pointless is required first because
lawless needs the reflections **sorted and on a consistent, symmetry-assigned
setting**.

### Step 7 (new): combined integrated data → unmerged Laue MTZ

New module `src/xia2/Modules/LaueTOF/laue_tof_mtz.py`, writing with `gemmi` (already a
DIALS dependency; `export_mtz.py` uses it). Input: the combined
`integrated.{expt,refl}` from the combine step. Output: `unmerged.mtz`.

Requirements, each of which has been observed to fail silently (from
`aimless/CLAUDE.md`, "What the input MTZ must contain"):

| item | requirement | failure if wrong |
|---|---|---|
| `LAMBDA` column | per-observation wavelength in Å from `wavelength_cal`; label `LAMBDA` (lawless accepts `LAMBDA`/`LAM`/`WAVELENGTH`, case-insensitive, first in that order wins) | lawless stops when `LAUE` is given |
| `M/ISYM` | must be present and **written as 1 throughout** for unreduced indices | pointless reads the file as *merged*, collapses observations and silently drops extra columns |
| batch records | one per orientation, with `LDTYPE = 3` (Laue), a cell, non-zero `PHIRANGE`, `BSCALE = 1`, and `UMAT` = the refined crystal orientation | no batch records ⇒ pointless reports one batch regardless of the `BATCH` column; missing `UMAT` ⇒ `SCALES SECONDARY`/`ABSORPTION` refines flat **with no warning** |
| datasets | `HKL_base` as dataset 0 (H,K,L,M/ISYM,BATCH), data in dataset ≥ 1, batch `LBSETID` pointing at the data set | pointless dies with `Dataset::pxdname, setid not found` |
| `ROT` column | per-batch φ; optional | absent ⇒ lawless warns and substitutes the batch number (harmless for stationary exposures) |
| `BATCH` column | one batch number per orientation, matching `imageset_id`/`id` as used by the SHELX export | mis-scaled or merged-across-settings data |
| `ALAMBD`/`DELAMB` in batch header | band centre and half-width from the per-batch λ range | cosmetic/diagnostic |

Also write `I`, `SIGI` (summation or profile intensities as integrated), and
`XDET`/`YDET` where available. Do not write `SCALEUSED` — nothing has scaled the
data yet.

Working reference implementations of exactly this header-writing (gemmi, MTZ batch
word offsets already worked out): `/Users/aaronfinke/LADI_files/TIM_2020/aaron/addbatches.py`
(laue-dials → pointless) and `mtz1tomtz.py` (LAUEGEN `.mtz1`), plus
`/Volumes/Finke_NMX/Mandi/CuZnSOD/aaron/lauenorm2mtz.py`. Copy the offset table
(`I_LDTYPE=14`, `F_UMAT=6`, `F_PHISTT=36`, `F_PHIRANGE=47`, `F_BSCALE=43`,
`F_ALAMBD=86`, `F_DELAMB=87`, …) rather than re-deriving it.

### Step 8 (new): POINTLESS — sort and set symmetry

```
pointless HKLIN unmerged.mtz HKLOUT sorted.mtz XMLOUT pointless.xml <<EOF
LAUEGROUP P 6/m m m
CHOOSE SPACEGROUP P 65 2 2
EOF
```

- **Give the symmetry explicitly.** Pointless cannot reliably determine the Laue
  group from *unnormalised* Laue data: the wavelength dependence inflates the
  disagreement between symmetry mates, so even the identity operation scores badly
  and the ranking is unreliable (`lawless.md`, "A note on POINTLESS"). Default
  `lauegroup`/`space_group` in xia2 to what the user gave to `dials.index` (or to
  the `refine_bravais_settings` recommendation), and let both be overridden.
- **Version matters**: use the local **1.16.1** build (`/Users/aaronfinke/pointless/pointless`).
  1.13.6 (CCP4 9) does not pass the wavelength column through. 1.16.1 prints
  *"Additional unrecognised columns passed unchanged to output file"* — assert this
  in the log, or verify `LAMBDA` is present in `sorted.mtz` before calling lawless,
  and fail with a clear message naming the version requirement.
- `COPY` is available if only sorting/reindexing to a given setting is wanted
  (`keywords.cpp:1125`); expose it as `pointless.mode = choose | copy`.

### Step 9 (new): LAWLESS — wavelength normalisation, scaling, merging

```
lawless HKLIN sorted.mtz HKLOUT scaled.mtz XMLOUT lawless.xml <<EOF
PROBE NEUTRON
RUN 1 BATCH 1 TO <nbatches>
SCALES BATCH BFACTOR ON SECONDARY 0
LAUE NORMGPR <lam_min> <lam_max>
LAUE NORMLAMREF <lam_ref>
SDCORRECTION NOREFINE 1.0 0.0 0.0
REJECT 4
OUTPUT MERGED UNMERGED
EOF
```
(the recommended starting point for a neutron TOF dataset, from `lawless.md`
"Workflow for Laue data"). Defaults xia2 should set and log:

- `PROBE NEUTRON` — switches off automatic `ANOMALOUS`, zeroes the polarisation
  factor, makes `SCALES BATCH` the default primary mode (stationary exposures, so φ
  is not a scaling variable), rewords radiation-damage output and reports
  rejections against 2θ. Every change is printed by the program.
- `LAUE NORMGPR lam_min lam_max` — GP normalisation is the better default: it
  follows structure a low-order polynomial cannot and is **fixed after its
  pre-pass**, so it cannot trade against the other scale terms. `NORMCHEBYSHEV
  <degree> <lam_min> <lam_max>` is the alternative (refined jointly; up to 5
  non-overlapping ranges). They are **mutually exclusive — giving both is fatal**.
  Derive the default `lam_min`/`lam_max` from the `wavelength_cal` range of the
  combined data (or pass ≤ 0 to let lawless take them from the data) and
  `NORMLAMREF` from the intensity-weighted centre of the band.
- `NORMGPRBINS` (default 50), `NORMGPRLENGTH`, `NORMGPRMATERN` exposed as
  overrides. `NORMGPRPERRUN` is experimental — exposed but off.
- `SDCORRECTION ... SDLAMBDA [<scale>]` exposed for the case where χ² trends with
  wavelength; the exponent is fitted, not given, and is self-limiting.
- `LAMBDAONLY` exposed as a diagnostic mode (wavelength term only, no primary/B/
  secondary scales, no outlier rejection) — useful for isolating `w(λ)`.

Outputs to register with `FileHandler`: `scaled.mtz` (merged) and the unmerged
output (which now carries a `LAMBDA` column), `lawless.xml`, the log, and
**`LAMBDANORM`** — a self-contained gnuplot script of the fitted `w(λ)` curve
(`gnuplot -p LAMBDANORM`; GP fit carries a 1σ band, Chebyshev does not). Parse the
normalisation table and the merging statistics out of the XML for the xia2 summary
and report.

Optional final SHELX file from the merged MTZ via the existing
`src/xia2/Wrappers/CCP4/Mtz2various.py` wrapper (HKLF4), alongside the unscaled
HKLF2 checkpoint written by `dials.export format=shelx`.

### Lorentz: exactly once

**DIALS records nothing about it.** `corrections.lorentz=True` multiplies `I`,
`B` and their variances by `L = sin²θ/λ⁴` inside the C++ integrator
(`tof_integration.h`) and writes no column and no flag; the only trace is
*"Adding Lorentz correction"* in `tof_integrate.log`. An `integrated.refl` from
someone else therefore cannot be asked whether it is corrected. xia2 adds a
`lorentz_applied` boolean column at the integrate step so the answer travels
with the data through combine and into the MTZ, and refuses to combine
orientations that disagree.


`dials.tof_integrate corrections.lorentz` defaults to **False**; lawless `LORENTZ`
defaults to **NONE** and *nothing infers it* — not even `PROBE NEUTRON`. So by
default the Laue Lorentz factor is applied **nowhere**, and applying it twice is as
wrong as not at all. Model this as a single xia2 knob:

```
lorentz.applied_by = *integrate | lawless | none
```

- `integrate` (default) → `dials.tof_integrate corrections.lorentz=True`, lawless
  gets no `LORENTZ` keyword. Keeps the intermediate SHELX/MTZ physically correct
  whatever is done downstream.
- `lawless` → integration left at `corrections.lorentz=False`, lawless gets
  `LORENTZ TOF` (`I' = I·sin²θ/λ⁴`; `LORENTZ LAUE` is `I'=I·sin²θ`, differing only
  in where the λ⁴ is booked — the normalisation curve absorbs whichever is left in).
- `none` → neither; for data already corrected upstream.

Whichever is chosen, log the decision explicitly. Expect the overall R-merge to
move a lot when the correction is on (re-weighting of resolution shells), while
per-shell values stay put — that is not a loss of quality.

## Resolution: observed now, calculated next

`integration_type=observed` integrates only the spots that were found and
indexed, which on MANDI is a couple of hundred reflections of the tens of
thousands the geometry can reach. `integration_type=calculated` integrates every
predicted reflection out to `calculated.dmin`, and the question is what that
should be. **Design (user, 2026-09-21: multi-step, starting from the absolute
d_min of the detector geometry):**

1. **Pass 1, `observed`** as now: index, refine, integrate the found spots. This
   is what fixes the geometry, and it is cheap.
2. **Work out where to stop.** The hard floor is geometric -
   `d = lambda / (2 sin theta)` at the shortest wavelength and the largest
   scattering angle any pixel sees - and nothing beyond it can have been
   recorded. `geometric_d_min` computes it and the import step logs it
   (**implemented 2026-09-21**):

   | | wavelengths | 2theta max | geometric d_min | predicted reflections | observed integrated |
   |---|---|---|---|---|---|
   | NMX config1, 50 bins | 1.87-3.54 Å | 134.0° | **1.016 Å** | 93 304 | 1 359 |
   | MANDI CuZnSOD, 50 bins | 1.92-4.08 Å | 139.3° | **1.027 Å** | 82 630 | ~230 |

   The wavelength range comes from the scan's `time_of_flight` and the flight
   path, not from `beam.get_wavelength_range()`: FormatESSNMX derives that (and
   the two agree to three decimals) but **FormatMANDI hard-codes 2.0-4.0 Å**
   where the histogram gives 1.92-4.08, which alone moves MANDI's limit from
   1.066 to 1.027 Å. Prediction is cheap - `TOFReflectionPredictor` gives those
   counts in 1.6 s - so the number can be reported before committing to a pass.
3. **Pass 2, `calculated`** to that d_min, so that everything measurable is
   measured. Two to three orders of magnitude more reflections than pass 1, most
   of them weak, which is the point: the cutoff is a merging decision, not an
   integration one.
4. **Cut afterwards on CC½** from lawless, applied at merging. Re-running lawless
   with `RESOLUTION` costs nothing, where re-integrating costs everything, so the
   pipeline integrates wide once and cuts as often as it likes.

Open points for the implementation:

- Phil shape: `integration.expand = True/False` plus
  `integration.expand_d_min = *geometric observed <float>`, or fold it into
  `integration_type` with the d_min derived when it is `calculated`.
- A cost guard (`integration.max_predicted`): 93 000 shoeboxes is a different
  proposition from 1 359, and a large cell with an optimistic d_min could run for
  hours. Predict first, log the count, refuse above the guard.
- `integration.wavelength_range` should default to the TOF-derived band, so that
  predictions outside what was recorded are dropped rather than integrated as
  noise - particularly on MANDI, where the stated band is not the measured one.
- Whether to keep the pass 1 output as a checkpoint (it is a strict subset) or
  overwrite it.

## Files to add / modify

Follow the modern `xia2.ssx` pattern (functional/CLI orchestration that bypasses
the legacy `Schema/XWavelength` monochromatic model). Multiple exposures ⇒ loop the
per-exposure steps over the N orientations (lighter than SSX batching — each
orientation is one experiment, not thousands of stills).

1. **`src/xia2/cli/laue_tof.py`** (new) — entry point, mirror `src/xia2/cli/ssx.py`:
   `ArgumentParser(phil=phil_scope, read_experiments=False, read_reflections=False,
   check_format=False)`, `setup_logging(logfile="xia2.laue_tof.log", …)`, clear the
   dials logger handler, `cleanup(cwd)` context, then call
   `run_xia2_laue_tof(cwd, params)`. Citations (add LSCALE/Arzt 1999 and the
   aimless/pointless references used by lawless) + `write_citations`.

2. **`src/xia2/Modules/LaueTOF/__init__.py`** (new).

3. **`src/xia2/Modules/LaueTOF/xia2_laue_tof.py`** (new) — the PHIL scope
   (`full_phil_str`) and `run_xia2_laue_tof(cwd, params)` orchestrator. PHIL groups:
   ```
   input   { image=/template=/directory=  (repeatable → one entry per exposure)
             vanadium_run=None  empty_run=None }
   binning { … }                  # see Step 0; Auto = bin only un-binned input
   spotfinding { min_spot_size=15          # not the ToF auto-default of 6
                 threshold_algorithm = *dispersion_extended radial_profile dispersion
                 phil=None }
   indexing { unit_cell=None space_group=None min_spots=None
              method=None                  # None = walk the retry ladder
              ladder = fft3d real_space_grid_search fft1d   # ordered, retried in turn
              max_cell=None                # None -> 1.3 x longest cell edge, not Auto
              d_min_start=4.0
              outlier { algorithm=mcd  iqr_multiplier=0.8  separate_panels=True }
              seed_from_first=True         # seed later orientations from the first
                                           # success (primitive-setting A matrix)
              max_strong=35000             # guard the ToF predictor assert
              strategy=None                # drop-in index.phil, wins over all above
              phil=None }
   bravais_settings { enabled=True phil=None }
   integration { method=profile_1d_ibix   # summation is written whatever is set
                 fallback_method=summation  # profile fitting can abort the run
                 min_profile_fraction=0.5   # or fit almost nothing and "succeed"
                 integration_type=observed
                 background_model=linear3d
                 mask=both                  # ellipse | seed_skewness | both
                 mask_metric=i_over_sigma   # how both is decided
                 lorentz=True
                 bbox_tof_padding= bbox_xy_padding= wavelength_range=
                 absorption{…} phil=None }
   absorption { … passthrough / corrections.absorption overrides (default off) }
   lorentz { applied_by = *integrate lawless none }
   scaling {
     enabled = True
     pointless {
       executable = None          # default: $CCP4/bin/pointless, else PATH;
                                  # user points this at /Users/aaronfinke/pointless/pointless
       mode = *choose copy
       lauegroup = None           # default from indexing.space_group / bravais result
       space_group = None
       keywords = None            # extra raw keyword lines appended
     }
     lawless {
       executable = None          # default: PATH 'lawless', then ~/aimless/build/aimless
       probe = *neutron xray
       normalisation = *gpr chebyshev none
       lam_min = None  lam_max = None  lam_ref = None   # default from wavelength_cal
       chebyshev_degree = 6
       gpr { bins = 50  length = None  matern = False  per_run = False }
       scales = "BATCH BFACTOR ON SECONDARY 0"
       sdcorrection = "NOREFINE 1.0 0.0 0.0"
       sdlambda = None
       reject = 4
       resolution = None
       lambda_only = False        # LAMBDAONLY diagnostic mode
       keywords = None            # extra raw keyword lines appended
     }
   }
   output  { shelx=True mtz=True intensity=auto composition=CH }
   workflow.steps = find_spots+index+bravais+refine+integrate+combine+export+
                    unmerged_mtz+pointless+lawless
   nproc=<auto>
   ```
   Orchestrator loops the per-exposure steps below over each orientation (in a
   per-exposure subdirectory), then runs a single combine + export + the scaling
   tail in a `scale/` subdirectory. Decorated with `@report_timing` (reuse
   `xia2.Modules.SSX.util.report_timing`).

4. **`src/xia2/Modules/LaueTOF/laue_tof_programs.py`** (new) — one function per
   step, run **per exposure** (loop over the N orientations). **Recommended: drive
   the validated DIALS CLI via `subprocess.run`** (exactly the pattern of
   `xia2.Modules.SSX.data_integration_standard.run_import`), which guarantees the
   TOF code paths are taken and is trivial to verify against the known-good recipe.
   Per-exposure steps:
   - `classify_input` → h5py probe of the NeXus entry: `nxlauetof` (ready) |
     `tofraw` | `nxsnsevent` | `unknown`. See
     [Step 0](#step-0-new-tof-binning-with-essnmx).
   - `bin_run` (only when the input is not `nxlauetof`, or when two binnings are
     requested) → `essnmx-reduce` (NMX) / `essmandi-reduce` (MANDI) with
     `--nbins`, writing the coarse (index) and fine (integrate) reductions
   - `import_run` → `dials.import <image> output.experiments=imported.expt`
   - `find_spots` → `dials.find_spots imported.expt`
   - `index` → `dials.index imported.expt strong.refl [unit_cell=][space_group=]
     [indexing.method=] [max_cell=] [refinement_protocol.d_min_start=]
     [refinement.reflections.outlier.*]`, driven by the retry ladder and acceptance
     test of [Indexing defaults](#indexing-defaults-from-the-nmx-study). After the
     first success, later exposures are seeded with that crystal's **primitive-setting**
     A matrix (`change_of_basis_op_to_primitive_setting` in the caller). Log every
     attempt with its method, indexed count, indexed fraction and RMSD
   - `refine_bravais_settings` (if `bravais_settings.enabled`) →
     `dials.refine_bravais_settings indexed.expt indexed.refl`, then parse
     `bravais_summary.json` and log the scoring table (reuse the parsing approach
     from the existing `RefineBravaisSettings` wrapper). Informational for a known
     space group; report table + recommended setting.
   - `refine` → `dials.refine indexed.expt indexed.refl` → `refined.{expt,refl}`
   - `tof_integrate` → `dials.tof_integrate refined.expt refined.refl
     [corrections.incident_run=<vanadium>] [corrections.empty_run=<empty>]
     [corrections.lorentz=True] method= background_model=` → `integrated.{expt,refl}`
   Then, across all orientations:
   - `combine` — concatenate the per-exposure `integrated.{expt,refl}` (experiments
     + reflection tables, assigning unique identifiers/batches per orientation —
     reuse `assign_unique_identifiers` / reflection-table concat as in
     multiplex/SSX). Resolve indexing consistency across orientations first
     (see Risks).
   - `export` → `dials.export format=shelx <combined>` (multi-batch HKLF2 with
     wavelengths). **Not** `format=mtz` — refused for TOF.
   - `write_unmerged_mtz` → `laue_tof_mtz.py` (below) → `unmerged.mtz`.
   - `run_pointless` → `pointless HKLIN unmerged.mtz HKLOUT sorted.mtz XMLOUT
     pointless.xml`, keywords on stdin (`subprocess.run(..., input=keywords)`),
     `LAUEGROUP`/`CHOOSE SPACEGROUP` always written explicitly. Verify afterwards
     that `sorted.mtz` still has the wavelength column; if not, fail with the
     "needs pointless > 1.13.6" message.
   - `run_lawless` → `lawless HKLIN sorted.mtz HKLOUT scaled.mtz XMLOUT lawless.xml`,
     keywords on stdin, built from the `scaling.lawless` PHIL. Parse `lawless.xml`
     for merging statistics and the `<WavelengthNormalisation>` block; register
     `scaled.mtz`, the unmerged output, `LAMBDANORM`, `ROGUES`, the log and the XML
     with `FileHandler`.
   Each function raises a clear `ValueError` on non-zero return / stderr, records
   outputs via `xia2.Handlers.Files.FileHandler`, and per-step user phil files are
   merged as `dials.<prog> <user.phil> …` (SSX precedence pattern). For the two
   CCP4 programs the equivalent escape hatch is `scaling.<prog>.keywords`, appended
   after the generated keyword lines so the user always wins.

5. **`src/xia2/Modules/LaueTOF/laue_tof_mtz.py`** (new) — the combined-integrated →
   unmerged Laue MTZ converter described in Step 7, plus a small
   `check_wavelength_column(mtz_path)` helper used to validate pointless output.
   Keep it importable and unit-testable on its own; it is the piece most likely to
   need instrument-specific tweaks.

6. **`setup.py`** — add console script:
   `"xia2.laue_tof=xia2.cli.laue_tof:run"` in `console_scripts`.

7. **`tests/regression/test_laue_tof.py`** (new) — `dials_data("isis_sxd_example_data")`
   (`sxd_nacl_run.nxs`, `sxd_vanadium_run.nxs`, `sxd_empty_run.nxs`). The public SXD
   data is a single NaCl orientation, so this exercises the pipeline at N=1. Run
   `xia2.laue_tof` via subprocess with
   `unit_cell=5.64,5.64,5.64,90,90,90 space_group=225` + vanadium/empty runs; assert
   `integrated.{expt,refl}` and the SHELX `.hkl` exist and reflection count is
   sane. Model on `tests/regression/test_ssx.py`. Mark network/slow.
   Scaling tail: gate the pointless/lawless assertions on the executables being
   found (`pytest.mark.skipif`), since neither is a DIALS dependency. Add a
   **converter-only unit test** that does not need them: build `unmerged.mtz` from
   the integrated output and assert `LAMBDA` present, `M/ISYM == 1`, one batch
   record per orientation with `LDTYPE == 3`, non-zero `PHIRANGE`, a non-identity
   `UMAT`, and base columns in dataset 0. (Add a multi-orientation test later
   if/when multi-setting SXD data is available.)

8. **Docs** — `doc/sphinx/laue_tof/index.rst` + `basic_usage.rst`, added to the main
   `doc/sphinx/index.rst` toctree. Model on `doc/sphinx/multiplex/` and
   `serial_crystallography.rst`. Document the multi-exposure / stationary-crystal
   model, the vanadium/empty-run inputs, the pointless/lawless requirement (with
   the version caveat and how to point xia2 at local builds), the Lorentz
   single-application rule, and how to read `LAMBDANORM`.

## Reuse (do not re-implement)

- `dials.util.options.ArgumentParser`, `iotbx.phil` — CLI/PHIL (as in `cli/ssx.py`).
- `xia2.Handlers.Streams.setup_logging`, `xia2.Handlers.Files.{FileHandler,cleanup}`,
  `xia2.Handlers.Citations`, `xia2.Applications.xia2_main.write_citations`.
- `xia2.Modules.SSX.util.report_timing` and the `subprocess.run` orchestration
  pattern in `xia2.Modules.SSX.data_integration_standard.run_import`.
- Multi-experiment combine helpers (`assign_unique_identifiers`, reflection-table
  concat) as used in `xia2.Modules.MultiCrystal` / SSX data reduction.
- `src/xia2/Wrappers/Dials/RefineBravaisSettings.py` — reference for running
  `dials.refine_bravais_settings` and parsing `bravais_summary.json`.
- `src/xia2/Wrappers/CCP4/Pointless.py` and `Aimless.py` — legacy Driver-based
  wrappers. **Reference, not reuse**: they assume monochromatic rotation data and
  the CCP4 `Driver` machinery. Take from them the keyword vocabulary, the XML
  parsing and the log-scraping idioms; drive the new programs with `subprocess.run`
  like the rest of this pipeline. `src/xia2/Wrappers/CCP4/Mtz2various.py` **is**
  directly reusable for the final merged-MTZ → SHELX HKLF4 conversion.
- MTZ batch-header writing: `/Users/aaronfinke/LADI_files/TIM_2020/aaron/addbatches.py`
  and `mtz1tomtz.py`, `/Volumes/Finke_NMX/Mandi/CuZnSOD/aaron/lauenorm2mtz.py` —
  known-good gemmi code with the batch word offsets already correct.
- Indexing parameter study: `/Users/aaronfinke/ppase_allconfigs/` — the reports
  (`NMX_indexing_recommendations.html`, `NMX_indexing_study.html`), the harness
  (`nmxopt/nmx_index_opt.py`, `residuals.py`, `panels.py`, `README.md`), the
  per-configuration recipes (`configN/opt_NNbins/best*/{find_spots,index}.phil`,
  `run.sh`) and every trial with its metrics and rejection reason
  (`results.jsonl`). `cco/` is the same harness on a large-cell P 21 21 21 case.
  Take the phil files as the defaults' provenance and as regression fixtures; do
  not re-derive them.
- TOF reduction: **essnmx** at `/Users/aaronfinke/ess/packages/essnmx`
  (`ess.nmx`, `ess.mandi`; console scripts `essnmx-reduce`,
  `essnmx_reduce_mcstas`, `essmandi-reduce`). xia2 shells out to these, exactly as
  it does to the DIALS and CCP4 programs — no scipp workflow is imported into xia2.
- Program documentation: `/Users/aaronfinke/aimless/lawless.md` (keywords, input
  requirements, recommended workflow), `/Users/aaronfinke/aimless/CLAUDE.md`
  (implementation notes, MTZ pitfalls, Lorentz discussion),
  `/Users/aaronfinke/aimless/test_data/run_aimless.sh` (invocation pattern).
- All heavy lifting stays in DIALS / pointless / lawless; xia2 only orchestrates +
  logs + handles files.

## Risks / open items

- **Consistent indexing across orientations**: each exposure is indexed
  independently, so alternative indexings / indexing ambiguities must be resolved
  before combining (otherwise equivalents merge incorrectly downstream). Options:
  index/refine jointly against a shared crystal model, or resolve ambiguities with
  `dials.cosym` / reindex to a common setting (the SSX/multiplex approach) before
  the combine step. Decide during implementation once multi-orientation data is on
  hand. The NMX study makes seeding the obvious first move: once one exposure
  indexes, use its crystal model for the rest (converting to the primitive setting
  in the caller — see the `known_symmetry.A_matrix` defect), which both enforces a
  common setting and rescues exposures that never autoindex. Note that **pointless
  cannot be relied on to catch a wrong choice for Laue data** — see the
  symmetry-determination point below.
- **The unmerged-MTZ converter is the main new risk surface.** Every requirement in
  the Step 7 table has been observed to fail *silently* (merged-file collapse,
  single-batch reporting, flat secondary correction). Validate the file before
  running pointless — column labels, `M/ISYM`, batch count, `LDTYPE`, `PHIRANGE`,
  `UMAT` non-identity, dataset layout — and fail loudly with the specific cause.
- **Pointless cannot determine the Laue group from unnormalised Laue data.** The
  wavelength dependence inflates disagreement between symmetry mates, so the
  ranking is unreliable and even the identity operation scores poorly. Always pass
  `LAUEGROUP` + `CHOOSE SPACEGROUP`; treat any pointless symmetry *suggestion* as
  informational only, and never let it silently reindex.
- **Program availability and versions**: pointless must be **> 1.13.6** (local
  1.16.1 at `/Users/aaronfinke/pointless`); lawless (0.0.1, `/Users/aaronfinke/aimless/build`)
  is pre-release and not in CCP4, so its keywords may still move. Both need the
  CCP4 environment sourced (`/Applications/ccp4-9/bin/ccp4.setup-sh`) for libraries
  and symmetry data. Detect both up front, report versions in the log, and degrade
  gracefully — if either is missing, stop cleanly after `unmerged.mtz` with a
  message saying what to install and how to point `scaling.*.executable` at it,
  rather than failing at the end of a long run.
- **Double Lorentz / double normalisation**: see
  [Lorentz: exactly once](#lorentz-exactly-once). Similarly, `dials.tof_integrate`
  already divides by the measured vanadium spectrum, so lawless's `w(λ)` fits only
  the *residual* wavelength dependence — expect a much flatter curve than for a
  LADI/LAUEGEN dataset, and be suspicious of a curve with structure finer than the
  physics justifies (usually too many GP training bins). Open question to settle on
  real data: whether it is better to skip the vanadium correction at integration
  and let lawless fit the whole spectrum. Both paths should remain available.
- **Environment**: must run against the tested **DIALS 3.30** build at
  `/Users/aaronfinke/dials/dials-v3-30-0` (currently `v3.dev-1508-gb78e3bffe`); this
  is also the build the NMX indexing study was made with. Its `modules/xia2` is
  upstream xia2 on `main`, **not** this fork — install the `xia2N` checkout editable
  into that environment (or point `modules/xia2` at it) before testing. The
  Phenix-2.1 bundled DIALS lacks `dials.tof_integrate` and the TOF formats — do not
  use it. The repo's `.venv` (DIALS conda python + `--system-site-packages`) carries
  the fork, DIALS 3.30 and **essnmx/essreduce editable plus the scipp stack**, so
  all three toolchains are importable and on PATH from one interpreter.
- **essnmx availability**: reduction pulls in the scipp stack (scipp, scippnexus,
  sciline, essreduce, tof, plopp, dask), which is a much heavier dependency than
  anything else xia2 needs. Keep it **optional and lazily invoked** — a missing
  `essnmx-reduce` must only fail runs whose input actually needs binning, with a
  message saying which file and why, never a run given NXlauetof input.
- **Indexing without a known cell**: works for many Laue-TOF cases in DIALS 3.30;
  keep `indexing.method` (e.g. `pink_indexer`) exposed as a fallback for hard cases.
  Supplying `unit_cell`+`space_group` remains the most robust path — and note that
  every difficult-layout success in the NMX study also needed a **fixed
  `max_cell`**, which without a known cell has to come from the user.
- **Some geometries never autoindex.** config7 and config9 produced nothing in 480
  and 330 de-novo trials, including transfers of the recipes that cracked their
  near-identical sibling — a property of what those layouts record, not of the
  tuning. The pipeline must therefore treat "did not index" as a routine outcome:
  seed from another exposure, say so in the log and the report, and carry on rather
  than aborting the run.
- `refine_bravais_settings` stays a toggleable step (`bravais_settings.enabled`)
  even though it works with TOF in 3.30 — informational when the space group is
  already known.
- **One batch per orientation** means the scale model has few parameters per batch;
  with only a handful of settings, `SCALES BATCH BFACTOR ON SECONDARY 0` may still
  be over-parameterised. Keep `scaling.lawless.scales` easily overridable and
  report the parameter/observation ratio.

## Environment and external programs

Everything runs from `.venv/` in this repo (gitignored), built from the DIALS
conda python with `--system-site-packages`:

```
source .venv/bin/activate
export PATH=/Users/aaronfinke/dials/dials-v3-30-0/conda_base/bin:$PATH
```

The second line matters: `.venv/bin` holds the `xia2.*` scripts but **not** the
DIALS dispatchers, and the steps shell out to `dials.import` etc., so without the
conda `bin` on `$PATH` the first step fails with *"Unable to find dials.import on
$PATH"*.

| what | where |
|---|---|
| DIALS 3.30 | `/Users/aaronfinke/dials/dials-v3-30-0` (`v3.dev-1508-gb78e3bffe`) |
| essnmx / essreduce | `/Users/aaronfinke/ess/packages/` (editable installs; scipp stack from PyPI, numpy pinned to the DIALS 2.4.6 so the compiled extensions keep working) |
| xia2 (this fork) | editable install, so `xia2.*` console scripts run this tree, not the DIALS tree's `modules/xia2` |
| pointless 1.16.1 | `/Users/aaronfinke/pointless/pointless` |
| lawless 0.0.1 | `/Users/aaronfinke/aimless/build/aimless` |
| NMX indexing study | `/Users/aaronfinke/ppase_allconfigs/` (+ `cco/` for the large-cell case) |

Two things in the venv are hand-made and will not survive recreating it:
`.venv/share/cctbx` is a symlink to the DIALS `conda_base/share/cctbx` (without
it `libtbx.load_env` raises `FileNotFoundError`), and `.venv/bin/pytest` is a
shim that runs this venv's interpreter (otherwise a `pytest` from another conda
env on `$PATH` wins). `ruff` is not installed in the venv; the one used so far is
`/Users/aaronfinke/miniforge3/envs/nmx_strategy/bin/ruff`.

## Verification

Run everything in the DIALS 3.30 conda env
(`/Users/aaronfinke/dials/dials-v3-30-0`), with xia2N installed editable into it,
and with the CCP4 environment sourced for pointless/lawless.

1. **Smoke**: `dials.tof_integrate -h` and `xia2.laue_tof -h` both succeed;
   `python -c "import xia2.cli.laue_tof"`; `/Users/aaronfinke/pointless/pointless -i`
   reports **1.16.1**; `/Users/aaronfinke/aimless/build/aimless` banner reports
   **LAWLESS 0.0.1**.
2. **Binning / input classification**: `classify_input` returns `nxlauetof` for
   `ppase_allconfigs/cco/config1_binned{50,200}.h5` and `tofraw` for
   `config1_sampling_output.h5` (both confirmed); `essnmx-reduce --help` and
   `essmandi-reduce --help` both run in the project venv. Reducing a raw NMX file at `--nbins 50` produces a file that
   `dials.import` accepts and on which `ExperimentList.all_tof()` is true.
3. **End-to-end on real SXD data** (single orientation, N=1):
   ```
   dials.data get isis_sxd_example_data
   xia2.laue_tof image=<…>/sxd_nacl_run.nxs \
       unit_cell=5.64,5.64,5.64,90,90,90 space_group=225 \
       vanadium_run=<…>/sxd_vanadium_run.nxs empty_run=<…>/sxd_empty_run.nxs \
       scaling.pointless.executable=/Users/aaronfinke/pointless/pointless \
       scaling.lawless.executable=/Users/aaronfinke/aimless/build/aimless
   ```
   Expect: `imported.expt` → `strong.refl` → `indexed.*` → `bravais_summary.json`
   (logged table) → `refined.*` → `integrated.{expt,refl}` → combined SHELX `.hkl`
   (with a wavelength column) → `unmerged.mtz` → `sorted.mtz` → `scaled.mtz` +
   `LAMBDANORM`, plus `xia2.laue_tof.log`. Confirm `ExperimentList.all_tof()` is
   true on `imported.expt` and integrated reflections carry `wavelength_cal`.
4. **Converter checks on `unmerged.mtz`** (gemmi, before pointless): `LAMBDA`
   present with a sane Å range; `M/ISYM` all 1; one batch record per orientation
   with `LDTYPE == 3`, non-zero `PHIRANGE`, `BSCALE == 1` and a `UMAT` matching the
   refined crystal orientation; H,K,L,M/ISYM,BATCH in dataset 0 and I,SIGI,LAMBDA in
   dataset 1 with `LBSETID` pointing at it.
5. **Pointless passthrough**: `sorted.mtz` still has `LAMBDA`, and the log contains
   *"Additional unrecognised columns passed unchanged to output file"*. Observation
   count must match `unmerged.mtz` (a large drop means it was read as merged).
6. **Lawless**: runs without the *"no wavelength column"* stop; the log reports
   *"Per-reflection wavelength taken from column LAMBDA"*, the `PROBE NEUTRON`
   change block, and a wavelength-normalisation table + ASCII plot. Inspect
   `LAMBDANORM` (`gnuplot -p LAMBDANORM`) — with vanadium normalisation already
   applied at integration the curve should be close to flat. Check χ² against
   resolution, intensity and batch; compare R-merge/CC½ with and without
   `LAUE NORM…` (run `LAMBDAONLY` for an isolated view of `w(λ)`).
7. **Lorentz applied exactly once**: run with `lorentz.applied_by=integrate` and
   `=lawless` and confirm merged intensities agree up to an overall wavelength-
   dependent scale that the normalisation curve absorbs (per-shell statistics should
   barely move; overall R-merge may move a lot — that is re-weighting, see above).
8. **Indexing defaults**: on an NMX dataset (`~/ppase_allconfigs/configN/`),
   confirm the pipeline's default recipe reproduces the study's balanced pick to
   within noise — indexed count, indexed fraction and positional RMSD in the same
   band as `configN/opt_50bins/best_by_balanced/best.json` — and that the retry
   ladder recovers a layout (config8/config10) that only `fft1d` indexes. Confirm a
   seeded C-centred crystal comes back with the *input* cell, not a doubled one.
9. **Multi-orientation** (when multi-setting data is available): pass several
   `image=` exposures; confirm each is processed in its own subdirectory, the
   combined export carries one batch per orientation, and lawless's `RUN 1 BATCH 1
   TO N` covers them all.
10. **Regression**: `pytest tests/regression/test_laue_tof.py`.
11. Sanity-check output intensities/resolution against the DIALS SXD NaCl reference.
