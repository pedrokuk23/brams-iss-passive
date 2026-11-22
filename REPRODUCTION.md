# Independent reproduction (9-10 Nov 2026)

**Scope.** This file records a historical reproduction, made on 9-10 Nov 2026 on the version of the analysis of that date. The table below refers to that version. It is followed by the state of the current version.

The reproduction showed that the search outputs (`p2_*` and `blind_*` shards, `blind_summary.json`) and the ISS detections are reproducible from the raw recordings. At that date the feasibility model predicted 1,306 detections; the ISS timing rows refer to the first interpolation of the NASA ephemeris (cubic Hermite); the blind-search table gave the wider of a bootstrap over recordings (`ci_cluster`) and one over satellites (`ci_object`).

The whole pipeline was rerun from an empty `results/` folder on a separate Linux machine (8 CPU cores, Python 3.13,
numpy 2.5, scipy 1.18, sgp4 2.27, astropy 8.0, matplotlib 3.11). The inputs were the code and auxiliary files of
this deposit and the 1,121 BRAMS recordings, which were downloaded again from the BRAMS archive. The SHA-256 checksums
of all three folders were verified (`raw_wav` 293, `raw_wav_envisat` 40, `raw_wav_phase2` 788). The new results were
compared with the published ones by `code/compare_results.py`, which matches records by their text fields and compares
numbers by relative difference, and by `code/key_numbers.py`.

## Headline numbers

| Quantity | Published | Reproduced |
|---|---|---|
| Phase 1 blind search (`blind_summary.json`, 3,188 values) | — | all values identical |
| Phase 2 per-receiver thresholds z_th (BEOTTI, BEOUDS, FRMAUB, BEHUMA, BEUCCE) | 5.443, 6.226, 6.854, 5.316, 5.346 | identical |
| Phase 2 detections | 1,751 | 1,751 |
| Phase 2 detected passes, ISS complex / DTC / V2-mini / AST / V1.x | 174/419, 438/3,035, 1,028/12,031, 46/318, 61/5,700 | 174/419, 438/3,035, 1,028/12,033, 46/318, 61/5,700 |
| Feasibility test: observed, predicted, 95% predictive interval | 1,592, 1,305.9, [1,129, 1,481] | 1,592, 1,305.8, [1,129, 1,481] |
| Passes with M < 1 detected, Phase 1 / Phase 2 | 7/6,401, 0/78 | 7/6,401, 0/78 |
| Mantel-Haenszel odds ratio shadow vs sunlit, V2-mini Phase 1 / Phase 2 | 47.58 / 16.41 | 47.58 / 16.41 |
| Same, DTC Phase 1 / Phase 2 | 4.92 / 3.76 | 4.92 / 3.77 |
| Hold-out false alarms of the ISS search | 1/602 | 1/602 |
| ISS detection tables (Table 3, without/with line removal) | — | all values within 1e-3 (relative) |
| BlueBird search and nulls | — | all values within 1e-3 |
| Ionospheric fits (`iono_fit_v2.json`) | — | all 142 values within 1e-2 (e.g. an offset of 7.956 vs 8.010 ms) |
| Along-track comparison with the NASA ephemeris (`oem_cmp_v2.json`) | — | 125 of 135 values within 1e-2; the largest relative difference is 7%, on a near-zero distance (-5.6 vs -5.2 m) |

The fresh Phase 2 pass list has two more V2-mini passes than the published one (12,033 vs 12,031); neither was
detected. (An earlier version of this table summed the per-receiver counts together with the class total and so
doubled every count; `code/key_numbers.py` now reads the class total only.) The small numerical differences are most likely round-off from different numpy and scipy versions.

## Differences that do not change the paper

- **Tobit antenna-map fit (`antenna_map_tobit.json`).** All elevation coefficients are shifted by about -0.12, a
  constant that trades off against the free offset of the model. The coefficient of the highest elevation bin, which
  has only 13 blocks, moved from 7.41 to 9.30. The bootstrap of the antenna map (`antenna_map_bootstrap*.json`) uses
  random resampling and differs in the last digits.
- **Synthetic echoes at Ottignies (`synth_BEOTTI.json`).** 11 of 3,818 records are matched differently; the summary
  differs in the last digits.
- **Catalog scan (`catalog_scan.json`).** 296 of 14,657 records change their "is the ISS" flag. That flag compares the
  peak with that of the ISS to 1e-6, so floating-point round-off changes it.
- **Phase 1 shard files (`blind_shard*.json`).** These hold the same searches distributed differently among the
  parallel shards. Their summary is identical.
- **Files the published folder had and a fresh run does not produce.** These are the intermediate files of earlier
  runs (`p2_dev*`, `p2_cloud*`, `blind_v3_blankonly/`, `oem_cmp.json`, `iono_fit_final.json`, ...). The fresh run
  writes its Phase 2 shards as `p2_rerun*.json`.

## Fixes this reproduction required (all in this version of `code_and_results.zip`)

1. `run_all.sh`: the steps now run in dependency order, so the script works from an empty `results/` folder.
   - The Doppler-residual tracks are computed before Table 3.
   - `iss_table.py` was added.
   - The first catalog scan now runs before the close-pass search.
   - `detect_theory.py` now runs before `close_echo.py`.
   - `air_mf.py` now runs before the figures.
2. `phase2_analyze.py` reads the WAV headers when `results/phase2_meta.json` is absent.
3. `cat3_eval.py` writes `catalog_scan.json` as the list of records that `close_scan.py` and `glint_close.py` read.
4. Four figure scripts (`fig_starlink.py`, `fig_detectability.py`, `fig_feasibility.py`, `fig_blind.py`) wrote to the
   manuscript folder outside the repository. They now write to `results/figures/`, or to `$FIGDIR`.
5. `antmap_run.py` and `air_uniform.py` took the file time from the full path. A folder name with an underscore broke
   them. They now use the file name only.
6. `bb2_run.py`, `bb2_null.py` and `bb_null_clean.py` used every recording in the folder. They now read the fixed lists
   in `bluebird_files.json`. With the full set of recordings, the search had stopped with an all-NaN error on one of
   the additional files.
7. `null_run.py` and `null_v2.py` used every recording in the folder. With all recordings present, they absorbed the
   hold-out sample, and `null_holdout.json` came out empty. They now read the fixed lists in `null_files.json`.
8. `air_mf.py`, which writes the input of Fig. 2, was missing from the package. It now writes
   `results/air_mf_rows.npy`.

The published results were produced before these fixes. They are unchanged, and the table above shows that the fixed
code reproduces them.

## Current version (22 Nov 2026)

The analysis was revised after the historical reproduction: echo energy of the detectability model at the exact search frequency, beacon power of 340 W, one transmitter position for the timing and echo-energy analyses with receiver offsets relative to Oudsbergen, a 10-point Lagrange interpolation of the full archived NASA files, intervals from a two-way bootstrap over satellites and time clusters (over physical passes for the ISS complex), and a joint admissible-offset test of the catalog association. The revised feasibility model predicts about 1,552 detections (`results/feasibility_validation.json`). The searches and their selection keep the configuration of the original runs (README, Configurations): the transmitter position written in each file and the planning sensitivity model (`data/planning/planning_model.json`).

Checks of the current version, made from the published deposit (raw recordings checked against `SHA256SUMS`):

| Check | Result |
|---|---|
| `repro_critical.sh` from the 293 recordings (time base, ISS detections, tracks, NASA comparison, joint fit, offsets, hold-out nulls) | `pps_audit.json`, `iss_detection_table_v2.json`, `iono_tracks_v2.json`, `oem_cmp_v2.json`, `iono_fit_v2.json`, `iss_checks2.json`, `null_holdout.json` identical to the published files |
| Reanalysis of the archived search outputs (`blind_analyze`, `phase2_analyze`, `phase2_compare`, `eclipse`, `feasibility_validate`, `stats_cluster`, `catalog_association`, `coincidence`, `maubeuge_separate`, `shadow_glm`) | identical to the published files |
| Phase 2 search (`phase2_search.py`, historical configuration) on four recordings (306 real and null searches) | identical to the published shards (z, tau, nu) |
| Phase 2 metadata (`phase2_meta.py`) | identical to `results/phase2_meta.json` |
| Phase 2 planning values with the frozen planning model | the 24,284 values of `sigma_min` of the frozen selection (`data/planning/phase2_byfile.json`) are reproduced exactly from their ranges and arcs |
