# Opportunistic multistatic sensing of objects in low Earth orbit with an open meteor radar network

Code and data for the paper *Opportunistic multistatic sensing of objects in low Earth orbit with an open meteor radar network* (P. Kukulka de Albuquerque, submitted to Advances in Space Research).

Repository: https://github.com/pedrokuk23/brams-iss-passive. Archived version with the raw recordings: Zenodo, doi:10.5281/zenodo.24050975.

The raw WAV files are not in this repository (they total about 4.1 GB); download them from the BRAMS portal or from the Zenodo deposit as described under *Setup*. Some large intermediate results that the code regenerates are also left out (see `.gitignore`).

The BRAMS network (Belgium) records a continuous 49.97 MHz beacon at many GPS-timed receivers and publishes the raw audio. This repository uses those recordings, without new hardware, to:

1. calibrate the timing of the network with aircraft tracked by ADS-B, and map the relative gain of each receiver with a full day of aircraft echoes;
2. remove narrowband interference lines from the SDR receivers and detect the International Space Station (ISS) with a matched filter built from its element set (TLE), at up to five receivers;
3. measure the along-track error of the TLE from the residual Doppler and compare it with the NASA ephemeris (predictions propagated less than one day agree to 36 ms rms, about 270 m, over 17 passes of two campaigns; the pass-to-pass pattern is reproduced by an SGP4 fit to the NASA ephemeris);
4. separate a small ionospheric term and compare it with the Dourbes ionosonde and with GNSS TEC maps;
5. scan the catalog and detect a BlueBird (AST SpaceMobile) satellite;
6. search every satellite blindly with one frozen procedure (`PROTOCOL.md`): Phase 1 on the recordings of 3-5 Nov 2026 (development) and Phase 2 on 788 new recordings of 19 Oct - 2 Nov 2026 (test), with detection probabilities per class against sign-reversed nulls;
7. decide in advance whether a target is worth observing (`code/brams_feasibility.py`), a rule fitted on Phase 1 and checked for consistency on Phase 2.

## Layout

```
code/      analysis scripts (run them from inside code/)
data/      inputs (see DATA_LICENSES.md); raw BRAMS WAV files are not included
results/   intermediate and final products used in the paper, and the figures
```

## Setup

```bash
pip install -r requirements.txt      # or requirements-lock.txt for the exact versions of the published run
```

Python 3.13 (tested; exact versions with their dependencies in `requirements-lock.txt`).

**Raw data.** Download the BRAMS WAV files listed in `data/raw_wav/FILES.txt` from the BRAMS data portal (https://brams.aeronomie.be, dataset doi:10.18758/x620ti37) and place them in `data/raw_wav/`. The same 293 files, unmodified, are also in the Zenodo deposit of this paper (doi:10.5281/zenodo.24050975), split into five archives (`brams_wav_20261103.zip`, `brams_wav_20261104_part1.zip` to `_part3.zip` and `brams_wav_20261105.zip`); extract them with `for z in brams_wav_2026110*.zip; do unzip -j -o "$z" -d data/raw_wav/; done`; the Envisat archive `brams_wav_envisat_20261103-05.zip` goes to `data/raw_wav_envisat/` and the six `brams_wav_phase2_part*.zip` to `data/raw_wav_phase2/` in the same way. The Phase 2 element sets are in `data/tle/phase2/` (also in `auxiliary_data_phase2.zip` on Zenodo). File names follow `RAD_BEDOUR_YYYYMMDD_HHMM_<station>_SYS001.wav`.

**Element sets.** All element sets used are in `data/tle/` (Space-Track.org, redistributed under its blanket approval for basic SSA data, except the catalog snapshot `leo.txt`, taken from CelesTrak; see `DATA_LICENSES.md`), and `data/tle/QUERIES.md` gives the exact Space-Track queries that retrieve them.

## Running everything

`run_all.sh` runs the whole pipeline in order once the raw WAV files are in `data/raw_wav/` (verify them with `sha256sum -c SHA256SUMS` inside that folder). `python -m pytest -q tests` checks the core formulas and the published numbers, including independent checks of the Earth-fixed transformation and the station coordinates against astropy, of the numerical Doppler against the analytic expression from the SGP4 velocity, of the recovery of an injected timing error and of the Chapman integral; when raw files are present it also requires all 293 of them and verifies their checksums. The 40 Envisat recordings (`data/raw_wav_envisat/FILES.txt`, `SHA256SUMS`) are in the Zenodo deposit as a separate archive; place them in `data/raw_wav_envisat/`. The 788 Phase 2 recordings (`data/raw_wav_phase2/FILES.txt`, `SHA256SUMS`) are in the deposit as six `brams_wav_phase2_part*.zip` archives; place them in `data/raw_wav_phase2/`. `repro_critical.sh` recomputes the ISS detections, the along-track offsets and the hold-out false alarms from the raw files; in a clean copy these reproduce the published JSON files exactly.

**Independent reproduction.** An earlier version of the pipeline was rerun from an empty `results/` folder on a separate machine, with all 1,121 recordings downloaded again from the BRAMS archive (checksums verified); it reproduced the search outputs (blind-search shards and ISS detections) within round-off. The later revisions of the analysis (including the echo energy at the exact search frequency, one transmitter position for every receiver and the two-way bootstrap intervals) were recomputed from those search outputs and from the raw recordings on the main machine. `REPRODUCTION.md` gives the exact scope, the outcome and the fixes that run required. `code/compare_results.py <published results> <new results>` compares two result folders record by record, and `code/key_numbers.py` prints the headline numbers of both side by side.

The null distributions, the BlueBird search and its nulls read fixed file lists (`code/null_files.json`, `code/bluebird_files.json`), so they do not depend on which other recordings are in `data/raw_wav/`. The manuscript figures are written to `results/figures/`, or to the folder given in the environment variable `FIGDIR`.

## Pipeline

All scripts are run from `code/`. The results of every step are already in `results/`, so any step can be rerun on its own.

| Step | Script | Output | Paper |
|---|---|---|---|
| WAV reader (BRAMS RIFF chunks `BRA1`, `data`, `BRA2`) | `bwav.py` | — | Sec. 2 |
| Common geometry, SGP4 propagation (TEME→ECEF with GMST), baseband conversion, direct-signal suppression | `station_run.py` | — | Secs. 2–3 |
| ADS-B extraction from the adsb.lol daily archive | `adsb_extract.py` | `data/adsb/*.jsonl` | Sec. 4 |
| Aircraft matched filter, uniform parameters for all stations and days | `run_uniform.sh` → `air_uniform.py` | `results/uni_<station>_<day>.json` | Sec. 4, Supplementary S4 |
| Aircraft statistics per receiver | `aircraft_stats_uniform.py` | `results/aircraft_uniform.json` | Supplementary Table S1 |
| Histogram of the time offset τ | `fig_tau_hist.py` | `results/figures/fig_tau_hist.*` | Supplementary Fig. S1 |
| Matched filter of every ADS-B aircraft near the link (Ottignies, 3 Nov 10:00-10:25 UTC) | `air_mf.py` | `results/air_mf_rows.npy` | Fig. 2 |
| Example spectrograms with ADS-B overlay | `paper_figs3.py` | `fig_aircraft_example.png` | Fig. 2 |
| Spread of the aircraft offsets vs. ADS-B sampling and SNR | `adsb_lat.py` | `results/adsb_latency.json` | Sec. 4, Supplementary S4 |
| Joint (τ, Δf) search for the ISS | `iss2d.py <wav files>` | printed | Sec. 5 |
| Geometry map and ISS detection figure | `paper_figs2.py` | `fig_geometry.pdf`, `fig_iss_detection.png`, `data/stations/coords.json` | Figs. 1, 5 |
| Time-resolved residual Doppler of the ISS echo | `iono_track.py <wav files>` | `results/iono_tracks.json` | Sec. 6 |
| Least-squares fit: TLE timing + Chapman ionosphere | `iono_fit_final.py` | `results/iono_fit_final.json` | Table 4 (first version) |
| Residual figure | `paper_figs4.py` | `fig_residual.pdf` | Fig. 6 (first version) |
| Comparison with the NASA ephemeris (OEM, 10-point Lagrange interpolation of the full archived files, GCRS→ITRS) | `oem_cmp.py` | `results/oem_cmp.json` (first version: nine Ottignies/Oudsbergen tracks) | Sec. 6 |
| Sensitivity to hmF2 and to the ionosonde values | `iono_sens.py` | printed | Sec. 7 |
| Receiver differences: aircraft elevation, sensitivity and ISS geometry | `ant_geom.py` | `results/aircraft_elev.json` | Sec. 11.2 |
| Time-base audit of every file (GPS PPS / block time stamps) | `pps_audit.py` | `results/pps_audit.json` | Sec. 2.3 |
| Time-base correction of the Oudsbergen offsets | applied from `pps_audit.json` | `results/dt_pps_corrected.json` | Table 4 |
| ISS detection statistic for every pass and receiver | `iss_table.py` | `results/iss_detection_table.json` | Table 3 |
| Null distribution: ISS search on recordings without the ISS (fixed file list `code/null_files.json`) | `null_run.py` | `results/null_iss.json` | Sec. 5.2 |
| Synthetic echoes injected into real recordings (`synth.py` holds the shared functions) | `synth_run.py` (env `ST=BEOTTI` or `BEOUDS`) | `results/synth_<station>.json` | Sec. 5.2 |
| Detection probability, false alarm and estimator error; figure | `synth_eval.py` | `results/synthetic_summary.json`, `fig_synthetic.*` | Fig. 5, Sec. 5.2 |
| BlueBird (AST SpaceMobile) search, wider grid, and its null distribution (need `data/tle/leo.txt`; fixed file lists `code/bluebird_files.json`) | `bb2_run.py`, `bb2_null.py` (functions in `bb2.py`) | `results/bb2_results.json`, `results/bb2_null.json`, `results/bluebird_0310.json` | Sec. 8 |
| Dechirped spectrograms of the BlueBird echo | `bb_fig.py` | `fig_bluebird.*` | Fig. 7 |
| Catalog scan of all LEO objects (needs `data/tle/leo.txt`) | `run_cat3.sh` → `catscan3.py`; `cat3_eval.py` | `results/catscan/`, `results/catalog_scan.json` | Supplementary S1 |
| **Version 2 (ASR manuscript)** | | | |
| Narrowband line detection and removal (SDR receivers) | `linefilt.py` | — | Sec. 3.1 |
| ISS detection statistic without/with line removal | `iss_table_v2.py` | `results/iss_detection_table_v2.json` | Table 3 |
| Time-resolved residual with line removal (`CLEAN=1`) | `iono_track.py` | `results/iono_tracks_clean_sdr.json`, `results/iono_tracks_v2.json` | Sec. 6 |
| NASA ephemeris comparison for all receivers (`TRACKS=...`) | `oem_cmp.py` | `results/oem_cmp_v2.json` | Sec. 6.1, Table 4 |
| Fit with time-base correction, all receivers | `iono_fit_v2.py` | `results/iono_fit_v2.json` | Table 4 |
| Sensitivity to hmF2 and ionosonde values, all receivers | `iono_sens_v2.py` | `results/iono_sens_v2.json` | Sec. 7 |
| Residual figure, all receivers | `paper_figs4_v2.py` | `fig_residual.*` | Fig. 6 |
| Null distribution with line removal (fixed file list `code/null_files.json`) | `null_v2.py` | `results/null_v2.json` | Sec. 5.2 |
| Synthetic echoes at Humain with line removal | `synth_v2.py` (env `ST=BEHUMA`) | `results/synth_v2_BEHUMA.json` | Sec. 5.2 |
| Aircraft echoes along each ADS-B track, 1-s blocks, full day | `antmap_run.py` → `antmap_blocks.py` | `results/antmap/am_*.json.gz` | Supplementary S5 |
| Censored (Tobit) fit of the receiver gain maps | `antmap_tobit.py` | `results/antenna_map_tobit.json` | Supplementary Fig. S2 |
| Gain maps applied to the ISS geometry | `antmap_iss_pred.py` | `results/antmap/iss_pred.json` | Sec. 11.2 |
| Antenna-map and interference figures | `fig_antenna_map.py`, `fig_interference.py` | `fig_antenna_map.*`, `fig_interference.*` | Figs. 3, 4 |
| BlueBird search with line removal | `bb_clean.py` | `results/bb_clean.json` | Sec. 8 |
| BlueBird null distribution with line removal | `bb_null_clean.py` | `results/bb_null_clean.json` | Sec. 8 |
| Comparison of the ionospheric term with GNSS TEC (ROB maps) | `gnss_tec_cmp.py` | `results/gnss_tec_cmp.json` | Sec. 7 |
| GNSS TEC figure | `fig_gnss_tec.py` | `fig_gnss_tec.*` | Supplementary Fig. S6 |
| Independent (hold-out) false-alarm test with the frozen criterion | `null_holdout.py` | `results/null_holdout.json` | Sec. 5.2 |
| Sensitivity of the ISS statistic to the coherent block length | `block_length_test.py` | `results/block_length_test.json` | Sec. 3.2 |
| Cluster bootstrap (by aircraft track) of the antenna maps | `antmap_boot.py` | `results/antenna_map_bootstrap.json` | Supplementary S5 |
| Linearized (one-step Newton) cluster bootstrap with 1,000 replicates, sandwich cross-check | `antmap_boot1s.py` | `results/antenna_map_bootstrap_1s.json` | Supplementary S5 |
| Ionospheric model checks: WGS84 height, conditioning, injection with other profiles | `iono_checks.py` | `results/iono_checks.json` | Sec. 7, Supplementary S8 |
| Close passes (< 400 km) repeated with line removal, nulls; check against the ISS echo | `close_scan.py`, `close_check.py` | `results/close_scan.json`, `results/close_check.json` | Supplementary S9 |
| Spectrograms of the close Starlink candidates | `fig_starlink.py` | `fig_starlink.png` | Supplementary Fig. S5 |
| Envisat: window statistic and glint-matched statistic, with nulls | `envisat.py`, `glint.py` | `results/envisat.json`, `results/glint_envisat.json` | Supplementary S2–S3 |
| Glint-matched search on the close passes at Oudsbergen (DTC and older Starlink, CubeSats), with nulls | `glint_close.py` | `results/glint_close.json` | Supplementary S9 |
| Measured C/N0 and E/N0 of the echoes; z versus energy | `detectability.py`, `close_echo.py` | `results/detectability_cn0.json`, `results/close_echo.json` | Sec. 9.4 |
| Radar-equation model: effective gains, cross sections, sensitivity, smallest facet, Envisat limits | `detect_theory.py`, `fig_detectability.py` | `results/detect_theory.json`, `fig_detectability.pdf` | Sec. 9, Table 6, Fig. 8 |
| Blind catalog-wide search: candidate passes, fast statistic, sign-reversed nulls, association, detection probability | `blind_enumerate.py`, `fastsearch.py`, `blind_search.py`, `blind_analyze.py`; protocol in `PROTOCOL.md` | `results/blind_passes.json`, `results/blind_*.json`, `results/blind_summary.json` | Sec. 10 |
| **Phase 2 and feasibility (9 Nov 2026)** | | | |
| Phase 2 plan: predicted passes of the target classes, cross section needed, physical limit; file selection | `phase2_plan.py`, `phase2_select.py` (needs `data/tle/phase2/`) | `results/phase2_plan_passes.json` (462 MB, not in the deposit: rerun `phase2_plan.py`), `results/phase2_candidate_files.json`, `results/phase2_rank.json`; the frozen selection read by the search is `data/planning/phase2_files.json` (the 788 retrieved recordings) and `data/planning/phase2_byfile.json` (passes searched in each) | Sec. 10.1 |
| Phase 2 search (frozen statistic and sign-reversed null for every selected pass) | `phase2_search.py <shard> <nshards>` (env `P2DIR`, `OUTTAG`) | `results/p2_cloud*.json`, `results/p2_dev*.json`, logs in `results/p2logs/` | Sec. 10 |
| Phase 2 analysis: thresholds from the Phase 2 nulls, association, detection probability; sensitivity with the Phase 1 thresholds frozen; out-of-sample and effective false-alarm rates, thresholds without null removal | `phase2_analyze.py` (WAV header fields in `results/phase2_meta.json`; env `P2GLOB` selects the shards of one run) and `far_checks.py` (also called by `blind_analyze.py`) | `results/phase2_summary.json` (keys `far_checks`, `no_null_removal`) | Sec. 10.2 |
| Phase 2 against the Phase 1 prediction; ISS by sigma_min | `phase2_compare.py` | `results/phase2_compare.json` | Sec. 10.2, Fig. 9 |
| Detection with the satellite in sunlight or in the Earth's shadow (exploratory, both phases) | `eclipse.py` | `results/eclipse.json` | Sec. 10.4 |
| Feasibility test of a target (library and command line) | `brams_feasibility.py` | — | Sec. 10.5 |
| Feasibility curve fitted on Phase 1, prediction of Phase 2 | `feasibility_validate.py`, `fig_feasibility.py` | `results/feasibility_curve.json`, `results/feasibility_validation.json`, `fig_feasibility_flow.pdf`, `fig_feasibility_valid.pdf` | Sec. 10.5, Supplementary S6–S7, Figs. S3–S4, Table S2 |
| Uncertainties by bootstrap over recordings; predictive interval and error decomposition of the feasibility test | `stats_cluster.py` | `results/cluster_stats.json` | Sec. 10, Table 7 |
| Shadow versus sunlight stratified by receiver, margin, elevation and hour; Mantel-Haenszel odds ratio (exploratory) | `shadow_strata.py` | `results/shadow_strata.json` | Sec. 10.4 |
| Blind-search figure (ISS versus sigma_min; shadow and sunlight) | `fig_blind.py` | `fig_blind.pdf` | Fig. 9 |
| Receiver timing offsets relative to Oudsbergen (version-5 receivers; Ottignies), shared by the timing scripts | `rx_offsets.py` | — | Sec. 6 |
| Graphical abstract | `fig_graphical_abstract.py` | `graphical_abstract.png` | — |
| **Further checks** | | | |
| Doppler-equivalent offset versus geometric along-track difference (OEM), receiver offsets (version-5 receivers, Ottignies), light time | `iss_checks2.py` | `results/iss_checks2.json` | Secs. 2.4, 6.1 |
| ISS timing fits with non-overlapping 2-s windows: four interleaved sets, formal errors, latency model per set | `iss_nonoverlap.py` | `results/iss_nonoverlap.json` | Sec. 6.1 |
| Phase 2 detection rates and pooled out-of-sample false-alarm rate with Maubeuge reported separately | `maubeuge_separate.py` | `results/maubeuge_separate.json` | Sec. 10.2 |
| Angle between the bistatic bisector and candidate reflector normals (nadir, Sun, Sun-tracking) for Starlink passes | `sun_geometry.py` | `results/sun_geometry.json` | Sec. 10.4 |
| Reflection versus unintended emission: statistic along the bistatic and the one-way Doppler | `uemr_test.py` (needs the Phase 2 recordings; 200 detections drawn from the recordings of 19-26 Oct listed in `data/raw_wav_phase2/SUBSET_CHECKS.txt`) | `results/uemr_test.json` | Sec. 10.3 |
| End-to-end injection into Phase 2 recordings through the frozen blind search, 96 recordings, clean tracks (no ISS module, no detected echo on the track) | `inject_test.py` (needs the Phase 2 recordings; files drawn from `data/raw_wav_phase2/SUBSET_CHECKS.txt`; `SUBSET=none` uses every file found; env `NFILE`) | `results/inject_test.json` | Sec. 10.3 |
| Search along the one-way Doppler (unintended emission) with its own nulls, a third of the Phase 2 recordings at hand | `uemr_search.py <shard> <nshards>` (frozen list `data/raw_wav_phase2/SUBSET_UEMR_SEARCH.txt`; env `OUTTAG`), `uemr_search_analyze.py` (env `UEMRGLOB`) | `results/uemr_search_*.json`, `results/uemr_search.json` | Sec. 10.3 |
| Alternative explanations of the sunlight effect: orbital altitude (crude, near-zenith and stratified), null exceedances by UTC window | `confound_checks.py` | `results/confound_checks.json` | Sec. 10.4 |
| ISS timing in Phase 2: residual Doppler, joint fit, time-base correction; comparison with the archived NASA OEM files (`data/oem_archive/`) and with later element sets | `phase2_iss_timing.py` (needs the Phase 2 recordings; all 45 recordings with an ISS detection), `phase2_iss_summary.py` (passes, latency, OEM comparison by age, pass order) | `results/phase2_iss_timing.json`, `results/phase2_iss_summary.json` | Sec. 6 |
| ISS timing of both campaigns against the OEM predictions with one rule (latest file created before the pass; age from its initial state and from its creation), per-pass table | `iss_oem_combined.py` | `results/iss_oem_combined.json` | Sec. 6.2, Table 5, Supplementary S10 |
| SGP4 element set fitted to three days of each OEM file: along-track residual, its spectrum, and comparison with delta t at the Phase 2 passes | `sgp4_oem_fit.py` (needs `data/tle/phase2/`) | `results/sgp4_oem_fit.json` | Sec. 6, Supplementary S10 |
| Cross-receiver coincidence of the Phase 2 detections and of the null exceedances | `coincidence.py` | `results/coincidence.json` | Sec. 10.3 |
| Other objects of the Phase 2 element sets (debris included) whose Doppler matches the echo track of each detection; test with a common tau at several receivers | `catalog_association.py` (needs `data/tle/phase2/`) | `results/catalog_association.json` | Sec. 10.3 |
| Sunlight effect of the V2-mini passes in one logistic model (receiver, margin, elevation, altitude, UTC window), errors clustered by satellite and by recording | `shadow_glm.py` | `results/shadow_glm.json` | Sec. 10.4 |
| Supplementary Table S2 (ISS passes against the OEM predictions) | `supp_table_s10.py` | `results/table_s10.tex` | Supplementary S10 |
| Every number quoted in the paper, from the result files; check of the manuscript against them | `paper_numbers.py`, `check_numbers.py` (needs the LaTeX sources of the manuscript, which are not part of this repository: place them in `paper/` or set `PAPERDIR`; without them the check is skipped) | `results/paper_numbers.json` | — |
| Comparison of two result folders; headline numbers | `compare_results.py`, `key_numbers.py` | — | `REPRODUCTION.md` |
| Development versions kept for reference, not used for the published numbers | `fastsearch_v1.py`, `fastsearch_v2.py`, `fast_validate.py` (`results/fast_validate.json`), `bb_track.py`, `fig_iss_detection.py` | — | — |

`iono_fit.py` and `air_overlay.py` are helper modules that other scripts load. `figstyle.py` holds the common figure style. `fig_graphical_abstract.py` draws `graphical_abstract.png`.

## Notes

- **Frame conversion.** SGP4 outputs positions in TEME. The rotation to the Earth-fixed frame must use the Greenwich mean sidereal time (`sgp4.propagation.gstime`), not the Earth rotation angle; the latter introduces an error of about 26 km and a spurious timing offset of several seconds.
- **Time base.** Each WAV file carries its own timing record (`BRA2` chunk): sample indices at each GPS pulse per second for the analog receiver of Ottignies, and time stamps of each 1008-sample block for the SDR receivers. The analysis uses the file start time in the metadata; `pps_audit.py` shows that this is exact to 0.1 ms at Ottignies but off by up to ±15 ms at the SDR receivers, whose start time inherits the jitter of the first block time stamp. The ISS results correct every SDR receiver for this (`iono_fit_v2.py`).
- **Configurations.** Two transmitter-position configurations are used, and each script selects its own. The searches and the products archived from the original runs (aircraft timing, ISS detection table without line removal, null distributions and hold-out nulls, synthetic echoes, block-length test, antenna maps, BlueBird, Envisat and close-pass searches, Phase 1 and Phase 2 blind searches and their selection, analysis and metadata, one-way emission search) use the transmitter position written in each file, as in the original runs; these scripts set `BRAMS_TX_HEADER=1` in their first line. The revised timing and echo-energy analyses (`iss_table_v2.py`, `iono_track.py`, `oem_cmp.py`, `iono_fit_v2.py`, `iono_sens_v2.py`, `iono_checks.py`, `iss_checks2.py`, `iss_nonoverlap.py`, `phase2_iss_timing.py`, `detectability.py`, `close_echo.py`, `inject_test.py` and the figure scripts) use one position for every receiver, 50.09756 N, 4.58852 E, 236 m (`bwav.py`). The headers differ from it by about 4 m horizontally and 37 m in altitude at Oudsbergen and by about 0.28 km at Ottignies. Scripts that load `blind_analyze.py` or `phase2_analyze.py` inherit the historical configuration. The selection of Phase 2 (`phase2_plan.py`) used the planning station coordinates of `data/stations_geom.json` and the planning sensitivity model (`data/planning/planning_model.json`); Fig. 9a and the selection caps are on that scale. For receivers without their own calibration, the Phase 1 analysis (`blind_analyze.py`) assigned the Ottignies gain to Humain and Uccle, and the Phase 2 planning assigned the Oudsbergen gain to Humain and the Ottignies gain to Uccle; Fig. 9a preserves these archived conventions and is not a comparison with the recalibrated cross sections of Section 9.
- **Recordings used as nulls.** `close_scan.py`, `glint.py`, `glint_close.py`, `envisat.py` and `null_holdout.py` choose their recordings without a target among all files in the folder. They reproduce the published results only when `data/raw_wav/` and `data/raw_wav_envisat/` hold exactly the 293 and 40 recordings of `FILES.txt`; the BlueBird search and the development nulls of the ISS search read fixed lists (`bluebird_files.json`, `null_files.json`).
- **Language.** Comments in the scripts are partly in Portuguese.
- **Structure.** Several scripts load shared code with `exec(open(...).read())` instead of imports; this keeps each step self-contained but makes unit testing harder, and is a known limitation.

## How to cite

Scripts that read the raw files accept `WAVDIR=<folder>` to point to the WAV files. If you use this code, please cite the paper (reference to be added after publication) and the data sources listed in `DATA_LICENSES.md`, in particular the BRAMS dataset (doi:10.18758/x620ti37).
