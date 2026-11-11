#!/bin/bash
# Reproduces the results of the ASR manuscript in order. Run from the repository root after placing the
# BRAMS WAV files listed in data/raw_wav/FILES.txt in data/raw_wav/ (check them with: cd data/raw_wav && sha256sum -c SHA256SUMS).
# The steps run in dependency order, so the script works from an empty results/ folder; a step can also be run alone
# when its inputs are already in results/.
set -euo pipefail
mkdir -p results/figures results/antmap
cd code
export WAVDIR=../data/raw_wav/ DAYDIR=../data/raw_wav/ BBDIR=../data/raw_wav/ P2DIR=../data/raw_wav_phase2/
# run N shards of a script in parallel and stop if any of them fails
shards() { local pids=() s; for s in 0 1 2 3; do "$@" $s 4 & pids+=($!); done; for s in "${pids[@]}"; do wait "$s"; done; }
python3 pps_audit.py                                   # time base of every file
sh run_uniform.sh && python3 aircraft_stats_uniform.py && python3 fig_tau_hist.py      # aircraft timing (Table 2, Supplementary Fig. S1)
ISS="20261103_0040 20261103_0220 20261103_2355 20261104_0130 20261104_2305"   # the five ISS passes
# Doppler-residual tracks without line removal (Ottignies, five passes; Oudsbergen, four passes), used by synth.py
python3 iono_track.py $(for p in $ISS; do ls ../data/raw_wav/RAD_BEDOUR_${p}_BEOTTI_SYS001.wav; done) $(for p in 20261103_0220 20261103_2355 20261104_0130 20261104_2305; do ls ../data/raw_wav/RAD_BEDOUR_${p}_BEOUDS_SYS001.wav; done)
python3 iss_table.py && python3 iss_table_v2.py       # ISS detection without/with line removal (Table 3)
CLEAN=1 OUT=../results/iono_tracks_clean_sdr.json python3 iono_track.py $(for p in $ISS; do ls ../data/raw_wav/RAD_BEDOUR_${p}_*_SYS001.wav; done | grep -v BEOTTI)
python3 -c "import json;a=json.load(open('../results/iono_tracks.json'));b=json.load(open('../results/iono_tracks_clean_sdr.json'));m={k:v for k,v in a.items() if k.endswith('BEOTTI')};m.update(b);json.dump(m,open('../results/iono_tracks_v2.json','w'))"
TRACKS=../results/iono_tracks_v2.json OUT=../results/oem_cmp_v2.json python3 oem_cmp.py
python3 iono_fit_v2.py && python3 receiver_delay.py && python3 iono_sens_v2.py && python3 gnss_tec_cmp.py   # Table 4, Secs. 6-7
python3 iono_checks.py                                 # WGS84 height, identifiability and injection tests (Sec. 7)
python3 iss_checks2.py                                 # Doppler-equivalent offset vs geometric along-track, latency model, light time (Sec. 2.5, 6)
sh run_cat3.sh && python3 cat3_eval.py                 # first catalog scan (needed by the close-pass search)
python3 close_scan.py && python3 close_check.py        # close passes (Starlink) with line removal; ISS ambiguity check (Sec. 8.1)
python3 fig_starlink.py                                # Supplementary Fig. S3
EVDIR=../data/raw_wav_envisat/ python3 envisat.py && EVDIR=../data/raw_wav_envisat/ python3 glint.py   # Envisat (Supplementary S2-S3)
python3 glint_close.py                                 # glint search on the close passes (supplementary)
python3 detectability.py && python3 detect_theory.py && python3 close_echo.py && python3 fig_detectability.py   # Sec. 9, Table 5, Fig. 9
python3 null_run.py && STS='BEHUMA BEUCCE BEOUDS FRMAUB' OUT=../results/null_v2.json python3 null_v2.py && python3 null_holdout.py   # false alarms
ST=BEOTTI python3 synth_run.py && ST=BEOUDS python3 synth_run.py && python3 synth_eval.py && ST=BEHUMA python3 synth_v2.py         # synthetic echoes
python3 block_length_test.py
# antenna maps (full day of 4 Nov; ADS-B extract in data/adsb/adsb_be_nov04_day.jsonl.gz)
gunzip -kf ../data/adsb/adsb_be_nov04_day.jsonl.gz
for s in BEOTTI BEOUDS FRMAUB BEHUMA BEUCCE; do python3 antmap_run.py $s ../data/adsb/adsb_be_nov04_day.jsonl ../data/raw_wav ../results/antmap/am_$s.json; done
python3 antmap_blocks.py ../results/antmap && python3 antmap_tobit.py ../results/antmap 6 && NBOOT=10 python3 antmap_boot.py ../results/antmap 6 && NBOOT=1000 python3 antmap_boot1s.py ../results/antmap 6 && python3 antmap_iss_pred.py
# BlueBird and catalog scan
python3 bb2_run.py && python3 bb2_null.py && python3 bb_clean.py && python3 bb_null_clean.py && python3 bb_tle_update.py
# figures
python3 air_mf.py && python3 paper_figs2.py && python3 paper_figs3.py && python3 paper_figs4_v2.py && python3 fig_antenna_map.py && python3 fig_interference.py && python3 fig_gnss_tec.py && python3 bb_fig.py
cd .. && python3 -m pytest -q tests

# Phase 1 blind search (PROTOCOL.md items 1-5); about 4 CPU-hours
cd code
python3 blind_enumerate.py
export OUTTAG=blind_rerun; shards python3 blind_search.py; unset OUTTAG   # new shards; the published ones are not read
# BLINDGLOB stays exported: eclipse.py and feasibility_validate.py below re-run blind_analyze.py and must read the same shards
export BLINDGLOB='../results/blind_rerun*.json'
python3 blind_analyze.py
cd ..

# Phase 2 (validation campaign) and feasibility test; the search takes about 5 CPU-hours (shards can run in parallel).
# The shards of this run are written as results/p2_rerun*.json; phase2_search.py skips only recordings already in
# files of this tag, and phase2_analyze.py reads only them (P2GLOB), so the published shards can stay in place.
cd code
python3 phase2_plan.py && python3 phase2_select.py
export OUTTAG=p2_rerun; shards python3 phase2_search.py; unset OUTTAG
export P2GLOB='../results/p2_rerun*.json'           # kept for every later script that reads the Phase 2 shards
python3 phase2_analyze.py && python3 phase2_compare.py && python3 eclipse.py
python3 feasibility_validate.py && python3 fig_feasibility.py && python3 stats_cluster.py && python3 shadow_strata.py && python3 fig_blind.py
# further checks: attitude geometry of the shadow effect, reflection vs unintended emission, end-to-end injection
python3 sun_geometry.py && python3 uemr_test.py && python3 inject_test.py && python3 inject_clean.py
python3 confound_checks.py
python3 phase2_iss_timing.py && python3 phase2_iss_summary.py   # ISS timing in Phase 2 against the archived NASA OEM files (Sec. 6)
python3 iss_oem_combined.py && python3 sgp4_oem_fit.py         # both campaigns with one rule; SGP4 fitted to the OEM (Sec. 6, Supplementary S10)
python3 coincidence.py                                 # cross-receiver coincidence of the Phase 2 detections (Sec. 10.2)
python3 catalog_association.py                         # other catalog objects compatible with each detection (Sec. 10.2)
python3 shadow_glm.py                                  # sunlight effect in one logistic model, clustered errors (Sec. 10.3)
python3 supp_table_s10.py                              # Supplementary Table S2
OUTTAG=_rerun shards python3 uemr_search.py; UEMRGLOB='../results/uemr_search_rerun_*.json' python3 uemr_search_analyze.py   # one-way search, frozen list of recordings
python3 fig_graphical_abstract.py                      # graphical abstract (needs the Phase 2 results above)
python3 paper_numbers.py && python3 check_numbers.py   # every number quoted in the paper, checked against the manuscript
cd ..
