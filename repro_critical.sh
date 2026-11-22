#!/bin/bash
# Recomputes from the raw WAV files the three critical results: ISS detections with line removal (Table 3), along-track offsets (Table 4) and hold-out false alarms.
# Run from the repository root with the 293 WAV files in data/raw_wav/; compare results/ with a saved copy afterwards.
set -e
cd code
export WAVDIR=../data/raw_wav/ DAYDIR=../data/raw_wav/ BBDIR=../data/raw_wav/
python3 pps_audit.py
python3 iss_table_v2.py
ISS="20261103_0040 20261103_0220 20261103_2355 20261104_0130 20261104_2305"
CLEAN=1 OUT=../results/iono_tracks_clean_sdr.json python3 iono_track.py $(for p in $ISS; do ls ../data/raw_wav/RAD_BEDOUR_${p}_*_SYS001.wav; done | grep -v BEOTTI)
OUT=../results/iono_tracks_otti.json python3 iono_track.py $(for p in $ISS; do ls ../data/raw_wav/RAD_BEDOUR_${p}_BEOTTI_SYS001.wav; done)   # separate file: results/iono_tracks.json also holds the Oudsbergen tracks
python3 -c "import json;a=json.load(open('../results/iono_tracks_otti.json'));b=json.load(open('../results/iono_tracks_clean_sdr.json'));m={k:v for k,v in a.items() if k.endswith('BEOTTI')};m.update(b);json.dump(m,open('../results/iono_tracks_v2.json','w'))"
TRACKS=../results/iono_tracks_v2.json OUT=../results/oem_cmp_v2.json python3 oem_cmp.py
python3 iono_fit_v2.py && python3 iss_checks2.py
python3 null_holdout.py
echo DONE
