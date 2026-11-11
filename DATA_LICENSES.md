# Data sources and licenses

The code is MIT-licensed (see `LICENSE`). The data in `data/` and the derived products in `results/` come from third parties and keep their own terms:

| Folder | Content | Source | License / terms |
|---|---|---|---|
| `data/raw_wav/` (not included) | BRAMS raw WAV files, 5-min, 3–5 Nov 2026 (list in `FILES.txt`) | H. Lamy, *Raw WAV data from the BRAMS network*, BIRA-IASB, doi:10.18758/x620ti37 | CC BY 4.0 — cite the dataset |
| `data/stations/coords.json` | Receiver and transmitter coordinates, read from the BRAMS WAV metadata | as above | CC BY 4.0 |
| `data/adsb/` | ADS-B tracks extracted for Belgium, 09:55–10:35 UTC on 3 and 4 Nov 2026, and 04:55–21:35 UTC on 4 Nov 2026 (`adsb_be_nov04_day.jsonl.gz`) | adsb.lol, `globe_history_2026` (https://github.com/adsblol/globe_history_2026) | Open Database License (ODbL) 1.0 — attribution and share-alike apply to these extracts |
| `data/oem/` | ISS state vectors (CCSDS OEM, EME2000), excerpts of the files of 2 and 4 Nov 2026 | NASA JSC Trajectory Operations and Planning (https://spotthestation.nasa.gov/trajectory_data.cfm) | U.S. Government work |
| `data/giro/` | Autoscaled foF2 and hmF2, Dourbes DB049 (DPS-4D) | GIRO DIDBase (https://giro.uml.edu), data provider: Royal Meteorological Institute of Belgium | CC BY-NC-SA 4.0 — acknowledge the DB049 provider; follow the GIRO Rules of the Road |
| `data/tle/iss.txt` | ISS element set used in the paper (epoch 2026-11-03 12:44 UTC) | USSPACECOM via Space-Track.org | Basic SSA data, used and redistributed under the Space-Track blanket approval, which requires citing Space-Track.org. The exact queries are in `data/tle/QUERIES.md` |
| `data/tle/leo.txt`, `data/tle/dtc_ids.txt` | 15,097 LEO element sets for the catalog scan and Phase 1 of the blind search (latest per object as of 3 Nov 2026, 22:21 UTC); list of the Starlink direct-to-cell satellites | USSPACECOM element sets through CelesTrak (downloaded 4 Nov 2026; see `data/tle/QUERIES.md`) | Basic SSA data; cite CelesTrak (T. S. Kelso, https://celestrak.org) and Space-Track.org |
| `data/tle/envisat_gp_history.txt`, `data/tle/phase2/` | Element sets of Envisat (2-7 Nov 2026) and of the Phase 2 period | USSPACECOM via Space-Track.org | as above |
| `data/raw_wav_envisat/`, `data/raw_wav_phase2/` (not included) | BRAMS raw WAV files of the Envisat set and of Phase 2 (lists and checksums included) | BRAMS dataset, doi:10.18758/x620ti37 | CC BY 4.0 — cite the dataset |
| `data/tle/bluebird_gp_history.txt` | Element sets of SPACEMOBILE-008, -009, -010, epochs 3–8 Nov 2026 | USSPACECOM via Space-Track.org | as above |
| `data/gnss_tec/rob_vtec_subset.txt` | Vertical TEC over Belgium (48.5–52.5°N, 2.5–7°E, 0.5° grid) at the times of the five ISS passes, extracted from the near-real-time IONEX maps of the Royal Observatory of Belgium | ROB/STCE, https://gnss.be/SpaceWeather/Products/IONEX/ (Bergeot et al., 2014, doi:10.1051/swsc/2014028) | Publicly available; cite Bergeot et al. (2014) and the ROB GNSS service |
| `results/` | Products derived from the above | this work | Derived from the sources above; the ADS-B-derived files remain under ODbL and the ionosonde-derived values under CC BY-NC-SA |
| `data/oem_archive/ISS_OEM_*.txt` | Nine complete ISS Orbit Ephemeris Message files issued 19 Oct - 6 Nov 2026 (archive of the NASA ISS trajectory data, `https://nasa-public-data.s3.amazonaws.com/iss-coords/<date>/ISS_OEM/ISS.OEM_J2K_EPH.txt`) | NASA JSC, Trajectory Operations and Planning | U.S. Government work, not subject to copyright in the United States |
