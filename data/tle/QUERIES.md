# Element sets: origin and how to retrieve them

All element sets are basic space situational awareness data from USSPACECOM, obtained from Space-Track.org
with the author's account, except `leo.txt`, taken from CelesTrak (see the table). They are redistributed here under the Space-Track "blanket approval for
transfer/redistribution of basic SSA data" (TLEs and OMMs), which requires citing Space-Track.org.
Anyone with a free Space-Track account can retrieve them again with the queries below (log in first,
then open each URL in the same browser session).

| File | Content | Query |
|---|---|---|
| `iss.txt` | ISS (NORAD 25544), element set of epoch 2026-11-03 12:44 UTC used for all ISS results | `https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/25544/EPOCH/2026-11-01--2026-11-06/orderby/EPOCH%20asc/format/3le` (the set used is the one of epoch 26307.53092409) |
| `leo.txt` | 15,097 objects in low Earth orbit, latest element set of each as of 2026-11-03 22:21 UTC (98% with epochs on 3 Nov), used for the catalog scan and Phase 1 of the blind search. Downloaded from CelesTrak (active group, https://celestrak.org/NORAD/elements/gp.php?GROUP=active, filtered to low Earth orbit) on 4 Nov 2026, 01:19 UTC; CelesTrak redistributes the USSPACECOM element sets and its names tag the Starlink direct-to-cell satellites with `[DTC]` (the CelesTrak format has no `0 ` before the name, unlike the Space-Track 3LE files) | to reproduce it later from Space-Track, query `gp_history` for each NORAD number with `EPOCH/<2026-11-03T22:21:00` and keep the latest set; the direct-to-cell satellites are listed in `dtc_ids.txt` |
| `dtc_ids.txt` | NORAD numbers and names of the 637 Starlink direct-to-cell satellites tagged `[DTC]` in `leo.txt`, so that the class can be assigned to element sets that lack the tag (Phase 2) | derived from `leo.txt` |
| `bluebird_gp_history.txt` | All element sets of SPACEMOBILE-008, -009 and -010 (NORAD 69589, 69590, 69591) with epochs from 3 to 8 Nov 2026 | `https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/69589,69590,69591/EPOCH/2026-11-03--2026-11-08/orderby/EPOCH%20asc/format/3le` |
| `envisat_gp_history.txt` | All element sets of Envisat (NORAD 27386) with epochs from 2 to 7 Nov 2026 | `https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/27386/EPOCH/2026-11-02--2026-11-07/orderby/EPOCH%20asc/format/3le` |
| `phase2/gp_2026-10-19_24.txt`, `phase2/gp_2026-10-24_29.txt`, `phase2/gp_2026-10-29_11-03.txt` | Phase 2 (19 Oct - 2 Nov 2026): all element sets of objects with mean motion above 14.85 rev/day (low Earth orbit below about 650 km) with epochs in each window, retrieved by the author from the `gp_history` class; the middle file holds about half of the sets per object per day, which still covers every day of its window. Each recording is processed with the set of epoch closest to it | `https://www.space-track.org/basicspacedata/query/class/gp_history/EPOCH/<start>--<end>/MEAN_MOTION/%3E14.85/orderby/NORAD_CAT_ID,EPOCH/format/3le` with the windows of the file names |

Cite as: U.S. Space Command, Space-Track.org, general perturbations element sets (gp and gp_history classes), accessed November 2026.

SHA-256 checksums of `leo.txt` and of the three Phase 2 history files (`phase2/`) are in `SHA256SUMS_catalog`.
