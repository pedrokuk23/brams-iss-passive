# Detection protocol

Version 3, frozen on 9 Nov 2026 after the development run on 3-5 November 2026 (Phase 1). Items 2b, 2c, 4 and 4b were
revised during Phase 1, as recorded below; Phase 2 (19 October - 2 November 2026) applies this version unchanged,
including the thresholds of item 4 re-estimated with the same rule on the Phase 2 nulls.

One procedure for every object and every recording.

1. Candidate passes. Every (object, recording) pair in which the predicted bistatic Doppler stays in the usable band
   (30 Hz < |f_D| < 900 Hz) for at least 20 s of the 5-min file and the smallest geometric mean range sqrt(R1 R2)
   during those seconds is at most 600 km. Objects: the full catalog snapshot (data/tle/leo.txt) plus the ISS,
   the BlueBird satellites and Envisat with their own element sets (epoch closest to the recording).
2. Statistic. Line removal at the SDR receivers; dechirp along the element-set Doppler; 0.5-s coherent blocks; mean
   block power on the grid tau in [-6, 6] s (0.5 s), nu in [-12, 12] Hz (0.5 Hz); z = (P - median)/(1.4826 MAD) over
   the grid; maximum of z (code/fastsearch.py, validated against code/synth.py:search within 6% on known echoes).
2b. Burst blanking (revision of 9 Nov 2026, decided on the development data of Phase 1 after inspecting the largest
   values of the first complete run, which were short broadband bursts crossing the predicted track): for each timing
   offset, a 0.5-s block is excluded when the median power at residual frequencies of 13-40 Hz on either side of the
   search band exceeds twice the median of that quantity over the blocks of the search. Applied identically to the
   candidate and null searches. Known echoes (ISS complex, DTC Starlink, SPACEMOBILE-008) keep 90-100% of their z.
2c. Limiter (same revision): before averaging, the power of each block in each cell is divided by that block's noise
   level (median of the 13-40 Hz ring divided by ln 2) and capped at 20 (13 dB), so that a narrowband signal that the
   fast-sweeping prediction crosses within one block near closest approach cannot dominate the mean. On the development
   data this lowers such crossings and the largest null values from 25-50 to 3-6, while real echoes keep z of 14-70.
3. Detection. z >= z_th(receiver) with the peak inside the grid (|tau| < 6 s and |nu| < 12 Hz).
4. Threshold. For each receiver, z_th is the value exceeded by a fraction 1e-3 of null searches inside the grid.
   Null search: for every candidate pass, the same Doppler history with its sign reversed, in the same recording.
   No orbiting object produces a Doppler history that rises from negative to positive, while the noise, interference
   and meteor echoes are those of the real search.
   (Revision of 8 Nov 2026, made after a partial run and before any result was used: nulls built by applying a pass
   to another recording of the same receiver contained real echoes, because satellites of the same orbital plane
   follow nearly the same Doppler history a few minutes apart. That design is discarded.)
4b. Association. Candidates are taken in decreasing z. A candidate whose strongest 2-s window lies on the echo track
   of a stronger candidate in the same recording (raw frequency within 3 Hz of that track within 1.5 s) is the same
   echo: it is counted as an ambiguous association, not a detection; otherwise it is a new echo attributed to it. The same rule is applied to the null searches: a null
   peak that falls on an echo attributed to a real object (within 5 s and 5 Hz) is a crossing of that echo, which the
   association removes in the real search as well, and the false-alarm rate is the rate of the remaining null peaks.
   The thresholds and the association are iterated, starting from z = 7, until they no longer change.
4c. The ISS modules and the vehicles docked to it, which have their own element sets, form one object (ISS complex).
5. Reporting. Detection probability per object class and receiver with Clopper-Pearson 95% intervals; expected
   number of false alarms = (number of searches) x (false-alarm rate); detection probability versus the echo energy
   predicted by the radar equation for a unit cross section.
6. Phase 2 (validation campaign, 19 Oct - 2 Nov 2026, element sets from Space-Track gp_history with mean motion > 14.85
   rev/day) uses items 2-5 unchanged. Item 1 is replaced, before any Phase 2 recording is retrieved, by a selection on
   predicted geometry only (code/phase2_plan.py, code/phase2_select.py):
   - target classes: ISS complex, AST SpaceMobile, Starlink DTC and V2-mini, and Starlink V1.x as a control class;
   - a pass is searched when it is physically detectable, sigma_min <= sigma_max = 4 pi A^2 / lambda^2 with A the total
     flat area of its type (ISS: no limit; BlueBird Block 2: 223 m^2; BlueWalker 3: 64 m^2; Starlink V2-mini/DTC: 116 m^2;
     V1.x: 29 m^2), and its Phase 1 detection probability is at least about 1%: sigma_min <= 300 m^2 (V2-mini, V1.x),
     1000 m^2 (DTC), 3000 m^2 (AST), any (ISS complex);
   - recordings: every slot with an ISS-complex pass at the five receivers, and at Oudsbergen and Humain the 200 slots per
     receiver with the largest expected number of detections (sum of the Phase 1 detection probabilities of their passes).
   In Phase 1, passes outside the physical limit were detected at 0.14% (12 of 8511), the false-alarm rate.
6b. Correction of the class assignment of Phase 2 (9 Nov 2026, after the Phase 2 search and before any class result was
   reported): objects named SPACEMOBILE-SSRnn DEB (launch debris, 245 passes) had been assigned to the AST SpaceMobile
   class by name prefix and are reported as 'Other'; BlueBird Block 1 (SPACEMOBILE-001 to -005) carries a 64-m^2 array
   like BlueWalker 3, so its passes with sigma_min above the 64-m^2 limit (22 passes) are reported as beyond the physical
   limit. Thresholds, association and all other classes are unaffected (code/phase2_analyze.py).
7. Phase 2 sensitivity: the same association with the thresholds frozen at their Phase 1 values
   (results/phase2_summary.json, key frozen_phase1_thresholds).
8. Exploratory, not part of the frozen procedure: detection probability with the satellite in sunlight or in the Earth's
   shadow at closest approach (code/eclipse.py), noticed in Phase 2 and checked on Phase 1.

## Notes added after the freeze (10-11 Nov 2026; the procedure above is unchanged)

- The freeze is declared by the author. This file was edited after Phase 2 (items 6b, 7, 8 and these notes), and its
  modification time is not evidence of the date of the freeze.
- Item 6: the 200 slots per receiver at Oudsbergen and Humain include 21 that are also ISS-complex slots, so the paper
  counts 179 further slots per receiver (788 recordings in total).
- Item 6, last line: the Phase 1 figure of 12 of 8,511 passes beyond the physical limit was computed with the first
  version of the detectability model. With the revised model of the paper (energy of the glint lobe, n_b = 2T, the
  receiver's own threshold), passes with margin M < 1 were detected in 13 of 9,978 cases in Phase 1 and 14 of 3,524 in
  Phase 2 (`results/feasibility_validation.json`, key `impossible`). The selection of Phase 2 itself used the first
  version and is not changed.
- Checks of the false-alarm rule (out of sample and with explained nulls removed from the denominator) and of the
  thresholds without null removal are reported by `code/far_checks.py`; they do not change the procedure.
- Everything below was done after the Phase 2 results were examined (11 Nov 2026) and is
  reported in the paper as such: the grouping by sunlight in the feasibility test; the revised detectability model
  (energy of the glint lobe, n_b = 2T, the receiver's own threshold, logistic P_d(M) with a slope prior, gain
  corrections fitted on Phase 1); the out-of-sample false-alarm checks; the end-to-end injection; the emission tests
  (`uemr_test.py`, `uemr_search.py`); the stratifications by hour, elevation and orbital altitude
  (`shadow_strata.py`, `confound_checks.py`); the geometry test (`sun_geometry.py`); the ISS timing of Phase 2 (`phase2_iss_timing.py`). None of them changes the search,
  the thresholds or the association of items 1-6.
- Item 6b: "before any class result was reported" — the reclassifications were made after the search output existed.
