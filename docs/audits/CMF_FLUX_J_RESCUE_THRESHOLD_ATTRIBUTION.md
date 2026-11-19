# CMF_FLUX J_RESCUE threshold attribution (Phase 2.3)

Audit of the 75,738 NRESC events per observer render on
`m15_2p82d_SOB_l14_lte_adp2_tauto` (stock OBS path: DO_NONMONO_SOB=F,
DO_NONMONO_KNOP=F, USE_LAM_ES=T, FG_J_CMF_V12 + LAM_ES_TC,
J_RESCUE=T, RESC_JMAX=1.0E-3, DAMP_FJ=1.0).

Frozen baseline executable:
`exe_frozen/844b17b0fb430c6e_79b0651256a8b508/cmf_flux.exe`
(SHA256 prefix `79b0651256a8b508`).
Diagnostic-extended executable (Phase 2.3 JEV fields, JDMP dump):
`exe_frozen/844b17b0fb430c6e_d0125473c474290c/cmf_flux.exe`
(SHA256 prefix `d0125473c474290c`). Event data written by this exe is
diagnostic-only; OBSFLUX byte-identity vs the baseline freeze is
verified in sec. 5.

## 1. Executive summary

The 75,738 NRESC events per render are fully deterministic
(byte-identical OBSFLUX and identical event counts across
single-thread and 8-thread renders on the frozen exe) and are
driven by the INNER-boundary velocity inversion this model carries
— not thread noise, not sporadic instability.

Attribution vs the phase objective:

- **(B) systematic formal-solver behavior in compressed cells —
  the primary mechanism.** Every event sits at depth indices
  59-70, where the reconstructed SIGMA_TRUE is -5.2 to -155.6
  (median -31.0). V12 clamps this input to -0.999. J climbs
  smoothly inward from ~1e-9 to ~1e-3 across the whole compressed
  column; the offending values live at the inner boundary.
- **(C) possibly-legitimate boundary J — partially open.** At
  idx 70 the local source scale (ZETA=ETA/CHI ~ 3e-6) is ~300x
  below computed J. A scattering-dominated diffusion ceiling
  (J ~ S/(1-omega)) could justify ~1e-3 only if omega ~ 0.997;
  cannot be fully ruled out on these fields, but the tail events
  (up to 94x cap), the uniform spike across all 44,924 frequencies,
  and 36% F_DID_NOT_CONV incidence tilt toward artifact.
- **(A) cap miscalibration — secondary.** The KNOP-derived "5e-5
  peak" calibration came from a different solver/model and does not
  transfer; but the sweep shows the cap is doing real work, not
  clipping noise: at 2e-3 the event count only falls to 64k and the
  spectrum diverges 78% in spectral mass — the field inflates into
  whatever headroom it gets.

Recommendation: **retain RESC_JMAX=1.0e-3** (see sec. 11). Raising
it admits larger unverified fields into the warm-start chain and
materially reshapes the emergent spectrum without resolving the
events. The real fix is a formal solution valid at sigma_true
<-1 (the KNOP arbitrary-velocity path, Knop et al. 2009), not
cap recalibration.

## 2. Baseline reproducibility

Three single-thread A/A renders (job 13337357) and a ten-render
ensemble at the production thread count of 8 (job 13337358), all on
the frozen `79b0651256a8b508` exe with identical model inputs.

| render | threads | NRESC | OBSFLUX sha256 (prefix) | termination |
|--------|---------|-------|-------------------------|-------------|
| aa1 | 1 | 75738 | bd78a10852347cde | GLOBAL_CMF_GRID_UNDERCOVERAGE |
| aa2 | 1 | 75738 | bd78a10852347cde | GLOBAL_CMF_GRID_UNDERCOVERAGE |
| aa3 | 1 | 75738 | bd78a10852347cde | GLOBAL_CMF_GRID_UNDERCOVERAGE |
| e01-e10 | 8 | 75738 | bd78a10852347cde | GLOBAL_CMF_GRID_UNDERCOVERAGE |

Measured facts:

- Single-threaded and 8-thread renders produce byte-identical OBSFLUX
  and identical NRESC=75,738. The event count is fully deterministic;
  it is not thread-scheduling or reduction-order noise.
- Every render terminates in `OBS_FRAME_SUB_V9` with
  `GLOBAL_CMF_GRID_UNDERCOVERAGE` - the same expected production
  termination as the Phase 2.2 render. Not a completed model.
- NES_SAN=0, NFATAL=0, F_DID_NOT_CONV=0, FEDD-trigger=0,
  first-frequency events=0 (from JTOT/JEND, not throttled prints).
- Event freq indices span 1964-49064 with 44,924 unique frequencies
  affected - essentially the entire continuum grid, not a few bad
  frequencies.
- Depth histogram: idx 70 (ND, inner boundary) carries 66,035 events
  (87%); idx 59-69 carry the remaining 9,703.
- rjmax/RESC_JMAX: min 1.0000, med 1.17, p95 12.7, p99 51.6,
  max 93.8.
- Restore delta djmax = max|RJ_before - JNU_PREV|: med 1.98e-4,
  max 9.3e-2. djmax/rjmax: med 0.17.

Orientation correction vs the Phase 2.2 memo: depth index increases
INWARD (R(59)=1.601e14 cm, R(70)=1.455e14 cm in 1e10-cm units; V
rises 2632 -> 7074 km/s inward). The events sit at the INNER
boundary region of the non-monotonic velocity field (V increases
inward, dV/dR < 0 => true sigma < -1), not the outer edge.

## 3. Active-cap attribution

Production cap (obs/COMP_JCONT_V4.INC:740-748):

```
RESC_CAP   = MIN(RESC_JMAX, 1.0E4*MAXVAL(ABS(JNU_PREV)))
RESC_CAP_IT= MIN(1.0E10, MAX(1.0E4*MAXVAL(ABS(JNU_PREV)),1.0E-3))
```

Measured on the extended-record run (rmax_1000, 75,738 events):

- **ABSOLUTE limiter: 75,738 / 75,738 events.** `CAP_PREV` ranges
  9.64-10.0 (median ~10.0) — the absolute cap is always the smaller
  term. Raising RESC_JMAX genuinely moves the bound; it will not be
  silently overridden by the relative cap.
- `prevj` = MAXVAL(ABS(JNU_PREV)) sits at ~9.89e-4 — the previous
  field pins at ~0.99x the cap. Each restore writes a field that
  peeks just under the cap, the next frequency recomputes ~1e-3+,
  triggers again: a stable pin, not a runaway.
- The restore's replacement field itself carries a boundary spike
  (JP at idx 70 = 9.89e-4 in the JDMP profile) — the pin is
  structural, inherited from the previous frequency's clipped spike.

## 4. True-sigma reconstruction (measured)

`SIGT` per event = SIGMA_TRUE_CORE run on a diagnostic temp array:

- **Median SIGMA_TRUE at event depths: -30.97.** Range: -155.6 to
  -5.24. p95 = -5.37. **All 75,738 events sit at true sigma < -1** —
  the worst at -155.6. Not marginal: the events are deep inside the
  compression regime V12's -0.999 clamp conceals.
- The physical field (JDMP): outer ejecta ~49 km/s at R~1.25e15 cm
  rising inward to 7,074 km/s at R~1.45e14 cm — a velocity inversion
  (fast inner ejecta under a stalled outer shell, the interacting-SN
  profile this model carries). dV/dR is strongly negative there.
- SIG=-0.999 in every event = the clamped solver input, a constant.
- Offline replica `sigma_true_check.py` reproduces the routine to
  ~1e-13 on power-law profiles including both boundaries; the render
  calls the production routine itself, so orientation and stencil are
  definitionally identical.

## 4a. Where the J inflation lives (JDMP profile)

Full radial profile for 200 dumped events (sample freq 1965+):

| idx | R (1e10 cm) | V (km/s) | SIGMA_TRUE | RJ (pre-restore) | JNU_PREV |
|-----|-------------|----------|------------|------------------|----------|
| 1 | 125,312 | 49 | -1.04 | 5.06e-9 | 5.21e-9 |
| 5 | 116,082 | 49 | -1.06 | 5.92e-9 | 6.10e-9 |
| 66 | 14,556 | 6,943 | -33.5 | 6.99e-4 | 6.83e-4 |
| 68 | 14,548 | 7,058 | -31.3 | 9.64e-4 | 9.43e-4 |
| 70 | 14,547 | 7,074 | -31.0 | 1.01e-3 | 9.89e-4 |

J climbs smoothly inward — ~1e-9 at the outer edge, ~1e-7 by
idx ~15, ~1e-4 by the low 60s — reaching ~1e-3 exactly at the inner
boundary (idx 70). The whole compressed column inflates; only the
inner boundary cells cross the cap. ZETA (= ETA/CHI, the local
source function, OPACITIES_V4.INC:171) at those cells is ~3e-6 —
computed J sits ~300x above local source scale. A scattering-
dominated diffusion ceiling J ~ S/(1-omega) could reach ~1e-3 only
at omega ~ 0.997; the tail events (up to 94x cap) exceed any such
ceiling. FEDD at the boundary converges to ~0.335-0.95 (finite,
physical-looking) — the Eddington factor gives no warning.

- F_DID_NOT_CONV=T on 27,433 events (36%) — on a third of affected
  frequencies the f iteration did not even converge, yet still
  produced finite over-cap J. The iteration's stopping criterion is
  blind to this failure class.

## 5. Diagnostic-exe integrity checks

- New-exe + diag-OFF (`off_newexe`): OBSFLUX sha256
  `bd78a10852347cde` — **byte-identical** to every baseline render
  (both exes). The diagnostic extension and rebuild are numerically
  inert end-to-end.
- New-exe + diag-ON, RESC_JMAX=1e-3 (`rmax_1000`): **NRESC=75,738 —
  identical to baseline.** Same freq range, same depth histogram.
  The extension does not perturb the rescue decision or count.

## 6. Dimensional and physical audit of the J ceiling

RJ is the comoving-frame mean intensity J_nu from the formal solve
(`RJ(1:ND)=TC(1:ND)`, obs/cmf_flux_sub_v5.f:1032), in the code's
intensity units (same units as ETA/CHI = the source function).
The field is per-frequency; there is no grid or model normalization
in the stored values.

The source comment (COMP_JCONT_V4.INC:720-738) calibrates
RESC_JMAX~1e-3 from "legitimate J peaks near 5e-5 ... measured on
the converged KNOP field". That calibration came from the KNOP
interaction-model run, a different formal solver on a different
model than this stock V12 render. The measured event RJ (median
1.13e-3 at idx 70) sits ~20x above the quoted peak - but the
adjacent-frequency field is smooth: djmax/rjmax med 0.17, i.e. the
J field varies only ~17% between neighbouring frequencies at the
same depth. A runaway would show J >> JNU_PREV; instead JNU_PREV
itself sits at the same ~1e-3 scale (CAP_PREV ~ 1e4*1e-3 = 10).

Dimensional verdict: RJ is a stored mean intensity in the code's
intensity units (the same units as ETA/CHI) with no normalization —
so a single absolute ceiling IS dimensionally consistent across the
grid. What a flat number cannot express is per-cell physics: the
meaningful bound at a boundary cell is set by that cell's source
scale and scattering fraction, not a global constant. The 1e-3
value is dimensionally valid but physically uncalibrated for this
model family — and the sweep data shows the artifact makes any
single number arbitrary in a different way: whatever ceiling is
chosen, the field inflates toward it.

## 7. Threshold sensitivity (RESC_JMAX sweep)

Same model, same diag-exe, 8 threads, J_RESCUE=T, DAMP_FJ=1.0,
only RESC_JMAX varied (job 13337580):

| case | RESC_JMAX | NRESC | uniq freqs | OBSFLUX | spec-mass diff |
|------|-----------|-------|-----------|---------|----------------|
| baseline | 1.0e-3 | 75,738 | 44,924 | bd78a108 | - |
| rmax_1100 | 1.1e-3 | 71,752 | ~44k | d81a2d4f | 10.4% |
| rmax_1250 | 1.25e-3 | 73,821 | 43,742 | a1bfe073 | 16.9% |
| rmax_2000 | 2.0e-3 | 64,032 | 37,808 | 34a7ed5b | 78.3% |

CAP_PREV never became the limiter (ABSOLUTE on 100% of events in
every case — cap_prev ~10 >> any swept value).

Readings:

- Event counts fall modestly and NON-monotonically (1.25e-3 >
  1.1e-3): each cap changes JNU_PREV's pin, so the field plays out
  differently. The baseline over-cap histogram was a floor on the
  artifact's amplitude — the restore was capping the MEASURED
  rjmax, so the true field values above the cap were never
  recorded at full amplitude.
- At 2e-3 still 64,032 interventions — the artifact does not
  saturate at low thresholds; the boundary field keeps inflating
  relative to any plausible cap in this range.
- All cases terminate at GLOBAL_CMF_GRID_UNDERCOVERAGE (expected
  production termination), no FATAL, no non-finite events.

## 8. Spectral effect of restoration

Baseline OBSFLUX vs each sweep case (aligned-array diff; OBSFLUX
layout = freq block + flux block + Luminosity cumulative tail):

| case | med rel | p95 rel | p99 rel | spec-mass diff | luminosity med | final lum |
|------|---------|---------|---------|----------------|----------------|-----------|
| 1.1e-3 | 0.0 | 0.22 | 18.8 | 10.4% | 14% | identical |
| 1.25e-3 | 0.0 | 0.44 | 16.8 | 16.9% | 30% | identical |
| 2.0e-3 | 0.0 | 4.29 | 41.1 | 78.3% | 141% | identical |

- Half the spectrum is untouched (med_rel=0) — differences
  concentrate in the affected frequency region, reaching 40x
  relative change at p99 for the widest cap.
- Total cumulative luminosity is conserved (final values
  identical across cases) — the field redistribution changes
  where the luminosity accrues, not the total.
- Restore mechanics (JDMP): the whole RJ vector is replaced by
  JNU_PREV, but since JP ~= RJ outside the spike cells the net
  effect is confined to the inner-boundary cells — the restore
  swaps this frequency's clipped spike for the previous
  frequency's clipped spike, then propagates via the JNU_PREV
  warm-start chain into the next frequency's solve.

So J_RESCUE is load-bearing for the emergent spectrum: the
interventions materially change what comes out (10-78%
spectral-mass divergence between cap settings).

## 9. J_RESCUE boundary-coupling audit

`J_RESCUE` currently gates two additional inner-boundary corrections
independent of the rescue/restore logic:

1. Negative-DBB floor, obs/COMP_JCONT_V4.INC:117-121:
   `(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP .OR. J_RESCUE) .AND. DBB<0`
   sets DBB=1e-99, DDBBDT=0. Physical limit of an inwards-decreasing
   Planck gradient is zero diffusion flux; the floor keeps RJ(ND)>0
   for JFEAU_IBC_V2.
2. Negative-DTDR replacement, obs/cmf_flux_sub_v5.f:1546-1554:
   same gate; DTDR<0 with valid ROSS_MEAN(ND) and STARS_LUM is
   replaced by the luminosity-derived DTDR. Without valid Rosseland
   data the case falls through to the DBB guard.

Both corrections act at depth ND = 70 - the same index where 87% of
the NRESC events fire. J_RESCUE=F would therefore silently change
the inner boundary condition in addition to disabling the restore;
it was correctly excluded from the sweep protocol.

Measured activation in this model: the DTDR replacement prints to
LUER when it fires; no such line appears in any render log, and the
two printed DBB values (4.6e-5, 6.2e-5) are positive. The DTDR
correction is armed but inactive. The DBB floor has no message of
its own; a silent firing at frequencies where DBB is not printed
cannot be excluded from the log alone, but no evidence requires it.

A separate control (`J_RESCUE_BC_FIX`) would decouple the boundary
corrections from the restore gate; not implemented in this phase.

## 10. Scientific interpretation

The LAM_ES_TC update is a variable-Eddington-factor / Lambda-type
iteration (Hubeny & Mihalas 2014, ch. 13); it is not ALI - no
Rybicki-Hummer approximate Lambda operator is used. FEDD
convergence certifies the f iteration stopped; it does not certify
the formal solution in the compressed regime. 36% of events carried
F_DID_NOT_CONV=T — the iteration literally did not converge on
a third of the offending frequencies, yet still emitted
finite-but-inflated J.

Knop, Hauschildt & Baron (2009, A&A 496, 295, DOI
10.1051/0004-6361:200810290) built their formal solver because
monotonic-flow solvers fail exactly here: velocity inversions,
stagnation, compression. The measured sigma_true range (-5 to -156,
median -31) is the regime their method exists for; V12's -0.999
clamp makes the same input look like mild expansion.

What the measurements establish:

- The events are systematic (44,924/49k frequencies), deterministic
  (13 byte-identical renders), localized (idx 59-70, the inversion
  cells), and marginal in the median (1.17x cap) with a real tail
  (94x).
- The J field at those cells is inflated vs local source scale by
  ~300-1000x; whether that is an artifact or a legitimate
  diffusion-limit asymptote (J ~ S/(1-omega)) is not fully
  decidable from J/ZETA/ETA alone — the discriminator is an
  independent formal solver on the same grid, i.e. the
  DO_NONMONO_KNOP path.
- The restore is a correct stopgap: it substitutes the previous
  frequency's clipped boundary field — same spike shape, pinned
  just under cap — so JNU_PREV remains smooth in frequency.

Unresolved physics questions:

1. Is the ~1e-3 boundary J a correct photon-trapping asymptote or
   an inflated artifact? Decide by running the same model through
   the KNOP/FSM formal path and comparing J at idx 59-70.
2. Does the restored boundary J bias the emergent spectrum
   systematically (it feeds the Eddington factor -> flux)? The
   10-78% spectral divergence shows the field is load-bearing;
   direction of bias is not established.
3. FEDD converges to 0.33-0.95 at the boundary cells — a
   plausible-looking f on an implausible J. The Eddington
   iteration is blind to this class.

## 11. Recommendation

**Retain RESC_JMAX=1.0e-3.** The measured evidence does not
support recalibrating upward:

- Raising the cap to 2e-3 removed only 15% of events while letting
  the boundary field grow into the headroom — the artifact's
  amplitude is cap-relative, not absolute.
- Higher caps admit larger unverified fields into the JNU_PREV
  warm-start chain and diverge the spectrum up to 78% in spectral
  mass — a material physics change with no evidence it is toward
  truth.
- The threshold's KNOP-derived calibration (5e-5 peaks) does not
  transfer to this model's inversion physics, but the cap's actual
  role here is different: it is an early-warning bound on a region
  where the formal solver is out-of-domain. It is doing that job.

Better long-term criterion (future phase, needs review): bound J
against the local source scale (e.g. flag when RJ exceeds
kappa*ZETA or the diffusion-limit J=S/(1-omega) at that cell) —
a physically motivated per-cell test rather than a flat number —
AND/OR route the compressed cells to the arbitrary-velocity formal
solve. Do not change the production default under this audit;
the sweep data, not opinion, should drive that review.

## Appendix: run inventory

Run root:
`/viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/m15_2p82d_SOB_l14_lte_adp2_tauto/cmfgen/obs_ensemble/`

- Baseline: `aa{1,2,3}` (1 thread, job 13337357), `e01..e10`
  (8 threads, job 13337358); exe 79b0651256a8b508; PINIT
  `CMF_FLUX_PARAM_INIT.diag` sha bf8dba2d…
- Sweep (job 13337580, diag-exe d0125473c474290c):
  `off_newexe` (diag off, identity check), `rmax_1000` (1.0e-3),
  `rmax_1100` (1.1e-3), `rmax_1250` (1.25e-3), `rmax_2000` (2.0e-3).
- `baseline_summary.json` — machine-readable 13-render summary
  (hashes, event stats, terminations).
- `<run>/manifest.json` — per-render exe sha, input shas, env,
  threads, exit code, termination, output shas, event stats.
- `<run>/J_RESCUE_EVENTS` — JEV + JTOT (JEND writes only on clean
  exit; this model always aborts at UNDERCOVERAGE, so JTOT is the
  durable record).
- `rmax_*/J_RESCUE_DUMP` — 200 bounded JDMP radial profiles per run.
- `figs/jres_phase23_overview.{png,pdf}`,
  `figs/jres_phase23_luminosity.png`.
- Analysis: `analyze_jres.py`, `sigma_true_check.py`,
  `make_audit_figs.py`.

Known harness caveat (fixed): `run_one.sh` originally `cp -f`'d the
per-run PINIT through a hardlinked inode, so early manifests
recorded the last-written init instead of the per-run one; fixed
with `--remove-destination` + hashing the source PINIT, and all
existing manifests were repaired to the true values. Render data
was unaffected — OBSFLUX hashes prove all 13 baseline renders ran
the same effective config.
