# EXPONX model diagnostics — m15_2p82d_FRESH_sob_dieoff

Model run dir: /ptmp/sbrennan/SN2024pba/cmfgen_validation/.../cmfgen/2.82days
Build commit 9c5d5a25; DO_NONMONO_SOB=T, DO_NONMONO_KNOP=F, DO_NONLOC_SOB=F,
NEG_OPAC_OPT=ESEC_CHK, CHECK_LINE_OPAC default TRUE.

## NEG_OPAC file (production's own log, LST_ITERATION writes)

- 159 flagged (group,depth) cells logged; ~85 unique transitions across
  species: HI n=9-18 (27), He2 (48), C IV (3), O IV 2s2_8p-8s (1 group),
  O V 6z-6g (2), O VI (2), S IV, Sk IV.
- Depths: d15-30 (Ne 1e8-1e9) and d63-70 (Ne 1e10-1e14).
- Most negative logged TAU(Sob) = -0.984 (angle-free, group-summed).
- Every logged cell was crushed: CHIL->1.0, CHIL_MAT->1/NUM_SIM_LINES
  before SOBJBAR_SIM — GAM>0 at every flagged cell.

## Independently counted inversions (audit 5/6 inventory)

- 1027 (line,depth) cells with per-line CHIL_MAT<0 = 148 unique
  transitions x depths (mostly O IV Rydberg/intercombination).
- Per-line negatives do NOT reach EXPONX (GAM built from group sum).
  They enter the apportionment ST_ON_SL at sobjbar_sim.f:204-218 only —
  finite even when CHIL_MAT<0.

## Distribution of EXPONX arguments actually realized

- Flagged cells: GAM>0 by construction (crushed).
- Unflagged cells with group-summed tau in (-0.5,0): GAM negative but
  >-0.5/T1. No diagnostic records per-cell GAM; model completed 727
  iterations with no NaN/overflow evidence in OUTGEN, so realized
  negative GAM stayed in the tame regime (|GAM| < ~few).
- No Inf/NaN signatures in converged outputs; final applied change
  0.049% — consistent with stable evaluation.

## Minimum additional diagnostics if deeper proof wanted

1. In SOBJBAR_SIM, log min(GAM) per (group,depth) once per iteration —
   ~5 lines, no physics change. NOT implemented (read-only audit).
2. Extend NEG_OPAC to record the most negative UNflagged group sum
   (window -0.5 < tau < 0), which is the true residual-risk cell.
3. No diagnostic restart performed or required — existing NEG_OPAC
   already demonstrates the guard is active.
