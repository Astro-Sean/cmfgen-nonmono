# EXPONX negative-opacity audit — master report

Commit audited: `a029fd37` (main). Read-only; no source/model changes.

## Mechanism summary

EXPONX has no negative-X guard, but negative GAM's only route to it is
the group-SUMMED CHIL through SOBJBAR_SIM/SOBEW. Production has a
deliberate, verified upstream guard (CHECK_LINE_OPAC, default TRUE):

  summed angle-free Sobolev tau < -0.5  =>  NEG_OPACITY(I)=T,
  CHIL->1.0, every CHIL_MAT member -> 1/NUM_SIM_LINES.

At flagged cells GAM>0 always. Per-line inversions inside non-flagged
groups never reach EXPONX individually — they only enter the
apportionment at sobjbar_sim.f:204-218 (finite for any sign).

## The residual latent window

A group whose summed angle-free tau lies in (-0.5, 0) is NOT crushed.
Then GAM=GAMH/T1 with T1>=1e-4 gives GAM in (-5000, 0) — crossing
EXPONX's -709.78 overflow point. In this model it is unrealized (all
159 flagged cells were below -0.5 and crushed; run completed), but the
structure permits it: one sufficiently negative group that does not trip
the -0.5 threshold... — precisely: a group with summed tau in (-0.5,0)
AND a near-singular direction T1~1e-4. That combination produces
GAM ~ -5e3 -> EXPONX=+Inf -> BETA=Inf -> likely NaN/Inf propagation.
No current diagnostic records whether such cells occur.

## Confirmed vs theoretical

CONFIRMED (P0): guard exists, is active (ESEC_CHK + CHECK_LINE_OPAC,
VADAT + defaults + NEG_OPAC log), prevents the flagged cells reaching
EXPONX negative; sob_rh_couple explicitly guarded; JFNK replay and obs
code duplicate the crush; EXPONX itself mathematically exact.

LATENT (P2, not a defect in this model): the (-0.5,0) window is an
unbounded negative-GAM channel into EXPONX; d_EXPONX_dX shares it.
Also the apportionment ST_ON_SL can go negative (CHIL_MAT<0, CHIL>0),
giving ZNET extrapolation of mixed sign — finite, physically odd,
never flagged.

PHYSICS GAP (P2): crushed cells lose real weak-maser gain (high-n
recombination inversions at d63-70, Ne~1e13); unflagged mild inversions
get unsaturated-gain EXPONX>1 with no saturation or multi-resonance
coupling (DO_NONLOC_SOB=F). Both are approximations, neither is a bug.

## Literature check

- DAH15 (MNRAS 449, 4304): floor on |1+sigma mu^2| — present verbatim;
  says nothing about negative CHIL.
- RH78 (ApJ 219, 654): Sobolev multi-resonance coupling — related to
  sob_rh_couple's double-resonance term; not an inversion treatment.
- Lamers, Cerruti-Sola & Perinotto 1987 (ApJS 71, 983): standard
  tau_S form — the GAMH/|1+sigma mu^2| convention matches.
- 1984 JQSRT transfer-with-inversion method / Gray (2012) maser
  modelling: the field-standard answer is unsaturated gain + explicit
  saturation limit; production has the gain expression but no
  saturation model — the crush is the de facto ceiling. Not a
  validation of the crush as physics.

## Defect classification

Not a demonstrated defect (P0-verified safe in this model). Classified
as a latent numerical vulnerability + unresolved physical approximation
(P2), with a P3 diagnostic gap (unflagged negative group sums are
invisible to current logging).
