# CMFGEN KNOP provisional-store population-consumption audit (Phase 2.16)

Status: source-level lifecycle audit complete; the population A/B
rate test and the spectral comparisons still need live model runs
(same limitation as Phase 2.15). Provenance: HEAD 9964bcd6.

## 1. Complete store lifecycle (data flow)

    FG_J_CMF_NONMONO_V1 (per-frequency ray solve)
        |
        +--> ray intensities I_F -> moments -> RJ (population J)
        |
        v
    KNOP_RAY_SOLVE_V1 per ray LS (refresh, NPASS passes)
        exit class: 0 exact / 1 tol / 2 bulk / 3 stall / 4 cap /
                    5 div / 6 plat / -1 unsolved / -2 build fail
        |
        +-- converged channel: WRITE_STORE
        |     STORE = T2 + INEXT_RELAX*(T1-T2)   (blended)
        +-- fail-close channel: WRITE_STORE_CAP (FG_COMMIT_CONV)
              V1 capped at FG_COMMIT_NC*max(V2, 1e-4*VMAX_F)
              (per-cell ratchet, no blend)
        |
        v
    STORE = I_P_DEPTH_NEXT / I_M_DEPTH_NEXT (module SAVE arrays)
        |
        +--> SEED_FLD at next refresh pass (same ray)
        +--> upwind coupling I_*_DEPTH_NEXT read at l+/-1 during
        |    every subsequent frequency solve (lines ~5163-5233)
        |    -> next solve's I_F -> moments -> RJ -> rates
        +--> FGINEXT_GS checkpoint file (fort.291 read/write,
             lines ~7131-7346), validated for NaN/negative

    NXIT -> NXT_LS -> NX_CNT histogram -> LUER print ONLY.
    No per-ray status travels with the store or out of the refresh.

Consumers: SEED_FLD (refresh seed), the upwind-coupling reads in
the frequency march (every ray's solve), the checkpoint
persistence. All three read raw store values with no
provisional/converged tag. **No consumer can distinguish a
provisional store from a converged one.**

## 2. Population consumption and certification

- RJ (population J) is produced inside the formal solve each
  frequency, from intensities whose compressing-branch upwind is
  I_M_DEPTH_NEXT (line 5163ff: "Compressing (GAM<0): upwind =
  next frequency = I_M_DEPTH_NEXT"). A provisional store value is
  consumed as upwind data exactly as a converged one is.
- Rates: RJ enters SOURCE (2083), NUST/VB/VC (2101-2106), STEQ_T
  (2683-2704) in cmfgen_sub.f.
- Population certification = MAXCH (population demand) under the
  mode controller; the field side is the GS maturity gate:
  `GS_MAX_ALL = MAX(GS_MAXCHG, GS_SW_MAXCHG, GS_RF_DCHG)` vs
  `GS_MATURE_TOL = 0.05` (solve_for_pops.f:213,1248-1285).
- There is no FG_RF_CONV flag. GET_FG_GS_MATURITY exports GS_OK
  (I_NEXT_GS_OK checkpoint flag), MAXCHG, counts, RF_DCHG -- but
  nothing equivalent to the MOM-side MOM_RF_CONV ("refresh
  converged this iteration"). NXIT classes die in the histogram.

## 3. The Phase 2.15 five-step sequence is possible

1. Ray exits bulk (IXIT=2) under KN_RAY_BULK_TOL>0.
2. WRITE_STORE_CAP commits the capped field to the store.
3. Next frequency solve reads it as upwind; moments->RJ->rates.
4. Populations iterate on that field.
5. Field maturity is judged by GS_MAX_ALL vs 0.05 -- and a
   bulk-converged refresh at residual ~0.03 is BELOW 0.05:
   `GS_FG_IMMATURE=.FALSE.` reads the field as MATURE.

The 5-step chain is structurally real: a provisional store feeds
rates, and the maturity gate cannot see it because the floor
(6e-2) and the observed bulk residual (~3e-2) sit inside the
maturity band.

## 4. The m15 non-convergence mechanism

Evidence chain from the FRESH_knop OUTGEN plus the source:

- Refresh ends with bulk exits at maxres 1.8e-2 - 3.1e-2
  (9-15x KN_RAY_TOL=2e-3).
- GS_MAX_ALL ~0.03 < GS_MATURE_TOL=0.05 -> field reads mature
  -> freeze releases.
- The released populations run on a field 15x above ray
  tolerance; demand (MAXCH) does not fall 10% below release
  demand, or the sweep drifts above GS_UNSAFE_TOL=0.5 ->
  release counted no-progress (N_REL_NOPROG++), refreeze.
- After GS_REL_NOPROG_LIM=10 such cycles: STOP
  "release/refreeze cycle never matured the lagged field;
  stopping as non-convergent".

So the observed death is a freeze/release limit cycle driven by
a maturity band that admits bulk-floor residuals. FG_COMMIT_NC
supplies the provisional field; the maturity band (not the cap
itself) is what let it certify. This is consistent with - but
not yet a proof of - causation; the alternative (a genuinely
non-contracting field) is not excluded by this evidence.

## 5. Provisional vs converged vs invalid

- FULL_CONVERGED_STORE: all rays ix it 0/1 (exact/tol); today
  indistinguishable from a bulk store once committed.
- PROVISIONAL_PROGRESS_STORE: bulk/stall/cap/plat exit; kept for
  refresh progress (WRITE_STORE_CAP ratchet). Must not be the
  final certified field.
- INVALID_STORE: NaN/<=0/failed build; fort.291 read-back
  validation exists (lines 7152-7155) and fails closed.

Smallest state mechanism: a per-refresh LOGICAL
FG_RF_ALL_CONVERGED (MOM_RF_CONV analogue), set =ALL(NXT_LS<=1)
over the production rays at the census block (~9125), exported
through GET_FG_GS_MATURITY. That is one flag, no store-level
tagging, and it survives to the maturity gate where it is
needed. A per-ray array is not required because maturity is
judged refresh-wide.

## 6. Final-convergence policy options

- A (provisional for init only): require all-ray ix it<=1 before
  the field can certify mature. Simplest; strictest; may park
  marginal rays (the parked-store tradeoff the cap exists for).
- B (final verification pass): allow provisional stores during
  iteration; before final population convergence is reported,
  run one refresh and require all-ray FULL; re-check MAXCH.
  Adds one refresh per certification attempt.
- C (hybrid): tighten the maturity band as MAXCH falls
  (GS_MATURE_TOL scaled to KN_RAY_TOL when MAXCH near tol).
  Cheapest; risk of oscillating classification near the band.

Recommendation (audit-level, not implemented): A or B. C is a
band rescale, not a contract; the clean fix is an explicit
all-ray convergence flag feeding the maturity gate plus a final
verified refresh before the run reports converged - i.e. B
with A's flag.

## 7. Velocity-inversion model search

No deck exists. All available models (m15_2p82d_FRESH_knop,
handoff_smoke/knop_mom, tmp/v13_audit) carry the same hydro:
monotonic V, SIGMA -1.04 to -156.9, all 70 cells compressing.
Needed: an SN_HYDRO_DATA whose V(r) contains a genuine local
non-monotonicity or a decelerating shelf that places seeded
cells on both sides of GAM=0 - e.g. a shock-flattened or
re-accelerating profile from a radiation-hydro interaction run.
Synthetic sign flips were already ruled out.

## 8. Safety/regression

No source changes; knray_p214 frozen; defaults, controller,
equations, observer, populations untouched. Diagnostics only.
