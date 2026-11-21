# DS-corner lag diagnostic proposal (NOT implemented - for review)

## Mechanism under test

The ordered P march at fg_j_cmf_nonmono_v1.f ~12172-12408 forms
slope/source terms for row J using the PRE-MARCH source of row
J-1 (DSD=KN_ORD_DS(J-1,L) taken before row J-1's update), the
documented lagged C-corner treatment of Knop et al. 2009 Paper II.
It is intentional, not a bug. The question is only how much of the
81P/82P residual this lag contributes. Rows 1-2 co-move with the
same sign in 96-98% of columns (phase 2.23-2.25), consistent with
coupling through J-1, but no causal fraction is measured.

## Proposed instrumentation (observational only)

In the P march loop, at the point where DSD is consumed for row J:

1. Record DSD_pre(J,L) = KN_ORD_DS(J-1,L) as used (already lagged).
2. After row J-1's update completes, record
   DSD_post(J,L) = KN_ORD_DS(J-1,L) refreshed.
3. Write delta = DSD_post - DSD_pre into a diagnostic array
   KR_DBG_DSL(NI,NCF) alongside existing KR_DBG_FLD* dumps.
4. Do NOT feed DSD_post back into the march - the solve must be
   bit-identical with the diagnostic on or off.
5. Optional shadow pass (separate flag, default F): recompute row
   J's slope/coefficient terms substituting DSD_post, apply the
   row-J update locally, and record the hypothetical intensity
   difference dI_shadow(J,L). This estimates the causal transfer
   change the lag would have produced. Shadow values must never
   enter FLD_P, KR_DBG_FLD*, residual counters, or any store.

## Quantities the diagnostic answers

- per-cell lag magnitude |DSD_post - DSD_pre| vs |DSD_pre|
- correlation between lag magnitude and per-pass defect sign/size
- |dI_shadow| vs the observed per-pass |dI| at rows 1-2 and
  propagation into rows 3-4
- whether the lag localizes to the grazing-ray boundary pair or
  appears on converged control rays (ray 84 has the same march)

## Controls

- Flag-gated (e.g. KN_DBG_DSLAG=T/F, default F); default path
  bit-identical to current source.
- Run the identical VADAT/model with flag on: production numbers
  (census, residuals, FG_RF_ALL_CONV) must be unchanged.
- Compile check on both flag settings; diff of solver objects on
  flag-off path must be nil apart from the added writes.

## Constraints

- No change to march order, C-corner treatment, tolerances, or
  any store/commit path.
- Diagnostic arrays deallocated/zeroed when flag off.
- Review required before implementation: the shadow pass touches
  solver internals even though it writes nothing back.
