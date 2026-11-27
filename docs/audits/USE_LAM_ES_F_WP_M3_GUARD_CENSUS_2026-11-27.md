# USE_LAM_ES=F WP-M3 - equation-modifying guards, censused

Date: 27-Nov-2026. Evidence for `checklists/use_lam_es_eq_f.md`
WP-M3: every guard that fires at scale on the
`MOM_J_CMF_NONMONO_V1` + refresh route, its file:line, what it
replaces, and whether it feeds or masks the WP-M2 mechanism
(the frozen-snapshot Jacobi map at spectral radius ~ census a).
Fixture numbers are from the WP-M1 `rfj` run
(docs/audits/USE_LAM_ES_F_WP_M2_STALL_MECHANISM_2026-11-27.md);
production counts are the diagnosis values from m15 r01/noali.

## Structural fact that frames the table

The census coefficient NU_CUP = `LAG_SCL_NM*PSI` (compressing)
or `LAG_SCL_NM*PSIPREV` (expanding) at
mom_j_cmf_nonmono_v1.f:3556-3599 is the *bounded* lagged weight.
Since `LAG_SCL = min(1, max(TB-|TA|-|TC|,0)/weight)` (:3127-3138)
and the weight contains the same PSI, the measured a-value is
capped near ~1 by construction. Rows with excess <= 0 (the
production `xs<=0` census) get LAG_SCL=0: their lagged terms are
nulled entirely (a=0, diagonal-only rows). The marginal tail the
stall lives on is the residual band - rows whose excess is small
but positive, scaled to as little as LAGB_MIN=0.0265 yet still
a~0.96.

## Census table

| Guard | Where (mom_j_cmf_nonmono_v1.f) | What it replaces | Census counters | Production | Fixture |
|---|---|---|---|---|---|
| MOM_LAG_BND (T) | :3106-3145 interior; :3285-3298, :3338-3352 row 1 | Scales the row's whole lagged source XLAG_NM by min(1, excess/weight); replaces unbounded amplification of rows where the lagged weight exceeds the M-matrix diagonal excess | MOM_LAGB_CNT, MOM_LAGB_MIN | ~245k xs<=0 rows/pass | 5489-5553 rows, min 0.0265 |
| MOM_EPS_LIM (T) | :3019-3050 interior; :3396-3401 boundary | Donor-balance cap on the +-T5 lagged EPS midpoint transfer between rows I/I+1; replaces uncapped debits that drove interior rows negative | MOM_EPS_CAP, MOM_EPS_CAP_SW | ~1E5 cells | 0 |
| MOM_XM_POS (T) | :3160-3168 | Clamps the assembled RHS XM_NM(I) to 0 when the lagged debit exceeds the local source balance; replaces a negative RHS that returns negative J from a well-conditioned solve | MOM_XMNEG, MOM_XMWORST | ~17k rows | 0 |
| Post-solve negJ rejection + runaway | :3869-3900 (significant -> XM_FB_NM, roundoff -> abs); :3953-3970 (runaway -> XM_FB_NM) | Rejects a failed solve back to the previous valid field instead of storing negative or exploded J | MOM_NEGJ_LRG, MOM_N_NEGJ_T, MOM_GS_NQUAR, MOM_QUAR_AGE, MOM_N_RUNAWAY | ~11k cells/blk call | 0 |
| MOM_GAM_FAC_MIN=1E-2 (parameter) | :597; applied :2297, :2313, :2418, :2437 | Floors |1+sigma*K| in the GAM factors; replaces the zero-GAM resonance where the implicit diagonal collapses | PASS_NRES ('gam-floor' per pass, :6334) | 10144 cells/pass | 1 of 4256 (last-freq slice) |
| Off-grid foot reads | JUP_NM/HUP_NM guards (:2847-2856, :2904-2921) + census 'offgrid' | JUP=0 when the characteristic foot lands outside 1..NCF | PASS_NOFF | F8 item (184 on m15) | 14308 reads |

## Feeding vs masking assessment

- **MOM_LAG_BND: masks the scale, not the existence.** The
  measured a (0.964 fixture, 0.995 production) is post-bound;
  the raw |PSI|/TB on the deepest rows is ~36x larger
  (a/LAG_SCL at the 0.0265 minimum). The bound caps the measured
  marginality at ~1 by construction and nulls the xs<=0 rows
  outright - it cannot remove the near-excess tail (a->1), which
  is where the stall lives. It feeds nothing into the map; it is
  the reason the map is only marginally rather than wildly
  non-contractive. WP-M5 flag: report unbounded a alongside the
  bounded value; do not remove the guard.
- **MOM_EPS_LIM: preconditioner, inactive in the reproduction.**
  Its pairwise-conservative cap bounds the same lagged weight
  LAG_BND measures. No feed or mask role in the fixture stall;
  production's ~1E5 count means it fires where debits are
  largest - plausibly the same marginal cells.
- **MOM_XM_POS: symptom-masking.** The negative RHS it clips is
  produced by the marginal lagged debits at absorbing rows; the
  clip discards the debit silently (row lands on XM=0). Inactive
  in the fixture (0 clips) because the fixture's absorbing
  outer cells keep positive local source. WP-M5 flag: census
  only.
- **Post-solve negJ/runaway rejection: downstream mask.**
  Reverting failed solves to XM_FB_NM reseeds the lagged store
  with stale values - the sweep-side comment (:4232-4243) notes
  the Jacobi write re-injected the refresh-removed error (~2M
  cells >50% between retained field and next sweep while refresh
  residuals stayed ~1E4). Inactive in the fixture because the
  bound keeps solves finite.
- **MOM_GAM_FAC_MIN: feeds the marginal band.** The floor keeps
  PSI nonzero at the 1+sigma*K -> 0 resonance; resonance cells
  are exactly the zero-excess rows LAG_BND then nulls. The two
  guards bracket the resonance from opposite sides. Production's
  10144/pass means the m15 CDS sits on the resonance; the
  fixture barely does (1 cell) - a coverage gap of the
  reproduction, not a defect.
- **MOM_RF_BULK_EXIT (early-exit, not a guard): masks the
  signature in the fixture.** The fixture's compressing drifts
  are sub-50% so the bulk gate fires at pass 4 before the stall
  floor (pass ~10). Production's higher >50% counts evidently
  kept it closed. Recorded as a fixture deviation, not a defect.

## Guards flagged for WP-M5 gating (census/reporting changes only)

1. Report unbounded a = |PSI|/NU_TB beside the bounded value in
   the nu-census line - the current amax under-reports raw
   marginality by up to ~1/LAG_SCL.
2. XM_POS clip count is already censused (MOM_XMNEG); keep it -
   it is a symptom meter, not the defect.
3. The negJ->XM_FB rejection's reseed behaviour is already
   documented in-code; no removal.

No guard is removed or weakened by this census; all findings are
reporting-level. The stall mechanism itself (WP-M2) persists
with every guard at production defaults.

## Verdict

WP-M3 VERIFIED: five equation-modifying guards censused with
file:line, replacement behaviour, production+fixture counts, and
a feeding/masking verdict each. LAG_BND is the only one that
reshapes the measured defect (caps a at ~1, nulls xs<=0 rows);
the marginal tail it leaves still exceeds the 60-pass budget.
