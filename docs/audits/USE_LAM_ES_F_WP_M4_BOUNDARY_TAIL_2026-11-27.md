# USE_LAM_ES=F WP-M4 - boundary and tail defects (F8, F10)

Date: 27-Nov-2026. Evidence for `checklists/use_lam_es_eq_f.md`
WP-M4: the two deferred diagnosis items, each labelled
VERIFIED-OK or defect with evidence. Fixture: the committed
WP-M1 `rfj` arm (NCF=32 build) plus a new driver arg21 sweeping
the midpoint N/H closure GF.

## F8 - off-grid characteristic-foot reads get JUP=0

Question: is the free-outflow frequency boundary physical for
the CDS foot reaching the red edge (production: 184 reads,
largest matched offset 8166 of 8331)?

Structure of the read (mom_j_cmf_nonmono_v1.f):

- Compressing rows (Q_NM<0) at frequency L read the lagged
  upwind field at bin L+MOM_NM_J(I); when L+MOM_NM_J > NCF the
  read is skipped and JUP_NM stays 0 (zeroed per call at :2813,
  guards :2847-2856, :2863-2874). Census: MOM_GS_NOFF :2882
  (sweep), PASS_NOFF :6324 (refresh).
- Under MOM_FREQ_MATCH (auto-set T on this route) the offset is
  capped: N_CAP = min(NCF_IN-L, CEILING(T4)) at :2456-2460. A
  foot landing beyond the grid is clamped INTO the last bin -
  reads stay in-grid by construction. The off-grid population is
  therefore exactly the L=NCF rows: the match block's guard
  (NCF_IN-L >= 1) fails at the last bin and every compressing
  row there gets JUP=0. The fixture nm dump confirms NM_J=0 at
  all 445 compressing cells for the last-frequency solve; the
  refresh-side augres census counted 14308 off-grid reads =
  ~445 cells x ~32 pass evals, i.e. all of them at the red edge.

So the boundary is two-layered: (a) feet beyond the band are
clamped to the edge bin's J as the upwind value (an
approximation: the edge bin stands in for all un-modelled
redward flux), and (b) the edge bin itself is solved with zero
compressing inflow - free outflow, the only consistent choice
without a band-edge model (persistence or extrapolation would
fabricate flux).

Verdict: VERIFIED-OK for the stall, with a coverage caveat.

- Mechanism sign: JUP=0 removes the lagged term from the row's
  RHS entirely - the off-grid row's a-value drops to 0, i.e. it
  acts OPPOSITE to the non-contraction, never feeding it. It
  cannot produce or amplify the WP-M2 rho~1 map; it is a pure
  source truncation that under-prices the compressing field at
  the redmost bin.
- Count: 184 production reads ~ one frequency's worth of CDS
  compressing cells - negligible against ~1.6M rows.
- Caveat flagged for WP-M5/M6: the clamp-into-edge rule means
  every compressing cell whose physical foot overshoots the
  band reads a J that was itself built with zero compressing
  inflow. If m15's CDS feet saturate at the edge across many
  frequencies (the 8166-bin offset suggests the CDS's Doppler
  reach exceeds the band), the band-edge J is systematically
  under-priced at exactly the marginal cells. Whether m15's
  ln-nu span is wide enough is a model-data question the
  fixture cannot answer; a production census of clamped
  (N_CAP=NCF-L) reads per frequency would settle it.

## F10 - closures from expansion-biased V13 at sigma<-1

Question: does the V13 closure misprice the moment rows where
the march stalls?

Provenance: DO_NONMONO_SOB + USE_LAM_ES=F keeps FG_J_CMF_V13 as
the FG solver (comp_j_blank.f:928) -> GET_MOMS_NON_REL ->
MOD_RAY_MOM_STORE; V13 is expansion-biased at sigma<-1
(comp_j_blank.f:891), and the degenerate |H|<1E-3|N| branch
fills NMID_ON_HMID = SIGN(1E10,N) (get_moms_non_rel.f:216). The
closure enters the row coefficients through 1+sigma*g, g =
NMID_ON_HMID.

Test (driver arg21, committed fixture): sweep the midpoint
closure GF over five decades and check the stall signature.

| GF | compressing | expanding | refresh | census amax |
|----|-------------|-----------|---------|-------------|
| 0.5 (committed) | STALL pass ~16 | survives | converged | 0.964 |
| 0.1 | STALL | survives | converged | 0.773 |
| 2.0 | STALL | survives | converged | 0.504 |
| 10.0 | STALL | survives | converged | 0.499 |
| -0.5 | no stall | survives | converged | 0.964 |
| 1E10 (degenerate branch) | no stall | survives | converged | 0.499 |

Read: the closure value modulates the marginality - amax falls
from 0.96 to 0.50 as |GF| grows (the sigma*g weight raises the
diagonal excess the same way MOM_H_IMPL does) - but the stall
signature is invariant over GF in [0.1, 10], the plausible
mispricing range. Only regime-changing values remove it: GF<0
flips the sign structure of the compressing rows and GF=1E10
saturates 1+sigma*g so the lagged share is negligible - both
are different physics, not a mispriced closure in the same
regime.

Verdict: VERIFIED-OK for the stall mechanism. The V13
expansion bias is real provenance damage (the closure feeding
the marginal rows is measured by a solver built for the
opposite regime) and it does move the row coefficients - but
the non-contractive channel is the frozen-snapshot lagged J
read (WP-M2), which the closure only modulates, not creates.
A mispriced g inside its plausible range cannot rescue the map
or explain the stall; it shifts the fixed point and the exact
a-values. Keep it flagged as wrong-regime input for the fixed
point's accuracy; it is not the gate the fix must pass.

## WP-M4 exit

- F8: VERIFIED-OK (correct free-outflow boundary; acts opposite
  to the defect; small count). Coverage caveat recorded for the
  band-edge census on the production arm (WP-M6).
- F10: VERIFIED-OK for the stall (signature invariant over
  GF 0.1-10; amax only modulated). Wrong-regime provenance
  stands as an accuracy concern, not the contraction defect.
