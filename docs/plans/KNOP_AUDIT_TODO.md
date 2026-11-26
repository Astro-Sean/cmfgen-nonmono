# KNOP_AUDIT_TODO.md - independent audit of DO_NONMONO_KNOP

Baseline: 53ee4e14 (test_knop_lagged_field pin commit)
Audit target: `cmfgen/plane/fg_j_cmf_nonmono_v1.f` and the refresh stack
(KNOP_RAY_SOLVE_V1 / BUILD_FREQ / MATCH_FREQ / SRC_PRIME / SLOPES /
MARCH / SWEEP_FLD / SWEEP_FLD_ORD / WRITE_STORE).

Priority order per the brief: P1 incorrect transfer physics or invalid
radiation fields; P2 solver dispatch, numerical robustness and
convergence integrity; P3 performance, diagnostics, maintainability.

Status key: FIXED / VERIFIED / OPEN / BLOCKED.

-------------------------------------------------------------------------------

## Pillar 1 - inward ordered-sweep fallback dispatch

Claimed bug: `ELSE` branch of the inward I_M dispatch calling
SWEEP_FLD_ORD instead of SWEEP_FLD.

Finding: NOT PRESENT at baseline. Commit 1b393a5a (17-Oct-2026) already
corrected both dispatches; at 53ee4e14 the I_M block reads
`IF(ORD_SWEEP) CALL SWEEP_FLD_ORD ELSE CALL SWEEP_FLD` and the I_P block
is identical in structure. git blame confirms the correct calls were
never reverted.

Disable conditions audited:

- `ORD_SWEEP = FG_ORD_SWEEP .AND. HAVE_PACK .AND. NI >= 3` (single
  expression covering flag, pack availability and slope width).
- SPW allocation failure: `IF(IOS.NE.0)ORD_SWEEP=.FALSE.` after
  `ALLOCATE(SPW(NCF,4))` - falls back correctly.
- HAVE_PACK=F: pack allocation failure is cleaned up (every PK_* array
  deallocated individually) and the alternating path never reads them:
  all `PK_*` references inside SWEEP_FLD sit within one
  `IF(HAVE_PACK)THEN ... ELSE CALL BUILD_FREQ` block.
- NI<3: slope endpoints need three rows; the `NI.GE.3` term in the
  ORD_SWEEP expression excludes it. NI=2 is legal for the alternating
  sweep (SLOPES and MARCH index only I and I+1); NI<2 returns early.
- SPW is referenced only inside SWEEP_FLD_ORD and its KN_ORD_* helpers;
  the fallback path contains no SPW access.

Status: VERIFIED at baseline; regression pin added (no code change
required): `test_ord_sweep_dispatch_fallback` in
unit_tests/python/test_knop_lagged_field.py asserts the dispatch
structure, all four disable conditions, and the no-unallocated-access
property. NI=2 is covered by the NI.GE.3 term pin.

-------------------------------------------------------------------------------

## Pillar 2 - missing I_NEXT donor physics (sink defect) - FIXED

Defect confirmed: a compressing ray point (GAM<0) whose matched upwind
donor was unavailable kept `+|GAM|/DLN` inside chi_hat while the
modified source dropped the `|Q|*I_NEXT` inflow. The implicit outflow
stayed in the operator with no matching explicit inflow: a sink the
discretised equation never wrote down. The missing-donor source branch
already intended to "drop the frequency coupling entirely"; the chi_hat
side was simply missed.

Missing-donor taxonomy (audited against mask construction and the
NM_RAY/NFM_L caps):

1. Genuine frequency boundary: `FREQ_IDX_ITER+NM_RAY > NCF_IN` in the
   sweep, `L+NFM_L(I) > NCF` in the refresh. Unreachable while the
   matched offsets are capped on-grid (both are), but kept as the
   boundary rule: outflow term stays, inflow is exactly zero.
2. Missing sparse-storage depth: GAM<0 at a ray point whose two
   bracketing model depths carry no COMP_DEPTH_IDX slot. This is the
   real mask-gap case - cubic SIGMA_RAY interpolation overshoots
   sigma=-1 between unflagged depths, or a marginal compression the
   +/-1 mask expansion did not reach.
3. `HAS_COMPRESSING=F`: no sparse store exists at all; every
   compressing overshoot cell previously kept the sink. Now every such
   cell drops the coupling.
4. Boundary ray points beyond the model depth grid (K<1 or K+1>ND):
   same treatment - the outer thick-boundary extension cells are
   documented as monotonic, so dropping the coupling there is
   consistent with their existing handling.

Fix implemented (both paths, one gate each):

- Sweep: new file-scope `LOGICAL FUNCTION KN_COUPLING_OK(I2,LS,NM2,
  NCF2,ND2)`. Returns T for expanding cells and for genuine red-edge
  boundaries; for interior compressing cells it requires a
  COMP_DEPTH_IDX slot on either bracketing depth. Applied to the +|GAM|/
  DLN chi_hat increment (INIT and general branches), the
  SIGN*dGAMdR/DLN chi_hat derivative, and the signed Q = GAM/(DLN*chi)
  term; dropped cells get Q=0 explicitly. Because the SP-build donor
  read is the same predicate, both sides of the coupling are now
  present or absent together; the FG_NOSLOT census still counts exactly
  the cells whose coupling was dropped (comment updated).
  NM_RAY is OMP-private in the ray loop, so the caller passes the
  private value; ND/NCF_IN are routine arguments, passed likewise.
- Refresh: new contained `LOGICAL FUNCTION KN_RF_COUPLING(I2,L2)`
  mirroring the same rule via the IKC sparse-slot mask (the refresh's
  donor convention: SRC_PRIME needs IKC1|IKC2 and L+NFM<=NCF). Applied
  identically inside BUILD_FREQ (L=1 and L>1 branches).

Consequence noted and accepted: when the whole coupling is dropped the
cell reverts to the plain chi operator - i.e. the field is taken
locally flat in frequency at that point. This is the only consistent
closure without a donor; the alternative (outflow-only) is exactly the
sink being removed. A self-donor (NM=0 self-slot at the last frequency)
produces the same algebra, which supports the choice.

Numerical verification: drv_knop_ray cases `gap` / `gapr` construct a
compressing band with the mask covering depths 3-4 only. In `gap`,
cell 5 keeps GAM<0 but has no donor, so its whole coupling is dropped;
in `gapr` cell 5 has GAM=0 exactly, producing the identical discrete
operator (plain chi, Q=0, no upwind term, dGAMdR=0). The two committed
stores must agree to solve tolerance - test T11 in
unit_tests/python/test_knop_ray.py.

FG_NOSLOT role: diagnostic only. It is summed into FG_NOSLOT_CNT once
per frequency and reported via GET_FG_NOSLOT -> KN_NOSLOT line in
cmfgen_sub.f; it does not feed the field-maturity decision
(GS_NBAD/GS_NTOT and the residual do).

Status: FIXED (this session) + static pins
(test_missing_donor_consistent_drop) + numerical pin (T11).

-------------------------------------------------------------------------------

## Pillar 3 - finite-value protection standardised + fail-closed build - FIXED

Defects found in BUILD_FREQ (refresh coefficient builder):

- The L>1 branch used NaN-only tests (`x.NE.x`) on CHI_L, ETA_L, GAM,
  GAM_MAG, dGAMdR, DCHR_L, Q_L, SR_L and on the DTAU fallback. +Inf
  passed every one of them: an Inf DTAU then marched a wrong-but-finite
  field through the A0-A4 branch (T1>=40 -> EEL=0, finite coefficients)
  and could be committed.
- The DLN spacing was never floored in the L>1 loops, unlike the
  in-sweep chi_hat loop which floors DLN_C at 1.0E-30 before division.
- The in-sweep dCHIdR and signed-Q loops divided by DLN_C without the
  floor their own chi_hat loop applied - a 0/negative DLN_RAY could
  poison the derivative while chi_hat itself was guarded.

Changes:

- All BUILD_FREQ guards (L=1 and L>1) now use the
  `.NOT.(ABS(x).LE.HUGE(x))` finite test, matching the in-sweep and the
  refresh L=1 branch convention. Sub-floor-but-finite values still take
  the documented clamps.
- DLN floors added to every division site in BUILD_FREQ (matched the
  in-sweep floor form) and to the two previously unfloored in-sweep
  loops (dCHIdR update, signed Q).
- Fail-closed propagation: any guard that detects a NON-FINITE value,
  or a non-positive matched span, sets the host-scope DIVERGED flag.
  The pack-build loop exits on it and the ray takes the diverged-ray
  exit (`NDIV=1`, GOTO 800 -> heal, RESID=1, NXIT=5). The no-pack
  SWEEP_FLD exits the frequency loop and the pass loop with IXIT=5.
  A failed coefficient construction therefore propagates as a failed
  ray solve and can never be accepted as a converged refresh; the
  maturity gate reads RESID=1 as unconverged. The floor/zero values are
  still written so the abandoned local solve stays finite.
- A post-loop finite check on A0L-A4L was added as the last belt: every
  march coefficient is a pure function of the (now finite) DTAU/E set,
  so a non-finite there can only come from an overflow the per-term
  guards missed.
- MARCH and WRITE_STORE/WRITE_STORE_CAP already handled non-finite
  field values (NaN/zero clamps on VNEW, VMAX cap path); no change
  needed there.

Note: guards that fire on legitimately sub-floor finite values (e.g.
CHI_L < 1.0E-30 from thin boundary cells) do NOT fail the ray - only
non-finite values do, matching the "invalid value = non-finite"
boundary the audit brief drew.

Numerical verification: drv_knop_ray cases `nanchi` (NaN in
KN_CHI_EXT(4,3)) and `infgam` (+Inf GAM(4)) must now exit NXIT=5,
NDIV=1, RESID=1.0 - tests T12/T13.

Status: FIXED (this session) + static pins
(test_build_freq_fail_closed) + numerical pins (T12, T13).

-------------------------------------------------------------------------------

## Pillar 4 - in-sweep / refresh operator equivalence - ANALYSED

Coefficient-by-coefficient comparison performed:

- chi_hat: identical structure - chi + |GAM|/DLN where DLN is dLOG_NU
  for expanding and the matched span for compressing cells. The sweep
  divides by DLN_RAY after flooring; the refresh now does the same.
- The compressing-branch replacements (redder upwind, |Q| inflow) use
  the same matched offset/span pair (NM_RAY/DLN_RAY vs NFM_L/DLNM_L)
  and the same bracket-depth interpolation rule (single-donor fallback
  included) in both paths.
- DCHR: both add SIGN(1,GAM)*dGAMdR/DLN to d(chi)/dR - sign verified:
  for GAM<0 the |GAM| term's derivative is -dGAMdR, matching the
  compressing form -dGAMdR/DLN in the INIT branch.
- DTAU: identical trapezoid + (dchi/dR)/6 correction, with the
  min/max clamp applied under the same rule (sweep skips it for
  INSERT and MON methods; refresh keys on the captured KN_METHOD_SAV/
  KN_INSERT_SAV - verified consistent).
- Q, SOURCE_PRIME, A0-A4 (Eddington factor block): same forms.
- Frequency indexing verified end-to-end: FREQ_IDX_ITER is 0-based per
  frequency, the sparse store is slot-per-frequency with slot s holding
  frequency s+1, SEED_FLD/WRITE_STORE agree, and the last (phantom)
  slot is never written - consistent on both paths.
- Store-slot conventions verified against SEED_FLD
  (FLD(:,LS2)=STORE(:,LS2-1)) and the in-sweep capture at
  I_*_DEPTH_NEXT(...,FREQ_IDX_ITER).
- Ordered-path donor choice: KN_ORD_SPV recomputes the matched offset
  inline (NNF2 = MAX(1,CEILING(DRF_L/KN_DLN_N)) capped at NCF-L)
  rather than reading NFM_L. The formula is identical to MATCH_FREQ,
  so the ordered and alternating paths select the same donor; PK_Q
  already carries the DLNM_L span. Duplication noted for P3
  maintainability - a single shared match helper would remove the
  second copy.
- Expanding-cell coupling verified consistent: sweep Q = GAM/
  (dLOG_NU*chi_hat) with SOURCE_PRIME += Q*I_M_PREV mirrors refresh
  Q_L = GAM/(DLNM_L*CHI_L) with SP_L += Q_L*FLD(I,L-1); DLNM_L for an
  expanding cell is the adjacent-bin span (uniform grid: dLOG_NU).

Remaining gap (documented, not yet exercised): the refresh rebuilds
CHI/ETA from the captured monotone cubics at the ray points, while the
sweep interpolates the live arrays - the two agree exactly at grid
nodes but can differ by the cubic residual between them. That is
inherent to the captured-store design, not a bug. A strict
coefficient-level numerical A/B (same inputs, same outputs) is still
worthwhile - see Pillar 7 test 10; the manufactured-driver cases
already exercise full refresh equivalence (T8) and both sweep
strategies (T3).

Status: VERIFIED analytically; numerical A/B test listed under P7.

-------------------------------------------------------------------------------

## Pillar 5 - convergence criteria - ANALYSED

Gate chain audited:

- Per-ray exits: NXIT 0 (exact, no compressing cells), 1 (residual
  tol), 2 (bulk bad-fraction), 3 (same-direction stall), 4 (pass cap),
  5 (diverged), 6 (hot-pass cap). NXIT>=0-2 with residual<=KN_RAY_TOL
  (or bulk exit under KN_RAY_BULK_TOL) is a commit candidate;
  stall/cap/diverged exits never commit (FG_COMMIT_CONV gate at
  commit).
- Diverged rays report RESID=1.0 so the max-norm cannot read them as
  converged; they are counted in FG_GS_NDIV and share-gated in
  SOLVE_FOR_POPS (W2 17-Oct-2026).
- The best-pass census (NBAD_BLK) keeps stall exits honest: the gate
  reads MIN(GS_NBAD, GS_NBAD_BLK), so a bulk exit cannot hide the bad
  cells a stall left above tolerance.
- Capped writes (FG_INEXT_CAPPED) and clamp events feed the census;
  missing donors now feed FG_NOSLOT (diagnostic) and no longer create
  sink cells.
- New (this audit): a failed coefficient build exits as diverged
  (RESID=1) - can never satisfy the maturity gate.
- Population-release hold: the solve holds on the GS/KNOP maturity
  verdict; a bulk exit therefore cannot "pass" the field - the bad-cell
  census still counts its leftovers.

Residual risk noted for the vault: a ray that exits NXIT=2 (bulk) with
a large-but-sub-tolerance residual commits its field; the hold gate
reads MIN of last/best bad census, so a bulk-only pass still counts
those cells - accepted behaviour by design.

Status: VERIFIED analytically; no code change needed.

-------------------------------------------------------------------------------

## Pillar 6 - transfer-equation physics - ANALYSED

Derivation check of the implemented operator:

- chi_hat = chi + |GAM|/dln(nu) is the standard frequency-advection
  term in the comoving-frame formal integral (Knop et al. 2009 sign
  convention). For GAM>=0 the upwind is the previous (bluer) frequency
  and the source carries +Q*I_PREV with signed Q>0; for GAM<0 the
  comoving direction reverses, the upwind becomes the next (redder)
  frequency, and the source must carry |Q|*I_NEXT - verified in both
  sweep directions (inward I_M and outward I_P) and in the refresh
  SRC_PRIME. Using signed Q there would subtract the donor a second
  time (the opacity term already carries +|GAM|/DLN); confirmed
  consistent after the P2 gate.
- Negative radial velocity and gradient reversals: GAM carries the
  sign of (1+sigma*mu^2); the donor/upwind choice keys on GAM's sign
  only, so a V<0 flow with positive GAM and an ordinary V>0 expansion
  take the same branch - verified against the rev driver case (two
  bands, two reversals).
- The deliberate omission of the 3*a_nu*I_nu proportional term is
  inherited from the CMFGEN V13 frequency-form convention (the S_tilde
  removal comment documents this). The omitted term couples the same
  frequency's intensity back into the source; in the discrete form used
  here it would modify the effective opacity by a uniform factor rather
  than change the upwind structure, so the approximation is
  self-consistent with the moment equations CMFGEN solves (same
  convention). Flagged as a verified approximation, not an
  inconsistency; a controlled A/B against the full term remains
  possible future work.

Status: VERIFIED (approximation documented); no correction required.

-------------------------------------------------------------------------------

## Pillar 7 - required regression tests

Implemented this session:

- Static pins (test_knop_lagged_field.py):
  test_ord_sweep_dispatch_fallback (P1: both directions, all disable
  conditions, no unallocated access), test_missing_donor_consistent_drop
  (P2: both gates, both chi/DCHR/Q terms, boundary rule), 
  test_build_freq_fail_closed (P3: no NaN-only guards, DIVERGED flag
  count, abort sites).
- Numerical pins (test_knop_ray.py + drv_knop_ray.f90):
  T11 gap/gapr dropped-coupling equivalence (P2 test 7 missing donor);
  T12 NaN opacity -> NXIT=5/RESID=1 (P3 failure handling);
  T13 Inf GAM -> NXIT=5/RESID=1 (P3 failure handling).
  Pre-existing coverage mapped: monotonic-limit equivalence (T1/ord
  cases), single compressing zone upwind (T1), multiple compression
  zones (T9 rev), ordered/fallback equivalence (T3), refresh
  equivalence end-to-end (T8), cap/failure exits (T5), checkpoint (T7).

Still OPEN (recommended next):
- Velocity reversal through GAM=0 (P7.4: partially covered by rev;
  a dedicated zero-crossing case would pin the sign transition).
- Blue/red frequency-boundary numerical pins (P7.5 - static pins only
  today).
- Discrete photon/energy balance (P7.11) - requires a conservation
  checker on the manufactured problem.
- Strict coefficient-level A/B between in-sweep and refresh on
  identical inputs (P7.10 / Pillar 4 residual item).

Status: PARTIAL - core numerical cases landed; remainder listed.

-------------------------------------------------------------------------------



===============================================================================
## PHASE 2 - NUMERICAL VALIDATION (tasks 1-3, 5)

New harness: drv_knop_ray 'ab'/'abx' cases (contained DRV_AB). The
driver feeds the identical manufactured problem (ND=8 core ray,
NCF=6, single compressing band, VDOP refinement dilated so NI=ND)
through the production in-sweep solver FG_J_CMF_NONMONO_V1 for six
outer iterations, then calls FG_KNOP_REFRESH_V1 on the captured KN_*
store. Debug dump arrays KN_DBG_* (sweep) and KR_DBG_* (refresh),
allocation-gated in FG_J_CMF_NONMONO_MOD_V1, expose the discrete
operator per (cell, frequency): CHI_hat, Q, DCHR, DTAU, S' (both
directions), W/WINX, DLNM, NM/NFM, KC/IKC brackets, A0-A4, and the
committed FLD fields. 'ab' runs the ordered refresh (pack + SPW);
'abx' forces alternating SWEEP_FLD. Tests: T14 (both), T15 (GAM
crossing). Numerical pins in test_knop_ray.py, 16/16 pass.

### Task 1 - coupling gates (DONE)

The sweep (KN_COUPLING_OK), ordered sweep (KN_ORD_SPV), and refresh
(KN_RF_COUPLING) make identical coupling-availability decisions:

- Frequency offset map: NM_RAY == NFM_L at every (cell, freq)
  including the red-edge cell (5,6): nm=1, nfm=1 (matched-offset
  identical; NNF2 in KN_ORD_SPV recomputes the same CEILING span -
  duplication retained, demonstrated equivalent).
- Donor availability: (IKC==1) == (KC>=1) everywhere. Sweep KC are
  sparse-store slots, refresh JPR are ray rows via J_PNT of the same
  bracket ends (K=RAY_PNT-ND_ADD, same code); on the production grid
  J_PNT is non-identity (depth d -> row d+1) so raw indices differ
  while resolving to identical model depths. WINX (interpolation
  weight) and DLNM (matched log span) agree bitwise.
- Red-edge closure: at L=NCF the sweep reads the force-zeroed
  phantom store slot with NM=1/DLN=last interval -> pure outflow;
  the refresh default NFM_L=1/DLNM=KN_DLN_P(NCF) (= same interval)
  and SRC_PRIME L+NFM>NCF -> SP=SR produce the identical operator.
  Verified to 3.7e-8.

### Task 2 - GAM zero crossings (DONE)

Six profiles (all-positive, all-negative, pos->neg, neg->pos,
adjacent-cell zero, continuous glide through zero): all converge,
NDIV=0, NNEG=0, NCAP=0, NSTORE=30. gpp/gpn0 (GAM(4)=+/-1e-9*GAMP
across the sign flip) commit BIT-IDENTICAL fields - the discrete
operator is exactly continuous through GAM=0; the sign flip only
selects the upwind direction. Exact-zero (gzz) differs from the
microstep at 3.8e-3, inside the fixed-point tolerance band (the
differing source term is O(1e-10)); not an intensity jump.

### Task 3 - frequency boundaries (DONE at the A/B level)

- Red edge: verified numerically above (phantom-donor pure outflow,
  sweep and refresh identical, 3.7e-8).
- Blue edge (L=1): expanding cells carry Q=0 and no bluer donor in
  both paths; all L=1 coefficients agree to the same <=2.4e-7 as
  the interior. Compressing cells at L=1 keep their redder donor -
  no coupling silently lost (identified in the DBS/DBR dumps).
- Multiple compression regions / reversals: 'rev' case (two bands,
  sign flips twice) converges (T9); 'gcr' resolved glide through
  zero converges (T15).

### Task 5 - formal/refresh operator equivalence (DONE)
Identical inputs both paths. Field-independent coefficients:
- chi_hat, Q: <=5.0e-8 relative (span evaluation order:
  LOG(KN_FREQ ratio) vs K2*DLN_N_SAV).
- DCHR: <=2.4e-7. DTAU: identical on all interior intervals; row
  NI carries the boundary value in the sweep only (bookkeeping
  slot, not part of the march).
- A0-A4 (march coefficients, both direction conventions verified
  identical roles): <=7e-8 except A0 <=2.0e-6 (EXP(-T1)
  amplification of the <=5e-8 DTAU difference).
- Modified source S': <=4.9e-5 - bounded by the lagged-iterate
  state difference (sweep store at outer-iteration 6 vs refresh
  converged to 8e-6), not an operator defect.
- Same numbers under ordered (pack/SPW) and alternating sweeps.

REGRESSION FOUND AND REVERTED: the earlier phase-2 red-edge edits
(MATCH_FREQ NFM_L=0/DLNM_L=0 'own-bin clamp', KN_ORD_SPV L<NCF
guard removal, NFM_L=0 DLN-floor exemptions) were wrong in both
directions. The A/B caught them: NFM_L=0 floored DLN to 1e-30 and
gave chi_hat ~1e23 / dtau ~4e26 at (5,6) - an artificial sink far
worse than the original code - while KN_ORD_SPV's NNF2=0 selected
an own-frequency inflow donor instead of pure outflow. The sweep's
actual convention is NM=1 + zeroed phantom slot + last-interval
DLN, which the ORIGINAL refresh defaults already reproduced. All
three edits reverted; T14 now pins the correct closure.

## Follow-up (validation-test correction, ref b8940d4c)

Task A - strengthened A/B (T14). The test now requires exactly one
finite non-sentinel record per (depth, frequency) on both paths for
the full printed grid (48 cells core ray, periapsis bookkeeping rows
must be inactive on BOTH sides or the run fails). Compared directly:
CHI_HAT, Q, DCHR, DTAU (interior), WINX, DLNM, A0-A4, SR (new dump).
Result: wx/dlnm bitwise, chih/q/sr <=1e-7, dchr <=2.4e-7, a0 <=2.1e-6
(EXP(-T1) span sensitivity), all interior dtau <=3.6e-8 -- on the
ordered (ab), alternating (abx) AND non-central (ab2) rays. ab2 runs
NPA=2 with KN_DBG_LS=2 (PA(2)=1.2E4, 7 model depths on-ray): 42
active records, same agreement scale.

Donor identity (A2): KN_DBG_KC1/KC2 are sparse-store slots mapped to
depths via COMP_DEPTH_IDX; KR_DBG_JPR1/JPR2 are ray rows mapped via
J_PNT (driver prints both as JMAP). The test requires identical
donor DEPTHS wherever both paths hold a bracket end, identical
NM/NFM offsets and availability -- a same-availability/
different-depth selection now fails. Verified identical on all
three rays.

SPM/SPP (field-dependent) are not cross-compared: sweep and refresh
use different iterate states. Instead the construction rule is
asserted: where no donor exists the modified source is SR bitwise.
Red-edge closure (A4) is pinned directly: at the compressing
last-frequency cell SPM==SR bitwise on both paths while Q stays
finite -- the phantom slot contributes exactly zero, no self-donor.

Task B - resolved GAM crossing (T15). New cases gz[pn][1-4] set
GAM(4)=+/-1e-2..1e-5*GAMP (the 1e-2 step shifts chi_hat ~1.4e-7 vs
CHI~4e-3 -- measurable at ES16.8). Coefficient-level assertions:
gzz cell 4 has Q==0 and SPM==SPP==SR bitwise at all frequencies;
|coef(eps)-coef(0)| contracts ~10x/decade toward the zero operator
for both signs (observed 1.39e-7 -> 1.3e-10 chi_hat); the chi_hat
increment is sign-blind (|GAM|/DLN); Q flips sign symmetrically
where the compressing donor exists; gzb keeps cell 4 quiet between
opposite-sign neighbours; gcr shows Q<0 left and Q>0 right of the
glide. The 2% field allowance is gone - micro-pair store identity
is a secondary check only.

### Still open

- Task 4: discrete conservation checker (MARCH recursion extracted:
  FLD(I+1)=FLD(I)*A0 + SP(I)*A1 + SP(I+1)*A2 + DS(I)*A3 +
  DS(I+1)*A4 inward; outward swaps A1/A2 and negates slope terms).
  Not started: needs a DS_L dump plus the controlled absorption/
  emission cases; deferred to a dedicated session.
- Task 3 residual: dedicated blue-edge compressing-contribution
  case at L=1 beyond coefficient agreement.
- Task 6: controlled monotonic-vs-V13 then non-monotonic spectrum.


-------------------------------------------------------------------------------

## Commit log

- ced96c12 fg_j_cmf_nonmono_v1.f: KN_DBG_*/KR_DBG_* allocation-gated
  operator dumps; FG_NOSLOT_RAY dealloc in the rebuild teardown.
- ee06b167 drv_knop_ray.f90: ab/abx A/B cases + gpos..gcr profiles.
- b8940d4c test_knop_ray.py: T14 operator equivalence (both sweep
  modes) + T15 GAM sign-flip continuity pins. 16/16 pass.
- <follow-up> fg_j_cmf_nonmono_v1.f: KN_DBG_SR/KR_DBG_SR local-source
  dumps (field-independent leg of SPM/SPP).
- <follow-up> drv_knop_ray.f90: SR columns, JMAP/GAMV donor-map and
  GAM prints, ab2 non-central-ray A/B, KR_DBG dumps + gz[pn][1-4]
  resolved perturbations and gzb zero-cell profiles on the shared
  path.
- 316c8df6 test_knop_ray.py: T14 mandatory coverage + donor-depth
  identity + red-edge SP==SR; T15 coefficient-level zero-crossing
  sequence. 17/17 pass.

## Second follow-up (regression-test gap closure, ref 316c8df6)

T1: DCHR is now in the T14 comparison set (was dumped but never
asserted). Max rel(abs) on ab/abx: chih 5e-8(6e-11), q 4.8e-8, dchr
2.4e-7(7.2e-15), dtau-interior 3.6e-8(1.0), wx/dlnm bitwise 0,
sr 9.4e-8, a0 2.1e-6(4.4e-34), a1 7.2e-8, a2 1e-9, a3 3.6e-8,
a4 6e-9. Denominator is the local reference magnitude with only a
1e-9*scale absolute floor.

T2: _db rejects duplicates, wrong widths, nonfinite values;
missing/unexpected keys fail through the coverage gate with a
diagnostic, not a KeyError.

T3: at every frequency where the compressing branch has a real
in-grid donor, Q(gzp)*Q(gzn) < 0 is asserted explicitly and the
magnitudes must match within 5% (4 resolved flips). Frequencies
where a donor is structurally absent are classified, not forced:
l=1 expansion (no l-1) and l=6 compression (l+NFM>NCF) -- the
latter must additionally show SPM==SR (beyond-grid donor consumed
nothing).

T4: continuity uses three dimensionless deviations -- |chih-chih0|/
chih0, |q|/q_scale with q_scale=|GAM1|/(dlnm*chih0) derived from the
dumped span and zero-case opacity, and |spm-sr|/sr within the same
run (SR drifts ~1e-8 run-to-run). Each sign must contract
monotonically and resolved pairs must sit in [0.03,0.5] -- measured
0.10 at every level for both signs (linear in |GAM|).

T5: donor identity is verified by reconstruction: (SPM-SR)/|Q|
reproduces FLD(depth4, l-1) for expansion and the JPR-bracket
interpolation (1-w)*FLD(d4,l+1)+w*FLD(d5,l+1) for compression,
within 5% (lagged donor = previous iterate of the same converged
field). 5 expansion + 4 compression donors verified; the wrong
stencil would fail the match.

T6: GAM=0 donor-access audit by inspection: both SRC_PRIME and
KN_ORD_SPV take the expansion branch at GAM>=0 where the only read
is the in-grid working-field FLD(R,L-1), guarded by L>1, weighted
by exactly-zero Q. No sparse-store slot or out-of-grid index is
reachable at GAM=0; contribution is bitwise zero (asserted). No
counter needed.
- Phase-2 red-edge regression reverted in tree before commit; net
  diff vs 4bb68a85 is additive only.
- 2af1c196 fg_j_cmf_nonmono_v1.f: KN_COUPLING_OK + KN_RF_COUPLING
  donor gates; consistent drop for missing donors; finite guards +
  fail-closed coefficient build; DLN floor parity.
- 0d38a516 drv_knop_ray.f90: gap/gapr/nanchi/infgam cases.
- df1f2c4e test_knop_lagged_field.py: P1-P3 static pins.
- 4bb68a85 test_knop_ray.py: T11-T13 numerical pins.
- KNOP_AUDIT_TODO.md and the vault note are gitignored local
  artefacts (/*.md and cmfgen_obsidian/ in .gitignore), kept on disk.
- Full build.sh clean rebuild after the edits: Build succeeded,
  31/31 unit test files passed.
- Baseline reference: 53ee4e14.

## CMF_FLUX / observer-frame spectrum audit (ref 5a602ce1)

Confirmed data path:
- Per-iteration OBSFLUX in cmfgen_sub.f (SECTION='CONTINUUM' ->
  CONT_VEL=T -> COMP_OBS_V2 with IPLUS from the current-iterate
  formal solve) is monitoring only.  The file header itself states
  the final spectrum is computed with a different, more accurate
  routine.
- The final spectrum is cmfgen/obs/cmf_flux_sub_v5.f: reruns the continuum
  formal solve on the converged ETA/CHI (its own COMP_JCONT_V4.INC),
  fills IPLUS for every ray at every frequency, then COMP_OBS_V3
  transforms to the observer frame and writes OBSFREQ/OBSFRAME.
- IPLUS = 2*CV_BOUND is the outgoing emergent CMF intensity at RMAX,
  ray-ordered consistently with MU_AT_RMAX/HQW_AT_RMAX.  All solver
  branches fill it: V13 (prod)/V12 (obs) for monotonic, KNOP ->
  FG_J_CMF_NONMONO_V1 in both.  IPLUS is not a moment reconstruction.
- CONT_VEL is structural (SECTION=='CONTINUUM'), not a user flag.
  NORDFLUX output is a CMF flux used only for non-continuum sections.

Defects found and fixed:
1. Trial-spectrum contamination: JFNK_RAD_SWEEP_V1 (contained in
   CMFGEN_SUB) called COMP_OBS_V2 with trial-population IPLUS,
   writing into the same SAVE'd IPLUS_STORE/NU_STORE state and
   advancing the shared OBS_INDX cursor.  Both call sites are now
   guarded with .NOT. JFNK_IN_EVAL; production (always F) unchanged.
2. Silent extrapolation/NaN in COMP_OBS_V2 and V3: the coverage gate
   uses MIN_FREQ_CONV_FAC only, which bounds the LOW end.  For a ray
   with FREQ_CONV_FAC>1 (mu~0) CMF_FREQ can exceed NU_STORE(1); the
   binary search then returned ML_ST=ML_END=1 and LIN_INT divided by
   zero -> silent NaN flux.  Now a per-ray range check writes a
   diagnostic and STOPs.  A non-finite NEW_IPLUS (NaN/Inf passed the
   .LE.0 clamp untouched) is likewise rejected.
3. SOB render solver: obs dispatched DO_NONMONO_SOB models to
   FG_J_CMF_V12 (expansion-biased, wrong/aborts for sigma<-1) while
   KNOP got FG_J_CMF_NONMONO_V1.  The INC branch and the
   SET_FG_GS_OBS_GRID registration now cover DO_NONMONO_SOB too;
   SIGMA_TRUE was already built under both flags.

Config: no new keyword needed.  REL_OBS (DO_FULL_REL_OBS) in VADAT
and DO_CMF_RELO (DO_CMF_REL_OBS) in obs controls are already enforced
for SN_MODEL by check_param_consistency/check_cmf_flux_param_consis.
CONT_VEL needs no input.

Tests: unit_tests/python/test_comp_obs.py + f90/drv_comp_obs.f90,
6/6 pass -- zero-V limit vs CMF (exact), relativistic FCF/ICF vs
independent evaluation (3.3e-16 rel), frequency-only transform
(1.8e-16), RMAX^2 Jy scaling (4.0 exact), non-finite IPLUS rejected,
insufficient top coverage rejected.

Open: Test B (legacy monotonic spectrum comparison) and Test C
(controlled non-monotonic spectrum) still need model runs; solver
independence at the render level now follows from the shared
FG_J_CMF_NONMONO_V1 + COMP_OBS_V3 path but has not been run
end-to-end on a converged model.  KNOP conservation checker and
monotonic V13 comparison remain the next milestones before claiming
non-monotonic observed-spectrum validation.

## Exact-zero GAM donor access + parser hardening (ref a176d338)

Production fixes in fg_j_cmf_nonmono_v1.f -- at GAM=0 the expanding
branches still evaluated the donor expression, and 0*NaN is NaN in
IEEE arithmetic, so a stale or invalid donor could contaminate the
source despite a mathematically zero coupling. Four sites now take an
explicit GAM==0 branch that returns the local source and reads no
frequency donor:
- minus-ray sweep: SOURCE_PRIME = SOURCE_RAY + Q*I_M_PREV (was .GE.0)
- plus-ray sweep: SOURCE_PRIME = SOURCE_RAY + Q*I_P_PREV (was .GE.0)
- SRC_PRIME (refresh): SP_L = SR_L + Q_L*FLD(I,L-1) (was .GE.0)
- KN_ORD_SPV (ordered refresh): SPV = PK_SR + PK_Q*FLD(R,L-1)
Nonzero-GAM expanding and compressing equations are unchanged; the
debug print sites that only display I_M_PREV/I_P_PREV were left.

Proof of no read: new KN_DBG_NDON(i,l) counter is incremented inside
each .GT.0 donor-read branch (allocation-gated; production pays
nothing). A GAM=0 cell can only trip it if the zero guard is
bypassed, so reverting .GT. to .GE. fails the test immediately.
Coverage: abz (new NCASE=33, ab profile with GAM(4)=0 through DRV_AB)
asserts zero reads on the sweep minus/plus and refresh paths, zero
reads on compressing cells (counter is branch-selective), and reads
present on every expanding cell at l>1 (counter liveness). gzz cells
4-5 assert zero reads on the refresh path.

Parser hardening in test_knop_ray.py: JMAP and GAMV now get the same
validation as the DB blocks -- duplicates, malformed rows (including
records shorter than three tokens carrying a known block name), wrong
widths, non-finite GAMV, and out-of-range depth mappings all fail
with diagnostics. JMAP/GAMV must form contiguous 1..N row sets.
DBSD/DBRD rows added to the grid-coverage gate.

Suite: 18/18 pass (T14 now includes abz). Full build.sh clean,
32/32 test files.

===============================================================================
## FIELD-DRIFT METRIC INVESTIGATION (09-Nov-2026) - REPORTED + DIAGNOSTIC FIXES

Question: are the large initial FG (~1.0) and MOM (~1E23) field drifts
a harmless consequence of restoring an older checkpoint?

Verdict: NO on the checkpoint-artifact claim, MIXED on the metric
semantics. Three independent findings:

1. knop_mom's own continuous run (no restart) produced the IDENTICAL
   iteration-13 diagnostics (fg=1.000E+00, mom=1.426E+23, same bad/quar
   counts) the restarted probes showed. The tag-12 checkpoint itself
   recorded the field as non-stationary at dump time (ray CHG=1.0
   pinned, mom CHG=3.278E+02). The drift is a chronic property of this
   model's lagged-field map, not a restart artifact.
2. The one-sided sweep normalization IS asymmetric (T3 built from the
   old column peak): a cell the store held at ~1E-47 taking a fresh
   ~1E-19 prints ~1E4. Verified on real consecutive-generation data
   (knop_mom tag-13 -> probe-cp tag-14, native ifort reader on the
   continuation-record checkpoint layout): max one-sided 7.4E3 (I_P) /
   5.0E3 (I_M) driven by floor cells.
3. But the churn is real at significant cells: of cells carrying >1%
   of column peak, 32.2M/77.2M I_P (42%) and 9.3M/27.9M I_M (33%)
   changed >50% under the symmetric metric between consecutive
   generations; column-peak ratio median 1.001 - amplitude stable,
   structure restructured. negJ/quar persist every iteration at ~1-2%
   of active cells, 100% of dumped failing cells at sigma_true < -1
   (compressing regime, frequencies ~156-200).

Diagnostic-only corrections applied (no gate/threshold/equation
change):

- cmfgen/plane/mom_j_cmf_nonmono_v1.f: new module metrics MOM_GS_MAXSYM
  (symmetric drift, |new-old|/max(|new|,|old|,1e-4*max(old pk,
  fresh pk,1e-28))) and MOM_GS_ARGV(3) (stored/fresh/col-peak at the
  MOM_GS_MAXCHG argmax cell), folded on the same counting gates as
  MOM_GS_MAXCHG in the J and H loops; refresh driver snapshots
  MOM_SW_MAXSYM and saves/restores ARGV alongside ARGI/ARGL.
  GET_MOM_GS_MATURITY returns SW_MAXSYM/ARGI/ARGL/ARGV (4 new trailing
  args; sole caller SOLVE_FOR_POPS updated).
- cmfgen/new_main/mod_subs/solve_for_pops.f:
  * ITER_SUM 'fg' row renamed 'fg_rf' and gained ' sw=': the printed
    number was the refresh residual KN_RF_RESID, mislabeled
    'lagged field drift' (GS_NBAD=-1 sentinel falls back to 'fg /
    sweep drift' when no refresh ran).
  * 'mom' row gains ' sym:' = MOM_SW_MAXSYM.
  * KNOP field block gains a 'sym' row: symmetric drift + argmax cell
    (dep,frq) + (stored,fresh,col-peak) triple.
  * The two 'largest drift' prints relabeled 'residual/sweep max'
    (GS_MAXCHG/MOM_MAXCHG hold the refresh residual after a refresh,
    the sweep drift otherwise).
  * GS_DRIFT_DIAG mom line gains /sym and the argmax triple.
  * ITER_HISTORY columns: fg_d -> fg_rfd (renamed, same value), new
    fg_swd + fg_rfdel, new mom_symd (14 -> 17 trailing ES columns).
- Vault updated: reference/Diagnostics and OUTGEN.md (ITER_SUM rows,
  ITER_HISTORY column names), control/Field maturity and GS freeze.md
  (MOM_SW_MAXSYM/ARGV observables, GS_DRIFT_DIAG sym channel).

Consequences still open (not part of this change):
- The lagged field never converged on m15_2p82d: refresh residual
  ~1.0 for 7+ iterations, mom sweep drift 1E11-1E31 every iteration.
  The maturity gates are the only thing standing between JFNK-FRZ and
  this field; JFNK-FRZ admitted at it=15 while fg still read 5-7%.
  Gate-threshold review belongs to the JFNK task chain (JT7/JT8),
  not to a diagnostic change.
- MOM_J_NONMONO_ERRORS is capped at 80k lines (20k cells); probes'
  copies are stale Nov-4 files (the running processes' own copies
  were never written - the REWIND on INIT did not fire). Any future
  per-cell forensics should raise or roll the cap.
