# Outstanding issues -- cmfgen_nonmono_dev_V7

Audit of `MOM_J_CMF_NONMONO_V1` for the TypeIIn_KNOP_TRIDIAG run
(22-Sep-2026).  Trigger: outer luminosity blow-up,
`RLUMST(1) = 2.85E8 -> 2.97E11 -> 1.36E15` over iterations 1-3.

Each claim is tagged VERIFIED (source + runtime evidence),
SUPPORTED (source evidence, mechanism inferred), or OPEN.

## Status update (30-Sep-2026, evening)

The marginal-release trial from this morning's section is
runtime-validated and its one structural gap is closed.

VERIFIED (jobs 12040589 + 12040591, Dessart m15 dif/dbg configs,
USE_LAM_ES=F and T): the marginal-release machinery pulsed correctly
-- it30/it32 trials released, resacc accepted both (nrej=0),
lum_ratio moved 0.16 -> 1.03 -- but N_FROZEN counted every frozen
iteration including accepted-release ones, and the raw-mature gate
(MOM_SW_MAXCHG <= 10; observed >= ~1.3e3, asymptotes ~1e3 on the
one-bin-upwind marginal cells) can never produce the 2 consecutive
mature sweeps needed to clear it.  Both runs stopped at
GS_FREEZE_LIM=60 while progressing: VERIFIED that un-refunded
counting makes the marginal tier self-defeating.

FIXED (solve_for_pops.f, rebuilt exe 15:56): accepted trials refund
N_FROZEN to MOM_TRIAL_MIN and increment N_TRIAL_ACC (persisted in
GS_STATE, survives requeue); rejected trials keep counting; new
VADAT key GS_TRIAL_LIM (default 15, <=0 off) bounds consecutive
accepted trials between mature sweeps.  MOM_TRIAL terms prints the
9th trlim term + ntr=.

Paper-fidelity fixes in the same build (mom_j_cmf_nonmono_v1.f):
dead midpoint CHI_HAT loop removed (read uninitialized/stale
GAM_MAG_NM); expanding-row PSIPREV now |GAM*GAM_FAC_PREV| so a sign
flip of 1+sigma*K_on_J cannot inject a negative upwind source
(Knop I Eq. 9 positivity).

3-iteration itest on interactive (job 12042808,
dessart_m15_2p82d_KNOP_MOM_dif_itest.yaml): COMPLETED 0:0, bordered
tables correct incl. `frozen= T` (a format-descriptor bug, 4(A,A)
-> 3(A,A), was caught by this run and fixed before resubmission).
Full model resubmitted: FINAL 12043642, OBS 12043643 on small.

STILL OPEN from below: the single-cell MOM_SW_MAXCHG asymptote
(~1e3) is unchanged by the refund -- it now only gates the strict
maturity tier, but a genuinely contracting field still cannot reach
raw-mature on this model, so N_FROZEN resets to 0 remain
unreachable.  If the new run converges only via trial pulses, decide
between a per-cell vs bulk MAXCHG metric or formally promoting the
marginal tier to maturity.

## Status update (23-Sep-2026, morning)

Second validation (`TypeIIn_KNOP_TRIDIAG_EPS2`, job 11942988, sigma*n
mask): REGRESSION.  it1 L(d=1)=2.95E40 on the FIRST iteration (vs
~3E8 healthy in both prior runs); negJ ~7E5/pass; capped writes
climbed to 3123 (non-finite XM in double precision, not the old
float32 channel).  argmax hopped over scattered depths/frequencies
-- the instability went diffuse, not localized.

The three-variant trend is monotonic and decisive:

| variant | it1 L(d=1) | it2 |
|---|---|---|
| all-backward EPS (original) | 2.85E8 | 3.0E11 |
| 1+sigma*K mask (few flips) | 3.04E8 | 7.5E10 |
| sigma*n mask (many flips) | 2.95E40 | crashed |

More forward-redirected midpoints = strictly worse.  Interpretation:
the backward J_prev read is a within-sweep Gauss-Seidel update
(fresh, contractive); the J_next read is a stale Jacobi read (adds
lag).  Each flipped midpoint trades fresh coupling for stale
coupling, so redirecting EPS expands the map.  The EPS redirect
concept is counterproductive; the divergence driver lies elsewhere
(|EPS|>HU non-M-matrix rows and/or the intrinsic map).

Also found: a format/variable-type crash (forrtl 61) in the new
diagnostic WRITE at dump_mom_jnext_gs -- a missing A descriptor put
' nonMM rows=' on an I10 slot.  This is what killed job 11942988.

FIX (09:57 build): EPS_SGN flip disabled (array stays +1, machinery
retained for a possible conservative-form retry), format bug fixed.
EPS restored to unconditional backward difference = original code.

New validation: TypeIIn_KNOP_TRIDIAG_DBG2 on the interactive
partition (submission: debug), FINAL 11947270, OBS 11947271.
Purpose: fast-iteration monitor of the baseline + nonMM/argmax
diagnostics to localize the true driver.  DBG (11947244) was
cancelled after 5 min -- it loaded the defective binary before the
relink.

## Status update (22-Sep-2026, late evening)

First validation run (`TypeIIn_KNOP_TRIDIAG_EPS`, job 11942275)
result: STILL DIVERGING, cancelled at iteration 4.

- L(d=1): 3.04E8 -> 7.55E10 -> 6.34E16 (vs old run 2.85E8 ->
  2.97E11 -> 1.36E15).  ~4x lower at it2 but the same explosive
  trend -- not a cold-start transient.
- `capped=0` throughout: the HRF float32->Inf channel is confirmed
  closed (VERIFIED).  No XMrhs=Inf/NaN rows.
- negJ counts fell ~10x (it2: ~2E4 vs old ~5E5) but refresh still
  aborts as diverging: residuals contract ~10 passes then explode
  (it2: 2.3E12 -> best 1.6E4 -> 8.9E12).

Source audit of the new run's error file produced a REAL defect in
the EPS_SGN mask:

- The dominant negJ site moved to depth ~474 (sigma~-2.9, mild
  compression) where `NMID_ON_J ~ -0.06` is NEGATIVE.  The strong
  block (depths 170-189, sigma~-38) also has `NMID_ON_J < 0`.
- NMID_ON_J is the n-moment (mu^3-type) Eddington factor and is
  legitimately SIGNED: n<0 under inward streaming.  The EPS term is
  the n-weighted part of the comoving frequency derivative, so its
  characteristic direction is set by sign(sigma*n), NOT by the
  J-transport criterion 1+sigma*K that the first implementation used.
- Consequence: at sigma*n>0 cells (e.g. the n<0 strong block) the
  Kmid mask wrongly applied the forward difference, sign-inverting
  the lagged EPS source exactly where the field was already unstable.

FIX (22:45 build): EPS_SGN mask now keys on
`AV_SIGMA_TRUE_NM(I)*NMID_ON_J(I) < 0` (i.e. EPS<0), which is the
correct upwind direction for the n-weighted frequency transport.
Note the matrix form collapses: with the sign-correspondent mask
s*EPS=-|EPS| always, so TA=-HL-|EPS|, TC=-HU+|EPS| regardless of
direction -- the mask only selects which lagged field the RHS reads
(J_prev vs J_next) and the ln-nu spacing.  The TC>0 (non-M-matrix)
exposure whenever |EPS|>HU is unchanged and remains Issue 3.

Also added: argmax diagnostics in MOM_KNOP_REFRESH_V1 -- each pass
line now reports the depth/frequency index of the max relative
change (`argI/argL`) to localize the unstable mode.

New validation: `TypeIIn_KNOP_TRIDIAG_EPS2` (fresh cold start),
FINAL job 11942988, OBS 11942989, running the 22:45 binary.

## Status update (22-Sep-2026, evening)

Implemented in `plane/mom_j_cmf_nonmono_v1.f` and rebuilt into
`exe/cmfgen_dev.exe` (21:08 build):

- Issue 1 (HRF precision): FIXED.  `HRF` is now `REAL(KIND=LDP)`;
  the float32 -> Inf channel is closed.
- Issue 2 (upwind coverage): PARTIALLY FIXED.  The `EPS*Jsum`
  (sigma*NMID_ON_J*J) coupling is now upwinded per midpoint:
  a new `EPS_SGN_NM` mask (`1+AV_SIGMA_TRUE*Kmid < 0`, the
  midpoint analogue of the `Q_NM` criterion) switches each
  midpoint's contribution to the forward difference
  `EPS*Jsum_next - EPS*Jsum_k`, with the lagged `J_next` read
  through `JUP_NM` (same store/semantics as `PSI_NM*J_next`,
  including the refresh `JNEXT_RF_READ` snapshot) and the
  forward `ln(nu_k/nu_{k+1})` spacing in the EPS coefficient.
  Applied at: interior tridiagonal TA/TB/TC, interior RHS,
  the `OUT_BC_TYPE>1` boundary row, the negative-J rebuild
  path, and the H reconstruction.  The `W*dH/dlnu` (`HS*H_prev`)
  term deliberately stays bluer-directed: under `N_ON_J`,
  `NMID_ON_HMID=0` makes it a pure V/R aberration term whose
  transport sign never reverses; all sigma-carrying transport
  is in the EPS term.  Under `G_ONLY`/`MIXED` that argument
  does not hold and an `H_next` store would be needed (OPEN).
- Validation: config `configs/TypeIIn_10d_DEV7_KNOP_TRIDIAG_EPS.yaml`
  (fresh `TypeIIn_KNOP_TRIDIAG_EPS` dir -- cold start required, a
  warm restart would restore the diverged `J_NEXT_FRQ` from the
  old run's `MOMJNXT_GS`).  Production submitted via
  run_model_grid.py: FINAL job 11942275, OBS 11942276, running
  `cmfgen_nonmono_dev_V7/exe/cmfgen_dev.exe` (21:08 build).
  Interactive test 11942188 (prep-only validation) pending on
  AssocMaxJobsLimit.  Watch OUTGEN for: XMrhs=Inf rows,
  MOM_KNOP_REFRESH residual growth, capped/neg-J counters, and
  RLUMST(1) vs the old run's 2.9E8 -> 3E11 -> 1.4E15 -> 2E40.
- Issues 3 (unclamped EPS in the matrix) and 4 (intrinsic
  non-contractive map) remain OPEN pending the validation run.

## Setup status (all VERIFIED, no action needed)

- `USE_LAM_ES=F` + `DO_NONMONO_KNOP=T` dispatches to
  `MOM_J_CMF_NONMONO_V1` per frequency (OUTGEN "Calling
  MOM_J_CMF_NONMONO_V1").  `SIGMA_TRUE`, `J_NEXT_FRQ` GS store,
  `SET_MOM_NU_NEXT` forward spacing, `MOM_KNOP_REFRESH` driver,
  `MOMJNXT_GS` checkpointing and the `N_LAM_RWARM`/immaturity gates
  are all engaged as designed.
- `SOL_METH=TRI` / `NUM_BNDS=3` do NOT reach the moment solver.
  The moment equation is assembled tridiagonal and solved with
  `THOMAS` regardless; TRI only selects the BA block solver in
  `SOLVEBA_V13`.  With `USE_LAM_ES=F` the BA off-diagonal depth
  blocks are zero (no `TX/TVX` coupling from `do_var_cont.f`), so
  TRI reduces to per-depth DIAG solves.  During forced-LAMBDA
  iterations `TEMP_CHAR='DIAG'` is used regardless.  TRI is a
  no-op for this failure.
- `RLUMST(1)` blow-up is a real divergent radiation field:
  `SOB(1)` is built from the moment `RSQHNU` anchored at
  `HFLUX_AT_OB = HBC_CMF(1)*RJ(1)` (comp_j_blank.f:845), so
  `L(d=1)` tracks the outer-boundary `RJ(1)` of the moment solve.
- Sibling comparison: `TypeIIn_KNOP` (same model, `USE_LAM_ES=T`,
  FG path) converges smoothly `L_outer 2.0E8 -> 1.1E9` over it1-5.
  This run is the first MOM-path attempt on the TypeIIn model.

## ~~Issue 1: HRF single-precision -> Inf injection~~ (FIXED)

`MOM_KNOP_REFRESH_V1` declares the lagged-H scratch as default
`REAL` (float32, max ~3.4E38):

- `plane/mom_j_cmf_nonmono_v1.f:1958` `REAL, ALLOCATABLE :: HRF(:,:)`
- `:2131` `HRF(1:ND_NM,L)=RSQHNU_NM(1:ND_NM)` -- LDP -> float32 store
- `:2113` `RSQHNU_NM(1:ND_NM)=HRF(1:ND_NM,L-1)` -- Inf read back

`RSQHNU_NM` is bounded only by the `|H|<=J` clip against `XM_NM`,
which the 1E40 cap (`:1420-1428`) allows to reach ~1E40.  Any
value above float32 max overflows to `Inf` on the store; the next
re-solve then builds `XMrhs=+/-Infinity` via `VB*Inf`/`VC*Inf`
while `TA/TB/TC` stay finite, and the whole frequency (all 724
depths) is zeroed.  This is exactly the OUTGEN signature at
refresh pass 9 of iteration 3 (`capped=1448 = 2*724`,
`XMrhs=Infinity` / `-Infinity` on consecutive frequencies).

The build has no `-r8`/`-autodouble` (viper_intel_nodebug flags:
`-extend-source -O3 $(AVX512) -align array64byte`), so default
`REAL` is real*4.  The "float32 is adequate" comment predates the
1E40 XM cap, which made the runaway range reachable.

**Fix (21:08 build, VERIFIED in run 11942275):** `HRF` promoted to
`REAL(KIND=LDP)`.  `capped=0` throughout the EPS validation run and
no `XMrhs=Inf/NaN` rows -- the Inf channel is confirmed closed.
This stops false whole-frequency zeroing and makes residuals
meaningful; it did not cure the underlying divergence (issues 3-4).

## ~~Issue 2: sign-aware upwinding only on the dJ/dlnu term~~
## (RESOLVED -- the redirect is harmful, reverted)

At compressing depths (`Q_NM<0`, physical frequency characteristics
run red->blue, upwind = k+1) only the J-equation `dJ/dlnu` coupling
reads the redder field (`PSI_NM*J_next`, explicit weight = implicit
diagonal weight -- correct, matches FG `S + |Q_k|*I_NEXT`).

Every other frequency coupling still reads the bluer field (k-1)
unconditionally:

- `HS_NM*RSQHNU_PREV` in the H reconstruction (`:1531-1539`) --
  the `W*dH/dlnu` term.
- `VB*H_prev + VC*H_prev` in the J-eq RHS (`:1226-1241`) -- the
  same term after elimination into the tridiagonal.
- `EPS_PREV*J_prev` terms (`:1249-1250`, `:1305-1306`, rebuild
  path `:1341-1343`) -- the sigma*nmid*J frequency coupling.
  Under `N_TYPE=N_ON_J`, `NMID_ON_HMID` is identically zero
  (`get_moms_non_rel.f:169`), so `W=GAMH*|1+0|=GAMH` carries no
  sigma dependence at all: the EPS terms are the ENTIRE
  sigma-carrying part of the H equation, and they are
  backward-differenced toward blue at every depth.

**Resolution (23-Sep-2026): TRIED AND REVERTED.**  The EPS*Jsum
coupling was implemented with a per-midpoint upwind redirect to the
lagged `J_next` field (`EPS_SGN_NM` machinery, forward difference
`EPS*Jsum_next - EPS*Jsum_k`, forward `ln(nu_k/nu_{k+1})` spacing,
applied at all 6 sites).  Two mask criteria were tested:

- `1+AV_SIGMA_TRUE*Kmid < 0` (job 11942275): it2 L=7.5E10 -- still
  diverging.
- `AV_SIGMA_TRUE*NMID_ON_J < 0` (job 11942988): it1 L=2.95E40 --
  catastrophic regression.

The trend is monotonic: more flipped midpoints = strictly worse.
The backward `J_prev` read is a within-sweep Gauss-Seidel update
(fresh field, contractive); `J_next` is a stale Jacobi read (last
sweep's field).  Redirecting EPS removes contractive coupling from
the map -- the original all-backward form is the more stable
choice.  EPS is restored to unconditional backward difference
(EPS_SGN machinery retained but disabled, 09:57 build).

The `HS*H` term was never flipped (correctly: under N_ON_J it is a
pure aberration term whose transport never reverses).

## Issue 3: unclamped EPS inside the tridiagonal matrix (SUPPORTED -- OPEN)

`EPS = GAMH*AV_SIGMA_TRUE*NMID_ON_J/(1+W)` (`:1077-1082`) uses
unclamped sigma_true (measured -5 to -63 at bad rows in
`MOM_J_NONMONO_ERRORS`) and enters the matrix:

    TA = -HL - EPS,   TC = -HU + EPS,   TB = ... - EPS(I-1) + EPS(I)

(`:1150-1153`).  The header claim "the tridiagonal system is always
well-conditioned" (`:294-299`) only covers `chi_hat`/`PSI`/`W`; EPS
is excluded from that argument.  With |EPS| exceeding HL/HU in
optically thick compressions the off-diagonals can flip sign
(`TC>0` seen in the first-bad-row dumps) and the M-matrix structure
is lost -- a plausible direct source of the ~1E6 significant
negative-J resets per sweep, which then seed the lagged map.

Additional data point: `NMID_ON_J` itself goes negative at the
strongest compressions (ray-solve N-moment sign flip, e.g. row 330
in the error file), so the EPS weight changes sign there as well.

**Next candidate fix:** with issue 2 resolved, this is now the
prime suspect.  A positivity-preserving option: where |EPS|>HU,
move that midpoint's implicit EPS contribution into the RHS on the
same-frequency lagged field (`J_NEXT_FRQ(:,k)` in the main sweep /
`JNEXT_RF_READ(:,L)` in refresh), keeping the tridiagonal a clean
M-matrix while preserving the fixed point.  The new `nonMM rows=`
and `max|EPS|/HU=` diagnostics (running in job 11947270) measure
how widespread the violation is.

## Issue 4: intrinsic non-contractive lagged map (documented, OPEN)

Already recorded in AGENTS.md (21-Sep-2026): the in-sweep solve
produces significant negative-J rows at sigma_true ~ -37..-55 even
with the `PSI_NM` weight fix, the 1E40 cap and best-pass restore.
On this model the refresh diverges every iteration (best residuals
1E2 -> 1.7E7 -> 7.7E9 and rising) and the restored field feeds the
next in-sweep solve, which amplifies it into `RJ(1)` -- hence the
outer-luminosity runaway.  The TypeIIn compression strengths
(sigma_true to -63, 7 velocity sign changes) are beyond the
containment demonstrated on MOM_TFIT (job 11889372, pinned at
~1.5E8).

Whether the residual instability survives a complete sign-aware
upwind treatment (issue 2) or needs an additional fix (e.g.
EPS-side limiting, issue 3) is OPEN until those are tried.

## Recommended order of work (updated 23-Sep-2026)

1. ~~`HRF` -> LDP (issue 1)~~ DONE, verified (capped=0).
2. ~~EPS forward-difference at compressing cells (issue 2)~~ TRIED,
   REVERTED -- the redirect expands the map; original backward form
   retained.
3. Baseline + diagnostics run in progress (DBG2, job 11947270):
   measure `nonMM rows=`, `max|EPS|/HU=`, `argmax(d,f)` to confirm
   issue 3 as the driver and localize it.
4. If confirmed: positivity-preserving EPS treatment (explicit EPS
   on the same-frequency lagged field at |EPS|>HU cells) -- keeps
   the matrix M-matrix, preserves the fixed point.
5. Until the MOM path holds `RLUMST(1)` stable over multiple
   iterations, `USE_LAM_ES=T` (FG path) is the only converged
   option for this model; MOM path results should not be trusted.

Validation note: a rebuild alone is not validation -- any change
must run through `scripts/run_model_grid.py` and show a stable
`RLUMST(1)` over multiple iterations, not just a clean first pass.

## Status update (24-Sep-2026, afternoon) -- F1/F3/F4 landed, new boundary finding

Commits (each rebuilt clean, separate commits):

1. `36b03eb` (F1): outer-BC upwind now keys on the row's own signed
   coefficient `A_BND = GAM(1)*(HBC+NBC*sigma_true(1))` instead of
   interior `Q_NM(1)`.  `PSI_NM(1)=|A_BND|` unchanged; the `.GT.1`
   branch keeps `Q_NM(1)` (that row is the negated interior row).
2. `d8507dc` (F3): `GET_MOMS_NON_REL` size-checks the 10
   MOD_RAY_MOM_STORE arrays (`SIZE(K_ON_J)/=ND` -> dealloc) so an
   ND change cannot write out of bounds.
3. `257ef3b` (F4): `KN_M_OUTBC` scalar -> `KN_M_OUTBC(NCF)`;
   captured at slot K, replayed as `KN_M_OUTBC(L)` in refresh.

Validation run: `M15_HeOCa_25d_DEV7_knop_MOM_bndfix`, FINAL
11960113 / OBS 11960114 (interactive), cancelled after 21
iterations -- diverging.

CAVEAT (audit correction): the 11960113 binary was built from a
DIRTY tree -- five functional changes were still uncommitted at
build time (solve_for_pops FORCE_LAMBDA reset, NG FORCE_LAMBDA
gate, GS_FREEZE smooth guard, do_ng_accel_band 0/0 guard,
check_ions_present init guard) plus format-split cleanups.  The
F1 mechanism is orthogonal to all of them so the finding stands,
but the run does not map to a single commit.  Those changes are
now committed (504e4c9, 5829cf7, 080fbc6, 026195c, 619dc3a,
fc5ed0e); subsequent builds are clean.

F1 verified at its target: ZERO depth-1 records in
MOM_J_NONMONO_ERRORS over 21 great iterations (47 records total,
all interior depths 409-412; baseline had 247 depth-1 records --
earlier "thousands" claim corrected by direct count).  The
boundary row no longer sits on the `|XM|/10` fixer floor.

NEW ISSUE (VERIFIED, HIGH) -- boundary row is non-contractive
when HBC < 0:

The boundary row is
    TB(1) = K_ON_J(1)*QF(1)/DTAU(1) + |A_bnd| + HBC
    XM(1) = |A_bnd| * J_next(1)          (compressing branch)
With HBC < 0 the diagonal shrinks below |A_bnd|, so the lagged
J_next weight exceeds 1 -> expansive GS map localized at depth 1.
OUTGEN confirms: `argmax(d,f)=(1,6471)`, max rel change 2.37.
J(1) then grows ~1.5-1.7x/iter; L(d=1) = HBC_CMF(1)*RJ(1) crossed
zero at it7-8 and ran -1.4E8 -> -1.5E11 by it21 (HBC_CMF(1)<0 =
negative emergent flux HNU_AT_OB from the FG solve; RJ(1) stays
positive, no fixer records).  Interior XM magnitudes inflated to
~-1E9 as the corrupted boundary field spread.

The old misclassification masked this: expanding branch gave
PSIPREV*J_prev < 0 -> XM(1) < 0 -> fixer pinned J(1) to the floor
every sweep, damping the feedback (baseline still exploded at
it70 via a different path).  Old depth-1 records show
PSIPREV(1) < 0, i.e. HBC < 0 was already present in the baseline
-- the negative boundary Eddington is a pre-existing formal-solve
artifact (HNU_AT_OB < 0, or JNU_STORE(1) on its 1E-30 floor making
the ratio garbage), not created by F1.

Interpretation: with A_bnd < 0 the boundary row is degenerate
under EITHER classification -- expanding gives a negative explicit
source (fixer), compressing gives an expansive lagged coupling.
The consistent fix likely needs the boundary Eddington in the row
to be physically sane (HBC < 0 at the outer boundary is net inward
flux -- unphysical for this atmosphere, a formal-solve artifact).
Candidate: floor the HBC used in the boundary-row diagonal /
classification only (not the HBC_CMF used for HFLUX_AT_OB
diagnostics), which restores diagonal >= |A_bnd| and makes the
compressing weight <= 1.  Alternatively diagnose why FG returns
HNU_AT_OB < 0 (N_NEG_INT counters, grazing-ray negatives, or the
JNU_STORE(1) 1E-30 floor inflating the H/J ratio at EUV).  Same
defect class as the F2 interior PSIPREV sign mismatch, at the
boundary.

## Status update (24-Sep-2026, evening) -- boundary HBC floor implemented

Audit refinements to the proposed fix, applied in `368fa76`:

- The floor must cover BOTH branches, not just the diagonal: with
  HBC < 0 the expanding branch still injects
  PSIPREV(1) = GAM*(HBC_PREV + NBC_PREV*sigma_true(1)) < 0
  (negative explicit source -> fixer).  Consistent formulation:
  row-local `HBC_ROW = max(HBC,0)` and
  `HBC_PREV_ROW = max(HBC_PREV,0)` (the .NOT.(.GE.0) test also
  catches NaN), used in all three places -- A_BND
  (classification), the TB(1) diagonal, and PSIPREV(1).
  Physically: when the measured H/J is negative the Eddington part
  of the closure is unusable; dropping it substitutes H(1)=0, the
  physical condition for an unilluminated boundary, and
  A_BND = GAM*NBC*sigma_true carries the pure transport sign.
- `HBC_CMF` is untouched: HFLUX_AT_OB keeps reporting the measured
  (possibly negative) flux, so the diagnostic stays honest.
- The `OUT_BC_TYPE>1` branch is unchanged (different row, no
  evidence it fails).
- Applies automatically in MOM_KN_REFRESH re-calls: they restore
  HBC/HBC_SAVE from the captured KN_MHBC store and pass through
  the same boundary-row code.

Instrumentation: each floored frequency writes a 3-line record
(FREQ_INDX, FREQ, HBC, NBC, HBC_PREV, NBC_PREV, sigma_true(1),
JNU_STORE(1), HNU_AT_OB) to `MOM_BND_HBC_FLOOR` (append mode, lazy
open so no empty file in runs that never floor).  A per-sweep
count appears in the GS sweep diagnostic as `bnd HBC floors=`.

Audit context (VERIFIED): argmax frequencies in the diverged run
localize to ~1230-1270 A (FREQ_INDX ~6400-6500) and ~22 um
(~23024) -- near-zero-flux regions where the thick-BC extension
can emit more inward than the star emits outward
(CV_BOUND = 0.5*(I_P - I_M) < 0 -> HNU_AT_OB < 0).  Zero
negative-intensity and zero JNU-clamp events over the whole run,
so this is the extension back-emission channel, not the I<0 or
degenerate-J channels.  The CV_BOUND/HNU_AT_OB/HBC construction
is identical to stock V13 -- an inherited latent vulnerability,
not a nonmono regression.

Validation runs (two; the first was confounded):

`M15_HeOCa_25d_DEV7_knop_MOM_bndflr` (floor + PRUNED MODEL_SPEC),
FINAL 11961324 / OBS 11961325 -- DIVERGED, cancelled at it5.
L(d=1): 3.3E8 -> 3.0E8 -> 2.0E11 -> 2.0E15 -> 2.2E30 (~4-15
dex/iter), argmax at inner depths 383-448, 15 non-finite XM
writes (bndfix had 0 in 21 its).  CONFOUNDED: cmfgen_base/
MODEL_SPEC was swapped at 13:04 between the bndfix submission
(11:28) and this one -- bndflr got the pruned spec (He2 removed,
O2 at 3,3,5, NCF=25131) while bndfix ran the full 18-species spec
(NCF=23723).  The interior explosion is SPEC-DRIVEN, not
floor-seeded: the pruned spec lacks the opacity/ionisation stages
that anchor the EUV radiation field at inner depths.  Treat the
pruned spec as its own open issue -- the same run with the full
spec contracts (below).

`M15_HeOCa_25d_DEV7_knop_MOM_bndflrF` (floor + the exact full
18-species MODEL_SPEC bndfix ran, via cmfgen_base_fullspec/),
FINAL 11961508 / OBS 11961509 -- CONVERGING (clean A/B):

- GS max rel. change: 1.0 -> 4.9 -> 0.87 -> 0.98 -> 0.73 -> 0.42
  -> 0.28 -> 0.23 -> 0.21 -> 0.19 -> 0.16 -> 0.14 -> 0.13 (it1-14)
  -- below 1.0 from it3 and monotonically contracting.  bndfix at
  the same iterations: 7.4, 15, 8.2, 4.1, 2.9, 2.4 -- always >1,
  boundary mode still active, diverged by it21.
- argmax: (200,2) -> (411,6089) -> (1,64xx) for it3-6 (still the
  largest mover but contracting, i.e. a converging fixed point not
  an expansive map) -> (2,10802) it7-13 -> (3,10802) it14.
  Off depth 1 from it7.
- L(d=1): 3.34E8 -> 3.11E8 -> 2.92E8 -> 2.76E8 -> 2.64E8 -> 2.55E8
  -> 2.48E8 -> ... -> 2.24E8 (it14): smooth positive decay, no
  zero crossing (bndfix crossed to -1.4E8 at it8).  HBC_CMF stays
  unfloored so this is the measured flux -- the boundary J is no
  longer corrupted.
- L(d=ND) = 0 throughout (ZERO_FLUX).
- MOM_J_NONMONO_ERRORS: ZERO depth-1 records; ~47 records at
  depths 409-412 (the marginal inner cells, bounded).
- Zero non-finite XM writes (bndflr/pruned had 15).
- bnd HBC floors: 9607 -> 6746 -> 5359 -> 5042 -> 4985 -> 4935 ->
  4915 -> 4861 -> 4834 -> 4788 -> 4695 -> 4664 -> 4622 -> 4283 ->
  4197 -> 4109 -> 4036 -> 3974 -> 3919 (it1-23), decreasing as the
  field converges.  Per-event log: ~6.3k of 81k events have
  HBC < -0.01 (the meaningful negative Eddingtons), ALL at
  FREQ_INDX 4000-6999 (~55-200 A EUV/FUV) -- the predicted window;
  ~44k events are the JNU_STORE(1)=1E-30 underflow tail
  (HBC ~ -1E-37, floor to 0 is exact).  Min HBC ~ -0.41.
  No optical/NIR floor events -- extension back-emission diagnosis
  SUPPORTED.
- it18-23: maxchg 0.091 -> 0.064 (contracting ~8%/iter), argmax
  (2,11257) -> (1,12484), L(d=1) 2.16E8 -> 2.08E8 smooth positive
  decay, no depth-1 records, zero non-finite writes.

VERDICT (it23): the floor eliminates the boundary expansive mode
and the run is on a sustained contraction -- the first MOM-path
trajectory on this model that is not diverging.  Full convergence
not yet reached (maxchg ~0.064 still decaying ~8%/iter); the job
continues past this writing via graceful-stop requeue.

====================================================================
UPDATE 24-Sep-2026 (pm): bndflrF CONVERGED; interior negJ cause
found and guarded
====================================================================

bndflrF final outcome: FINAL 11961508 exited cleanly at 61 min,
"Convergence assessment CONVERGED" at iteration 74 (single_final
mode).  L(d=1)=1.1509E8 erg/s, L(d=ND)=0, last max % increase
2.06% at depth 15.  OBS 11961509 completed.  The floor run is a
converged reference; it converged on pure LAMBDA stepping because
the persistent interior negJ events held GS_FREEZE the whole run.

Interior negJ mechanism -- CONFIRMED (was inferred from records):

- Diagnostics commit 6723811 added a per-record dump of
  TB,TA,TC,absorp=(1-COH)*DTAUONQ,COH,DTAUONQ to
  MOM_J_NONMONO_ERRORS.
- Dedicated run M15_HeOCa_25d_DEV7_knop_MOM_negdiag (FINAL
  11962966, same physics as bndflrF) captured 49 records in the
  first sweep before cancellation: 33 had COH=ESEC/CHI>1, 20 had
  outright TB<0, 22 had TB<|TA|+|TC| (non-dominant).  Deepest
  cell: XM=-5.7E9 at TB=0.004 vs off-diagonal sum 0.055 --
  near-singular; RHS positive throughout.
- Cause: ESEC and CHI are interpolated onto the expanded grid
  independently, so ESEC can exceed CHI at near pure-scattering
  cells (thick inner ejecta at FUV freqs, sigma~0).  (1-COH)*DTAUONQ
  then subtracts up to ~0.34 from a ~0.34 diagonal, collapsing it;
  the row solve flips to a huge negative J that drags neighbours
  down through -TC*XM coupling.
- These events recur identically every sweep (the crossing is a
  fixed property of the grid), so MOM_NEGJ>0 forever -> permanent
  GS_FREEZE (F5 symptom) -> pure-LAMBDA convergence.

Guard commit fec6fed: absorption fraction floored at the pure
scattering limit, TB uses MAX(0,1-COH_VEC)*DTAUONQ in both
interior branches; the .GT.1 outer-boundary row (negated interior
equation, dead under OUT_BC_TYPE=1) gets the symmetric
MIN(COH,1)-1 term.  COH_VEC left raw so the record still shows
the ESEC>CHI crossing.

Validation run M15_HeOCa_25d_DEV7_knop_MOM_cohgd (FINAL 11963093,
identical config to negdiag/bndflrF, fresh dir): ZERO negJ records
through the first 3 iterations (vs ~49 in sweep 1 unguarded).
Watch for: MOM_NEGJ=0 -> GS_FREEZE release -> first BA steps and
faster convergence than bndflrF's 74 pure-LAMBDA iterations.

TypeIIn remains a DIFFERENT mechanism (do not conflate): its negJ
flood sits at compressing cells (sigma ~ -3 to -98), only 1.8% of
records have ESEC/CHI>1, nonMM rows=0 and max|EPS|/HU~25 -- the
negative bands are the non-diagonally-dominant Z-matrix /
advection-dominated class (F7), not the diagonal-collapse defect.
Its requeue onto the guard+diagnostic binary will print TB/TA/TC
at those rows for direct confirmation.

Open after this change:

- Whether negJ=0 is enough for GS_FREEZE release: other gate terms
  (MOM_MAXCHG<0.05, capped counts, FG maturity) already passed in
  bndflrF, so cohgd should show BA steps once its MOM lagged field
  matures (~it15-50).  If records still appear, the TB diagnostics
  will show the surviving mechanism directly.
- The guard removes a symptom of independent ESEC/CHI interpolation;
  the interpolation itself is untouched (physically correct since
  absorption cannot be negative, but ESEC>CHI cells still have a
  ~1-4% opacity inconsistency worth understanding).
- TypeIIn |EPS|/HU~25 nondominant-row negatives: needs a separate
  treatment (implicit EPS, or dominance-restoring weighting) --
  unchanged, unguarded.
- F5 gate design stands: with this guard M15 reaches negJ=0, but a
  model whose interior negatives are intrinsic (TypeIIn) still
  cannot release the freeze; a bulk/class-aware gate may still be
  wanted.

====================================================================
UPDATE 24-Sep-2026 (late pm) -- strong-compression issue review
(TypeIIn bndflr run, guard+diagnostic binary, jobs 11962018/19)
====================================================================

Row diagnostics now live in TypeIIn records (TB,TA,TC,absorp,COH,
DTAUONQ per negJ event).  1837 records reviewed; findings:

ISSUE 3 (unclamped EPS in the matrix) -- CLOSED / SUPERSEDED.
The depth-upwind EPS form (23-Sep) is verified working at full
compression strength: nonMM rows=0 every sweep, TC<0 in 100% of
records, TB>0 in 100% (vs the |EPS|>HU -> TC>0 channel the issue
predicted).  The matrix is no longer a suspect; max|EPS|/HU~25
rows exist but are sign-definite.  What remains is covered by
issue 4/4a below.

ISSUE 4 (non-contractive lagged map) -- CONFIRMED, refined:
- MOM_KNOP_REFRESH diverges EVERY iteration: best residuals
  13.9 -> 1.3E9 -> 6.5E10 -> 27.7 -> 1.35E9 -> 1.06E3 -> 36 -> 32
  -> 21 (it1-9), mid-pass explosions to ~1E20 before the abort
  guard; ~850-940 negJ/pass; unconverged -> pre-refresh field
  retained each time.  The sign-aware EPS work did not cure it.
- Mechanism is NOT the matrix: all rows are M-matrix and ~98%
  diagonally dominant.  The failures are faithful negative
  solutions to NEGATIVE RHS (XM_SAVE<0 at seed cells from the
  explicit EPS_PREV/PSIPREV/H-coupling terms), which then
  propagate through non-dominant rows.

NEW SUB-ISSUE 4a (VERIFIED): near-singular transparent rows.
The shallow compressing band (expanded depths 18-26,
sigma~-1.006, K_ON_J~0.99, chi~1E-11 i.e. transparent) has
TB ~ |TA|+|TC| -- the J-equation there is nearly a pure
recursion.  A small negative RHS seed amplifies across the
block into a saturated -285 plateau.  The deep shocked shell
(d457-473, sigma~-4.7, K_ON_J~0.2) is a separate band fed by
coupling rather than singularity.  Two morphologies, same
map-level cause.

DIAGNOSTIC ARTIFACT (VERIFIED, cosmetic): the GS sweep summary
counters (nonMM/negEPS/max|EPS|/HU/capped/bnd floors) are wiped
whenever the last solver call of the iteration is FREQ_INDX=1 --
the reset at FREQ_INDX==1 runs inside that call AFTER the counters
accumulate.  Pre-requeue refresh aborted on odd (ascending)
passes -> last call NCF -> real counts (23M negEPS).  Post-requeue
it ran all 60 passes -> ends on freq 1 -> counters show only the
freq-1 call's own values (argmax (262,1) confirms).  EPS coupling
did NOT vanish; the sweep-summary line is simply unreliable.
MOM_NEGJ_LRG/MOM_GS_MAXCHG (which drive GS_FREEZE) are restored
from pass totals at :2447-2449, so the gate itself is not fooled.

F5 consequence confirmed: MOM_MAXCHG (refresh best residual)
never nears the 0.05 maturity threshold -> GS_FREEZE holds ->
zero BA iterations for TypeIIn, permanently.  The gate is doing
its job here -- the field genuinely is immature -- but it means
this model family cannot converge via the MOM path at all.

Standing verdict (unchanged, now stronger): for TypeIIn-class
strong compression the MOM path must not be trusted; USE_LAM_ES=T
(FG path) remains the only converged option.  A real fix needs to
target the map: damp the compressing-band iterate (lower
JNEXT_RELAX locally), subcycle the shell cells, or solve the
coupled frequency block simultaneously instead of lagging it.
Row-level guards cannot help -- the rows are already healthy.

====================================================================
AUDIT 24-Sep-2026 (late pm) -- is USE_LAM_ES=T (FG path) suitable
for the strong-compression model?
====================================================================

What the switch does (VERIFIED in source):

- comp_j_blank.f:770-772: RJ(1:ND)=TC(1:ND) -- J comes straight from
  the ray moments; no moment solve runs.  Under DO_NONMONO_KNOP the
  rays are solved by FG_J_CMF_NONMONO_V1 with the correct compressing
  treatment (chi_hat split, sign-aware I_PREV/I_NEXT upwind, lagged
  I_NEXT store + FG_KNOP_REFRESH_V1).
- do_var_cont.f:332-333: `ELSE IF (USE_LAM_ES)THEN` is an EMPTY
  branch (upstream design) -- TX/TVX stay zero, so the SE/BA update
  gets NO J-variation: a pure Lambda iteration for populations.
- do_var_cont.f:333-357: the DO_NONMONO_KNOP branch that fills the
  diagonal ALO (TX(K,K,1)=dJ/dCHI, TX(K,K,2)=dJ/dETA from
  GET_FG_ALO_DIAG) sits AFTER the USE_LAM_ES test and is unreachable
  under USE_LAM_ES=T.  The FG solver computes ALO_DJDCHI/ALO_DJDETA
  on every call regardless (fg_j_cmf_nonmono_v1.f:3460-3527,
  ALO_VALID) -- the acceleration exists but is discarded.
- comp_j_blank.f:794-798: HBC_CMF(1) is stale under LAM_ES (moment
  solver never fills it); HFLUX_AT_OB falls back to HNU_AT_OB only
  under ZERO_FLUX + KNOP.  Other inner-boundary methods would keep
  the stale value.
- fg_j_cmf_nonmono_v1.f:88-95: the I_nu-proportional comoving term
  is omitted by convention (same as stock V13/V11) -- documented,
  but at sigma~-76 it is a caveat to carry, not an error.

Runtime evidence (VERIFIED):

- The FG lagged map IS contractive on TypeIIn where the MOM map is
  not: TypeIIn_KNOP refresh converged in ~6 passes (residual 2.7E-2,
  0 diverged rays, 0 clamps); KNOP_ES decaying 0.15->0.057;
  ES_ION 0.11->0.055.  The ray-level fixed point converges.
- No FG-path run has yet reached the convergence criterion on this
  model: TypeIIn_KNOP stopped at it21 (external NUM_ITS=0, still
  contracting at MAXCH~10%); KNOP_ES (job 11962582) and
  KNOP_ES_ION (job 11961300) in progress, all iterations so far
  healthy (BA computed at it10/it20).
- ~1000-1500 s/iter at 4 CPUs; refresh is ~70% of that (R1 in
  [[future_fixes]]: cpus_per_task 16 is the cheap win).

Assessment:

- Suitable in physics: J from the rays is the correct comoving mean
  intensity including compression; better founded than the MOM path
  while the J_NEXT map fails to converge.
- Limitation 1 (biggest): pure Lambda SE update -- unaccelerated,
  so thick-region convergence is slow and could stall or
  limit-cycle.  The fix is already in the code but unreachable:
  reorder do_var_cont so DO_NONMONO_KNOP fills the ALO diagonal
  before/inside the USE_LAM_ES branch (~10 lines, no new machinery;
  same fixed point, just an accelerated iterate).  Candidate F8.
- Limitation 2: no row-level J diagnostics on this path -- only
  N_NEG_INT/N_NAN_JNU clamp counters; interior field damage is less
  visible than on the MOM path.
- Limitation 3: convergence to tolerance not yet demonstrated --
  needs one run taken to completion before it is called the
  production path (project rule: stability is not convergence).

F8 implemented (24-Sep-2026, commit 2e6fcf0): the
DO_NONMONO_KNOP/DO_NONMONO_SOB branch now precedes the empty
USE_LAM_ES branch in do_var_cont.f, so KNOP+LAM_ES=T fills
TX(K,K,1)=dJ/dCHI, TX(K,K,2)=dJ/dETA from GET_FG_ALO_DIAG -- the
SE update becomes a diagonal ALI rather than a pure Lambda
iteration.  SOB unaffected (no ALO store; explicit zeroing ==
stale-zero state).  Plain LAM_ES unchanged.  Rebuilt into
exe/cmfgen_dev.exe (18:19 build; the plt_jh.exe link failure in
spec_plt is a pre-existing gap -- libsubs now pulls
control_variable_mod via the mom solver's USE_LAM_ES import, and
spec_plt does not link -lmod; unrelated to this change).

Validation: TypeIIn_10d_DEV7_KNOP_ES_ALO (config committed
f28c684, identical physics to TypeIIn_KNOP_ES -- COH_ES=F, FIX_T,
FIX_NE, DO_NG=F -- fresh independent model dir), FINAL 11964585 /
OBS 11964586 on small (24h).  Success criterion: the ALO-coupled
SE update reaches the EPS_TERM=1.0 convergence assessment where
pure-Lambda ES runs have not yet been taken to completion.  Side
effect: the running KNOP_ES (11962582) and KNOP_ES_ION (11961300)
jobs pick up the reorder binary at their next requeue.

F8 integration audit (24-Sep-2026, commit 18e498a adds the
runtime confirmation line): the ALO wiring was verified against
the Knop ALI requirements --

- Operator consistency: ALO_DJDETA/ALO_DJDCHI/LAMSTAR_DIAG are
  accumulated in fg_j_cmf_nonmono_v1.f:3496-3527 from the same
  JQW ray weights and the Knop-split chi_hat the formal solution
  uses.  LAMSTAR_DIAG = dJ/dS' is sanity-checked in (0,1] on
  INIT; ALO_VALID is set only when the INTEGRAL method ran, so
  DIFFERENCE-method frequencies cannot inject a mismatched
  operator.
- Per-frequency correspondence: ALO arrays are depth-only and
  re-zeroed per formal call.  DO_VAR_CONT runs in the same
  per-frequency block immediately after COMP_J_BLANK
  (cmfgen_sub.f:1874->1975 and :2365->2565), so GET_FG_ALO_DIAG
  returns the ALO of the frequency being varied -- not a stale
  or last-frequency value.
- Correct slot: TX(K,K,1/2) are the per-frequency dCHI/dETA
  temporaries folded into slots 3/4 (VCHI/VETA chain rule),
  slot 5 (electron scattering), the line-variation matrices
  and slot 6 (T) at do_var_cont.f:443-534 before zeroing at
  :542 -- the identical consumption path as VAR_MOM_J_CMF_V12.
- Ordering: DO_VAR_CONT is gated on COMPUTE_BA .AND. .NOT.
  LAMBDA_ITERATION, so the fill only runs when the BA matrix is
  actually assembled.  FG_KNOP_REFRESH_V1 runs after the sweep
  (cmfgen_sub.f:2937); the ALO reflects the pre-refresh field,
  which is acceptable since Lambda* is opacity/geometry-derived
  (captured opacities) and the lagged upwind-source term is an
  outer-iteration quantity by construction.
- Approximation level: diagonal in depth and frequency.  The
  frequency coupling (compressing upwind source) is iterated,
  not in the operator -- the standard diagonal-ALI convention.
  DO_SRCE_VAR_ONLY=F in this config, so the fill is not
  overwritten.
- Runtime proof: 'KNOP ALO: dJ/dCHI diag [...] dJ/dETA diag
  [...]' prints once per BA sweep (18e498a).  Job 11964585 was
  requeued at 18:3x so it execs the diagnostic binary; first
  BA expected after the warm-up once MAXCH<50.

24-Sep-2026 (evening): per user request ALL cmfgen jobs were
cancelled (KNOP_ES 11962582, KNOP_ES_ION 11961300, TRIDIAG
11964466, TRIDIAG_bndflr 11962018, M15 cohgd 11963093, ALO
11964585, plus their OBS/dependent jobs) and a single fresh
USE_LAM_ES=T run was submitted on TypeIIn_KNOP: FINAL 11964791
/ OBS 11964792 (interactive, config submission=debug, 30G
clamp -- same conditions the previous it1-21 run survived).
The previous run dir was archived to
cmfgen/10.00days_prev_it21, so this is a true cold start
(no POINT/FGINEXT_GS/MOMJNXT_GS state).  It runs the 18:33
cmfgen_dev.exe carrying both the F8 reorder and the KNOP ALO
diagnostic line.  This is now the sole validation vehicle:
first BA (KNOP ALO line + BA computed) expected after the
warm-up iterations.

A second identical-physics run was submitted on the normal
partition per user request: TypeIIn_10d_DEV7_KNOP_LN (config
TypeIIn_10d_DEV7_KNOP_LN.yaml, model dir TypeIIn_KNOP_LN,
stella inputs symlinked from TypeIIn_KNOP/stella), FINAL
11964798 / OBS 11964799 on small (23:59 limit, 95G -- no
OOM risk unlike the debug run's 30G clamp).  It runs
uninterrupted in 24h chunks vs the debug run's 2h requeue
cycle; both share the same config physics and the F8
binary.

25-Sep-2026: MOM refresh optimization round (jobs 11966433 gsdesc,
11967320 keepneg, 11967617 gatebulk; commits 28b0cdd, b5ff320,
7b1e570).  Findings and changes:

* Frozen JNEXT_RF_READ on descending passes was tested live
  (MOM_RF_GS_DESC=T): residual floor ~40 at it1, tracking baseline
  (13.9/10.6/5.5 vs 40.8/12.5/5.6) -- roughly neutral, not the
  limit-cycle driver.  Kept as VADAT option, default F.
* MOM_RF_KEEPNEG=T quarantines repaired cells (negJ reset,
  non-finite, cap) -- the stored J_NEXT_FRQ keeps the prior iterate.
  First implementation inflated the reported residual (quarantined
  cells pinned the metric at ~3.5E3); fixed by counting them
  separately as MOM_GS_NQUAR.
* Same-direction stall exit added to MOM_KNOP_REFRESH (two
  consecutive same-direction residual ratios >0.90 after pass 8):
  fired at pass 9-16 on the limit cycle, cutting ~420s of wasted
  passes to ~70-100s per iteration.  Unconditional improvement.
* Bulk metrics MOM_GS_NBAD5/NBAD50/NTOT/NQUAR now recorded per
  sweep and per refresh pass (cells>5%, cells>50%, quar= lines).
* MOM_GATE_BULK=T (VADAT, default F): BA maturity judged by bad-cell
  fraction (<=0.1% cells >5%, zero capped, maxchg<0.5, negJ/quar
  <10% of cells) instead of max-norm<5% + negJ>0.  The default gate
  cannot open on this model: negJ_LRG ~750/pass in refresh and ~9E5
  in-sweep keep it permanently frozen even though the bad fraction
  is ~0.002-3% -- M15 bndflrF showed the same (converged it74 in
  Lambda with the negJ gate frozen).
* Baseline behavior across great iterations on TypeIIn_10d:
  MOM best residual 13.9 -> 10.6 -> 5.5 -> 4.9 -> 5.3 -> 3.1 ->
  3.9 -> 3.1 -> 2.8 (bndflr it1-9): slow contraction, plateauing
  near ~2.7-3.9, far above the 5% max-norm gate.  The great-
  iteration loop contracts; the inner refresh oscillates.

25-Sep-2026 (later): sweep-drift maturity gate (commit 069ed05,
job 11969828 gatebulk resubmit).  The first GATE_BULK version still
gated on the refresh-pass residual counters (MOM_GS_*), which measure
lagged-map self-consistency: nbad5 decayed ~7-10%/iter but
asymptoted, and refresh maxchg plateaus at ~1.3-1.7 on marginally
contractive compressing cells.  Fix: MOM_SW_* snapshot of the
sweep-level counters at refresh entry (fresh XM vs stored J_NEXT --
the drift the populations actually see).  Gate now uses sweep drift:
freeze while >10% of cells drift >5%, unsafe on capped writes,
>0.1% of cells drifting >50%, any cell >10x stored, negJ/quar
fractions >10%, or refresh residual >1e2.  NOTE: a resubmit on an
immature checkpoint replays the 4-iteration warm-up and rebuilds
J_NEXT -- field maturity progress is lost (populations carry over);
avoid resubmitting immature runs.

Result: gatebulk opened the freeze at it25 (nbad5_sw crossed 10%);
it25-27 ran clean Lambda under DO_LAM_AUTO with corrections
17->15->14->13%, no immature-field returns, nbad50 ~300, negJ/quar
~2000.  First time the USE_LAM_ES=F path has exited the forced-
Lambda freeze on this model; the iteration sequence now matches the
T path's (warm-up -> field maturation -> DO_LAM_AUTO mode choice).
Winning keys ported to TypeIIn_10d_DEV7_KNOP_TRIDIAG.yaml
(MOM_RF_KEEPNEG=T, MOM_GATE_BULK=T; GS_DESC left F -- neutral).

Post-gate observation (11969828, it28-36): first BA computed at
it28 (kick +3.4e3/-1.3e17%); thereafter the model cycles
Lambda<->BA like the T path -- each BA kick trips the sweep-drift
unsafe bound (maxchg 33-54, nbad50 ~2.4e4) for 1-2 frozen Lambda
iterations, then re-opens.  Kick sizes so far 3.4e3 -> 1.5e4 ->
2.7e4 (it28/31/35) -- bounded and absorbed each cycle but not yet
shrinking; whether the BA cycling nets faster convergence than the
pure-Lambda baseline (F-full 11965659 contracting ~10%/iter at
it40+) is the open question for the rest of this leg.

## Status update (30-Sep-2026) -- stall-point audit + freeze escalation

A sixth external review audited the KNOP stall paths.  Fixes landed
(commits a9ab6b7, be1213e, 72641f2, d57d105, 5e7d5cd, 7b07fdf; shell
watchdog ed48071 in the main repo): exit-side hysteresis on the
compressing-depth set (rebuilds were re-zeroing I_NEXT on a flapping
sigma=-1 boundary), in-loop bounds on the five ray-bracket searches,
persistent-quarantine age counters (MOM_GS_NPQUAR), a BA_NREJ_LIM cap
on consecutive rejected BA steps, the GS_STATE freeze-state file plus
KNOP_STALL_MAX=60 watchdog in the generated slurms, and the
MOM_TRIAL_MIN=30 marginal-release BA trial (resacc-gated) for fields
that asymptote but never reach the strict maturity bound.

Full detail: [[convergence_issues]] section 11.8 and the 2026-09-30
dated entry.  The marginal-release path has not yet taken a BA step
on a live model -- runtime validation pending.

## Status update (30-Sep-2026, late) -- direction-ordered refresh + EPS limiter

Seventh round, targeted at the m15_2p82d MOM stall (job 12045050, the
fully-compressing m15 model: sigma<-1 at all 70 depths, inner jump
sigma ~ -157).  Literature check: Knop et al. Paper II (A&A 509,
arXiv:0909.4264) Sect. 3.2 solves the wavelength coupling as a
boundary-value problem with a physically-guided GS pass order --
follow the direction information propagates.  Ceulemans et al. 2025
(JQSRT 343) put a stability criterion on the per-bin Doppler shift
(matched by DO_NONMONO_FREQ_MATCH).

Changes (commit 5b814ae, on top of c02664c which carries the PSIPREV
positivity fix + GS trial-release machinery that job 12045050 ran):

1. MOM_RF_DIR = AUTO/ALT/ASC/DESC (default AUTO).  Per-sweep census
   of compressing J-rows (Q_NM<0, the actual upwind selector);
   >=90% one-sided -> refresh runs only the dominant direction.
   This grid measured cmp_frac=0.765 -> alternating stays (a mixed
   field genuinely needs both directions: A passes freshen expanding
   reads, D passes freshen compressing reads).
2. Compressing-row write suppression in the main sweep while a
   dominant-direction refresh owns the lagged store (needs
   MOM_RF_ENABLED + KN_M_READY; inert on ALT/mixed grids).
3. MOM_EPS_LIM (default T): donor-cell positivity limiter on the
   lagged EPS_PREV*J_up midpoint scatter -- a row may not be debited
   past its accumulated positive RHS; recipient credited only what
   was debited (conservative pair).  Applied to interior scatter,
   negated boundary-row share, and the catastrophic-fallback rebuild.
   Target: the negJ population at advection-dominated transparent
   cells (|EPS|/HU ~ 20, chi << chi_hat).

Job 12045050 cancelled (stalled, n_frozen 38/60); resubmitted as
12046369 warm-restart (OBS chain 12046373).  Binary =
ptmp2 exe/cmfgen_dev.exe built 22:17 from ptmp1@5b814ae file-sync
(ptmp2 git clone cannot push/fetch origin from the login node --
commit objects live in the ptmp1 tree only; see
m15_2p82d_knop_MOM_apu/cmfgen/2.82days/RUN_NOTES_12046369.md).

First refresh of the resubmit: eps_cap ~ 2.7-2.9M capped transfers
per pass, pass negJ ~2.0e4 (was ~7.5e4-class), sweep negJ 330k
(was 2.4M), DGTSV/TB<=0/nonfinite all zero, HBC floors 0.  Residual
still pinned ~3.6e5 and sweep maxchg 8.1e17 (single bluest-edge
cell) -- the limiter fixed the sign-defect class, the coupling
contraction itself is still the open problem for this model.
Verdict pending: watch negJ/quar/residual trend to it50 before the
watchdog fires (n_frozen=60).
