# KNOP remediation WP-7 -- trace levels and the ill-conditioned SE block

Date: 26/27-Nov-2026. Status: LANDED pending push (flags default F;
armed arms are the human A/B). Three independent defects, one package.

## WP-7a -- POP_TRACE_MAXFAC inside the RES_ACC/LAMBDA gate

The per-row fractional-correction cap POP_TRACE_MAXFAC sat inside
`IF(DO_BA_RES_ACC .AND. .NOT. LAMBDA_IT)` in SOLVEBA_V13, so on the
all-LAMBDA KNOP stall (every call is a LAMBDA call while the warmup
runs) a runaway trace row was never clamped: the demanded correction
kept growing and pinned MAXCH at the trust cap. The block now runs
unconditionally, immediately after the STEQ_RAW snapshot and
POP_FLOOR_STEP, so the pre-clamp demand stays visible in the
STEQ_RAW-based channels (BA_DEM_EXM) exactly as before. Structural
check in test_knop_lagged_field.py asserts the block sits outside
the gate and that no clamp remains inside it.

## WP-7b -- the DIAG solve contract (DGEEQU INFO + DGECON)

Two holes in CMF_BLKBAND_V3, both closing the "solve contract" gap
already logged open in docs/CHECKLIST.md:

- DGEEQU's IFAIL was ignored on both solve paths. On TRI it now
  routes to the 9999 soft-failure path (the TRI fallback machinery
  retries with the diagonal blocks only). On the DIAG compute path
  it counts the depth failed, zeroes that column and moves on --
  the same handling DGETRF already had.
- A failed DIAG depth previously returned FLAG=.TRUE. with a zeroed
  STEQ column, so a caller could not distinguish "solved" from
  "partly failed". NDFAIL now counts failed depths and the routine
  returns FLAG=(NDFAIL==0); SOLVEBA_V13's `!SUCCESS -> STOP` path
  then halts instead of certifying an incomplete grid. This is a
  non-flag-gated contract fix: failed solves were always wrong, they
  are now reported wrong.

DGECON also runs on every factored depth (was POP_AUDIT-only) so a
real reciprocal condition estimate exists as solve-contract evidence;
CN_ANORM is recorded unconditionally for the same reason.

Exposed while testing: CN_WORK was allocated 2*N but DGECON's WORK
must be 4*N -- the overflow smashed the heap and crashed the process
at exit once DGECON ran unconditionally (the POP_AUDIT gate had kept
it rare enough to miss). Fixed to 4*N.

## WP-7c -- symmetric elimination + conservation at the most
##         populated row

- New standalone routine CMF_BLKBAND_PIN_SYM_V1(C_MAT, RHS, XBAR,
  PIN, N, N_PIN): standard symmetric elimination for a prescribed
  unknown. All column carries RHS <- RHS - C_MAT(:,J)*XBAR(J) run
  first against the unmodified matrix (order-independent for multiple
  pins), then each pinned row/column is zeroed, the diagonal set to 1
  and RHS(J)=XBAR(J). Under [BA_SYM_PIN] (default F) the DIAG solve
  applies it to each level row J with POPS(J,K) <= MAX(1E-30,
  BA_DEM_PFLOOR) -- the same census and J=1..N-2 range as the
  POP_AUDIT floor count, so ED/T and the conservation rows are never
  pinned. The pinned correction solves as exactly zero, matching the
  post-solve POP_MASK_FLOOR result while leaving that mask untouched
  in SOLVEBA_V13. This is the correct form of the matrix pin whose
  row-only prototype crashed the coupled block (exit 174 on the m15
  fast KNOP run, job 12868544): removing the column as well keeps the
  near-null direction out of the factorisation.
- Under [BA_NC_MAXPOP] (default F), GENERATE_FULL_MATRIX_V3 places
  each species' number-conservation equation at the row of its most
  populated level (argmax POPS at that depth), per Arramy et al. 2024
  Sec. 2.1.1 -- conditioning improves when the conserved equation
  sits on a well-populated unknown. The displaced rate equation
  relocates to EQ_SPECIES (the top pseudo-ion's ground row, whose
  rate equation was already sacrificed to NC in the default path), so
  the equation count is unchanged. Ground-state rows (ATM(ID)%EQXzV)
  are excluded from the candidate set: the ionization-conservation
  replacement below the insertion point may overwrite them and would
  silently drop the NC row. When the argmax row IS EQ_SPECIES or no
  candidate exists (single-level ions), the placement reduces to the
  default.

## Verification (all actual outputs)

- test_cmf_blkband_pin.py: 13 sym-pin checks added to the existing
  driver -- manufactured near-null trace column (coupling 1E-43):
  the pinned block solve stays finite, the three majors match a
  reference 3x3 solve to <1E-12, the pinned correction returns the
  prescribed value, and a nonzero XBAR case verifies the column
  carry keeps majors exact. Null-pin checks unchanged.
- test_mat5pen.py: 9/0 including the new drv_cmf_diag battery --
  'ok' returns FLAG=T with finite corrections on both depths;
  'badrow' (one rate row identically zero, nonzero residual) trips
  DGEEQU IFAIL=3, reports 'DIAG solve failed at 1 of 2 depths' and
  returns FLAG=F. Depth-1 corrections identical to 'ok'.
- test_gfm.py: 21/0 including drv_gfm_nc -- 'ncdef' keeps NC at
  EQ_SPECIES (maxerr 0), 'ncmax' moves NC to the argmax row and
  relocates the displaced rate eq to EQ_SPECIES (maxerr 0 on every
  row and STEQ element), ED/T rows untouched in both.
- test_knop_lagged_field.py: WP-7a structural checks pass; the file
  also needed its L_CONT_READY pin updated to KN_RAY_TOL_EFF (the
  WP-6 rename had silently broken that check -- the semantic intent
  is unchanged). File still ends at the documented pre-existing red
  (KN-SEED report strings renamed upstream).
- test_vadat_keys.py: 136/0 (+2): both flags default F absent and
  land on CVM when armed.
- Flag-off battery: test_knop_frc 35/0, test_knop_ordering 50/0,
  test_knop_osgn 14/0, test_iter_conv_gate 18/0, test_conv_exit
  53/0, test_sigma_true 16/0 -- all identical.
- test_alteration_coverage.py: cmf_blkband_pin_sym_v1 registered
  under test_cmf_blkband_pin. The remaining manifest FAIL is the
  parallel session's uncommitted subs/two/Makefile work, not this
  change.
- build.sh / build_debug.sh: see the build log verdicts.

## Flag-off identity

WP-7a moves an existing block -- identical arithmetic, now run on
every call (the all-LAMBDA stall path gains the clamp it always
needed; non-LAMBDA paths are unchanged). WP-7b's contract change is
a bug fix on the stock path (a failed depth reporting success was
never a valid result); the only new overhead is one DGECON per
depth per call. WP-7c's two flags are F by default: with BA_SYM_PIN
off the DIAG path is byte-identical; with BA_NC_MAXPOP off the NC
insertion is the same two assignments merged into one loop.

## Owed measurement

- BA_SYM_PIN=T on a floor-heavy model (m15 stall): the DIAG solve
  should complete with fewer DGETRF/GECON failures and no change to
  committed corrections versus the post-solve mask.
- BA_NC_MAXPOP=T: watch the same model's DIAG solve rcond (the
  POP_AUDIT RUNAWAY line) -- the conservation row on a populated
  level should lift the block's reciprocal condition estimate.
