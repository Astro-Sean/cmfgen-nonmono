# Convergence strategy audit: dev7 non-monotonic solver vs the Knop/PHOENIX literature

Date: 2026-10-31.  Scope: what convergence criteria and iteration
ordering the non-monotonic CMF literature actually uses, how dev7's
machinery (`FG_KNOP_REFRESH_V1`, `MOM_KNOP_REFRESH_V1`,
`SOLVE_FOR_POPS` gates, `SOLVEBA_V13` resacc) maps onto it, and the
adapted strategy implemented in this tree.  Claims tagged
VERIFIED / SUPPORTED / UNVERIFIED / CONTRADICTED / BLOCKED per repo
policy.

Sources reviewed:

- Baron & Hauschildt 2004, A&A 427, 987 (astro-ph/0408212) -- the
  original non-monotonic CMF BVP algorithm ("Paper I" in Knop-II's
  numbering).
- Hauschildt & Baron 2004, A&A 417, 317 -- wavelength-derivative
  discretizations (first vs second order) and the xi-mixing
  stabilization.
- Knop, Hauschildt & Baron 2009, A&A 496, 295 (aa10290-08) --
  Paper I: positive-generalized-opacity formal solution.
- Knop, Hauschildt & Baron 2009, A&A 501, 813 (aa11793-09) --
  Paper II: large-scale implementation; iteration ordering and
  convergence criteria.
- Ceulemans, De Ceuster & Decin 2025, JQSRT 343, 109470
  (arXiv:2505.08309) -- oscillatory-eigenmode analysis and the
  frequency-matching stabilization.

## 1. What the papers converge, in what order

### 1.1 Baron & Hauschildt 2004 (the reference algorithm)

The non-monotonic wavelength derivative turns the CMF RTE from an
initial-value problem into a boundary-value problem in (i,l).  Their
algorithm:

1. At fixed source function (frozen chi, eta), assemble the FULL
   wavelength-space discretization on each characteristic: for a
   second-order-in-s derivative the row couples (i-1,i,i+1) x
   (l-1,l,l+1) -- a large banded sparse system.
2. Solve it EXACTLY with a general sparse solver (SuperLU).  There is
   no inner iteration at all; the BVP is converged by construction
   before any operator-splitting step.
3. ALI step: `(1 - eps*Lambda*) J = J^FS - eps*Lambda* J^old`, with
   Lambda* tridiagonal in wavelength and full spatial bandwidth,
   iterated to convergence.
4. Repeat the outer (great) iteration until the model converges.

VERIFIED (astro-ph/0408212 text): "Since the equation is first order
in wavelength at each spatial point there is only one boundary
condition, whose wavelength sense depends on the local sign of the
coefficient"; local upwind differencing is mandatory for stability.

### 1.2 Hauschildt & Baron 2004 (xi-mixing)

Two wavelength discretizations (first-order upwind; second order) are
mixed by xi.  VERIFIED from Paper II Sec. 2: fully implicit xi=1 is
required for strong wavelength couplings -- "the interpretation of
the entire wavelength discretization as a source term is not valid
for all optical depths".  dev7 already uses xi=1 throughout.

### 1.3 Knop et al. 2009 Paper II (the method dev7 mirrors)

Because the full matrix solve does not scale, Paper II replaces it
with a PHYSICALLY ORDERED Gauss-Seidel march (Sec. 2):

- March along the characteristic in space (the direction photons
  travel).
- At each spatial point march in wavelength in the LOCAL upwind
  direction, i.e. in the sign of a_i (expanding: short->long
  wavelength; compressing: long->short).
- With a linear source-function interpolant the march is
  "quasi-analytic": every term on the RHS is already computed, so
  each cell needs exactly one update -- a direct solve, not an
  iteration.
- With the general (parabolic) scheme they iterate Eqs. 31/32; the
  only lagged terms are the C-corner coefficients C^+ I_{i-1,l+1},
  C^- I_{i+1,l-1} -- the SPATIALLY DOWNWIND x WAVELENGTH-UPWIND
  corner.  Everything on the (i,l) path the information actually
  travels is fresh.

SUPPORTED (Sec. 4.3, + A&A 496 formal-solution paper): a converged
field is then fed to the ALI step; the ALI itself iterates
`(1-eps*Lambda*)J = J^FS - eps*Lambda*J^old` with Lambda* tridiagonal
in wavelength x full spatial bandwidth, solved by blockwise GS/Jacobi
in wavelength with a one-point overlap, and converged to the internal
relative tolerance 1e-8.  "ALL wavelength points have to be converged
simultaneously" -- the old recursive scheme that freezes the previous
wavelength and never revisits it is explicitly superseded.

The NLTE-relevant ordering the paper guarantees and dev7 must
reproduce:

    frozen (chi, eta) -> CONVERGED wavelength-coupled BVP field
                      -> source/SE update (Lambda* / populations)
                      -> repeat

The populations never see a non-fixed-point field.  That ordering --
not any particular pass count or tolerance -- is what makes the
scheme safe to couple to statistical equilibrium.

### 1.4 Ceulemans et al. 2025 (why it still stalls)

VERIFIED (Secs. 2.3, 3): the standard discretization admits
spurious eigenvalues in -1 < lambda < 0 (oscillatory modes) when the
per-cell comoving Doppler drift exceeds the bin spacing, and those
modes survive iteration ordering alone.  Their fix: match the
frequency grid along the ray so the upwind sample sits on the
characteristic foot -- dev7's `DO_NONMONO_FREQ_MATCH` (currently
default F) is the same idea applied to the compressing slot.  Their
analysis also confirms the marginal-link mechanism already documented
in [[MOM_KNOP_STIFFNESS_AUDIT]]: links with
`|Q|/(chi+|Q|) -> 1` have gain -> 1 under any ordering that keeps
them lagged.

## 2. What dev7 currently does, mapped onto the papers

| Paper element | dev7 equivalent | Status |
|---|---|---|
| BVP solved exactly (B&H04) / quasi-analytic ordered march (Knop II) before any SE step | In-sweep: expanding cells GS-live, compressing cells Jacobi-lagged on the I_NEXT depth store (one lag per great iteration).  Post-sweep `FG_KNOP_REFRESH_V1` iterates the lagged field with ALTERNATING GLOBAL frequency passes (ascending: expanding links fresh; descending: compressing links fresh). | CONTRADICTED in part: the alternating l-major pass order leaves one coupling direction stale EVERY pass; the corner that Knop lags (spatial-downwind x upwind-frequency slope term) is NOT the link dev7 lags. |
| Ordered march converges in a few iterations | `KN_RAY_MAXPASS=60` alternating passes, adaptive omega, stall exits | Contradicted on compressing rays: marginal |Q|-links keep the pass-residual pinned; refresh routinely exits on stall rather than tolerance. |
| Converged BVP before SE update | `SOLVE_FOR_POPS` gate: GS_FREEZE forces LAMBDA while the field is immature; marginal tiers allow bounded resacc-gated BA trials | SUPPORTED in intent: the gate is the right NLTE-safety structure.  It fails open only if the "maturity" metric passes while the field is not self-consistent. |
| Lambda* tridiagonal in wavelength | Diagonal ALO only (TX(K,K,1/2) from FG or MOM diagonal in `do_var_cont.f`) | SUPPORTED as a documented deviation: cheaper, weaker; the upwind-frequency sensitivity is dropped.  Acceptable for step 1; revisiting means cross-frequency Jacobian blocks -- out of scope. |
| All wavelengths converge simultaneously | Refresh iterates all NCF frequencies together | VERIFIED consistent. |
| Frequency matching (Ceulemans) | `DO_NONMONO_FREQ_MATCH` (default F) | Implemented, off by default. |

### 2.1 The structural deviation that matters

dev7's refresh iterates `DO L=LB,LE,LDI` over frequency with the
whole ray's depth march inside each frequency (l-major, alternating
l-direction per pass).  On an ascending pass every compressing
upwind read `FLD(JPR, l+NFM)` is stale; on a descending pass every
expanding read `FLD(i, l-1)` is stale.  Strong links therefore
carry a half-lagged GS map whose convergence factor approaches the
stronger of the two link weights -- exactly the observed stall at
|Q|/(chi+|Q|) -> 1.

Knop's ordering is the transpose: march i along the characteristic;
at each row i march l in the row's upwind direction.  Under dev7's
stencil that gives:

- Expanding reads `FLD(i,l-1)` are SAME-ROW: fresh under an
  ascending in-row l-march, unconditionally.
- Compressing reads `FLD(JPR(i), l+NFM)` are CROSS-ROW through the
  depth-grid interpolant: fresh whenever the donor rows
  JPR1/JPR2(i) have already been written this pass; one of the two
  bracketing donors is typically still unmarched -> this stays a
  corner term, with weight <= |Q|*(1-WINX) or |Q|*WINX.
- The parabolic-slope terms DS(j), DS(j+1) need the source of the
  NOT-YET-MARCHED neighbor row (SP(j+1) for the inward march;
  SP(j-1) for the outward march).  That is precisely Knop's lagged
  C-corner: the only mandatory lag, weighted by the slope
  coefficients A3/A4 and killed entirely where the monotonicity
  clamp drives DS -> 0 (i.e. in the steep-gradient compression zones
  where it matters most).

So the ordered march replaces "half of ALL links stale every pass"
with "one weak spatial corner lagged".  Expected behaviour:
no-compressing rays exact in one pass (unchanged); compressing rays
converge in a few passes at a rate set by the corner weight instead
of by the full |Q| link.  SUPPORTED by the derivation above;
UNVERIFIED until runtime-tested.

### 2.2 NLTE-safety requirements vs current state

What "NLTE safe" requires, translated from the papers' ordering:

1. The SE/BA update must act on a field that is the fixed point of
   the transfer map at the CURRENT chi, eta -- not merely a field
   whose drift is small.  dev7 approximates this with the maturity
   gate (NBAD/NTOT fractions, MAXCHG, NCONV store churn).  The gate
   semantics are correct; the weakness is that the refresh rarely
   reaches the fixed point, so "mature" measures seed drift rather
   than self-consistency.  The ordered march closes that gap for FG.
2. Low-intensity floor cells must not dominate convergence metrics.
   VERIFIED already handled: the 1e-4 column-peak significance floor
   in MARCH/WRITE_STORE.
3. Catastrophic cells (non-finite, capped, negative-J repairs,
   quarantines, boundary-floor events) stay visible to the gate.
   VERIFIED in the MOM gate tiers; the FG gate reads NCONV/NCAP.
4. A bounded number of BA trials on a marginal field, each
   individually guarded by residual acceptance + trust radius.
   VERIFIED: `DO_BA_RES_ACC`, `BA_TRUST_CUR`, `BA_NREJ_LIM`,
   `GS_TRIAL_LIM`, `MOM_TRIAL_MIN` refund on acceptance.
5. Stale-Jacobian corrections bounded.  PARTIAL: the trust-radius
   zeroing in SOLVEBA_V13 only runs under `DO_BA_RES_ACC`; with
   resacc off the "BA not computed" path can emit unbounded
   corrections from a stored inverse (the 1.1e12 event in
   [[convergence_issues]]).  Recommendation: keep DO_BA_RES_ACC=T on
   all KNOP runs; do not run it0 scale-ups without it.
6. Iteration-ordering caveat: the refresh improves the LAGGED STORE
   that the next sweep reads; the J used by this iteration's SE
   update comes from the sweep.  A converged store therefore means
   "next sweep's J is a fixed point" -- the gate must test store
   drift (NCONV) AND sweep residual together.  Current gate does
   both.  SUPPORTED.

### 2.3 The MOM path

The per-row ordering does not transfer: `MOM_J_CMF_NONMONO_V1` solves
a GLOBAL tridiagonal system in depth per frequency -- there is no
depth-march to reorder.  Its Jacobi lag is inherent to the
per-frequency decomposition.  Options consistent with the papers:

- Keep direction-ordered passes + Anderson acceleration (current
  machinery) and the stricter gate tiers -- the MOM path is a
  moment analogue, not the paper's operator, and OBS only supports
  `USE_LAM_ES=T` anyway.
- If MOM convergence is needed, the paper-conform direction is the
  block/frequency-range iteration of [[BLOCK_GS_DESIGN]]
  (isolate the strongly coupled frequency runs), not more global
  passes.

Recommendation: `USE_LAM_ES=T` (the FG/ray path) is the
paper-conform route and the one OBS supports; the ordered refresh
makes it converge; MOM stays on the gated machinery as documented.

## 3. Adopted strategy (this tree)

1. `SWEEP_FLD_ORD` -- a new pass driver in `KNOP_RAY_SOLVE_V1` that
   marches rows in the characteristic direction (j ascending for
   I_M, descending for I_P) and marches frequency inside each row
   in the row's upwind direction (sign(GAM): ascending on expanding
   rows so the same-row `FLD(j,l-1)` read is fresh; descending on
   compressing rows so a self-mapped `FLD(j,l+NNF)` read is fresh).
   A four-deep ring `SPW(l,mod(r-1,4)+1)` holds the modified source
   S' of the four rows a march step touches; the not-yet-marched
   neighbour's row is previewed once per step by `KN_ORD_SPW_FILL`
   from its pre-march field -- the Knop C-corner.  `KN_ORD_SPV`
   reproduces `SRC_PRIME` elementwise including the
   `DO_NONMONO_FREQ_MATCH` remap of `MATCH_FREQ` (inline NNF2 =
   ceil(DRF_L/KN_DLN_N), capped on-grid); `KN_ORD_SEC`/`KN_ORD_DS`
   reproduce `SLOPES` on the ring including the endpoint
   extrapolations and monotonicity clamps.  Same fixed point, same
   residual metric, same guards (INEXT cap, NaN floor, omega
   damping, divergence abort, stall exit).  Gated on VADAT
   `KN_RAY_ORDERED` (default T); F restores the legacy alternating
   passes.  Requires the coefficient pack and NI>=3; the
   alternating scheme is the automatic fallback otherwise.
   Corner case: a ray with no compressing cells still runs two
   ordered passes, because pass 1 carries the seed's lagged slope
   preview (the alternating scheme's single-pass exact-exit no
   longer applies).
2. Gate semantics unchanged: converged refresh + bounded store churn
   still gate BA; resacc + trust radius still guard each BA step.
   With the ordered march the FG gate's "RES<TOL" becomes a real
   fixed-point criterion rather than a stall exit.
3. `DO_NONMONO_FREQ_MATCH` stays opt-in; recommend enabling it on
   compressing models (it weakens the |Q|->1 links at the source).
4. No change to the MOM path beyond what is already in place; the
   ordered march does not apply to a global tridiagonal solve.

### Validation plan (before claiming PASS)

- `make cmfgen` clean.
- Model run through `scripts/run_model_grid.py` on the non-monotonic
  test model with `DO_NONMONO_KNOP=T, USE_LAM_ES=T,
  KN_RAY_ORDERED=T`: expect refresh NPASS to fall from the stall
  limit to a handful; GS_NBAD/NCONV to reach the mature tier so BA
  releases without trial pulses; luminosity drift absent.
- A/B against `KN_RAY_ORDERED=F` on the same setup to isolate the
  ordering effect.
