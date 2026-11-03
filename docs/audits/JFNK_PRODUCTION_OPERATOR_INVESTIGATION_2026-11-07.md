# JFNK production-operator correctness investigation (Phase 1)

Date: 2026-11-07. Tree: `cmfgen_jfnk`, HEAD `22d63df9` + the working-tree
changes listed in section H. Scope: establish whether the production
`JFNK_RESIDUAL_V1`/`JFNK_JVP_V1` pair is a sufficiently consistent linear
operator for the existing restarted GMRES, without introducing a new
preconditioner, Schur complement, FGMRES, or frequency solver.

## A. Production residual architecture

`JFNK_RESIDUAL_V1(N, X, Y)` (new_main/mod_subs/jfnk_cmfgen_interface_v1.f,
~line 1664) maps the correction vector `x` to the masked statistical-
equilibrium residual. Depth-major ordering, `idx = (k-1)*NT + j`;
`j <= NT-2` are population corrections in fractional units,
`j = NT-1` is the electron-density correction, `j = NT` is the
log-temperature correction (`x_frac(NT) = exp(x_T) - 1`).

Frozen-J residual (the documented operator):

    F(x) = STEQ_coll(n0*(1+x_frac)) + STEQ_rad(n0) + CMAT_RAD * x_frac

Collisional/non-radiative terms are recomputed at the perturbed
populations; the stored radiative term `STEQ_RAD_FLAT` is added
unchanged; `JFNK_BA_RAD_MATVEC_V1` contributes the linearized radiative
increment through the stored band matrix. Radiation quantities held
fixed: all J/photo-rate-weighted radiative terms (they enter only via
the stored term and the linearized band product). Dependencies still
active: collisional rates, auto-ionization, charge exchange, advection
(`STEQ_ADVEC_V4` or `STEQ_CO_MOV_DERIV_V3` under `DO_CO_MOV_DDT`),
electron-density terms, two-photon rates, level dissolution.

Coupled-J residual: `JFNK_COUPLED_LOC=T` routes the same residual
through `JFNK_SWEEP_CB`, registered by `cmfgen_sub.f` as
`JFNK_RAD_SWEEP_V1` (~line 5755). Each evaluation replays the continuum
transfer (and the line section under the Sobolev-coupled option) at the
perturbed populations, then assembles the same SE equations the outer
iteration uses, with `FREEZE_I_NEXT` armed so the lagged Knop
`I_NEXT`/`J_NEXT` stores are not modified. Frozen vs coupled are
therefore different nonlinear operators: frozen is a fixed-radiation
subproblem with a linearized radiative increment; coupled probes the
full SE+radiation map. They are not interchangeable, and the driver
already treats an initial `F(0)` mismatch differently per mode
(coupled: adopt evaluated residual as the Newton baseline; frozen:
skip the Newton solve).

State handling in the residual: `POPS_WORK` is perturbed, `POPS_BASE`
is the anchor; `ED`, `T`, module-level `STEQ_T`/`STEQ_ED`/`STEQ_EHB`,
per-ion `SE%STEQ*`, two-photon rates, and the derived-atmosphere state
(`POPTOION` -> `SUP_TO_FULL_V4` -> `COMP_LEV_DIS_BLK`) are all saved
and restored per call by code inspection. Non-finite components are
capped to signed `1D20` and counted (`jfnk_resid_nbad`); broken
coupled sweeps set `JFNK_EVAL_BROKEN` and return sentinel rows.
`JFNK_MASK_FIXED_UNKNOWNS` zeroes fixed `NE`/`T` rows; PTC adds
`(x - x0)*s/dtau` anchored at the baseline.

`JFNK_JVP_V1`: forward difference `[F(x+h v) - F(x)]/h` with
`H_STEP = eps*(XNORM+0.1)/VNORM` subject to a constraint-margin cap
(`H_CAP` from population-clip and `exp` headroom), `JFNK_PROBE_HSCALE`,
`JFNK_HADAPT_MULT`, and up to `JFNK_CLIP_RETRY` h/10 retries when the
perturbed evaluation newly crosses a constraint (`JFNK_N_NEW_CLIP`).
The effective step is exposed via `JFNK_PROBE_H_LAST`; this change adds
`JFNK_PROBE_CAPPED` to record when the cap binds.

## B. Residual repeatability

Instrumentation: `JFNK_OP_PROBE_V1` emits `OPPROBE REP` (4 baseline
re-evals) and `OPPROBE REPP` (perturbed-point re-eval) lines on the
live residual at admission under `[JFNK_HSCAN]`.

Measured (surrogate `se_f`, deterministic): `df_inf = 0` bitwise on
all of mild/hard/shock, baseline and perturbed. This verifies the
measurement path; it does not establish production repeatability.

Static determinism audit of the production residual: the evaluator is
written as a pure function of its declared state -- every mutated
array has a matched save/restore by inspection, and the coupled sweep
has explicit bracket lists for radiation, line, moment/Knop refresh,
Eddington-record (`ACCESS_F`) and file-cursor state. The residual
itself restores `POPS_WORK`, `ED`, `T`, STEQ arrays and re-syncs
atmospheric state from restored populations before returning. No
counter increments inside the residual enter the residual value; clip
counting is a function of `(x, state)` only.

Caveat not yet resolved by measurement: production repeatability in
both modes requires a live model. The battery emits exactly the data
needed (REP/REPP rows per admission); a real run is listed as pending
in section I.

## C. JVP accuracy

New code: `new_main/mod_subs/jfnk_op_probe_v1.f90`, module
`JFNK_OP_PROBE_MOD`, routine `JFNK_OP_PROBE_V1`. It is callback-based
(resid/jvpa/bamv/dotp + nt/nd layout) so the same measurement code
runs (a) on the surrogate via `se_driver seprobe` and (b) on the
production residual at admission via the adapters
`JFNK_OP_JVPA_V1`/`JFNK_OP_BAMV_V1` added to the interface file and the
call site inserted after the deferred `JFNK_FACTOR_PRECOND_V1(.FALSE.)`
in `jfnk_solve_corr_v1.f` (fenced from `jfnk_resid_nbad`,
`JFNK_N_NEW_CLIP`, `JFNK_EVAL_BROKEN`).

Battery per activation: 4 baseline REP evals; 13 direction classes
(`rnd1/2`, `pops`, `ed`, `t`, `k1`, `kmid`, `kworst`, `blk3`, `basis`,
`edt`, `spike`, `xdep`); per direction the production JVP call (with
`heff`, cap/retry/clip bookkeeping), a BA action, a central difference
at `heff`, and a 5-point h-sweep `heff * {1e-2..1e2}` of forward FD vs
the reference; then ADD (unnormalized `u+v`, policy defect `d_pol` vs
shared-h defect `d_op`) and HOM (`J(2u)` vs `2 J(u)`).

Measured on the surrogate (se_jexact as exact reference, eps=1e-6):

- Repeatability: bitwise 0 (surrogate is deterministic).
- Population directions: eFD ~2.5e-7 at heff=1e-6, eFC ~1e-10;
  h-sweep is textbook: ef proportional to h above the noise floor,
  plateau at emin ~1e-8 near h* ~1e-8..1e-7. FD JVP resolves the exact
  surrogate operator to truncation accuracy.
- ED/T-only directions: eFD ~4e-2..6e-2 on mild/hard, hstar at the top
  of the sweep (heff*1e2). Cause: the weighted norm on a tiny-norm
  image (njv ~2-8e4 vs 1e14 for pop dirs) inflates the relative error;
  the absolute error is still truncation-scaled. Flagged as the class
  to watch on production.
- ADD: d_pol ~1.9e-7, d_op ~9.0e-8 at eps=1e-6 (both scale with eps)
  -- additivity defect sits at the FD truncation floor; the norm-scaled
  h policy adds no measurable extra nonlinearity on the surrogate.
- HOM: d_pol = 0 exactly (norm-scaled h makes `J(a v) = a J(v)`
  bitwise); d_op ~1e-7 = truncation.

Production numbers are pending a live admission (section I).

## D. GMRES consistency

From source (jfnk_solver_v1.f90): the GMRES residual norms are weighted
(`dotprd` = `JFNK_DOT_V1` with `JFNK_RSCALE`). `JFNK_TRES_FAC` (default
0=off) measures `true_res = ||b - A x||_w` with one extra matvec at
every done-exit; a claimed convergence is accepted only when
`true_res <= TRES_FAC*tol`, else `tres_rej` restarts the cycle;
non-finite true_res exits `gmres_tres_nonfinite` rather than restarting
on a bad operator. Restart-cap, stagnation-soft and hookstep
trust-radius rescales all re-measure on the returned vector
(`ltres`/`jfnk_gmres_true_res(_ok)`), and the caller gate in
`newton_hook_v1` accepts the step only when
`true_res < 0.9*ref_norm`; `true_res > |F|` skips a useless re-solve.
The Hessenberg LS residual is printed against `true_res` at every
restart (`GMRES RESTART DIAG`) -- audit notes recorded hess_res
15x-450x below true_res on exits, which is why the explicit
measurement gates acceptance.

Coverage assessment: ordinary convergence, restart, restart-cap,
stagnation-soft, and trust-radius-rescale exits all produce an explicit
true-residual measurement on the returned vector. The forcing term is
computed in the preconditioned (weighted) norm, consistent with the
left-preconditioned formulation used. Remaining gap for the production
run to confirm: that `TRES_FAC` is actually armed in the runs being
audited (default off), and that `jvp_soft` exits record `ltres` before
returning (they do, lines ~2896-2918).

## E. BA comparison

`JFNK_OP_BAMV_V1` forms `BA_A*v` (per-depth NT x NT blocks) plus the
stored `CMAT_RAD_BANDS` off-diagonals when `JFNK_BA_RAD_BANDS_READY`,
plus the live PTC diagonal `s/dtau`, plus the fixed-row mask. This is
the JV_ANAL quantity of the JVP diagnostic extended with PTC, so it is
the derivative content of the preconditioner as GMRES sees it.

By construction the BA action is an approximation, not an exact
Jacobian: it omits (i) residual nonlinearity (it is the frozen linear
assembly), (ii) cross-depth coupling beyond the stored bands
(advection/J-coupling terms not in `BA_A`), (iii) under coupled-J,
essentially the whole radiation-response derivative -- the sweep's
recomputed J fields have no representation in BA_A. Under frozen-J the
FD-vs-BA comparison isolates exactly terms (i)+(ii); the JVP
diagnostic already prints `|FD - BA*v|` per direction and this
battery's `eFD` column is the same quantity over 13 direction classes.

Surrogate cross-check: with `se_jexact` in the bamv slot the FD JVP
matches to truncation (eFD ~2.5e-7); the same OPPROBE records on a
production admission give the frozen-J BA error directly. `JFNK_BA_JVP`
remains a diagnostic, not promoted.

## F. LU investigation (T21)

Reproduced the census directly on the stored surrogate blocks
(J0 = A diag(n0) + C, the se_factor_lu replica: ridge ladder
{0,1e-12,1e-8,1e-4,1e-2}, reg rows -> identity, DGEEQU, DGETRF).

Fortran MZDK measurements (backward error ||Mz-f||/||f|| per depth):

- mild: k1 4.0e7, k2 4.7e-2, k3 3.9e8, k4 8.1e-7
- hard: k1 2.8e10, k2 1.0e4, k3 1.1e8, k4 6.9e11
- shock: k1 4.6e13, k2 9.1e-7, k3 1.4e-1, k4 1.3e2

Block analysis (numpy replica of the identical pipeline):

- Effective-M condition numbers 1e15..1e27 per depth (well past
  1/eps_mach ~ 4.5e15), so no double-precision LU can reach a backward
  error below ~cond*eps: the measured 1e7..1e11 residuals sit inside
  that bound -- the failures are structural, not factorization bugs.
- Null-space attribution: the smallest right-singular vectors are
  concentrated on population columns 30-32 (mild/hard) and 25-33
  (shock); the corresponding columns of J0 have norm ~1e-8 against
  block max 1e7..1e13 (relative 1e-16..1e-22). These columns belong to
  population unknowns sitting at the n0 = 1e-30 floor: their Jacobian
  columns carry only the ~1e-8 conservation/decay coupling, so the
  blocks have an approximate 3-dimensional null space. Column
  equilibration rescales the column but cannot recover information the
  physics does not contain.
- Zero rows are handled on the row side by reg->identity (nreg up to
  19/block on mild); the residual null space is on the column
  (unknown) side, which no row fix touches.
- Left-null vectors concentrate on population/conservation rows (34,
  33, 30-31) -- near-dependent SE equations, consistent with a
  conservation-constraint duplicate, which is what puts the f_rhs
  component outside the effective range.

Conclusion: T21 XFAILs are retained and justified. The defect is
variable-space scaling (population-floor unknowns decouple at ~30
decades below block scale), not incorrect assembly, not pivot
strategy, not insufficient ridge. Ridge escalation cannot fix a column
whose coupling is absent by construction; candidate remedies are
architectural (mask floor-pinned unknowns, or solve in a rescaled
unknown basis) and are out of scope here.

## G. Recommended implementation

Justified by measurements so far:

1. Keep the probe battery as the production measurement channel
   (done, under JFNK_HSCAN, opt-in). Do not add VADAT keys.
2. On the first production admissions, read `OPPROBE` rows: REP/REPP
   for state contamination, per-tag eFD/eFC/emin/hstar for the h
   operating point, caps/retr/clip for the constraint regime, ADD d_pol
   vs d_op for h-policy nonlinearity, before touching H_STEP policy.
3. T21: keep XFAIL; the fix class is unknown-space rescaling or
   floor-unknown masking, to be designed separately -- not a ridge or
   LSQ tweak.
4. If production ADD shows d_pol >> d_op (policy-driven nonlinearity),
   only then evaluate the cycle-consistent-h alternative; the surrogate
   shows the current norm-scaled rule is exact under homogeneity and
   truncation-limited under additivity, so there is no evidence yet for
   changing the default h policy.

No changes to defaults, physical equations, floors/clipping, GMRES
method, or preconditioner were made or are recommended now.

## H. Test report

- Compiler: Intel Fortran 2024.0.2 (`/mpcdf/.../intel/2024.0.2/bin/ifort`).
- Build of the touched production files, ifort,
  `-qopenmp -mkl -nbs -extend-source -traceback -O3 -I mod -module mod`:
  `jfnk_cmfgen_int_mod.f`, `jfnk_op_probe_v1.f90`,
  `jfnk_cmfgen_interface_v1.f`, `jfnk_solve_corr_v1.f` all compile clean
  (one pre-existing ES-descriptor remark). A full `./build.sh` was not
  re-run; only the four objects were rebuilt against the existing
  `mod/` interface set.
- Source revision: `22d63df9` + this working-tree diff.
- Commands run:
  - `bash unit_tests/build.sh` (rebuild se_driver incl. probe module)
  - `./bin/se_driver seprobe systems/{mild,hard,shock} 1e-6`
  - `python3 unit_tests/python/test_jfnk_opint.py` -> 36/36 checks pass
  - `python3 unit_tests/run_all_tests.py` -> 21/21 test files pass
    (incl. new test_jfnk_opint.py; T21 XFAILs retained)
  - `sepc` MZDK census + numpy block replica for section F.

Changed files: `new_main/mod_subs/jfnk_op_probe_v1.f90` (new),
`jfnk_cmfgen_interface_v1.f` (cap flag set + 2 adapters),
`jfnk_cmfgen_int_mod.f` (JFNK_PROBE_CAPPED), `jfnk_solve_corr_v1.f`
(HSCAN call site + fence), `control_variable_mod.f` (doc only),
`mod_subs/Makefile` (object + rule), `unit_tests/build.sh`,
`unit_tests/se_driver.f90` (seprobe + 2 adapters),
`unit_tests/python/test_jfnk_opint.py` (new), `unit_tests/README.md`
(inventory line). Not committed per workflow; the tree remains at the
checked-out HEAD plus this diff.

## I. Remaining limitations

1. No complete CMFGEN model has been run for this investigation. Every
   production-mode number (frozen-J and coupled-J REP/DIR/HS/ADD/HOM,
   real cap/retry/clip incidence, true-residual ratios per GMRES exit)
   requires one live admission with `JFNK_HSCAN=T`; the instrumentation
   is in place and verified on the surrogate, the measurements are
   pending.
2. Repeatability on production is asserted by inspection (save/restore
   audit), not by measurement; the REP/REPP rows will confirm or
   localize contamination on the first HSCAN admission.
3. The BA-vs-FD comparison on production (section E) is wired but
   unmeasured; the coupled-J BA limitation is established by
   construction (no radiation response in BA_A), not yet quantified.
4. The T21 null-column analysis is on the surrogate systems; the same
   census on production `JFNK_BA_A` blocks is the natural next probe.
5. The probe battery costs ~100 residual evals per admission (~a full
   radiation sweep each under coupled mode) -- deliberately bounded and
   opt-in, but not free.

## J. Production readiness

Current state, stated exactly:

- Component level: the production JFNK operator is instrumented for
  direct measurement at the admission baseline; surrogate-level checks
  pass (repeatability, FD-vs-exact truncation behavior, additivity,
  homogeneity). T21 LU defects are explained structurally (population-
  floor null columns, cond ~1e20-1e27), not masked by an XFAIL label.
- Validated on surrogate systems: yes, where the reference is exact.
- Validated on a complete CMFGEN model: not yet -- pending one
  instrumented admission per mode.
- Not conflating categories: surrogate verification of the measurement
  code path is not production verification of the operator.

Proposed next step for review: run one existing KNOP (coupled) and one
SOB (frozen) admission each with `JFNK_HSCAN=T` + `JFNK_TRES_FAC>0`,
collect the OPPROBE/GMRES-TRES rows, then decide whether any
correctness fix is warranted before considering h-policy or solver
changes.
