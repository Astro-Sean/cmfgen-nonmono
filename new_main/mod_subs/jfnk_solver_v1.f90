!----------------------------------------------------------------------
! JFNK_SOLVER_V1: Jacobian-Free Newton-Krylov solver with GMRES(m)
! and Hookstep/trust-region globalization.
!
! Based on the methodology of:
!   Arramy, de la Cruz Rodriguez & Leenaarts (2024), A&A 690, A12
!   https://doi.org/10.1051/0004-6361/202449963
! and the reference implementation:
!   apwillis1/JFNK-Hookstep (https://github.com/apwillis1/JFNK-Hookstep)
!
! This module provides:
!   - GMRES(m) with LEFT preconditioning (default) and restart
!   - Hookstep (dogleg trust-region) globalization
!   - Finite-difference Jacobian-vector products (forward difference)
!   - Newton iteration with convergence diagnostics
!   - Eisenstat-Walker forcing terms for adaptive GMRES tolerance
!
! The solver is Jacobian-free: it never constructs the full Jacobian
! matrix.  The Jacobian-vector product J*v is approximated by:
!   J*v ≈ [F(x + eps*v) - F(x)] / eps
! where F is the residual function and eps is a finite-difference step.
!
! NOTE (review §5): the JVP captures the exact directional derivative
! of the implemented frozen-J residual, including nonlinear collisional,
! electron-density and temperature dependencies where those dependencies
! are explicitly evaluated by the residual.  It does NOT capture
! dF_rad/dT or other radiation-dependent nonlinearities because those
! have been replaced by the frozen linear CMAT_RAD representation.
! For fixed J, fixed T, this is the correct Jacobian.  For free T,
! the Jacobian is approximate (see review §4).
!
! The preconditioner M is applied as z = M^{-1} * v.
! Both LEFT and RIGHT preconditioning are supported; the default is
! LEFT (jfnk_use_left_pc = .true.), following Arramy et al. (2024).
!   LEFT:  Krylov subspace K(M^{-1}A, M^{-1}r), initial residual IS
!          preconditioned.  The GMRES stopping tolerance must be
!          measured in ||M^{-1}r|| units (see EW forcing-term fix in
!          newton_hook_v1).
!   RIGHT: Krylov subspace K(A*M^{-1}, r), initial residual is NOT
!          preconditioned.  The preconditioner is applied to each
!          Arnoldi basis vector before the matvec, and to the Krylov
!          solution before adding to x.
! For CMFGEN non-monotonic V4, M is the fixed-J BA matrix
! (equilibrated depth-block LU with local adaptive ridge).
! Fix #9 (review 2026-09-11): the default is LEFT preconditioning
! (jfnk_use_left_pc = .true.).  The production V1 configuration is:
!   GMRES:
!     preconditioning = LEFT
!     inner product    = EUCLIDEAN (RSCALE-weighted via JFNK_DOT_V1)
!     operator         = J (frozen, finite-difference JVP)
!     preconditioner   = M^{-1} (equilibrated depth-block BA LU)
!
! editted by Sean Brennan (V4, 2026-09)
!----------------------------------------------------------------------
module jfnk_solver_v1_mod
   implicit none
   save
   !
   ! Module-level state (mirrors the `newton` module in the reference)
   !
   double precision, allocatable :: jfnk_x(:)       ! current solution vector
   double precision, allocatable :: jfnk_fx(:)      ! current residual F(x)
   double precision              :: jfnk_tol        ! current |F(x)|
   ! Estimated max relative population correction still required,
   ! max_i |F_i|/|BA_ii| over population rows, refreshed by each
   ! JFNK_RESIDUAL_V1 call.  JFNK_MAXCH_TOL>0 adds it as a second
   ! convergence test so a residual norm dominated by already-
   ! converged rows cannot stall an otherwise-finished solve.
   double precision              :: jfnk_maxch_est = -1d0
   double precision              :: jfnk_maxch_tol = 0d0
   double precision              :: jfnk_del        ! trust-region radius
   integer                       :: jfnk_nits       ! Newton iterations done
   integer                       :: jfnk_gits       ! GMRES iterations (last)
   character(len=96)             :: jfnk_exit_reason = 'not set'
   !
   ! Eisenstat-Walker forcing term state (added 2026-09-07):
   ! The forcing term eta_k controls how accurately GMRES solves the
   ! inner linear system J*s = -F at each Newton step.  Instead of a
   ! fixed tolerance, eta_k adapts based on the Newton convergence rate:
   !
   !   eta_k = gamma * (|F_k| / |F_{k-1}|)^alpha
   !
   ! with safeguards:
   !   eta_k = max(eta_k, 0.5 * |F_k| / |F_{k-1}|)   (prevent over-solving)
   !   eta_k = min(eta_k, eta_max)                    (cap at user max)
   !   eta_k = max(eta_k, eta_min)                    (floor for stability)
   !
   ! When |F| is decreasing rapidly (far from root), eta_k is large and
   ! GMRES stops early (inexact Newton is fine).  When |F| is decreasing
   ! slowly (near root), eta_k is small and GMRES solves more accurately.
   !
   ! References:
   !   Eisenstat & Walker (1994), SIAM J. Optim. 4, 393
   !   Kelley (2003), "Solving Nonlinear Equations with Newton's Method"
   !
   ! Note: the module already has a blanket SAVE statement (line above),
   ! so module variables are implicitly saved.  The per-variable SAVE
   ! attributes that were here are redundant and rejected by
   ! gfortran -std=f2008; removed per ChatGPT audit fix #3 (2026-09-09).
   double precision :: jfnk_tol_prev = 0d0    ! |F| at previous Newton step
   double precision :: jfnk_eta = 0d0          ! current Eisenstat-Walker eta
   ! jfnk_ew_adapt (VADAT [JFNK_FADAPT], default T): T keeps the
   ! adaptive EW Choice-2 forcing term below; F restores the fixed
   ! user gtol each Newton step -- the fixed path stays one branch
   ! away for A/B comparison.  (Devin, 17-Oct-2026)
   logical :: jfnk_ew_adapt = .true.
   ! eta_prev stores the previous Newton step's forcing term: the
   ! Choice-2 safeguard eta_k >= gamma*eta_{k-1}^alpha needs it.
   ! (Oct-2026 EW fix)
   double precision :: jfnk_eta_prev = 0d0     ! eta at previous Newton step
   double precision :: jfnk_eta_gamma = 0.9d0  ! EW gamma parameter
   ! EW Choice 2 requires alpha in (1,2] for superlinear local
   ! convergence (Eisenstat & Walker 1994, Thm 3.3); the previous
   ! alpha=0.5 contracted too slowly for the theory to apply.  The
   ! golden ratio is the reference value used by EW and Kelley.
   double precision :: jfnk_eta_alpha = 1.618033988749895d0
   ! Floor matches the bottom of the range Arramy et al. (2024,
   ! A&A 690, A12) found useful for multilevel NLTE.  Below this,
   ! GMRES spends formal solutions without improving the Newton step.
   ! eta_max=0.9 is the standard EW cap; 0.1 forced GMRES to solve
   ! every Newton step far more accurately than inexact-Newton
   ! theory requires.
   double precision :: jfnk_eta_min = 3.0d-4   ! minimum forcing term
   double precision :: jfnk_eta_max = 0.9d0    ! maximum forcing term
   !
   ! Diagnostics counters
   !
   integer :: jfnk_total_f_evals = 0   ! total residual evaluations
   integer :: jfnk_total_jvp_evals = 0 ! total JVP evaluations
   ! jfnk_fx_fresh (03-Oct-2026): set by the driver/baseline rebuild
   ! when jfnk_fx already equals F(jfnk_x) -- lets newton_hook skip
   ! the entry residual re-eval (a full formal sweep in coupled mode).
   logical :: jfnk_fx_fresh = .false.
   integer :: jfnk_precond_calls = 0   ! total preconditioner applications
   integer :: jfnk_precond_fails = 0   ! preconditioner failures
   integer :: jfnk_nan_rejects = 0     ! NaN/Inf rejected steps
   !
   ! Stabilization usage counters (review §6, 2026-09-10):
   ! Track how often each artificial safeguard fires.
   !
   integer :: jfnk_n_step_scaled = 0      ! |s|_inf > 10 scaling triggered
   integer :: jfnk_n_zmax_hit = 0         ! preconditioner |z| > ZMAX clamped
   integer :: jfnk_n_reg_proj = 0         ! regularized null rows projected to 0
   integer :: jfnk_n_cauchy_used = 0      ! Cauchy fallback accepted
   integer :: jfnk_n_ls_used = 0          ! line search accepted
   integer :: jfnk_n_tr_rejected = 0       ! trust-region rejected step
   integer :: jfnk_n_descent_retry = 0     ! non-descent direction retry
   !
   ! Maximum GMRES iterations per Newton step (added 2026-09-07).
   ! When > 0, overrides the default limit of m (one restart cycle).
   ! When 0, the default is used (jfnk_gits = m, or 9999 if del == 0).
   ! Set from JFNK_MAX_GMRES control variable via solve_for_pops_jfnk_v1.f.
   !
   integer :: jfnk_max_gmres_its = 0
   !
   ! Maximum GMRES restarts per solve (added 22-Oct-2026).  Guards
   ! against restart storms on an invariant Krylov space; default 1
   ! preserves the historical cap.  The total iteration count is still
   ! bounded by imx regardless of this cap.  Set from the JFNK_GMRES_RST
   ! control variable via JFNK_SOLVE_CORR_V1.
   !
   integer :: jfnk_gmres_max_restart = 1
   !
   ! True-residual acceptance factor (added 22-Oct-2026, numerical-
   ! methods review Phase 1).  When > 0, a claimed LSQ convergence
   ! (res <= tol) in gmres_m_v1 must survive an explicit
   ! ||b - A x||_w measurement: the Hessenberg residual is only a
   ! model of the true residual, and with a finite-difference matvec
   ! it has been observed 150-450x low at restart (the Arnoldi basis
   ! does not compose a fixed linear operator).  Set from the
   ! JFNK_TRES_FAC control variable via JFNK_SOLVE_CORR_V1;
   ! 0 = off preserves the historical accept-on-hess_res path.
   !
   double precision :: jfnk_gmres_tres_fac = 0d0
   !
   ! Left vs right GMRES preconditioning (Arramy 2024 prefers left).
   ! Set from JFNK_LEFT_PRECOND before NEWTON_HOOK_V1.
   ! editted by Sean Brennan (08-Sep-2026)
   !
   logical :: jfnk_use_left_pc = .true.
   !
   ! Preconditioner refresh bookkeeping (Issue 3).
   !
   integer :: jfnk_last_gmres_its = 0
   ! 2026-09-13 TR-units fix: norm of the step returned by the last
   ! gmres_m_v1 call in the units the hookstep constrains.  Under left
   ! PC that is ||y|| = ||s||_w (V orthonormal in dotprd); under right
   ! PC it is sqrt(y' G y) = ||s||_w with G = Z' W Z, so the
   ! trust-region radius is in physical-step units either way.
   double precision :: jfnk_last_ynrm = 0d0
   !
   ! editted by Devin, 12-Nov-2026: true preconditioned linear residual
   ! ||M^{-1}(b - A x)||_w measured inside gmres_m_v1, exported so the
   ! newton_hook_v1 stagnation gate can test the real residual rather
   ! than the Hessenberg least-squares model residual.  With a finite-
   ! difference matvec the Krylov basis need not span a fixed linear
   ! operator; OUTGEN showed hess_res 15x-450x below true_res on exits
   ! this gate used to accept.  _ok marks whether the last gmres_m_v1
   ! exit carried a measurement.
   double precision :: jfnk_gmres_true_res = 0d0
   logical :: jfnk_gmres_true_res_ok = .false.
   integer :: jfnk_need_precond_refresh = 0
   integer :: jfnk_gmres_m_dim = 10
   double precision :: jfnk_pc_ref_frac = 0.8d0
   !
   ! Ridge floor: the relaxation path in jfnk_cmfgen_interface_v1.f will
   ! not relax JFNK_BA_RIDGE_FAC below this value.  Captured at first
   ! factorization so a user-supplied ridge (e.g. 0.1) is preserved
   ! instead of being relaxed back to the historical 1e-2 floor.
   ! Added 2026-09-10.
   !
   double precision :: jfnk_ridge_floor = 1.0d-2
   logical :: jfnk_ridge_floor_set = .false.
   !
   ! Preconditioner audit flags (jfnk-precond-fixes branch, Oct-2026).
   ! All default to the historical behaviour; each maps to a VADAT key
   ! that is assigned by the interface layer (JFNK_INIT_V1 /
   ! JFNK_FACTOR_PRECOND_V1) so the extracted block routines below stay
   ! free of CONTROL_VARIABLE_MOD dependencies and unit-testable.
   !
   ! jfnk_pc_signed: [JFNK_PTC_SIGNED] fold/ridge/rescale use
   !    sign(d_i)*sigma_i so |d_i+shift| = |d_i|+sigma_i on
   !    loss-dominated (d<0) rows instead of moving the diagonal
   !    toward zero.
   ! jfnk_pc_reg_zerocol: [JFNK_REG_ZERO_COL] regularized rows also
   !    get their column zeroed (the frozen T/Ne path already does
   !    this); otherwise a substituted row still couples to other
   !    equations through its column entries.
   ! jfnk_pc_tiny_equil: [JFNK_ROW_TINY_EQUIL] classify near-null rows
   !    on the column-equilibrated block instead of the raw block, so
   !    a row in different units is not mistaken for a null row.
   ! jfnk_rcond_gate/rcond_min/itry_cap: [JFNK_PC_RCOND_GATE,
   !    JFNK_PC_RCOND_MIN] replace the x=1 forward-error micro-test by
   !    an explicit reciprocal-condition threshold and cap the ridge
   !    ladder attempts per depth.
   !
   logical :: jfnk_pc_signed = .false.
   logical :: jfnk_pc_reg_zerocol = .false.
   logical :: jfnk_pc_tiny_equil = .false.
   logical :: jfnk_rcond_gate = .false.
   double precision :: jfnk_rcond_min = 1.0d-13
   integer :: jfnk_itry_cap = 3
   !
   ! Persistent workspace for newton_hook_v1 and gmres_m_v1 (audit §3,
   ! 09-Sep-2026): replaces large automatic arrays that were on the stack.
   ! For N=40000 (NT=200, ND=200), these accounted for ~2.5MB of stack.
   ! Allocated on first call and resized if N changes.
   !
   double precision, allocatable :: jfnk_w_x(:)       ! x_ in Newton
   double precision, allocatable :: jfnk_w_fx(:)      ! fx_ in Newton
   double precision, allocatable :: jfnk_w_x2(:)      ! x__ in Newton
   double precision, allocatable :: jfnk_w_fx2(:)     ! fx__ in Newton
   double precision, allocatable :: jfnk_w_s(:)       ! s in Newton
   double precision, allocatable :: jfnk_w_fxpc(:)    ! fx_pc in Newton
   double precision, allocatable :: jfnk_w_gmres_w(:) ! w in GMRES
   double precision, allocatable :: jfnk_w_gmres_z(:) ! z in GMRES
   integer :: jfnk_w_alloc_n = 0                       ! current allocation size
   !
   ! Line-search globalization control (added 2026-09-10, review §3):
   ! When the full Newton step (α=1) increases the residual, backtrack
   ! with α = β, β^2, β^3, ... until the Armijo sufficient-decrease
   ! condition is met.  This avoids re-solving GMRES when the Newton
   ! direction is good but the step is too large.
   !
   ! References:
   !   Nocedal & Wright (2006), Numerical Optimization, ch.3.1–3.2
   !   Eisenstat & Walker (1994), SIAM J. Optim. 4, 393
   !
   logical          :: jfnk_do_line_search = .true.   ! enable line search
   integer          :: jfnk_ls_max_steps = 20          ! max backtracking steps
   double precision :: jfnk_ls_beta_fac = 0.5d0      ! backtracking factor
   double precision :: jfnk_ls_armijo_c = 1.0d-4      ! Armijo parameter
   !
   ! Non-monotone line search (Li & Fukushima 2000, Appl. Math.
   ! Comput.; Grippo-Lampariello-Lucidi reference norm): when the
   ! coupled residual carries evaluation noise, a monotone Armijo test
   ! rejects steps whose true decrease is masked by sweep noise.  The
   ! reference norm is the max of the last jfnk_ls_nmwin_s ACCEPTED
   ! |F| values (ring buffer below), so a bounded temporary increase
   ! still satisfies the descent condition.
   !
   logical          :: jfnk_ls_nonmono = .true.
   integer          :: jfnk_ls_nmwin_s = 4
   double precision :: jfnk_ls_nmhist(8) = -1d0
   integer          :: jfnk_ls_nmidx = 0
   !
   ! Cauchy (preconditioned-residual) fallback direction (review §3):
   ! When the Newton line search fails, try d_C = M^{-1}*F as a
   ! fallback direction with its own simple sufficient-decrease line
   ! search.  This gives "Newton when Newton works; stabilized
   ! residual correction when it doesn't."
   !
   logical          :: jfnk_do_cauchy = .true.  ! enable Cauchy fallback
   !
   ! Step-safeguard caps (parameterized 2026-10 for unit testing):
   ! jfnk_step_inf_cap is the maximum allowed |s|_inf per Newton step
   ! (= MAX_JFNK_COR semantics: factor-10 fractional population change);
   ! jfnk_runaway_fac is the catastrophic weighted-step norm, measured
   ! relative to ref_norm = ||M^{-1} F||_w, above which the step is
   ! rejected and the trust radius shrunk.  The reference must be
   ! relative: under a non-Euclidean JFNK_RSCALE the absolute norm
   ! scale is arbitrary (|F|_w ~ 1e10 under DIAG), so an absolute cap
   ! either always fires or never does.  A Newton step satisfies
   ! |s|_w ~ |M^{-1} F|_w; 100x that is genuine singular-inverse
   ! garbage (the 1e31 class this check was added for).
   !
   double precision :: jfnk_step_inf_cap = 10d0
   double precision :: jfnk_runaway_fac  = 1.0d2
   !
   ! Logarithmic unknown and pseudo-transient continuation.
   ! Set from VADAT [JFNK_LOGX] and [JFNK_PTC] in RD_CONTROL_VARIABLES
   ! (both default T).  See cmfgen_obsidian/solvers/JFNK_SOLVER.md.
   !
   ! jfnk_log_unknown: x = ln(n/n0).  The step cap is then ln(MAX_JFNK),
   ! written into jfnk_step_inf_cap by the driver.
   ! jfnk_ptc_active: G(x) = F(x) + (x - x0)*s/dtau, s_i =
   ! jfnk_ptc_scale(i).  dtau grows when |F| falls and shrinks when it
   ! does not.
   !
   ! edited by Devin, 06-Oct-2026 (SER anchor fix): the pseudo-time term
   ! is now anchored at jfnk_ptc_x0 = x_k, the current Newton iterate --
   ! switched evolution relaxation, Kelley & Keyes (1998, SIAM J. Numer.
   ! Anal. 35, 508) and Knoll & Keyes (2004, JCP 193, 357 Sec 4.2).  The
   ! previous form anchored at the solve baseline x0 = 0, so the term
   ! pulled later iterates back toward x = 0 and a converged |G| could
   ! hide F(x) = -(x)*s/dtau /= 0.  newton_hook_v1 re-anchors x0 at each
   ! accepted iterate and strips the term from every residual it reads,
   ! so jfnk_fx/jfnk_tol and all Newton logic see F and |F|.
   !
   logical :: jfnk_log_unknown = .true.
   logical :: jfnk_ptc_active = .true.
   double precision :: jfnk_ptc_dtau = 1d0
   double precision :: jfnk_ptc_dtau_max = 1d8
   !
   ! FIX_T_AUTO adaptive damped-T step guard (generic: the driver names
   ! a variable block [dst,dend] of the row index jt inside a repeating
   ! record of length nt).  When on, every Newton direction leaving
   ! GMRES and every Cauchy direction has its T components scaled by
   ! jfnk_tauto_alpha and capped at +-jfnk_tauto_dmax.  The damping acts
   ! on the step only; the residual, JVP and preconditioner still see
   ! the live T columns, so the Krylov operator keeps the true
   ! temperature subspace (directive: never zero the T columns when T
   ! is merely damped).
   !
   logical          :: jfnk_tauto_on   = .false.
   double precision :: jfnk_tauto_alpha = 1d0
   double precision :: jfnk_tauto_dmax  = 0.5d0
   integer          :: jfnk_tauto_nt    = 0
   integer          :: jfnk_tauto_jt    = 0
   integer          :: jfnk_tauto_dst   = 0
   integer          :: jfnk_tauto_dend  = 0
   double precision, allocatable :: jfnk_ptc_scale(:)
   logical :: jfnk_ptc_scale_ready = .false.
   ! Anchor x_k for the pseudo-time term (above).  JFNK_SET_PTC_SCALE_V1
   ! allocates it at 0 (the baseline x = 0 is the first anchor);
   ! newton_hook_v1 re-anchors at each accepted iterate via
   ! jfnk_ptc_set_anchor_v1.
   double precision, allocatable :: jfnk_ptc_x0(:)
   logical :: jfnk_ptc_anchor_ready = .false.
   !
   ! Set by JFNK_RESIDUAL_V1 on every call: number of residual
   ! components that were NaN or |F|>1D30 before capping.  The F(0)
   ! identity test in JFNK_SOLVE_CORR_V1 uses it to tell a broken
   ! eval (non-finite entries) from an honest lagged-field offset,
   ! which the uniform |F|>magnitude comparison cannot do.
   !
   integer :: jfnk_resid_nbad = 0
   !
   ! Explicit Krylov-failure signalling (audit 17-Sep-2026).
   ! jfnk_jvp_failed    : set by JFNK_JVP_V1 when the finite-difference
   !                      JVP produced too many NaN/Inf entries to trust.
   !                      Replaces the old "JV=0" convention, which GMRES
   !                      misread as a genuine zero-action direction
   !                      (artificial happy breakdown / stagnation).
   ! jfnk_pc_clip_event : set by JFNK_PRECOND_V1 when the |z|>1D8 hard
   !                      clip or a NaN->0 write fires.  Clipping makes
   !                      M^{-1} nonlinear, which invalidates the GMRES
   !                      Krylov construction -- the solve is aborted
   !                      (info=3) and a preconditioner refresh requested
   !                      instead of silently solving a different problem.
   ! Both are reset at the top of each Newton iteration and checked
   ! after every matvec/psolve call inside GMRES.
   !
   logical :: jfnk_jvp_failed = .false.
   logical :: jfnk_pc_clip_event = .false.
   !
   ! GMRES hookstep-recalculation flag (edited by Devin, 05-Oct-2026,
   ! recompute-pass descent-skip fix): set by gmres_m_v1 when it returns
   ! from its info==2 path, where hook_step_v1 recalculated the Krylov
   ! coefficients for a new trust radius on the unchanged Arnoldi basis.
   ! newton_hook_v1 tests this flag -- not ginfo -- to skip the descent
   ! check on recompute passes, because gmres_m_v1 resets info=0 before
   ! returning so ginfo==2 never survives to the check.  Cleared at the
   ! top of each Newton iteration, like the failure flags above.
   !
   logical :: jfnk_gmres_recomputed = .false.
   !
   ! Last completed Arnoldi cycle of gmres_m_v1 (basis dimension jj and
   ! residual norm beta).  The info==2 hookstep recompute combines these
   ! with the caller's h/v arrays, so they must describe that caller's
   ! basis.  They live at module scope -- not as gmres_m_v1 SAVE
   ! variables -- because diagnostic callers
   ! (jfnk_gmres_dim_experiment_v1 runs gmres_m_v1 with m_test up to
   ! 40 > production m) would otherwise leave a stale pair behind: a
   ! recompute pass would then index y(1:jj), h(1:jj+1,1:jj) and
   ! v(:,1:jj) past the caller's m-sized arrays.  Any external
   ! gmres_m_v1 caller must save and restore this pair.
   ! (numerical audit fix, 2026-10)
   !
   integer          :: jfnk_gmres_jj    = 0
   double precision :: jfnk_gmres_beta  = 0d0
   ! The opt-in [USE_JFNK_SOLVER] path needs this checkpoint to replace
   ! only the last restart cycle during an info=2 hookstep recompute.
   double precision, allocatable :: jfnk_gmres_x_cycle_start(:)
   !
   ! Right-preconditioning step map (edited by Devin, 2026-10,
   ! right-PC TR-units fix): under right PC the GMRES iterate y maps
   ! to the physical step via s = Z y with Z = M^{-1} V, so
   ! ||s||_w^2 = y' G y with G = Z' W Z (dotprd inner product).
   ! hook_step_g_v1 constrains this G-norm instead of ||y||, keeping
   ! the trust radius in physical-step units.  Stored per completed
   ! Arnoldi cycle like the jj/beta pair above -- any external
   ! gmres_m_v1 caller must save/restore them the same way.
   !
   double precision, allocatable :: jfnk_gmres_z(:,:)   ! M^{-1} V basis
   double precision, allocatable :: jfnk_gmres_g(:,:)   ! Z' W Z Gram
   !
   ! Cycle-start residual of the previous restart cycle (23-Oct-2026,
   ! LGMRES diagnosis): Baker-Jessup-Manteuffel stagnation shows up as
   ! a repeating/alternating residual direction across cycles, so the
   ! restart diag reports |<r_c, r_{c-1}>|/(|r_c||r_{c-1}|).  Purely
   ! diagnostic -- never enters the solve.
   !
   double precision, allocatable :: jfnk_gmres_rprev(:)
   !
   ! LGMRES augmentation state (23-Oct-2026; Baker, Jessup &
   ! Manteuffel 2005, "A Technique for Accelerating the Convergence of
   ! Restarted GMRES").  The coupled-SOB solve shows exactly the
   ! cyclic restart stagnation that paper targets: the cycle-1
   ! residual drops ~12% while cycles 2+ make <1% each and the
   ! restart-boundary residual direction returns cos(r_c,r_{c-1}) ~
   ! 0.99 -- every restart rebuilds the same Krylov space because the
   ! converged directions were discarded with the basis.  LGMRES keeps
   ! the last K per-cycle solution corrections z_i = x_e - x_s and
   ! solves each restart's least-squares over V_m + span(z_i) jointly.
   ! The image A z_i = r_s - r_e is a residual difference already in
   ! hand at cycle start, so augmentation costs zero extra matvecs.
   !
   ! Right-preconditioning only: the physical step s = Z y + D c
   ! appends d_i in step space, equivalent to t-space vectors M d_i --
   ! no forward M action is needed.  With JFNK_LEFTPC=T the flag is
   ! ignored and a note prints once per solve.
   !
   ! [JFNK_LGMRES] default F: plain restarted GMRES unchanged.
   ! [JFNK_LGMRES_K] default 2 (clamped 1..8): corrections retained.
   ! Vectors are kept only within a single linear solve -- no
   ! recycling across Newton steps.
   !
   logical :: jfnk_lgmres_on = .false.
   integer :: jfnk_lgmres_kmax = 2
   integer :: jfnk_lgmres_n = 0
   double precision, allocatable :: jfnk_lgmres_d(:,:)    ! z_i
   double precision, allocatable :: jfnk_lgmres_ad(:,:)   ! op z_i (A z or M^{-1}A z)
   double precision, allocatable :: jfnk_lgmres_dlast(:)  ! completed-cycle step
   double precision, allocatable :: jfnk_lgmres_q(:,:)    ! orth image complements
   !
   ! [JFNK_GDEFL] (default F): GMRES-DR deflated restarting (Morgan
   ! 2002, SISC 24(1) 20-37; Parks et al. 2006, SISC 28(5) 1651-1674
   ! for the GCRO-DR sequence version).  At each restart the smallest-
   ! |theta| harmonic Ritz pairs of the completed cycle's Hessenberg
   ! block are extracted and their Arnoldi images written into the
   ! same augmentation slots the LGMRES retention uses, so the next
   ! cycle's least-squares keeps the eigenmodes restarted GMRES
   ! forgets -- the restart-angle ~0.99 stall class.  A complex pair
   ! keeps both real vectors, so the live count can exceed
   ! [JFNK_GDEFL_K] (default 4, clamped 1..8) by one.  jfnk_aug_kmax
   ! is the combined augmentation capacity (LGMRES + GDEFL slots),
   ! hard bound 8 matching the lg_kmax workspace bound in gmres_m_v1.
   ! [JFNK_RECYC] (default F) additionally carries the retained set
   ! across gmres_m_v1 calls within one activation -- the operator
   ! has moved one Newton step, but augmentation only needs
   ! approximately-correct images.
   logical :: jfnk_gdefl_on = .false.
   integer :: jfnk_gdefl_kmax = 4
   integer :: jfnk_aug_kmax = 0
   logical :: jfnk_recyc_on = .false.
   !
   ! [JFNK_CYCSTALL] (default F): a completed restart cycle that trims
   ! the LS residual by less than [JFNK_CYCSTALL_FAC] (default 0.02)
   ! is rebuilding the same invariant space; [JFNK_CYCSTALL_N]
   ! (default 2) consecutive such cycles exit info=1 and hand the
   ! partial direction to the Newton inexact-direction path instead
   ! of grinding the iteration cap out.
   logical :: jfnk_cycstall_on = .false.
   double precision :: jfnk_cycstall_frac = 2.0d-2
   integer :: jfnk_cycstall_nmax = 2
   !
   ! GMRES stagnation criterion (audit 17-Sep-2026): an iteration counts
   ! as non-improving unless it reduces the residual by more than this
   ! relative factor.  The previous 1D-14 threshold only fired on
   ! machine-level stalls; with finite-difference JVP noise, gains
   ! below ~1D-3 are within the noise floor and not worth continuing
   ! for (3 consecutive non-improving iterations still trigger the
   ! stagnation exit).
   !
   double precision :: jfnk_gmres_stag_tol = 1.0d-3
   !
   ! One-shot diagnostic flag (Finding 3, review 2026-09-10):
   ! When .true., jfnk_newton_diagnostic_v1 and jfnk_gmres_dim_experiment_v1
   ! are called on the first Newton step.  Set from JFNK_DIAGNOSTIC VADAT
   ! key via solve_for_pops_jfnk_v1.f.  Default .false. for production.
   !
   logical          :: jfnk_diagnostic_flag = .false.
   !
   ! Persistent workspace for Cauchy direction (avoids stack allocation)
   !
   double precision, allocatable :: jfnk_w_cauchy(:)  ! M^{-1}*F direction
   !
contains
   !
   ! jfnk_pc_set_options_v1: copy the parsed VADAT values into the
   ! solver-module mirrors.  rd_control_variables.f cannot import the
   ! mirrors directly: their names collide with the
   ! CONTROL_VARIABLE_MOD declarations (ifort #6405), so the sync goes
   ! through this routine.
   !
   subroutine jfnk_pc_set_options_v1(ptc_signed, reg_zerocol, &
        tiny_equil, rcond_gate, rcond_min, itry_cap)
      implicit none
      logical, intent(in) :: ptc_signed, reg_zerocol, tiny_equil, &
           rcond_gate
      double precision, intent(in) :: rcond_min
      integer, intent(in) :: itry_cap
      jfnk_pc_signed = ptc_signed
      jfnk_pc_reg_zerocol = reg_zerocol
      jfnk_pc_tiny_equil = tiny_equil
      jfnk_rcond_gate = rcond_gate
      jfnk_rcond_min = rcond_min
      jfnk_itry_cap = itry_cap
   end subroutine jfnk_pc_set_options_v1
!
   ! jfnk_eval_cap_v1: effective per-activation residual-eval cap.
   ! JFNK_EVAL_LIM (int_mod) is the activation-local budget the driver
   ! trims per rung; JFNK_MAX_EVALS is the VADAT control.  The
   ! sentinel in JFNK_RESIDUAL_V1 must see the tighter of the two,
   ! with 0 disabling a cap independently.  Pure so the unit tests
   ! can lock the semantics directly.
   !
   pure integer function jfnk_eval_cap_v1(eval_lim, max_evals) &
         result(cap)
      implicit none
      integer, intent(in) :: eval_lim, max_evals
      cap = max_evals
      if(eval_lim > 0 .and. (cap <= 0 .or. eval_lim < cap)) &
         cap = eval_lim
   end function jfnk_eval_cap_v1
!
   ! jfnk_cpl_sob_cap_v1: effective residual-eval cap for a coupled-SOB
   ! activation ([JFNK_CPL_SOB_MAX_EVALS], default 0 = inherit).  Layers
   ! the coupled cap over the usual eval_lim/max_evals pair; 0 disables
   ! a cap independently at each level.  (23-Oct-2026)
   !
   pure integer function jfnk_cpl_sob_cap_v1(eval_lim, max_evals, &
        cpl_max) result(cap)
      implicit none
      integer, intent(in) :: eval_lim, max_evals, cpl_max
      cap = jfnk_eval_cap_v1(eval_lim, max_evals)
      if(cpl_max > 0 .and. (cap <= 0 .or. cpl_max < cap)) cap = cpl_max
   end function jfnk_cpl_sob_cap_v1
!
   ! jfnk_cpl_sob_admit_v1: reach-fail fallback admission predicate
   ! ([JFNK_COUPLED_SOB], default F; 23-Oct-2026).  The coupled tier
   ! fires under SOB only when the frozen tier just certified the
   ! residual unreachable (reach_skip), the activation is not already
   ! coupled, no fallback fired yet this activation, the admission is
   ! not a dirchk diagnostic call, the coupled-SOB config is armed,
   ! the consistency check has not latched dead, and the sweep
   ! callback is registered.  Extracted pure so the unit tests lock
   ! the gate semantics; JFNK_SOLVE_CORR_V1 calls it at the skip site.
   !
   pure logical function jfnk_cpl_sob_admit_v1(reach_skip, &
        coupled_loc, fired, dirchk, sob, coupled, cpl_sob, dead, &
        cb_assoc) result(arm)
      implicit none
      logical, intent(in) :: reach_skip, coupled_loc, fired, dirchk
      logical, intent(in) :: sob, coupled, cpl_sob, dead, cb_assoc
      arm = reach_skip .and. .not. coupled_loc .and. .not. fired &
         .and. .not. dirchk .and. sob .and. coupled .and. cpl_sob &
         .and. .not. dead .and. cb_assoc
   end function jfnk_cpl_sob_admit_v1
!
end module jfnk_solver_v1_mod


!----------------------------------------------------------------------
! JFNK_PTC_SET_ANCHOR_V1: re-anchor the pseudo-time term at x (the
! current Newton iterate x_k).  SER form (Kelley & Keyes 1998; Knoll &
! Keyes 2004 Sec 4.2): G(x) = F(x) + (x - x_k)*s/dtau vanishes at the
! anchor, so the stored residual at the linearization point is F, not
! G -- all Newton norms, ared/pred, EW forcing and the SER dtau update
! then act on |F|.  No-op unless the PTC machinery is armed
! ([JFNK_PTC], default T; USE_JFNK=F never reaches this).
! edited by Devin, 06-Oct-2026
!----------------------------------------------------------------------
subroutine jfnk_ptc_set_anchor_v1(n, x)
   use jfnk_solver_v1_mod
   implicit none
   integer,          intent(in) :: n
   double precision, intent(in) :: x(n)
   integer :: ios
   !
   if(.not. (jfnk_ptc_active .and. jfnk_ptc_scale_ready)) return
   if(allocated(jfnk_ptc_x0)) then
      if(size(jfnk_ptc_x0) /= n) deallocate(jfnk_ptc_x0)
   end if
   if(.not. allocated(jfnk_ptc_x0)) then
      allocate(jfnk_ptc_x0(n), stat=ios)
      if(ios /= 0) then
         jfnk_ptc_anchor_ready = .false.
         return
      end if
   end if
   jfnk_ptc_x0 = x
   jfnk_ptc_anchor_ready = .true.
end subroutine jfnk_ptc_set_anchor_v1


!----------------------------------------------------------------------
! JFNK_PTC_STRIP_V1: g := g - (x - x0)*s/dtau, so that g = F(x) when the
! residual callback returned g = G(x).  newton_hook_v1 calls it after
! every residual evaluation so jfnk_fx/jfnk_tol are F and |F|.  No-op
! unless the PTC term was actually included (same readiness gate as
! JFNK_RESIDUAL_V1, plus a matching anchor and sizes).
! edited by Devin, 06-Oct-2026
!----------------------------------------------------------------------
subroutine jfnk_ptc_strip_v1(n, x, g)
   use jfnk_solver_v1_mod
   implicit none
   integer,          intent(in)    :: n
   double precision, intent(in)    :: x(n)
   double precision, intent(inout) :: g(n)
   integer :: i
   !
   if(.not. (jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
             jfnk_ptc_anchor_ready)) return
   if(.not. (allocated(jfnk_ptc_scale) .and. allocated(jfnk_ptc_x0))) &
      return
   if(size(jfnk_ptc_scale) /= n .or. size(jfnk_ptc_x0) /= n) return
   do i = 1, n
      g(i) = g(i) - (x(i) - jfnk_ptc_x0(i)) * jfnk_ptc_scale(i) / &
         max(jfnk_ptc_dtau, 1d-8)
   end do
end subroutine jfnk_ptc_strip_v1


!----------------------------------------------------------------------
! JFNK_PTC_STRIP_DIR_V1: js := js - s*scale/dtau.  The finite-difference
! JVP of G is J_F*s + s*scale/dtau (the anchor offset cancels in the
! difference); stripping recovers the F-Jacobian action so the physical
! predicted residual ||F - J_F*s|| and the descent test <F, J_F*s> are
! consistent with the stripped ared.  The GMRES operator itself keeps
! the diagonal -- it correctly solves (J_F + D/dtau) s = F(x_k).
! Same readiness gate as jfnk_ptc_strip_v1.
! edited by Devin, 06-Oct-2026
!----------------------------------------------------------------------
subroutine jfnk_ptc_strip_dir_v1(n, s, js)
   use jfnk_solver_v1_mod
   implicit none
   integer,          intent(in)    :: n
   double precision, intent(in)    :: s(n)
   double precision, intent(inout) :: js(n)
   integer :: i
   !
   if(.not. (jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
             jfnk_ptc_anchor_ready)) return
   if(.not. allocated(jfnk_ptc_scale)) return
   if(size(jfnk_ptc_scale) /= n) return
   do i = 1, n
      js(i) = js(i) - s(i) * jfnk_ptc_scale(i) / &
         max(jfnk_ptc_dtau, 1d-8)
   end do
end subroutine jfnk_ptc_strip_dir_v1


!----------------------------------------------------------------------
! JFNK_ALLOC_WORKSPACE: Allocate or resize persistent workspace arrays
! for newton_hook_v1 and gmres_m_v1.  Called at the start of
! newton_hook_v1 to ensure workspace is available.
! audit §3, 09-Sep-2026
!----------------------------------------------------------------------
subroutine jfnk_alloc_workspace(n, luer)
   use jfnk_solver_v1_mod
   implicit none
   integer, intent(in) :: n, luer
   integer :: ios
   !
   if(jfnk_w_alloc_n == n) return  ! already allocated at this size
   !
   if(allocated(jfnk_w_x)) deallocate(jfnk_w_x, jfnk_w_fx, jfnk_w_x2, &
      jfnk_w_fx2, jfnk_w_s, jfnk_w_fxpc, jfnk_w_gmres_w, jfnk_w_gmres_z, &
      jfnk_w_cauchy)
   !
   allocate(jfnk_w_x(n), jfnk_w_fx(n), jfnk_w_x2(n), jfnk_w_fx2(n), &
      jfnk_w_s(n), jfnk_w_fxpc(n), jfnk_w_gmres_w(n), jfnk_w_gmres_z(n), &
      jfnk_w_cauchy(n), stat=ios)
   if(ios /= 0) then
      if(luer > 0) write(luer,*) 'jfnk_alloc_workspace: allocation failed'
      stop 'jfnk_alloc_workspace: alloc failed'
   end if
   jfnk_w_alloc_n = n
end subroutine jfnk_alloc_workspace


!----------------------------------------------------------------------
! NEWTON_HOOK_V1: Newton iteration with GMRES-hookstep globalization.
!
! Solves F(x) = 0 using:
!   - GMRES(m) for the Krylov inner solve
!   - Hookstep (trust-region) for globalization
!   - Finite-difference Jacobian-vector products
!
! Interface (matches the JFNK-Hookstep reference pattern):
!   f(n,x, y)       — residual: y = F(x)
!   df(n,x, y)      — JVP: y = J(x)*x_vec  (finite difference internally)
!   dfp(n, x)       — preconditioner: x = M^{-1} * x  (in-place)
!   sub()           — called at end of each Newton iteration (diagnostics)
!   dp(n,a,b)       — dot product: d = a.b
!
! editted by Sean Brennan (V4, 2026-09)
!----------------------------------------------------------------------
subroutine newton_hook_v1(f, df, dfp, sub, dp, m, n, gtol, &
                            tol, del, mndl, mxdl, nits, info, luer)
   use jfnk_solver_v1_mod
   use jfnk_cmfgen_int_mod, only: JFNK_NT, JFNK_COUPLED_LOC
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   external                        :: f, df, dfp, sub
   double precision, external      :: dp
   integer,          intent(in)    :: m, n
   double precision, intent(in)    :: gtol, tol, del, mndl, mxdl
   integer,          intent(in)    :: nits
   integer,          intent(inout) :: info
   integer,          intent(in)    :: luer
   !
   ! W25: SAVE'd GMRES basis buffer -- one heap allocation per
   ! (n,m) instead of one per Newton step (audit 17-Oct-2026).
   double precision, allocatable, save :: v(:,:)
   double precision :: tol_, del__, tol__, tol_initial
   double precision :: mxdl_, ared, pred, snrm, s_inf
   double precision :: snrm_tr   ! 2026-09-13 TR-units fix: ||y|| norm the hookstep constrains
   double precision :: gres, gdel, h((m+1)*m)
   double precision :: del_try, scale_fac, s_inf_orig
   double precision :: eta_sav    ! EW Choice-2 safeguard gamma*eta_{k-1}^alpha
   double precision :: ref_norm
   double precision :: pred_res
   double precision :: pred_res_physical  ! ||F - J*s||_w (physical norm)
   double precision :: pred_res_gmres      ! ||M^{-1}(F-J*s)|| (preconditioned, diag only)
   double precision :: pred_floor          ! floor for safe ared/pred division
   double precision :: rho                 ! ared/pred ratio (dimensionless)
   logical          :: pred_negligible     ! |pred| <= pred_floor
   logical          :: maxch_converged
   integer :: ginfo, n_reject, ios_nw, i
   integer, parameter :: n_reject_max = 8
   !
   ! Line-search and Cauchy fallback variables (review §3, 2026-09-10)
   !
   logical          :: ls_accepted
   double precision :: phi_0, dphi_0, phi_alpha, alpha_ls, dtau_was
   double precision :: phi_ref   ! non-monotone LS reference (Oct-2026)
   integer          :: ls_step
   !
   ! Cauchy direction safeguards (audit 17-Sep-2026): the
   ! preconditioned-residual direction d_C = M^{-1}F now gets the same
   ! finiteness / inf-norm / runaway checks as the Newton direction.
   !
   logical          :: cauchy_ok
   double precision :: cauchy_nrm, s_inf_c
   !
   ! V4 A/B production diagnostic variables (2026-09-09)
   !
   double precision :: dir_deriv, cos_angle, norm_Js
   double precision :: max_pop_corr, max_t_corr
   integer :: n_gmres_restarts
   character(len=30) :: step_status
   !
   ! V4 descent check (2026-09-09):
   ! If the GMRES direction is not a descent direction (<F,Js>_w <= 0),
   ! tighten the forcing term and re-solve.  This prevents the
   ! Eisenstat-Walker formula from producing too-loose tolerances
   ! (eta ~ 0.9) that yield non-descent directions with poor
   ! preconditioners.
   !
   integer :: descent_retry_count
   integer, parameter :: descent_retry_max = 3
   !
   ! Heap-allocated workspace (audit §3, 09-Sep-2026):
   ! x_, fx_, x__, fx__, s, fx_pc were automatic arrays of size N,
   ! contributing ~480KB (N=10000) to ~1.9MB (N=40000) of stack.
   ! Making them ALLOCATABLE,SAVE moves them to the heap.
   !
   double precision, allocatable, save :: x_(:), fx_(:)
   double precision, allocatable, save :: x__(:), fx__(:)
   double precision, allocatable, save :: s(:), fx_pc(:)
   !
   ! V4 audit fix (2026-09-09): Physical predicted residual workspace.
   ! J_s = J*s (one extra JVP per candidate step), F_pred = F(x) - J*s.
   ! Used to compute pred_res_physical = ||F - J*s||_w in the SAME
   ! weighted physical norm as tol_ = ||F||_w, so that ared/pred is
   ! dimensionally consistent.  The Hessenberg least-squares residual
   ! from GMRES (pred_res) lives in the preconditioned space
   ! ||M^{-1}(F - J*s)|| and must NOT be used in ared/pred.
   !
   double precision, allocatable, save :: J_s(:), F_pred(:)
   !
   ! (Re)allocate if N changed.  The stat= variables must be reset
   ! before each conditional allocate: on a re-entry where the arrays
   ! already hold the right shape, no allocate runs and a stale nonzero
   ! ios_nw/ginfo would otherwise fabricate an alloc failure.
   if(allocated(x_)) then
      if(size(x_) .ne. n) deallocate(x_, fx_, x__, fx__, s, fx_pc, J_s, F_pred)
   end if
   ios_nw = 0
   if(.not. allocated(x_)) &
      allocate(x_(n), fx_(n), x__(n), fx__(n), s(n), fx_pc(n), &
               J_s(n), F_pred(n), stat=ios_nw)
   if(ios_nw /= 0) then
      if(info == 1) write(luer,*) 'newton_hook_v1: workspace alloc failed'
      info = 4
      jfnk_exit_reason = 'alloc_fail_workspace'
      return
   end if
   !
   if(allocated(v)) then
      if(size(v,1) .ne. n .or. size(v,2) .ne. m+1) deallocate(v)
   end if
   ginfo = 0
   if(.not. allocated(v)) allocate(v(n,m+1), stat=ginfo)
   if(ginfo /= 0) then
      if(info == 1) write(luer,*) 'newton_hook_v1: unable to allocate memory'
      info = 4
      jfnk_exit_reason = 'alloc_fail_v'
      return
   end if
   !
   jfnk_nits = 0
   jfnk_gits = 0
   jfnk_del  = del
   mxdl_     = mxdl
   ginfo     = info
   n_reject  = 0
   !
   ! Right preconditioning is supported (2026-10, right-PC TR-units
   ! fix): gmres_m_v1 stores the preconditioned basis Z = M^{-1} V and
   ! its Gram G = Z' W Z, and hook_step_g_v1 constrains ||s||_w =
   ! sqrt(y' G y).  This supersedes the 17-Sep-2026 audit disable
   ! (issue 16), which forced LEFT because the hookstep then bounded
   ! ||y|| = ||M s||_w instead of the physical step norm.
   !
   !
   ! Initial residual.  The driver's F(0) test and the baseline
   ! rebuild both leave jfnk_fx = F(jfnk_x) already -- when they did,
   ! skip the re-eval (a full formal sweep in coupled mode).
   !
   if(jfnk_fx_fresh) then
      jfnk_fx_fresh = .false.
   else
      call f(n, jfnk_x, jfnk_fx)
      jfnk_total_f_evals = jfnk_total_f_evals + 1
      ! edited by Devin, 06-Oct-2026 (SER anchor fix): the entry eval
      ! carries the pseudo-time term whenever it was added by the
      ! residual (anchor = 0 from JFNK_SET_PTC_SCALE_V1); strip it so
      ! jfnk_fx/jfnk_tol are F and |F|.  No-op when PTC is off.
      call jfnk_ptc_strip_v1(n, jfnk_x, jfnk_fx)
   end if
   jfnk_tol = dsqrt(dp(n, jfnk_fx, jfnk_fx))
   tol_initial = jfnk_tol
   ! Fresh non-monotone LS history per activation (Oct-2026).
   jfnk_ls_nmhist = -1d0
   jfnk_ls_nmidx  = 0
   ! EW Choice-2 memory is per-Newton-iteration, not per-activation:
   ! a stale eta from a previous solve would safeguard nothing and
   ! make the forcing path stateful across calls.
   jfnk_eta_prev  = 0d0
   !
   ! 2026-09-13 TR-units fix: for M≈J the full Newton step has
   ! ||y|| ≈ ||F||_w, so |F|/10 capped every first step to a 10%
   ! residual reduction.  Start at 2*|F| instead.
   if(del < 0d0)  jfnk_del = 2d0 * jfnk_tol
   ! Auto trust-region (del<0) must not discard a user-supplied mxdl:
   ! only expand the cap when the caller did not set one (mxdl<=0).
   ! (Oct-2026 JFNK_MXDL fix)
   if(del < 0d0 .and. mxdl <= 0d0)  mxdl_ = 1d99
   !
   if(info == 1) write(luer,'(A,I3,A,ES12.4,A,ES12.4)') &
      '  JFNK: nits=', jfnk_nits, ' res=', real(jfnk_tol), &
      ' maxch_est=', real(jfnk_maxch_est)
   if(info == 1) flush(luer)
   !
   call sub()
   x_    = jfnk_x
   fx_   = jfnk_fx
   tol_  = jfnk_tol
   tol__ = 1d99
   !
   ! edited by Devin, 06-Oct-2026 (SER anchor fix): anchor the
   ! pseudo-time term at the entry iterate x_k.  The entry eval above
   ! already ran with anchor 0 (jfnk_x = 0 there in production, so
   ! G = F anyway); this makes the anchor explicit for later evals.
   call jfnk_ptc_set_anchor_v1(n, x_)
   !
   if(jfnk_tol < tol) then
      if(info == 1) write(luer,*) '  JFNK: input already converged'
      info = 0
      jfnk_exit_reason = 'converged_at_entry'
      return
   end if
   !
   ! --- Main Newton loop ---
   !
   ! All loops in this subroutine are explicitly named to prevent
   ! bare cycle/exit from binding to the wrong DO construct.
   ! (Finding 1, review 2026-09-10: a bare `cycle` inside descent_retry
   !  was silently binding to the inner loop instead of the outer
   !  Newton loop, causing stale trust-region radii to be reused.)
   newton_loop: do while(.true.)
      !
      ! V4 descent check: reset retry counter at the start of each
      ! Newton iteration.
      !
      descent_retry_count = 0
      !
      ! Audit 17-Sep-2026: per-iteration failure flags.  Any matvec
      ! (JVP) or psolve (preconditioner) call below can set these;
      ! GMRES aborts with info=3 on either, and the post-solve handler
      ! rejects the direction and requests a preconditioner refresh.
      !
      jfnk_jvp_failed    = .false.
      jfnk_pc_clip_event = .false.
      ! edited by Devin, 05-Oct-2026: same per-iteration reset for the
      ! hookstep-recalculation flag consumed by the descent check below.
      jfnk_gmres_recomputed = .false.
      !
      ! editted by Devin, 04-Oct-2026: also catch jfnk_del == 0.  The
      ! shrink producers (0.5*min(del_try, snrm_tr)) yield exactly 0
      ! whenever GMRES returned no step (snrm_tr = 0), and a zero
      ! radius slipped the old jfnk_del > 0 test: hookstep then saw
      ! del == 0 (the no-trust-region sentinel), produced s = 0, and
      ! the pred_negligible branch looped with no reject counter.
      ! del /= 0 keeps the (unreachable) explicit no-TR input mode.
      if(del /= 0d0 .and. jfnk_del < mndl) then
         if(info == 1) write(luer,*) '  JFNK: trust region too small'
         info = 3
         jfnk_exit_reason = 'trust_region_too_small'
            return
      end if
      !
      ! Solve for hookstep s subject to |s| < del
      !
      ! Eisenstat-Walker forcing term (added 2026-09-07):
      ! Instead of a fixed GMRES tolerance gres = gtol * |F|, we compute
      ! an adaptive forcing term eta_k that prevents over-solving the
      ! inner GMRES system when the Newton iteration is far from the root.
      !
      ! On the first Newton step (jfnk_nits == 0), we use the user-specified
      ! gtol as the initial forcing term.  On subsequent steps, we compute:
      !
      !   eta_k = gamma * (|F_k| / |F_{k-1}|)^alpha
      !
      ! with the safeguard eta_k = max(eta_k, 0.5 * |F_k| / |F_{k-1}|)
      ! to prevent the "oversolving" phenomenon where GMRES solves too
      ! accurately relative to the Newton convergence rate.
      !
      ! The forcing term is clamped to [eta_min, eta_max] for stability.
      ! When del == 0 (no trust region, pure Newton), we use a very loose
      ! tolerance to take as many GMRES iterations as needed.
      !
      s        = 0d0
      del_try  = jfnk_del
      if(jfnk_nits == 0 .or. jfnk_tol_prev .le. 0d0 .or. &
           .not. jfnk_ew_adapt) then
         ! First Newton step (or fixed-forcing mode): user-specified
         ! gtol as initial eta, but clamped to [eta_min, eta_max] like
         ! every later step (audit 17-Sep-2026 -- an unclamped gtol
         ! could bypass the forcing-term bounds entirely).
         jfnk_eta = min(max(gtol, jfnk_eta_min), jfnk_eta_max)
      else
         ! Eisenstat-Walker Choice 2 (Eisenstat & Walker 1994, SIAM J.
         ! Optim. 4, 393, Eq 4.2): eta_k = gamma*(|F_k|/|F_{k-1}|)^alpha.
         ! Their safeguard against over-solving is on eta_{k-1}, not on
         ! the residual ratio: when gamma*eta_{k-1}^alpha > 0.1 the
         ! contraction theory still holds with
         ! eta_k = max(eta_k, gamma*eta_{k-1}^alpha).  The previous code
         ! used max(eta_k, 0.5*|F_k|/|F_{k-1}|) -- a ratio safeguard that
         ! EW show is NOT covered by the Choice-2 convergence proof and
         ! that pins eta near the contraction rate, defeating the
         ! adaptive tolerance.  (Oct-2026 EW fix)
         jfnk_eta = jfnk_eta_gamma * &
            (jfnk_tol / jfnk_tol_prev)**jfnk_eta_alpha
         eta_sav = jfnk_eta_gamma * jfnk_eta_prev**jfnk_eta_alpha
         if(eta_sav > 0.1d0) jfnk_eta = max(jfnk_eta, eta_sav)
         ! Clamp to [eta_min, eta_max] for numerical stability
         jfnk_eta = max(jfnk_eta, jfnk_eta_min)
         jfnk_eta = min(jfnk_eta, jfnk_eta_max)
      end if
      jfnk_eta_prev = jfnk_eta
      !
      ! Bug fix (09-Sep-2026): EW forcing-term tolerance must be in the
      ! same units as the GMRES residual.  Under left preconditioning
      ! (the default, jfnk_use_left_pc=.true.), GMRES measures
      ! ||M^{-1}(b - Ax)||, not ||b - Ax||.  Using jfnk_tol (=||F||,
      ! physical units) as the reference makes gres off by the
      ! conditioning of M, so GMRES can quit far too early or grind
      ! far too long relative to the intended inexact-Newton accuracy.
      !
      ! Fix: compute ref_norm = ||M^{-1} F(x)|| via one extra psolve
      ! call when left-pc is active.  This is cheap (one preconditioner
      ! application per Newton step, vs m applications per GMRES cycle).
      !
      if(jfnk_use_left_pc) then
         fx_pc = fx_
         call dfp(n, fx_pc)
         !
         ! Audit 17-Sep-2026: if the preconditioner clipped during this
         ! forcing-norm evaluation, M^{-1} is nonlinear this step.
         ! Reject the step, shrink the trust region, request a refresh.
         !
         if(jfnk_pc_clip_event) then
            ! editted by Devin, 04-Oct-2026: abort instead of retry.
            ! The clip is a pure function of the (unchanged) LU factors
            ! and the input vector; shrinking del cannot change it
            ! (del is only consulted after GMRES completes), and the
            ! refresh flag is only consumed via sub() on accepted
            ! steps, so every "refresh + retry" replayed a bit-
            ! identical psolve and burned n_reject_max for nothing.
            if(info == 1) write(luer,'(A)') &
               '  JFNK: preconditioner clip in forcing norm -- &
               &deterministic, aborting Newton'
            jfnk_x   = x_
            jfnk_fx  = fx_
            jfnk_tol = tol_
            info = 3
               jfnk_exit_reason = 'precond_clip_forcing_norm'
            return
         end if
         ref_norm = dsqrt(dp(n, fx_pc, fx_pc))
      else
         ref_norm = jfnk_tol
      end if
      gres     = jfnk_eta * ref_norm
      gdel     = jfnk_del
      ! jfnk_gits is both the iteration budget in and the used count
      ! out of gmres_m_v1.  A ginfo==2 hookstep recompute performs no
      ! Arnoldi work, so the prior solve's count is deliberately
      ! retained rather than re-armed to m.
      if(ginfo /= 2) jfnk_gits = m
      if(del == 0d0) jfnk_gits = 9999
      !
      ! Override with user-specified max GMRES iterations if set
      ! (added 2026-09-07: JFNK_MAX_GMRES control variable).
      ! When jfnk_max_gmres_its > 0, it caps the total GMRES iterations
      ! per Newton step, allowing multiple restart cycles up to the limit.
      !
      if(jfnk_max_gmres_its > 0 .and. ginfo /= 2) jfnk_gits = jfnk_max_gmres_its
      !
      if(info == 1 .and. ginfo /= 2) then
         if(jfnk_ptc_active) then
            write(luer,'(A)') ' '
            write(luer,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
               '  JFNK: eta=', real(jfnk_eta), &
               '   gres=', real(gres), &
               '   dtau=', real(jfnk_ptc_dtau)
         else
            write(luer,'(A)') ' '
            write(luer,'(A,ES10.3,A,ES10.3)') &
               '  JFNK: eta=', real(jfnk_eta), '   gres=', real(gres)
         end if
      end if
      !
      ! V4 descent check (2026-09-09):
      ! Wrap the GMRES solve in a retry loop.  If the resulting direction
      ! is not a descent direction (<F,Js>_w <= 0), tighten the forcing
      ! term and re-solve.  This prevents the Eisenstat-Walker formula
      ! from producing too-loose tolerances (eta ~ 0.9) that yield
      ! non-descent directions with poor preconditioners.
      !
      descent_retry: do
      !
      jfnk_gmres_m_dim = m
      call gmres_m_v1(m, n, s, fx_, df, dfp, dp, h, v, gres, gdel, &
                       pred_res, jfnk_gits, ginfo, luer)
      !
      ! Issue 3: request LU refresh when GMRES used most of its budget.
      !
      jfnk_last_gmres_its = jfnk_gits
      if(jfnk_gmres_m_dim > 0 .and. &
         dble(jfnk_last_gmres_its) >= &
         jfnk_pc_ref_frac * dble(jfnk_gmres_m_dim)) then
         jfnk_need_precond_refresh = 1
      end if
      !
      ! Audit 17-Sep-2026: Krylov-solve failure (info=3 from gmres_m_v1).
      ! A matvec produced a non-finite JVP action, or the preconditioner
      ! clipped mid-solve making M^{-1} nonlinear.  The partially built
      ! Krylov space (and any step s derived from it) is unusable -- do
      ! NOT let it reach the step-cap/runaway/hookstep code, and do NOT
      ! reuse the space (ginfo=1 forces a fresh GMRES next pass).
      !
      if(ginfo == 3) then
         ! editted by Devin, 04-Oct-2026: abort instead of retry.  The
         ! failure (JVP non-finite or precond clip) is a deterministic
         ! function of x_, fx_ and the unchanged preconditioner --
         ! neither a halved del (only consulted inside the completed
         ! GMRES) nor the refresh flag (never consumed on reject
         ! paths) can alter the next GMRES pass, so the old loop
         ! replayed an identical abort n_reject_max times per solve.
         if(info == 1) write(luer,'(A)') &
            '  JFNK: Krylov solve failed (JVP/preconditioner) -- &
            &deterministic, aborting Newton'
         jfnk_x   = x_
         jfnk_fx  = fx_
         jfnk_tol = tol_
         info = 3
            jfnk_exit_reason = 'krylov_jvp_or_pc_failure'
         return
      end if
      !
      ! Bug fix: do NOT overwrite ginfo (GMRES return status) with info
      ! (Newton verbose flag).  ginfo must retain the GMRES status (0,
      ! 1, or 2) for the hookstep recalculation logic at the top of the
      ! next Newton iteration.  The old code had "ginfo = info" which
      ! always set ginfo=1, preventing the hookstep recalculation path
      ! (ginfo=2) from ever being taken.
      !
      ! |s| in the RSCALE-weighted inner product.  The hookstep bounds
      ! Krylov coefficients, not the physical step after M^{-1}.  A
      ! singular DIAG BA inverse can yield |s| ~ 1e31 with pred = 0,
      ! which previously set del = |s|/2 and looped forever (nits stuck
      ! at 0, F evaluated at POPS*(1+1e31)).
      ! editted by Sean Brennan (08-Sep-2026)
      snrm = dsqrt(dp(n, s, s))
      s_inf = maxval(abs(s))
      ! 2026-09-13 TR-units fix: the trust-region radius must be updated
      ! in the norm hook_step_v1 actually constrains — the Krylov
      ! coefficient norm ||y|| (= ||M*s||_w under right preconditioning),
      ! not the physical step norm ||s||_w.
      snrm_tr = jfnk_last_ynrm
      !
      ! One-time step diagnostic (09-Sep-2026)
      if(jfnk_nits == 0 .and. info == 1) then
         write(luer,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
            '  JFNK STEP DIAG: |s|_w=', real(snrm), &
            ' |s|_inf=', real(s_inf), ' |s|_max=', real(maxval(s))
         write(luer,'(A,ES14.6,A,ES14.6)') &
            '    s_min=', real(minval(s)), ' s_max=', real(maxval(s))
      end if
      !
      ! 2026-09-13 TR-units fix: GMRES stagnation (info=1) after a large
      ! linear-residual reduction is usually FD noise, not failure.
      ! If the predicted residual still beats 90% of the RHS norm, keep
      ! the inexact direction instead of aborting the whole Newton solve.
      ! edited by Devin, 05-Oct-2026 (stagnation-gate norm fix):
      ! pred_res is the Hessenberg least-squares residual in the units
      ! GMRES actually solves -- ||M^{-1}(F - J*s)|| under the forced
      ! left preconditioning -- so it must be compared against
      ! ref_norm = ||M^{-1}F||, not tol_ = ||F||_w.  The two differ by
      ! the conditioning of M, which made the old gate fire (or not
      ! fire) on the wrong threshold.
      ! editted by Devin, 12-Nov-2026: the Hessenberg residual is only a
      ! Krylov-space model; under a finite-difference matvec the basis
      ! need not span a fixed operator, and OUTGEN showed hess_res
      ! 15x-450x below the measured true_res on exits this gate used to
      ! accept.  Require the measured value below the same threshold,
      ! and treat a missing/non-finite measurement as unverifiable --
      ! abort rather than commit the solve to a bad direction.
      if(ginfo == 1 .and. pred_res < 0.9d0 * ref_norm .and. &
         jfnk_gmres_true_res_ok .and. &
         ieee_is_finite(jfnk_gmres_true_res) .and. &
         jfnk_gmres_true_res < 0.9d0 * ref_norm) then
         if(info == 1) write(luer,'(A,ES12.4,A,ES12.4)') &
            '  JFNK: GMRES stagnated; using inexact direction &
            &(pred_res/ref=', real(pred_res / max(ref_norm, 1d-30)), &
            ' true_res/ref=', &
            real(jfnk_gmres_true_res / max(ref_norm, 1d-30))
         ginfo = 0
      else if(ginfo == 1) then
         if(info == 1) write(luer,'(A)') &
            '  JFNK: GMRES stagnated -- aborting Newton (Lambda fallback)'
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         info = 3
            jfnk_exit_reason = 'gmres_stagnation_abort'
         return
      end if
      !
      ! Direction-preserving inf-norm scaling (09-Sep-2026):
      ! The DIAG BA preconditioner with ridge produces Newton steps with
      ! huge |s_inf| (up to 400+) for species with small BA diagonal.
      ! The RSCALE-weighted trust region doesn't penalize these because
      ! RSCALE is small for trace species.  Previously the runaway check
      ! rejected the ENTIRE step, causing trust-region collapse.
      !
      ! Fix: scale the ENTIRE step by min(1, S_MAX_SCALE/s_inf) instead
      ! of rejecting it.  This preserves the Newton direction (unlike
      ! per-element clipping which changes the direction) while limiting
      ! the maximum fractional change to S_MAX_SCALE per component.
      ! A small step in an approximately-correct direction can still
      ! reduce the residual (damped Newton principle).
      !
      s_inf_orig = s_inf
      if(s_inf > jfnk_step_inf_cap) then
! jfnk_step_inf_cap (default 10.0d0 = MAX_JFNK_COR, 1.0D1): steps larger
! than a factor-of-10 population correction are unphysical for the
! clipping that is applied after the Newton update.  Scale the whole
! step back so the inf-norm equals the maximum allowed correction
! (audit 2026-09-10).
         scale_fac = jfnk_step_inf_cap / s_inf
         s = s * scale_fac
         snrm = dsqrt(dp(n, s, s))
         snrm_tr = snrm_tr * scale_fac   ! 2026-09-13 TR-units fix
         s_inf = jfnk_step_inf_cap
         jfnk_n_step_scaled = jfnk_n_step_scaled + 1
         if(info == 1) write(luer,'(A,ES10.3,A,ES10.3,A,ES10.3)') &
            '  JFNK: scaled step by ', real(scale_fac), &
            ' |s|=', real(snrm), ' |s|_inf=', real(s_inf)
      end if
      !
      ! FIX_T_AUTO: damp/cap the T components of the finished Newton
      ! step BEFORE the JVP prediction and the trial eval, so every
      ! candidate state the line search sees is already bounded.
      !
      call jfnk_tauto_damp_v1(n, s)
      !
      ! jfnk_tauto_damp_v1 modifies the T components of s in place;
      ! the norms below must describe the step actually taken, not the
      ! pre-damping candidate.  snrm_tr = snrm: the damped s IS the
      ! trust-region-constrained quantity from here on (the Cauchy
      ! path already recomputes its own cauchy_nrm).
      !
      snrm = dsqrt(dp(n, s, s))
      s_inf = maxval(abs(s))
      snrm_tr = snrm
      !
      ! Weighted |s| can still explode if RSCALE is large.  Reject if
      ! the overall weighted norm is truly catastrophic -- i.e. far
      ! above the ||M^{-1} F||_w scale the step is built at.  An
      ! absolute cap is meaningless under non-Euclidean RSCALE (this
      ! model: |F|_w ~ 5e10, every step read as "runaway" and the
      ! solver aborted before a single trial F evaluation).  Not
      ! flag-gated: the relative check reproduces the intended
      ! semantics under every JFNK_RSCALE mode.
      if(snrm > jfnk_runaway_fac * max(ref_norm, 1.0d-30)) then
         !
         ! Fix (2026-09-11): the retry must produce a DIFFERENT, smaller
         ! trial step.  Previously this halved del_try — the trust-region
         ! radius, initialised to ~|F|/10 — which is typically orders of
         ! magnitude larger than |s|.  The hookstep constraint in
         ! gmres_m_v1/hook_step_v1 only binds when del < |s|, so halving
         ! del_try changed nothing; and with ginfo=1 the next pass
         ! re-solved the SAME unconstrained GMRES system, returning a
         ! bit-identical step that was rejected again.  The loop was
         ! therefore guaranteed to burn all n_reject_max retries and
         ! abort without ever testing a constrained step.
         !
         ! Now: shrink the trust radius BELOW the rejected step norm
         ! (same convention as the ared<0 branch: del = 0.5*|s|) and set
         ! ginfo=2 to request a hookstep recalculation, which reuses the
         ! existing Krylov space (h, v, jj, beta) to solve the
         ! trust-region-constrained least-squares problem |s| <= del
         ! at no extra JVP cost.
         !
         jfnk_del = 0.5d0 * min(del_try, snrm_tr)   ! 2026-09-13 TR-units fix
         if(info == 1) write(luer,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
            '  JFNK: runaway |s|=', real(snrm), ' |s|_inf=', real(s_inf), &
            ' del=', real(del_try), ' -> hookstep del=', real(jfnk_del)
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         n_reject = n_reject + 1
         if(n_reject >= n_reject_max .or. jfnk_del < mndl) then
            if(info == 1) write(luer,'(A)') &
               '  JFNK: too many rejected steps — aborting Newton'
            info = 3
               jfnk_exit_reason = 'too_many_rejected_steps'
            return
         end if
         ginfo = 2
         cycle newton_loop
      end if
      !
      ! V4 audit fix (2026-09-09): Compute the PHYSICAL predicted residual
      ! ||F - J*s||_w using one extra JVP, instead of using the Hessenberg
      ! least-squares residual from GMRES.  Under left preconditioning, the
      ! Hessenberg residual is ||M^{-1}(F - J*s)||, which is in a different
      ! norm/space than ||F||_w.  Using it as pred_res in ared/pred mixes
      ! physical and preconditioned norms, causing the trust region to
      ! collapse even when the Newton direction is good.
      !
      ! The linearized model for x_new = x - s is:
      !   F(x - s) ≈ F(x) - J*s
      ! so the predicted residual is ||F(x) - J*s||_w in the same weighted
      ! physical norm as ||F(x)||_w = tol_.
      !
      ! jfnk_x and jfnk_fx are still x_ and fx_ at this point (GMRES does
      ! not modify them), so the JVP evaluates J at the correct point.
      !
      call df(n, s, J_s)
      jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
      !
      ! Audit 17-Sep-2026: if this JVP failed, J_s is unreliable and so
      ! is pred_res_physical / dir_deriv.  Same failure path as a
      ! Krylov abort: reject, shrink, fresh GMRES next pass.
      !
      if(jfnk_jvp_failed) then
         ! editted by Devin, 04-Oct-2026: abort instead of retry --
         ! same determinism argument as the ginfo==3 path above.
         if(info == 1) write(luer,'(A)') &
            '  JFNK: JVP failure in pred-residual eval -- &
            &deterministic, aborting Newton'
         jfnk_x   = x_
         jfnk_fx  = fx_
         jfnk_tol = tol_
         info = 3
            jfnk_exit_reason = 'jvp_failure_pred_eval'
         return
      end if
      ! edited by Devin, 06-Oct-2026 (SER anchor fix): df returns the
      ! FD action of G, J_G*s = J_F*s + s*scale/dtau.  Strip the
      ! diagonal so F_pred, pred_res_physical, norm_Js and dir_deriv
      ! are in F terms, consistent with the stripped ared.
      call jfnk_ptc_strip_dir_v1(n, s, J_s)
      do i = 1, n
         F_pred(i) = fx_(i) - J_s(i)
      end do
      pred_res_physical = dsqrt(dp(n, F_pred, F_pred))
      !
      ! Keep the GMRES Hessenberg residual for diagnostics only.
      !
      pred_res_gmres = pred_res
      !
      ! V4 A/B diagnostic: compute ||J*s||_w and GMRES restart count
      !
      norm_Js = dsqrt(dp(n, J_s, J_s))
      n_gmres_restarts = jfnk_gits / max(m, 1)
      if(mod(jfnk_gits, max(m,1)) > 0) n_gmres_restarts = n_gmres_restarts + 1
      n_gmres_restarts = max(0, n_gmres_restarts - 1)
      !
      ! V4 descent check (2026-09-09):
      ! If the GMRES direction is not a descent direction (<F,Js>_w <= 0),
      ! tighten the forcing term and re-solve.  This is a standard
      ! inexact-Newton technique that prevents the Eisenstat-Walker
      ! formula from producing too-loose tolerances (eta ~ 0.9) that
      ! yield non-descent directions with poor preconditioners.
      !
      ! The descent check uses the PHYSICAL inner product <F, J*s>_w,
      ! not the preconditioned one.  A direction can be a descent
      ! direction in the preconditioned norm but not in the physical
      ! norm when M is a poor approximation of J.
      !
      ! Skip the check for hookstep recalculations because the Krylov
      ! basis is fixed: tightening eta would discard the recalculated
      ! step and force a full fresh solve for no gain.
      ! edited by Devin, 05-Oct-2026 (recompute-pass descent-skip fix):
      ! ginfo cannot identify those passes -- gmres_m_v1 resets info=0
      ! before returning from its info==2 recompute path, so the old
      ! ginfo /= 2 test never fired and recompute passes were needlessly
      ! re-solved.  The module flag jfnk_gmres_recomputed survives the
      ! return and is used instead.
      !
      dir_deriv = dp(n, fx_, J_s)
      !
      ! Bail when the measured exit residual already shows the Krylov
      ! basis amplifies the operator: jfnk_gmres_true_res is the FD
      ! residual of the assembled step, set on non-converged exits, so
      ! true_res > |F| means the returned vector cannot be improved by
      ! a tighter forcing term -- re-solving burns a full GMRES budget
      ! on the same broken map (observed in OUTGEN: hess_res 1.6e7 vs
      ! |F-Js| 5.0e15; the eta-tightened re-solve failed the same way).
      ! A converged solve with |F-Js|>|F| is ordinary overshoot on a
      ! curved residual -- the retry shortens the step and recovers
      ! (T2 Rosenbrock regressed to an infinite stall when this gate
      ! used pred_res_physical).  The step still reaches the trial
      ! eval below, where the ared <= 0 rejection shrinks del and
      ! requests a hookstep recompute on the existing basis at no eval
      ! cost.
      !
      if(dir_deriv <= 0d0 .and. .not. jfnk_gmres_recomputed .and. &
         jfnk_eta > jfnk_eta_min .and. &
         descent_retry_count < descent_retry_max .and. &
         jfnk_gmres_true_res_ok .and. &
         ieee_is_finite(jfnk_gmres_true_res) .and. &
         jfnk_gmres_true_res > ref_norm) then
         if(info == 1) write(luer,'(A,ES10.3,A,ES10.3,A)') &
            '  JFNK: non-descent dir (<F,Js>=', real(dir_deriv), &
            ') with true_res/|F|=', &
            real(jfnk_gmres_true_res / max(ref_norm, 1d-30)), &
            ' -- keeping direction, not retrying'
      else if(dir_deriv <= 0d0 .and. .not. jfnk_gmres_recomputed .and. &
         jfnk_eta > jfnk_eta_min .and. &
         descent_retry_count < descent_retry_max) then
         !
         descent_retry_count = descent_retry_count + 1
         jfnk_n_descent_retry = jfnk_n_descent_retry + 1
         jfnk_eta = max(jfnk_eta_min, jfnk_eta * 0.1d0)
         gres = jfnk_eta * ref_norm
         ginfo = 0  ! silent fresh re-solve
         if(jfnk_max_gmres_its > 0) then
            jfnk_gits = jfnk_max_gmres_its
         else
            jfnk_gits = m
         end if
         if(info == 1) write(luer,'(A,ES10.3,A,I2,A,ES12.4)') &
            '  JFNK: non-descent dir (<F,Js>=', real(dir_deriv), &
            '), tightening eta to ', descent_retry_count, &
            ' eta=', real(jfnk_eta)
         cycle descent_retry
      end if
      !
      ! Exit the descent_retry loop — the direction is accepted.
      ! Without this exit, the do loop is infinite and the trust-region
      ! logic below is unreachable.
      exit descent_retry
      !
      end do descent_retry
      !
      ! V4 diagnostic (2026-09-09): Comprehensive Newton direction diagnosis.
      ! Printed once on the first Newton step to understand why the
      ! computed direction does not provide a useful physical reduction.
      ! Covers: sign verification, directional derivative, component-wise
      ! analysis, alpha-sweep, RSCALE audit, JVP-vs-BA comparison.
      !
      ! DISABLED for clean A/B production test (2026-09-09).
      ! The GMRES dimension experiment and Newton direction diagnostic
      ! are diagnostic-only and must not run in production comparisons.
      ! Re-enable by setting JFNK_DIAGNOSTIC = T in VADAT.
      ! (Finding 3, review 2026-09-10: implemented as a real runtime
      !  flag read from JFNK_DIAG VADAT key, not a source-level comment.)
      !
      if(jfnk_diagnostic_flag .and. jfnk_nits == 0 .and. info == 1) then
         call jfnk_newton_diagnostic_v1(n, s, fx_, J_s, F_pred, &
              dp, df, luer)
         call jfnk_gmres_dim_experiment_v1(n, x_, fx_, df, dfp, dp, &
              gres, jfnk_del, luer)
      end if
      !
      jfnk_x = x_ - s
      !
      ! Evaluate new residual
      !
      call f(n, jfnk_x, jfnk_fx)
      jfnk_total_f_evals = jfnk_total_f_evals + 1
      ! edited by Devin, 06-Oct-2026 (SER anchor fix): strip the
      ! pseudo-time term so jfnk_fx/jfnk_tol are F and |F|.
      call jfnk_ptc_strip_v1(n, jfnk_x, jfnk_fx)
      !
      jfnk_tol = dsqrt(dp(n, jfnk_fx, jfnk_fx))
      ared = tol_ - jfnk_tol
      pred = tol_ - pred_res_physical
      !
      ! V4 audit fix (2026-09-09): Safe ared/pred ratio.  If pred is
      ! negligibly small (|pred| <= pred_floor), the linearized model
      ! predicts essentially zero reduction.  Do NOT invent a large
      ! rho — instead classify as "negligible predicted reduction" and
      ! use a conservative trust-region response (reject step, shrink
      ! trust radius, reuse Krylov space via ginfo=2).
      !
      pred_floor = 1.0d-12 * max(tol_, 1.0d0)
      pred_negligible = (abs(pred) <= pred_floor)
      if(.not. pred_negligible) then
         rho = ared / pred
      else
         rho = 0.0d0   ! forces conservative trust-region response below
      end if
      !
      if(info == 1) then
         !
         ! V4 A/B production diagnostic (2026-09-09):
         ! Comprehensive per-Newton-iteration output for clean m=10 vs
         ! m=20 comparison.  Reports all quantities needed to determine
         ! whether the nonlinear solve progresses or stalls.
         !
         dir_deriv = dp(n, fx_, J_s)
         cos_angle = dir_deriv / max(tol_ * norm_Js, 1.0d-30)
         !
         ! Max population and temperature corrections from s
         ! JFNK_NT = unknowns per depth point; last 2 are Ne, T
         !
         if(JFNK_NT > 2) then
            max_pop_corr = 0.0d0
            max_t_corr = 0.0d0
            do i = 1, n
               if(mod(i-1, JFNK_NT) <= JFNK_NT - 3) then
                  max_pop_corr = max(max_pop_corr, abs(s(i)))
               else if(mod(i-1, JFNK_NT) == JFNK_NT - 1) then
                  max_t_corr = max(max_t_corr, abs(s(i)))
               end if
            end do
         else
            max_pop_corr = s_inf
            max_t_corr = 0.0d0
         end if
         !
         write(luer,'(A,I4,A)') '  --------- JFNK ITER ', jfnk_nits, &
            ' -----'
         write(luer,'(A,ES12.4)') '    weighted residual norm       :', &
            tol_
         write(luer,'(A,ES12.4)') '    weighted Newton step norm    :', &
            snrm
         write(luer,'(A,ES12.4)') '    linear model residual |F-Js| :', &
            pred_res_physical
         write(luer,'(A,ES12.4)') '    weighted |Js|                :', &
            norm_Js
         write(luer,'(A,ES12.4)') '    directional derivative <F,Js>:', &
            dir_deriv
         write(luer,'(A,ES12.4)') '    cos(F,Js), -1 is descent     :', &
            cos_angle
         if(pred_negligible) then
            write(luer,'(A)') '    acceptance ratio : predicted'// &
               ' decrease negligible'
         else
            write(luer,'(A,ES12.4)') '    actual residual decrease     :', &
               ared
            write(luer,'(A,ES12.4)') '    model-predicted decrease     :', &
               pred
            write(luer,'(A,ES12.4)') '    acceptance ratio ared/pred   :', &
               rho
         end if
         write(luer,'(A,ES12.4)') '    trust-region radius          :', &
            jfnk_del
         write(luer,'(A,I4)')    '    GMRES iterations             :', &
            jfnk_gits
         write(luer,'(A,I3)')    '    GMRES restarts               :', &
            n_gmres_restarts
         write(luer,'(A,ES12.4)') '    largest pop correction       :', &
            max_pop_corr
         write(luer,'(A,ES12.4)') '    largest T correction         :', &
            max_t_corr
         !
         ! Step disposition
         !
         if(jfnk_tol < tol_) then
            step_status = 'ACCEPTED'
         else if(del == 0d0) then
            step_status = 'ACCEPTED (full Newton)'
         else if(jfnk_tol > tol__) then
            step_status = 'ACCEPTED (previous)'
         else if(pred_negligible) then
            step_status = 'REJECTED (pred negligible)'
         else if(ared <= 0d0) then
            step_status = 'REJECTED (no decrease)'
         else if(rho < 0.75d0) then
            step_status = 'ACCEPTED (rho<0.75)'
         else if(snrm_tr < jfnk_del * 0.9d0) then   ! 2026-09-13 TR-units fix
            step_status = 'ACCEPTED (full Newton)'
         else
            step_status = 'ACCEPTED (enlarge TR)'
         end if
         write(luer,'(A,A)') '    step disposition              : ', &
            trim(step_status)
         flush(luer)
      end if
      !
      ! NaN/Inf guard (moved before descent check, 09-Sep-2026):
      ! If the residual is NaN/Inf, the descent check (pred <= 0) may
      ! fire first (when pred is finite but <= 0), taking the generic
      ! "shrink and retry" path without ever reaching the explicit NaN
      ! branch below.  This undercounts jfnk_nan_rejects.  Fix: check
      ! for NaN/Inf here, ahead of the descent check.
      !
      if(jfnk_tol /= jfnk_tol .or. jfnk_tol > 1d30) then
         jfnk_nan_rejects = jfnk_nan_rejects + 1
         if(info == 1) write(luer,*) '  JFNK: NaN/Inf in residual, &
            &rejecting step'
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         ! Constrain below the rejected step (same fix as the runaway
         ! branch, 2026-09-11): halving del_try does not bind when
         ! del >> |s|, and ginfo=1 re-solved the identical unconstrained
         ! step.  ginfo=2 recomputes a constrained hookstep from the
         ! existing Krylov space.
         jfnk_del = 0.5d0 * min(del_try, snrm_tr)   ! 2026-09-13 TR-units fix
         n_reject = n_reject + 1
         if(n_reject >= n_reject_max .or. jfnk_del < mndl) then
            if(info == 1) write(luer,'(A)') &
               '  JFNK: too many rejected steps — aborting Newton'
            info = 3
               jfnk_exit_reason = 'too_many_rejected_steps'
            return
         end if
         ginfo = 2
         cycle newton_loop
      end if
      !
      ! V4 audit fix (2026-09-09): The `pred <= 0` descent test was
      ! removed earlier because it was based on the (wrong) preconditioned
      ! Hessenberg residual.  Now that pred uses the physical predicted
      ! residual, pred <= 0 has real mathematical meaning: the linearized
      ! model does not predict a reduction in the physical residual norm.
      ! However, we let the standard trust-region ratio logic (ared < 0
      ! branch below) handle this case, because it uses ginfo=2 (hookstep
      ! recalculation, reusing the Krylov space) which is more efficient
      ! than a dedicated pred<=0 branch with ginfo=1 (full GMRES restart).
      !
      ! V4 audit fix (2026-09-09): Removed the duplicated NaN/Inf check
      ! that was here.  The primary NaN check above (before the descent
      ! comment block) is sufficient.  The second copy was a stale
      ! leftover from before the descent test was removed.
      !
      !================================================================
      ! Line-search globalization (review §3, 2026-09-10)
      !================================================================
      ! When the full Newton step (α=1) increases the residual (ared < 0)
      ! and the direction is a descent direction (dir_deriv > 0), try
      ! Armijo backtracking: α = β, β^2, β^3, ... before falling through
      ! to the trust-region rejection logic.
      !
      ! The Armijo condition for φ(α) = ½||F(x − αs)||_w² is:
      !   φ(α) ≤ φ(0) + c·α·φ'(0)
      ! where φ(0) = ½·tol_² and φ'(0) = −<F, Js>_w = −dir_deriv.
      !
      ! This avoids re-solving GMRES when the Newton direction is good
      ! but the full step is too large — exactly the situation observed
      ! in the A/B tests (optimal α ~ 10⁻³, full step increases ||F||
      ! by 26×).  Each backtracking step costs one residual evaluation,
      ! vs m JVPs for a GMRES re-solve.
      !
      ! If the line search fails, try the Cauchy (preconditioned-residual)
      ! direction d_C = M⁻¹F as a fallback (review §3, item 3).
      !
      ls_accepted = .false.
      !
      if(jfnk_do_line_search .and. ared < 0d0 .and. dir_deriv > 0d0) then
         !
         phi_0  = 0.5d0 * tol_ * tol_
         dphi_0 = -dir_deriv           ! < 0 for descent direction
         !
         ! Non-monotone reference (Oct-2026): phi_ref replaces phi_0 in
         ! the sufficient-decrease test.  With jfnk_ls_nonmono=F this
         ! collapses to the monotone Armijo test.
         !
         phi_ref = phi_0
         if(jfnk_ls_nonmono) then
            do i = 1, min(jfnk_ls_nmwin_s, size(jfnk_ls_nmhist))
               if(jfnk_ls_nmhist(i) > 0d0) &
                  phi_ref = max(phi_ref, 0.5d0 * jfnk_ls_nmhist(i)**2)
            end do
            ! Cap a stale history reference so a trial cannot claim an
            ! unbounded increase. [USE_JFNK_SOLVER] defaults to F.
            if(jfnk_del > 0d0) phi_ref = min(phi_ref, 1.21d0 * phi_0)
         end if
         !
         ls_loop: do ls_step = 1, jfnk_ls_max_steps
            alpha_ls = jfnk_ls_beta_fac ** dble(ls_step)
            jfnk_x = x_ - alpha_ls * s
            call f(n, jfnk_x, jfnk_fx)
            jfnk_total_f_evals = jfnk_total_f_evals + 1
            ! SER anchor fix (06-Oct-2026): strip to F before |F|.
            call jfnk_ptc_strip_v1(n, jfnk_x, jfnk_fx)
            jfnk_tol = dsqrt(dp(n, jfnk_fx, jfnk_fx))
            phi_alpha = 0.5d0 * jfnk_tol * jfnk_tol
            if(phi_alpha <= phi_ref + jfnk_ls_armijo_c * alpha_ls * dphi_0) then
               ! Armijo sufficient decrease satisfied
               ls_accepted = .true.
               s    = alpha_ls * s          ! scale the step
               snrm = dsqrt(dp(n, s, s))     ! update weighted norm
               s_inf = maxval(abs(s))        ! update inf norm
               snrm_tr = alpha_ls * snrm_tr  ! 2026-09-13 TR-units fix
               if(jfnk_del > 0d0 .and. jfnk_tol <= tol_) then
                  jfnk_del = min(mxdl_, snrm_tr * 2d0)
               end if
               ared = tol_ - jfnk_tol
               if(info == 1) write(luer,'(A,ES10.3,A,ES14.4,A,ES14.4)') &
                  '  JFNK: line search accepted α=', real(alpha_ls), &
                  ' ||F||=', real(jfnk_tol), ' ||s||=', real(snrm)
               jfnk_n_ls_used = jfnk_n_ls_used + 1
               exit ls_loop
            end if
         end do ls_loop
         !
         ! Audit 17-Sep-2026 (critical): a failed line search left
         ! jfnk_x/jfnk_fx/jfnk_tol holding the LAST failed Armijo trial.
         ! The Cauchy fallback and the trust-region acceptance tests
         ! below would then evaluate that stale trial state instead of
         ! the Newton point.  Restore the Newton state unconditionally.
         !
         if(.not. ls_accepted) then
            jfnk_x   = x_
            jfnk_fx  = fx_
            jfnk_tol = tol_
            if(info == 1) write(luer,'(A,I3)') &
               '  JFNK: line search failed (', jfnk_ls_max_steps, &
               ' steps) -- restored Newton state'
         end if
      end if
      !
      !================================================================
      ! Cauchy (preconditioned-residual) fallback (review §3, item 3)
      !================================================================
      ! When the Newton line search fails (or the Newton direction is
      ! not a descent direction), try d_C = M⁻¹F as a fallback.  This
      ! is a damped Jacobi / preconditioned steepest-descent step that
      ! uses the existing preconditioner and does NOT require GMRES.
      !
      ! The step is x_new = x − α·d_C where d_C = M⁻¹F.  We use a simple
      ! sufficient-decrease test ||F(x−α·d_C)|| < ||F(x)|| (no Armijo,
      ! since computing <F, J·d_C> would need an extra JVP).
      !
      if(jfnk_do_cauchy .and. .not. ls_accepted .and. ared < 0d0) then
         !
         ! Compute Cauchy direction: d_C = M⁻¹ * F(x)
         !
         ! edited by Devin, 05-Oct-2026 (Cauchy workspace fix):
         ! jfnk_alloc_workspace is never called, so jfnk_w_cauchy was
         ! only sized implicitly by Fortran 2003 realloc-on-assignment.
         ! Allocate it explicitly instead -- the implicit form also
         ! silently kept a wrong-length array if N ever shrank.
         if(allocated(jfnk_w_cauchy)) then
            if(size(jfnk_w_cauchy) /= n) deallocate(jfnk_w_cauchy)
         end if
         ! ios_nw must be reset here too: with the workspace already
         ! allocated the conditional allocate is skipped and a stale
         ! stat value would fake an alloc failure (seen under gfortran
         ! on the second setauto admission).
         ios_nw = 0
         if(.not. allocated(jfnk_w_cauchy)) &
            allocate(jfnk_w_cauchy(n), stat=ios_nw)
         if(ios_nw /= 0) then
            if(info == 1) write(luer,'(A)') &
               '  JFNK: cauchy workspace alloc failed -- skipping step'
            jfnk_exit_reason = 'cauchy_alloc_failure'
            info = 3
            return
         end if
         jfnk_w_cauchy = fx_          ! start from F(x)
         call dfp(n, jfnk_w_cauchy)   ! d_C = M⁻¹ * F
         !
         ! FIX_T_AUTO: the Cauchy direction is a step candidate too --
         ! its T components get the same damping as the Newton step.
         !
         call jfnk_tauto_damp_v1(n, jfnk_w_cauchy)
         !
         ! Audit 17-Sep-2026 (critical #2): the Cauchy direction must
         ! pass the SAME safeguards as the Newton direction -- it comes
         ! straight from the (possibly poorly conditioned) BA
         ! preconditioner, and previously bypassed the finiteness,
         ! inf-norm-cap, and runaway-norm checks entirely.
         !
         cauchy_ok = .true.
         if(jfnk_pc_clip_event) then
            ! M^{-1} clipped during this evaluation -- nonlinear operator.
            cauchy_ok = .false.
            jfnk_need_precond_refresh = 1
            if(info == 1) write(luer,'(A)') &
               '  JFNK: preconditioner clipped in Cauchy direction -- &
               &skipped'
         else if(any(jfnk_w_cauchy /= jfnk_w_cauchy) .or. &
                 any(abs(jfnk_w_cauchy) > huge(1d0))) then
            cauchy_ok = .false.
            if(info == 1) write(luer,'(A)') &
               '  JFNK: non-finite Cauchy direction -- skipped'
         end if
         if(cauchy_ok) then
            cauchy_nrm = dsqrt(dp(n, jfnk_w_cauchy, jfnk_w_cauchy))
            ! Same relative cap as the Newton step: d_C = M^{-1} F is
            ! built at the ref_norm scale by construction.
            if(cauchy_nrm > jfnk_runaway_fac * max(ref_norm, 1.0d-30)) then
               cauchy_ok = .false.
               if(info == 1) write(luer,'(A,ES12.4)') &
                  '  JFNK: runaway Cauchy norm ', real(cauchy_nrm), &
                  ' -- skipped'
            else
               ! Same direction-preserving inf-norm cap as the Newton step.
               s_inf_c = maxval(abs(jfnk_w_cauchy))
               if(s_inf_c > jfnk_step_inf_cap) then
                  jfnk_w_cauchy = jfnk_w_cauchy * &
                     (jfnk_step_inf_cap / s_inf_c)
                  cauchy_nrm    = cauchy_nrm * (jfnk_step_inf_cap / s_inf_c)
                  if(info == 1) write(luer,'(A,ES10.3)') &
                     '  JFNK: Cauchy direction scaled by ', &
                     real(jfnk_step_inf_cap / s_inf_c)
               end if
            end if
         end if
         !
         ! Backtracking line search on Cauchy direction
         ! (start at α=1 and backtrack)
         !
         if(cauchy_ok) then
         cauchy_loop: do ls_step = 0, jfnk_ls_max_steps
            if(ls_step == 0) then
               alpha_ls = 1.0d0
            else
               alpha_ls = jfnk_ls_beta_fac ** dble(ls_step)
            end if
            jfnk_x = x_ - alpha_ls * jfnk_w_cauchy
            call f(n, jfnk_x, jfnk_fx)
            jfnk_total_f_evals = jfnk_total_f_evals + 1
            ! SER anchor fix (06-Oct-2026): strip to F before |F|.
            call jfnk_ptc_strip_v1(n, jfnk_x, jfnk_fx)
            jfnk_tol = dsqrt(dp(n, jfnk_fx, jfnk_fx))
            if(jfnk_tol < tol_) then
               ! Simple sufficient decrease
               ls_accepted = .true.
               s    = alpha_ls * jfnk_w_cauchy  ! step in Newton convention
               snrm = dsqrt(dp(n, s, s))
               s_inf = maxval(abs(s))
               ! Audit 17-Sep-2026 (critical #18): the trust-region
               ! radius must be the PHYSICAL weighted step norm
               ! alpha*||d_C||_w, not alpha*||F||_w -- ||M^{-1}F||_w is
               ! generally not equal to ||F||_w.
               snrm_tr = alpha_ls * cauchy_nrm
               if(jfnk_del > 0d0) jfnk_del = min(mxdl_, snrm_tr * 2d0)
               ared = tol_ - jfnk_tol
               if(info == 1) write(luer,'(A,ES10.3,A,ES14.4)') &
                  '  JFNK: Cauchy step accepted α=', real(alpha_ls), &
                  ' ||F||=', real(jfnk_tol)
               jfnk_n_cauchy_used = jfnk_n_cauchy_used + 1
               exit cauchy_loop
            end if
         end do cauchy_loop
         end if
         !
         ! Audit 17-Sep-2026 (critical #1): a failed Cauchy search left
         ! jfnk_x/jfnk_fx/jfnk_tol holding the LAST failed trial.  The
         ! trust-region tests below would then evaluate that stale
         ! state instead of the Newton point -- restore it and report
         ! ared=0 so the rejection logic sees the true situation.
         !
         if(.not. ls_accepted) then
            jfnk_x   = x_
            jfnk_fx  = fx_
            jfnk_tol = tol_
            ared     = 0d0
            if(info == 1) write(luer,'(A)') &
               '  JFNK: Cauchy fallback also failed -- &
               &restored Newton state'
         end if
      end if
      !
      !================================================================
      ! Trust-region acceptance / rejection (existing logic)
      !================================================================
      ! If the line search or Cauchy fallback accepted a step, skip
      ! the trust-region rejection branches and go straight to the
      ! step-accepted path (ginfo = 0).
      !
      if(ls_accepted .and. ared >= 0d0) then
         ginfo = 0  ! monotone Armijo decrease
      else if(ls_accepted .and. del == 0d0) then
         ginfo = 0  ! no trust region; keep the non-monotone LS outcome
      else if(ls_accepted) then
         ! [USE_JFNK_SOLVER] defaults to F. A non-monotone increase in
         ! trust-region mode must restore the base and take the normal
         ! rejection path so it cannot grow del or skip its counter.
         jfnk_n_tr_rejected = jfnk_n_tr_rejected + 1
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         jfnk_del = 0.5d0 * snrm_tr
         ginfo = 2
      else if(del == 0d0) then
         if(info == 1) write(luer,*) '  JFNK: took full newton step'
         ginfo = 0  ! step accepted, clear hookstep flag
      else if(jfnk_tol > tol__) then
         if(info == 1) write(luer,*) '  JFNK: accepting previous step'
         jfnk_x   = x__
         jfnk_fx  = fx__
         jfnk_tol = tol__
         jfnk_del = del__
         ginfo = 0  ! step accepted, clear hookstep flag
      else if(pred_negligible) then
         ! V4 audit fix (2026-09-09): Predicted reduction is negligible.
         ! The linearized model does not predict a reduction in the
         ! physical residual norm.  Reject the step and shrink the trust
         ! radius.  Reuse the Krylov space (ginfo=2) for efficiency.
         if(info == 1) write(luer,*) '  JFNK: pred negligible, &
            &shrinking trust region'
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         ! Shrink below the rejected step (2026-09-11): halving del_try
         ! alone does not bind when del >> |s|, so the hookstep would
         ! return the identical step and loop.
         jfnk_del = 0.5d0 * min(del_try, snrm_tr)   ! 2026-09-13 TR-units fix
         ginfo   = 2
      else if(ared <= 0d0) then
         ! editted by Devin, 12-Nov-2026: a trial that fails to reduce
         ! |F| is a rejection, not a retry.  ared == 0 used to land in
         ! the rho < 0.75 branch below, which retries without touching
         ! n_reject -- OUTGEN showed it halving del from ~6 to 1e-8
         ! over hundreds of zero-progress trials until
         ! trust_region_too_small.  Count it as a rejection so
         ! no-progress ends the solve in n_reject_max passes instead of
         ! crawling the del floor; ared < 0 (norm increase) takes the
         ! same path it always did.
         if(info == 1) write(luer,*) '  JFNK: norm did not decrease, &
            &try smaller step'
         n_reject = n_reject + 1
         jfnk_n_tr_rejected = jfnk_n_tr_rejected + 1
         if(n_reject >= n_reject_max .or. jfnk_del < mndl) then
            if(info == 1) write(luer,'(A)') &
               '  JFNK: too many rejected steps -- aborting Newton'
            jfnk_x = x_
            jfnk_fx = fx_
            jfnk_tol = tol_
            info = 3
            jfnk_exit_reason = 'too_many_rejected_steps'
            return
         end if
         jfnk_x = x_
         jfnk_fx = fx_
         jfnk_tol = tol_
         jfnk_del = 0.5d0 * snrm_tr   ! 2026-09-13 TR-units fix
         ginfo   = 2
      else if(rho < 0.75d0) then
         if(info == 1) write(luer,*) '  JFNK: step ok, trying smaller &
            &step'
         x__     = jfnk_x
         fx__    = jfnk_fx
         tol__   = jfnk_tol
         ! 2026-09-13 TR-units fix: radii in hookstep (||y||) units
         if(rho > 0.1d0)  del__ = snrm_tr
         if(rho <= 0.1d0) del__ = 0.5d0 * snrm_tr
         jfnk_del = 0.7d0 * snrm_tr
         ! edited by Devin, 05-Oct-2026 (retry-pass base-point fix):
         ! restore the module state to the linearization point before
         ! the recompute pass.  JFNK_JVP_V1 differentiates F at jfnk_x
         ! against jfnk_fx; leaving the trial here made the next pass's
         ! pred-residual JVP (and any descent-retry GMRES) evaluate J at
         ! the wrong base point against b = fx_.  The trial itself is
         ! kept in x__/fx__/tol__/del__ for the "accept previous" branch.
         jfnk_x   = x_
         jfnk_fx  = fx_
         jfnk_tol = tol_
         ginfo   = 2
      else if(snrm_tr < jfnk_del * 0.9d0) then
         if(info == 1) write(luer,*) '  JFNK: step good, took full &
            &newton step'
         if(snrm < jfnk_runaway_fac * max(ref_norm, 1.0d-30)) &
            jfnk_del = min(mxdl_, snrm_tr * 2d0)   ! 2026-09-13 TR-units fix
         ginfo = 0  ! step accepted, clear hookstep flag
      else if(jfnk_del < mxdl_ * 0.9d0) then
         if(info == 1) write(luer,*) '  JFNK: step good, trying larger &
            &step'
         x__     = jfnk_x
         fx__    = jfnk_fx
         tol__   = jfnk_tol
         del__   = jfnk_del
         if(snrm < jfnk_runaway_fac * max(ref_norm, 1.0d-30)) &
            jfnk_del = min(mxdl_, snrm_tr * 2d0)   ! 2026-09-13 TR-units fix
         ! edited by Devin, 05-Oct-2026 (retry-pass base-point fix):
         ! same restore as the rho<0.75 branch above -- the module
         ! state must sit at the linearization point (x_, fx_) because
         ! JFNK_JVP_V1 differentiates F at jfnk_x against jfnk_fx.
         jfnk_x   = x_
         jfnk_fx  = fx_
         jfnk_tol = tol_
         ginfo   = 2
      end if
      !
      ! Check if need to try another s
      !
      if(ginfo == 1 .or. ginfo == 2) then
         cycle newton_loop
      end if
      !
      ! End of iteration
      !
      ! Save |F| for the next Eisenstat-Walker forcing term computation.
      ! jfnk_tol_prev is used at the top of the next Newton step to compute
      ! eta_{k+1} = gamma * (|F_{k+1}| / |F_k|)^alpha.
      jfnk_tol_prev = tol_
      !
      ! Pseudo-time (Knoll & Keyes 2004; SER update of Kelley & Keyes
      ! 1998).  Grow dtau when |F| fell, so the added term
      ! (x - x_k)*s/dtau fades and the iteration becomes Newton.
      ! Shrink it when |F| did not fall.  jfnk_tol/jfnk_tol_prev are
      ! |F| -- the strip above removed the pseudo-time term from every
      ! residual this routine reads.  One line per accepted iteration;
      ! the sweep chatter stays in JFNK_EVAL.log.
      !
      if(jfnk_ptc_active) then
         dtau_was = jfnk_ptc_dtau
         if(jfnk_tol < jfnk_tol_prev) then
            jfnk_ptc_dtau = min(jfnk_ptc_dtau_max, &
               jfnk_ptc_dtau * jfnk_tol_prev / max(jfnk_tol, 1d-30))
         else
            jfnk_ptc_dtau = max(1d-2, 0.5d0 * jfnk_ptc_dtau)
         end if
         if(info == 1) then
            write(luer,'(A)') ' '
            write(luer,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
               '  JFNK: dtau ', real(dtau_was), ' -> ', &
               real(jfnk_ptc_dtau), &
               '    |F| ', real(jfnk_tol_prev), ' -> ', real(jfnk_tol)
         end if
      end if
      !
      ! Push the accepted |F| into the non-monotone LS history ring
      ! (Oct-2026, Li-Fukushima reference set).
      !
      jfnk_ls_nmidx = mod(jfnk_ls_nmidx, size(jfnk_ls_nmhist)) + 1
      jfnk_ls_nmhist(jfnk_ls_nmidx) = jfnk_tol
      !
      jfnk_nits = jfnk_nits + 1
      n_reject = 0
      call sub()
      x_    = jfnk_x
      fx_   = jfnk_fx
      tol_  = jfnk_tol
      tol__ = 1d99
      !
      ! edited by Devin, 06-Oct-2026 (SER anchor fix): the accepted
      ! iterate is the new anchor x_k for the pseudo-time term, so the
      ! next GMRES solves (J_F + D/dtau) s = F(x_k) with the
      ! linearization-point residual equal to F.
      call jfnk_ptc_set_anchor_v1(n, x_)
      !
      ! A coupled residual solve can carry sweep noise, so maxch alone
      ! cannot certify its root. [JFNK_COUPLED] defaults to F.
      !
      maxch_converged = jfnk_maxch_tol > 0d0 .and. &
         jfnk_maxch_est >= 0d0 .and. jfnk_maxch_est < jfnk_maxch_tol
      if(maxch_converged .and. jfnk_tol >= tol .and. info == 1) then
         if(JFNK_COUPLED_LOC) then
            write(luer,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
               '  JFNK: maxch-only convergence blocked in coupled mode; maxch_est=', &
               real(jfnk_maxch_est), ' |F|=', &
               real(jfnk_tol), ' |F|/|F0|=', &
               real(jfnk_tol/max(tol_initial,1d-30))
         else
            write(luer,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
               '  JFNK: maxch-only convergence; maxch_est=', &
               real(jfnk_maxch_est), ' |F|=', real(jfnk_tol), &
               ' |F|/|F0|=', real(jfnk_tol/max(tol_initial,1d-30))
         end if
      end if
      if(jfnk_tol < tol .or. &
         (maxch_converged .and. .not. JFNK_COUPLED_LOC)) then
         if(info == 1) write(luer,*) '  JFNK: NONLINEAR CONVERGED', &
            ' (res=', real(jfnk_tol), ' maxch_est=', real(jfnk_maxch_est), ')'
         info = 0
         jfnk_exit_reason = 'converged'
            return
      else if(jfnk_nits == nits) then
         if(info == 1) write(luer,*) '  JFNK: MAXIMUM NONLINEAR ITERATIONS REACHED'
         if(info == 1) write(luer,*) '  JFNK: NONLINEAR SOLVE NOT YET CONVERGED'
         info = 2
         jfnk_exit_reason = 'max_newton_iterations'
            return
      end if
      !
   end do newton_loop
   !
end subroutine newton_hook_v1


!----------------------------------------------------------------------
! GMRES_M_V1: GMRES(m) with left or right preconditioning and hookstep.
!
! Solves A*x = b for x, subject to |x| < delta (hookstep).
!
! Right preconditioning (legacy): A*M^{-1}*u = b, x = M^{-1}*u
!   - Initial residual: r = b - A*x  (not preconditioned)
!   - Arnoldi: z = M^{-1}*v_j, w = A*z
!   - Update: x += M^{-1}*V*y
!
! Left preconditioning (Arramy 2024 default, Issue 6):
!   M^{-1}*A*x = M^{-1}*b
!   - Initial residual: r = M^{-1}*(b - A*x)
!   - Arnoldi: w = A*v_j, w = M^{-1}*w
!   - Update: x += V*y  (no final psolve)
!
! Bug fix (09-Sep-2026): del (trust radius) is now INTENT(IN) and the
! predicted residual norm is returned separately as pred_res.  Previously
! del was INTENT(INOUT) and was overwritten by hook_step_v1 with the
! predicted residual, which:
!   (a) corrupted the trust radius on GMRES restart (goto 1)
!   (b) gave wrong pred = tol_ - 0 = tol_ when del==0 (no trust region)
!       instead of pred = tol_ - res (the GMRES residual)
!
! editted by Sean Brennan (08-Sep-2026)
!----------------------------------------------------------------------
subroutine gmres_m_v1(m, n, x, b, matvec, psolve, dotprd, h, v, &
                       res, del, pred_res, its, info, luer)
   use, intrinsic :: ieee_arithmetic
   use jfnk_solver_v1_mod, only: jfnk_precond_fails, &
                                  jfnk_total_jvp_evals, jfnk_use_left_pc, &
                                  jfnk_last_ynrm, jfnk_jvp_failed, &
                                  jfnk_pc_clip_event, jfnk_gmres_stag_tol, &
                                  jfnk_gmres_recomputed, &
                                  jfnk_gmres_jj, jfnk_gmres_beta, &
                                  jfnk_gmres_z, jfnk_gmres_g, &
                                  jfnk_gmres_x_cycle_start, &
                                  jfnk_gmres_max_restart, jfnk_gmres_tres_fac, &
                                  jfnk_gmres_true_res, jfnk_gmres_true_res_ok, &
                                  jfnk_gmres_rprev, jfnk_exit_reason, &
                                  jfnk_lgmres_on, jfnk_lgmres_kmax, &
                                  jfnk_lgmres_n, jfnk_lgmres_d, &
                                  jfnk_lgmres_ad, jfnk_lgmres_dlast, &
                                  jfnk_lgmres_q, jfnk_gdefl_on, &
                                  jfnk_gdefl_kmax, jfnk_aug_kmax, &
                                  jfnk_recyc_on, jfnk_cycstall_on, &
                                  jfnk_cycstall_frac, jfnk_cycstall_nmax
   ! Fix #7 (review 2026-09-11): jfnk_precond_calls removed from GMRES.
   ! JFNK_PRECOND_V1 owns the sole preconditioner application counter
   ! (JFNK_N_PRECOND_APPLIES).  Incrementing here double-counted every
   ! psolve call.
   implicit none
   integer,          intent(in)    :: m, n
   double precision, intent(inout) :: x(n)
   double precision, intent(in)    :: b(n)
   external                        :: matvec, psolve
   double precision, external      :: dotprd
   double precision, intent(inout) :: h(m+1, m)
   double precision, intent(inout) :: v(n, m+1)
   ! res is dual-use: on entry it carries the caller's weighted
   ! residual target (gres = eta*ref_norm), captured into tol below;
   ! on exit it returns the final Hessenberg LS residual of the last
   ! cycle.  newton_hook_v1 reads pred_res, not res, for its logic.
   double precision, intent(inout) :: res
   double precision, intent(in)    :: del
   double precision, intent(out)   :: pred_res
   integer,          intent(inout) :: its, info
   integer,          intent(in)    :: luer
   !
   double precision :: tol
   double precision :: res_best   ! 2026-09-13 TR-units fix: best LS residual this restart cycle
   integer :: n_stag              ! 2026-09-13 TR-units fix: consecutive non-improving residuals
   double precision :: h_(m+1, m), y(m+1), p(m+1), work(6*m+16)
   double precision :: s_tmp  ! re-orthogonalization scratch
   ! LGMRES work (23-Oct-2026): hard bound KMAX=8 matches the VADAT
   ! clamp on JFNK_LGMRES_K.  The augmented LS is (j+1+k) x (j+k) --
   ! dense, because the appended image columns are full -- so it goes
   ! through DGELSY like the Hessenberg system, and the hookstep takes
   ! an augmented Gram built once per cycle end.
   integer, parameter :: lg_kmax = 8
   double precision :: h_aug(m+1+lg_kmax, m+lg_kmax)
   double precision :: hd_aug(m+1+lg_kmax, m+lg_kmax)
   double precision :: p_aug(m+1+lg_kmax), y_aug(m+lg_kmax)
   double precision :: g_aug(m+lg_kmax, m+lg_kmax)
   double precision :: work_aug(6*(m+lg_kmax)+16)
   integer :: piv_aug(m+lg_kmax)
   integer :: lg_k, lg_i, lg_l, lg_nv, lg_rank
   double precision :: lg_rho, lg_scale
   logical :: lg_active, lg_have_last, gd_active
   ! aug_cap: actual slot capacity this solve (LGMRES + GDEFL combined,
   ! bounded by the lg_kmax workspace bound).  jvp_soft: a JVP failure
   ! killed a Krylov column mid-cycle -- the cycle closes on the valid
   ! basis and exits through the soft-exit path.  res_prev_end /
   ! n_cycstall back the [JFNK_CYCSTALL] per-cycle contraction guard.
   integer :: aug_cap
   logical :: jvp_soft
   double precision :: res_prev_end
   integer :: n_cycstall
   double precision :: true_res  ! V4 audit: true residual at restart
   double precision :: del_cycle, x_cycle_norm
   integer :: imx, piv(m), rank, i, j, ios_gm
   ! beta and jj are no longer SAVE variables: the cross-call basis
   ! state the info==2 hookstep recompute needs lives in the module
   ! pair jfnk_gmres_jj / jfnk_gmres_beta (which external diagnostic
   ! callers can save/restore), synced at the cycle checkpoint below.
   ! (numerical audit fix, 2026-10)
   double precision :: beta
   integer :: jj
   logical :: done, left_pc
   ! ltres: a JFNK_TRES_FAC true-residual measurement was taken this
   ! cycle -- the restart-diag block reuses it instead of spending a
   ! second eval on the same vector.
   logical :: ltres, tres_rej
   ! Lucky-breakdown cycle end: h(j+1,j)~0 means the Krylov space is
   ! invariant, so every later Arnoldi step is built on the zeroed
   ! v(:,j+1) and can only spin until n_stag aborts the solve.  End the
   ! cycle, take the LSQ step, and let the j==m restart path re-seed
   ! the basis from the true residual.
   logical :: lb_hit
   integer :: n_restart
   !
   !
   ! Heap-allocated workspace (audit §3, 09-Sep-2026):
   ! w, z were automatic arrays of size N (~160KB for N=10000).
   !
   double precision, allocatable, save :: w(:), z(:)
   ! edited by Devin, 05-Oct-2026 (alloc-failure stale-output fix):
   ! initialize the outputs before the workspace allocation below -- an
   ! alloc failure returns early, and newton_hook_v1 reads pred_res
   ! (stagnation gate) and jfnk_last_ynrm (snrm_tr) unconditionally.
   pred_res = 0d0
   jfnk_last_ynrm = 0d0
   ! Alloc-failure returns below precede the per-exit resets; without
   ! this a stale measured residual from a previous solve could pass
   ! the caller's true-residual stagnation gate.
   jfnk_gmres_true_res_ok = .false.
   jfnk_gmres_true_res = 0d0
   if(allocated(w)) then
      if(size(w) .ne. n) deallocate(w, z)
   end if
   if(.not. allocated(w)) then
      allocate(w(n), z(n), stat=ios_gm)
      if(ios_gm /= 0) then
         if(info == 1) write(luer,*) 'gmres_m_v1: workspace alloc failed'
         info = 1
         jfnk_exit_reason = 'gmres_alloc_fail'
         return
      end if
   end if
   if(allocated(jfnk_gmres_x_cycle_start)) then
      if(size(jfnk_gmres_x_cycle_start) /= n) &
         deallocate(jfnk_gmres_x_cycle_start)
   end if
   if(.not. allocated(jfnk_gmres_x_cycle_start)) then
      allocate(jfnk_gmres_x_cycle_start(n))
      jfnk_gmres_x_cycle_start = x
   end if
   !
   left_pc = jfnk_use_left_pc
   ! The correction vectors live in x-space under both conventions;
   ! only the image map differs (A z vs M^{-1} A z) -- handled at pair
   ! formation below, one extra psolve per restart for left-PC.
   ! Augmentation runs when either retention scheme is armed:
   ! [JFNK_LGMRES] fills slots from completed-cycle steps, [JFNK_GDEFL]
   ! from harmonic Ritz extraction at each restart.  The capacity is
   ! the driver-set combined bound, rebuilt here when a direct caller
   ! (unit tests) left jfnk_aug_kmax at its zero default.
   lg_active = jfnk_lgmres_on .or. jfnk_gdefl_on
   gd_active = jfnk_gdefl_on
   aug_cap = jfnk_aug_kmax
   if(jfnk_lgmres_on .and. .not. jfnk_gdefl_on) &
      aug_cap = max(aug_cap, jfnk_lgmres_kmax)
   if(jfnk_gdefl_on .and. .not. jfnk_lgmres_on) &
      aug_cap = max(aug_cap, jfnk_gdefl_kmax)
   if(jfnk_lgmres_on .and. jfnk_gdefl_on) &
      aug_cap = max(aug_cap, jfnk_lgmres_kmax + jfnk_gdefl_kmax)
   aug_cap = min(aug_cap, lg_kmax)
   if(aug_cap < 1) lg_active = .false.
   lg_k = 0
   pred_res = 0d0
   !
   ! Right-PC step map (2026-10 TR-units fix): Z = M^{-1} V and its
   ! Gram G = Z' W Z are stored so the hookstep constrains ||s||_w
   ! and s = Z y is recovered without a second (possibly nonlinear)
   ! psolve application.
   !
   if(.not. left_pc) then
      if(allocated(jfnk_gmres_z)) then
         done = allocated(jfnk_gmres_g)
         if(done) done = (size(jfnk_gmres_g,1) >= m)
         if(size(jfnk_gmres_z,1) /= n .or. size(jfnk_gmres_z,2) < m &
            .or. .not. done) then
            if(allocated(jfnk_gmres_g)) deallocate(jfnk_gmres_g)
            deallocate(jfnk_gmres_z)
         end if
      end if
      if(.not. allocated(jfnk_gmres_z)) then
         allocate(jfnk_gmres_z(n,m), jfnk_gmres_g(m,m), stat=ios_gm)
         if(ios_gm /= 0) then
            if(info == 1) write(luer,*) 'gmres_m_v1: z/g alloc failed'
            info = 1
            jfnk_exit_reason = 'gmres_alloc_fail_zg'
            return
         end if
      end if
   end if
   !
   ! Augmentation storage is independent of the preconditioning side.
   if(lg_active) then
      if(allocated(jfnk_lgmres_d)) then
         if(size(jfnk_lgmres_d,1) /= n .or. &
            size(jfnk_lgmres_d,2) /= aug_cap) then
            deallocate(jfnk_lgmres_d, jfnk_lgmres_ad, &
                       jfnk_lgmres_dlast, jfnk_lgmres_q)
         end if
      end if
      if(.not. allocated(jfnk_lgmres_d)) then
         allocate(jfnk_lgmres_d(n,aug_cap), &
                  jfnk_lgmres_ad(n,aug_cap), &
                  jfnk_lgmres_dlast(n), &
                  jfnk_lgmres_q(n,aug_cap), stat=ios_gm)
         if(ios_gm /= 0) then
            if(info == 1) write(luer,*) 'gmres_m_v1: lgmres alloc'// &
               ' failed -- running plain GMRES'
            lg_active = .false.
         end if
      end if
      if(lg_active) then
         ! [JFNK_RECYC] keeps the retained set across solves; stale
         ! images are still useful augmentation constraints.
         if(.not. jfnk_recyc_on) jfnk_lgmres_n = 0
         jfnk_lgmres_n = min(jfnk_lgmres_n, aug_cap)
      end if
   end if
   !
   ! Hookstep recalculation with new del
   !
   ! editted by Devin, 04-Oct-2026: only valid when a completed Arnoldi
   ! cycle left a basis (jj > 0).  Early returns that never reach the
   ! jj assignment -- DGELSY error, near-zero beta, JVP/clip abort --
   ! used to reach this branch with a stale or uninitialized jj and
   ! index h(1:j+1,1:j) over a Krylov space that does not belong to
   ! this point; jj = 0 on every fresh-solve entry invalidates it.
   if(info == 2 .and. jfnk_gmres_jj > 0 .and. &
      jfnk_gmres_jj <= m) then
      jj   = jfnk_gmres_jj
      beta = jfnk_gmres_beta
      if(left_pc) then
         x_cycle_norm = dsqrt(max(0d0, &
              dotprd(n, jfnk_gmres_x_cycle_start, &
                     jfnk_gmres_x_cycle_start)))
         del_cycle = del
         if(del_cycle > 1d-30) then
            call hook_step_v1(jj, h, m, beta, del_cycle, y, pred_res)
         else
            y = 0d0
            pred_res = beta
         end if
         z = matmul(v(:,1:jj), y(1:jj))
      else
         ! Right-PC: constrain ||s||_w = sqrt(y' G y), recover s = Z y
         ! from the stored preconditioned basis.
         x_cycle_norm = dsqrt(max(0d0, &
              dotprd(n, jfnk_gmres_x_cycle_start, &
                     jfnk_gmres_x_cycle_start)))
         del_cycle = del
         if(del_cycle > 1d-30) then
            call hook_step_g_v1(jj, h, m, beta, del_cycle, &
                 jfnk_gmres_g, size(jfnk_gmres_g,1), y, pred_res)
         else
            y = 0d0
            pred_res = beta
         end if
         jfnk_last_ynrm = dsqrt(max(0d0, dot_product(y(1:jj), &
              matmul(jfnk_gmres_g(1:jj,1:jj), y(1:jj)))))
         z = matmul(jfnk_gmres_z(:,1:jj), y(1:jj))
      end if
      x = jfnk_gmres_x_cycle_start + z
      if(del > 0d0) then
         x_cycle_norm = dsqrt(max(0d0, dotprd(n, x, x)))
         if(x_cycle_norm > del) x = x * (del / x_cycle_norm)
      end if
      jfnk_last_ynrm = dsqrt(max(0d0, dotprd(n, x, x)))
      info = 0
      ! edited by Devin, 05-Oct-2026 (recompute-pass descent-skip fix):
      ! tell newton_hook_v1 this step came from a recalculation on the
      ! unchanged Krylov basis, so its descent check skips the pass --
      ! info=0 above hides that fact from the caller.
      jfnk_gmres_recomputed = .true.
      jfnk_gmres_true_res_ok = .false.
      jfnk_exit_reason = 'gmres_recompute_pass'
      return
   end if
   jj = 0
   jfnk_gmres_jj = 0   ! fresh solve: no basis survives for recompute
   if(allocated(jfnk_gmres_rprev)) deallocate(jfnk_gmres_rprev)   ! fresh solve: no cross-solve residual history
   ! [JFNK_RECYC] carries the augmentation set over from the previous
   ! solve; otherwise each solve starts with an empty set.
   if(.not. jfnk_recyc_on) jfnk_lgmres_n = 0
   lg_have_last = .false.   ! no completed cycle yet to correct from
   jvp_soft = .false.
   res_prev_end = -1d0
   n_cycstall = 0
   !
   tol = res
   imx = its
   its = 0
   v   = 0d0
   jfnk_last_ynrm = 0d0   ! 2026-09-13 TR-units fix: fresh GMRES solve
   jfnk_gmres_true_res_ok = .false.   ! no exit measurement yet
   n_restart = 0
   !
1 continue
   n_stag = 0
   res_best = 1d99
   lb_hit = .false.
   ltres = .false.
   jfnk_gmres_x_cycle_start = x
   del_cycle = del
   !
   ! Initial residual: r = b - A*x  (then left-precondition if requested)
   !
   beta = dsqrt(dotprd(n, x, x))
   if(beta == 0d0)  w = 0d0
   if(beta /= 0d0) then
      call matvec(n, x, w)
      jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
      ! Audit 17-Sep-2026: abort the Krylov solve on an explicit JVP
      ! failure instead of interpreting a zeroed JV as zero action.
      ! Softened (gdefl work): on a restart the accumulated x already
      ! holds the completed cycles' step -- return it through the
      ! inexact-direction gate instead of discarding the solve.  A
      ! failure on the very first residual still aborts (x is the
      ! input guess; nothing to salvage).
      if(jfnk_jvp_failed) then
         if(its > 0) then
            if(info == 1) write(luer,*) '  GMRES: JVP failure at'// &
               ' restart residual -- returning partial step'
            info = 1
            jfnk_exit_reason = 'gmres_jvp_soft_exit'
            return
         end if
         if(info == 1) write(luer,*) '  GMRES: JVP failure -- aborting'
         info = 3
         jfnk_exit_reason = 'gmres_jvp_failure'
         return
      end if
   end if
   w = b - w
   !
   ! LGMRES diagnosis (23-Oct-2026): the physical residual at cycle
   ! start is what the next Arnoldi basis is built on.  If successive
   ! cycles restart on nearly the same direction (cos ~ +-1), the
   ! Baker-Jessup-Manteuffel augmentation mechanism is the relevant
   ! cure; a small cosine means each cycle explores fresh space and
   ! retained correction vectors would not repair the stall.
   !
   if(allocated(jfnk_gmres_rprev)) then
      if(size(jfnk_gmres_rprev) /= n) deallocate(jfnk_gmres_rprev)
   end if
   if(.not. allocated(jfnk_gmres_rprev)) allocate(jfnk_gmres_rprev(n))
   if(n_restart > 0 .and. luer > 0) then
      s_tmp = dotprd(n, w, jfnk_gmres_rprev)
      write(luer,'(A,I3,A,ES12.4)') &
         '  GMRES RESTART ANGLE: cyc=', n_restart, &
         '  |<r_c,r_prev>|/(|r||r_prev|)=', &
         real(abs(s_tmp) / max(dsqrt(dotprd(n, w, w)) * &
              dsqrt(dotprd(n, jfnk_gmres_rprev, jfnk_gmres_rprev)), &
              1d-300))
   end if
   ! [JFNK_GDEFL] deflated restart (Morgan 2002, GMRES-DR): extract
   ! the smallest-|theta| harmonic Ritz directions of the completed
   ! cycle into the augmentation slots.  Must run BEFORE the LGMRES
   ! append below (its vectors occupy the leading slots) and before
   ! v(:,1)/h are overwritten by the new cycle -- the completed
   ! basis and (jj+1) x jj Hessenberg block are still intact here.
   ! Extraction costs one jj x jj DGEEV; no JVPs.
   !
   if(gd_active .and. n_restart > 0 .and. jj >= 2 .and. aug_cap > 0) then
      call jfnk_gdefl_extract_v1(n, m, jj, h, m+1, v, n, dotprd, &
                                 left_pc, jfnk_lgmres_on, luer)
   end if
   !
   ! LGMRES: pair the completed cycle's step d = x_e - x_s (saved at
   ! the update below) with its operator image, recovered from the
   ! residual difference r_s - r_e = rprev - w -- no extra matvec.
   ! rprev is still the previous cycle's start residual here.  A
   ! non-finite or empty pair is dropped; at capacity the oldest
   ! vector leaves (FIFO, the recent corrections carry the direction
   ! GMRES is circling).  Under left PC the Arnoldi operator is
   ! B = M^{-1}A, so the stored image needs one psolve on the true
   ! residual difference (cheap vs a JVP; still before the r0 psolve).
   if(jfnk_lgmres_on .and. lg_have_last .and. n_restart > 0) then
      if(all(ieee_is_finite(jfnk_lgmres_dlast)) .and. &
         all(ieee_is_finite(jfnk_gmres_rprev - w))) then
         if(jfnk_lgmres_n == aug_cap) then
            jfnk_lgmres_d(:,1:aug_cap-1) = &
               jfnk_lgmres_d(:,2:aug_cap)
            jfnk_lgmres_ad(:,1:aug_cap-1) = &
               jfnk_lgmres_ad(:,2:aug_cap)
            jfnk_lgmres_n = aug_cap - 1
         end if
         jfnk_lgmres_n = jfnk_lgmres_n + 1
         jfnk_lgmres_d(:,jfnk_lgmres_n) = jfnk_lgmres_dlast
         if(left_pc) then
            z = jfnk_gmres_rprev - w
            call psolve(n, z)
            if(jfnk_pc_clip_event) then
               if(info == 1) write(luer,*) &
                  '  GMRES: precond clip -- aborting'
               info = 3
               jfnk_exit_reason = 'gmres_pc_clip'
               return
            end if
            jfnk_lgmres_ad(:,jfnk_lgmres_n) = z
         else
            jfnk_lgmres_ad(:,jfnk_lgmres_n) = jfnk_gmres_rprev - w
         end if
      end if
      lg_have_last = .false.
   end if
   jfnk_gmres_rprev = w
   !
   if(left_pc) then
      call psolve(n, w)
      ! Audit 17-Sep-2026: a clip during psolve makes M^{-1} nonlinear.
      if(jfnk_pc_clip_event) then
         if(info == 1) write(luer,*) '  GMRES: precond clip -- aborting'
         info = 3
         jfnk_exit_reason = 'gmres_pc_clip'
         return
      end if
   end if
   !
   beta = dsqrt(dotprd(n, w, w))
   if(beta < 1d-30) then
      ! Bug fix (09-Sep-2026): beta near zero means the current x (or s)
      ! already solves the system — that's convergence, not a
      ! preconditioner failure.  Previously this incremented
      ! jfnk_precond_fails, polluting diagnostics every time GMRES got
      ! a trivially-easy solve.
      if(info == 1) write(luer,*) &
         '  GMRES: near-zero initial residual (converged)'
      info = 0
      pred_res = 0d0
      jfnk_exit_reason = 'gmres_near_zero_rhs'
      return
   end if
   v(:,1) = w / beta
   !
   h = 0d0
   do j = 1, m
      its = its + 1
      !
      if(left_pc) then
         ! Left: w = M^{-1} * A * v_j
         call matvec(n, v(:,j), w)
         jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
         ! JVP soft-fail: only this direction is unevaluable (budget
         ! sentinel, clip floor, non-finite perturbed eval); the
         ! accumulated step is still valid.  The failed matvec left
         ! w = JV = 0, so the MGS below leaves column j dead and
         ! DGELSY ignores it -- the cycle closes on the completed
         ! basis and exits through the soft-exit path rather than
         ! discarding the solve.
         if(jfnk_jvp_failed) then
            if(its == 1) then
               ! First Krylov column is unevaluable: no completed
               ! basis and x is still the input guess, so there is
               ! nothing to salvage -- keep the hard fail.
               if(info == 1) write(luer,*) &
                  '  GMRES: JVP failure -- aborting'
               info = 3
               jfnk_exit_reason = 'gmres_jvp_failure'
               return
            end if
            if(info == 1) write(luer,*) '  GMRES: JVP failure --'// &
               ' closing cycle on completed basis'
            jvp_soft = .true.
         end if
         call psolve(n, w)
         if(jfnk_pc_clip_event) then
            if(info == 1) write(luer,*) '  GMRES: precond clip -- aborting'
            info = 3
            jfnk_exit_reason = 'gmres_pc_clip'
            return
         end if
      else
         ! Right: w = A * M^{-1} * v_j.  Keep z_j = M^{-1} v_j and the
         ! j-th Gram column G(:,j) = Z' W z_j -- the hookstep bounds
         ! ||s||_w = sqrt(y' G y) and the step is s = Z y, so no second
         ! psolve is needed at step formation.
         z = v(:,j)
         call psolve(n, z)
         if(jfnk_pc_clip_event) then
            if(info == 1) write(luer,*) '  GMRES: precond clip -- aborting'
            info = 3
            jfnk_exit_reason = 'gmres_pc_clip'
            return
         end if
         jfnk_gmres_z(:,j) = z
         do i = 1, j
            jfnk_gmres_g(i,j) = dotprd(n, jfnk_gmres_z(:,i), z)
            jfnk_gmres_g(j,i) = jfnk_gmres_g(i,j)
         end do
         call matvec(n, z, w)
         jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
         ! Same soft-fail as the left-PC site: the failed matvec left
         ! w = 0, so column j goes dead in the least-squares and the
         ! cycle closes on the completed basis.  z_j and its Gram
         ! column were stored before the matvec and stay valid.
         if(jfnk_jvp_failed) then
            if(its == 1) then
               ! First Krylov column is unevaluable: no completed
               ! basis and x is still the input guess, so there is
               ! nothing to salvage -- keep the hard fail.
               if(info == 1) write(luer,*) &
                  '  GMRES: JVP failure -- aborting'
               info = 3
               jfnk_exit_reason = 'gmres_jvp_failure'
               return
            end if
            if(info == 1) write(luer,*) '  GMRES: JVP failure --'// &
               ' closing cycle on completed basis'
            jvp_soft = .true.
         end if
      end if
      !
      ! Modified Gram-Schmidt with re-orthogonalization (twice MGS)
      ! Re-orthogonalization prevents loss of orthogonality in the
      ! Krylov basis for ill-conditioned problems, which causes the
      ! least-squares residual to underestimate the true residual
      ! and leads to poor GMRES restart behaviour.
      !
      do i = 1, j
         h(i,j) = dotprd(n, w, v(1,i))
         w = w - h(i,j) * v(:,i)
      end do
      ! Second MGS pass (re-orthogonalization)
      do i = 1, j
         s_tmp = dotprd(n, w, v(1,i))
         h(i,j) = h(i,j) + s_tmp
         w = w - s_tmp * v(:,i)
      end do
      h(j+1,j) = dsqrt(dotprd(n, w, w))
      !
      ! Lucky-breakdown guard (09-Sep-2026): if h(j+1,j) is numerically
      ! zero, the exact solution already lies in the current Krylov
      ! subspace (a "happy breakdown").  Dividing by zero injects NaNs
      ! into v(:,j+1), silently corrupting the rest of the solve.
      ! Fix: skip the divide, set v(:,j+1)=0, and let the least-squares
      ! residual (res) naturally collapse to ~0, triggering done.
      !
      if(h(j+1,j) < 1d-14 * beta) then
         if(info == 1) write(luer,'(A,I3,A,ES12.4)') &
            '  GMRES: lucky breakdown at it=', its, &
            ' h(j+1,j)=', real(h(j+1,j)), ' -- cycle end'
         v(:,j+1) = 0d0
         lb_hit = .true.
      else
         v(:,j+1) = w / h(j+1,j)
      end if
      !
      ! Solve least-squares for y
      !
      if(lg_active .and. jfnk_lgmres_n > 0) then
         !
         ! LGMRES (Baker-Jessup-Manteuffel 2005): append the retained
         ! correction images A z_i to the LS as dense columns.  The
         ! basis is [V_{j+1}, q_1..q_k] with q_i the orthonormalized
         ! complement of A z_i; column j+i holds [V' A z_i ; R_i]
         ! where R is the (lower-triangular) q-expansion.  A z_i is
         ! stored from the restart pairing, so no matvec is needed.
         ! The complement must be re-orthogonalized each iteration
         ! because V_{j+1} grows -- n-dot cost is trivial vs a JVP.
         ! A dependent image (already inside the span) zeroes its
         ! column: DGELSY's rank-revealing solve then keeps its
         ! coefficient ~0 and the step ignores that direction.
         !
         lg_k = jfnk_lgmres_n
         h_aug = 0d0
         h_aug(1:j+1,1:j) = h(1:j+1,1:j)
         do lg_i = 1, lg_k
            do lg_l = 1, j+1
               h_aug(lg_l, j+lg_i) = &
                  dotprd(n, v(1,lg_l), jfnk_lgmres_ad(1,lg_i))
            end do
            jfnk_lgmres_q(:,lg_i) = jfnk_lgmres_ad(:,lg_i)
            do lg_l = 1, j+1
               jfnk_lgmres_q(:,lg_i) = jfnk_lgmres_q(:,lg_i) - &
                  h_aug(lg_l, j+lg_i) * v(:,lg_l)
            end do
            do lg_l = 1, lg_i-1
               s_tmp = dotprd(n, jfnk_lgmres_q(1,lg_i), &
                              jfnk_lgmres_q(1,lg_l))
               h_aug(j+1+lg_l, j+lg_i) = s_tmp
               jfnk_lgmres_q(:,lg_i) = jfnk_lgmres_q(:,lg_i) - &
                  s_tmp * jfnk_lgmres_q(:,lg_l)
            end do
            lg_scale = dsqrt(max(0d0, dotprd(n, jfnk_lgmres_ad(1,lg_i),&
                                             jfnk_lgmres_ad(1,lg_i))))
            lg_rho = dsqrt(max(0d0, dotprd(n, jfnk_lgmres_q(1,lg_i), &
                                           jfnk_lgmres_q(1,lg_i))))
            if(lg_rho > 1d-8 * lg_scale) then
               jfnk_lgmres_q(:,lg_i) = jfnk_lgmres_q(:,lg_i) / lg_rho
               h_aug(j+1+lg_i, j+lg_i) = lg_rho
            else
               h_aug(:, j+lg_i) = 0d0
               jfnk_lgmres_q(:,lg_i) = 0d0
            end if
         end do
         p_aug = 0d0
         p_aug(1) = beta
         piv_aug = 0
         ! DGELSY overwrites A with its QR factors, so it gets a copy;
         ! h_aug must stay intact for the residual reconstruction and
         ! the augmented hookstep below (same convention as h_/h).
         hd_aug = h_aug
         call dgelsy(j+1+lg_k, j+lg_k, 1, hd_aug, m+1+lg_kmax, p_aug, &
                     m+1+lg_kmax, piv_aug, 1.0d-14, lg_rank, &
                     work_aug, 6*(m+lg_kmax)+16, i)
         if(i /= 0) then
            if(info == 1) write(luer,*) '  GMRES: dgelsy error ', i
            info = 1
            jfnk_exit_reason = 'gmres_dgelsy_error'
            return
         end if
         y_aug(1:j+lg_k) = p_aug(1:j+lg_k)
         y(1:j) = y_aug(1:j)
         p_aug(1:j+1+lg_k) = &
            -matmul(h_aug(1:j+1+lg_k,1:j+lg_k), y_aug(1:j+lg_k))
         p_aug(1) = p_aug(1) + beta
         res = dsqrt(dot_product(p_aug(1:j+1+lg_k), &
                                 p_aug(1:j+1+lg_k)))
      else
         p(1) = beta
         p(2:j+1) = 0d0
         h_(1:j+1,1:j) = h(1:j+1,1:j)
         piv = 0
         !
         ! DGELSY workspace (ChatGPT audit fix #1, 2026-09-09): the
         ! least-squares solve for the GMRES y-coefficients needs ~5*M
         ! workspace for this rank-deficient Hessenberg problem.  The
         ! original 4*M+1 was undersized and could overwrite adjacent
         ! stack memory inside GMRES, producing otherwise-inexplicable
         ! crashes around the Newton-loop return/deallocation.  Use a
         ! generous 6*M+16 (matches the declaration above).
         !
         call dgelsy(j+1, j, 1, h_, m+1, p, m+1, piv, 1.0d-14, rank, work, &
                      6*m+16, i)
         if(i /= 0) then
            if(info == 1) write(luer,*) '  GMRES: dgelsy error ', i
            info = 1
            jfnk_exit_reason = 'gmres_dgelsy_error'
            return
         end if
         y = p
         !
         ! Residual of the (preconditioned) least-squares problem
         !
         p(1:j+1) = -matmul(h(1:j+1,1:j), y(1:j))
         p(1) = p(1) + beta
         res = dsqrt(dot_product(p(1:j+1), p(1:j+1)))
      end if
      ! 2026-09-13 TR-units fix: tolerate small FD-noise residual rises.
      ! The old `res > res_` test aborted the whole Newton solve on any
      ! non-monotone wiggle; now require 3 consecutive non-improvements.
      ! Audit 17-Sep-2026: the improvement threshold is jfnk_gmres_stag_tol
      ! (default 1D-3), tied to the finite-difference JVP noise floor --
      ! the previous 1D-14 value let GMRES grind out sub-noise gains for
      ! many iterations before declaring stagnation.
      if(res > res_best*(1d0 - jfnk_gmres_stag_tol)) then
         n_stag = n_stag + 1
      else
         n_stag = 0
         res_best = min(res_best, res)
      end if
      done = (res <= tol .or. its == imx .or. n_stag >= 3)
      ! Throttled: print first it, every 5th, and on completion —
      ! not every GMRES iteration (was flooding OUTGEN). 2026-09-11
      if(info == 1 .and. (its == 1 .or. mod(its,5) == 0 .or. done)) &
         write(luer,'(A,I3,A,ES12.4)') &
         '    GMRES: it=', its, ' res=', real(res)
      !
      if(done .or. j == m .or. lb_hit) then
         if(del > 0d0) then
            if(lg_active .and. jfnk_lgmres_n > 0) then
               ! LGMRES hookstep: the step is V y + D c (left-PC) or
               ! Z y + D c (right-PC), so the ||s||_w constraint needs
               ! the Gram of the combined basis.  Under left PC the
               ! basis is V itself (W-orthonormal, Gram = I); under
               ! right PC it is the stored Z'WZ.  Cross and aug-aug
               ! blocks are built once per cycle end -- k*j dotprd
               ! calls, trivial vs a JVP.
               lg_k = jfnk_lgmres_n
               g_aug = 0d0
               if(left_pc) then
                  do i = 1, j
                     g_aug(i,i) = 1d0
                  end do
               else
                  g_aug(1:j,1:j) = jfnk_gmres_g(1:j,1:j)
               end if
               do lg_i = 1, lg_k
                  do lg_l = 1, j
                     if(left_pc) then
                        s_tmp = dotprd(n, v(1,lg_l), &
                                       jfnk_lgmres_d(1,lg_i))
                     else
                        s_tmp = dotprd(n, jfnk_gmres_z(1,lg_l), &
                                       jfnk_lgmres_d(1,lg_i))
                     end if
                     g_aug(lg_l, j+lg_i) = s_tmp
                     g_aug(j+lg_i, lg_l) = s_tmp
                  end do
                  do lg_l = 1, lg_i
                     s_tmp = dotprd(n, jfnk_lgmres_d(1,lg_i), &
                                    jfnk_lgmres_d(1,lg_l))
                     g_aug(j+lg_i, j+lg_l) = s_tmp
                     g_aug(j+lg_l, j+lg_i) = s_tmp
                  end do
               end do
               call hook_step_g_v1(j+lg_k, h_aug, m+lg_kmax, beta, &
                    del_cycle, g_aug, m+lg_kmax, y_aug, pred_res)
               y(1:j) = y_aug(1:j)
            else if(left_pc) then
               call hook_step_v1(j, h, m, beta, del_cycle, y, pred_res)
            else
               call hook_step_g_v1(j, h, m, beta, del_cycle, &
                    jfnk_gmres_g, size(jfnk_gmres_g,1), y, pred_res)
            end if
         else
            ! No trust region (del==0): use the unconstrained
            ! least-squares residual as the predicted residual.
            pred_res = res
         end if
         ! Restart updates may temporarily exceed del. Project only the
         ! returned total step so a single-cycle solve keeps its original
         ! hookstep while restart accumulation remains bounded on return.
         if(lg_active .and. jfnk_lgmres_n > 0) then
            ! Step = basis*y_aug(1:j) + D*y_aug(j+1:j+k): the basis is
            ! V under left PC (Arnoldi already runs in x-space) and
            ! the stored Z = M^{-1}V under right PC.
            if(left_pc) then
               z = matmul(v(:,1:j), y_aug(1:j)) + &
                   matmul(jfnk_lgmres_d(:,1:jfnk_lgmres_n), &
                          y_aug(j+1:j+jfnk_lgmres_n))
            else
               z = matmul(jfnk_gmres_z(:,1:j), y_aug(1:j)) + &
                   matmul(jfnk_lgmres_d(:,1:jfnk_lgmres_n), &
                          y_aug(j+1:j+jfnk_lgmres_n))
            end if
         else if(left_pc) then
            z = matmul(v(:,1:j), y(1:j))
         else
            z = matmul(jfnk_gmres_z(:,1:j), y(1:j))
         end if
         if(lg_active) then
            ! Retain this cycle's step for the next restart's pair:
            ! its image is recovered from the residual difference at
            ! the next cycle start, so no storage of A z is needed.
            jfnk_lgmres_dlast = z
            lg_have_last = .true.
         end if
         x = x + z
         if(del > 0d0 .and. done) then
            x_cycle_norm = dsqrt(max(0d0, dotprd(n, x, x)))
            if(x_cycle_norm > del) x = x * (del / x_cycle_norm)
         end if
         jfnk_last_ynrm = dsqrt(max(0d0, dotprd(n, x, x)))
         if(its == imx) info = 2
         if(n_stag >= 3) info = 1   ! 2026-09-13 TR-units fix: stagnation, not any rise
         if(res <= tol) info = 0
         jj = j
         ! Checkpoint the completed cycle for the info==2 hookstep
         ! recompute: caller's h/v hold this basis and beta is this
         ! cycle's starting residual norm.  (2026-10 audit fix)
         jfnk_gmres_jj   = j
         jfnk_gmres_beta = beta
         if(done .or. jvp_soft) then
            !
            ! JVP soft-fail exit (gdefl work): a killed column closed
            ! the cycle above and x already holds the step built from
            ! the completed basis.  Measure the true residual when the
            ! operator still evaluates -- an eval that fails leaves the
            ! direction unverifiable and the caller's stagnation gate
            ! falls back to the same abort the old hard exit produced.
            !
            if(jvp_soft) then
               ! Steps accumulate across cycles, so x can exceed the
               ! trust radius even though each cycle's step was
               ! hookstep-bounded.  Rescale like the restart-cap exit
               ! before the caller sees the direction.
               if(del > 0d0) then
                  x_cycle_norm = dsqrt(max(0d0, dotprd(n, x, x)))
                  if(x_cycle_norm > del) x = x * (del / x_cycle_norm)
                  jfnk_last_ynrm = dsqrt(max(0d0, dotprd(n, x, x)))
               end if
               if(.not. ltres) then
                  ! The flag still marks the killed column; this eval
                  ! is a fresh direction, so the stale flag must not
                  ! veto a measurement that succeeds.
                  jfnk_jvp_failed = .false.
                  call matvec(n, x, z)
                  jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
                  if(.not. jfnk_jvp_failed) then
                     w = b - z
                     if(left_pc) call psolve(n, w)
                     if(.not. jfnk_pc_clip_event) then
                        true_res = dsqrt(dotprd(n, w, w))
                        ltres = .true.
                        jfnk_gmres_true_res = true_res
                        jfnk_gmres_true_res_ok = .true.
                     end if
                  end if
               end if
               info = 1
               jfnk_exit_reason = 'gmres_jvp_soft_exit'
               return
            end if
            ! [JFNK_TRES_FAC] (default 0 = off): measure the true
            ! weighted linear residual ||b - A x||_w at every done-exit
            ! with one extra matvec.  The Hessenberg LS residual is only
            ! a model of that quantity; with a finite-difference matvec
            ! the Arnoldi vectors need not span a fixed linear operator,
            ! and hess_res has been observed 150-450x below the true
            ! residual at restart (RESTART DIAG, npq_dbg run).  A
            ! claimed convergence (res <= tol) is accepted only when
            ! true_res <= JFNK_TRES_FAC*tol; a rejected claim restarts
            ! the cycle on the measured residual (bounded by the
            ! restart cap and imx), or exits 'gmres_tres_reject' when no
            ! restart headroom remains.  A non-finite true_res at ANY
            ! done-exit exits 'gmres_tres_nonfinite' immediately --
            ! non-converged exits otherwise only record true_res.
            tres_rej = .false.
            if(jfnk_gmres_tres_fac > 0d0) then
               jfnk_jvp_failed = .false.
               call matvec(n, x, z)
               jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
               ! Soft fail: the step is formed; an unverifiable
               ! verification eval hands info=1 to the caller gate
               ! (which aborts identically) instead of discarding x.
               if(jfnk_jvp_failed) then
                  if(luer > 0) write(luer,*) '  GMRES: JVP failure'// &
                     ' in TRES check -- soft exit'
                  info = 1
                  jfnk_exit_reason = 'gmres_jvp_soft_exit'
                  return
               end if
               w = b - z
               if(left_pc) then
                  call psolve(n, w)
                  if(jfnk_pc_clip_event) then
                     if(luer > 0) write(luer,*) '  GMRES: precond'// &
                        ' clip in TRES check -- aborting'
                     info = 3
                     jfnk_exit_reason = 'gmres_pc_clip_tres'
                     return
                  end if
               end if
               true_res = dsqrt(dotprd(n, w, w))
               ltres = .true.
               jfnk_gmres_true_res = true_res
               jfnk_gmres_true_res_ok = .true.
               ! info already holds the exit code here (0 on a claimed
               ! convergence), so it cannot double as the verbose gate:
               ! this line exists to audit exactly those claims.
               if(luer > 0) &
                  write(luer,'(A,I3,3(A,ES12.4))') &
                  '  GMRES TRES: it=', its, ' hess_res=', real(res), &
                  ' true_res=', real(true_res), ' acc_tol=', &
                  real(jfnk_gmres_tres_fac * tol)
               ! A non-finite true_res cannot verify ANY done-exit, and
               ! the measured residual vector is itself non-finite: the
               ! ordinary paths would restart the cycle and re-run
               ! matvecs on an operator that just produced NaN/Inf (the
               ! lb_hit stagnation bypass reaches the restart at label 2
               ! even without a convergence claim).  Exit non-converged
               ! with a dedicated reason instead of restarting
               ! (07/10-Nov-2026 audits).  ieee_is_finite covers NaN and
               ! Inf -- Inf already fails the > test below but belongs
               ! to this unverifiable class, not the "too big" one.
               if(.not. ieee_is_finite(true_res)) then
                  if(luer > 0) write(luer,'(A)') &
                     '  GMRES: non-finite true residual in TRES '// &
                     'check -- exit unverifiable, not restarting'
                  info = 1
                  jfnk_exit_reason = 'gmres_tres_nonfinite'
                  return
               end if
               if(res <= tol .and. &
                  true_res > jfnk_gmres_tres_fac * tol) tres_rej = .true.
            end if
            if(tres_rej) then
               if(its == imx .or. &
                  n_restart >= jfnk_gmres_max_restart) then
                  if(luer > 0) write(luer,'(A)') &
                     '  GMRES: convergence claim failed the true-'// &
                     'residual check; no restart headroom left'
                  info = 1
                  jfnk_exit_reason = 'gmres_tres_reject'
                  return
               end if
               if(luer > 0) write(luer,'(A)') &
                  '  GMRES: convergence claim failed the true-'// &
                  'residual check -- restarting on the measured residual'
               goto 2
            end if
            ! Convergence and the total-iteration cap always exit.
            ! Stagnation exit (n_stag>=3) is bypassed when lb_hit ended
            ! the cycle: the stall was caused by the invariant Krylov
            ! space, not slow progress, so a fresh basis is worth it.
            ! Restarts are bounded by imx (its keeps accumulating).
            if(res <= tol .or. its == imx .or. .not. lb_hit) then
               ! editted by Devin, 12-Nov-2026: an info=1/2 exit hands
               ! the caller only the Hessenberg model residual, which a
               ! finite-difference matvec can leave orders below the
               ! true ||b - A x||_w.  Measure once here so the Newton
               ! stagnation gate tests the real residual.  Converged
               ! exits (info=0) and cycles already measured by the TRES
               ! check (ltres) add no extra matvec.
               if(info /= 0 .and. .not. ltres) then
                  jfnk_jvp_failed = .false.
                  call matvec(n, x, z)
                  jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
                  if(jfnk_jvp_failed) then
                     if(luer > 0) write(luer,*) '  GMRES: JVP failure'// &
                        ' in exit true-res measure -- soft exit'
                     info = 1
                     jfnk_exit_reason = 'gmres_jvp_soft_exit'
                     return
                  end if
                  w = b - z
                  if(left_pc) then
                     call psolve(n, w)
                     if(jfnk_pc_clip_event) then
                        if(luer > 0) write(luer,*) '  GMRES: precond'// &
                           ' clip in exit true-res measure -- aborting'
                        info = 3
                        jfnk_exit_reason = 'gmres_pc_clip_tres_exit'
                        return
                     end if
                  end if
                  true_res = dsqrt(dotprd(n, w, w))
                  ltres = .true.
                  jfnk_gmres_true_res = true_res
                  jfnk_gmres_true_res_ok = .true.
               end if
               if(info == 0) then
                  jfnk_exit_reason = 'gmres_converged'
               else if(its == imx) then
                  jfnk_exit_reason = 'gmres_max_iterations'
               else
                  jfnk_exit_reason = 'gmres_stagnation'
               end if
               if(luer > 0) then
                  if(info == 0) then
                     write(luer,'(A,I4,A,ES12.4)') &
                        '  GMRES: CONVERGED  its=', its, ' res=', real(res)
                  else if(its == imx) then
                     write(luer,'(A,I4,A,ES12.4)') &
                        '  GMRES: MAXIMUM ITERATIONS REACHED  its=', its, &
                        ' res=', real(res)
                  end if
               end if
               return
            end if
         end if
2        continue
         !
         ! V4 audit (2026-10): True residual at restart.
         ! The Hessenberg least-squares residual `res` is computed with
         ! dot_product (Euclidean) on the coefficient vector, while the
         ! Arnoldi basis V is orthonormal in the weighted (dotprd) inner
         ! product.  If V loses orthogonality, res underestimates the
         ! true weighted physical residual ||F - J*s||_w.
         !
         ! Compute the true residual explicitly and compare.
         !
         if(info == 1) then
            ! When the TRES check above already measured ||b - A x||_w
            ! this cycle, reuse it -- a second eval on the same vector
            ! buys nothing and doubles the diagnostic cost.
            if(.not. ltres) then
               w = b
               jfnk_jvp_failed = .false.
               call matvec(n, x, z)
               jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
               if(jfnk_jvp_failed) then
                  write(luer,*) '  GMRES: JVP failure in restart diag -- soft exit'
                  info = 1
                  jfnk_exit_reason = 'gmres_jvp_soft_exit'
                  return
               end if
               w = b - z
               if(left_pc) then
                  call psolve(n, w)
                  if(jfnk_pc_clip_event) then
                     write(luer,*) '  GMRES: precond clip in restart diag -- aborting'
                     info = 3
                     jfnk_exit_reason = 'gmres_pc_clip_restart_diag'
                     return
                  end if
               end if
               true_res = dsqrt(dotprd(n, w, w))
               ltres = .true.
               jfnk_gmres_true_res = true_res
               jfnk_gmres_true_res_ok = .true.
            end if
            write(luer,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)') &
               '  GMRES RESTART DIAG: it=', its, &
               ' hess_res=', real(res), &
               ' true_res=', real(true_res), &
               ' ratio=', real(true_res / max(res, 1d-30))
         end if
         !
         ! [JFNK_CYCSTALL] (default F): per-cycle contraction guard.
         ! A cycle that trims the LS residual by less than
         ! jfnk_cycstall_frac is re-solving a nearly invariant space --
         ! the JVP budget is better spent on the inexact direction
         ! than on jfnk_cycstall_nmax more identical cycles.  Exit via
         ! info=1 like the other non-converged exits; the caller's
         ! true-residual gate decides whether the partial direction
         ! is usable.
         if(jfnk_cycstall_on) then
            if(res_prev_end > 0d0) then
               if(res > res_prev_end * (1d0 - jfnk_cycstall_frac)) then
                  n_cycstall = n_cycstall + 1
               else
                  n_cycstall = 0
               end if
               if(n_cycstall >= jfnk_cycstall_nmax) then
                  if(luer > 0) write(luer,'(A,I3,A,ES12.4)') &
                     '  GMRES: cycle stall (', n_cycstall, &
                     ' flat cycles) -- exiting, res=', real(res)
                  info = 1
                  jfnk_exit_reason = 'gmres_cycle_stall'
                  return
               end if
            end if
            res_prev_end = res
         end if
         !
         ! Restart cap: a Krylov space that goes invariant because the
         ! preconditioned operator cannot reach the residual tail stays
         ! invariant on restart.  Cycles then die in ~2 its and restart
         ! storm accumulates bad LSQ steps until imx.  The cap is a
         ! storm guard only -- set [JFNK_GMRES_RST] above 1 when the
         ! goal is to exercise the full JFNK_MAXGMRES budget; the
         ! true-residual diagnostic above prints at every restart so
         ! a persistently flat cycle is recorded, not hidden.  On
         ! exceeding the cap, fall back to the caller's
         ! inexact-direction path via info=1.
         n_restart = n_restart + 1
         if(n_restart > jfnk_gmres_max_restart) then
            if(del > 0d0) then
               x_cycle_norm = dsqrt(max(0d0, dotprd(n, x, x)))
               if(x_cycle_norm > del) then
                  x = x * (del / x_cycle_norm)
                  ltres = .false.
                  jfnk_gmres_true_res_ok = .false.
               end if
               jfnk_last_ynrm = dsqrt(max(0d0, dotprd(n, x, x)))
            end if
            ! Same unverifiable-direction hazard as the stagnation exit
            ! above: x is re-scaled to the trust radius here, so a TRES
            ! measurement taken on the unscaled x is stale -- measure
            ! after the rescale when no current-cycle value exists.
            if(.not. ltres) then
               jfnk_jvp_failed = .false.
               call matvec(n, x, z)
               jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
               if(jfnk_jvp_failed) then
                  if(luer > 0) write(luer,*) '  GMRES: JVP failure'// &
                     ' in exit true-res measure -- soft exit'
                  info = 1
                  jfnk_exit_reason = 'gmres_jvp_soft_exit'
                  return
               end if
               w = b - z
               if(left_pc) then
                  call psolve(n, w)
                  if(jfnk_pc_clip_event) then
                     if(luer > 0) write(luer,*) '  GMRES: precond'// &
                        ' clip in exit true-res measure -- aborting'
                     info = 3
                     jfnk_exit_reason = 'gmres_pc_clip_tres_exit'
                     return
                  end if
               end if
               true_res = dsqrt(dotprd(n, w, w))
               jfnk_gmres_true_res = true_res
               jfnk_gmres_true_res_ok = .true.
            end if
            if(luer > 0) write(luer,'(A,I3,A)') &
               '  GMRES: restart cap (', n_restart, &
               ') -- treating as stagnation'
            info = 1
            jfnk_exit_reason = 'gmres_restart_cap_stagnation'
            return
         end if
         ! edited by Devin, 05-Oct-2026: luer>0 guard, same convention
         ! as every other write in this routine.
         if(del > 0d0 .and. luer > 0)  write(luer,*) '  GMRES: warning! &
            &restart affects hookstep'
         goto 1       ! (j==m) restart
      end if
      !
   end do
   !
end subroutine gmres_m_v1


!----------------------------------------------------------------------
! JFNK_GDEFL_EXTRACT_V1: harmonic-Ritz deflation for GMRES-DR.
!
! At a GMRES restart, the k smallest-|theta| harmonic Ritz pairs of
! the completed cycle's Hessenberg block approximate the eigenmodes
! of the preconditioned operator B that restarted GMRES discards and
! re-discovers each cycle (Morgan 2002, SISC 24(1) 20-37; the
! Parks et al. 2006 GCRO-DR uses the same machinery for sequences of
! systems).  The retained directions fill the augmentation slots the
! LGMRES machinery already consumes, so the next cycle's dense
! least-squares sees span[V_jj g_i] jointly with the new Krylov
! space.
!
! Harmonic Ritz problem (Morgan 2002 eq. 4.1): the pairs (theta, g)
! are eigenpairs of Hhat = H + h_{jj+1,jj}^2 f e_jj^T with
! f = H^{-T} e_jj.  The retained step-space direction is V_jj g under
! left PC and Z g = M^{-1} V_jj g under right PC -- matching the
! [V D]/[Z D] step the augmented LS builds.  Its operator image is
! V_{jj+1} H_bar g by the Arnoldi relation B V_jj = V_{jj+1} H_bar,
! so extraction needs no extra matvecs.  A complex-conjugate pair
! contributes its real 2-d invariant span (two vectors).
!
! [JFNK_GDEFL] (default F) gates the caller site; [JFNK_GDEFL_K]
! bounds the mode count.  All failures (singular H, dgeev error,
! non-finite vectors) degrade to "keep the existing set" -- the
! augmented LS only treats the set as constraint directions.
!----------------------------------------------------------------------
subroutine jfnk_gdefl_extract_v1(n, m, jj, h, ldh, v, ldv, dotprd, &
                                 left_pc, keep_lg, luer)
   use, intrinsic :: ieee_arithmetic
   use jfnk_solver_v1_mod, only: jfnk_gmres_z, jfnk_gdefl_kmax, &
                                jfnk_lgmres_d, jfnk_lgmres_ad, jfnk_lgmres_n
   implicit none
   integer,          intent(in) :: n, m, jj, ldh, ldv, luer
   double precision, intent(in) :: h(ldh,m), v(ldv,m+1)
   double precision, external   :: dotprd
   logical,          intent(in) :: left_pc, keep_lg
   !
   integer :: i, l, k_g, k_cap, info_e, lwk, is, ns
   double precision :: th, dnorm
   double precision, allocatable :: hh(:,:), ht(:,:), f(:), wr(:), &
                                    wi(:), vr(:,:), wk(:), gs(:,:), &
                                    hg(:), dv(:), av(:)
   integer, allocatable :: ipiv(:), ord(:)
   logical, allocatable :: used(:)
   !
   if(jj < 2) return
   if(.not. left_pc .and. .not. allocated(jfnk_gmres_z)) return
   if(.not. allocated(jfnk_lgmres_d)) return
   !
   ! Slot budget: reserve one for the LGMRES last-step pair when that
   ! retention is also armed.
   k_cap = size(jfnk_lgmres_d,2)
   if(keep_lg) k_cap = k_cap - 1
   k_cap = min(k_cap, jfnk_gdefl_kmax)
   if(k_cap < 1) return
   !
   lwk = 8*jj
   allocate(hh(jj,jj), ht(jj,jj), f(jj), wr(jj), wi(jj), vr(jj,jj), &
            wk(lwk), gs(jj,k_cap), hg(jj+1), dv(n), av(n), ipiv(jj), &
            ord(jj), used(jj), stat=info_e)
   if(info_e /= 0) return
   !
   ! f = H^{-T} e_jj, then the rank-1 update on the last column.
   !
   ht = transpose(h(1:jj,1:jj))
   f = 0d0
   f(jj) = 1d0
   call dgesv(jj, 1, ht, jj, ipiv, f, jj, info_e)
   if(info_e /= 0) then
      if(luer > 0) write(luer,'(A,I3)') &
         '  GDEFL: H^T solve failed, info=', info_e
      goto 90
   end if
   hh = h(1:jj,1:jj)
   hh(:,jj) = hh(:,jj) + h(jj+1,jj)**2 * f
   !
   call dgeev('N', 'V', jj, hh, jj, wr, wi, ht, 1, vr, jj, wk, lwk, &
              info_e)
   if(info_e /= 0) then
      if(luer > 0) write(luer,'(A,I3)') &
         '  GDEFL: dgeev failed, info=', info_e
      goto 90
   end if
   !
   ! Ascending-|theta| order (insertion sort; jj <= m).
   !
   do i = 1, jj
      ord(i) = i
   end do
   do i = 2, jj
      is = ord(i)
      th = wr(is)**2 + wi(is)**2
      do l = i-1, 1, -1
         if(wr(ord(l))**2 + wi(ord(l))**2 <= th) exit
         ord(l+1) = ord(l)
         ord(l) = is
      end do
   end do
   !
   ! Collect real coefficient vectors: one per real eigenvalue, two
   ! per complex-conjugate pair (the real invariant span).  A partner
   ! is dropped only when the slot budget is spent.
   !
   used = .false.
   ns = 0
   do i = 1, jj
      if(ns >= k_cap) exit
      is = ord(i)
      if(used(is)) cycle
      if(wi(is) == 0d0) then
         ns = ns + 1
         gs(:,ns) = vr(:,is)
         used(is) = .true.
      else if(wi(is) > 0d0) then
         ns = ns + 1
         gs(:,ns) = vr(:,is)
         used(is) = .true.
         if(is+1 <= jj) then
            used(is+1) = .true.
            if(ns < k_cap) then
               ns = ns + 1
               gs(:,ns) = vr(:,is+1)
            end if
         end if
      else
         ! wi < 0: conjugate already taken with its wi > 0 partner.
         used(is) = .true.
      end if
   end do
   !
   ! Map coefficients to (step-space direction, operator image) pairs
   ! and unit-normalize the direction -- the image is scaled by the
   ! same factor so the pair stays consistent.
   !
   k_g = 0
   do i = 1, ns
      hg = matmul(h(1:jj+1,1:jj), gs(:,i))
      av = matmul(v(:,1:jj+1), hg)
      if(left_pc) then
         dv = matmul(v(:,1:jj), gs(:,i))
      else
         dv = matmul(jfnk_gmres_z(:,1:jj), gs(:,i))
      end if
      if(.not. all(ieee_is_finite(dv)) .or. &
         .not. all(ieee_is_finite(av))) cycle
      dnorm = dsqrt(max(0d0, dotprd(n, dv, dv)))
      if(dnorm < 1d-30) cycle
      k_g = k_g + 1
      jfnk_lgmres_d(:,k_g)  = dv / dnorm
      jfnk_lgmres_ad(:,k_g) = av / dnorm
   end do
   jfnk_lgmres_n = k_g
   !
   if(luer > 0) write(luer,'(A,I2,A,ES12.4)') &
      '  GDEFL: retained ', k_g, ' Ritz dirs, |theta|_min=', &
      real(sqrt(wr(ord(1))**2 + wi(ord(1))**2))
   !
90 continue
   if(allocated(hh)) deallocate(hh, ht, f, wr, wi, vr, wk, gs, hg, &
                                dv, av, ipiv, ord, used)
   return
   !
end subroutine jfnk_gdefl_extract_v1


!----------------------------------------------------------------------
! HOOK_STEP_V1: Compute the hookstep (trust-region constrained step)
! using SVD of the upper Hessenberg matrix.
!
! Based on Viswanath (2008) arXiv:0809.1498 and the JFNK-Hookstep
! reference implementation.
!
! Bug fix (09-Sep-2026): del (trust radius) is now INTENT(IN) and
! the predicted residual norm is returned separately as pred_res.
! Previously del was INTENT(INOUT) and was overwritten with the
! predicted residual, which corrupted the trust radius on GMRES
! restart and gave wrong pred when del==0 (no trust region).
!
! editted by Sean Brennan (V4, 2026-09)
!----------------------------------------------------------------------
subroutine hook_step_v1(j, h, m, beta, del, y, pred_res)
   implicit none
   integer,          intent(in)    :: j, m
   double precision, intent(in)    :: h(m+1, j)
   double precision, intent(in)    :: beta
   double precision, intent(in)    :: del
   double precision, intent(out)   :: y(j)
   double precision, intent(out)   :: pred_res
   !
   double precision :: a(j+1, j), s(j), u(j+1, j+1), vt(j, j)
   double precision :: work(8*(j+1)+16)
   double precision :: p(j+1), q(j), mu, qn, s_tol
   double precision :: dq_dmu, dmu
   integer :: info, i, i_ms
   !
   ! Guard: j=0 means no GMRES iterations were done.  There is no
   ! least-squares problem to solve.  Return a zero step to let the
   ! Newton loop evaluate the residual and decide what to do next.
   ! Without this guard, DGESVD is called with LDVT=0 < max(1,N)=1,
   ! causing "Intel MKL ERROR: Parameter 11 was incorrect on entry
   ! to DGESVD."
   !
   if(j < 1) then
      y = 0d0
      pred_res = beta   ! zero step → predicted residual = beta (no improvement)
      return
   end if
   !
   a = h(1:j+1, 1:j)
   !
   call dgesvd('A', 'A', j+1, j, a, j+1, s, u, j+1, vt, j, work, &
                8*(j+1)+16, info)
   if(info /= 0) then
      ! SVD failed — fall back to zero step.
      y = 0d0
      pred_res = beta
      return
   end if
   !
   p(1:j) = beta * u(1, 1:j)
   !
   ! Singular-value threshold (ChatGPT audit fix #2, 2026-09-09):
   ! The Hessenberg least-squares system can be rank-deficient when
   ! fixed T/Ne variables and masked equations remove unknowns.  A zero
   ! (or near-zero) singular value makes q = p/s inject Inf/NaN, which
   ! then corrupts the trust-region step and the predicted residual.
   ! DGESVD returns singular values in descending order, so s(1) is the
   ! largest; threshold relative to it and explicitly remove null
   ! directions instead of dividing by them.
   !
   s_tol = max(s(1) * 1.0d-14, 1.0d-300)
   q = 0d0
   do i = 1, j
      if(s(i) > s_tol) q(i) = p(i) / s(i)
   end do
   qn = dsqrt(dot_product(q, q))
   if(qn <= del) then
      mu = 0d0
   else
      ! Moré-Sorensen (1983): Newton iteration on the secular equation
      !   phi(mu) = 1/||q(mu)|| - 1/del = 0,   q_i = s_i p_i/(mu + s_i^2)
      ! with
      !   d||q||/dmu = -(1/||q||) * sum_i q_i^2 / (mu + s_i^2)
      !   mu_new = mu + ||q||*(1 - ||q||/del) / (d||q||/dmu)
      ! Both factors are negative while ||q|| > del, so mu strictly
      ! increases and the iteration converges quadratically (the
      ! previous mu <- 1.1*mu search converged only linearly and
      ! could overshoot the required trust-region radius).
      !
      mu = max(s(j)*s(j)*1d-6, 1d-99)
      do i_ms = 1, 50
         qn = 0d0
         dq_dmu = 0d0
         do i = 1, j
            if(s(i) > s_tol) then
               q(i) = p(i) * s(i) / (mu + s(i)*s(i))
               qn = qn + q(i)*q(i)
               dq_dmu = dq_dmu + q(i)*q(i) / (mu + s(i)*s(i))
            else
               q(i) = 0d0
            end if
         end do
         qn = dsqrt(qn)
         if(qn <= del) exit
         dq_dmu = -dq_dmu / qn
         if(dq_dmu == 0d0) exit
         dmu = qn * (1.0d0 - qn/del) / dq_dmu
         mu = mu + dmu
         if(mu <= 0d0) mu = 1.0d-99
         if(abs(dmu) <= 1.0d-12 * mu) exit
      end do
      qn = dsqrt(dot_product(q, q))
      if(qn > del .and. qn > 0d0) q = q * (del / qn)
   end if
   !
   y = matmul(q, vt)
   !
   ! Predicted residual for the hookstep
   !
   p = -matmul(h(1:j+1, 1:j), y(1:j))
   p(1) = p(1) + beta
   pred_res = dsqrt(dot_product(p, p))
   !
end subroutine hook_step_v1


!----------------------------------------------------------------------
! HOOK_STEP_G_V1: hookstep under a general G-norm trust constraint
!   y' G y <= del^2.
!
! Used under right preconditioning: the physical step is s = Z y with
! Z = M^{-1} V, and G = Z' W Z makes the constraint ||s||_w rather
! than ||y|| = ||M s||_w (the "TR norm invalid" issue noted in the
! 17-Sep-2026 audit).  Cholesky G = R' R maps to the Euclidean
! problem in yhat = R y:
!   min ||beta e1 - H R^{-1} yhat||   s.t.  ||yhat|| <= del,
! solved by the standard hook_step_v1, then y = R^{-1} yhat.
!
! A numerically singular G (preconditioner collapse, e.g. a psolve
! that maps a basis vector to ~0) falls back to the Euclidean
! hookstep -- the radius is then in ||y|| units, which is wrong but
! still bounded.
!
! edited by Devin, 2026-10
!----------------------------------------------------------------------
subroutine hook_step_g_v1(j, h, m, beta, del, g, ldg, y, pred_res)
   implicit none
   integer,          intent(in)  :: j, m, ldg
   double precision, intent(in)  :: h(m+1,j), beta, del, g(ldg,m)
   double precision, intent(out) :: y(j), pred_res
   !
   double precision :: r(j,j), ht(m+1,j), gmax
   integer :: i, info
   !
   if(j < 1) then
      y = 0d0
      pred_res = beta
      return
   end if
   !
   r = g(1:j,1:j)
   ! Ridge relative to the Gram scale: psolve output can make Z
   ! columns nearly dependent without being exactly singular.
   gmax = max(maxval(abs(r)), 1d-300)
   do i = 1, j
      r(i,i) = r(i,i) + 1.0d-12 * gmax
   end do
   call dpotrf('U', j, r, j, info)
   if(info /= 0) then
      ! Fallback: Euclidean hookstep in y-units.  This is not the
      ! physical ||s||_w norm the G-factorization was built to
      ! constrain -- the step is still bounded by del, but in the
      ! wrong metric.  Rare: the ridge above already absorbs
      ! near-dependent Z columns.
      call hook_step_v1(j, h, m, beta, del, y, pred_res)
      return
   end if
   ! Ht = H R^{-1} via a right-side triangular solve Ht * R = H.
   ht(1:j+1,1:j) = h(1:j+1,1:j)
   call dtrsm('R', 'U', 'N', 'N', j+1, j, 1.0d0, r, j, ht, m+1)
   call hook_step_v1(j, ht, m, beta, del, y, pred_res)
   ! y currently holds yhat; overwrite with R^{-1} yhat.
   call dtrsv('U', 'N', 'N', j, r, j, y, 1)
end subroutine hook_step_g_v1


!----------------------------------------------------------------------
! JFNK_NEWTON_DIAGNOSTIC_V1: Comprehensive Newton direction diagnosis.
!
! Printed once on the first Newton step to understand why the computed
! direction does not provide a useful physical reduction.
!
! Diagnostics:
!   1. Sign verification: ||F - J*s|| vs ||F + J*s||
!   2. Directional derivative: <F, J*s>_w (descent direction test)
!   3. Component-wise top-20 analysis of |s(i)|
!   4. Alpha-sweep of linearized model ||F - alpha*J*s||_w
!   5. RSCALE audit (min/max/ratio)
!   6. GMRES convergence in both norms
!   7. JVP-vs-BA comparison for the actual Newton direction
!
! editted by Sean Brennan (V4, 2026-09)
!----------------------------------------------------------------------
subroutine jfnk_newton_diagnostic_v1(n, s, fx, J_s, F_pred, &
                                      dp, df, luer)
   use jfnk_solver_v1_mod, only: jfnk_total_jvp_evals
   use jfnk_cmfgen_int_mod, only: JFNK_NT, JFNK_ND, JFNK_RSCALE, &
                                  JFNK_RSCALE_READY, JFNK_BA_DIAG, &
                                  JFNK_BA_DIAG_READY, JFNK_BA_A, &
                                  JFNK_LU_READY, JFNK_BA_RIDGE_ADD
   implicit none
   integer,          intent(in) :: n, luer
   double precision, intent(in) :: s(n), fx(n), J_s(n), F_pred(n)
   double precision, external   :: dp
   external                      :: df
   !
   double precision :: Fp_minus(n), Fp_plus(n)
   double precision :: norm_F, norm_Js, norm_Fp_minus, norm_Fp_plus
   double precision :: directional_deriv
   double precision :: alpha_list(6), alpha, norm_model
   double precision :: rs_min, rs_max, rs_ratio
   double precision :: BA_s(n), diff(n), diff_norm, ba_norm
   double precision :: s_arr(n), s_temp(n)
   integer :: i, k, kk, ii, jj, irow, icol, ntop, ntop_actual
   ! W25: SAVE'd sort workspace -- one allocation per n.
   integer, allocatable, save :: idx_sorted(:)
   double precision, allocatable, save :: s_abs(:)
   !
   ! Index mapping: i = (k-1)*JFNK_NT + j, where k=depth, j=level
   ! j=1..NT-2: populations, j=NT-1: electron density, j=NT: temperature
   !
   integer :: depth, level
   character(len=20) :: var_name
   !
   if(luer <= 0) return
   !
   write(luer,'(A)') ' '
   write(luer,'(A)') '=== JFNK Newton Direction Diagnostic ==='
   write(luer,'(A)') ' '
   !
   ! --- 1. Sign verification ---
   !
   ! The code uses: jfnk_x = x_ - s, F_pred = fx - J_s
   ! If GMRES solves J*s = F, then F_pred = F - J*s should be small.
   ! If GMRES solves J*s = -F, then F + J*s should be small.
   !
   norm_F = dsqrt(dp(n, fx, fx))
   norm_Js = dsqrt(dp(n, J_s, J_s))
   !
   Fp_minus = fx - J_s   ! F - J*s (current convention)
   Fp_plus = fx + J_s    ! F + J*s (opposite sign)
   norm_Fp_minus = dsqrt(dp(n, Fp_minus, Fp_minus))
   norm_Fp_plus = dsqrt(dp(n, Fp_plus, Fp_plus))
   !
   write(luer,'(A)') '--- 1. Sign Verification ---'
   write(luer,'(A,ES14.4)') '  ||F||_w           = ', norm_F
   write(luer,'(A,ES14.4)') '  ||J*s||_w         = ', norm_Js
   write(luer,'(A,ES14.4)') '  ||F - J*s||_w     = ', norm_Fp_minus
   write(luer,'(A,ES14.4)') '  ||F + J*s||_w     = ', norm_Fp_plus
   write(luer,'(A,ES14.4)') '  ||F-J*s||/||F||   = ', norm_Fp_minus / max(norm_F, 1.0d-30)
   write(luer,'(A,ES14.4)') '  ||F+J*s||/||F||   = ', norm_Fp_plus / max(norm_F, 1.0d-30)
   if(norm_Fp_minus < norm_Fp_plus) then
      write(luer,'(A)') '  => F - J*s is smaller: GMRES solves J*s = F (sign correct)'
   else
      write(luer,'(A)') '  => F + J*s is smaller: GMRES may solve J*s = -F (sign WRONG?)'
   end if
   write(luer,'(A)') ' '
   !
   ! --- 2. Directional derivative <F, J*s>_w ---
   !
   ! For phi = 0.5*||F||_w^2, the derivative along x(alpha) = x - alpha*s is:
   !   dphi/dalpha = -<F, J*s>_w
   ! A descent direction requires <F, J*s>_w > 0.
   !
   directional_deriv = dp(n, fx, J_s)
   write(luer,'(A)') '--- 2. Directional Derivative ---'
   write(luer,'(A,ES14.4)') '  <F, J*s>_w        = ', directional_deriv
   write(luer,'(A,ES14.4)') '  ||F||_w           = ', norm_F
   write(luer,'(A,ES14.4)') '  ||J*s||_w         = ', norm_Js
   write(luer,'(A,ES14.4)') '  cos(F, J*s)_w     = ', &
      directional_deriv / max(norm_F * norm_Js, 1.0d-30)
   if(directional_deriv > 0.0d0) then
      write(luer,'(A)') '  => <F,J*s>_w > 0: descent direction (good)'
   else if(directional_deriv < 0.0d0) then
      write(luer,'(A)') '  => <F,J*s>_w < 0: NOT a descent direction (bad)'
   else
      write(luer,'(A)') '  => <F,J*s>_w ~ 0: direction is orthogonal to F'
   end if
   write(luer,'(A)') ' '
   !
   ! --- 3. Component-wise top-20 analysis ---
   !
   write(luer,'(A)') '--- 3. Top-20 Components of |s(i)| ---'
   if(allocated(s_abs)) then
      if(size(s_abs) .ne. n) deallocate(s_abs, idx_sorted)
   end if
   if(.not. allocated(s_abs)) allocate(s_abs(n), idx_sorted(n))
   do i = 1, n
      s_abs(i) = abs(s(i))
      idx_sorted(i) = i
   end do
   ! Simple selection sort for top 20 (n can be large, but 20 is small)
   ntop = min(20, n)
   do i = 1, ntop
      do k = i+1, n
         if(s_abs(k) > s_abs(i)) then
            ! Swap
            alpha = s_abs(i); s_abs(i) = s_abs(k); s_abs(k) = alpha
            jj = idx_sorted(i); idx_sorted(i) = idx_sorted(k); idx_sorted(k) = jj
         end if
      end do
   end do
   write(luer,'(A)') '  rank  idx    depth  level  var          s(i)         F(i)' // &
      '        J*s(i)      F-J*s(i)    RSCALE(i)   BA_diag(i)  |Js|/|F|'
   ntop_actual = min(ntop, n)
   do i = 1, ntop_actual
      ii = idx_sorted(i)
      if(JFNK_NT > 0 .and. JFNK_ND > 0) then
         depth = (ii - 1) / JFNK_NT + 1
         level = mod(ii - 1, JFNK_NT) + 1
         if(level <= JFNK_NT - 2) then
            write(var_name,'(A,I3)') 'pop_l', level
         else if(level == JFNK_NT - 1) then
            var_name = 'Ne'
         else
            var_name = 'T'
         end if
      else
         depth = 0
         level = 0
         write(var_name,'(A)') 'unknown'
      end if
      if(JFNK_RSCALE_READY .and. allocated(JFNK_RSCALE)) then
         write(luer,'(I4,I8,I6,I5,2X,A8,5ES12.3,2ES10.2,ES10.2)') &
            i, ii, depth, level, var_name, &
            s(ii), fx(ii), J_s(ii), F_pred(ii), &
            DBLE(JFNK_RSCALE(ii)), &
            merge(DBLE(JFNK_BA_DIAG(level, depth)), 0.0d0, &
                  JFNK_BA_DIAG_READY .and. allocated(JFNK_BA_DIAG)), &
            abs(J_s(ii)) / max(abs(fx(ii)), 1.0d-30)
      else
         write(luer,'(I4,I8,I6,I5,2X,A8,4ES12.3,ES10.2)') &
            i, ii, depth, level, var_name, &
            s(ii), fx(ii), J_s(ii), F_pred(ii), &
            abs(J_s(ii)) / max(abs(fx(ii)), 1.0d-30)
      end if
   end do
   continue  ! W25: workspace now SAVE'd, freed with the process
   write(luer,'(A)') ' '
   !
   ! --- 4. Alpha-sweep of linearized model ---
   !
   write(luer,'(A)') '--- 4. Linearized Model Alpha-Sweep ---'
   write(luer,'(A)') '  alpha      ||F - alpha*J*s||_w   ||F-alpha*J*s||/||F||'
   alpha_list = (/1.0d0, 0.5d0, 0.25d0, 0.1d0, 0.01d0, 0.001d0/)
   do i = 1, 6
      alpha = alpha_list(i)
      do k = 1, n
         Fp_minus(k) = fx(k) - alpha * J_s(k)
      end do
      norm_model = dsqrt(dp(n, Fp_minus, Fp_minus))
      write(luer,'(ES10.3,2ES20.4)') alpha, norm_model, &
         norm_model / max(norm_F, 1.0d-30)
   end do
   write(luer,'(A)') ' '
   !
   ! --- 5. RSCALE audit ---
   !
   if(JFNK_RSCALE_READY .and. allocated(JFNK_RSCALE)) then
      rs_min = JFNK_RSCALE(1)
      rs_max = JFNK_RSCALE(1)
      do i = 2, n
         if(JFNK_RSCALE(i) < rs_min) rs_min = JFNK_RSCALE(i)
         if(JFNK_RSCALE(i) > rs_max) rs_max = JFNK_RSCALE(i)
      end do
      rs_ratio = rs_max / max(rs_min, 1.0d-30)
      write(luer,'(A)') '--- 5. RSCALE Audit ---'
      write(luer,'(A,ES14.4)') '  min RSCALE  = ', rs_min
      write(luer,'(A,ES14.4)') '  max RSCALE  = ', rs_max
      write(luer,'(A,ES14.4)') '  max/min     = ', rs_ratio
      write(luer,'(A)') ' '
   end if
   !
   ! --- 6. JVP-vs-BA comparison for actual Newton direction ---
   !
   if(JFNK_LU_READY .and. allocated(JFNK_BA_A)) then
      write(luer,'(A)') '--- 6. JVP vs BA_full*s (actual Newton direction) ---'
      BA_s = 0.0d0
      do kk = 1, JFNK_ND
         do ii = 1, JFNK_NT
            do jj = 1, JFNK_NT
               BA_s((kk-1)*JFNK_NT + ii) = &
                  BA_s((kk-1)*JFNK_NT + ii) + &
                  DBLE(JFNK_BA_A(ii, jj, kk)) * s((kk-1)*JFNK_NT + jj)
            end do
         end do
      end do
      ! Subtract ridge (BA_A includes ridge on diagonal)
      do kk = 1, JFNK_ND
         do ii = 1, JFNK_NT
            BA_s((kk-1)*JFNK_NT + ii) = BA_s((kk-1)*JFNK_NT + ii) - &
               DBLE(JFNK_BA_RIDGE_ADD(ii, kk)) * s((kk-1)*JFNK_NT + ii)
         end do
      end do
      diff = J_s - BA_s
      diff_norm = 0.0d0
      ba_norm = 0.0d0
      do i = 1, n
         diff_norm = max(diff_norm, abs(diff(i)))
         ba_norm = max(ba_norm, abs(BA_s(i)))
      end do
      write(luer,'(A,ES14.4)') '  ||JFD*s||_inf     = ', maxval(abs(J_s))
      write(luer,'(A,ES14.4)') '  ||BA*s||_inf      = ', ba_norm
      write(luer,'(A,ES14.4)') '  ||JFD*s - BA*s||  = ', diff_norm
      write(luer,'(A,ES14.4)') '  diff/||BA*s||     = ', &
         diff_norm / max(ba_norm, 1.0d-30)
      write(luer,'(A)') ' '
   end if
   !
   write(luer,'(A)') '=== End Newton Direction Diagnostic ==='
   write(luer,'(A)') ' '
   flush(luer)
   !
end subroutine jfnk_newton_diagnostic_v1


!----------------------------------------------------------------------
! JFNK_GMRES_DIM_EXPERIMENT_V1: Controlled GMRES dimension experiment.
!
! Runs GMRES with m = 5, 10, 20, 40 on the SAME physical state
! (same F, JVP, preconditioner, RSCALE, trust-region radius).
! For each m, reports:
!   - GMRES iterations and restart status
!   - ||s||_w, ||s||_inf
!   - ||F-Js||_w (physical linear residual)
!   - ||M^{-1}(F-Js)|| (preconditioned residual)
!   - <F,Js>_w (directional derivative)
!   - cos(F,Js) (alignment with residual)
!   - alpha_opt (optimal step length from linearized model)
!   - Trust-region scaling info
!   - Alpha sweep for m=5 and m=40
!
! The experiment uses 40 total GMRES iterations for all m values
! (fair comparison: same total work, different subspace sizes).
!
! editted by Sean Brennan (V4, 2026-09)
!----------------------------------------------------------------------
subroutine jfnk_gmres_dim_experiment_v1(n, x_cur, fx_cur, df, dfp, dp, &
                                          gres, del, luer)
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx, jfnk_total_jvp_evals, &
                                  jfnk_max_gmres_its, jfnk_gmres_jj, &
                                  jfnk_gmres_beta, jfnk_gmres_z, &
                                  jfnk_gmres_g, jfnk_gmres_x_cycle_start
   use jfnk_cmfgen_int_mod, only: JFNK_RSCALE, JFNK_RSCALE_READY
   implicit none
   integer,          intent(in) :: n, luer
   double precision, intent(in) :: x_cur(n), fx_cur(n), gres, del
   external                      :: df, dfp
   double precision, external   :: dp
   !
   integer, parameter :: n_tests = 4
   integer :: m_tests(n_tests)
   integer :: i_test, m_test, i, j, info_test, its_test, n_restarts
   ! W25: SAVE'd test/save workspace -- one heap allocation per n
   ! instead of one per experiment call (audit 17-Oct-2026).
   double precision, allocatable, save :: v_test(:,:), h_test(:)
   double precision, allocatable, save :: s_test(:), J_s_test(:), &
     &   F_pred_test(:)
   double precision :: pred_res_test, res_test
   double precision :: norm_F, norm_Js, norm_Fp, dir_deriv, cos_angle
   double precision :: alpha_list(5), alpha, norm_model, alpha_opt
   double precision :: snrm_test, s_inf_test, snrm_scaled, scale_fac_test
   integer :: save_max_gmres, save_gmres_jj
   double precision :: save_gmres_beta
   double precision, allocatable, save :: save_jfnk_x(:), &
     &   save_jfnk_fx(:)
   ! Right-PC basis/Gram checkpoint: the test solves below clobber
   ! jfnk_gmres_z/g, so save/restore them with the jj/beta pair.
   double precision, allocatable, save :: save_gmres_z(:,:), &
     &   save_gmres_g(:,:)
   double precision, allocatable, save :: save_gmres_x_cycle_start(:)
   logical :: have_gmres_zg
   logical :: have_gmres_x_cycle_start
   double precision :: rs_min, rs_max, rs_ratio
   character(len=12) :: status_str
   !
   if(luer <= 0) return
   !
   m_tests = (/5, 10, 20, 40/)
   !
   ! Save state that GMRES might modify.  The test solves below leave
   ! jfnk_gmres_jj/jfnk_gmres_beta describing the LAST test's basis in
   ! h_test/v_test (with m_test up to 40 > production m); without the
   ! restore, a later info==2 hookstep recompute would read/write past
   ! the caller's m-sized h/v/y arrays.  (2026-10 audit fix)
   if(allocated(save_jfnk_x)) then
      if(size(save_jfnk_x) .ne. n) deallocate(save_jfnk_x, &
     &   save_jfnk_fx)
   end if
   if(.not. allocated(save_jfnk_x)) allocate(save_jfnk_x(n), &
     &   save_jfnk_fx(n))
   save_jfnk_x = jfnk_x
   save_jfnk_fx = jfnk_fx
   save_max_gmres = jfnk_max_gmres_its
   save_gmres_jj   = jfnk_gmres_jj
   save_gmres_beta = jfnk_gmres_beta
   have_gmres_x_cycle_start = allocated(jfnk_gmres_x_cycle_start)
   if(have_gmres_x_cycle_start) then
      if(allocated(save_gmres_x_cycle_start)) then
         if(size(save_gmres_x_cycle_start) /= &
            size(jfnk_gmres_x_cycle_start)) &
            deallocate(save_gmres_x_cycle_start)
      end if
      if(.not. allocated(save_gmres_x_cycle_start)) &
         allocate(save_gmres_x_cycle_start( &
            size(jfnk_gmres_x_cycle_start)))
      save_gmres_x_cycle_start = jfnk_gmres_x_cycle_start
   end if
   have_gmres_zg = allocated(jfnk_gmres_z)
   if(have_gmres_zg) then
      if(allocated(save_gmres_z)) then
         if(size(save_gmres_z,1) .ne. size(jfnk_gmres_z,1) .or. &
            size(save_gmres_z,2) .ne. size(jfnk_gmres_z,2)) &
            deallocate(save_gmres_z, save_gmres_g)
      end if
      if(.not. allocated(save_gmres_z)) &
         allocate(save_gmres_z(size(jfnk_gmres_z,1), &
            size(jfnk_gmres_z,2)), &
            save_gmres_g(size(jfnk_gmres_g,1), &
            size(jfnk_gmres_g,2)))
      save_gmres_z = jfnk_gmres_z
      save_gmres_g = jfnk_gmres_g
   end if
   jfnk_max_gmres_its = 0   ! use default (m) for each test
   !
   if(allocated(s_test)) then
      if(size(s_test) .ne. n) deallocate(s_test, J_s_test, &
     &   F_pred_test)
   end if
   if(.not. allocated(s_test)) allocate(s_test(n), J_s_test(n), &
     &   F_pred_test(n))
   !
   norm_F = dsqrt(dp(n, fx_cur, fx_cur))
   !
   write(luer,'(A)') ' '
   write(luer,'(A)') '=== GMRES Dimension Experiment ==='
   write(luer,'(A)') ' '
   write(luer,'(A,ES14.4)') '  ||F||_w  = ', norm_F
   write(luer,'(A,ES14.4)') '  gres     = ', gres
   write(luer,'(A,ES14.4)') '  del      = ', del
   !
   if(JFNK_RSCALE_READY .and. allocated(JFNK_RSCALE)) then
      rs_min = JFNK_RSCALE(1)
      rs_max = JFNK_RSCALE(1)
      do i = 2, n
         if(JFNK_RSCALE(i) < rs_min) rs_min = JFNK_RSCALE(i)
         if(JFNK_RSCALE(i) > rs_max) rs_max = JFNK_RSCALE(i)
      end do
      rs_ratio = rs_max / max(rs_min, 1.0d-30)
      write(luer,'(A,ES14.4,A,ES14.4,A,ES14.4)') &
         '  RSCALE: min=', rs_min, ' max=', rs_max, ' ratio=', rs_ratio
   end if
   write(luer,'(A)') ' '
   !
   do i_test = 1, n_tests
      m_test = m_tests(i_test)
      !
      if(allocated(v_test)) then
         if(size(v_test,1) .ne. n .or. size(v_test,2) .ne. m_test+1) &
     &      deallocate(v_test, h_test)
      end if
      if(.not. allocated(v_test)) allocate(v_test(n, m_test+1), &
     &   h_test((m_test+1)*m_test))
      !
      ! Restore JVP state for each test
      jfnk_x = x_cur
      jfnk_fx = fx_cur
      s_test = 0.0d0
      !
      ! Run GMRES: 40 total iterations for fair comparison
      its_test = 40
      info_test = 0    ! silent, fresh solve
      res_test = gres
      !
      call gmres_m_v1(m_test, n, s_test, fx_cur, df, dfp, dp, &
                       h_test, v_test, res_test, del, pred_res_test, &
                       its_test, info_test, luer)
      !
      ! Compute J*s_test (restore JVP state first)
      jfnk_x = x_cur
      jfnk_fx = fx_cur
      call df(n, s_test, J_s_test)
      jfnk_total_jvp_evals = jfnk_total_jvp_evals + 1
      !
      ! Diagnostics
      norm_Js = dsqrt(dp(n, J_s_test, J_s_test))
      dir_deriv = dp(n, fx_cur, J_s_test)
      cos_angle = dir_deriv / max(norm_F * norm_Js, 1.0d-30)
      F_pred_test = fx_cur - J_s_test
      norm_Fp = dsqrt(dp(n, F_pred_test, F_pred_test))
      snrm_test = dsqrt(dp(n, s_test, s_test))
      s_inf_test = maxval(abs(s_test))
      alpha_opt = dir_deriv / max(norm_Js * norm_Js, 1.0d-30)
      !
      n_restarts = its_test / m_test
      if(mod(its_test, m_test) > 0) n_restarts = n_restarts + 1
      n_restarts = max(0, n_restarts - 1)
      !
      if(info_test == 0) then
         status_str = 'converged'
      else if(info_test == 1) then
         status_str = 'stagnated'
      else if(info_test == 2) then
         status_str = 'max_its'
      else
         status_str = 'unknown'
      end if
      !
      write(luer,'(A)') ' '
      write(luer,'(A,I3)') '  --- m = ', m_test
      write(luer,'(A,I4,A,A,A,I3)') '  GMRES its=', its_test, &
         '  status=', trim(status_str), '  restarts=', n_restarts
      write(luer,'(A,ES14.4,A,ES14.4)') '  ||s||_w=', snrm_test, &
         '  ||s||_inf=', s_inf_test
      write(luer,'(A,ES14.4,A,ES14.4)') '  ||F-Js||_w=', norm_Fp, &
         '  ||M^-1(F-Js)||=', pred_res_test
      write(luer,'(A,ES14.4,A,ES14.4)') '  <F,Js>_w=', dir_deriv, &
         '  cos(F,Js)=', cos_angle
      write(luer,'(A,ES14.4)') '  alpha_opt=', alpha_opt
      !
      ! Trust-region scaling info
      if(s_inf_test > 1.0d0) then
         scale_fac_test = 1.0d0 / s_inf_test
         snrm_scaled = snrm_test * scale_fac_test
      else
         scale_fac_test = 1.0d0
         snrm_scaled = snrm_test
      end if
      write(luer,'(A,ES14.4,A,ES14.4)') '  del=', del, &
         '  ||s||_w/del=', snrm_test / max(del, 1.0d-30)
      write(luer,'(A,ES14.4,A,ES14.4)') '  inf_scale_fac=', scale_fac_test, &
         '  ||s_scaled||_w=', snrm_scaled
      if(snrm_test < del * 0.9d0) then
         write(luer,'(A)') '  => Trust region NOT active (full Newton step)'
      else
         write(luer,'(A)') '  => Trust region ACTIVE (step constrained)'
      end if
      !
      ! Alpha sweep for m=5 (worst) and m=40 (best)
      if(m_test == 5 .or. m_test == 40) then
         write(luer,'(A)') '  Alpha sweep:'
         alpha_list = (/1.0d0, 0.1d0, 0.01d0, 0.001d0, 0.0001d0/)
         do i = 1, 5
            alpha = alpha_list(i)
            do j = 1, n
               F_pred_test(j) = fx_cur(j) - alpha * J_s_test(j)
            end do
            norm_model = dsqrt(dp(n, F_pred_test, F_pred_test))
            write(luer,'(A,ES10.3,A,ES14.4,A,ES10.3)') &
               '    alpha=', alpha, ' ||F-aJs||=', norm_model, &
               ' ratio=', norm_model / max(norm_F, 1.0d-30)
         end do
      end if
      !
      continue  ! W25: SAVE'd workspace
   end do
   !
   ! Restore state (see entry: the gmres basis pair must describe the
   ! caller's h/v again, not the last test solve's h_test/v_test)
   jfnk_x = save_jfnk_x
   jfnk_fx = save_jfnk_fx
   jfnk_max_gmres_its = save_max_gmres
   jfnk_gmres_jj      = save_gmres_jj
   jfnk_gmres_beta    = save_gmres_beta
   if(have_gmres_zg) then
      if(allocated(jfnk_gmres_z)) deallocate(jfnk_gmres_z, jfnk_gmres_g)
      allocate(jfnk_gmres_z(size(save_gmres_z,1), &
         size(save_gmres_z,2)), jfnk_gmres_g(size(save_gmres_g,1), &
         size(save_gmres_g,2)))
      jfnk_gmres_z = save_gmres_z
      jfnk_gmres_g = save_gmres_g
      continue  ! W25: SAVE'd workspace
   end if
   if(have_gmres_x_cycle_start) then
      if(allocated(jfnk_gmres_x_cycle_start)) &
         deallocate(jfnk_gmres_x_cycle_start)
      allocate(jfnk_gmres_x_cycle_start( &
         size(save_gmres_x_cycle_start)))
      jfnk_gmres_x_cycle_start = save_gmres_x_cycle_start
   else if(allocated(jfnk_gmres_x_cycle_start)) then
      deallocate(jfnk_gmres_x_cycle_start)
   end if
   !
   continue  ! W25: SAVE'd workspace
   !
   write(luer,'(A)') ' '
   write(luer,'(A)') '=== End GMRES Dimension Experiment ==='
   write(luer,'(A)') ' '
   flush(luer)
   !
end subroutine jfnk_gmres_dim_experiment_v1


!----------------------------------------------------------------------
! JFNK_TAUTO_DAMP_V1: FIX_T_AUTO step guard.  Scales the temperature
! components of a Newton/Cauchy direction by jfnk_tauto_alpha and caps
! each at +-jfnk_tauto_dmax.  T components sit at flat index
! i = (k-1)*jfnk_tauto_nt + jfnk_tauto_jt for k in [dst,dend]; all
! other components (populations, Ne) pass through unchanged.  Called
! on the finished physical step AFTER the Krylov solve and before the
! trial eval, so the line search, JVP prediction and commit all see
! the same bounded candidate.
!----------------------------------------------------------------------
subroutine jfnk_tauto_damp_v1(n, s)
   use jfnk_solver_v1_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: s(n)
   integer :: k, i
   !
   if(.not. jfnk_tauto_on) return
   if(jfnk_tauto_nt <= 0 .or. jfnk_tauto_jt <= 0) return
   if(jfnk_tauto_dst <= 0 .or. jfnk_tauto_dend < jfnk_tauto_dst) return
   !
   do k = jfnk_tauto_dst, jfnk_tauto_dend
      i = (k-1)*jfnk_tauto_nt + jfnk_tauto_jt
      if(i >= 1 .and. i <= n) then
         s(i) = s(i) * jfnk_tauto_alpha
         if(s(i) >  jfnk_tauto_dmax) s(i) =  jfnk_tauto_dmax
         if(s(i) < -jfnk_tauto_dmax) s(i) = -jfnk_tauto_dmax
      end if
   end do
   return
end subroutine jfnk_tauto_damp_v1


!----------------------------------------------------------------------
! JFNK_BLK_LSQ_SOLVE: apply one LSQ-factored block to an np-vector.
! Replicates the JFNK_PRECOND_V1 'LSQ' apply on a single np x np block:
!   w = rsf^2 * b
!   u = csf * (B_raw^T w + ridge*w + reg-row substitution correction)
!   solve (B_eq^T B_eq + lam^2 I) u = rhs via the DPOTRF factor chol
!   z = csf * u
! Shared by the production SCHUR preconditioner and the unit-test
! driver so both exercise identical arithmetic.
! info_out: 0 on success, LAPACK info otherwise.
!----------------------------------------------------------------------
subroutine jfnk_blk_lsq_solve(np, chol, braw, ridge, reg, rsf, csf, &
     b, z, info_out)
   implicit none
   integer, intent(in) :: np
   double precision, intent(in) :: chol(np,np), braw(np,np)
   double precision, intent(in) :: ridge(np), rsf(np), csf(np), b(np)
   logical, intent(in) :: reg(np)
   double precision, intent(out) :: z(np)
   integer, intent(out) :: info_out
   integer :: i, j
   double precision :: w(np), u(np)
   !
   do i = 1, np
      w(i) = b(i) * rsf(i) * rsf(i)
   end do
   do j = 1, np
      u(j) = 0d0
      do i = 1, np
         u(j) = u(j) + braw(i,j) * w(i)
      end do
      u(j) = u(j) + ridge(j) * w(j)
   end do
   do i = 1, np
      if(reg(i)) then
         do j = 1, np
            u(j) = u(j) - w(i) * braw(i,j)
         end do
         u(i) = u(i) + w(i) - w(i) * ridge(i)
      end if
   end do
   do j = 1, np
      u(j) = u(j) * csf(j)
   end do
   call dpotrs('L', np, 1, chol, np, u, np, info_out)
   if(info_out /= 0) then
      z = 0d0
      return
   end if
   do j = 1, np
      z(j) = u(j) * csf(j)
   end do
   return
end subroutine jfnk_blk_lsq_solve


!----------------------------------------------------------------------
! JFNK_SCHUR_FORM_V1: assemble and factor the dense temperature Schur
! complement
!   S(k,k) = JTT(k) - JTn(k) . Jnn(k)^{-1} JnT(k)
! where Jnn^{-1} JnT is evaluated with the same per-depth LSQ block
! solve used at apply time (chol/braw/ridge/reg/rsf/csf as above).
! np_ = nt-2 population rows per depth; T sits at flat row nt,
! Ne at nt-1.  S is nd x nd; factored in place by DGETRF.
!----------------------------------------------------------------------
subroutine jfnk_schur_form_v1(nd, nt, np_, jt, chol, braw, ridge, &
     reg, rsf, csf, s, sipiv, info)
   implicit none
   integer, intent(in) :: nd, nt, np_, jt
   double precision, intent(in) :: chol(nt,nt,nd), braw(nt,nt,nd)
   double precision, intent(in) :: ridge(nt,nd), rsf(nt,nd), csf(nt,nd)
   logical, intent(in) :: reg(nt,nd)
   double precision, intent(out) :: s(nd,nd)
   integer, intent(out) :: sipiv(nd), info
   integer :: k, i, ierr
   double precision :: w(np_), b(np_), jtt, chl(np_,np_)
   !
   s = 0d0
   info = 0
   do k = 1, nd
      do i = 1, np_
         b(i) = braw(i,jt,k)
      end do
      chl = chol(1:np_,1:np_,k)
      call jfnk_blk_lsq_solve(np_, chl, braw(1:np_,1:np_,k), &
           ridge(1:np_,k), reg(1:np_,k), rsf(1:np_,k), csf(1:np_,k), &
           b, w, ierr)
      if(ierr /= 0) w = 0d0
      jtt = braw(jt,jt,k)
      s(k,k) = jtt
      do i = 1, np_
         s(k,k) = s(k,k) - braw(jt,i,k) * w(i)
      end do
   end do
   call dgetrf(nd, nd, s, nd, sipiv, info)
   return
end subroutine jfnk_schur_form_v1


!----------------------------------------------------------------------
! JFNK_SCHUR_APPLY_V1: two-level (n,T) Schur-complement preconditioner
! apply z = M^{-1} v on the flat state vector v(nd*nt).
!
!   z_n(k) = Jnn(k)^{-1} v_n(k)          (LSQ block solve, pops only)
!   t(k)   = v_T(k) - JTn(k) . z_n(k)
!   z_T    = S^{-1} t                    (dense nd x nd DGETRS)
!   z_n(k) = z_n(k) - Jnn(k)^{-1} (JnT(k) z_T(k))
!
! T is written to row nt of each depth block, Ne (nt-1) zeroed; reg
! population rows are not projected here -- the caller's post-chain
! does that, exactly as for the LU/LSQ paths.
!----------------------------------------------------------------------
subroutine jfnk_schur_apply_v1(nd, nt, np_, jt, jne, chol, braw, &
     ridge, reg, rsf, csf, s, sipiv, v, info)
   implicit none
   integer, intent(in) :: nd, nt, np_, jt, jne
   double precision, intent(in) :: chol(nt,nt,nd), braw(nt,nt,nd)
   double precision, intent(in) :: ridge(nt,nd), rsf(nt,nd), csf(nt,nd)
   logical, intent(in) :: reg(nt,nd)
   double precision, intent(in) :: s(nd,nd)
   integer, intent(in) :: sipiv(nd)
   double precision, intent(inout) :: v(nd*nt)
   integer, intent(out) :: info
   integer :: k, i, b, ierr
   double precision :: zn(np_,nd), cc(np_), trhs(nd), zt(nd), bb(np_)
   double precision :: chl(np_,np_)
   !
   info = 0
   do k = 1, nd
      b = (k-1) * nt
      chl = chol(1:np_,1:np_,k)
      call jfnk_blk_lsq_solve(np_, chl, braw(1:np_,1:np_,k), &
           ridge(1:np_,k), reg(1:np_,k), rsf(1:np_,k), csf(1:np_,k), &
           v(b+1:b+np_), zn(1,k), ierr)
      if(ierr /= 0) then
         info = ierr
         return
      end if
      trhs(k) = v(b+jt)
      do i = 1, np_
         trhs(k) = trhs(k) - braw(jt,i,k) * zn(i,k)
      end do
   end do
   zt = trhs
   call dgetrs('N', nd, 1, s, nd, sipiv, zt, nd, ierr)
   if(ierr /= 0) then
      info = ierr
      return
   end if
   do k = 1, nd
      b = (k-1) * nt
      do i = 1, np_
         bb(i) = braw(i,jt,k) * zt(k)
      end do
      chl = chol(1:np_,1:np_,k)
      call jfnk_blk_lsq_solve(np_, chl, braw(1:np_,1:np_,k), &
           ridge(1:np_,k), reg(1:np_,k), rsf(1:np_,k), csf(1:np_,k), &
           bb, cc, ierr)
      if(ierr /= 0) then
         info = ierr
         return
      end if
      do i = 1, np_
         v(b+i) = zn(i,k) - cc(i)
      end do
      if(jne >= 1 .and. jne <= nt) v(b+jne) = 0d0
      v(b+jt) = zt(k)
   end do
   return
end subroutine jfnk_schur_apply_v1

! ----------------------------------------------------------------------
! jfnk_pc_factor_block_v1: factorization of ONE depth block of the BA
! preconditioner, extracted verbatim from the per-depth loop of
! JFNK_FACTOR_PRECOND_V1 (jfnk_cmfgen_interface_v1.f) so the unit tests
! run the identical production code (precedent: jfnk_schur_apply_v1).
!
! a(n,n) arrives already carrying the PTC fold and the frozen T/Ne
! identity substitution: those write the persistent JFNK_BA_A /
! JFNK_BA_DIAG which the banded path, the Jacobi fallback and the
! admission gates also read, so they stay in the caller.  pre_reg is
! the caller-built population-floor mask (needs POPS_BASE/ATM); this
! routine treats its rows exactly like detected near-null rows.
!
! All jfnk_pc_* module flags default to the historical behaviour; each
! flag branch is a VADAT-gated audit fix (JFNK_PRECOND_AUDIT.md).
! ----------------------------------------------------------------------
   subroutine jfnk_pc_factor_block_v1(a, n, pre_reg, row_tiny, &
        ridge_piv, lsq_lam, lsq_mode, luer, &
        lu, ipiv, rsf, csf, reg, reg_tiny, ridge_rows, &
        rcond, pivrat, itry_used, ridge_used, lu_ok)
      use jfnk_solver_v1_mod, only: jfnk_ridge_floor, &
           jfnk_pc_signed, jfnk_pc_reg_zerocol, jfnk_pc_tiny_equil, &
           jfnk_rcond_gate, jfnk_rcond_min, jfnk_itry_cap
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: n, luer
      double precision, intent(in) :: a(n,n)
      logical, intent(in) :: pre_reg(n), lsq_mode
      double precision, intent(in) :: row_tiny, ridge_piv, lsq_lam
      double precision, intent(out) :: lu(n,n), rsf(n), csf(n), &
           ridge_rows(n), rcond, pivrat, ridge_used
      integer, intent(out) :: ipiv(n), itry_used
      logical, intent(out) :: reg(n), reg_tiny(n), lu_ok
!
      integer :: i, j, itry, ios, info_k, info_g, itry_max
      double precision :: eps, ridge_add, blkmax, scale_k, rowmax_i
      double precision :: row_cnd, col_cnd, amax, lsq_l2, anorm_k
      double precision :: s_tmp, piv_min, piv_max, rcond_k
      double precision, allocatable :: mtst(:,:), b_tst(:), x_tst(:)
      double precision, allocatable :: gn_work(:,:)
      double precision, allocatable :: gecon_work(:)
      integer, allocatable :: gecon_iwork(:)
      logical, allocatable :: reg_raw(:), reg_eq(:)
      double precision, allocatable :: csf_pre(:)
      integer :: n_tiny_chg
!
      lu_ok = .false.
      ridge_used = 0.0d0
      itry_used = -1
      rcond = -1.0d0
      pivrat = -1.0d0
      ridge_rows(1:n) = 0.0d0
      reg(1:n) = .false.
      reg_tiny(1:n) = .false.
!
      scale_k = 0.0d0
      do i = 1, n
         scale_k = max(scale_k, abs(a(i,i)))
      end do
      scale_k = max(scale_k, 1.0d0)
!
      blkmax = 0.0d0
      do i = 1, n
         do j = 1, n
            blkmax = max(blkmax, abs(a(i,j)))
         end do
      end do
!
      allocate(mtst(n,n), b_tst(n), x_tst(n), stat=ios)
      if(ios /= 0) return
      if(lsq_mode)then
         allocate(gn_work(n,n), stat=ios)
         if(ios /= 0) return
      end if
!
! The rcond gate (F5) replaces the x=1 forward-error micro-test by an
! explicit reciprocal-condition threshold and caps the retry ladder;
! the default path keeps the historical 0..8 ladder unchanged.
!
      if(jfnk_rcond_gate)then
         itry_max = max(jfnk_itry_cap, 1) - 1
      else
         itry_max = 8
      end if
!
      do itry = 0, itry_max
         reg(1:n) = .false.
         reg_tiny(1:n) = .false.
         if(itry == 0)then
            eps = jfnk_ridge_floor
         else if(itry == 1)then
            eps = 1.0d-16
         else if(itry == 2)then
            eps = 1.0d-14
         else if(itry == 3)then
            eps = 1.0d-12
         else if(itry == 4)then
            eps = 1.0d-10
         else if(itry == 5)then
            eps = 1.0d-8
         else if(itry == 6)then
            eps = 1.0d-6
         else if(itry == 7)then
            eps = 1.0d-4
         else
            eps = 1.0d-2
         end if
         ridge_add = max(eps, jfnk_ridge_floor)
         lu(1:n,1:n) = a(1:n,1:n)
         if(jfnk_pc_signed)then
! [JFNK_PTC_SIGNED]: the ridge must increase |d|, so it follows the
! sign of the post-fold diagonal (d_i*(1+ridge) form).  A positive
! ridge added to a loss-dominated d<0 row instead shrinks |d| and can
! push the pivot through zero.  For |d_i| below the 1e-30*scale_k
! floor the sign is ambiguous by definition; sign(x,0)=x keeps +ridge
! on an exactly-zero diagonal, matching the unsigned path there.
            do i = 1, n
               lu(i,i) = lu(i,i) + ridge_add * &
                    sign(max(abs(a(i,i)), 1.0d-30*scale_k), a(i,i))
            end do
         else
            do i = 1, n
               lu(i,i) = lu(i,i) + ridge_add * &
                    max(abs(a(i,i)), 1.0d-30*scale_k)
            end do
         end if
!
! Near-null row detection.  Default (historical): test the raw row
! against the block maximum, which is dominated by the stiffest rows
! so a full-rank row in different units can be misclassified.
! [JFNK_ROW_TINY_EQUIL] (F3): scale-free test -- a row counts as
! null only when every entry is negligible relative to its column
! maximum, share_i = max_j |a_ij|/colmax_j < row_tiny.  A row in
! small units that still owns its column keeps share ~1; a genuinely
! null row keeps share ~0.  DGEEQU factors cannot express this: its
! column scaling normalizes the row-scaled matrix, leaving a small
! row's absolute magnitude dominant anyway.
!
         if(jfnk_pc_tiny_equil)then
            allocate(reg_raw(n), reg_eq(n), csf_pre(n), stat=ios)
            if(ios == 0)then
               do j = 1, n
                  csf_pre(j) = 0.0d0
                  do i = 1, n
                     csf_pre(j) = max(csf_pre(j), abs(a(i,j)))
                  end do
               end do
               do i = 1, n
                  rowmax_i = 0.0d0
                  do j = 1, n
                     if(csf_pre(j) > 0.0d0) rowmax_i = &
                          max(rowmax_i, abs(a(i,j))/csf_pre(j))
                  end do
                  reg_eq(i) = (rowmax_i < row_tiny)
                  rowmax_i = 0.0d0
                  do j = 1, n
                     rowmax_i = max(rowmax_i, abs(a(i,j)))
                  end do
                  reg_raw(i) = (rowmax_i < row_tiny*blkmax)
               end do
               n_tiny_chg = count(reg_raw .neqv. reg_eq)
               if(n_tiny_chg > 0 .and. luer > 0)then
                  write(luer,'(A,I6,A)') &
                    '  JFNK: ROW_TINY_EQUIL reclassifies ', &
                    n_tiny_chg, ' rows vs raw test'
               end if
               reg_tiny(1:n) = reg_eq(1:n)
               deallocate(reg_raw, reg_eq, csf_pre)
            else
! allocation failure must not silently disable regularization;
! fall back to the historical raw-block test.
               if(blkmax > 0.0d0)then
                  do i = 1, n
                     rowmax_i = 0.0d0
                     do j = 1, n
                        rowmax_i = max(rowmax_i, abs(a(i,j)))
                     end do
                     if(rowmax_i < row_tiny*blkmax) &
                          reg_tiny(i) = .true.
                  end do
               end if
            end if
         else
            if(blkmax > 0.0d0)then
               do i = 1, n
                  rowmax_i = 0.0d0
                  do j = 1, n
                     rowmax_i = max(rowmax_i, abs(a(i,j)))
                  end do
                  if(rowmax_i < row_tiny*blkmax) reg_tiny(i) = .true.
               end do
            end if
         end if
!
! Substitute identity for regularized rows (and, under
! [JFNK_REG_ZERO_COL], columns) in the factored block only -- the
! residual and JVP never see the substitution.  Column zeroing
! removes the stale coupling of the regularized unknown into the
! other equations; without it the projected solve error grows
! linearly with the column magnitude (unit test T-regcol).
!
         do i = 1, n
            if(reg_tiny(i) .or. pre_reg(i))then
               reg(i) = .true.
               do j = 1, n
                  lu(i,j) = 0.0d0
               end do
               if(jfnk_pc_reg_zerocol)then
                  do j = 1, n
                     lu(j,i) = 0.0d0
                  end do
               end if
               lu(i,i) = 1.0d0
            end if
         end do
!
         call dgeequ(n, n, lu, n, rsf, csf, row_cnd, col_cnd, amax, &
              info_k)
         if(info_k /= 0) cycle
         do j = 1, n
            do i = 1, n
               lu(i,j) = lu(i,j) * rsf(i) * csf(j)
            end do
         end do
!
         if(lsq_mode)then
            lsq_l2 = (lsq_lam * max(amax, 1.0d0))**2
            call dsyrk('L', 'T', n, n, 1.0d0, lu, n, 0.0d0, &
                 gn_work, n)
            do i = 1, n
               gn_work(i,i) = gn_work(i,i) + lsq_l2
            end do
            lu(1:n,1:n) = gn_work(1:n,1:n)
            mtst(1:n,1:n) = lu(1:n,1:n)
            call dpotrf('L', n, lu, n, info_k)
         else
! F4: DGECON below is called with NORM='I' (infinity norm), so ANORM
! must be the maximum ROW sum; the pre-fix code supplied the maximum
! column sum (1-norm), which on unsymmetric blocks underestimates
! rcond by an order-dependent factor (~9x at n=50, ~240x at n=800 in
! the unit test).  rcond is diagnostic only today, so this changes
! logged values, not control flow.
            anorm_k = 0.0d0
            do i = 1, n
               s_tmp = 0.0d0
               do j = 1, n
                  s_tmp = s_tmp + abs(lu(i,j))
               end do
               anorm_k = max(anorm_k, s_tmp)
            end do
            mtst(1:n,1:n) = lu(1:n,1:n)
            call dgetrf(n, n, lu, n, ipiv, info_k)
         end if
         if(info_k == 0)then
!
! Acceptance test.  Default: synthetic-RHS forward-error micro-test
! (b = M*1 must return x=1 to 1e-4).  [JFNK_PC_RCOND_GATE] (F5): the
! micro-test is a hidden condition-number gate (~1e13) since the true
! backward error stays ~1e-17; the gate instead compares the DGECON
! reciprocal condition estimate against JFNK_PC_RCOND_MIN directly.
!
            if(jfnk_rcond_gate .and. .not. lsq_mode)then
               if(.not. allocated(gecon_work))then
                  allocate(gecon_work(4*n), gecon_iwork(n), stat=ios)
                  if(ios /= 0) return
               end if
               call dgecon('I', n, lu, n, anorm_k, rcond_k, &
                    gecon_work, gecon_iwork, info_g)
               if(info_g == 0) rcond = rcond_k
               ! Fail closed on a non-finite estimate (10-Nov-2026
               ! audit): a NaN rcond makes '< min' .false. and would
               ! accept a NaN-contaminated LU as a valid factorization.
               if(info_g /= 0 .or. .not. ieee_is_finite(rcond_k) .or. &
                    rcond_k < jfnk_rcond_min)then
                  ridge_used = 0.0d0
                  cycle
               end if
            else
               b_tst(1:n) = 0.0d0
               do i = 1, n
                  do j = 1, n
                     b_tst(i) = b_tst(i) + mtst(i,j)
                  end do
               end do
               x_tst(1:n) = b_tst(1:n)
               if(lsq_mode)then
                  call dpotrs('L', n, 1, lu, n, x_tst, n, info_g)
               else
                  call dgetrs('N', n, 1, lu, n, ipiv, x_tst, n, info_g)
               end if
               if(info_g /= 0 .or. &
                    maxval(abs(x_tst(1:n)-1.0d0)) > 1.0d-4)then
                  ridge_used = 0.0d0
                  cycle
               end if
               if(lsq_mode)then
                  rcond = -1.0d0
               else
                  if(.not. allocated(gecon_work))then
                     allocate(gecon_work(4*n), gecon_iwork(n), stat=ios)
                     if(ios /= 0) return
                  end if
                  call dgecon('I', n, lu, n, anorm_k, rcond_k, &
                       gecon_work, gecon_iwork, info_g)
                  if(info_g == 0) rcond = rcond_k
               end if
            end if
            lu_ok = .true.
            ridge_used = ridge_add
            itry_used = itry
            piv_min = huge(1.0d0)
            piv_max = 0.0d0
            do i = 1, n
               piv_min = min(piv_min, abs(lu(i,i)))
               piv_max = max(piv_max, abs(lu(i,i)))
            end do
            pivrat = piv_min / max(piv_max, 1.0d-300)
            if(.not. lsq_mode .and. &
                 ridge_piv > 0.0d0 .and. pivrat < ridge_piv)then
! DGETRF can return INFO=0 on a near-singular block whose smallest
! pivot still avoids an exact zero; the solve inherits the full
! amplification.  With JFNK_RIDGE_PIV>0 the successful-but-ill-
! conditioned factorization is rejected so the ladder escalates.
               lu_ok = .false.
               ridge_used = 0.0d0
               cycle
            end if
            if(jfnk_pc_signed)then
               do i = 1, n
                  ridge_rows(i) = ridge_add * &
                       sign(max(abs(a(i,i)), 1.0d-30*scale_k), a(i,i))
               end do
            else
               do i = 1, n
                  ridge_rows(i) = ridge_add * &
                       max(abs(a(i,i)), 1.0d-30*scale_k)
               end do
            end if
            exit
         end if
      end do
!
      if(allocated(gecon_work)) deallocate(gecon_work, gecon_iwork)
      if(allocated(gn_work)) deallocate(gn_work)
      deallocate(mtst, b_tst, x_tst)
      return
   end subroutine jfnk_pc_factor_block_v1
!
! ----------------------------------------------------------------------
! jfnk_pc_apply_block_v1: apply of ONE depth block, extracted verbatim
! from JFNK_PRECOND_V1 so the unit tests run the identical production
! code.  v(n) is the input right-hand side on entry and the
! preconditioned output on return; rhs(n) is an unmodified copy of the
! input used by the Jacobi fallback.  fixed(n) is the caller-built
! fixed-unknown mask (frozen T/Ne rows), reg(n) the factorization's
! regularized-row mask; both project z_i = 0, which is linear and keeps
! M^{-1} a consistent operator for GMRES.  Nonfinite outputs are zeroed
! and counted in nzclip -- the only nonlinearity, and the caller aborts
! the GMRES cycle when nzclip grows (jfnk_pc_clip_event).
!
! The LSQ B^T w expansion assumes the row-only regularization of the
! default path; under [JFNK_REG_ZERO_COL] the substituted block also
! has zeroed columns, for which (B'^T w)_i = w_i exactly on reg rows.
! ----------------------------------------------------------------------
   subroutine jfnk_pc_apply_block_v1(v, rhs, n, lu_ok, lu, ipiv, &
        rsf, csf, a_raw, ridge_rows, ba_diag, reg, fixed, lsq_mode, &
        zmax, nproj, nzero, nzclip, nzmax_hit, nsolved)
      use jfnk_solver_v1_mod, only: jfnk_pc_reg_zerocol
      implicit none
      integer, intent(in) :: n
      double precision, intent(inout) :: v(n)
      double precision, intent(in) :: rhs(n)
      logical, intent(in) :: lu_ok, lsq_mode, reg(n), fixed(n)
      double precision, intent(in) :: lu(n,n), rsf(n), csf(n), &
           a_raw(n,n), ridge_rows(n), ba_diag(n), zmax
      integer, intent(in) :: ipiv(n)
      integer, intent(inout) :: nproj, nzero, nzclip, nzmax_hit, &
           nsolved
!
      integer :: i, j, info
      double precision :: diag_use, wi_lsq
      double precision, allocatable :: w_lsq(:)
!
      if(.not. lu_ok)then
         do j = 1, n
            diag_use = ba_diag(j) + ridge_rows(j)
            if(reg(j))then
               v(j) = 0.0d0
               nproj = nproj + 1
            else if(abs(diag_use) > 1.0d-30)then
               v(j) = rhs(j) / diag_use
            else
               v(j) = 0.0d0
               nzero = nzero + 1
            end if
            nsolved = nsolved + 1
         end do
         return
      end if
!
      do j = 1, n
         v(j) = v(j) * rsf(j)
      end do
      if(lsq_mode)then
         allocate(w_lsq(n))
! w = ROW_SF v' and rhs = B^T w.  B = A + ridge + reg-row substitution
! expanded as B^T w = A^T w + D_r w + sum_{i reg} w_i (e_i - A(i,:)^T)
! - w_i D_r(i) e_i so no second n x n block is needed.
         do j = 1, n
            w_lsq(j) = v(j) * rsf(j)
         end do
         call dgemv('T', n, n, 1.0d0, a_raw, n, w_lsq, 1, 0.0d0, &
              v, 1)
         do j = 1, n
            v(j) = v(j) + ridge_rows(j) * w_lsq(j)
         end do
         do i = 1, n
            if(reg(i))then
               wi_lsq = w_lsq(i)
               do j = 1, n
                  v(j) = v(j) - wi_lsq * a_raw(i,j)
               end do
               if(jfnk_pc_reg_zerocol)then
! With the column also zeroed, (B'^T w)_i = w_i exactly: the raw
! DGEMV left the other rows' column-i contributions in component i.
                  v(i) = wi_lsq
               else
                  v(i) = v(i) + wi_lsq - wi_lsq * ridge_rows(i)
               end if
            end if
         end do
         do j = 1, n
            v(j) = v(j) * csf(j)
         end do
         call dpotrs('L', n, 1, lu, n, v, n, info)
      else
         call dgetrs('N', n, 1, lu, n, ipiv, v, n, info)
      end if
      if(info /= 0)then
         do j = 1, n
            diag_use = ba_diag(j) + ridge_rows(j)
            if(reg(j))then
               v(j) = 0.0d0
               nproj = nproj + 1
            else if(abs(diag_use) > 1.0d-30)then
               v(j) = rhs(j) / diag_use
            else
               v(j) = 0.0d0
               nzero = nzero + 1
            end if
         end do
      else
         do j = 1, n
            v(j) = v(j) * csf(j)
         end do
         do j = 1, n
            if(reg(j))then
               v(j) = 0.0d0
               nproj = nproj + 1
            else if(fixed(j))then
               v(j) = 0.0d0
            else if(.not. (abs(v(j)) <= huge(1.0d0)))then
               v(j) = 0.0d0
               nzero = nzero + 1
               nzclip = nzclip + 1
            else if(abs(v(j)) > zmax)then
               nzmax_hit = nzmax_hit + 1
            end if
            nsolved = nsolved + 1
         end do
      end if
      if(allocated(w_lsq)) deallocate(w_lsq)
      return
   end subroutine jfnk_pc_apply_block_v1
!
! ----------------------------------------------------------------------
! jfnk_rscale_branch_v1: JFNK_RSCALE mode dispatch, extracted as a
! pure function so the unit tests exercise the exact production
! branch order.  F7: the pre-fix chain tested mode(1:4)=='DIAG'
! before 'DIAG_RDG', which made DIAG_RDG unreachable.  The specific
! 8-character match must be tested before the 4-character prefix.
! Return codes: 0 EUCLIDEAN (default), 1 SQRT_STE, 2 DIAG, 3 DIAG_RDG.
! ----------------------------------------------------------------------
   integer function jfnk_rscale_branch_v1(mode)
      implicit none
      character(len=*), intent(in) :: mode
      character(len=8) :: m8
      m8 = mode
      if(m8 == 'SQRT_STE')then
         jfnk_rscale_branch_v1 = 1
      else if(m8 == 'DIAG_RDG')then
         jfnk_rscale_branch_v1 = 3
      else if(m8(1:4) == 'DIAG')then
         jfnk_rscale_branch_v1 = 2
      else
         jfnk_rscale_branch_v1 = 0
      end if
      return
   end function jfnk_rscale_branch_v1
!
! ----------------------------------------------------------------------
! jfnk_pc_sign_vec_v1: per-row sign vector for [JFNK_PTC_SIGNED].
! sgn_i = sign(d_i) when |d_i| clears the floor (the caller supplies
! the ROW_TINY-style floor), else +1: a numerically zero diagonal
! cannot decide a sign, and such rows are regularized anyway, so +1
! reproduces the unsigned shift on them.  Shared by the residual PTC
! term, the factorization fold and the apply rescale so all three use
! the same convention.
! ----------------------------------------------------------------------
   subroutine jfnk_pc_sign_vec_v1(d, n, floor_abs, sgn)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: d(n), floor_abs
      double precision, intent(out) :: sgn(n)
      integer :: i
      do i = 1, n
         if(abs(d(i)) >= floor_abs)then
            sgn(i) = sign(1.0d0, d(i))
         else
            sgn(i) = 1.0d0
         end if
      end do
      return
   end subroutine jfnk_pc_sign_vec_v1
!
! ----------------------------------------------------------------------
! jfnk_pc_sign_stats_v1: [JFNK_PC_SIGN_DIAG] report data.  For rows
! with d < 0 and sigma > 0 (the unsigned-shift cancellation class),
! collect the ratio sigma/|d| and count the rows where the unsigned
! fold nearly annihilates the diagonal (|d+sigma| < 0.1|d|).  The
! caller prints; this routine only reduces the arrays.
! ----------------------------------------------------------------------
   subroutine jfnk_pc_sign_stats_v1(d, sig, n, n_neg, n_cancel, &
        rmin, rmed, rp99)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: d(n), sig(n)
      integer, intent(out) :: n_neg, n_cancel
      double precision, intent(out) :: rmin, rmed, rp99
      integer :: i, m
      double precision, allocatable :: r(:)
!
      n_neg = 0
      n_cancel = 0
      allocate(r(n))
      do i = 1, n
         if(d(i) < 0.0d0 .and. sig(i) > 0.0d0)then
            n_neg = n_neg + 1
            r(n_neg) = sig(i) / abs(d(i))
            if(abs(d(i)+sig(i)) < 0.1d0*abs(d(i))) &
                 n_cancel = n_cancel + 1
         end if
      end do
      rmin = 0.0d0; rmed = 0.0d0; rp99 = 0.0d0
      if(n_neg > 0)then
         call jfnk_pc_hpsort(r, n_neg)
         rmin = r(1)
         m = max(1, int(0.5d0*n_neg + 0.5d0))
         rmed = r(min(m, n_neg))
         m = max(1, int(0.99d0*n_neg + 0.5d0))
         rp99 = r(min(m, n_neg))
      end if
      deallocate(r)
      return
   end subroutine jfnk_pc_sign_stats_v1
!
! Heapsort for the quantile report (LAPACK has no sort).
   subroutine jfnk_pc_hpsort(a, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(inout) :: a(n)
      integer :: i, j, l, ir
      double precision :: t
      if(n < 2) return
      l = n/2 + 1
      ir = n
      do
         if(l > 1)then
            l = l - 1
            t = a(l)
         else
            t = a(ir)
            a(ir) = a(1)
            ir = ir - 1
            if(ir == 1)then
               a(1) = t
               return
            end if
         end if
         i = l
         j = 2*l
         do while(j <= ir)
            if(j < ir)then
               if(a(j) < a(j+1)) j = j + 1
            end if
            if(t < a(j))then
               a(i) = a(j)
               i = j
               j = 2*j
            else
               j = ir + 1
            end if
         end do
         a(i) = t
      end do
   end subroutine jfnk_pc_hpsort


!----------------------------------------------------------------------
! JFNK_PC_FINAL_GUARD_V1: last-chance finite-output scan on the fully
! post-processed preconditioner result (10-Nov-2026 audit).  The
! per-depth apply guards its own intermediate values, but the
! log-unknown rescale, PTC diagonal rescale, T-ladder scaling and the
! enrichment backsolve run afterwards and can still overflow or return
! NaN on a block the factor stage accepted.  Non-finite entries are
! zeroed -- the only safe value inside a matvec -- and counted so the
! caller can raise jfnk_pc_clip_event and abort the Krylov solve
! rather than silently deliver a nonlinear M^{-1}.
! Extracted as a pure routine so the unit tests run the identical
! production code.
!----------------------------------------------------------------------
   subroutine jfnk_pc_final_guard_v1(n, v, nbad)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      integer, intent(in) :: n
      double precision, intent(inout) :: v(n)
      integer, intent(out) :: nbad
      integer :: i
      nbad = 0
      do i = 1, n
         if(.not. ieee_is_finite(v(i))) then
            v(i) = 0.0d0
            nbad = nbad + 1
         end if
      end do
   end subroutine jfnk_pc_final_guard_v1


!----------------------------------------------------------------------
! JFNK_H_SELECT_V1: plateau selector for [JFNK_HADAPT] (10-Nov-2026
! audit).  jv_sw holds the JVPs sampled by the h-sweep (one column per
! multiplier; ok_sw marks evals that survived the residual checks).
! A forward difference inside the stable region between truncation
! error (large h) and residual noise (small h) changes little across a
! decade; noise-dominated points fluctuate randomly and truncation-
! dominated points drift systematically, so the best adjacent-pair
! agreement locates the plateau.  imult returns the upper member of
! the best pair -- the more noise-immune endpoint, with agreement
! bounding the truncation drift -- or 0 when no adjacent pair agrees
! within tol, meaning the direction is not FD-resolvable at any
! scanned step.
!----------------------------------------------------------------------
   subroutine jfnk_h_select_v1(jv_sw, ok_sw, n, nh, tol, &
        imult, rel_min)
      implicit none
      integer, intent(in) :: n, nh
      double precision, intent(in) :: jv_sw(n, nh), tol
      logical, intent(in) :: ok_sw(nh)
      integer, intent(out) :: imult
      double precision, intent(out) :: rel_min
      integer :: i, k
      double precision :: rel_k, na, nb
      imult = 0
      rel_min = huge(1.0d0)
      do k = 2, nh
         if(.not. ok_sw(k) .or. .not. ok_sw(k-1)) cycle
         na = 0.0d0
         nb = 0.0d0
         rel_k = 0.0d0
         do i = 1, n
            na = max(na, abs(jv_sw(i,k-1)))
            nb = max(nb, abs(jv_sw(i,k)))
            rel_k = max(rel_k, abs(jv_sw(i,k) - jv_sw(i,k-1)))
         end do
         rel_k = rel_k / max(max(na, nb), tiny(1.0d0))
         if(rel_k < rel_min)then
            rel_min = rel_k
            imult = k
         end if
      end do
      if(rel_min > tol) imult = 0
   end subroutine jfnk_h_select_v1


!----------------------------------------------------------------------
! JFNK_BTD_FACTOR_V1: block-Thomas factorization of a scaled
! block-tridiagonal matrix for [JFNK_BANDED_TDMA].  The DGBTRF path on
! the assembled banded matrix costs O(N*KL*KU) ~ 3E15 flops at
! NT~2E3/ND~70 -- infeasible per JFNK activation.  The block
! recurrence costs ND*(4/3)*NT^3 (~seconds once per activation) and
! each apply is O(ND*NT^2).  lu_out arrives holding the diagonal
! blocks D_k and returns the pivoted LU factors of the modified
! diagonals D~_k = D_k - al_k c_{k-1}; c_out arrives holding the
! uppers U_k and returns c_k = D~_k^{-1} U_k; al holds the lowers
! (read-only).  nfail returns the failing depth index, 0 on success
! -- the caller retries with a ridge or falls back to the per-depth
! LU.
!----------------------------------------------------------------------
   subroutine jfnk_btd_factor_v1(nt, nd, al, lu_out, c_out, &
        ipiv, nfail, luer)
      implicit none
      integer, intent(in) :: nt, nd, luer
      double precision, intent(in) :: al(nt,nt,nd)
      double precision, intent(inout) :: lu_out(nt,nt,nd), &
           c_out(nt,nt,nd)
      integer, intent(out) :: ipiv(nt,nd), nfail
      integer :: k, info
      nfail = 0
      do k = 1, nd
         if(k > 1)then
            call dgemm('N','N',nt,nt,nt,-1.0d0,al(1,1,k),nt, &
                 c_out(1,1,k-1),nt,1.0d0,lu_out(1,1,k),nt)
         end if
         call dgetrf(nt,nt,lu_out(1,1,k),nt,ipiv(1,k),info)
         if(info /= 0)then
            nfail = k
            if(luer > 0) write(luer,'(A,I6,A,I4)') &
                 '  JFNK-BTD: dgetrf failed at depth ',k,' info=',info
            return
         end if
         if(k < nd)then
            call dgetrs('N',nt,nt,lu_out(1,1,k),nt,ipiv(1,k), &
                 c_out(1,1,k),nt,info)
            if(info /= 0)then
               nfail = k
               if(luer > 0) write(luer,'(A,I6,A,I4)') &
                    '  JFNK-BTD: dgetrs failed at depth ',k, &
                    ' info=',info
               return
            end if
         end if
      end do
   end subroutine jfnk_btd_factor_v1


!----------------------------------------------------------------------
! JFNK_BTD_APPLY_V1: in-place block-Thomas solve M x = b using the
! factors from jfnk_btd_factor_v1.  Forward sweep overwrites v_k with
! D~_k^{-1}(v_k - L_k v_{k-1}); backward sweep completes
! v_k = w_k - c_k v_{k+1}.  In-place is safe because each block is
! consumed only by its immediate neighbour in sweep order.
!----------------------------------------------------------------------
   subroutine jfnk_btd_apply_v1(nt, nd, al, lu_f, c_f, ipiv, v, &
        nfail)
      implicit none
      integer, intent(in) :: nt, nd
      double precision, intent(in) :: al(nt,nt,nd), lu_f(nt,nt,nd), &
           c_f(nt,nt,nd)
      integer, intent(in) :: ipiv(nt,nd)
      double precision, intent(inout) :: v(nt,nd)
      integer, intent(out) :: nfail
      integer :: k, info
      nfail = 0
      do k = 1, nd
         if(k > 1)then
            call dgemv('N',nt,nt,-1.0d0,al(1,1,k),nt,v(1,k-1),1, &
                 1.0d0,v(1,k),1)
         end if
         call dgetrs('N',nt,1,lu_f(1,1,k),nt,ipiv(1,k),v(1,k),nt, &
              info)
         if(info /= 0)then
            nfail = k
            return
         end if
      end do
      do k = nd-1, 1, -1
         call dgemv('N',nt,nt,-1.0d0,c_f(1,1,k),nt,v(1,k+1),1, &
              1.0d0,v(1,k),1)
      end do
   end subroutine jfnk_btd_apply_v1
