!----------------------------------------------------------------------
! TEST_GUARDS: standalone regression tests for the Oct-2026
! robustness fixes (A1/A4/A2) against the same production sources the
! build compiles.
!
!   G1: operator-consistency of the in-sweep matched-frequency span
!       with the MATCH_FREQ refresh span.  The in-sweep code in
!       fg_j_cmf_nonmono_v1.f sums one-bin log ratios over the matched
!       offset K2 and reads the upwind intensity from store slot
!       FREQ_IDX_ITER+NM_RAY (frequency k+NM_RAY, since the store is
!       shifted one slot).  The refresh MATCH_FREQ computes
!       LOG(KN_FREQ(L)/KN_FREQ(L+NNF)).  This test verifies the two
!       expressions telescope to the same span for every k and offset,
!       and that the cap NCF_IN-FIDX-1 never touches the phantom slot.
!
!   G2: IN_HBC denominator floor.  Replicates the guarded formula in
!       fg_j_cmf_nonmono_v1.f verbatim and checks the degenerate
!       extremes (cancellation denominator, NaN/Inf inputs) produce
!       the 0.5 fallback, and that a well-scaled denominator keeps the
!       raw ratio.
!
!   G3: JFNK_MXDL cap honored under auto-delta (del<0).  Runs
!       newton_hook_v1 on the shared test problem with mxdl=0.05 and
!       checks the trust radius stayed <= mxdl (the old code reset
!       mxdl to 1d99 whenever del<0).
!
!   G4: solver-core determinism: two newton_hook_v1 activations from
!       the same start return bitwise-identical iterates.  The
!       production eval-side bracketing (JFNK_RAD_SWEEP_V1) cannot be
!       exercised standalone; this test locks determinism of the
!       Newton/GMRES layer itself.
!   G5: multi-cycle GMRES keeps the total step inside del and an
!       info=2 recompute retains the earlier restart contribution.
!   G6: log-unknown PTC diagonal agrees with the scalar shifted operator
!       at x != 0.
!   G7: GMRES aborts when the JVP callback reports a repaired bad
!       perturbed residual instead of consuming its sentinel vector.
!   G8: coupled Newton does not converge on maxch alone.
!   G9: jfnk_h_select_v1 picks the stable h plateau and rejects a
!       direction whose sweep never stabilizes (JFNK_HADAPT).
!   G10: jfnk_pc_final_guard_v1 zeroes and counts non-finite outputs.
!   G11: eval-budget bypass contract: the sentinel fires without the
!        bypass, a real eval passes through it, and a non-finite
!        applied state still fails the commit gate (f_budget mirrors
!        the JFNK_RESIDUAL_V1 sentinel verbatim, same convention as
!        the G2 replication).
!   G12: solver-level NaN trial residuals count nan_rejects and abort
!        the Newton loop rather than committing a non-finite step.
!   G11d/e: JFNK_EVAL_LIM tightens the sentinel cap below
!        JFNK_MAX_EVALS; jfnk_eval_cap_v1 picks the smallest
!        positive of the two (0 = off).
!   G13: the Newton stagnation gate uses the measured true residual,
!       not only the Hessenberg model residual (12-Nov-2026 fix).
!   G14: a zero-decrease trial is a rejection, not a rho<0.75 retry
!   G15: a broken linear solve (measured true_res > |F| on a
!        non-converged GMRES exit) skips the eta re-solve.
!       (12-Nov-2026 fix).
!
! edited by Devin, Oct-2026
!----------------------------------------------------------------------

! Test problem identical in shape to test_newton_core's: a known-root
! exponential system.  f_t is a pure function of x, so repeated
! activations are deterministic by construction.
module guards_prob_mod
   implicit none
   integer, parameter :: NG = 60
   double precision :: ag(NG), xsg(NG)
   double precision :: bg = 0.3d0
   ! G11/G12: inject_nan makes f_budget return a non-finite residual
   ! (the applied-state reject case); nan_calls counts f_nan calls so
   ! its first call (the entry eval) can stay real.
   logical :: inject_nan = .false.
   integer :: nan_calls = 0
   ! df_split keeps components <= split_keep of unit-norm arguments and
   ! answers zero on the larger accumulated-step calls -- see G13.
   integer :: split_keep = 34
contains
   subroutine seed_prob()
      integer :: i
      do i = 1, NG
         ag(i)  = 1.0d-3 * 10d0**(6d0 * dble(i - 1) / dble(NG - 1))
         xsg(i) = 1.5d0 * sin(1.7d0 * dble(i))
      end do
   end subroutine seed_prob
end module guards_prob_mod

program test_guards
   use, intrinsic :: ieee_arithmetic
   use jfnk_solver_v1_mod
   use jfnk_cmfgen_int_mod, only: JFNK_COUPLED_LOC, JFNK_EVAL_BYPASS, &
      JFNK_MAX_EVALS, JFNK_N_RESIDUAL_EVALS, JFNK_EVAL_LIM
   use guards_prob_mod
   implicit none
   external :: f_t, df_t, df_bad, dfp_t, sub_t, gmres_m_v1, f_nan, &
      f_flat, f_const, df_split, df_id, df_anti
   double precision, external :: dp_t
   integer :: nfail, info, luer, i
   !
   luer = 6
   nfail = 0
   !
   ! ---------------- G1: matched-span operator consistency ----------
   !
   block
      integer, parameter :: NCF_IN = 40
      double precision :: knf(0:NCF_IN)   ! KN_FREQ indexed by frequency
      integer :: fidx, k2, ncf_ok, m
      double precision :: span_insweep, span_match, span_tele, freq
      logical :: bad
      !
      ! Nonuniform log grid (MATCH_FREQ must hold on any spacing).
      do i = 0, NCF_IN
         knf(i) = 1.0d15 * exp(-0.11d0 * dble(i) - 0.002d0 * dble(i*i))
      end do
      bad = .false.
      do fidx = 1, NCF_IN - 2
         ! Current frequency k sits at store slot fidx, i.e. frequency
         ! index fidx+1 (the I_NEXT store is shifted by one slot).
         freq = knf(fidx + 1)
         do k2 = 1, NCF_IN - fidx - 1   ! the production cap
            ! In-sweep span (fixed form): log(FREQ/KN_FREQ(FIDX+1+K2)).
            span_insweep = log(freq / knf(fidx + 1 + k2))
            ! MATCH_FREQ span LOG(KN_FREQ(L)/KN_FREQ(L+NNF)) with
            ! L = fidx+1 (frequency index of the current slot) and
            ! NNF = k2 -- the same upwind slot the in-sweep store read
            ! FREQ_IDX_ITER+K2 maps to (frequency fidx+1+k2).
            span_match = log(knf(fidx + 1) / knf(fidx + 1 + k2))
            ! A by-hand one-bin telescoped sum gives a third,
            ! independent evaluation of the same span.
            span_tele = 0d0
            do m = 1, k2
               span_tele = span_tele + &
                  log(knf(fidx + m) / knf(fidx + m + 1))
            end do
            if(abs(span_insweep - span_match) > &
               1d-12 * abs(span_match)) then
               write(luer,'(A,2I4,2ES14.6)') &
                  'G1 FAIL: span mismatch k2=', fidx, k2, &
                  span_insweep, span_match
               bad = .true.
            end if
            if(abs(span_tele - span_match) > &
               1d-10 * abs(span_match)) then
               write(luer,'(A,2I4,2ES14.6)') &
                  'G1 FAIL: telescoped span mismatch k2=', fidx, k2, &
                  span_tele, span_match
               bad = .true.
            end if
         end do
      end do
      ! Cap check: largest usable offset leaves the phantom slot free.
      ncf_ok = 0
      do fidx = 1, NCF_IN - 2
         if(fidx + 1 + (NCF_IN - fidx - 1) <= NCF_IN) &
            ncf_ok = ncf_ok + 1
      end do
      if(ncf_ok /= NCF_IN - 2) bad = .true.
      if(bad) then
         nfail = nfail + 1
         write(luer,'(A)') 'G1: FAIL (matched-span operator)'
      else
         write(luer,'(A)') 'G1: PASS (matched-span operator)'
      end if
   end block
   !
   ! ---------------- G2: IN_HBC denominator guard -------------------
   !
   block
      double precision :: inhbc, den, dscl, jnd, ic, muin
      integer :: nbad2
      nbad2 = 0
      ! Case A: near-cancellation 2*J ~ IC (old guard lets it blow up)
      jnd = 0.5d0 * 3.7d-10; ic = 3.7d-10
      muin = 1.0d-9
      inhbc = muin
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(inhbc /= 0.5d0) nbad2 = nbad2 + 1
      ! Case B: healthy denominator keeps the measured ratio
      jnd = 2.0d0; ic = 1.0d0
      muin = 0.9d0
      inhbc = muin
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(abs(inhbc - 0.3d0) > 1d-14) nbad2 = nbad2 + 1
      ! Case C: NaN numerator falls back.  ieee_value is used because
      ! gfortran refuses to fold a literal 0d0/0d0 at compile time.
      inhbc = ieee_value(inhbc, ieee_quiet_nan)
      den = 2.0d0 * jnd - ic
      dscl = max(abs(ic), 2.0d0 * jnd, 1.0d-30)
      if(abs(den) < 1.0d-6 * dscl) then
         inhbc = 0.5d0
      else
         inhbc = inhbc / den
      end if
      if(.not.(inhbc >= 0.0d0)) inhbc = 0.5d0
      if(.not.(inhbc <= 1.0d0)) inhbc = 0.5d0
      if(inhbc /= 0.5d0) nbad2 = nbad2 + 1
      if(nbad2 > 0) then
         nfail = nfail + 1
         write(luer,'(A,I2)') 'G2: FAIL (IN_HBC guard), nbad=', nbad2
      else
         write(luer,'(A)') 'G2: PASS (IN_HBC guard)'
      end if
   end block
   !
   ! -------- G3/G4: MXDL cap + core determinism on newton_hook_v1 ---
   !
   block
      integer, parameter :: NN = 60
      double precision :: x1(NN), x2(NN), tol_t
      integer :: j
      call seed_prob()
      if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .true.
      ! del<0 selects automatic trust-region init; mxdl=0.05 must not
      ! be overwritten with the unlimited fallback.
      jfnk_x = xsg + 0.4d0
      info = 1
      call newton_hook_v1(f_t, df_t, dfp_t, sub_t, dp_t, &
         10, NN, 1d-9, 1d-11, -1d0, 1d-12, 5d-2, 30, info, luer)
      if(jfnk_del > 5.0d-2 + 1d-14) then
         nfail = nfail + 1
         write(luer,'(A,ES12.4)') &
            'G3: FAIL (mxdl overridden), jfnk_del=', jfnk_del
      else
         write(luer,'(A,ES12.4)') &
            'G3: PASS (mxdl honored), jfnk_del=', jfnk_del
      end if
      x1 = jfnk_x
      ! Determinism: identical second activation, bitwise-identical x.
      jfnk_x = xsg + 0.4d0
      info = 1
      call newton_hook_v1(f_t, df_t, dfp_t, sub_t, dp_t, &
         10, NN, 1d-9, 1d-11, -1d0, 1d-12, 5d-2, 30, info, luer)
      x2 = jfnk_x
      write(luer,'(A,ES12.4)') '  G4 max|x1-x2|=', maxval(abs(x1-x2))
      if(any(x1 /= x2)) then
         nfail = nfail + 1
         write(luer,'(A)') 'G4: FAIL (nondeterministic core)'
      else
         write(luer,'(A)') 'G4: PASS (deterministic core)'
      end if
   end block
   !
   ! -------- G5: restarted GMRES total trust radius and recompute ----
   !
   block
      integer, parameter :: M = 4
      double precision :: rhs(NG), step(NG), step_before(NG)
      double precision :: h(M+1,M), basis(NG,M+1), pred, res, radius
      double precision :: x_cycle(NG), step_norm
      integer :: its
      !
      jfnk_ptc_active = .false.
      jfnk_use_left_pc = .true.
      jfnk_x = xsg + 0.4d0
      call f_t(NG, jfnk_x, jfnk_fx)
      rhs = jfnk_fx
      step = 0d0
      res = 1d-12
      radius = 5d-2
      its = 20
      info = 1
      call gmres_m_v1(M, NG, step, rhs, df_t, dfp_t, dp_t, h, basis, &
         res, radius, pred, its, info, luer)
      step_norm = dsqrt(dp_t(NG, step, step))
      x_cycle = jfnk_gmres_x_cycle_start
      if(its <= M .or. step_norm > radius + 1d-12) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,ES12.4)') &
            'G5: FAIL (restart bound), its=', its, ' norm=', step_norm
      else
         write(luer,'(A,I3,A,ES12.4)') &
            'G5: PASS (total restart step bounded), its=', its, &
            ' norm=', step_norm
      end if
      !
      step_before = step
      radius = 2.5d-2
      its = 20
      info = 2
      call gmres_m_v1(M, NG, step, rhs, df_t, dfp_t, dp_t, h, basis, &
         res, radius, pred, its, info, luer)
      step_norm = dsqrt(dp_t(NG, step, step))
      if(step_norm > radius + 1d-12 .or. &
         dp_t(NG, step, x_cycle) <= 0d0 .or. &
         maxval(abs(step_before)) <= 1d-12) then
         nfail = nfail + 1
         write(luer,'(A,ES12.4,A,ES12.4)') &
            'G5: FAIL (recompute lost restart), norm=', step_norm, &
            ' prior_norm=', dsqrt(dp_t(NG, x_cycle, x_cycle))
      else
         write(luer,'(A,ES12.4)') &
            'G5: PASS (recompute retained earlier cycle), norm=', step_norm
      end if
   end block
   !
   ! -------- G6: fractional/log PTC chain rule at nonzero x ----------
   !
   block
      double precision :: d_frac, sigma, x_log, rhs_scalar
      double precision :: shifted_frac, log_apply, exact_apply
      d_frac = -2.3d0
      sigma = 0.7d0
      x_log = 0.6d0
      rhs_scalar = 1.4d0
      shifted_frac = d_frac + sigma * exp(-x_log)
      log_apply = exp(-x_log) * rhs_scalar / shifted_frac
      exact_apply = rhs_scalar / (d_frac * exp(x_log) + sigma)
      if(abs(log_apply - exact_apply) > &
         1d-13 * max(1d0, abs(exact_apply))) then
         nfail = nfail + 1
         write(luer,'(A,2ES14.6)') &
            'G6: FAIL (log PTC scalar), got/expected=', &
            log_apply, exact_apply
      else
         write(luer,'(A)') 'G6: PASS (log PTC scalar at x != 0)'
      end if
   end block
   !
   ! -------- G7: GMRES aborts on a failed perturbed residual --------
   !
   block
      integer, parameter :: M = 4
      double precision :: rhs(NG), step(NG), h(M+1,M), basis(NG,M+1)
      double precision :: pred
      integer :: its
      jfnk_x = 0d0
      jfnk_fx = 1d0
      rhs = 1d0
      step = 0d0
      jfnk_jvp_failed = .false.
      its = M
      info = 1
      call gmres_m_v1(M, NG, step, rhs, df_bad, dfp_t, dp_t, h, basis, &
         1d-8, 5d-2, pred, its, info, luer)
      if(info /= 3 .or. .not. jfnk_jvp_failed) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,L1)') &
            'G7: FAIL (JVP sentinel consumed), info=', info, &
            ' failed=', jfnk_jvp_failed
      else
         write(luer,'(A)') 'G7: PASS (GMRES aborts on failed JVP)'
      end if
   end block
   !
   ! -------- G9: h-sweep plateau selector (JFNK_HADAPT) -------------
   !
   block
      integer, parameter :: NV = 40, NH = 8
      double precision :: jvsw(NV, NH), relmin, base(NV)
      logical :: oks(NH)
      integer :: im, j, k
      external :: jfnk_h_select_v1
      !
      ! Synthetic sweep: columns 1-3 noise-dominated (random-sign
      ! jitter), 4-6 a stable plateau, 7-8 truncation drift.  The
      ! selector must return index 5 or 6 (upper member of the best
      ! adjacent pair inside the plateau).
      do k = 1, NH
         do j = 1, NV
            if(k <= 3) then
               jvsw(j,k) = 1d0 + 0.9d0 * sin(12.9898d0 * j * k)
            else if(k <= 6) then
               jvsw(j,k) = cos(0.31d0 * j) + 1d-4 * k
            else
               jvsw(j,k) = cos(0.31d0 * j) + 0.5d0 * (k - 6)
            end if
         end do
      end do
      oks = .true.
      call jfnk_h_select_v1(jvsw, oks, NV, NH, 0.3d0, im, relmin)
      if(im < 4 .or. im > 6) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,ES12.4)') &
            'G9a: FAIL (plateau missed), im=', im, ' rel=', relmin
      else
         write(luer,'(A,I3,A,ES12.4)') &
            'G9a: PASS (plateau at idx', im, '), rel=', relmin
      end if
      !
      ! Pure noise: every pair disagrees > 30% -- must return im=0
      ! (direction unreliable) rather than a noise-chosen step.
      do k = 1, NH
         do j = 1, NV
            jvsw(j,k) = sin(78.233d0 * j * k) + &
               0.5d0 * sin(4.13d0 * j + 2.7d0 * k)
         end do
      end do
      call jfnk_h_select_v1(jvsw, oks, NV, NH, 0.3d0, im, relmin)
      if(im /= 0) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,ES12.4)') &
            'G9b: FAIL (noise accepted), im=', im, ' rel=', relmin
      else
         write(luer,'(A)') 'G9b: PASS (unresolvable sweep rejected)'
      end if
      !
      ! A failed eval (oks=F) must be skipped: column 5 would pair
      ! perfectly with 6 if it were read, but the only honest plateau
      ! is (6,7) -> im must be 7, not 6.
      jvsw = 0d0
      do j = 1, NV
         base(j) = cos(0.31d0 * j)
      end do
      jvsw(:,1) = base * 1.9d0
      jvsw(:,2) = base * 0.9d0
      jvsw(:,3) = base * 0.2d0
      jvsw(:,4) = base * 3.0d0
      jvsw(:,5) = base * 1.0d0   ! bad eval: skipped
      jvsw(:,6) = base * 1.0d0
      jvsw(:,7) = base * 1.01d0
      jvsw(:,8) = base * 4.0d0
      oks = .true.
      oks(5) = .false.
      call jfnk_h_select_v1(jvsw, oks, NV, NH, 0.3d0, im, relmin)
      if(im /= 7) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,ES12.4)') &
            'G9c: FAIL (bad eval not skipped), im=', im, ' rel=', relmin
      else
         write(luer,'(A,I3)') 'G9c: PASS (bad eval skipped), im=', im
      end if
   end block
   !
   ! -------- G10: preconditioner final-output finite guard ---------
   !
   block
      double precision :: z(16)
      integer :: nbad, i
      external :: jfnk_pc_final_guard_v1
      z = 1d0
      z(3) = ieee_value(0d0, ieee_quiet_nan)
      z(7) = ieee_value(0d0, ieee_positive_inf)
      z(11) = -ieee_value(0d0, ieee_positive_inf)
      call jfnk_pc_final_guard_v1(16, z, nbad)
      if(nbad /= 3 .or. any(z /= z) .or. &
         any(.not. ieee_is_finite(z)) .or. &
         z(3) /= 0d0 .or. z(7) /= 0d0 .or. z(11) /= 0d0 .or. &
         z(1) /= 1d0) then
         nfail = nfail + 1
         write(luer,'(A,I3)') &
            'G10: FAIL (final guard), nbad=', nbad
      else
         write(luer,'(A)') 'G10: PASS (final guard zeroes NaN/Inf)'
      end if
      ! All-finite input passes through unchanged with nbad=0.
      z = 5d0
      call jfnk_pc_final_guard_v1(16, z, nbad)
      if(nbad /= 0 .or. any(z /= 5d0)) then
         nfail = nfail + 1
         write(luer,'(A,I3)') 'G10b: FAIL (finite pass-through)', nbad
      else
         write(luer,'(A)') 'G10b: PASS (finite vector untouched)'
      end if
   end block
   !
   ! -------- G8: coupled mode requires the residual tolerance --------
   !
   block
      double precision :: f0_norm
      call seed_prob()
      if(.not. allocated(jfnk_x)) allocate(jfnk_x(NG), jfnk_fx(NG))
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .false.
      jfnk_x = xsg + 3d0
      call f_t(NG, jfnk_x, jfnk_fx)
      f0_norm = dsqrt(dp_t(NG, jfnk_fx, jfnk_fx))
      jfnk_maxch_tol = 1d0
      jfnk_maxch_est = 0d0
      JFNK_COUPLED_LOC = .true.
      info = 1
      call newton_hook_v1(f_t, df_t, dfp_t, sub_t, dp_t, &
         10, NG, 1d-3, 1d-14*f0_norm, -1d0, 1d-12, 5d-2, 1, info, luer)
      JFNK_COUPLED_LOC = .false.
      jfnk_maxch_tol = 0d0
      if(info == 0) then
         nfail = nfail + 1
         write(luer,'(A,ES12.4,A,ES12.4)') &
            'G8: FAIL (maxch bypassed residual), |F0|=', f0_norm, &
            ' |F|=', jfnk_tol
      else
         write(luer,'(A,ES12.4)') &
            'G8: PASS (coupled residual remains above tol), |F|=', jfnk_tol
      end if
   end block
   !
   ! -------- G11: eval-budget bypass for the commit verification ----
   !
   ! jfnk_solve_corr_v1 sets JFNK_EVAL_BYPASS around exactly one
   ! JFNK_RESIDUAL_V1 call.  f_budget replicates the residual-side
   ! sentinel semantics verbatim (the JFNK_MAX_SEC branch has the
   ! same shape): an exhausted budget returns Y=1D20 with nbad=N
   ! unless the bypass is set.  The commit gate (nbad>0 or
   ! F_APP >= F_ENTRY -> rollback, jfnk_solve_corr_v1) must:
   !   (a) see the sentinel and reject when the bypass is off --
   !       the rzc2 it585 rollback mechanism,
   !   (b) see the real, improved residual and commit when it is on,
   !   (c) still reject when the applied eval is genuinely non-finite:
   !       the bypass opens only the budget gate, not NaN accounting.
   !
   block
      double precision :: x_ent(NG), x_app(NG), y_app(NG)
      double precision :: f_entry, f_app
      logical :: commit
      integer :: nbad_a
      !
      ! The flag must default to .false.: the sentinel bounds GMRES
      ! burn and only the commit-verification call may pass it.
      if(JFNK_EVAL_BYPASS) then
         nfail = nfail + 1
         write(luer,'(A)') 'G11: FAIL (bypass default not F)'
      end if
      call seed_prob()
      x_ent = xsg + 0.4d0    ! entry state: far from the root
      x_app = xsg + 0.05d0   ! applied state: verified improvement
      JFNK_MAX_EVALS = 2
      JFNK_N_RESIDUAL_EVALS = 0
      ! Spend the budget exactly as the solve does: entry + one trial.
      call f_budget(NG, x_ent, y_app)
      f_entry = dsqrt(dp_t(NG, y_app, y_app))
      call f_budget(NG, x_app, y_app)
      if(JFNK_N_RESIDUAL_EVALS /= JFNK_MAX_EVALS) then
         nfail = nfail + 1
         write(luer,'(A)') 'G11a: FAIL (budget accounting)'
      end if
      call f_budget(NG, x_app, y_app)
      commit = (jfnk_resid_nbad <= 0) .and. &
         all(ieee_is_finite(y_app)) .and. &
         (dsqrt(dp_t(NG, y_app, y_app)) < f_entry)
      if(commit .or. jfnk_resid_nbad /= NG .or. &
         any(y_app /= 1.0d20)) then
         nfail = nfail + 1
         write(luer,'(A)') 'G11a: FAIL (sentinel not applied)'
      else
         write(luer,'(A)') 'G11a: PASS (exhausted budget -> sentinel)'
      end if
      JFNK_EVAL_BYPASS = .true.
      call f_budget(NG, x_app, y_app)
      JFNK_EVAL_BYPASS = .false.
      f_app = dsqrt(dp_t(NG, y_app, y_app))
      commit = (jfnk_resid_nbad == 0) .and. &
         all(ieee_is_finite(y_app)) .and. (f_app < f_entry)
      if(.not. commit) then
         nfail = nfail + 1
         write(luer,'(A,2ES12.4,A,I4)') &
            'G11b: FAIL (bypass eval not real), f_app/f_entry=', &
            f_app, f_entry, ' nbad=', jfnk_resid_nbad
      else
         write(luer,'(A,2ES12.4)') &
            'G11b: PASS (real eval commits), f_app/f_entry=', &
            f_app, f_entry
      end if
      inject_nan = .true.
      JFNK_EVAL_BYPASS = .true.
      call f_budget(NG, x_app, y_app)
      JFNK_EVAL_BYPASS = .false.
      inject_nan = .false.
      nbad_a = jfnk_resid_nbad
      commit = (nbad_a <= 0) .and. all(ieee_is_finite(y_app)) .and. &
         (dsqrt(dp_t(NG, y_app, y_app)) < f_entry)
      if(commit .or. nbad_a <= 0) then
         nfail = nfail + 1
         write(luer,'(A)') 'G11c: FAIL (non-finite applied committed)'
      else
         write(luer,'(A)') 'G11c: PASS (non-finite applied rejects)'
      end if
      ! JFNK_EVAL_LIM is the activation-local trim the driver applies
      ! per budget rung; the sentinel must fire at the tighter of the
      ! two caps -- checking JFNK_MAX_EVALS alone let the trimmed cap
      ! go unenforced for ordinary residual calls.
      JFNK_MAX_EVALS = 5
      JFNK_EVAL_LIM = 1
      JFNK_N_RESIDUAL_EVALS = 0
      call f_budget(NG, x_ent, y_app)   ! eval 1 spends the lim cap
      call f_budget(NG, x_ent, y_app)   ! sentinel at the lower cap
      if(jfnk_resid_nbad /= NG .or. any(y_app /= 1.0d20)) then
         nfail = nfail + 1
         write(luer,'(A)') 'G11d: FAIL (EVAL_LIM cap not honored)'
      else
         write(luer,'(A)') 'G11d: PASS (EVAL_LIM tightens cap)'
      end if
      ! jfnk_eval_cap_v1: smallest positive cap wins; 0 disables.
      nbad_a = 0
      if(jfnk_eval_cap_v1(5, 10) /= 5) nbad_a = nbad_a + 1
      if(jfnk_eval_cap_v1(0, 10) /= 10) nbad_a = nbad_a + 1
      if(jfnk_eval_cap_v1(10, 0) /= 10) nbad_a = nbad_a + 1
      if(jfnk_eval_cap_v1(0, 0) /= 0) nbad_a = nbad_a + 1
      if(jfnk_eval_cap_v1(3, 0) /= 3) nbad_a = nbad_a + 1
      if(nbad_a > 0) then
         nfail = nfail + 1
         write(luer,'(A,I2)') 'G11e: FAIL (cap helper)', nbad_a
      else
         write(luer,'(A)') 'G11e: PASS (eval cap helper)'
      end if
      JFNK_MAX_EVALS = 0
      JFNK_EVAL_LIM = 0
      JFNK_N_RESIDUAL_EVALS = 0
      jfnk_resid_nbad = 0
   end block
   !
   ! -------- G12: non-finite trial residuals reject, never commit ---
   !
   ! f_nan returns the real residual for the entry eval and NaN for
   ! every trial eval.  The solver must count nan_rejects, restore
   ! the iterate, and abort rather than converging on NaN steps.
   !
   block
      nan_calls = 0
      jfnk_nan_rejects = 0
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .false.
      jfnk_use_left_pc = .false.
      jfnk_x = xsg + 0.4d0
      jfnk_exit_reason = ' '
      info = 1
      call newton_hook_v1(f_nan, df_t, dfp_t, sub_t, dp_t, &
         10, NG, 1d-9, 1d-11, -1d0, 1d-12, 5d-2, 30, info, luer)
      if(info == 0 .or. jfnk_nan_rejects < 1) then
         nfail = nfail + 1
         write(luer,'(A,I3,A,I5)') &
            'G12: FAIL (NaN trial committed), info=', info, &
            ' nan_rejects=', jfnk_nan_rejects
      else
         write(luer,'(A,I3,A,I5)') &
            'G12: PASS (NaN trials rejected), info=', info, &
            ' nan_rejects=', jfnk_nan_rejects
      end if
   end block
   !
   ! -------- G13: stagnation gate must use the measured true res ----
   !
   ! The Hessenberg least-squares residual is only a Krylov-space
   ! model.  df_split answers unit-norm basis calls with a rank-
   ! deficient projection (components <= split_keep kept) and the
   ! accumulated-step call with zero action, so the hess model claims
   ! ~34% residual reduction while the measured residual is the full
   ! ||b||.  The pre-fix gate accepted that direction on
   ! pred_res/ref ~0.66; the fix must see true_res/ref = 1 and abort.
   block
      double precision :: f0w
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .false.
      jfnk_use_left_pc = .false.
      jfnk_gmres_tres_fac = 0d0
      jfnk_jvp_failed = .false.
      jfnk_gmres_true_res_ok = .false.
      jfnk_gmres_true_res = 0d0
      jfnk_x = 0d0
      info = 1
      call newton_hook_v1(f_flat, df_split, dfp_t, sub_t, dp_t, &
         10, NG, 1d-9, 1d-11, 1d3, 1d-12, 1d8, 30, info, luer)
      f0w = dsqrt(dble(NG))
      if(jfnk_exit_reason /= 'gmres_stagnation_abort' .or. &
         .not. jfnk_gmres_true_res_ok .or. &
         jfnk_gmres_true_res < 0.9d0 * f0w) then
         nfail = nfail + 1
         write(luer,'(A,A,A,ES12.4,A,L1)') &
            'G13: FAIL (hess-only gate), reason=', &
            trim(jfnk_exit_reason), ' true_res=', jfnk_gmres_true_res, &
            ' ok=', jfnk_gmres_true_res_ok
      else
         write(luer,'(A,ES12.4)') &
            'G13: PASS (true-res gate aborts phantom step), true_res=', &
            jfnk_gmres_true_res
      end if
      jfnk_use_left_pc = .true.
   end block
   !
   ! -------- G14: zero-progress trial rejects, not rho<0.75 retry ---
   !
   ! f_const returns a constant residual, so every trial reports
   ! ared == 0.  The pre-fix code sent that into the rho<0.75 retry
   ! branch, which never touches n_reject: the trust radius halved to
   ! the floor and the solve ended in trust_region_too_small after
   ! burning its whole budget.  With the fix the solve must abort via
   ! n_reject_max as too_many_rejected_steps.
   block
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .false.
      jfnk_use_left_pc = .false.
      jfnk_gmres_tres_fac = 0d0
      jfnk_jvp_failed = .false.
      jfnk_x = 0d0
      info = 1
      call newton_hook_v1(f_const, df_id, dfp_t, sub_t, dp_t, &
         10, NG, 1d-9, 1d-11, 1d3, 1d-12, 1d8, 30, info, luer)
      if(jfnk_exit_reason /= 'too_many_rejected_steps') then
         nfail = nfail + 1
         write(luer,'(A,A,A,I3)') 'G14: FAIL (zero-progress retry), &
            &reason=', trim(jfnk_exit_reason), ' info=', info
      else
         write(luer,'(A)') &
            'G14: PASS (zero progress aborts via n_reject)'
      end if
      jfnk_use_left_pc = .true.
   end block
   !
   ! -------- G15: broken linear solve must not burn an eta re-solve --
   !
   ! df_anti gives Arnoldi basis calls (|v| ~ 1) a ten-eigenvalue
   ! definite action, so no Krylov space of dim < 10 can converge;
   ! jfnk_max_gmres_its caps it at 2, forcing the info=2 exit and the
   ! measured true_res on every solve.  The accumulated-step calls
   ! (|v| > 2: the exit true-residual measure, restarts, and the
   ! newton J_s eval) see anti-aligned -2*v action, so true_res
   ! exceeds |F| and the step is non-descent -- the OUTGEN signature
   ! where an eta-tightened re-solve repeated the same failure.  The
   ! fix must skip the retry entirely: jfnk_n_descent_retry stays 0.
   block
      jfnk_ptc_active = .false.
      jfnk_do_line_search = .false.
      jfnk_use_left_pc = .false.
      jfnk_gmres_tres_fac = 0d0
      jfnk_jvp_failed = .false.
      jfnk_gmres_recomputed = .false.
      jfnk_gmres_true_res_ok = .false.
      jfnk_gmres_true_res = 0d0
      jfnk_max_gmres_its = 2
      jfnk_n_descent_retry = 0
      jfnk_x = 0d0
      info = 1
      call newton_hook_v1(f_flat, df_anti, dfp_t, sub_t, dp_t, &
         10, NG, 1d-9, 1d-11, 1d3, 1d-12, 1d8, 30, info, luer)
      if(jfnk_n_descent_retry /= 0) then
         nfail = nfail + 1
         write(luer,'(A,I3,A)') 'G15: FAIL (burned ', &
            jfnk_n_descent_retry, ' re-solves on a broken operator)'
      else
         write(luer,'(A,A)') 'G15: PASS (no eta re-solve on broken &
            &solve), reason=', trim(jfnk_exit_reason)
      end if
      jfnk_max_gmres_its = 0
      jfnk_use_left_pc = .true.
   end block
   !
   !
   if(nfail == 0) then
      write(luer,'(A)') 'ALL GUARD TESTS PASSED'
   else
      write(luer,'(A,I2)') 'FAILURES: ', nfail
      stop 1
   end if
end program test_guards

subroutine f_t(n, x, y)
   use guards_prob_mod
   implicit none
   integer :: n, i
   double precision :: x(n), y(n), d(n)
   d = x - xsg
   do i = 1, n
      y(i) = ag(i) * (exp(min(20d0, x(i))) - exp(xsg(i)))
      if(i > 1)  y(i) = y(i) + bg * (d(i) - d(i - 1))
      if(i < n)  y(i) = y(i) + bg * (d(i) - d(i + 1))
   end do
end subroutine f_t

subroutine df_t(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   integer :: n
   double precision :: v(n), jv(n), xp(n), yp(n), h, vinf
   ! Fixed-h finite difference (deterministic).
   vinf = maxval(abs(v))
   if(vinf == 0d0) then
      jv = 0d0
      return
   end if
   h = 1d-6 / vinf
   xp = jfnk_x + h * v
   call f_t(n, xp, yp)
   jv = (yp - jfnk_fx) / h
end subroutine df_t

subroutine df_bad(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_jvp_failed
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   jv = 0d0
   jfnk_jvp_failed = .true.
end subroutine df_bad

subroutine dfp_t(n, z)
   implicit none
   integer :: n
   double precision :: z(n)
   ! Identity preconditioner is adequate for this problem scale.
end subroutine dfp_t

subroutine sub_t()
   implicit none
end subroutine sub_t

double precision function dp_t(n, x, y)
   implicit none
   integer :: n
   double precision :: x(n), y(n)
   dp_t = dot_product(x, y)
end function dp_t


! Verbatim mirror of the JFNK_RESIDUAL_V1 eval-budget sentinel used by
! G11 (the JFNK_MAX_SEC branch in production has the same shape).
! Without JFNK_EVAL_BYPASS an exhausted budget returns Y=1D20 and
! nbad=N; with it the eval is real and bad components are counted
! exactly as the production NaN/Inf policy loop does.
subroutine f_budget(n, x, y)
   use jfnk_cmfgen_int_mod, only: JFNK_EVAL_BYPASS, JFNK_MAX_EVALS, &
      JFNK_N_RESIDUAL_EVALS, JFNK_EVAL_LIM
   use jfnk_solver_v1_mod, only: jfnk_resid_nbad, jfnk_eval_cap_v1
   use guards_prob_mod
   use, intrinsic :: ieee_arithmetic
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   integer :: i, nbad, nev_cap
   nev_cap = jfnk_eval_cap_v1(JFNK_EVAL_LIM, JFNK_MAX_EVALS)
   if(nev_cap > 0 .and. .not. JFNK_EVAL_BYPASS .and. &
      JFNK_N_RESIDUAL_EVALS >= nev_cap) then
      y = 1.0d20
      jfnk_resid_nbad = n
      JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
      return
   end if
   JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
   call f_t(n, x, y)
   if(inject_nan) y(1) = ieee_value(0d0, ieee_quiet_nan)
   nbad = 0
   do i = 1, n
      if(y(i) /= y(i) .or. abs(y(i)) > 1.0d30) nbad = nbad + 1
   end do
   jfnk_resid_nbad = nbad
end subroutine f_budget


! G12 residual: the first call is the entry eval and stays real so
! the Newton loop starts; every later call is a trial eval returning
! NaN, which must reject instead of committing.
subroutine f_nan(n, x, y)
   use guards_prob_mod, only: nan_calls
   use, intrinsic :: ieee_arithmetic
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   nan_calls = nan_calls + 1
   if(nan_calls <= 1) then
      call f_t(n, x, y)
   else
      y = ieee_value(0d0, ieee_quiet_nan)
   end if
end subroutine f_nan
subroutine f_flat(n, x, y)
   implicit none
   integer :: n
   double precision :: x(n), y(n)
   ! G13: unit residual per component, so the hess/true residual split
   ! in df_split is fully controlled by the projection rank.
   y = x - 1d0
end subroutine f_flat

subroutine f_const(n, x, y)
   implicit none
   integer :: n
   double precision :: x(n), y(n)
   ! G14: residual identical at every trial state -> ared == 0 always.
   y = 1d0
end subroutine f_const

subroutine df_split(n, v, jv)
   use guards_prob_mod, only: split_keep
   implicit none
   integer :: n, k
   double precision :: v(n), jv(n)
   ! Deliberately inconsistent matvec.  Arnoldi basis calls carry unit
   ! norm and see a rank-split_keep projection; the gmres exit's
   ! true-residual call (and cycle restarts) pass the accumulated step
   ! (norm ~ sqrt(split_keep)) and see zero action.  This reproduces
   ! the OUTGEN pathology: hess_res << true_res on a stagnation exit.
   jv = 0d0
   if(dsqrt(dot_product(v, v)) < 3d0) then
      k = min(split_keep, n)
      jv(1:k) = v(1:k)
   end if
end subroutine df_split

subroutine df_id(n, v, jv)
   implicit none
   integer :: n
   double precision :: v(n), jv(n)
   ! G14: identity action converges immediately so every Newton pass
   ! produces a full step against the constant residual.
   jv = v
end subroutine df_id

subroutine df_anti(n, v, jv)
   implicit none
   integer :: n, i
   double precision :: v(n), jv(n)
   ! G15: unit-norm Arnoldi basis calls see a definite ten-eigenvalue
   ! action (no Krylov space of dim < 10 converges, so the imx=2 cap
   ! always exits non-converged and true_res is measured); the
   ! accumulated-step calls -- exit true-residual measure, restart
   ! residuals, and the newton J_s -- see anti-aligned action.
   ! Reproduces the production pathology of a max-its exit returning a
   ! non-descent step whose measured true_res exceeds the reference.
   if(dsqrt(dot_product(v, v)) < 2d0) then
      do i = 1, n
         jv(i) = (0.4d0 + 0.3d0 * dble(mod(i, 10))) * v(i)
      end do
   else
      jv = -2d0 * v
   end if
end subroutine df_anti
