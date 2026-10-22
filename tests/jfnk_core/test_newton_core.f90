!----------------------------------------------------------------------
! TEST_NEWTON_CORE: standalone unit test of newton_hook_v1 (the
! production Newton/GMRES/hookstep core in
! new_main/mod_subs/jfnk_solver_v1.f90, compiled unchanged).
!
! Synthetic problem, N = 60:
!   F_i(x) = a_i*(exp(x_i) - exp(xs_i))
!            + b*(d_i - d_{i-1}) + b*(d_i - d_{i+1}),  d = x - xs
! (Dirichlet-style ends: the missing neighbour term is dropped).
! The coupling acts on x - xs so that xs is the exact known root;
! applied to x directly it would leave b*L*xs as a residual offset.
! a_i is log-spaced 1e-3 .. 1e-3*1e6 (badly scaled, deterministic),
! b = 0.3, xs_i = 1.5*sin(1.7*i), so xs is a known root.
!
! The residual callback adds the pseudo-transient term exactly as
! JFNK_RESIDUAL_V1 does -- (x - jfnk_ptc_x0)*jfnk_ptc_scale/dtau --
! so the tests exercise the SER anchor fix (06-Oct-2026): with the
! old baseline-anchored term the Newton norms measured |G| and the
! returned iterate could satisfy |G| < tol while F /= 0.
!
! Tests:
!   T1: PTC off -- plain Newton from x0 = xs + 0.5.
!   T2: PTC on  -- same start, verify the returned iterate has
!       |F_true| < tol (fails under the old x*s/dtau form).
!   T3: PTC on, far start xs + 3.0 -- exp nonlinearity overshoots the
!       full Newton step, exercising trust region / line search /
!       Cauchy fallback.
!
! edited by Devin, 06-Oct-2026
!----------------------------------------------------------------------
module test_prob_mod
   implicit none
   integer, parameter :: NN = 60
   double precision :: a(NN), xs(NN), xbase(NN), diagJ(NN)
   double precision :: b
   ! T4: fixed linear operator for the gmres_m_v1 tres tests, plus an
   ! optional constant offset that makes the matvec affine (v -> Av+c0)
   ! -- the Arnoldi model then converges while the explicitly measured
   ! ||b - A x|| keeps an O(c0) floor.
   double precision :: amat(NN,NN), c0v(NN)
   logical :: matvec_bias = .false.
   ! T4f: inject a NaN into the matvec output on a chosen call index so
   ! the tres measurement lands on a non-finite residual.  -1 disables.
   integer :: matvec_ct = 0, matvec_nan_at = -1
end module test_prob_mod


program test_newton_core
   use jfnk_solver_v1_mod
   use test_prob_mod
   implicit none
   external :: f_test, df_test, dfp_test, sub_test
   external :: df_lin, dfp_id, dfp_jac
   double precision, external :: dp_test
   integer :: i, info, nfail, luer
   integer :: its_g
   integer :: nevals_b, nevals_c, ntres_call
   double precision :: x0(NN), f0(NN), ftrue(NN), xsol(NN)
   double precision :: x_g(NN), b_g(NN), res_g, pred_g
   double precision :: h_g(11,10), v_g(NN,11)
   double precision :: fnrm0, tol_t, ftru_nrm, maxerr
   !
   luer = 6
   nfail = 0
   b = 0.3d0
   do i = 1, NN
      a(i)  = 1.0d-3 * 10d0**(6d0 * dble(i - 1) / dble(NN - 1))
      xs(i) = 1.5d0 * sin(1.7d0 * dble(i))
   end do
   !
   if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
   !
   ! ---------------- T1: PTC off ----------------
   !
   jfnk_ptc_active = .false.
   xbase = 0d0
   x0 = xs + 0.5d0
   jfnk_x = x0
   call f_true(NN, jfnk_x, f0)
   fnrm0 = dsqrt(dp_test(NN, f0, f0))
   tol_t = 1.0d-10 * fnrm0
   diagJ = a * exp(x0) + 2d0 * b          ! lagged BA-style diagonal
   info = 1
   call newton_hook_v1(f_test, df_test, dfp_test, sub_test, dp_test, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       50, info, luer)
   xsol = jfnk_x
   maxerr = maxval(abs(xsol - xs))
   write(luer,'(A,I3,A,ES12.4,A,ES12.4)') &
      'T1: info=', info, '  max|x-xs|=', maxerr, '  |F0|=', fnrm0
   if(info == 0 .and. maxerr < 1.0d-6) then
      write(luer,'(A)') 'T1: PASS'
   else
      write(luer,'(A)') 'T1: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---------------- T2: PTC on, unknown relative to baseline --------
   !
   ! Production x is ln(n/n0) -- relative to the solve baseline.  Shift
   ! the problem the same way: unknown u with absolute point xbase + u,
   ! xbase = xs + 0.5, and solve starts at u = 0.
   !
   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 1d0
   xbase = xs + 0.5d0
   jfnk_x = 0d0
   call f_true(NN, xbase, f0)
   call set_ptc_like_production(NN, f0)   ! scale + anchor = 0
   fnrm0 = dsqrt(dp_test(NN, f0, f0))
   tol_t = 1.0d-10 * fnrm0
   diagJ = a * exp(xbase) + 2d0 * b
   info = 1
   call newton_hook_v1(f_test, df_test, dfp_test, sub_test, dp_test, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       50, info, luer)
   xsol = xbase + jfnk_x
   ! The assertion that distinguishes the SER-anchored form: evaluate
   ! the TRUE F (no pseudo-time term) at the returned point.
   call f_true(NN, xsol, ftrue)
   ftru_nrm = dsqrt(dp_test(NN, ftrue, ftrue))
   write(luer,'(A,I3,A,ES12.4,A,ES12.4,A,ES12.4)') &
      'T2: info=', info, '  |F_true|=', ftru_nrm, &
      '  |F0|=', fnrm0, '  tol=', tol_t
   if(info == 0 .and. ftru_nrm < tol_t) then
      write(luer,'(A)') 'T2: PASS'
   else
      write(luer,'(A)') 'T2: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---------------- T3: PTC on, far start ---------------------------
   !
   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 1d0
   jfnk_n_ls_used = 0
   jfnk_n_tr_rejected = 0
   jfnk_n_cauchy_used = 0
   xbase = 0d0
   x0 = xs + 3.0d0
   jfnk_x = x0
   call f_true(NN, x0, f0)
   call set_ptc_like_production(NN, f0)
   fnrm0 = dsqrt(dp_test(NN, f0, f0))
   tol_t = 1.0d-10 * fnrm0
   diagJ = a * exp(x0) + 2d0 * b
   info = 1
   call newton_hook_v1(f_test, df_test, dfp_test, sub_test, dp_test, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       50, info, luer)
   xsol = jfnk_x
   maxerr = maxval(abs(xsol - xs))
   write(luer,'(A,I3,A,ES12.4,A,I4,A,I4,A,I4)') &
      'T3: info=', info, '  max|x-xs|=', maxerr, &
      '  ls_used=', jfnk_n_ls_used, '  tr_rejected=', &
      jfnk_n_tr_rejected, '  cauchy_used=', jfnk_n_cauchy_used
   if(info == 0 .and. maxerr < 1.0d-5) then
      write(luer,'(A)') 'T3: PASS'
   else
      write(luer,'(A)') 'T3: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---------------- T4: JFNK_TRES_FAC acceptance ---------------------
   !
   ! gmres_m_v1 driven directly on A x = b, identity preconditioner,
   ! no trust region.  T4a uses the consistent matvec jv = A v: the
   ! hessenberg claim is accepted once true_res <= fac*tol.
   !
   amat = 0d0
   do i = 1, NN
      amat(i,i) = 2d0
      if(i > 1)  amat(i,i-1) = -0.5d0
      if(i < NN) amat(i,i+1) = -0.5d0
   end do
   b_g = matmul(amat, xs)
   c0v = 1.0d-6
   jfnk_use_left_pc = .false.
   jfnk_jvp_failed = .false.
   jfnk_pc_clip_event = .false.
   jfnk_ptc_active = .false.
   jfnk_gmres_max_restart = 1
   jfnk_gmres_stag_tol = 1.0d-3
   jfnk_total_jvp_evals = 0
   !
   matvec_bias = .false.
   matvec_ct = 0
   jfnk_gmres_tres_fac = 1.5d0
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-10
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   maxerr = maxval(abs(x_g - xs))
   ! The tres measurement is the last matvec call of an accepted solve;
   ! T4f reuses this call index to land a NaN on it.
   ntres_call = matvec_ct
   write(luer,'(A,I3,A,ES12.4,A)') &
      'T4a: info=', info, '  max|x-xs|=', maxerr, &
      '  exit='//trim(jfnk_exit_reason)
   if(info == 0 .and. maxerr < 1.0d-6 .and. &
      jfnk_exit_reason == 'gmres_converged') then
      write(luer,'(A)') 'T4a: PASS'
   else
      write(luer,'(A)') 'T4a: FAIL'
      nfail = nfail + 1
   end if
   !
   ! T4b: same system but the matvec adds a fixed offset.  The LS
   ! residual converges while ||b - A x - c0|| keeps the offset floor,
   ! so every claimed convergence must fail the tres gate; with
   ! max_restart=1 the first claim restarts on the measured residual
   ! and the second claim exits 'gmres_tres_reject' at the cap.
   !
   matvec_bias = .true.
   jfnk_exit_reason = ' '
   jfnk_total_jvp_evals = 0
   x_g = 0d0
   res_g = 1.0d-10
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   nevals_b = jfnk_total_jvp_evals
   write(luer,'(A,I3,A,I4,A,A)') &
      'T4b: info=', info, '  evals=', nevals_b, &
      '  exit='//trim(jfnk_exit_reason)
   if(info /= 0 .and. jfnk_exit_reason == 'gmres_tres_reject') then
      write(luer,'(A)') 'T4b: PASS'
   else
      write(luer,'(A)') 'T4b: FAIL'
      nfail = nfail + 1
   end if
   !
   ! T4c: same biased claim with no restart headroom -- the tres gate
   ! must reject on the FIRST claim, spending strictly fewer evals than
   ! T4b (which burned a second cycle before hitting the cap).  tol is
   ! loosened to 1e-6 so the claim lands inside cycle 1; the c0 floor
   ! in the true residual (~5e-5 there) stays far above 1.5*tol.
   !
   jfnk_gmres_max_restart = 0
   jfnk_exit_reason = ' '
   jfnk_total_jvp_evals = 0
   x_g = 0d0
   res_g = 1.0d-6
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   nevals_c = jfnk_total_jvp_evals
   write(luer,'(A,I3,A,I4,A,A)') &
      'T4c: info=', info, '  evals=', nevals_c, &
      '  exit='//trim(jfnk_exit_reason)
   if(info /= 0 .and. jfnk_exit_reason == 'gmres_tres_reject' .and. &
      nevals_c < nevals_b) then
      write(luer,'(A)') 'T4c: PASS'
   else
      write(luer,'(A)') 'T4c: FAIL'
      nfail = nfail + 1
   end if
   jfnk_gmres_max_restart = 1
   !
   ! T4d/T4e: left preconditioning -- the production default.  The tres
   ! check then applies psolve to the measured residual, so these run
   ! with a nontrivial (non-cI) diagonal M to exercise the M^{-1} path.
   !
   jfnk_use_left_pc = .true.
   !
   matvec_bias = .false.
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-10
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_jac, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   maxerr = maxval(abs(x_g - xs))
   write(luer,'(A,I3,A,ES12.4,A)') &
      'T4d: info=', info, '  max|x-xs|=', maxerr, &
      '  exit='//trim(jfnk_exit_reason)
   if(info == 0 .and. maxerr < 1.0d-6 .and. &
      jfnk_exit_reason == 'gmres_converged') then
      write(luer,'(A)') 'T4d: PASS'
   else
      write(luer,'(A)') 'T4d: FAIL'
      nfail = nfail + 1
   end if
   !
   matvec_bias = .true.
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-10
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_jac, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   write(luer,'(A,I3,A,A)') &
      'T4e: info=', info, '  exit='//trim(jfnk_exit_reason)
   if(info /= 0 .and. jfnk_exit_reason == 'gmres_tres_reject') then
      write(luer,'(A)') 'T4e: PASS'
   else
      write(luer,'(A)') 'T4e: FAIL'
      nfail = nfail + 1
   end if
   jfnk_use_left_pc = .false.
   matvec_bias = .false.
   !
   ! T4f: non-finite true residual.  A calibration pass (tol loosened
   ! to 1e-6 so the claim lands inside cycle 1) learns which matvec
   ! call the tres measurement uses; the second pass injects a NaN on
   ! exactly that call.  true_res=NaN cannot verify the claim --
   ! without the fail-closed guard, NaN > acc_tol is .false. and this
   ! would return 'gmres_converged' on an unverifiable residual.
   !
   matvec_bias = .false.
   jfnk_gmres_max_restart = 0
   matvec_ct = 0
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-6
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   ntres_call = matvec_ct
   matvec_nan_at = ntres_call
   matvec_ct = 0
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-6
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   write(luer,'(A,I3,A,A)') &
      'T4f: info=', info, '  exit='//trim(jfnk_exit_reason)
   if(info /= 0 .and. jfnk_exit_reason == 'gmres_tres_nonfinite') then
      write(luer,'(A)') 'T4f: PASS'
   else
      write(luer,'(A)') 'T4f: FAIL'
      nfail = nfail + 1
   end if
   !
   ! T4g: same non-finite claim but WITH restart headroom -- the
   ! nonfinite exit must still fire immediately, with no further
   ! matvec (a restart would run calls beyond the injected index).
   ! matvec_ct == matvec_nan_at proves nothing ran after the tres
   ! measurement.
   !
   matvec_ct = 0
   jfnk_gmres_max_restart = 1
   jfnk_exit_reason = ' '
   x_g = 0d0
   res_g = 1.0d-6
   its_g = 60
   info = 1
   call gmres_m_v1(10, NN, x_g, b_g, df_lin, dfp_id, dp_test, &
                   h_g, v_g, res_g, 0d0, pred_g, its_g, info, luer)
   write(luer,'(A,I3,A,I4,A,A)') &
      'T4g: info=', info, '  ct=', matvec_ct, &
      '  exit='//trim(jfnk_exit_reason)
   if(info /= 0 .and. jfnk_exit_reason == 'gmres_tres_nonfinite' &
      .and. matvec_ct == matvec_nan_at) then
      write(luer,'(A)') 'T4g: PASS'
   else
      write(luer,'(A)') 'T4g: FAIL'
      nfail = nfail + 1
   end if
   jfnk_gmres_max_restart = 1
   matvec_nan_at = -1
   jfnk_gmres_tres_fac = 0d0
   !
   if(nfail == 0) then
      write(luer,'(A)') 'ALL TESTS PASSED'
   else
      write(luer,'(A,I2)') 'FAILURES: ', nfail
      call exit(1)
   end if

contains

   ! F without the pseudo-time term, at absolute point x.
   subroutine f_true(n, x, y)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: y(n)
      integer :: i
      double precision :: d(n)
      d = x - xs
      do i = 1, n
         y(i) = a(i) * (exp(x(i)) - exp(xs(i)))
         if(i > 1) y(i) = y(i) + b * (d(i) - d(i-1))
         if(i < n) y(i) = y(i) + b * (d(i) - d(i+1))
      end do
   end subroutine f_true

   ! Mirror JFNK_SET_PTC_SCALE_V1: freeze s_i = max(|F0_i|, floor) from
   ! the start-point residual and set the first anchor to x0 = 0 (the
   ! relative-unknown baseline).
   subroutine set_ptc_like_production(n, f0v)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: f0v(n)
      double precision :: fmax, floor_use
      integer :: i
      if(allocated(jfnk_ptc_scale)) deallocate(jfnk_ptc_scale)
      if(allocated(jfnk_ptc_x0)) deallocate(jfnk_ptc_x0)
      allocate(jfnk_ptc_scale(n), jfnk_ptc_x0(n))
      fmax = maxval(abs(f0v))
      floor_use = max(1.0d-30, 1.0d-8 * fmax)
      do i = 1, n
         jfnk_ptc_scale(i) = max(abs(f0v(i)), floor_use)
      end do
      jfnk_ptc_scale_ready = .true.
      jfnk_ptc_x0 = 0d0
      jfnk_ptc_anchor_ready = .true.
   end subroutine set_ptc_like_production

end program test_newton_core


! Residual callback: y = F(xbase + x) plus the pseudo-time term,
! added exactly as JFNK_RESIDUAL_V1 adds it.
subroutine f_test(n, x, y)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   use test_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   double precision :: xa(n), d(n)
   integer :: i
   xa = xbase + x
   d = xa - xs
   do i = 1, n
      y(i) = a(i) * (exp(xa(i)) - exp(xs(i)))
      if(i > 1) y(i) = y(i) + b * (d(i) - d(i-1))
      if(i < n) y(i) = y(i) + b * (d(i) - d(i+1))
   end do
   if(jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
      jfnk_ptc_anchor_ready) then
      if(allocated(jfnk_ptc_x0) .and. allocated(jfnk_ptc_scale)) then
         if(size(jfnk_ptc_x0) == n .and. size(jfnk_ptc_scale) == n) then
            do i = 1, n
               y(i) = y(i) + (x(i) - jfnk_ptc_x0(i)) * &
                  jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
            end do
         end if
      end if
   end if
end subroutine f_test


! JVP callback: forward finite difference of f at the module base
! point (jfnk_x, jfnk_fx), mirroring JFNK_JVP_V1.
subroutine df_test(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: h, xpert(n), fpert(n), vinf
   vinf = maxval(abs(v))
   if(vinf == 0d0) then
      jv = 0d0
      return
   end if
   h = sqrt(epsilon(1d0)) * max(1d0, maxval(abs(jfnk_x))) / vinf
   xpert = jfnk_x + h * v
   call f_test(n, xpert, fpert)
   jv = (fpert - jfnk_fx) / h
end subroutine df_test


! Preconditioner: v := v / (diagJ + ptcdiag).  diagJ is the exact
! diagonal a_i*exp(x_i) + 2b lagged at the solve start (like the BA
! LU); ptcdiag = s/dtau while PTC is active.
subroutine dfp_test(n, v)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready
   use test_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: i
   double precision :: ptc
   do i = 1, n
      ptc = 0d0
      if(jfnk_ptc_active .and. jfnk_ptc_scale_ready) &
         ptc = jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
      v(i) = v(i) / (diagJ(i) + ptc)
   end do
end subroutine dfp_test


subroutine sub_test
   implicit none
end subroutine sub_test


! Linear matvec for the T4 tres tests: jv = A v, plus the constant
! offset c0v when matvec_bias is set.  The offset makes the action
! affine -- each Arnoldi product carries c0 while the accumulated
! step does not, so the least-squares residual under-reports the
! explicit ||b - matvec(x)|| the tres gate measures.
subroutine df_lin(n, v, jv)
   use test_prob_mod
   use, intrinsic :: ieee_arithmetic
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   matvec_ct = matvec_ct + 1
   jv = matmul(amat, v)
   if(matvec_bias) jv = jv + c0v
   ! T4f: land a NaN on the tres measurement call.  Deliberately does
   ! NOT raise jfnk_jvp_failed -- the guard under test is the solver's
   ! own non-finite check on the measured residual.
   if(matvec_ct == matvec_nan_at) &
      jv(1) = ieee_value(0d0, ieee_quiet_nan)
end subroutine df_lin


subroutine dfp_id(n, v)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
end subroutine dfp_id


! Nontrivial left preconditioner for T4d/T4e: a per-row varying
! diagonal (not a multiple of the identity), so the tres psolve branch
! is exercised the same way the production left-pc path uses it.
subroutine dfp_jac(n, v)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: i
   do i = 1, n
      v(i) = v(i) / (2d0 + 0.5d0 * dble(i) / dble(n))
   end do
end subroutine dfp_jac


double precision function dp_test(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_test = dot_product(a, b)
end function dp_test
