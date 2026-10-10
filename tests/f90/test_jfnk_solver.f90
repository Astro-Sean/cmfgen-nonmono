! Standalone unit tests for the production Newton-Krylov solver
! (newton_hook_v1 + gmres_m_v1 in jfnk_solver_v1.f90).
!
! The module is compiled against a stub jfnk_cmfgen_int_mod so the
! real production code path (trust region, line search, EW forcing,
! PTC, Cauchy fallback) is exercised, not a copy.
!
! Test problems are synthetic with known roots so convergence can be
! verified exactly.

module jfnk_test_prob
   implicit none
   integer :: prob_kind = 1
   ! 1 = linear F = A(x-x*)   (exact J and exact PC)
   ! 2 = nonlinear F_i = x_i + sin(x_i) - c_i  (diag approx PC)
   ! 3 = residual returns NaN (failure handling)
   ! 4 = linear, Jacobi (diagonal) PC only
   double precision, allocatable :: amat(:,:), ainv(:,:)
   double precision, allocatable :: xstar(:), cvec(:)
   integer :: n_eta_calls = 0
   double precision :: eta_hist(200) = 0d0
   integer :: nits_hist(200) = 0
   ! Per-iterate residual/dtau history, recorded in tsub (which fires
   ! once at entry and once per accepted Newton iterate) and dumped
   ! as CSV for tests/python/plot_solver_tests.py.
   double precision, allocatable :: hist_res(:), hist_dtau(:)
   integer :: hist_n = 0
   logical :: hist_on = .false.
end module jfnk_test_prob


subroutine tf(n, x, fx)
   use, intrinsic :: ieee_arithmetic
   use jfnk_test_prob
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_scale, &
        jfnk_ptc_scale_ready, jfnk_ptc_x0, jfnk_ptc_dtau, &
        jfnk_ptc_anchor_ready
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: fx(n)
   integer :: i
   select case (prob_kind)
   case (1, 4)
      fx = matmul(amat, x - xstar)
   case (2)
      fx = x + sin(x) - cvec
   case default
      fx = ieee_value(0d0, ieee_quiet_nan)
   end select
   if (jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
       jfnk_ptc_anchor_ready) then
      do i = 1, n
         fx(i) = fx(i) + (x(i) - jfnk_ptc_x0(i)) * &
                jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
      end do
   end if
end subroutine tf


subroutine tdf(n, v, jv)
   use jfnk_test_prob
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_ptc_active, &
        jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_dtau
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   integer :: i
   select case (prob_kind)
   case (1, 4)
      jv = matmul(amat, v)
   case (2)
      do i = 1, n
         jv(i) = (1d0 + cos(jfnk_x(i))) * v(i)
      end do
   case default
      jv = 0d0
   end select
   ! The production JVP finite-differences G = F + ptc, so the exact
   ! JVP here must include the ptc diagonal to be consistent.
   if (jfnk_ptc_active .and. jfnk_ptc_scale_ready) then
      do i = 1, n
         jv(i) = jv(i) + v(i) * jfnk_ptc_scale(i) / &
                 max(jfnk_ptc_dtau, 1d-8)
      end do
   end if
end subroutine tdf


subroutine tdf_fd(n, v, jv)
   ! Finite-difference JVP replicating JFNK_JVP_V1's step choice:
   ! h = eps * max(1,|x|inf) / |v|inf, forward difference.
   use jfnk_test_prob
   use jfnk_solver_v1_mod, only: jfnk_x, jfnk_fx
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: h, xinf, vinf
   double precision, allocatable :: xp(:), fp(:)
   external :: tf
   vinf = maxval(abs(v))
   if (vinf < 1d-300) then
      jv = 0d0
      return
   end if
   xinf = max(1d0, maxval(abs(jfnk_x)))
   h = 1d-8 * xinf / vinf
   allocate(xp(n), fp(n))
   xp = jfnk_x + h * v
   call tf(n, xp, fp)
   jv = (fp - jfnk_fx) / h
   deallocate(xp, fp)
end subroutine tdf_fd


subroutine tdfp(n, v)
   use jfnk_test_prob
   use jfnk_solver_v1_mod, only: jfnk_x
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: i
   select case (prob_kind)
   case (1)
      v = matmul(ainv, v)            ! exact inverse
   case (2, 4)
      do i = 1, n
         v(i) = v(i) / (1d0 + cos(jfnk_x(i)))
      end do
   case default
      ! identity
   end select
end subroutine tdfp


double precision function tdp(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   tdp = dot_product(a, b)
end function tdp


subroutine tsub()
   use jfnk_test_prob
   use jfnk_solver_v1_mod, only: jfnk_eta, jfnk_nits, jfnk_tol, &
                               jfnk_ptc_dtau
   implicit none
   n_eta_calls = n_eta_calls + 1
   if (n_eta_calls <= 200) then
      eta_hist(n_eta_calls) = jfnk_eta
      nits_hist(n_eta_calls) = jfnk_nits
   end if
   if (hist_on) then
      hist_n = hist_n + 1
      if (hist_n > size(hist_res)) call grow_hist()
      hist_res(hist_n) = jfnk_tol
      hist_dtau(hist_n) = jfnk_ptc_dtau
   end if
end subroutine tsub


subroutine grow_hist()
   use jfnk_test_prob
   implicit none
   double precision, allocatable :: tr(:), td(:)
   call move_alloc(hist_res, tr)
   call move_alloc(hist_dtau, td)
   allocate(hist_res(2 * size(tr)), hist_dtau(2 * size(td)))
   hist_res(1:size(tr)) = tr
   hist_dtau(1:size(td)) = td
end subroutine grow_hist


subroutine dump_hist(fname, info_)
   use jfnk_test_prob
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: info_
   integer :: k, lu
   open(newunit=lu, file=fname, status='replace')
   write(lu,'(A,I0)') '# solver_info=', info_
   write(lu,'(A)') '# iter,res,dtau'
   do k = 1, hist_n
      write(lu,'(I0,A,ES16.8,A,ES16.8)') k - 1, ',', &
         hist_res(k), ',', hist_dtau(k)
   end do
   close(lu)
end subroutine dump_hist


program test_jfnk_solver
   use jfnk_solver_v1_mod
   use jfnk_test_prob
   implicit none
   external :: tf, tdf, tdfp, tsub, tdf_fd
   double precision, external :: tdp
   integer :: n, info, nits_max, luer
   integer :: npass, nfail
   double precision :: tol, gtol, del, mndl, mxdl
   double precision, allocatable :: x0(:)
   integer :: i, j, seed(8)
   double precision :: err, fnrm

   luer = 6
   npass = 0
   nfail = 0
   hist_on = .false.
   hist_n = 0
   seed = (/1,2,3,4,5,6,7,8/)
   call random_seed(put=seed)

   write(luer,'(A)') '================================================'
   write(luer,'(A)') ' JFNK solver standalone unit tests'
   write(luer,'(A)') '================================================'

   !----------------------------------------------------------------
   ! T1: linear problem, exact Jacobian + exact PC -> ~1 Newton step
   !----------------------------------------------------------------
   n = 8
   call build_linear(n)
   prob_kind = 1
   call reset_state(n)
   allocate(x0(n)); x0 = 0d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .false., info, 0)
   err = maxval(abs(jfnk_x - xstar))
   call check('T1 linear Newton converged (info=0)', info == 0)
   call check('T1 linear |x-x*| < 1e-8', err < 1d-8)
   call check('T1 linear <= 3 Newton its', jfnk_nits <= 3)

   !----------------------------------------------------------------
   ! T2: nonlinear problem, diag-approx PC -> converges
   !----------------------------------------------------------------
   n = 6
   call build_nonlinear(n)
   prob_kind = 2
   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0.3d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .false., info, 0)
   fnrm = sqrt(dot_product(jfnk_fx, jfnk_fx))
   call check('T2 nonlinear converged (info=0)', info == 0)
   call check('T2 nonlinear |F| < 1e-7', fnrm < 1d-7)

   !----------------------------------------------------------------
   ! T3: trust-region cap respected every Newton step
   !    (verified indirectly: del=0.1 forces multi-step contraction
   !     and solver must still converge)
   !----------------------------------------------------------------
   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0.3d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0.1d0, .false., info, 0)
   fnrm = sqrt(dot_product(jfnk_fx, jfnk_fx))
   call check('T3 TR-limited converged (info=0)', info == 0)
   call check('T3 TR-limited |F| < 1e-7', fnrm < 1d-7)
   call check('T3 TR forced >1 Newton it', jfnk_nits > 1)

   !----------------------------------------------------------------
   ! T4: EW forcing eta stays inside [eta_min, eta_max]
   !----------------------------------------------------------------
   call check('T4 eta >= eta_min on all its', &
              minval(eta_hist(1:n_eta_calls)) >= jfnk_eta_min - 1d-12)
   call check('T4 eta <= eta_max on all its', &
              maxval(eta_hist(1:n_eta_calls)) <= jfnk_eta_max + 1d-12)

   !----------------------------------------------------------------
   ! T5: PTC active -> converges, dtau grows above dtau0
   !----------------------------------------------------------------
   prob_kind = 2
   call reset_state(n)
   jfnk_ptc_active = .true.
   allocate(jfnk_ptc_scale(n)); jfnk_ptc_scale = 1d0
   jfnk_ptc_scale_ready = .true.
   jfnk_ptc_anchor_ready = .false.
   jfnk_ptc_dtau = 1d0
   deallocate(x0); allocate(x0(n)); x0 = 0.3d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .true., info, 0)
   fnrm = sqrt(dot_product(jfnk_fx, jfnk_fx))
   call check('T5 PTC converged (info=0)', info == 0)
   call check('T5 PTC |F| < 1e-7', fnrm < 1d-7)
   call check('T5 PTC dtau grew', jfnk_ptc_dtau > 1d0)
   jfnk_ptc_active = .false.
   jfnk_ptc_scale_ready = .false.

   !----------------------------------------------------------------
   ! T6: ptc strip helpers return plain F exactly
   !----------------------------------------------------------------
   block
      double precision :: xs(4), gs(4)
      jfnk_ptc_active = .true.
      jfnk_ptc_scale_ready = .true.
      if (allocated(jfnk_ptc_scale)) deallocate(jfnk_ptc_scale)
      allocate(jfnk_ptc_scale(4)); jfnk_ptc_scale = 2d0
      jfnk_ptc_dtau = 4d0
      call jfnk_ptc_set_anchor_v1(4, (/0d0,1d0,2d0,3d0/))
      xs = (/0.5d0, 1.5d0, 2.5d0, 3.5d0/)
      gs = xs**2 + (xs - jfnk_ptc_x0) * jfnk_ptc_scale / jfnk_ptc_dtau
      call jfnk_ptc_strip_v1(4, xs, gs)
      call check('T6 ptc_strip returns F', &
                 maxval(abs(gs - xs**2)) < 1d-14)
      jfnk_ptc_active = .false.
      jfnk_ptc_scale_ready = .false.
   end block

   !----------------------------------------------------------------
   ! T7: NaN residual -> solver reports failure, does not hang
   !----------------------------------------------------------------
   n = 4
   call build_linear(n)
   prob_kind = 3
   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .false., info, 0)
   call check('T7 NaN residual -> info /= 0', info /= 0)

   !----------------------------------------------------------------
   ! T8: linear problem, Jacobi PC only -> converges (slower)
   !----------------------------------------------------------------
   n = 8
   call build_linear(n)
   prob_kind = 4
   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0d0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .false., info, 0)
   err = maxval(abs(jfnk_x - xstar))
   call check('T8 Jacobi-PC converged (info=0)', info == 0)
   write(luer,'(A,ES12.4)') '   T8 |x-x*| =', err
   call check('T8 Jacobi-PC |x-x*| < 1e-6', err < 1d-6)

   !----------------------------------------------------------------
   ! T9: eval counters incremented
   !----------------------------------------------------------------
   call check('T9 f-eval counter > 0', jfnk_total_f_evals > 0)
   call check('T9 jvp counter > 0', jfnk_total_jvp_evals > 0)

   !----------------------------------------------------------------
   ! T10: nonlinear solve with FD-JVP (production path uses finite
   ! differences, not an exact matvec) -> still converges.
   !----------------------------------------------------------------
   n = 6
   prob_kind = 2
   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0.3d0
   call run_newton(tf, tdf_fd, tdfp, tsub, tdp, n, x0, 0d0, .false., &
                   info, 0)
   fnrm = sqrt(dot_product(jfnk_fx, jfnk_fx))
   call check('T10 FD-JVP converged (info=0)', info == 0)
   call check('T10 FD-JVP |F| < 1e-7', fnrm < 1d-7)

   !----------------------------------------------------------------
   ! Convergence-history capture: sub() fires once at entry and once
   ! per accepted Newton iterate, so hist_res is the |F| trajectory.
   ! Written as CSV for tests/python/plot_solver_tests.py.
   !----------------------------------------------------------------
   allocate(hist_res(64), hist_dtau(64))
   n = 6
   prob_kind = 2

   call reset_state(n)
   deallocate(x0); allocate(x0(n)); x0 = 0.3d0
   hist_on = .true.; hist_n = 0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0d0, .false., info, 0)
   hist_on = .false.
   call dump_hist('solver_hist_near.csv', info)

   call reset_state(n)
   jfnk_ptc_active = .true.
   if (allocated(jfnk_ptc_scale)) deallocate(jfnk_ptc_scale)
   allocate(jfnk_ptc_scale(n)); jfnk_ptc_scale = 1d0
   jfnk_ptc_scale_ready = .true.
   jfnk_ptc_anchor_ready = .false.
   jfnk_ptc_dtau = 1d0
   x0 = 3.0d0
   hist_on = .true.; hist_n = 0
   call run_newton(tf, tdf, tdfp, tsub, tdp, n, x0, 0.1d0, .true., info, 0)
   hist_on = .false.
   call dump_hist('solver_hist_far_ptc.csv', info)
   jfnk_ptc_active = .false.
   jfnk_ptc_scale_ready = .false.

   write(luer,'(A)') '================================================'
   write(luer,'(A,I3,A,I3)') ' RESULTS: ', npass, ' passed, ', nfail, ' failed'
   write(luer,'(A)') '================================================'
   if (nfail > 0) stop 1

contains

   subroutine check(name, ok)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      if (ok) then
         npass = npass + 1
         write(luer,'(A,A)') '  PASS  ', name
      else
         nfail = nfail + 1
         write(luer,'(A,A)') '  FAIL  ', name
      end if
   end subroutine check

   subroutine reset_state(nn)
      integer, intent(in) :: nn
      if (allocated(jfnk_x)) deallocate(jfnk_x)
      if (allocated(jfnk_fx)) deallocate(jfnk_fx)
      allocate(jfnk_x(nn), jfnk_fx(nn))
      jfnk_x = 0d0
      jfnk_fx = 0d0
      jfnk_fx_fresh = .false.
      jfnk_tol = 0d0
      jfnk_tol_prev = 0d0
      jfnk_nits = 0
      jfnk_gits = 0
      jfnk_del = 0d0
      jfnk_use_left_pc = .true.
      jfnk_max_gmres_its = 0
      jfnk_ptc_active = .false.
      jfnk_ptc_scale_ready = .false.
      jfnk_ptc_anchor_ready = .false.
      jfnk_ptc_dtau = 1d0
      jfnk_jvp_failed = .false.
      jfnk_pc_clip_event = .false.
      jfnk_need_precond_refresh = 0
      n_eta_calls = 0
   end subroutine reset_state

   subroutine run_newton(f_, df_, dfp_, sub_, dp_, nn, xinit, del_in, ptc, info_, ilvl)
      implicit none
      external :: f_, df_, dfp_, sub_
      double precision, external :: dp_
      integer, intent(in) :: nn, ilvl
      double precision, intent(in) :: xinit(nn), del_in
      logical, intent(in) :: ptc
      integer, intent(inout) :: info_
      jfnk_x = xinit
      info_ = ilvl
      call newton_hook_v1(f_, df_, dfp_, sub_, dp_, 10, nn, &
           1d-3, 1d-10, del_in, 1d-20, 1d20, 30, info_, luer)
   end subroutine run_newton

   subroutine build_linear(nn)
      integer, intent(in) :: nn
      double precision :: r
      if (allocated(amat)) deallocate(amat, ainv, xstar)
      allocate(amat(nn,nn), ainv(nn,nn), xstar(nn))
      call random_number(amat)
      amat = 0.2d0*amat
      do i = 1, nn
         amat(i,i) = 2d0 + amat(i,i)
      end do
      call random_number(xstar)
      xstar = 0.5d0 + xstar
      ! exact inverse via LU-free path: solve A*Ainv = I by elimination
      ainv = 0d0
      do i = 1, nn
         ainv(i,i) = 1d0
      end do
      call invert(nn)
   end subroutine build_linear

   subroutine invert(nn)
      integer, intent(in) :: nn
      integer :: k
      double precision :: aug(nn, 2*nn)
      aug(:,1:nn) = amat
      aug(:,nn+1:2*nn) = ainv
      do k = 1, nn
         aug(k,:) = aug(k,:) / aug(k,k)
         do j = 1, nn
            if (j /= k) aug(j,:) = aug(j,:) - aug(j,k)*aug(k,:)
         end do
      end do
      ainv = aug(:,nn+1:2*nn)
   end subroutine invert

   subroutine build_nonlinear(nn)
      integer, intent(in) :: nn
      if (allocated(cvec)) deallocate(cvec)
      allocate(cvec(nn))
      call random_number(cvec)
      cvec = 0.4d0 + cvec   ! roots near x in [0.4,1.4]
   end subroutine build_nonlinear

end program test_jfnk_solver
