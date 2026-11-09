!----------------------------------------------------------------------
! TEST_GMRES_EDGE: edge-path coverage for the GMRES(m) inner solve and
! the More-Sorensen hookstep inside jfnk_solver_v1.f90 (production .o).
!
! The routines under test are externals in the production object, so
! this driver calls them directly with manufactured Hessenberg systems
! and callback triples -- no Newton driver, no model state.
!
!   E1: hook_step_v1 with j=0 returns the documented zero step
!       (y=0, pred_res=beta) instead of calling DGESVD on an empty
!       problem.
!   E2: hook_step_v1 on a well-conditioned H reduces the model
!       residual (pred_res < beta) and returns a finite step.
!   E3: hook_step_v1 on a rank-deficient H (a zero column -> a zero
!       singular value) must not inject Inf/NaN into the step.
!   E4: hook_step_v1 honours the trust radius: ||y|| <= del.
!   E5: hook_step_g_v1 with G = I is bitwise identical to the
!       Euclidean hook_step_v1 (the G-metric machinery reduces to
!       the plain solve when the metric is the identity).
!   E6: hook_step_g_v1 with an indefinite G drives DPOTRF to
!       info /= 0; the documented fallback is the Euclidean
!       hook_step_v1, which the output must match bitwise.
!   E7: gmres_m_v1 on the identity operator with an identity
!       preconditioner converges in one Arnoldi step; the
!       accumulated step x equals b.
!   E8: gmres_m_v1 is deterministic: the same problem twice gives
!       bitwise-identical iterates.
!   E9: gmres_m_v1 on a cyclic-shift operator (minimal polynomial
!       degree n) cannot converge inside m=4 with its=6 -- the
!       solve must return inside the iteration budget with
!       info = 2 (budget exhausted) or 1 (stagnation), a finite
!       step, and a set jfnk_exit_reason.
!----------------------------------------------------------------------
module gm_edge_op_mod
   implicit none
   integer, parameter :: NN = 32
   ! Cyclic-shift operator storage is trivial (index map only); the
   ! identity case needs no state.
   logical :: use_shift = .false.
end module gm_edge_op_mod

program test_gmres_edge
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use jfnk_solver_v1_mod
   use gm_edge_op_mod
   implicit none

   external :: mv_edge, ps_edge
   double precision, external :: dp_euclid

   integer, parameter :: M = 4
   double precision :: h(M+1, M), v(NN, M+1)
   double precision :: x(NN), b(NN), x2(NN)
   double precision :: res, del, pred_res, res_stall
   integer :: its, info
   integer :: i, nfail
   logical :: lok

   nfail = 0

   ! Production scalar defaults the GMRES body reads (VADAT-set values
   ! in the real run; fixed here so the driver is self-contained).
   jfnk_use_left_pc     = .true.
   jfnk_lgmres_on       = .false.
   jfnk_gdefl_on        = .false.
   jfnk_recyc_on        = .false.
   jfnk_cycstall_on     = .false.
   jfnk_gmres_max_restart = 1
   jfnk_gmres_stag_tol  = 1.0d-3
   jfnk_gmres_tres_fac  = 1.0d0
   jfnk_precond_fails   = .false.
   jfnk_jvp_failed      = .false.
   jfnk_pc_clip_event   = .false.

   call run_hook_checks(nfail)

   ! ---- E7: identity operator converges in one step -----------------
   use_shift = .false.
   do i = 1, NN
      b(i) = sin(dble(i) * 0.37d0) + 0.5d0
   end do
   x = 0d0
   res = 1.0d-10
   del = 1.0d6
   its = 20
   info = 0
   call gmres_m_v1(M, NN, x, b, mv_edge, ps_edge, dp_euclid, &
                   h, v, res, del, pred_res, its, info, -1)
   lok = (info == 0) .and. &
         maxval(abs(x - b)) < 1.0d-9 .and. &
         maxval(abs(x)) < huge(1d0)
   call report('E7 identity operator solves x=b', lok, nfail)

   ! ---- E8: determinism ---------------------------------------------
   x2 = 0d0
   res = 1.0d-10
   its = 20
   info = 0
   call gmres_m_v1(M, NN, x2, b, mv_edge, ps_edge, dp_euclid, &
                   h, v, res, del, pred_res, its, info, -1)
   lok = all(x2 == x)
   call report('E8 GMRES is bitwise deterministic', lok, nfail)

   ! ---- E9: non-converging operator inside the budget ---------------
   use_shift = .true.
   x = 0d0
   res = 1.0d-12
   its = 6
   info = 0
   jfnk_exit_reason = 'unset'
   call gmres_m_v1(M, NN, x, b, mv_edge, ps_edge, dp_euclid, &
                   h, v, res, del, pred_res, its, info, -1)
   lok = (info == 1 .or. info == 2) .and. &
         all(ieee_is_finite(x)) .and. &
         jfnk_exit_reason /= 'unset'
   call report('E9 non-converging operator exits inside budget', &
               lok, nfail)
   res_stall = res

   ! ---- E10: LGMRES augmentation engages and beats the stall --------
   ! LGMRES is right-preconditioning only (module doc: under
   ! JFNK_LEFTPC=T the flag is ignored), so flip the PC side here.
   ! The cyclic shift is a pathological stall operator (eigenvalues
   ! uniform on the unit circle): augmentation cannot fully fix it,
   ! but the retained correction slots must populate and the
   ! augmented least-squares must beat plain restarted GMRES on the
   ! same operator. Measured: res 0.555 vs 0.930, lg_n = 8.
   jfnk_use_left_pc   = .false.
   jfnk_lgmres_on     = .true.
   jfnk_lgmres_kmax   = 8
   jfnk_aug_kmax      = 8
   jfnk_gmres_max_restart = 20
   x = 0d0
   res = 1.0d-12
   its = 80
   info = 0
   jfnk_exit_reason = 'unset'
   call gmres_m_v1(M, NN, x, b, mv_edge, ps_edge, dp_euclid, &
                   h, v, res, del, pred_res, its, info, -1)
   lok = all(ieee_is_finite(x)) .and. &
         jfnk_lgmres_n == 8 .and. &
         (info == 0 .or. res < 0.65d0 * res_stall)
   call report('E10 LGMRES slots fill, beats restart stall', &
               lok, nfail)

   ! ---- E11: GMRES-DR deflated restart on the same problem ----------
   jfnk_lgmres_on   = .false.
   jfnk_gdefl_on    = .true.
   jfnk_gdefl_kmax  = 4
   x = 0d0
   res = 1.0d-12
   its = 80
   info = 0
   call gmres_m_v1(M, NN, x, b, mv_edge, ps_edge, dp_euclid, &
                   h, v, res, del, pred_res, its, info, -1)
   lok = all(ieee_is_finite(x)) .and. &
         (info == 0 .or. res < 0.5d0 * res_stall)
   call report('E11 GDEFL deflation improves stalled solve', &
               lok, nfail)

   if(nfail == 0) then
      write(*,*) 'TEST_GMRES_EDGE: all checks passed'
   else
      write(*,*) 'TEST_GMRES_EDGE:', nfail, 'FAILURES'
      stop 1
   end if

contains

   subroutine report(name, ok, nfail)
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      integer, intent(inout) :: nfail
      if(ok) then
         write(*,'(A,A)') '  PASS  ', name
      else
         write(*,'(A,A)') '  FAIL  ', name
         nfail = nfail + 1
      end if
   end subroutine report

   subroutine run_hook_checks(nfail)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(inout) :: nfail
      integer, parameter :: JM = 3, JP = 4
      double precision :: hh(JP, JM), gg(JM, JM)
      double precision :: y(JM), y_ref(JM), pr, pr_ref, beta, del
      integer :: i, j
      logical :: ok

      beta = 2.5d0
      ! Well-conditioned upper Hessenberg from a fake Arnoldi step.
      ! hook_step_v1's m is the workspace bound: h is (m+1, j), so with
      ! a (4,3) array m = JM = 3.
      hh = 0d0
      hh(1,1) = 1.0d0; hh(2,1) = 0.5d0
      hh(1,2) = 0.3d0; hh(2,2) = 0.8d0; hh(3,2) = 0.4d0
      hh(1,3) = 0.1d0; hh(2,3) = 0.2d0; hh(3,3) = 0.6d0; hh(4,3) = 0.7d0

      ! E1: j=0 guard. The y dummy has extent j=0 here, so only the
      ! pred_res = beta contract (no least-squares solve attempted)
      ! is observable from the caller side.
      del = 1.0d0
      call hook_step_v1(0, hh, JM, beta, del, y, pr)
      ok = pr == beta
      call report('E1 hook_step_v1 j=0 returns zero step', ok, nfail)

      ! E2: unconstrained LSQ reduces the model residual.
      del = 1.0d6
      call hook_step_v1(JM, hh, JM, beta, del, y, pr)
      ok = all(ieee_is_finite(y)) .and. pr < beta
      call report('E2 hook_step_v1 solves LSQ, residual drops', &
                  ok, nfail)

      ! E3: rank-deficient H (zero column 2) stays finite.
      hh(:,2) = 0d0
      call hook_step_v1(JM, hh, JM, beta, del, y, pr)
      ok = all(ieee_is_finite(y)) .and. ieee_is_finite(pr)
      call report('E3 rank-deficient H gives finite step', ok, nfail)

      ! E4: trust-radius bound respected.
      hh(2,2) = 0.8d0; hh(3,2) = 0.4d0
      del = 0.05d0
      call hook_step_v1(JM, hh, JM, beta, del, y, pr)
      ok = all(ieee_is_finite(y)) .and. &
           sqrt(sum(y*y)) <= del * (1d0 + 1d-8)
      call report('E4 hookstep honours ||y|| <= del', ok, nfail)

      ! E5/E6: G-metric variant vs Euclidean reference.
      del = 0.3d0
      call hook_step_v1(JM, hh, JM, beta, del, y_ref, pr_ref)

      ! G = I: the 1e-12 Gram ridge perturbs the Cholesky factor, so
      ! the G-metric step can only agree to ~1e-11, not bitwise.
      gg = 0d0
      do i = 1, JM
         gg(i,i) = 1d0
      end do
      call hook_step_g_v1(JM, hh, JM, beta, del, gg, JM, y, pr)
      ok = all(abs(y - y_ref) <= 1d-9 * (1d0 + abs(y_ref))) .and. &
           abs(pr - pr_ref) <= 1d-9 * (1d0 + abs(pr_ref))
      call report('E5 G=I hookstep matches Euclidean', ok, nfail)

      ! Indefinite G: DPOTRF fails on the first non-positive pivot,
      ! and the documented fallback calls hook_step_v1 on the
      ! original H -- bitwise identical to the reference above.
      gg = 0d0
      do i = 1, JM
         gg(i,i) = -1d0
      end do
      call hook_step_g_v1(JM, hh, JM, beta, del, gg, JM, y, pr)
      ok = all(y == y_ref) .and. pr == pr_ref
      call report('E6 indefinite G falls back to Euclidean hook', &
                  ok, nfail)
   end subroutine run_hook_checks

end program test_gmres_edge

!----------------------------------------------------------------------
! Callbacks: matvec applies either the identity or a cyclic shift;
! psolve is the identity; dp is the Euclidean inner product.
!----------------------------------------------------------------------
subroutine mv_edge(n, vin, jv)
   use gm_edge_op_mod, only: use_shift
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: vin(n)
   double precision, intent(out) :: jv(n)
   if(use_shift) then
      jv(1) = vin(n)
      jv(2:n) = vin(1:n-1)
   else
      jv = vin
   end if
end subroutine mv_edge

subroutine ps_edge(n, vin)
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: vin(n)
end subroutine ps_edge

double precision function dp_euclid(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_euclid = dot_product(a, b)
end function dp_euclid
