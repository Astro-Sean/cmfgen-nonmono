!----------------------------------------------------------------------
! TEST_STAG_GRADED: white-box check of WP-2 (jfnk.md) - the graded
! Krylov-stagnation response under [JFNK_STAG_POLICY='GRADED'].
!
! Problem: F(x) = A (x - x*), A = I + 20*S with S the cyclic
! permutation, x* = 1. The map is linear, J = A exactly, and the
! trial F_pred is exact, so every step the solver keeps is
! commitable when it reduces |F|. A's spectrum is a circle of
! radius 20 about 1, so a starved solve (m = 3, restart cap 0)
! contracts the residual to ~0.99*ref_norm before exiting ginfo=1:
! inside the WP-2.1(a) band [0.9,1)*ref that the stock 0.9-tier
! arm misses, which is exactly the ABORT/GRADED separation this
! test needs.
!
!   run A: jfnk_do_stag_graded = F (default ABORT) -> Newton aborts
!          on the first stall (info = 3, no iterate, no refresh
!          consumption). The "before" behaviour.
!   run B: GRADED -> the inexact direction is kept (ginfo -> 0),
!          the descent gate passes (<F,As> = |As|^2 > 0 for this
!          min-residual step), and at least one iterate commits.
!          JFNK_STUB_NFACTOR stays 0: arm (a) must not spend a
!          refresh.
!   run C: GRADED with a matvec that fails on every 4th call -- the
!          exit true-res measurement (call 4 of each 4-call solve:
!          3 Arnoldi + 1 exit measure) never lands, so true_res is
!          unverifiable and arm (b) fires: consume the refresh flag
!          on the rejection path via JFNK_FACTOR_PRECOND_V1, tighten
!          eta 10x, re-solve. Two retries then the abort stands:
!          NFACTOR must be exactly 2.
!----------------------------------------------------------------------
module stag_prob_mod
   implicit none
   integer, parameter :: NN = 12
   integer :: n_calls = 0
   logical :: fail_on_measure = .false.
end module stag_prob_mod


program test_stag_graded
   use jfnk_solver_v1_mod
   use jfnk_cmfgen_int_mod, only: JFNK_STUB_NFACTOR
   use stag_prob_mod
   implicit none
   external :: f_sg, df_sg, df_sg_fail, dfp_sg, sub_sg
   double precision, external :: dp_sg
   integer :: info, nfail, luer, i, nits_a, nits_b
   double precision :: x0(NN), tol_t, fnrm0
   !
   luer = 6
   nfail = 0
   do i = 1, NN
      x0(i) = 1d0 + 0.05d0 * dble(i)
   end do
   if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
   !
   jfnk_ptc_active = .false.
   jfnk_do_line_search = .true.
   jfnk_ls_nonmono = .false.
   jfnk_do_cauchy = .false.
   jfnk_gmres_max_restart = 0
   ! Budget = 2*m so j==m is reached with iterations to spare: the
   ! exit is the restart cap (ginfo=1), not the its-budget path
   ! (which returns a usable direction without a stall flag).
   jfnk_max_gmres_its = 6
   jfnk_do_stag_graded = .false.
   call f_sg(NN, x0, jfnk_fx)
   fnrm0 = dsqrt(dp_sg(NN, jfnk_fx, jfnk_fx))
   tol_t = 1.0d-14 * fnrm0
   !
   ! ---- run A: default ABORT ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   JFNK_STUB_NFACTOR = 0
   n_calls = 0
   fail_on_measure = .false.
   jfnk_x = x0
   info = 1
   call newton_hook_v1(f_sg, df_sg, dfp_sg, sub_sg, dp_sg, &
                       3, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   nits_a = jfnk_nits
   write(luer,'(A,I3,A,I4,A,I3)') 'runA: info=', info, &
      '  nits=', nits_a, '  nfactor=', JFNK_STUB_NFACTOR
   if(info == 3 .and. JFNK_STUB_NFACTOR == 0) then
      write(luer,'(A)') 'A1: PASS  ABORT policy aborts on an ' // &
         'in-band stall with no refresh spend'
   else
      write(luer,'(A)') 'A1: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---- run B: GRADED, verifiable direction (arm a) ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   jfnk_do_stag_graded = .true.
   JFNK_STUB_NFACTOR = 0
   n_calls = 0
   jfnk_x = x0
   info = 1
   call newton_hook_v1(f_sg, df_sg, dfp_sg, sub_sg, dp_sg, &
                       3, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   nits_b = jfnk_nits
   write(luer,'(A,I3,A,I4,A,I3)') 'runB: info=', info, &
      '  nits=', nits_b, '  nfactor=', JFNK_STUB_NFACTOR
   if(nits_b > nits_a) then
      write(luer,'(A,I4,A,I4)') 'B1: PASS  GRADED keeps stalling ' // &
         'directions that ABORT discards (nits B=', nits_b, &
         ' vs A=', nits_a, ')'
   else
      write(luer,'(A)') 'B1: FAIL  GRADED gained no iterates'
      nfail = nfail + 1
   end if
   if(JFNK_STUB_NFACTOR == 0) then
      write(luer,'(A)') 'B2: PASS  arm (a) spends no precond refresh'
   else
      write(luer,'(A)') 'B2: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---- run C: GRADED, unverifiable direction (arm b) ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   JFNK_STUB_NFACTOR = 0
   n_calls = 0
   fail_on_measure = .true.
   jfnk_x = x0
   info = 1
   call newton_hook_v1(f_sg, df_sg_fail, dfp_sg, sub_sg, dp_sg, &
                       3, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   write(luer,'(A,I3,A,I4,A,I3)') 'runC: info=', info, &
      '  nits=', jfnk_nits, '  nfactor=', JFNK_STUB_NFACTOR
   if(JFNK_STUB_NFACTOR == 2) then
      write(luer,'(A)') 'C1: PASS  refresh consumed on the ' // &
         'rejection path, capped at 2 retries'
   else
      write(luer,'(A)') 'C1: FAIL'
      nfail = nfail + 1
   end if
   if(info == 3) then
      write(luer,'(A)') 'C2: PASS  abort stands after two ' // &
         'failed re-solves'
   else
      write(luer,'(A)') 'C2: FAIL'
      nfail = nfail + 1
   end if
   !
   if(nfail == 0) then
      write(luer,'(A)') 'TEST_STAG_GRADED: ALL PASS'
   else
      write(luer,'(A,I3)') 'TEST_STAG_GRADED: FAILURES=', nfail
   end if
   if(nfail > 0) stop 1
end program test_stag_graded


subroutine f_sg(n, x, y)
   ! F(x) = (I + 20 S)(x - 1), S the forward cyclic permutation.
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   integer :: i
   do i = 1, n
      y(i) = (x(i) - 1d0) + 20d0 * (x(mod(i, n) + 1) - 1d0)
   end do
end subroutine f_sg


subroutine df_sg(n, v, jv)
   use stag_prob_mod, only: n_calls
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   integer :: i
   n_calls = n_calls + 1
   do i = 1, n
      jv(i) = v(i) + 20d0 * v(mod(i, n) + 1)
   end do
end subroutine df_sg


subroutine df_sg_fail(n, v, jv)
   ! Every 4th call is a solve's exit true-res measurement (3
   ! Arnoldi + 1 measure per starved solve); failing it leaves the
   ! direction unverifiable, which is what arm (b) keys on.
   use jfnk_solver_v1_mod, only: jfnk_jvp_failed
   use stag_prob_mod, only: n_calls, fail_on_measure
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   integer :: i
   n_calls = n_calls + 1
   do i = 1, n
      jv(i) = v(i) + 20d0 * v(mod(i, n) + 1)
   end do
   if(fail_on_measure .and. mod(n_calls, 4) == 0) &
      jfnk_jvp_failed = .true.
end subroutine df_sg_fail


subroutine dfp_sg(n, v)
   ! Identity preconditioner: M = I, so ref_norm = ||F|| and the
   ! measured true_res is the physical residual.
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
end subroutine dfp_sg


subroutine sub_sg
   implicit none
end subroutine sub_sg


double precision function dp_sg(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_sg = dot_product(a, b)
end function dp_sg
