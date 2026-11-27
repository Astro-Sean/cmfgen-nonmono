!----------------------------------------------------------------------
! TEST_BUDGET_STOP: white-box check of WP-3 (jfnk.md) - the
! budget-exhaustion graceful stop. The residual callback returns real
! F values for the first BS_BUDGET calls, then flips
! JFNK_EVAL_BUDGET_HIT and returns the production sentinel
! (Y = 1D20) -- the same pair the interface writes when the eval/wall
! budget runs out inside JFNK_RESIDUAL_V1.
!
!   run A: budget dies on call BS_BUDGET > 1 -- mid-solve, inside the
!          first GMRES matvec sequence. The solver must exit info=2
!          with jfnk_exit_reason='budget_stop' and leave jfnk_x/jfnk_fx
!          holding the last accepted (or entry) iterate, NOT the
!          sentinel trial. The driver maps info=2 to the best-iterate
!          commit path, so info=3 here would be a hard abort.
!   run B: budget dies on call 1 -- the entry eval itself is the
!          sentinel. Same graceful contract: info=2, entry state kept.
!
! With the flag untested the sentinel's |Y|~1e22 trial would land in
! the NaN/contraction guards or be trusted as a real F -- the old
! silent misroute this WP removes.
!----------------------------------------------------------------------
module bs_prob_mod
   implicit none
   integer, parameter :: NN = 12
   integer :: n_calls = 0
   integer :: bs_budget = 0
end module bs_prob_mod

subroutine f_bs(n, x, y)
   use bs_prob_mod
   use jfnk_cmfgen_int_mod, only: JFNK_EVAL_BUDGET_HIT
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   integer :: i
   n_calls = n_calls + 1
   ! bs_budget = number of real calls allowed; every later call is
   ! the sentinel pair (flag + 1D20), matching JFNK_RESIDUAL_V1.
   if(n_calls > bs_budget) then
      JFNK_EVAL_BUDGET_HIT = .true.
      y = 1.0d20
      return
   end if
   ! Same contracting map as test_stag_graded: F = A(x - 1),
   ! A = I + 20*S with S the cyclic permutation.
   y(1)    = x(1) - 1.0d0 + 20.0d0 * (x(NN) - 1.0d0)
   do i = 2, n
      y(i) = x(i) - 1.0d0 + 20.0d0 * (x(i-1) - 1.0d0)
   end do
end subroutine f_bs

subroutine df_bs(n, v, jv)
   use bs_prob_mod
   use jfnk_cmfgen_int_mod, only: JFNK_EVAL_BUDGET_HIT
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   integer :: i
   ! The JVP re-evals the residual in production, so a matvec can be
   ! the call that hits the budget -- count it on the same counter.
   n_calls = n_calls + 1
   if(n_calls > bs_budget) then
      JFNK_EVAL_BUDGET_HIT = .true.
      jv = v
      return
   end if
   jv(1) = v(1) + 20.0d0 * v(NN)
   do i = 2, n
      jv(i) = v(i) + 20.0d0 * v(i-1)
   end do
end subroutine df_bs

subroutine dfp_bs(n, v)
   use bs_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   ! Preconditioner in place: defect-correction iterate
   ! w <- w + (v - A w) starting at w = v. Not an inverse of
   ! I + 20S -- only a usable left preconditioner so GMRES sees a
   ! nontrivial operator.
   integer :: i, it
   double precision :: w(NN), r(NN), rhs(NN)
   rhs = v
   w = v
   do it = 1, 60
      r(1) = rhs(1) - (w(1) + 20.0d0 * w(NN))
      do i = 2, n
         r(i) = rhs(i) - (w(i) + 20.0d0 * w(i-1))
      end do
      if(maxval(abs(r)) < 1d-12 * max(1d0, maxval(abs(rhs)))) exit
      w = w + r
   end do
   v = w
end subroutine dfp_bs

subroutine sub_bs
   implicit none
   ! commit callback -- nothing to persist
end subroutine sub_bs

double precision function dp_bs(n, a, b)
   use bs_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_bs = dot_product(a, b)
end function dp_bs

program test_budget_stop
   use jfnk_solver_v1_mod
   use jfnk_cmfgen_int_mod, only: JFNK_EVAL_BUDGET_HIT, JFNK_EVAL_CAP
   use bs_prob_mod
   implicit none
   external :: f_bs, df_bs, dfp_bs, sub_bs
   double precision, external :: dp_bs
   integer :: info, nfail, luer, i
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
   jfnk_gmres_max_restart = 4
   jfnk_do_stag_graded = .false.
   !
   n_calls = 0
   bs_budget = 1000      ! measuring call must be a real eval
   call f_bs(NN, x0, jfnk_fx)
   fnrm0 = dsqrt(dp_bs(NN, jfnk_fx, jfnk_fx))
   tol_t = 1.0d-14 * fnrm0
   !
   ! ---- run A: sentinel inside the first GMRES solve ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   JFNK_EVAL_BUDGET_HIT = .false.
   JFNK_EVAL_CAP = 0
   n_calls = 0
   ! entry eval is call 1, first GMRES matvec is call 2 -- budget
   ! dies on the fourth call, inside the first solve.
   bs_budget = 4
   jfnk_x = x0
   info = 1
   write(luer,'(A,ES12.4,A,ES12.4,A,L2)') 'preA: fnrm0=', fnrm0, &
      ' tol_t=', tol_t, ' fresh=', jfnk_fx_fresh
   call newton_hook_v1(f_bs, df_bs, dfp_bs, sub_bs, dp_bs, &
                       3, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   write(luer,'(A,I3,A,L2,A,A)') 'runA: info=', info, &
      '  hit=', JFNK_EVAL_BUDGET_HIT, '  reason=', jfnk_exit_reason
   if(info == 2) then
      write(luer,'(A)') 'A1: PASS  budget hit exits info=2 ' // &
         '(driver best-iterate path)'
   else
      write(luer,'(A,I3)') 'A1: FAIL  info=', info
      nfail = nfail + 1
   end if
   if(jfnk_exit_reason == 'budget_stop') then
      write(luer,'(A)') 'A2: PASS  exit_reason=budget_stop'
   else
      write(luer,'(A,A)') 'A2: FAIL  reason=', jfnk_exit_reason
      nfail = nfail + 1
   end if
   if(maxval(abs(jfnk_fx)) < 1d19) then
      write(luer,'(A)') 'A3: PASS  jfnk_fx holds a real iterate, ' // &
         'not the sentinel'
   else
      write(luer,'(A)') 'A3: FAIL  jfnk_fx carries the sentinel'
      nfail = nfail + 1
   end if
   !
   ! ---- run B: sentinel on the entry eval ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   JFNK_EVAL_BUDGET_HIT = .false.
   n_calls = 0
   bs_budget = 0          ! every call is the sentinel
   jfnk_x = x0
   info = 1
   call newton_hook_v1(f_bs, df_bs, dfp_bs, sub_bs, dp_bs, &
                       3, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   write(luer,'(A,I3,A,A)') 'runB: info=', info, &
      '  reason=', jfnk_exit_reason
   if(info == 2 .and. jfnk_exit_reason == 'budget_stop') then
      write(luer,'(A)') 'B1: PASS  entry-eval sentinel exits ' // &
         'info=2 budget_stop'
   else
      write(luer,'(A,I3,A,A)') 'B1: FAIL  info=', info, &
         '  reason=', jfnk_exit_reason
      nfail = nfail + 1
   end if
   if(maxval(abs(jfnk_x - x0)) < 1d-12) then
      write(luer,'(A)') 'B2: PASS  entry state returned untouched'
   else
      write(luer,'(A)') 'B2: FAIL  jfnk_x moved off the entry state'
      nfail = nfail + 1
   end if
   !
   if(nfail == 0) then
      write(luer,'(A)') 'TEST_BUDGET_STOP: ALL PASS'
   else
      write(luer,'(A,I3)') 'TEST_BUDGET_STOP: FAILURES =', nfail
   end if
end program test_budget_stop
