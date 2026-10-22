!----------------------------------------------------------------------
! TEST_PTC_SIGN: Newton-Krylov comparison of the signed vs unsigned
! PTC convention on a negative-diagonal variant of the stiff
! exponential problem (audit F1).  Same machinery as
! test_newton_core but the residual is negated
!   F_i = -a_i*(exp(x_i) - exp(xs_i)) - b*(d_i - d_{i-1})
!                                  - b*(d_i - d_{i+1})
! so J_ii = -a_i*exp(x_i) - 2b < 0 everywhere, the loss-dominated
! rate-row regime of the real BA.
!
! Two conventions, mirroring the production contract:
!   unsigned: G = F + sigma*(x - x0),  M_diag = d - (fold) + sigma
!   signed:   G = F - sigma*(x - x0),  M_diag = d - sigma
! (sgn = sign(d_i) = -1 uniformly here).  The signed convention is
! what JFNK_PTC_SIGNED implements: the residual PTC term, the
! preconditioner fold, and the apply rescale all carry the same
! sign so |d + shift| = |d| + sigma.  The toy preconditioner is a
! lagged Jacobi diagonal folded the same way the LU fold is.
!
! 10 deterministic starts; per start report Newton calls and
! success for each convention.  Prints a summary table -- the
! verdict (helps/hurts/no difference) is reported in the audit doc,
! not asserted, since either outcome is information.
!
! edited by Devin, 22-Oct-2026
!----------------------------------------------------------------------
module test_ptc_sign_mod
   implicit none
   integer, parameter :: NN = 60
   double precision :: a(NN), xs(NN), xbase(NN), diagJ(NN)
   double precision :: b
   logical :: ptc_signed = .false.
   integer :: n_f_evals = 0
end module test_ptc_sign_mod


program test_ptc_sign
   use jfnk_solver_v1_mod
   use test_ptc_sign_mod
   implicit none
   external :: f_neg, df_neg, dfp_neg, sub_neg
   double precision, external :: dp_neg
   integer :: i, info, it, ic, luer
   double precision :: x0(NN), f0(NN), ftrue(NN), xsol(NN), fmax
   double precision :: fnrm0, tol_t, ftru_nrm, maxerr
   integer :: evals(10,2), infos(10,2)
   double precision :: ferr(10,2), xerr(10,2)
   double precision :: pert(10)
   !
   luer = 6
   open(99, status='scratch')   ! solver's own log lines
   b = 0.3d0
   do i = 1, NN
      a(i)  = 1.0d-3 * 10d0**(6d0 * dble(i - 1) / dble(NN - 1))
      xs(i) = 1.5d0 * sin(1.7d0 * dble(i))
   end do
   if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
   ! deterministic perturbations covering near and far starts
   pert = [ 0.2d0, 0.5d0, -0.4d0, 1.0d0, -0.8d0, &
            1.5d0, -1.2d0, 2.0d0, -1.6d0, 2.5d0 ]
   !
   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 1d0
   !
   write(luer,'(A)') ' trial  conv     info  evals   |F_true|'// &
      '     max|x-xs|'
   do it = 1, 10
      do ic = 1, 2
         ptc_signed = (ic == 2)
         xbase = 0d0
         x0 = xs + pert(it)
         jfnk_x = x0
         n_f_evals = 0
         call f_neg(NN, jfnk_x, f0)
         fnrm0 = dsqrt(dp_neg(NN, f0, f0))
         tol_t = 1.0d-10 * fnrm0
         ! lagged diagonal of J_F (negative)
         diagJ = -a * exp(x0) - 2d0 * b
         ! freeze s_i = max(|F0_i|, floor) like JFNK_SET_PTC_SCALE_V1
         if(allocated(jfnk_ptc_scale)) deallocate(jfnk_ptc_scale)
         if(allocated(jfnk_ptc_x0)) deallocate(jfnk_ptc_x0)
         allocate(jfnk_ptc_scale(NN), jfnk_ptc_x0(NN))
         fmax = maxval(abs(f0))
         do i = 1, NN
            jfnk_ptc_scale(i) = max(abs(f0(i)), 1.0d-8 * fmax, 1.0d-30)
         end do
         jfnk_ptc_scale_ready = .true.
         jfnk_ptc_x0 = 0d0
         jfnk_ptc_anchor_ready = .true.
         info = 1
         call newton_hook_v1(f_neg, df_neg, dfp_neg, sub_neg, dp_neg, &
                            20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, &
                            1.0d99, 50, info, 99)
         xsol = jfnk_x
         call f_true_neg(NN, xsol, ftrue)
         ftru_nrm = dsqrt(dp_neg(NN, ftrue, ftrue))
         evals(it,ic) = n_f_evals
         infos(it,ic) = info
         ferr(it,ic) = ftru_nrm
         xerr(it,ic) = maxval(abs(xsol - xs))
         write(luer,'(I5,A,I8,I8,ES12.4,ES12.4)') it, &
            merge(' signed  ',' unsigned', ptc_signed), info, &
            n_f_evals, ftru_nrm, xerr(it,ic)
      end do
   end do
   !
   write(luer,'(A)') ' summary: evals unsigned / signed per trial'
   do it = 1, 10
      write(luer,'(I5,2I8)') it, evals(it,1), evals(it,2)
   end do
   write(luer,'(A,I6,A,I6)') ' total evals unsigned=', &
      sum(evals(:,1)), '  signed=', sum(evals(:,2))
   write(luer,'(A,I2,A,I2)') ' converged unsigned=', &
      count(infos(:,1) == 0), '  signed=', count(infos(:,2) == 0)
   write(luer,'(A)') 'TEST_PTC_SIGN DONE'

contains

   ! F without the PTC term, at absolute point x (negated problem).
   subroutine f_true_neg(n, x, y)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: y(n)
      integer :: i
      double precision :: d(n)
      d = x - xs
      do i = 1, n
         y(i) = -a(i) * (exp(x(i)) - exp(xs(i)))
         if(i > 1) y(i) = y(i) - b * (d(i) - d(i-1))
         if(i < n) y(i) = y(i) - b * (d(i) - d(i+1))
      end do
   end subroutine f_true_neg

end program test_ptc_sign


! Residual callback: G = F(x) + sgn*sigma*(x - x0), sgn = -1 under
! the signed convention, +1 unsigned.
subroutine f_neg(n, x, y)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   use test_ptc_sign_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   double precision :: xa(n), d(n), sgn
   integer :: i
   n_f_evals = n_f_evals + 1
   xa = xbase + x
   d = xa - xs
   do i = 1, n
      y(i) = -a(i) * (exp(xa(i)) - exp(xs(i)))
      if(i > 1) y(i) = y(i) - b * (d(i) - d(i-1))
      if(i < n) y(i) = y(i) - b * (d(i) - d(i+1))
   end do
   sgn = merge(-1.0d0, 1.0d0, ptc_signed)
   if(jfnk_ptc_active .and. jfnk_ptc_scale_ready .and. &
      jfnk_ptc_anchor_ready) then
      if(allocated(jfnk_ptc_x0) .and. allocated(jfnk_ptc_scale)) then
         if(size(jfnk_ptc_x0) == n .and. size(jfnk_ptc_scale) == n) then
            do i = 1, n
               y(i) = y(i) + sgn * (x(i) - jfnk_ptc_x0(i)) * &
                  jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
            end do
         end if
      end if
   end if
end subroutine f_neg


subroutine df_neg(n, v, jv)
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
   call f_neg(n, xpert, fpert)
   jv = (fpert - jfnk_fx) / h
end subroutine df_neg


! Jacobi preconditioner with the PTC shift folded the same way as
! the production LU fold: d + sgn*sigma.
subroutine dfp_neg(n, v)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready
   use test_ptc_sign_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: i
   double precision :: ptc, sgn
   sgn = merge(-1.0d0, 1.0d0, ptc_signed)
   do i = 1, n
      ptc = 0d0
      if(jfnk_ptc_active .and. jfnk_ptc_scale_ready) &
         ptc = jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
      v(i) = v(i) / (diagJ(i) + sgn * ptc)
   end do
end subroutine dfp_neg


subroutine sub_neg
   implicit none
end subroutine sub_neg


double precision function dp_neg(n, aa, bb)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: aa(n), bb(n)
   dp_neg = dot_product(aa, bb)
end function dp_neg
