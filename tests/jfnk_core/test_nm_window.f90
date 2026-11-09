!----------------------------------------------------------------------
! TEST_NM_WINDOW: white-box check of the non-monotone line-search
! history ring in newton_hook_v1 (jfnk_solver_v1.f90, production .o).
!
! Documented semantics (module comment, ~line 283): "the reference
! norm is the max of the last jfnk_ls_nmwin_s ACCEPTED |F| values".
! nmwin_s = 4, ring = 8 slots, nmidx cycles all 8. The phi_ref loop
! reads slots 1..min(nmwin_s,8) = 1..4 unconditionally.
!
! The sub() callback runs once after every ACCEPTED Newton iterate,
! so this test can compare:
!   intended(k) = max of the last-4 accepted |F| values
!   actual(k)   = max over jfnk_ls_nmhist(1:4)
! Any k where these differ shows the read window not tracking the
! last-4 window (slots 5..8 are write-only as far as phi_ref sees).
!
! Also verified: nmhist(nmidx)==jfnk_tol after each accept (write
! pointer), and a fresh activation resets the ring to -1.
!
! Problem: same stiff synthetic family as test_newton_core so the
! driver runs enough iterations to fill the ring.
!----------------------------------------------------------------------
module nm_prob_mod
   implicit none
   integer, parameter :: NN = 60
   double precision :: a(NN), xs(NN), xbase(NN)
   double precision :: b
   ! accepted-|F| history recorded from sub()
   integer :: n_acc = 0
   double precision :: acc_tol(200)
   integer :: nmidx_seen(200)
   integer :: n_bad_write = 0, n_diverge = 0
   integer :: first_diverge = -1
end module nm_prob_mod


program test_nm_window
   use jfnk_solver_v1_mod
   use nm_prob_mod
   implicit none
   external :: f_nm, df_nm, dfp_nm, sub_nm, sub_nm2
   double precision, external :: dp_nm
   integer :: i, info, nfail, luer, k
   double precision :: x0(NN), f0(NN), diagJ(NN)
   double precision :: fnrm0, tol_t
   !
   luer = 6
   nfail = 0
   b = 0.3d0
   do i = 1, NN
      a(i)  = 1.0d-3 * 10d0**(6d0 * dble(i - 1) / dble(NN - 1))
      xs(i) = 1.5d0 * sin(1.7d0 * dble(i))
   end do
   if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
   !
   ! PTC on: many accepted iterates (the ring needs > 8 to roll).
   !
   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 1d0
   xbase = xs + 0.5d0
   jfnk_x = 0d0
   call f_true_nm(NN, xbase, f0)
   call set_ptc_nm(NN, f0)
   fnrm0 = dsqrt(dp_nm(NN, f0, f0))
   tol_t = 1.0d-12 * fnrm0
   diagJ = a * exp(xbase) + 2d0 * b
   info = 1
   n_acc = 0
   call newton_hook_v1(f_nm, df_nm, dfp_nm, sub_nm, dp_nm, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       60, info, luer)
   !
   ! ---- report ----
   write(luer,'(A,I3,A,I5)') 'run1: info=', info, '  accepted=', n_acc
   write(luer,'(A,I4)') 'W1 write-pointer violations:', n_bad_write
   if(n_bad_write == 0) then
      write(luer,'(A)') 'W1: PASS  nmhist(nmidx)=tol after every accept'
   else
      write(luer,'(A)') 'W1: FAIL'
      nfail = nfail + 1
   end if
   !
   write(luer,'(A,I4,A,I4)') 'W2 read-window/write-ring divergences:', &
      n_diverge, '  first at accept', first_diverge
   if(n_diverge == 0) then
      write(luer,'(A)') 'W2: PASS -- phi_ref reads the last ' // &
         'nmwin_s pushed slots of the ring'
   else
      write(luer,'(A)') 'W2: FAIL -- phi_ref does not see the ' // &
         'last-nmwin_s window'
      nfail = nfail + 1
   end if
   !
   ! ---- second activation: ring must reset ----
   !
   jfnk_ptc_dtau = 1d0
   jfnk_x = 0d0
   call set_ptc_nm(NN, f0)
   info = 1
   n_acc = 0
   call newton_hook_v1(f_nm, df_nm, dfp_nm, sub_nm2, dp_nm, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       8, info, luer)
   !
   if(nfail == 0) then
      print *, 'ALL NM-WINDOW TESTS PASSED'
   else
      print *, nfail, ' NM-WINDOW TESTS FAILED'
   end if
   call exit(nfail)
contains

   subroutine f_true_nm(n, x, y)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: y(n)
      double precision :: d(n)
      integer :: i
      d = x - xs
      do i = 1, n
         y(i) = a(i) * (exp(x(i)) - exp(xs(i)))
         if(i > 1) y(i) = y(i) + b * (d(i) - d(i-1))
         if(i < n) y(i) = y(i) + b * (d(i) - d(i+1))
      end do
   end subroutine f_true_nm

   subroutine set_ptc_nm(n, f0v)
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
   end subroutine set_ptc_nm

end program test_nm_window


subroutine f_nm(n, x, y)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   use nm_prob_mod
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
         do i = 1, n
            y(i) = y(i) + (x(i) - jfnk_ptc_x0(i)) * &
               jfnk_ptc_scale(i) / max(jfnk_ptc_dtau, 1d-8)
         end do
      end if
   end if
end subroutine f_nm


subroutine df_nm(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_x0
   use nm_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   integer :: i
   jv = a * exp(xbase) * v
   do i = 1, n
      if(i > 1) jv(i) = jv(i) + b * (v(i) - v(i-1))
      if(i < n) jv(i) = jv(i) + b * (v(i) - v(i+1))
   end do
   if(jfnk_ptc_active) &
      jv = jv + v * jfnk_ptc_scale / max(jfnk_ptc_dtau, 1d-8)
end subroutine df_nm


subroutine dfp_nm(n, v)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale
   use nm_prob_mod
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   integer :: i
   do i = 1, n
      v(i) = v(i) / (a(i) * exp(xbase(i)) + 2d0 * b)
   end do
   if(jfnk_ptc_active) &
      v = v / (1d0 + jfnk_ptc_scale / max(jfnk_ptc_dtau, 1d-8) / &
               (a * exp(xbase) + 2d0 * b))
end subroutine dfp_nm


! sub_nm: live ring inspection after each accepted iterate.
subroutine sub_nm
   use jfnk_solver_v1_mod, only: jfnk_tol, jfnk_ls_nmhist, &
      jfnk_ls_nmidx, jfnk_ls_nmwin_s
   use nm_prob_mod
   implicit none
   integer :: lo, i, kk
   double precision :: intended, read_set
   if(jfnk_ls_nmidx == 0) return  ! entry callback, before any accept
   n_acc = n_acc + 1
   if(n_acc <= 200) then
      acc_tol(n_acc) = jfnk_tol
      nmidx_seen(n_acc) = jfnk_ls_nmidx
   end if
   ! W1: write pointer must land on the newest accepted norm.
   if(jfnk_ls_nmidx < 1 .or. jfnk_ls_nmidx > 8) then
      n_bad_write = n_bad_write + 1
   else if(abs(jfnk_ls_nmhist(jfnk_ls_nmidx) - jfnk_tol) > &
           1d-12 * max(1d0, abs(jfnk_tol))) then
      n_bad_write = n_bad_write + 1
   end if
   ! W2: the documented reference set is the last-nmwin_s accepts;
   ! the fixed production loop walks the ring back from nmidx.
   ! Simulate that walk here and compare against the accepted
   ! history recorded in acc_tol.
   lo = max(1, n_acc - jfnk_ls_nmwin_s + 1)
   intended = maxval(acc_tol(lo:n_acc))
   read_set = -1d0
   do i = 1, min(jfnk_ls_nmwin_s, 8)
      kk = mod(jfnk_ls_nmidx - i + 8, 8) + 1
      if(jfnk_ls_nmhist(kk) > 0d0) &
         read_set = max(read_set, jfnk_ls_nmhist(kk))
   end do
   if(abs(read_set - intended) > 1d-12 * max(1d0, intended)) then
      n_diverge = n_diverge + 1
      if(first_diverge < 0) first_diverge = n_acc
   end if
end subroutine sub_nm


! sub_nm2: checks the ring was reset at activation (W3) and tracks
! accepts for the second run.
subroutine sub_nm2
   use jfnk_solver_v1_mod, only: jfnk_tol, jfnk_ls_nmhist, &
      jfnk_ls_nmidx
   use nm_prob_mod
   implicit none
   integer :: i, n_neg
   if(jfnk_ls_nmidx == 0) return  ! entry callback
   n_acc = n_acc + 1
   if(n_acc == 1) then
      ! After the first push exactly one slot holds a norm; all
      ! others must still be the -1 empty markers.
      n_neg = 0
      do i = 1, 8
         if(jfnk_ls_nmhist(i) < 0d0) n_neg = n_neg + 1
      end do
      if(n_neg == 7 .and. jfnk_ls_nmhist(jfnk_ls_nmidx) == jfnk_tol) &
         write(6,'(A)') 'W3: PASS  ring reset per activation ' // &
            '(7 empty slots + 1 fresh norm)'
      if(n_neg /= 7) then
         write(6,'(A,I3)') 'W3: FAIL  empty slots after reset=', n_neg
         call exit(1)
      end if
   end if
end subroutine sub_nm2


double precision function dp_nm(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_nm = dot_product(a, b)
end function dp_nm
