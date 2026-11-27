!----------------------------------------------------------------------
! TEST_NM_ACCEPT: white-box check of WP-1 (jfnk.md) - a line-search
! step that passes the GLL reference test with |F| still rising must
! be commitable in trust-region mode under [JFNK_NM_ACCEPT], instead
! of being discarded as nm_increase_rejected.
!
! Problem: F(x) = x * r**p, p = -0.9, r = ||x|| (root at x = 0).
! The map is homogeneous, so the exact Newton step is s = x/(1+p) =
! 10x along +x and |F| scales as r**(p+1) = r**0.1. Along the step,
!   |F(x - a s)| / |F(x)| = |10a - 1|**0.1:
!     a = 1     -> 1.246  (+25%, ared < 0, LS fires)
!     a = 0.5   -> 1.149  (+15%, fails the 1.21*phi_0 cap)
!     a = 0.25  -> 1.041  (+4.1%, inside the cap once the nm window
!                          holds an older, larger accepted norm)
!     a = 0.125 -> 0.87   (decrease; the first two passes land here
!                          and seed the window)
! J and J^{-1} are analytic (Sherman-Morrison), so the GMRES
! direction is the exact Newton direction and the LS ladder, not
! solver noise, decides the outcome.
!
! Runs the same activation twice on the production .o:
!   run A: jfnk_do_nm_accept = F (default) -> nm_increase_rejected rows,
!          no nm_armijo_increase commits (the "before" behaviour);
!   run B: jfnk_do_nm_accept = T, credit cap 2 -> GLL-passed increases
!          commit, the credit bound holds, and rejects count only
!          genuine rejections.
!----------------------------------------------------------------------
module nma_prob_mod
   implicit none
   integer, parameter :: NN = 3
   double precision, parameter :: PEXP = -0.9d0
   integer :: n_acc = 0          ! accepted-iterate callbacks seen
end module nma_prob_mod


program test_nm_accept
   use jfnk_solver_v1_mod
   use nma_prob_mod
   implicit none
   external :: f_nma, df_nma, dfp_nma, sub_nma
   double precision, external :: dp_nma
   integer :: info, nfail, luer, ios, run_max
   integer :: nrej_a, ninc_a, ninc_b, nrej_b, lsacc_b, tr_b0, tr_b1
   integer :: nrej_rows_b
   double precision :: x0(NN), tol_t, fnrm0
   character(len=512) :: line
   character(len=32) :: tok(32), lab
   !
   luer = 6
   nfail = 0
   x0 = (/ 0.06d0, 0.08d0, 0.10d0 /)
   if(.not. allocated(jfnk_x)) allocate(jfnk_x(NN), jfnk_fx(NN))
   !
   jfnk_ptc_active = .false.
   jfnk_do_line_search = .true.
   jfnk_ls_nonmono = .true.
   jfnk_ls_nmwin_s = 4
   jfnk_do_cauchy = .false.     ! keep the reject path clean
   jfnk_x = x0
   call f_nma(NN, jfnk_x, jfnk_fx)
   fnrm0 = dsqrt(dp_nma(NN, jfnk_fx, jfnk_fx))
   tol_t = 1.0d-14 * fnrm0
   !
   ! ---- run A: default (JFNK_NM_ACCEPT off) ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   jfnk_do_nm_accept = .false.
   n_acc = 0
   info = 1
   call newton_hook_v1(f_nma, df_nma, dfp_nma, sub_nma, dp_nma, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   write(luer,'(A,I3,A,I5,A,I5)') 'runA: info=', info, &
      '  accepted=', n_acc, '  n_ls_used=', jfnk_n_ls_used
   call scan_ledger(ninc_a, nrej_a, lsacc_b, run_max)
   write(luer,'(A,I4,A,I4)') 'runA: nm_armijo_increase rows=', &
      ninc_a, '  nm_increase_rejected rows=', nrej_a
   if(ninc_a == 0 .and. nrej_a >= 1) then
      write(luer,'(A)') 'A1: PASS  default path discards the ' // &
         'GLL-passed increase (nm_increase_rejected)'
   else
      write(luer,'(A)') 'A1: FAIL'
      nfail = nfail + 1
   end if
   !
   ! ---- run B: JFNK_NM_ACCEPT on, credit cap 2 ----
   !
   call execute_command_line('rm -f JFNK_LEDGER')
   jfnk_do_nm_accept = .true.
   jfnk_nm_credit_cap = 2
   jfnk_n_ls_used = 0
   tr_b0 = jfnk_n_tr_rejected
   jfnk_x = x0
   n_acc = 0
   info = 1
   call newton_hook_v1(f_nma, df_nma, dfp_nma, sub_nma, dp_nma, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       40, info, luer)
   tr_b1 = jfnk_n_tr_rejected
   write(luer,'(A,I3,A,I5,A,I5)') 'runB: info=', info, &
      '  accepted=', n_acc, '  n_ls_used=', jfnk_n_ls_used
   call scan_ledger(ninc_b, nrej_b, lsacc_b, run_max)
   nrej_rows_b = nrej_b
   write(luer,'(A,I4,A,I4,A,I4)') 'runB: nm_armijo_increase rows=', &
      ninc_b, '  nm_increase_rejected rows=', nrej_b, &
      '  ls_accepted rows=', lsacc_b
   !
   ! B1: at least one GLL-passed increase commits.
   if(ninc_b >= 1) then
      write(luer,'(A)') 'B1: PASS  nm_armijo_increase commits ' // &
         'under JFNK_NM_ACCEPT'
   else
      write(luer,'(A)') 'B1: FAIL  no nm_armijo_increase row'
      nfail = nfail + 1
   end if
   ! B2: every nm-armijo-increase row came through the line search.
   if(lsacc_b >= ninc_b) then
      write(luer,'(A)') 'B2: PASS  increase accepts are LS trials'
   else
      write(luer,'(A)') 'B2: FAIL  non-LS increase accept'
      nfail = nfail + 1
   end if
   ! B3: credit bound - no run of consecutive increase accepts may
   ! exceed the cap.
   if(run_max <= jfnk_nm_credit_cap) then
      write(luer,'(A,I2)') 'B3: PASS  consecutive increase ' // &
         'accepts bounded, max run=', run_max
   else
      write(luer,'(A,I2)') 'B3: FAIL  run of increase accepts=', &
         run_max
      nfail = nfail + 1
   end if
   ! B4: rejects counted only genuine rejections: the tr_rejected
   ! delta must equal the nm_increase_rejected row count.
   if(tr_b1 - tr_b0 == nrej_rows_b) then
      write(luer,'(A)') 'B4: PASS  n_tr_rejected counts only ' // &
         'genuine rejections'
   else
      write(luer,'(A,2I5)') 'B4: FAIL  tr_rejected delta / ' // &
         'reject rows=', tr_b1 - tr_b0, nrej_rows_b
      nfail = nfail + 1
   end if
   ! B5: the iterate advanced - accepts reached sub().
   if(n_acc >= 1) then
      write(luer,'(A,I4)') 'B5: PASS  accepted iterates ' // &
         'reached sub(), n_acc=', n_acc
   else
      write(luer,'(A)') 'B5: FAIL  iterate never advanced'
      nfail = nfail + 1
   end if
   !
   if(nfail == 0) then
      print *, 'ALL NM-ACCEPT TESTS PASSED'
   else
      print *, nfail, ' NM-ACCEPT TESTS FAILED'
   end if
   call exit(nfail)
contains

   ! Scan JFNK_LEDGER for this run: count nm_armijo_increase rows,
   ! nm_increase_rejected rows, ls_accepted=1 rows, and the longest
   ! run of consecutive nm_armijo_increase rows.
   subroutine scan_ledger(n_inc, n_rej, n_lsacc, run_max_out)
      implicit none
      integer, intent(out) :: n_inc, n_rej, n_lsacc, run_max_out
      integer :: uu, rl
      character(len=64) :: word
      n_inc = 0
      n_rej = 0
      n_lsacc = 0
      run_max_out = 0
      rl = 0
      open(newunit=uu, file='JFNK_LEDGER', status='old', &
           iostat=ios)
      if(ios /= 0) return
      read(uu,'(A)',iostat=ios) line          ! header
      do
         read(uu,'(A)',iostat=ios) line
         if(ios /= 0) exit
         ! cols: 17 ls_accepted, 24 taken_branch
         call field_at(line, 17, word)
         if(index(adjustl(word), '1') == 1) n_lsacc = n_lsacc + 1
         call field_at(line, 24, lab)
         if(trim(lab) == 'nm_armijo_increase') then
            n_inc = n_inc + 1
            rl = rl + 1
            if(rl > run_max_out) run_max_out = rl
         else
            rl = 0
         end if
         if(trim(lab) == 'nm_increase_rejected') n_rej = n_rej + 1
      end do
      close(uu)
   end subroutine scan_ledger

   subroutine field_at(str, want, out)
      implicit none
      character(len=*), intent(in) :: str
      integer, intent(in) :: want
      character(len=*), intent(out) :: out
      integer :: ii, kk, st, en
      out = ' '
      kk = 1
      st = 1
      do ii = 1, len_trim(str)
         if(str(ii:ii) == ',') then
            en = ii - 1
            if(kk == want) then
               out = adjustl(str(st:en))
               return
            end if
            kk = kk + 1
            st = ii + 1
         end if
      end do
      if(kk == want) out = adjustl(str(st:len_trim(str)))
   end subroutine field_at

end program test_nm_accept


subroutine f_nma(n, x, y)
   use nma_prob_mod, only: PEXP
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: x(n)
   double precision, intent(out) :: y(n)
   double precision :: r
   r = max(dsqrt(dot_product(x, x)), 1d-30)
   y = x * r**PEXP
end subroutine f_nma


subroutine df_nma(n, v, jv)
   ! J v = r^p v + p r^{p-2} (x . v) x  at the current jfnk_x.
   use jfnk_solver_v1_mod, only: jfnk_x
   use nma_prob_mod, only: PEXP
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: v(n)
   double precision, intent(out) :: jv(n)
   double precision :: r, xv
   r = max(dsqrt(dot_product(jfnk_x, jfnk_x)), 1d-30)
   xv = dot_product(jfnk_x, v)
   jv = r**PEXP * v + PEXP * r**(PEXP - 2d0) * xv * jfnk_x
end subroutine df_nma


subroutine dfp_nma(n, v)
   ! Exact inverse via Sherman-Morrison:
   ! J^{-1} = r^{-p} (I - p/(1+p) u u^T), u = x/r.
   use jfnk_solver_v1_mod, only: jfnk_x
   use nma_prob_mod, only: PEXP
   implicit none
   integer, intent(in) :: n
   double precision, intent(inout) :: v(n)
   double precision :: r, xv
   r = max(dsqrt(dot_product(jfnk_x, jfnk_x)), 1d-30)
   xv = dot_product(jfnk_x, v) / (r * r)
   v = r**(-PEXP) * (v - (PEXP / (1d0 + PEXP)) * xv * jfnk_x)
end subroutine dfp_nma


subroutine sub_nma
   use jfnk_solver_v1_mod, only: jfnk_ls_nmidx
   use nma_prob_mod, only: n_acc
   implicit none
   if(jfnk_ls_nmidx == 0) return   ! entry callback, before any accept
   n_acc = n_acc + 1
end subroutine sub_nma


double precision function dp_nma(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_nma = dot_product(a, b)
end function dp_nma
