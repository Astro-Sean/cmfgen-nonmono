!----------------------------------------------------------------------
! TEST_STEP_STATUS: white-box check of the WP-0 step ledger in
! newton_hook_v1 (production jfnk_solver_v1.o). Drives a solve on the
! same stiff synthetic family as test_nm_window, then parses
! JFNK_LEDGER (one CSV row per decided Newton pass + a terminal row).
!
! Assertions (WP-0.1/0.3/0.4):
!   T1  header matches the documented column list.
!   T2  every taken_branch is a known label; 'none' never written.
!   T3  nits nondecreasing; the k-th accept-labelled row carries
!       nits = k-1; accept rows total the sub()-counted accepts; a
!       'converged' terminal row has nits = #accepts.
!   T4  ls_accepted rows have 0 < alpha <= 1; cauchy_used rows are
!       also ls_accepted.
!   T5  eval accounting: the driver's f increments the stub's
!       JFNK_N_RESIDUAL_EVALS, so the last row's n_res_evals must
!       equal the driver-side call count AND the solver-side
!       jfnk_ev_entry + trial + ls + cauchy sum.
!   T6  a terminal label appears only in the last row.
!----------------------------------------------------------------------
module ss_prob_mod
   implicit none
   integer, parameter :: NN = 60
   double precision :: a(NN), xs(NN), xbase(NN)
   double precision :: b
   integer :: n_acc = 0
   integer :: n_f_calls = 0
end module ss_prob_mod


program test_step_status
   use jfnk_solver_v1_mod
   use jfnk_cmfgen_int_mod, only: JFNK_N_RESIDUAL_EVALS
   use ss_prob_mod
   implicit none
   external :: f_ss, df_ss, dfp_ss, sub_ss
   double precision, external :: dp_ss
   integer :: i, info, nfail, luer, lu, ios
   double precision :: x0(NN), f0(NN), diagJ(NN)
   double precision :: fnrm0, tol_t
   character(len=2000) :: line, hdr
   character(len=64) :: tok(32)
   integer :: ntok, nrows, iacc, nacc_rows, nterm, bad_lab
   integer :: nits_prev, nits_row, lseq
   integer :: v_lsacc, v_cauchy, v_nev
   double precision :: v_alpha
   logical :: last_was_term, any_ls
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
   ! Fresh ledger for this run: delete rather than truncate so the
   ! solver writes its header row.
   !
   open(newunit=lu, file='JFNK_LEDGER', status='old', iostat=ios)
   if(ios == 0) close(lu, status='delete')
   !
   jfnk_ptc_active = .true.
   jfnk_ptc_dtau = 1d0
   xbase = xs + 0.5d0
   jfnk_x = 0d0
   JFNK_N_RESIDUAL_EVALS = 0
   n_f_calls = 0
   call f_true_ss(NN, xbase, f0)
   call set_ptc_ss(NN, f0)
   fnrm0 = dsqrt(dp_ss(NN, f0, f0))
   tol_t = 1.0d-12 * fnrm0
   diagJ = a * exp(xbase) + 2d0 * b
   info = 1
   n_acc = 0
   call newton_hook_v1(f_ss, df_ss, dfp_ss, sub_ss, dp_ss, &
                       20, NN, 1.0d-2, tol_t, -1d0, 1.0d-12, 1.0d99, &
                       60, info, luer)
   write(luer,'(A,I3,A,I5,A,I8)') 'run: info=', info, &
      '  accepted=', n_acc, '  f_calls=', n_f_calls
   !
   ! ---- parse JFNK_LEDGER ----
   !
   open(newunit=lu, file='JFNK_LEDGER', status='old', iostat=ios)
   if(ios /= 0) then
      write(luer,'(A)') 'T1: FAIL  JFNK_LEDGER not written'
      call exit(1)
   end if
   read(lu,'(A)') hdr
   if(hdr == 'nits,fw,eta,gres,ref_norm,gmres_its,gmres_restarts,'// &
      'gmres_info,hess_res,true_res,true_res_ratio,dir_deriv_F,'// &
      'dir_deriv_G,cos_angle,alpha_ls,ls_steps,ls_accepted,'// &
      'nm_phi_ref,nm_phi_alpha,cauchy_used,del_before,del_after,'// &
      'n_reject,taken_branch,exit_reason,n_res_evals,wall_s') then
      write(luer,'(A)') 'T1: PASS  header matches column list'
   else
      write(luer,'(A)') 'T1: FAIL  unexpected header'
      nfail = nfail + 1
   end if
   !
   nrows = 0
   nacc_rows = 0
   nits_prev = -1
   bad_lab = 0
   nterm = 0
   last_was_term = .false.
   any_ls = .false.
   v_nev = -1
   do
      read(lu,'(A)',iostat=ios) line
      if(ios /= 0) exit
      if(len_trim(line) == 0) cycle
      if(last_was_term) then
         ! a data row after a terminal label: order violation
         nterm = nterm + 1
      end if
      nrows = nrows + 1
      call csv_split(line, tok, ntok)
      if(ntok < 27) then
         write(luer,'(A,I3,A)') 'T?: row with <27 fields: ', ntok, &
            ' -- treated as bad label'
         bad_lab = bad_lab + 1
         cycle
      end if
      read(tok(1),*) nits_row
      read(tok(15),*) v_alpha
      read(tok(17),*) v_lsacc
      read(tok(20),*) v_cauchy
      read(tok(26),*) v_nev
      ! T3: monotonic nits
      if(nits_row < nits_prev) then
         write(luer,'(A,I6,A,I6)') 'T3: nits regression ', &
            nits_prev, ' -> ', nits_row
         nfail = nfail + 1
      end if
      nits_prev = nits_row
      ! T2/T3: label classes
      select case(trim(tok(24)))
      case('nm_armijo_decrease','ls_accept_noTR','full_newton_noTR', &
           'accept_previous','accept_full_newton','accept_at_maxdel')
         nacc_rows = nacc_rows + 1
         if(nits_row /= nacc_rows - 1) then
            write(luer,'(A,A,A,I6)') 'T3: accept row ', &
               trim(tok(24)), ' carries nits=', nits_row
            nfail = nfail + 1
         end if
      case('reject_runaway_step','reject_nan_inf', &
           'reject_pred_negligible','reject_no_decrease', &
           'retry_smaller_tr','retry_larger_tr','nm_increase_rejected')
         continue
      case('converged','max_newton_iterations', &
           'converged_at_entry','trust_region_too_small', &
           'precond_clip_forcing_norm','krylov_jvp_or_pc_failure', &
           'gmres_stagnation_abort','jvp_failure_pred_eval', &
           'cauchy_alloc_failure','too_many_rejected_steps', &
           'alloc_fail_workspace','alloc_fail_v')
         last_was_term = .true.
         if(trim(tok(24)) == 'converged' .and. &
            nits_row /= n_acc) then
            write(luer,'(A)') 'T3: converged row nits /= #accepts'
            nfail = nfail + 1
         end if
      case default
         bad_lab = bad_lab + 1
         write(luer,'(A,A)') 'T2: unknown taken_branch=', trim(tok(24))
      end select
      ! T4: alpha/lsacc/cauchy consistency
      if(v_lsacc == 1) then
         any_ls = .true.
         if(v_alpha <= 0d0 .or. v_alpha > 1d0) then
            write(luer,'(A,ES12.4)') 'T4: ls_accepted alpha=', v_alpha
            nfail = nfail + 1
         end if
      end if
      if(v_cauchy == 1 .and. v_lsacc /= 1) then
         write(luer,'(A)') 'T4: cauchy_used without ls_accepted'
         nfail = nfail + 1
      end if
   end do
   close(lu)
   !
   if(bad_lab == 0) then
      write(luer,'(A,I5)') 'T2: PASS  all taken_branch labels known, rows=', &
         nrows
   else
      write(luer,'(A,I5)') 'T2: FAIL  unknown labels=', bad_lab
      nfail = nfail + 1
   end if
   if(nacc_rows == n_acc) then
      write(luer,'(A,I5)') 'T3: PASS  accept rows == sub() accepts=', &
         nacc_rows
   else
      write(luer,'(A,I5,A,I5)') 'T3: FAIL  accept rows=', nacc_rows, &
         ' sub() accepts=', n_acc
      nfail = nfail + 1
   end if
   if(nterm == 0) then
      write(luer,'(A)') 'T6: PASS  terminal labels only in last row'
   else
      write(luer,'(A,I5)') 'T6: FAIL  rows after terminal=', nterm
      nfail = nfail + 1
   end if
   ! T5: eval accounting -- last data row's n_res_evals must equal
   ! the driver's own f-call count and the solver's ev_* sum.
   iacc = jfnk_ev_entry + jfnk_ev_trial + jfnk_ev_ls + jfnk_ev_cauchy
   if(v_nev == n_f_calls .and. iacc == n_f_calls) then
      write(luer,'(A,I8)') 'T5: PASS  eval accounting reconciles at', &
         n_f_calls
   else
      write(luer,'(A,I8,A,I8,A,I8)') 'T5: FAIL  ledger=', v_nev, &
         ' solver ev sum=', iacc, ' driver f calls=', n_f_calls
      nfail = nfail + 1
   end if
   if(.not. any_ls) &
      write(luer,'(A)') 'note: no ls_accepted rows this run ' // &
         '(T4 vacuous)'
   if(jfnk_ev_ls_wasted > jfnk_ev_ls) then
      write(luer,'(A)') 'T5b: FAIL  wasted > ls evals'
      nfail = nfail + 1
   end if
   !
   if(nfail == 0) then
      print *, 'ALL STEP-STATUS TESTS PASSED'
   else
      print *, nfail, ' STEP-STATUS TESTS FAILED'
   end if
   call exit(nfail)
contains

   subroutine csv_split(lin, tok, ntok)
      implicit none
      character(len=*), intent(in) :: lin
      character(len=64), intent(out) :: tok(32)
      integer, intent(out) :: ntok
      integer :: p0, p1, k
      ntok = 0
      p0 = 1
      do k = 1, 32
         p1 = index(lin(p0:), ',')
         if(p1 == 0) then
            ntok = ntok + 1
            tok(ntok) = lin(p0:)
            exit
         end if
         ntok = ntok + 1
         tok(ntok) = lin(p0:p0 + p1 - 2)
         p0 = p0 + p1
      end do
   end subroutine csv_split

   subroutine f_true_ss(n, x, y)
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
   end subroutine f_true_ss

   subroutine set_ptc_ss(n, f0v)
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
   end subroutine set_ptc_ss

end program test_step_status


subroutine f_ss(n, x, y)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_scale_ready, jfnk_ptc_x0, &
      jfnk_ptc_anchor_ready
   use jfnk_cmfgen_int_mod, only: JFNK_N_RESIDUAL_EVALS
   use ss_prob_mod
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
   n_f_calls = n_f_calls + 1
   JFNK_N_RESIDUAL_EVALS = JFNK_N_RESIDUAL_EVALS + 1
end subroutine f_ss


subroutine df_ss(n, v, jv)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale, jfnk_ptc_x0
   use ss_prob_mod
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
end subroutine df_ss


subroutine dfp_ss(n, v)
   use jfnk_solver_v1_mod, only: jfnk_ptc_active, jfnk_ptc_dtau, &
      jfnk_ptc_scale
   use ss_prob_mod
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
end subroutine dfp_ss


subroutine sub_ss
   use jfnk_solver_v1_mod, only: jfnk_ls_nmidx
   use ss_prob_mod
   implicit none
   if(jfnk_ls_nmidx == 0) return
   n_acc = n_acc + 1
end subroutine sub_ss


double precision function dp_ss(n, a, b)
   implicit none
   integer, intent(in) :: n
   double precision, intent(in) :: a(n), b(n)
   dp_ss = dot_product(a, b)
end function dp_ss
