! drv_lam_handoff -- CLI driver for LAM_HANDOFF_MOD (the extracted
! RD_LAMBDA -> full-iteration handoff predicate and stability-count
! step). Runs a fixed battery and prints one NAME=VALUE line per
! check; test_lam_handoff.py parses and asserts on them.
!
! The driver calls the production module routines -- not a copy of
! the intended logic -- so a regression in the gate, the counter,
! or the predicate fails here exactly as it would in OUTGEN.
!
program drv_lam_handoff
  use set_kind_module
  use lam_handoff_mod
  implicit none
!
  integer cnt
  logical proven
  integer i
!
! ---- LAM_FIELD_GATE(has_lagged, held, freeze, raw_immature,
!                     n_mature, rel_stable) -------------------------
!
! No lagged field (stock, plain SOB+LAM_ES=T): held/frozen flags
! alone decide; the lagged-field inputs are unjudged.
  write(*,'(A,L1)') 'G_nolag_clean=', &
    lam_field_gate(.false.,.false.,.false.,.false.,0,1)
  write(*,'(A,L1)') 'G_nolag_held=', &
    lam_field_gate(.false.,.true.,.false.,.false.,0,1)
  write(*,'(A,L1)') 'G_nolag_freeze=', &
    lam_field_gate(.false.,.false.,.true.,.false.,0,1)
  write(*,'(A,L1)') 'G_nolag_rawimm=', &
    lam_field_gate(.false.,.false.,.false.,.true.,0,1)
!
! Lagged-field routes: any term blocks; a mature field clears.
  write(*,'(A,L1)') 'G_lag_held=', &
    lam_field_gate(.true.,.true.,.false.,.false.,5,2)
  write(*,'(A,L1)') 'G_lag_freeze=', &
    lam_field_gate(.true.,.false.,.true.,.false.,5,2)
  write(*,'(A,L1)') 'G_lag_rawimm=', &
    lam_field_gate(.true.,.false.,.false.,.true.,5,2)
  write(*,'(A,L1)') 'G_lag_window0of2=', &
    lam_field_gate(.true.,.false.,.false.,.false.,0,2)
  write(*,'(A,L1)') 'G_lag_window1of2=', &
    lam_field_gate(.true.,.false.,.false.,.false.,1,2)
  write(*,'(A,L1)') 'G_lag_window2of2=', &
    lam_field_gate(.true.,.false.,.false.,.false.,2,2)
  write(*,'(A,L1)') 'G_lag_window0of1=', &
    lam_field_gate(.true.,.false.,.false.,.false.,0,1)
  write(*,'(A,L1)') 'G_lag_window1of1=', &
    lam_field_gate(.true.,.false.,.false.,.false.,1,1)
!
! ---- LAM_STABLE_STEP(force,lum,frozen,field_nrdy,maxch,val_exit,
!                      n_stable,rdlam_acc -> cnt,proven) ------------
!
! Earning run: 5 consecutive quiet calls (maxch=2 < val_exit=5)
! with n_stable=5 -- proven latches only on the fifth.
  cnt=0; proven=.false.
  do i=1,5
    call lam_stable_step(.false.,.false.,.false.,.false., &
      2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  end do
  write(*,'(A,I0,A,L1)') 'S_earn5=',cnt,' proven=',proven
!
  cnt=0; proven=.false.
  do i=1,4
    call lam_stable_step(.false.,.false.,.false.,.false., &
      2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  end do
  write(*,'(A,I0,A,L1)') 'S_earn4=',cnt,' proven=',proven
!
! One demand at/above val_exit after 4 earns resets the streak.
  call lam_stable_step(.false.,.false.,.false.,.false., &
    6.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_reset_demand=',cnt,' proven=',proven
!
  cnt=3; proven=.false.
  call lam_stable_step(.true.,.false.,.false.,.false., &
    2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_reset_force=',cnt,' proven=',proven
!
  cnt=3; proven=.false.
  call lam_stable_step(.false.,.true.,.false.,.false., &
    2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_reset_lum=',cnt,' proven=',proven
!
! Field-not-ready resets: the fix for the SOB-mom gap -- a quiet
! demand against an immature lagged field cannot earn or keep the
! stability count.
  cnt=3; proven=.false.
  call lam_stable_step(.false.,.false.,.false.,.true., &
    2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_reset_field=',cnt,' proven=',proven
!
! A proven latch is also revoked by a field-not-ready call.
  cnt=5; proven=.true.
  call lam_stable_step(.false.,.false.,.false.,.true., &
    2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_unlatch_field=',cnt,' proven=',proven
!
! Freeze-neutral semantics: a frozen call with the field gate clean
! earns nothing but does not reset (the earn arm requires
! .NOT. frozen). In production a frozen call always carries
! field_nrdy=T and resets through that leg instead.
  cnt=3; proven=.false.
  call lam_stable_step(.false.,.false.,.true.,.false., &
    2.0_LDP,5.0_LDP,5,.false.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_frozen_neutral=',cnt,' proven=',proven
!
! Resacc-accepted probe pins the count at n_stable while the field
! reports ready.
  cnt=0; proven=.false.
  call lam_stable_step(.false.,.false.,.false.,.false., &
    20.0_LDP,5.0_LDP,5,.true.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_pin_clean=',cnt,' proven=',proven
!
! The pin must not fire against a not-ready field (the old latch
! survived a rejection and re-pinned the count on the same call).
  cnt=0; proven=.false.
  call lam_stable_step(.false.,.false.,.false.,.true., &
    20.0_LDP,5.0_LDP,5,.true.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_pin_field_blocked=',cnt,' proven=',proven
!
  cnt=0; proven=.false.
  call lam_stable_step(.true.,.false.,.false.,.false., &
    20.0_LDP,5.0_LDP,5,.true.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_pin_force_blocked=',cnt,' proven=',proven
!
  cnt=0; proven=.false.
  call lam_stable_step(.false.,.true.,.false.,.false., &
    20.0_LDP,5.0_LDP,5,.true.,cnt,proven)
  write(*,'(A,I0,A,L1)') 'S_pin_lum_blocked=',cnt,' proven=',proven
!
! ---- LAM_HANDOFF_OK(do_lam_auto,rd_lam,maxch,jfnk_unconv,
!                     nonmono,stable_ok,field_nrdy) -----------------
!
! Stock path (nonmono=F): MAXCH_CTRL < 50 alone decides; stable_ok
! and field_nrdy are not consulted.
  write(*,'(A,L1)') 'H_stock_clean=', &
    lam_handoff_ok(.true.,.true.,49.9_LDP,.false.,.false., &
                   .false.,.false.)
  write(*,'(A,L1)') 'H_stock_at50=', &
    lam_handoff_ok(.true.,.true.,50.0_LDP,.false.,.false., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_stock_jfnk=', &
    lam_handoff_ok(.true.,.true.,10.0_LDP,.true.,.false., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_stock_ignores=', &
    lam_handoff_ok(.true.,.true.,10.0_LDP,.false.,.false., &
                   .false.,.true.)
  write(*,'(A,L1)') 'H_stock_nordlam=', &
    lam_handoff_ok(.true.,.false.,10.0_LDP,.false.,.false., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_stock_noauto=', &
    lam_handoff_ok(.false.,.true.,10.0_LDP,.false.,.false., &
                   .true.,.false.)
!
! Non-monotonic path: stable_ok AND field ready are both required
! on top of the stock terms.
  write(*,'(A,L1)') 'H_nm_clean=', &
    lam_handoff_ok(.true.,.true.,20.0_LDP,.false.,.true., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_nm_unstable=', &
    lam_handoff_ok(.true.,.true.,20.0_LDP,.false.,.true., &
                   .false.,.false.)
  write(*,'(A,L1)') 'H_nm_field_nrdy=', &
    lam_handoff_ok(.true.,.true.,20.0_LDP,.false.,.true., &
                   .true.,.true.)
  write(*,'(A,L1)') 'H_nm_at50=', &
    lam_handoff_ok(.true.,.true.,50.0_LDP,.false.,.true., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_nm_jfnk=', &
    lam_handoff_ok(.true.,.true.,20.0_LDP,.true.,.true., &
                   .true.,.false.)
  write(*,'(A,L1)') 'H_nm_quiet_dem=', &
    lam_handoff_ok(.true.,.true.,1.0_LDP,.false.,.true., &
                   .true.,.false.)
!
  write(*,'(A)') 'DONE'
end program drv_lam_handoff
