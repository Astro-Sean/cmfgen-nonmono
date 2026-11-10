! drv_mom_aa_state -- CLI driver for the audit-D2 GET_MOM_AA_STATE
! census getter (plane/mom_j_cmf_nonmono_v1.f). Sets the production
! module state of MOD_MOM_J_NONMONO_V1 for each scenario and calls
! the real getter; the Python test asserts on the printed state.
! The point under test is the state-vs-activity separation:
! enabled flags and the alloc-fail latch are configuration state
! and must be reported live, while the AA counters are per-refresh
! activity and must read zero whenever MOM_RF_RAN is false (the
! refresh counter reset sits after the self-guard early returns,
! so a guarded-out or never-ran call would otherwise leak the
! previous refresh's tallies).
!
!   argv1:
!     'a' refresh ran, inner+outer AA enabled, counters nonzero
!     'b' refresh ran, inner AA disabled
!     'c' refresh ran, outer AA disabled
!     'd' refresh did NOT run, stale counters from a prior refresh
!     'e' formal-J path: refresh never ran, counters zero
!     'f' refresh did NOT run, outer alloc-fail latch set, stale
!         counters (AOF must stay distinguishable from 'disabled')
!     'g' refresh ran, both enabled, zero applications (enabled
!         state must not be inferred from counters)
!
program drv_mom_aa_state
  use set_kind_module
  use mod_mom_j_nonmono_v1
  implicit none
!
  character(len=4) mode
  logical aa_on,aa_out_on,afail,freqm,rf_ran
  integer*8 aapp,askip,asng,arst,aapp_o,arej_o,asng_o
!
  call getarg(1,mode)
!
! Baseline: every case writes the full module state it depends on,
! so scenario order cannot leak through the SAVE'd module vars.
!
  mom_rf_aa    = .true.
  mom_aa_out   = .true.
  mom_aa_out_fail = .false.
  mom_freq_match  = .false.
  mom_rf_ran   = .true.
  mom_aa_app   = 0_8
  mom_aa_skip  = 0_8
  mom_aa_sng   = 0_8
  mom_aa_rst   = 0_8
  mom_aa_out_app = 0_8
  mom_aa_out_rej = 0_8
  mom_aa_out_sng = 0_8
!
  select case (mode)
  case ('a')
    mom_aa_app=3_8; mom_aa_skip=1_8; mom_aa_rst=2_8
    mom_aa_out_app=2_8; mom_aa_out_rej=1_8
  case ('b')
    mom_rf_aa=.false.
    mom_aa_app=4_8
    mom_aa_out_app=1_8
  case ('c')
    mom_aa_out=.false.
    mom_aa_app=5_8
  case ('d')
    mom_rf_ran=.false.
    mom_aa_app=7_8; mom_aa_skip=2_8; mom_aa_sng=1_8; mom_aa_rst=3_8
    mom_aa_out_app=4_8; mom_aa_out_rej=2_8; mom_aa_out_sng=1_8
  case ('e')
    mom_rf_ran=.false.
  case ('f')
    mom_rf_ran=.false.
    mom_aa_out_fail=.true.
    mom_aa_app=6_8
    mom_aa_out_app=2_8
  case ('g')
    continue
  case default
    write(*,'(A)') 'BAD_MODE'
    stop 1
  end select
!
  call get_mom_aa_state(aa_on,aapp,askip,asng,arst, &
       aa_out_on,aapp_o,arej_o,asng_o,afail,freqm,rf_ran)
!
  write(*,'(A,L1)')      'AA_ON='    ,aa_on
  write(*,'(A,L1)')      'AA_OUT_ON=',aa_out_on
  write(*,'(A,L1)')      'AFAIL='    ,afail
  write(*,'(A,L1)')      'FREQM='    ,freqm
  write(*,'(A,L1)')      'RF_RAN='   ,rf_ran
  write(*,'(A,I0)')      'AAPP='     ,aapp
  write(*,'(A,I0)')      'ASKIP='    ,askip
  write(*,'(A,I0)')      'ASNG='     ,asng
  write(*,'(A,I0)')      'ARST='     ,arst
  write(*,'(A,I0)')      'AAPP_O='   ,aapp_o
  write(*,'(A,I0)')      'AREJ_O='   ,arej_o
  write(*,'(A,I0)')      'ASNG_O='   ,asng_o
!
end program drv_mom_aa_state
