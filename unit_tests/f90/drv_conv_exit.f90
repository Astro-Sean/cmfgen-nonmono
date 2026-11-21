! drv_conv_exit -- CLI driver for the nine-gate converged-exit
! replica CONV_EXIT_EVAL (LAM_HANDOFF_MOD) and the binding-margin
! classifier CONV_MARGIN_EVAL. Prints one NAME=VALUE line per
! check; test_conv_exit.py parses and asserts on them.
!
! The driver calls the production module routines -- not a copy of
! the intended logic -- so a drifted replica conjunct, label or
! margin classification fails here exactly as it would in OUTGEN.
!
program drv_conv_exit
  use set_kind_module
  use lam_handoff_mod
  implicit none
!
  logical lam_ok,all_met
  integer npass
  character(len=160) fails
  integer i,nbad,npop
  logical bits(9)
  integer mc_lam,last_ng
  real(kind=ldp) ex_maxch
  real(kind=ldp) act_m,rel_m,trb_m
  character(len=8) act_c,rel_c,rel_w
  real(kind=ldp), parameter :: EPS=1.0_LDP
!
! ---- exhaustive truth table over the nine gates ------------------
! Bit meaning (1 = gate passes):
!   1 lambda-state ok, 2 not jfnk-unconverged, 3 ex_maxch<eps,
!   4 dem_ok, 5 not hydro-revision, 6 not pops-held,
!   7 not resacc-pending, 8 not on-ng-iteration,
!   9 not field-not-ready.
! lam_ok is driven through last_lam=0,rd_lam=F,do_lam_auto=F and
! main_counter = 1 (pass) or 0 (fail); the on-ng gate drives
! last_ng = -1 (pass) or = main_counter (fail).
!
  nbad=0
  do i=0,511
    bits(1)=btest(i,0); bits(2)=btest(i,1); bits(3)=btest(i,2)
    bits(4)=btest(i,3); bits(5)=btest(i,4); bits(6)=btest(i,5)
    bits(7)=btest(i,6); bits(8)=btest(i,7); bits(9)=btest(i,8)
    mc_lam=merge(1,0,bits(1))
    last_ng=merge(-1,mc_lam,bits(8))
    ex_maxch=merge(0.5_LDP,2.0_LDP,bits(3))
    call conv_exit_eval(.false.,.false.,0,mc_lam, &
      .not.bits(2),ex_maxch,EPS,bits(4),.not.bits(5), &
      .not.bits(6),.not.bits(7),last_ng,.not.bits(9), &
      lam_ok,all_met,npass,fails)
    npop=count(bits)
    if (all_met .neqv. (npop .eq. 9) .or. npass .ne. npop) &
      nbad=nbad+1
    if (lam_ok .neqv. bits(1)) nbad=nbad+1
  end do
  write(*,'(A,I0)') 'EXH_nbad=',nbad
!
! ---- named cases -------------------------------------------------
!
! Case A: the eight pre-fix gates pass, field not ready -> HOLD.
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.true.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,I0,A,A)') 'A_holdfield=',all_met, &
    ' npass=',npass,' fails=',trim(fails)
!
! Case B: all nine pass -> MET.
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,I0,A,A)') 'B_allpass=',all_met, &
    ' npass=',npass,' fails=',trim(fails)
!
! Case C: each single failure in turn.
  call conv_exit_eval(.false.,.false.,0,0, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_lamstate=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .true.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_jfnk=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,2.0_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_exm=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.false.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_dem=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.true.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_hydro=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.false.,.true.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_held=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.true., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_pend=',all_met,' fails=',trim(fails)
!
  call conv_exit_eval(.false.,.false.,0,1, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,A)') 'C_onng=',all_met,' fails=',trim(fails)
!
! Case D: several gates fail together -> all labels present.
  call conv_exit_eval(.false.,.false.,0,0, &
    .true.,2.0_LDP,EPS,.false.,.false.,.true.,.true., &
    0,.true.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,I0,A,A)') 'D_multi=',all_met, &
    ' npass=',npass,' fails=',trim(fails)
!
! Forced-lambda leg: rd_lam=T,do_lam_auto=F passes the lambda gate
! even with last_lam=main_counter.
  call conv_exit_eval(.true.,.false.,0,0, &
    .false.,0.5_LDP,EPS,.true.,.false.,.false.,.false., &
    -1,.false.,lam_ok,all_met,npass,fails)
  write(*,'(A,L1,A,L1)') 'E_forcedlam=',all_met,' lam_ok=',lam_ok
!
! ---- CONV_MARGIN_EVAL ---------------------------------------------
! eps=0.05 -> small=1.0, sat_raw=2.0e3, ex_maxch=0.03 (passing).
! Print per case: act margin/channel, released margin/channel/why,
! stage-B margin.
!
! SAT-released exm must not win the binding slot: the run-class
! worst-margin defect.
  call conv_margin_eval(0.05_LDP,0.03_LDP,1.0e15_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_sat=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Unsaturated material exm binds and beats ex_maxch.
  call conv_margin_eval(0.05_LDP,0.03_LDP,500.0_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_bind=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Small exempted demand is released, never binding.
  call conv_margin_eval(0.05_LDP,5.0_LDP,0.4_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_small=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Stage-B needed releases the channel but is not a SAT release.
  call conv_margin_eval(0.05_LDP,0.03_LDP,1500.0_LDP,0.0_LDP, &
    2.0e3_LDP,.true.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_trcneed=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Stage-B complete (release) releases the channel too.
  call conv_margin_eval(0.05_LDP,0.03_LDP,1500.0_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.true.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_trcrel=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! SAT outranks a simultaneous trc-need in the reason string.
  call conv_margin_eval(0.05_LDP,0.03_LDP,5.0e3_LDP,0.0_LDP, &
    2.0e3_LDP,.true.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_sat_trc=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! exr binding: unsaturated relative-floor demand wins the active
! slot over a small-margin ex_maxch.
  call conv_margin_eval(0.05_LDP,0.03_LDP,0.0_LDP,800.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_exr_bind=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! exr SAT-released: the relative-floor channel reports on the
! released row.
  call conv_margin_eval(0.05_LDP,0.03_LDP,0.4_LDP,3.0e3_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_exr_sat=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! BA_EX_RELFRAC disabled -> exr=0 contributes nothing anywhere.
  call conv_margin_eval(0.05_LDP,0.03_LDP,0.4_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_exr_off=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! BA_DEM_PFLOOR disabled -> exm=0 contributes nothing anywhere.
  call conv_margin_eval(0.05_LDP,0.03_LDP,0.0_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,0.0_LDP,0.0_LDP,0, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_pfl_off=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Stage-B demand armed: trb margin reads vs absmax.
  call conv_margin_eval(0.05_LDP,0.03_LDP,0.4_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,100.0_LDP,10.0_LDP,3, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_trabs=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
! Stage-B armed is not a binding gate: a live trb demand cannot
! take the active slot from the exempted-demand channels.
  call conv_margin_eval(0.05_LDP,50.0_LDP,0.4_LDP,0.0_LDP, &
    2.0e3_LDP,.false.,.false.,1.0e6_LDP,10.0_LDP,3, &
    act_m,act_c,rel_m,rel_c,rel_w,trb_m)
  write(*,'(A,ES10.3,A,A,A,ES10.3,A,A,A,A,A,ES10.3)') &
    'M_trabs_hi=',act_m,' ch=',trim(act_c),' rel=',rel_m,' rch=', &
    trim(rel_c),' why=',trim(rel_w),' trb=',trb_m
!
  write(*,'(A)') 'DONE'
end program drv_conv_exit
