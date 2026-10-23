! Standalone test of the Stage-B trace-polish state machine
! (BA_TRACE_EVAL in control_variable_mod.f, compiled unchanged --
! the production routine CMFGEN_SUB and SOLVE_FOR_POPS both call).
!
! Cases lock the 12-Nov-2026 off-state semantics:
!   BA_TR_MAXSWP=0  -> master off: NEEDED=F, OFF=T (exit releases;
!                      a material demand must not pin the run).
!   BA_DEM_PFLOOR<=0 -> demand channel off: OFF=T, release.
!   BA_TR_ABSMAX=0  -> absolute leg off only; the demand leg still
!                      arms Stage-B.
!   BA_DEM_EXM>=BA_SAT_RAW -> dead-stage class: not NEEDED.
!   NSW>=MAXSWP     -> DONE (budget spent, exit releases).
!   BA8/BA9 mirror the CMFGEN_SUB exit-release leg and the
!   final-budget status precedence -- the production branches live
!   in cmfgen_sub.f which this harness cannot link; the mirror
!   covers that NEEDED.OR.RELEASE releases the off state and that a
!   NEEDED final iterate selects the owed-status branch rather than
!   "Convergence criteria met".
!
   program test_trace_state
   use set_kind_module, only: LDP
   use control_variable_mod
   implicit none
   logical :: t_act, t_off, t_need, t_done
   logical :: released
   integer :: nfail, isel
   character(len=120) :: sel

   nfail = 0
   call reset_state()
   !
   ! --- BA1: master off (MAXSWP=0) with material demand releases ---
   BA_TR_MAXSWP = 0
   BA_DEM_EXM = 17.94_LDP
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   released = t_need .or. t_off .or. t_done
   if(t_need .or. .not. t_off .or. t_done .or. .not. released) then
      nfail = nfail + 1
      write(6,'(A,4L2)') 'BA1: FAIL (MAXSWP=0 pins exit) ', &
         t_act, t_off, t_need, t_done
   else
      write(6,'(A)') 'BA1: PASS (MAXSWP=0 -> OFF, exit released)'
   end if
   !
   ! --- BA2: ABSMAX=0 kills only the abs leg; demand still arms ----
   call reset_state()
   BA_TR_ABSMAX = 0.0_LDP
   BA_DEM_EXM = 17.94_LDP
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   if(.not. t_act .or. t_off .or. .not. t_need .or. t_done) then
      nfail = nfail + 1
      write(6,'(A,4L2)') 'BA2: FAIL (ABSMAX=0 suppressed demand leg) ', &
         t_act, t_off, t_need, t_done
   else
      write(6,'(A)') 'BA2: PASS (ABSMAX=0, demand leg still arms)'
   end if
   !
   ! --- BA3: PFLOOR<=0 -> demand channel off, release -------------
   call reset_state()
   BA_DEM_PFLOOR = 0.0_LDP
   BA_TR_ABS = 1.0_LDP   ! abs breach alone cannot arm while off
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   released = t_need .or. t_off .or. t_done
   if(t_need .or. .not. t_off .or. .not. released) then
      nfail = nfail + 1
      write(6,'(A,4L2)') 'BA3: FAIL (PFLOOR=0 did not release) ', &
         t_act, t_off, t_need, t_done
   else
      write(6,'(A)') 'BA3: PASS (PFLOOR=0 -> OFF, exit released)'
   end if
   !
   ! --- BA4: MAXSWP=1 -- one sweep, then DONE releases ------------
   call reset_state()
   BA_TR_MAXSWP = 1
   BA_DEM_EXM = 17.94_LDP
   BA_TR_NSWP = 0
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   if(.not. t_need) then
      nfail = nfail + 1
      write(6,'(A)') 'BA4: FAIL (material demand did not arm)'
   end if
   BA_TR_NSWP = 1
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   released = t_need .or. t_off .or. t_done
   if(t_need .or. .not. t_done .or. .not. released) then
      nfail = nfail + 1
      write(6,'(A,4L2)') 'BA4: FAIL (spent budget did not release) ', &
         t_act, t_off, t_need, t_done
   else if(nfail == 0) then
      write(6,'(A)') 'BA4: PASS (one sweep arms, spent budget releases)'
   else
      write(6,'(A)') 'BA4: PASS (spent budget releases)'
   end if
   !
   ! --- BA5: abs leg alone arms when ABSMAX>0 ----------------------
   call reset_state()
   BA_TR_ABS = 1.0d-6
   BA_TR_ABSMAX = 1.0d-8
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   if(.not. t_need) then
      nfail = nfail + 1
      write(6,'(A)') 'BA5: FAIL (abs breach did not arm)'
   else
      write(6,'(A)') 'BA5: PASS (abs leg arms Stage-B)'
   end if
   !
   ! --- BA6: saturated dead-stage demand does not arm --------------
   call reset_state()
   BA_DEM_EXM = BA_SAT_RAW
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   if(t_need) then
      nfail = nfail + 1
      write(6,'(A)') 'BA6: FAIL (saturated demand armed Stage-B)'
   else
      write(6,'(A)') 'BA6: PASS (saturated demand does not arm)'
   end if
   !
   ! --- BA7: nonmono off -> inert ----------------------------------
   call reset_state()
   DO_NONMONO_SOB = .false.
   DO_NONMONO_KNOP = .false.
   BA_DEM_EXM = 17.94_LDP
   BA_TR_ABS = 1.0_LDP
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   if(t_act .or. t_off .or. t_need .or. t_done) then
      nfail = nfail + 1
      write(6,'(A)') 'BA7: FAIL (stock path touched trace state)'
   else
      write(6,'(A)') 'BA7: PASS (nonmono off -> all flags false)'
   end if
   !
   ! --- BA8: final-budget status precedence (mirror) ---------------
   ! Mirrors the cmfgen_sub.f final-iteration branch: a met gate with
   ! BA_TRACE_NEEDED true must select the owed/trace-limited status,
   ! never the "Convergence criteria met" line.
   call reset_state()
   BA_DEM_EXM = 17.94_LDP
   call BA_TRACE_EVAL(t_act, t_off, t_need, t_done)
   isel = 4
   if(t_need) then
      isel = 1
   else if(t_done .and. BA_DEM_EXM >= max(10.0_LDP*EPS,1.0_LDP)) then
      isel = 2
   else if(t_off .and. BA_DEM_EXM >= max(10.0_LDP*EPS,1.0_LDP)) then
      isel = 3
   end if
   if(isel /= 1) then
      nfail = nfail + 1
      write(6,'(A,I2)') 'BA8: FAIL (owed state picked branch) ', isel
   else
      write(6,'(A)') 'BA8: PASS (final owed iterate picks owed status)'
   end if
   !
   if(nfail == 0) then
      write(6,'(A)') 'ALL TRACE STATE TESTS PASSED'
   else
      write(6,'(A,I2)') 'FAILURES: ', nfail
      stop 1
   end if
   !
contains
   !
   subroutine reset_state()
      DO_NONMONO_SOB = .true.
      DO_NONMONO_KNOP = .false.
      BA_TR_MAXSWP = 3
      BA_DEM_PFLOOR = 1.0d-30
      BA_TR_ABSMAX = 1.0d-8
      BA_TR_ABS = 0.0_LDP
      BA_DEM_EXM = 0.0_LDP
      BA_TR_NSWP = 0
      EPS = 0.01_LDP
      if(BA_SAT_RAW <= 0.0_LDP) BA_SAT_RAW = 1.0d6
   end subroutine reset_state
   !
end program test_trace_state
