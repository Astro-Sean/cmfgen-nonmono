! knop_spec_gate.f
!
! P8.3 KNOP spectrum-output gate (17-Nov-2026).  Classifies the
! model's lagged-field and solver state from the files the
! population run already leaves in the model directory, before
! CMF_FLUX_SUB_V5 writes OBSFLUX or OBSFRAME.  No new convergence
! counter is introduced; the predicate consumes existing state.
!
! Evidence, searched in the parent directory (batobs.sh runs this
! executable in MODEL_DIR/obs, so '..' is the model run dir) with a
! local-directory fallback for hand runs:
!   GS_STATE               - written every SOLVE_FOR_POPS call under
!                            DO_NONMONO_KNOP (solve_for_pops.f); keys
!                            n_frozen, frozen, mode, n_mature.
!   .CMFGEN_KNOP_STALLED   - watchdog marker touched by the generated
!                            slurm when n_frozen reaches KNOP_STALL_MAX
!                            (submit_cmfgen.sh, default 60).
!   manifest.json          - run manifest (cmf_run_manifest.py); the
!                            last recorded exit_code / termination
!                            catches a killed or interrupted model
!                            that a stale-but-healthy GS_STATE cannot.
!
! States (see DO_KNOP_SPEC_GATE in cmf_flux_cntrl_var_mod.f):
!   1 KNOP_CONVERGED            - frozen=F, n_frozen=0, mode not
!                                 LAMBDA_GS
!   2 KNOP_MATURE_STALLED       - stall evidence and n_mature>0
!   3 KNOP_IMMATURE             - still frozen / freeze streak active
!   4 KNOP_EXECUTION_INCOMPLETE - GS_STATE absent or malformed, or
!                                 manifest records an abnormal end
! States 2-4 block production output; state 0 means the gate did not
! apply (non-KNOP model) and is handled by the caller's flag check.
!
! mom_refresh_conv, ba_rej_consec and lum_ratio are written to the
! KNOP_SPEC_STATUS evidence file but do not gate: mom_refresh_conv
! stays F on the FG-refresh KNOP route, and no lum_ratio threshold
! has been calibrated against a population runaway yet.
!
! edited by Sean Brennan
!
	SUBROUTINE KNOP_SPEC_GATE_EVAL(LUER)
	USE SET_KIND_MODULE
	USE CMF_FLUX_CNTRL_VAR_MOD
	IMPLICIT NONE
	INTEGER LUER
!
	INTEGER LUG,IOER,IP,IP2,IP3
	INTEGER N_FRZ,N_MAT,M_EX
	LOGICAL FEX,FRZ,STALL_MK,GOT_NF,MAN_BAD
	CHARACTER(LEN=256) LINE
	CHARACTER(LEN=64) MODE,MAN_TERM
!
! Threshold mirrors the submit-script watchdog default
! (KNOP_STALL_MAX=60, cmfgen_base/submit_cmfgen.sh).  The stall
! marker check catches runs whose wrapper used a custom value;
! this literal only matters for hand reruns outside the wrapper.
!
	INTEGER, PARAMETER :: NF_STALL=60
!
	KNOP_SPEC_STATE=0
	KNOP_SPEC_CLASS='NOT_APPLICABLE'
	KNOP_SPEC_ALLOW=.TRUE.
	IF(.NOT. DO_NONMONO_KNOP)RETURN
!
	N_FRZ=-1
	N_MAT=0
	FRZ=.FALSE.
	GOT_NF=.FALSE.
	MODE=''
	STALL_MK=.FALSE.
	MAN_BAD=.FALSE.
	M_EX=0
	MAN_TERM=''
!
	INQUIRE(FILE='../GS_STATE',EXIST=FEX)
	IF(FEX)THEN
	  OPEN(NEWUNIT=LUG,FILE='../GS_STATE',STATUS='OLD',
     1      ACTION='READ',IOSTAT=IOER)
	ELSE
	  INQUIRE(FILE='GS_STATE',EXIST=FEX)
	  IF(FEX)OPEN(NEWUNIT=LUG,FILE='GS_STATE',STATUS='OLD',
     1      ACTION='READ',IOSTAT=IOER)
	END IF
	IF(FEX)THEN
	  IF(IOER .EQ. 0)THEN
	    IOER=0
	    DO WHILE(IOER .EQ. 0)
	      READ(LUG,'(A)',IOSTAT=IOER)LINE
	      IF(IOER .NE. 0)EXIT
! 'n_frozen=' is an exact 9-char head match so 'n_frozen_tot='
! cannot alias it; 'frozen=' rides inside the 'iter=... frozen= X'
! line, so it is matched by position preceded by a space, which the
! 'n_frozen=' head cannot produce.
	      IF(LINE(1:9) .EQ. 'n_frozen=')THEN
	        READ(LINE(10:),*,IOSTAT=IP)N_FRZ
	        IF(IP .EQ. 0)GOT_NF=.TRUE.
	      ELSE IF(LINE(1:9) .EQ. 'n_mature=')THEN
	        READ(LINE(10:),*,IOSTAT=IP)N_MAT
	      ELSE IF(LINE(1:5) .EQ. 'mode=')THEN
	        MODE=ADJUSTL(LINE(6:))
	      ELSE
	        IP=INDEX(LINE,'frozen=')
	        IF(IP .GT. 1 .AND. LINE(IP-1:IP-1) .EQ. ' ')
     1          READ(LINE(IP+7:),*,IOSTAT=IP2)FRZ
	      END IF
	    END DO
	    CLOSE(LUG)
	  ELSE
	    GOT_NF=.FALSE.
	  END IF
	END IF
!
	INQUIRE(FILE='../.CMFGEN_KNOP_STALLED',EXIST=STALL_MK)
	IF(.NOT. STALL_MK)
     1    INQUIRE(FILE='.CMFGEN_KNOP_STALLED',EXIST=STALL_MK)
!
! Manifest: the last recorded segment decides.  A 'knop_stall'
! termination is left for the GS_STATE checks, which carry the
! richer stall classification; other non-normal terminations and
! any nonzero exit_code flag an incomplete model.
!
	INQUIRE(FILE='../manifest.json',EXIST=FEX)
	IF(.NOT. FEX)INQUIRE(FILE='manifest.json',EXIST=FEX)
	IF(FEX)THEN
	  OPEN(NEWUNIT=LUG,FILE='../manifest.json',STATUS='OLD',
     1      ACTION='READ',IOSTAT=IOER)
	  IF(IOER .NE. 0)
     1      OPEN(NEWUNIT=LUG,FILE='manifest.json',STATUS='OLD',
     2        ACTION='READ',IOSTAT=IOER)
	  IF(IOER .EQ. 0)THEN
	    DO WHILE(IOER .EQ. 0)
	      READ(LUG,'(A)',IOSTAT=IOER)LINE
	      IF(IOER .NE. 0)EXIT
	      IP=INDEX(LINE,'"exit_code"')
	      IF(IP .GT. 0)THEN
	        IP2=INDEX(LINE(IP:),':')
	        IF(IP2 .GT. 0)READ(LINE(IP+IP2:),*,IOSTAT=IP3)M_EX
	      END IF
	      IP=INDEX(LINE,'"termination"')
	      IF(IP .GT. 0)THEN
	        IP2=INDEX(LINE(IP:),':')
	        IF(IP2 .GT. 0)THEN
	          IP2=IP+IP2-1
	          IP3=INDEX(LINE(IP2:),'"')
	          IF(IP3 .GT. 0)THEN
	            IP3=IP2+IP3-1
	            IP2=INDEX(LINE(IP3+1:),'"')
	            IF(IP2 .GT. 0)MAN_TERM=LINE(IP3+1:IP3+IP2-1)
	          END IF
	        END IF
	      END IF
	    END DO
	    CLOSE(LUG)
	  END IF
	END IF
	IF(M_EX .NE. 0 .AND. MAN_TERM .NE. 'knop_stall')
     1    MAN_BAD=.TRUE.
	IF(MAN_TERM .NE. '' .AND. MAN_TERM .NE. 'normal_exit' .AND.
     1    MAN_TERM .NE. 'graceful_stop' .AND.
     2    MAN_TERM .NE. 'knop_stall')MAN_BAD=.TRUE.
!
	IF(.NOT. GOT_NF)THEN
	  KNOP_SPEC_STATE=4
	  KNOP_SPEC_CLASS='KNOP_EXECUTION_INCOMPLETE'
	ELSE IF(STALL_MK .OR. N_FRZ .GE. NF_STALL)THEN
	  IF(N_MAT .GT. 0)THEN
	    KNOP_SPEC_STATE=2
	    KNOP_SPEC_CLASS='KNOP_MATURE_STALLED'
	  ELSE
	    KNOP_SPEC_STATE=3
	    KNOP_SPEC_CLASS='KNOP_IMMATURE'
	  END IF
	ELSE IF(FRZ .OR. N_FRZ .GT. 0 .OR. MODE .EQ. 'LAMBDA_GS')THEN
	  KNOP_SPEC_STATE=3
	  KNOP_SPEC_CLASS='KNOP_IMMATURE'
	ELSE IF(MAN_BAD)THEN
	  KNOP_SPEC_STATE=4
	  KNOP_SPEC_CLASS='KNOP_EXECUTION_INCOMPLETE'
	ELSE
	  KNOP_SPEC_STATE=1
	  KNOP_SPEC_CLASS='KNOP_CONVERGED'
	END IF
	KNOP_SPEC_ALLOW=(KNOP_SPEC_STATE .LE. 1)
!
! Evidence record next to the spectrum files: enough to audit why a
! render was allowed or refused without rerunning the gate.
!
	OPEN(NEWUNIT=LUG,FILE='KNOP_SPEC_STATUS',STATUS='REPLACE',
     1    IOSTAT=IOER)
	IF(IOER .EQ. 0)THEN
	  WRITE(LUG,'(A,A)')'class=',TRIM(KNOP_SPEC_CLASS)
	  WRITE(LUG,'(A,I0)')'n_frozen=',N_FRZ
	  WRITE(LUG,'(A,L2)')'frozen=',FRZ
	  WRITE(LUG,'(A,A)')'mode=',TRIM(MODE)
	  WRITE(LUG,'(A,I0)')'n_mature=',N_MAT
	  WRITE(LUG,'(A,L2)')'stall_marker=',STALL_MK
	  WRITE(LUG,'(A,I0)')'manifest_exit=',M_EX
	  WRITE(LUG,'(A,A)')'manifest_term=',TRIM(MAN_TERM)
	  CLOSE(LUG)
	END IF
!
	WRITE(LUER,'(A,A,A,I0,A,L1)')' KNOP spectrum gate: ',
     1    TRIM(KNOP_SPEC_CLASS),' (state ',KNOP_SPEC_STATE,
     2    ', allowed=',KNOP_SPEC_ALLOW,')'
	RETURN
	END SUBROUTINE KNOP_SPEC_GATE_EVAL
