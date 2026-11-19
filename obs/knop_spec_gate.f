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
!                            n_frozen, run_id, frozen, mode, n_mature.
!   .CMFGEN_KNOP_STALLED   - watchdog marker written by the generated
!                            slurm when n_frozen reaches KNOP_STALL_MAX
!                            (submit_cmfgen.sh, default 60); its content
!                            is the slurm job id that stalled, so a
!                            stale marker from an earlier execution
!                            cannot stall-classify a later run.
!   manifest.json          - run manifest (cmf_run_manifest.py).  The
!                            authoritative model-execution record is the
!                            last runs[] segment labelled "FINAL"; an
!                            open OBS segment (null exit_code) is the
!                            render in progress, not the model result.
!
! States (see DO_KNOP_SPEC_GATE in cmf_flux_cntrl_var_mod.f):
!   1 KNOP_CONVERGED            - frozen=F, n_frozen=0, mode not
!                                 LAMBDA_GS, provenance verified
!   2 KNOP_MATURE_STALLED       - stall evidence and n_mature>0
!   3 KNOP_IMMATURE             - still frozen / freeze streak active
!   4 KNOP_EXECUTION_INCOMPLETE - GS_STATE absent or malformed, or the
!                                 FINAL segment is open or ended
!                                 abnormally
!   5 KNOP_PROVENANCE_UNVERIFIED- solver state looks converged but the
!                                 record cannot be bound to the
!                                 execution being rendered: no manifest,
!                                 no FINAL segment, or run_id mismatch
!                                 (e.g. a stale GS_STATE left by a warm
!                                 restart that died before iterating)
! Precedence: solver-state evidence is trusted only after execution
! integrity is established, and a would-be CONVERGED is downgraded to
! 5 whenever provenance cannot be verified.  States 2-4 set
! KNOP_SPEC_ALLOW=F (blocked under DO_KNOP_SPEC_GATE=T); state 5 sets
! ALLOW=T but the state>1 marker rules force every written spectrum
! to carry the '! knop_spec_status:' label.  State 0 means the gate
! did not apply (non-KNOP model).
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
	INTEGER N_FRZ,N_MAT,M_EX,SK
	LOGICAL FEX,FRZ,GOT_NF,MAN_FEX
	LOGICAL STALL_MK,MK_APP,FN_SEEN,FN_OPEN,FN_BAD,PRV_OK
	CHARACTER(LEN=256) LINE
	CHARACTER(LEN=64) MODE,MAN_TERM,SEG_LBL,FN_TERM,FN_JOB
	CHARACTER(LEN=64) GS_RUNID,MKID
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
	GS_RUNID=''
	STALL_MK=.FALSE.
	MK_APP=.FALSE.
	MAN_FEX=.FALSE.
	SEG_LBL=''
	FN_SEEN=.FALSE.
	FN_OPEN=.FALSE.
	FN_BAD=.FALSE.
	PRV_OK=.FALSE.
	M_EX=0
	MAN_TERM=''
	FN_TERM=''
	FN_JOB=''
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
	      ELSE IF(LINE(1:7) .EQ. 'run_id=')THEN
	        GS_RUNID=ADJUSTL(LINE(8:))
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
! Stall marker: content is the slurm job id of the execution that
! stalled (post 17-Nov-2026 wrapper); an empty file is a legacy
! marker with no binding.  A bound marker applies only when its id
! equals the GS_STATE run_id; a legacy marker applies only when
! the state carries no run_id to contradict it (fail closed where
! nothing can be disambiguated).
!
	INQUIRE(FILE='../.CMFGEN_KNOP_STALLED',EXIST=STALL_MK)
	IF(.NOT. STALL_MK)
     1    INQUIRE(FILE='.CMFGEN_KNOP_STALLED',EXIST=STALL_MK)
	MKID=''
	IF(STALL_MK)THEN
	  OPEN(NEWUNIT=LUG,FILE='../.CMFGEN_KNOP_STALLED',
     1      STATUS='OLD',ACTION='READ',IOSTAT=IOER)
	  IF(IOER .NE. 0)OPEN(NEWUNIT=LUG,
     1      FILE='.CMFGEN_KNOP_STALLED',STATUS='OLD',
     2      ACTION='READ',IOSTAT=IOER)
	  IF(IOER .EQ. 0)THEN
	    READ(LUG,'(A)',IOSTAT=IOER)MKID
	    MKID=ADJUSTL(MKID)
	    CLOSE(LUG)
	  END IF
	  IF(MKID .NE. '' .AND. MKID .NE. 'local')THEN
	    MK_APP=(GS_RUNID .EQ. MKID)
	  ELSE
	    MK_APP=(GS_RUNID .EQ. '' .OR. GS_RUNID .EQ. 'local')
	  END IF
	END IF
!
! Manifest: the last runs[] segment labelled "FINAL" is the model
! execution record.  Open segments carry null exit_code/termination
! and are skipped for classification except that an open FINAL
! means the model is still executing.  'knop_stall' is left for the
! solver-state checks, which carry the richer classification.
!
	INQUIRE(FILE='../manifest.json',EXIST=MAN_FEX)
	IF(.NOT. MAN_FEX)INQUIRE(FILE='manifest.json',EXIST=MAN_FEX)
	IF(MAN_FEX)THEN
	  OPEN(NEWUNIT=LUG,FILE='../manifest.json',STATUS='OLD',
     1      ACTION='READ',IOSTAT=IOER)
	  IF(IOER .NE. 0)
     1      OPEN(NEWUNIT=LUG,FILE='manifest.json',STATUS='OLD',
     2        ACTION='READ',IOSTAT=IOER)
	  IF(IOER .EQ. 0)THEN
	    DO WHILE(IOER .EQ. 0)
	      READ(LUG,'(A)',IOSTAT=IOER)LINE
	      IF(IOER .NE. 0)EXIT
	      IP=INDEX(LINE,'"label"')
	      IF(IP .GT. 0)THEN
	        IP2=INDEX(LINE,':')
	        IF(IP2 .GT. 0)THEN
	          IP3=INDEX(LINE(IP2+1:),'"')
	          IF(IP3 .GT. 0)THEN
	            IP3=IP2+IP3
	            IP2=INDEX(LINE(IP3+1:),'"')
	            IF(IP2 .GT. 0)SEG_LBL=LINE(IP3+1:IP3+IP2-1)
	          END IF
	        END IF
	      END IF
	      IF(SEG_LBL .EQ. 'FINAL')THEN
	        IP=INDEX(LINE,'"slurm_job_id"')
	        IF(IP .GT. 0)THEN
	          IP2=INDEX(LINE,':')
	          IF(IP2 .GT. 0)THEN
	            IP3=INDEX(LINE(IP2+1:),'"')
	            IF(IP3 .GT. 0)THEN
	              IP3=IP2+IP3
	              IP2=INDEX(LINE(IP3+1:),'"')
	              IF(IP2 .GT. 0)FN_JOB=LINE(IP3+1:IP3+IP2-1)
	            END IF
	          END IF
	        END IF
	        IP=INDEX(LINE,'"exit_code"')
	        IF(IP .GT. 0)THEN
	          FN_SEEN=.TRUE.
	          IF(INDEX(LINE,'null') .GT. 0)THEN
	            FN_OPEN=.TRUE.
	          ELSE
	            IP2=INDEX(LINE(IP:),':')
	            IF(IP2 .GT. 0)READ(LINE(IP+IP2:),*,
     1            IOSTAT=IP3)M_EX
	          END IF
	        END IF
	        IP=INDEX(LINE,'"termination"')
	        IF(IP .GT. 0)THEN
	          IF(INDEX(LINE,'null') .GT. 0)THEN
	            FN_OPEN=.TRUE.
	            FN_TERM=''
	          ELSE
	            IP2=INDEX(LINE(IP:),':')
	            IF(IP2 .GT. 0)THEN
	              IP2=IP+IP2-1
	              IP3=INDEX(LINE(IP2:),'"')
	              IF(IP3 .GT. 0)THEN
	                IP3=IP2+IP3-1
	                IP2=INDEX(LINE(IP3+1:),'"')
	                IF(IP2 .GT. 0)THEN
	                  FN_TERM=LINE(IP3+1:IP3+IP2-1)
	                  MAN_TERM=FN_TERM
	                END IF
	              END IF
	            END IF
	          END IF
	        END IF
	      END IF
	    END DO
	    CLOSE(LUG)
	  END IF
	END IF
! A FINAL segment that ended normally is good execution evidence;
! an open or abnormally-ended one is incomplete.  'knop_stall'
! defers to the solver-state stall classification below.
	FN_BAD=FN_SEEN .AND. .NOT. FN_OPEN .AND.
     1    ((M_EX .NE. 0 .AND. FN_TERM .NE. 'knop_stall') .OR.
     2     (FN_TERM .NE. '' .AND. FN_TERM .NE. 'normal_exit' .AND.
     3      FN_TERM .NE. 'graceful_stop' .AND.
     4      FN_TERM .NE. 'knop_stall'))
!
	IF(.NOT. GOT_NF)THEN
	  KNOP_SPEC_STATE=4
	  KNOP_SPEC_CLASS='KNOP_EXECUTION_INCOMPLETE'
	ELSE IF(FN_OPEN .OR. FN_BAD)THEN
	  KNOP_SPEC_STATE=4
	  KNOP_SPEC_CLASS='KNOP_EXECUTION_INCOMPLETE'
	ELSE
	  SK=1
	  IF(MK_APP .OR. N_FRZ .GE. NF_STALL .OR.
     1      FN_TERM .EQ. 'knop_stall')THEN
	    IF(N_MAT .GT. 0)THEN
	      SK=2
	    ELSE
	      SK=3
	    END IF
	  ELSE IF(FRZ .OR. N_FRZ .GT. 0 .OR. MODE .EQ. 'LAMBDA_GS')THEN
	    SK=3
	  END IF
! Provenance: a converged-looking state must be bound to the model
! execution being rendered.  run_id comes from the executing job's
! SLURM_JOB_ID; the FINAL segment records the same id.  Absent
! manifest, absent run_id ('local' non-slurm runs), or a mismatch
! all downgrade CONVERGED to PROVENANCE_UNVERIFIED -- rendered but
! marked, never silently clean.
	  IF(SK .EQ. 1)THEN
	    PRV_OK=MAN_FEX .AND. FN_SEEN .AND. .NOT. FN_OPEN .AND.
     1        FN_TERM .EQ. 'normal_exit' .AND.
     2        GS_RUNID .NE. '' .AND. GS_RUNID .NE. 'local' .AND.
     3        FN_JOB .EQ. GS_RUNID
	    IF(.NOT. PRV_OK)SK=5
	  END IF
	  KNOP_SPEC_STATE=SK
	  IF(SK .EQ. 2)KNOP_SPEC_CLASS='KNOP_MATURE_STALLED'
	  IF(SK .EQ. 3)KNOP_SPEC_CLASS='KNOP_IMMATURE'
	  IF(SK .EQ. 5)KNOP_SPEC_CLASS='KNOP_PROVENANCE_UNVERIFIED'
	  IF(SK .EQ. 1)KNOP_SPEC_CLASS='KNOP_CONVERGED'
	END IF
	KNOP_SPEC_ALLOW=(KNOP_SPEC_STATE .LE. 1 .OR.
     1    KNOP_SPEC_STATE .EQ. 5)
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
	  WRITE(LUG,'(A,A)')'run_id=',TRIM(GS_RUNID)
	  WRITE(LUG,'(A,L2)')'stall_marker=',STALL_MK
	  WRITE(LUG,'(A,L2)')'marker_applies=',MK_APP
	  WRITE(LUG,'(A,I0)')'final_exit=',M_EX
	  WRITE(LUG,'(A,A)')'final_term=',TRIM(FN_TERM)
	  WRITE(LUG,'(A,A)')'final_job=',TRIM(FN_JOB)
	  WRITE(LUG,'(A,L2)')'final_open=',FN_OPEN
	  WRITE(LUG,'(A,L2)')'provenance_ok=',PRV_OK
	  CLOSE(LUG)
	END IF
!
	WRITE(LUER,'(A,A,A,I0,A,L1)')' KNOP spectrum gate: ',
     1    TRIM(KNOP_SPEC_CLASS),' (state ',KNOP_SPEC_STATE,
     2    ', allowed=',KNOP_SPEC_ALLOW,')'
	RETURN
	END SUBROUTINE KNOP_SPEC_GATE_EVAL
