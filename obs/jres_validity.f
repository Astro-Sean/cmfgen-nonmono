! jres_validity.f
!
! Elementwise field-validity helpers for the J_RESCUE machinery in
! COMP_JCONT_V4.INC.  The previous checks used
! T=MAXVAL(ABS(field)) followed by T.NE.T to catch NaNs; a MAXVAL
! reduction is not required to propagate NaN on every compiler
! lane order, so a non-finite element could be masked by a larger
! finite one.  These routines test every element with
! IEEE_IS_FINITE and classify the failure so the caller can report
! the trigger, not just that something failed.
!
! JRES_TRIG returns:
!   0 - field valid (all elements finite, max|a| <= CAP)
!   1 - NaN present (BAD = first non-finite index)
!   2 - Inf present (BAD = first non-finite index)
!   3 - finite but max|a| > CAP (BAD = argmax, TMAX = max|a|)
! A non-finite element takes precedence over a magnitude
! violation: the argmax of a poisoned field is meaningless.
!
! edited by Sean Brennan
!
	INTEGER FUNCTION JRES_TRIG(A,N,CAP,BAD,TMAX)
	USE SET_KIND_MODULE
	USE IEEE_ARITHMETIC
	IMPLICIT NONE
	INTEGER N,BAD,I
	REAL(KIND=LDP) A(N),CAP,TMAX
	JRES_TRIG=0
	BAD=0
	TMAX=0.0_LDP
	DO I=1,N
	  IF(.NOT. IEEE_IS_FINITE(A(I)))THEN
	    BAD=I
	    IF(IEEE_IS_NAN(A(I)))THEN
	      JRES_TRIG=1
	    ELSE
	      JRES_TRIG=2
	    END IF
	    RETURN
	  END IF
	END DO
	DO I=1,N
	  IF(ABS(A(I)) .GT. TMAX)THEN
	    TMAX=ABS(A(I))
	    BAD=I
	  END IF
	END DO
	IF(TMAX .GT. CAP)JRES_TRIG=3
	RETURN
	END
!
! JRES_FIN: finite-only check for fields whose validity has no
! magnitude bound (e.g. GEDD, RSQHNU snapshots).
!
	LOGICAL FUNCTION JRES_FIN(A,N)
	USE SET_KIND_MODULE
	USE IEEE_ARITHMETIC
	IMPLICIT NONE
	INTEGER N,I
	REAL(KIND=LDP) A(N)
	JRES_FIN=.TRUE.
	DO I=1,N
	  IF(.NOT. IEEE_IS_FINITE(A(I)))THEN
	    JRES_FIN=.FALSE.
	    RETURN
	  END IF
	END DO
	RETURN
	END
!
! JRES_DIAG_STATE: counters and file state shared by JRESC_EV and
! JRESC_TOT (SAVE in separate routines does not share).
!
	MODULE JRES_DIAG_STATE
	INTEGER :: LUE=0
	INTEGER :: NC_CNT=0,NL_CNT=0,NE_CNT=0,NEL_CNT=0
	CHARACTER(LEN=12) :: LSECT=' '
	END MODULE
!
! JRESC_EV: one machine-readable J_RESCUE_EVENTS record per rescue
! intervention or fatal abort (DO_J_RES_DIAG, default F).
! Positional fields so the file parses with a whitespace split;
! counters kept per section so the totals do not depend on the
! print throttle.  STATUS='UNKNOWN': OUT_FLUX accumulates renders
! in APPEND, the event file must not - one render per file.
!
! Record layout (JEV):
!   type indx ml fl sect first fg mom L fdnc trig ftrig bad
!   rjmax fedmax rjesmax rcap rcapit snapok djmax dfmax
!   R V SIG SIGT ltc1 luedd
!
! type:   NRESC (previous-frequency restore), NES_SAN (RJ_ES
!         substitution), or FATAL (invalid field, no fallback).
! trig:   RJ/RJ_ES trigger (0 ok, 1 NaN, 2 Inf, 3 over-cap).
! ftrig:  FEDD trigger on the same scale (cap = 2).
! snapok: previous-frequency snapshot passed validation.
! djmax/dfmax: max|field_before - field_after| of the restore.
! ltc1:   SIGMA_TRUE < -1 at the offending depth.
! luedd:  field reached the LU_EDD record (F only for FATAL).
!
	SUBROUTINE JRESC_EV(ETYPE,FREQ_INDX,ML,FL,SECTION,FIRST_FREQ,
	1    FG_NM,MOM_NM,L,FDNC,TRIG,FTRIG,BAD,RJMAX,FEDMAX,RJESMX,
	2    RCAP,RCAPIT,SNAPOK,DJMAX,DFMAX,RB,VB,SIGB,SIGTB,LUER)
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	CHARACTER(LEN=*) ETYPE,SECTION,FG_NM,MOM_NM
	INTEGER FREQ_INDX,ML,L,TRIG,FTRIG,BAD,LUER
	REAL(KIND=LDP) FL,RJMAX,FEDMAX,RJESMX,RCAP,RCAPIT
	REAL(KIND=LDP) DJMAX,DFMAX,RB,VB,SIGB,SIGTB
	LOGICAL FIRST_FREQ,FDNC,SNAPOK,LTC1
	CHARACTER(LEN=1) LUF
	IF(LUE .EQ. 0)THEN
	  OPEN(NEWUNIT=LUE,FILE='J_RESCUE_EVENTS',STATUS='UNKNOWN',
	1    ACTION='WRITE')
	  WRITE(LUE,'(A)')'! JEV type indx ml fl sect first fg mom L'//
	1    ' fdnc trig ftrig bad rjmax fedmax rjesmax rcap rcapit'//
	2    ' snapok djmax dfmax R V SIG SIGT ltc1 luedd'
	END IF
	LTC1=(SIGTB .LT. -1.0_LDP)
	LUF='T'
	IF(ETYPE .EQ. 'FATAL')LUF='F'
	IF(LSECT .NE. '' .AND. LSECT .NE. SECTION)
	1    CALL JRESC_SECTOT(LSECT)
	LSECT=SECTION
	IF(ETYPE .NE. 'FATAL')THEN
	  IF(SECTION .EQ. 'CONTINUUM')THEN
	    IF(ETYPE .EQ. 'NRESC')THEN
	      NC_CNT=NC_CNT+1
	    ELSE
	      NE_CNT=NE_CNT+1
	    END IF
	  ELSE
	    IF(ETYPE .EQ. 'NRESC')THEN
	      NL_CNT=NL_CNT+1
	    ELSE
	      NEL_CNT=NEL_CNT+1
	    END IF
	  END IF
	END IF
	WRITE(LUE,'(A,1X,A,2(1X,I6),1X,ES16.8,1X,A,1X,L1,1X,A,1X,A,'//
	1   '1X,I3,1X,L1,2(1X,I2),1X,I4,5(1X,ES12.4),1X,L1,'//
	2   '2(1X,ES12.4),4(1X,ES12.4),1X,L1,1X,A)')
	3   'JEV',ETYPE,FREQ_INDX,ML,FL,TRIM(SECTION),FIRST_FREQ,
	4   TRIM(FG_NM),TRIM(MOM_NM),L,FDNC,TRIG,FTRIG,BAD,RJMAX,
	5   FEDMAX,RJESMX,RCAP,RCAPIT,SNAPOK,DJMAX,DFMAX,RB,VB,SIGB,
	6   SIGTB,LTC1,LUF
	RETURN
	END
!
! JRESC_SECTOT: pending-section JTOT flush; called on a section
! transition and by JRESC_TOT so the last section is not lost when
! it ends without a following event.
!
	SUBROUTINE JRESC_SECTOT(SECTION)
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	CHARACTER(LEN=*) SECTION
	IF(LUE .EQ. 0)RETURN
	IF(SECTION .EQ. 'CONTINUUM')THEN
	  WRITE(LUE,'(A,1X,A,2(1X,A,I8))')'JTOT',SECTION,
	1    'nresc=',NC_CNT,'nes_san=',NE_CNT
	ELSE
	  WRITE(LUE,'(A,1X,A,2(1X,A,I8))')'JTOT',SECTION,
	1    'nresc=',NL_CNT,'nes_san=',NEL_CNT
	END IF
	RETURN
	END
!
! JRESC_TOT: end-of-render record.  Flushes the last section's JTOT
! and writes cumulative per-section totals; opens the file even
! with zero events so a clean render still proves the flag was on.
!
	SUBROUTINE JRESC_TOT()
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	IF(LUE .EQ. 0)THEN
	  OPEN(NEWUNIT=LUE,FILE='J_RESCUE_EVENTS',STATUS='UNKNOWN',
	1    ACTION='WRITE')
	  WRITE(LUE,'(A)')'! JEV type indx ml fl sect first fg mom L'//
	1    ' fdnc trig ftrig bad rjmax fedmax rjesmax rcap rcapit'//
	2    ' snapok djmax dfmax R V SIG SIGT ltc1 luedd'
	END IF
	IF(LSECT .NE. '')CALL JRESC_SECTOT(LSECT)
	WRITE(LUE,'(A,4(1X,A,I8))')'JEND',
	1    'nresc_cont=',NC_CNT,'nesan_cont=',NE_CNT,
	2    'nresc_line=',NL_CNT,'nesan_line=',NEL_CNT
	RETURN
	END
