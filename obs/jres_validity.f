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
! Sealed flags make JTOT emission idempotent: the continuum totals are
! flushed at the CONTINUUM->LINE transition inside CMF_FLUX_SUB_V5 and
! again from JRESC_EV on the first LINE event and from JRESC_TOT at
! the RETURNs; a render that aborts mid-section (e.g. the OBS_FRAME
! undercoverage STOP) still carries the section JTOT already emitted.
	LOGICAL :: CONT_SEALED=.FALSE.,LINE_SEALED=.FALSE.
! J_RESCUE_DUMP unit + counter for JRESC_DMP (bounded full-field
! dump; JEV stays one line per event so parsers are not flooded).
	INTEGER :: LUD=0,NDMP=0
! J_RESCUE_FREQ unit + NSCND cap for cadence ('SCAN') profile
! dumps: separate from NDMP so trigger dumps are not displaced.
	INTEGER :: LUFQ=0,NSCND=0
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
! prevj/prevf: max|JNU_PREV| and max|FEDD_PREV| - the post-restore
!         maxima, so RJMAX_AFTER_RESTORE and the relative cap
!         1E4*prevj are recoverable without another pass (0 when no
!         snapshot applies, e.g. NES_SAN or FIRST_FREQ).
! etab/chib/zetab: ETA, CHI, ZETA at the offending depth - the
!         local emissivity/opacity/thermal-source scale the cap is
!         measured against.
!
	SUBROUTINE JRESC_EV(ETYPE,FREQ_INDX,ML,FL,SECTION,FIRST_FREQ,
	1    FG_NM,MOM_NM,L,FDNC,TRIG,FTRIG,BAD,RJMAX,FEDMAX,RJESMX,
	2    RCAP,RCAPIT,SNAPOK,DJMAX,DFMAX,RB,VB,SIGB,SIGTB,
	3    PREVJ,PREVF,ETAB,CHIB,ZETAB,LUER)
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	CHARACTER(LEN=*) ETYPE,SECTION,FG_NM,MOM_NM
	INTEGER FREQ_INDX,ML,L,TRIG,FTRIG,BAD,LUER
	REAL(KIND=LDP) FL,RJMAX,FEDMAX,RJESMX,RCAP,RCAPIT
	REAL(KIND=LDP) DJMAX,DFMAX,RB,VB,SIGB,SIGTB
	REAL(KIND=LDP) PREVJ,PREVF,ETAB,CHIB,ZETAB
	LOGICAL FIRST_FREQ,FDNC,SNAPOK,LTC1
	CHARACTER(LEN=1) LUF
	IF(LUE .EQ. 0)THEN
	  OPEN(NEWUNIT=LUE,FILE='J_RESCUE_EVENTS',STATUS='UNKNOWN',
	1    ACTION='WRITE')
	  WRITE(LUE,'(A)')'! JEV type indx ml fl sect first fg mom L'//
	1    ' fdnc trig ftrig bad rjmax fedmax rjesmax rcap rcapit'//
	2    ' snapok djmax dfmax R V SIG SIGT ltc1 luedd'//
	3    ' prevj prevf etab chib zetab'
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
	2   '2(1X,ES12.4),4(1X,ES12.4),1X,L1,1X,A,5(1X,ES12.4))')
	3   'JEV',ETYPE,FREQ_INDX,ML,FL,TRIM(SECTION),FIRST_FREQ,
	4   TRIM(FG_NM),TRIM(MOM_NM),L,FDNC,TRIG,FTRIG,BAD,RJMAX,
	5   FEDMAX,RJESMX,RCAP,RCAPIT,SNAPOK,DJMAX,DFMAX,RB,VB,SIGB,
	6   SIGTB,LTC1,LUF,PREVJ,PREVF,ETAB,CHIB,ZETAB
	RETURN
	END
!
! JRESC_DMP: bounded full-field dump for the first 200 NRESC events.
! The JEV line carries maxima only; the restore's depth structure
! (isolated outer-cell excess vs whole-profile change) needs the
! actual vectors.  J_RESCUE_DUMP keeps them out of the event file so
! JEV stays one-line-per-event.  200 blocks is enough to characterize
! the profile shape without an I/O flood (~3 MB).
!
! Per block: a JDMP header (indx fl sect bad nd) then ND rows of
!   i R V SIG SIGT RJ_before JNU_PREV FEDD FEDD_PREV ETA CHI ZETA
! RJ_before is the pre-restore formal field; the restore has not run
! at this call site.
!
	SUBROUTINE JRESC_DMP(ETYPE,FREQ_INDX,FL,SECTION,BAD,ND,
	1    RJ,JP,FEDD,FP,ETA,CHI,ZETA,R,V,SIG,SIGT,LUER)
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	CHARACTER(LEN=*) ETYPE,SECTION
	INTEGER FREQ_INDX,BAD,ND,LUER,I
	REAL(KIND=LDP) FL
	REAL(KIND=LDP) RJ(ND),JP(ND),FEDD(ND),FP(ND)
	REAL(KIND=LDP) ETA(ND),CHI(ND),ZETA(ND),R(ND),V(ND)
	REAL(KIND=LDP) SIG(ND),SIGT(ND)
	IF(ETYPE .EQ. 'SCAN')THEN
!	  Cadence dumps (every 500th frequency at the call site) give an
!	  aligned radial grid across renders for the V12-vs-KNOP field
!	  comparison; separate cap so they cannot displace trigger dumps.
	  IF(NSCND .GE. 150)RETURN
	  NSCND=NSCND+1
	ELSE IF(ETYPE .EQ. 'NRESC')THEN
	  IF(NDMP .GE. 200)RETURN
	  NDMP=NDMP+1
	ELSE
	  RETURN
	END IF
	IF(LUD .EQ. 0)THEN
	  OPEN(NEWUNIT=LUD,FILE='J_RESCUE_DUMP',STATUS='UNKNOWN',
	1    ACTION='WRITE')
	  WRITE(LUD,'(A)')'! JDMP indx fl sect bad nd then ND rows:'
	1    //' i R V SIG SIGT rj jprev fedd feddprev eta chi zeta'
	END IF
	IF(ETYPE .EQ. 'SCAN')THEN
	  WRITE(LUD,'(A,1X,I6,1X,ES16.8,1X,A,1X,I4,1X,I4)')
	1    'JSCN',FREQ_INDX,FL,TRIM(SECTION),BAD,ND
	ELSE
	  WRITE(LUD,'(A,1X,I6,1X,ES16.8,1X,A,1X,I4,1X,I4)')
	1    'JDMP',FREQ_INDX,FL,TRIM(SECTION),BAD,ND
	END IF
	DO I=1,ND
	  WRITE(LUD,'(A,1X,I4,12(1X,ES13.5))')'JDJ',I,R(I),V(I),
	1    SIG(I),SIGT(I),RJ(I),JP(I),FEDD(I),FP(I),ETA(I),CHI(I),
	2    ZETA(I)
	END DO
	FLUSH(LUD)
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
	  IF(CONT_SEALED)RETURN
	  CONT_SEALED=.TRUE.
	  WRITE(LUE,'(A,1X,A,2(1X,A,I8))')'JTOT',SECTION,
	1    'nresc=',NC_CNT,'nes_san=',NE_CNT
	ELSE
	  IF(LINE_SEALED)RETURN
	  LINE_SEALED=.TRUE.
	  WRITE(LUE,'(A,1X,A,2(1X,A,I8))')'JTOT',SECTION,
	1    'nresc=',NL_CNT,'nes_san=',NEL_CNT
	END IF
	FLUSH(LUE)
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
	FLUSH(LUE)
	RETURN
	END
!
! JRESC_FRQ: one machine-readable J_RESCUE_FREQ record per
! continuum frequency (DO_J_RES_DIAG, default F).  The JEV file
! only covers trigger frequencies; the V12-vs-KNOP comparison
! needs the field state at EVERY frequency, so this separate file
! keeps the full-sweep record without flooding JEV.  ~NCF lines.
!
! Record layout (JFRQ):
!   indx ml fl sect fg fdnc rjmax jdep rjnd fedmax gedmax
!   trig ftrig fired restored prevj prevf sig sigt eta chi zeta
! rjmax/jdep: max|RJ| and its depth index post-restore (or post-
! bypass under J_RESCUE_RESTORE=F - the record states which).
! rjnd: |RJ(ND)| - the inner-boundary value where the Phase 2.3
! inversion artifact lives.  trig/ftrig: this frequency's rescue
! triggers.  fired: the trigger armed an intervention (J_RESCUE
! block entered).  restored: the *_PREV copies actually ran -
! F under J_RESCUE_RESTORE=F or when fired=F.  prevj/prevf:
! max|JNU_PREV|/max|FEDD_PREV| snapshots (0 on FIRST_FREQ).
! sig/sigt/eta/chi/zeta: at depth jdep, the local source scales.
!
	SUBROUTINE JRESC_FRQ(FREQ_INDX,ML,FL,SECTION,FG_NM,FDNC,
	1    RJ,ND,JDEP,RJND,FEDMAX,GEDMAX,TRIG,FTRIG,FIRED,RESTOR,
	2    PREVJ,PREVF,SIGB,SIGTB,ETAB,CHIB,ZETAB,LUER)
	USE SET_KIND_MODULE
	USE JRES_DIAG_STATE
	IMPLICIT NONE
	CHARACTER(LEN=*) SECTION,FG_NM
	INTEGER FREQ_INDX,ML,ND,JDEP,TRIG,FTRIG,LUER
	REAL(KIND=LDP) FL,RJ(ND),RJND,FEDMAX,GEDMAX
	REAL(KIND=LDP) RJMAX,PREVJ,PREVF,SIGB,SIGTB,ETAB,CHIB,ZETAB
	LOGICAL FDNC,FIRED,RESTOR
	RJMAX=MAXVAL(ABS(RJ(1:ND)))
	IF(LUFQ .EQ. 0)THEN
	  OPEN(NEWUNIT=LUFQ,FILE='J_RESCUE_FREQ',STATUS='UNKNOWN',
	1    ACTION='WRITE')
	  WRITE(LUFQ,'(A)')'! JFRQ indx ml fl sect fg fdnc rjmax jdep'//
	1    ' rjnd fedmax gedmax trig ftrig fired restored prevj'//
	2    ' prevf sig sigt eta chi zeta'
	END IF
	WRITE(LUFQ,'(A,1X,I6,1X,I6,1X,ES16.8,1X,A,1X,A,1X,L1,'//
	1   '1X,ES12.4,1X,I4,3(1X,ES12.4),2(1X,I2),2(1X,L1),'//
	2   '7(1X,ES12.4))')
	3   'JFRQ',FREQ_INDX,ML,FL,TRIM(SECTION),TRIM(FG_NM),FDNC,
	4   RJMAX,JDEP,RJND,FEDMAX,GEDMAX,TRIG,FTRIG,FIRED,RESTOR,
	5   PREVJ,PREVF,SIGB,SIGTB,ETAB,CHIB,ZETAB
	RETURN
	END
