	SUBROUTINE GET_PROFILE_STORAGE_LIMITS_V2(
	1       NPROF_MAX,NFREQ_MAX,NPROF_MAX_VOIGT,NFREQ_MAX_VOIGT,MAX_SIM,
	1       LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,
	1       PROF_TYPE,VEC_TRANS_TYPE,NLINES,NCF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 19-Mar-2023  Now initialize MAX_SIM, NPROF_CNT_VOIGT, and NFREQ_AV_VOIGT
! Altered: 05-Jul-2022 	Call changed to V2
!			VOIGT variables and MAX_SIM added to call.
!			Extra diagnostic output to LUER.
!			
! These variables are returned.
!
	INTEGER NFREQ_MAX		!Maximum number of frequencies for  profile.
	INTEGER NPROF_MAX		!Maximum number of Stark profiles to be stored at any one time
	INTEGER NPROF_MAX_VOIGT		!Maximum number of Voigt profiles to be stored at any one time
	INTEGER NFREQ_MAX_VOIGT		!Maximum number of frequencies for Voigt profile.
	INTEGER MAX_SIM			!Maximum number of profiles (any type except SOB) treated at the same time.
!
! Addtional diagnostic variables (not returned).
!
	INTEGER NFREQ_AV_VOIGT		!Average number of frequencies for Voigt  profile.
	INTEGER NFREQ_AV		!Average number of frequencies per profile.
	INTEGER NPROF_CNT
	INTEGER NPROF_CNT_VOIGT
!
! Required as input:
!
	INTEGER NLINES
	INTEGER NCF
	INTEGER LINE_ST_INDX_IN_NU(NLINES)
	INTEGER LINE_END_INDX_IN_NU(NLINES)
	CHARACTER(LEN=*) VEC_TRANS_TYPE(NLINES)
	CHARACTER(LEN=*) PROF_TYPE(NLINES)
!
! Local variables
!
	INTEGER K
	INTEGER NL, NL_SAVE
	INTEGER ML
	INTEGER LN_INDX
	INTEGER NX
	INTEGER NV
	INTEGER NSIM
!
! Get maximum number of frequencies required by Voigt and Stark intrinsic line profiles.
! Potentially used by other routines to set profile storage.
!
	MAX_SIM=0
	NFREQ_MAX=0;  NFREQ_MAX_VOIGT=0
	NPROF_CNT=0;  NPROF_CNT_VOIGT=0
	NPROF_CNT_VOIGT=0; NFREQ_AV_VOIGT=0
	DO NL=1,NLINES
	  ML=LINE_END_INDX_IN_NU(NL)-LINE_ST_INDX_IN_NU(NL)+1
	  IF(VEC_TRANS_TYPE(NL) .EQ. 'SOB')THEN
	  ELSE IF(PROF_TYPE(NL) .EQ. 'VOIGT')THEN
	    NFREQ_MAX_VOIGT=MAX(NFREQ_MAX_VOIGT,ML)
	    NFREQ_AV_VOIGT=NFREQ_AV_VOIGT+ML
	    NPROF_CNT_VOIGT=NPROF_CNT_VOIGT+1
	  ELSE IF(PROF_TYPE(NL)(1:3) .NE. 'DOP')THEN
	    NFREQ_MAX=MAX(NFREQ_MAX,ML)
	    NFREQ_AV=NFREQ_AV+ML
	    NPROF_CNT=NPROF_CNT+1
	  END IF
	END DO
!
! Get the maximum number of lines that overlap, and that require storage of their intrinsic profile.
! Line computed with Doppler or Voigt profiles do not require storage.
!
	NPROF_MAX=0;  NPROF_MAX_VOIGT=0
	LN_INDX=1
	DO ML=1,NCF
	  NX=0; NV=0; NSIM=0
	  K=LN_INDX
	  DO NL=K,NLINES
	    NL_SAVE=NL
	    IF(VEC_TRANS_TYPE(NL) .EQ. 'SOB')THEN
	    ELSE IF(ML .GE. LINE_ST_INDX_IN_NU(NL) .AND.
	1           ML .LE. LINE_END_INDX_IN_NU(NL) )THEN
	      NSIM=NSIM+1
	      IF(PROF_TYPE(NL) .EQ. 'VOIGT')THEN
	        NV=NV+1
	      ELSE IF(PROF_TYPE(NL)(1:3) .NE. 'DOP')THEN
	        NX=NX+1
	      END IF
	    END IF
	    IF(LINE_END_INDX_IN_NU(NL) .LE. ML .AND. NL .EQ. LN_INDX)LN_INDX=LN_INDX+1
	    IF(LINE_ST_INDX_IN_NU(NL) .GT. ML)EXIT
	  END DO
	  NPROF_MAX=MAX(NX,NPROF_MAX)
	  NPROF_MAX_VOIGT=MAX(NV,NPROF_MAX_VOIGT)
	  MAX_SIM=MAX(MAX_SIM,NSIM)
	  IF(LN_INDX .EQ. NLINES)EXIT
	END DO
!
	IF(NPROF_MAX .NE. 0 .OR. NFREQ_MAX .NE. 0)THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'Total number of stark profiles is           ',NPROF_CNT
	  WRITE(6,*)'Maximum number of profiles to be stored is  ',NPROF_MAX
	  WRITE(6,*)'Maximum number of frequencies per profile is',NFREQ_MAX
	  WRITE(6,*)'Average number of frequencies per profile is',NFREQ_AV/NPROF_CNT
	END IF
	IF(NPROF_MAX_VOIGT .NE. 0 .OR. NFREQ_MAX_VOIGT .NE. 0)THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'Total number of Voigt profiles is                 ',NPROF_CNT_VOIGT
	  WRITE(6,*)'Maximum number of Voigt profiles to be stored is  ',NPROF_MAX_VOIGT
	  WRITE(6,*)'Maximum number of Voigt frequencies per profile is',NFREQ_MAX_VOIGT
	  WRITE(6,*)'Average number of frequncies/Voigt profile is     ',NFREQ_AV_VOIGT/NPROF_CNT_VOIGT
	  WRITE(6,*)' '
	END IF
	WRITE(6,*)'Maximum number of overlapping profiles is ',MAX_SIM
	WRITE(6,*)' '
!
	RETURN
	END
