!
! Program to create a NEW_R_GRID which is equally spaced in LOG(Tau) where
! TAU is based on the FLUX mean opacity.
!
	SUBROUTINE ADJUST_R_GRID_V4(POPS,ESEC,MAIN_COUNTER,DONE_R_REV,ND,NT)
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD, ONLY: DO_NONMONO_SOB,DO_NONMONO_KNOP
	USE UPDATE_KEYWORD_INTERFACE
	IMPLICIT NONE
!
! Altered 05-Oct-2026 : VEL_SPACE rejected under DO_NONMONO_SOB /
!                       DO_NONMONO_KNOP -- its R(V) interpolation
!                       assumes a monotonic velocity field (Devin).
! Aleterd 20-Aug-2016 : Introduced OB_ONLY option
! Altered 19-Aug-2015 : Check on makig sure SIGMA > -1 (cur_hmi, 23-Jun-2105)
! Altered 02-Dec-2012 : Changed calls to DO_TAU_REGRID and DO_VEL_REGRID.
!                       Now open R_REGRIDDING_LOG in this routine.
!                       Use monotonic interpolaton for V and SIGMA.
!
! Altered 09-Feb-2011 : Major rewrite of ADJUST_R_GRID_V3 which was  a version
!                       still under development. This routine is controlled by
!	                parameters read in from the file ADJUST_R_DEFAULTS. This
!                       will give greater flexibility, and allow the routine to
!                       be asily modified. Call is diferent from V2.
!
	INTEGER ND,NT
	INTEGER MAIN_COUNTER
!
	REAL(KIND=LDP) POPS(NT,ND)
	REAL(KIND=LDP) ESEC(ND)			!Electron scattering opacity
!
	LOGICAL DONE_R_REV
!
! For specifying grid.
!
	CHARACTER(LEN=10) REGRIDDING_METHOD
	CHARACTER(LEN=10) GRID_TYPE
	CHARACTER(LEN=10) OUT_BND_OPT			!Outer boundary option
	CHARACTER(LEN=10) IN_BND_OPT			!Inner boundary option
!
! Local variables.
!
	REAL(KIND=LDP) R_OLD(ND)
	REAL(KIND=LDP) V_OLD(ND)
	REAL(KIND=LDP) SIGMA_OLD(ND)
!
	REAL(KIND=LDP) LOG_R_OLD(ND)
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) dTAU_OLD(ND)
	REAL(KIND=LDP) TAU_OLD(ND)
	REAL(KIND=LDP) TAU(ND)
!
	REAL(KIND=LDP) TA(ND)			!Work vectors
	REAL(KIND=LDP) TB(ND)
	REAL(KIND=LDP) COEF(ND,4)
!
! The fine grid (FG) is chosen to cover the ionization front. The default values are
! -2.0 to 1.0D0 in log(TAU) space.
!
	REAL(KIND=LDP) FG_MIN			!Min Tau for FG
	REAL(KIND=LDP) FG_MAX			!Max Tau for FG
	REAL(KIND=LDP) FG_RANGE
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) DTAU2_ON_DTAU1
	REAL(KIND=LDP) LOG_TAU
	REAL(KIND=LDP) STRETCH_POW		!Power law exponent to stretch tau scale about 1
!
	REAL(KIND=LDP) OBND_PARAMS(5)		!Parameters specifying grid placement at outer boundary.
	REAL(KIND=LDP) IBND_PARAMS(5)		!Parameters specifying grid placement at nner boundary.
!
	INTEGER NUM_IBND_PARAMS		!Number of points inserted near inner boundary.
	INTEGER NUM_OBND_PARAMS		!Number of points inserted near outer boundary.
!
! The parameters describe on what iteration we should do the first grid adjustement,
! how often we should do  a grid adjustment, and how many grid adjustemnts should  be
! done. NO_R_REV is adjusted downwards after each iteration, and is output to ADJUST_R_DEFAULTS.
!
	INTEGER STRT_R_REV
	INTEGER FREQ_R_REV
	INTEGER NO_R_REV
!
	INTEGER NX
	INTEGER I,I1,I2,J
	INTEGER IOS
	INTEGER LU
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUSCR=8
	LOGICAL ERROR
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL ONLY_OB_DONE
	CHARACTER(LEN=80) STRING
!
	WRITE(T_OUT,*)'Entering ADJUST_R_GRID_V4'
!
! Set default parameters
!
        DONE_R_REV=.FALSE.
	STRT_R_REV=1; FREQ_R_REV=10
	OUT_BND_OPT='DEFAULT'
	IN_BND_OPT='DEFAULT'
!
! Read in parameters describing the new model.
!
	CALL GEN_ASCI_OPEN(LUIN,'ADJUST_R_DEFAULTS','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening ADJUST_R_DEFAULTS in ADJUST_R_GRID_V4, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
	CALL RD_STORE_INT(NO_R_REV,'N_ITS',L_TRUE,'Number of iterations remaining')
	CALL RD_STORE_INT(STRT_R_REV,'STRT_ITS',L_FALSE,'Iteration to start first R-grid revision')
	CALL RD_STORE_INT(FREQ_R_REV,'FREQ_ITS',L_FALSE,'Frequency for R-grid revisions')
	I=10; CALL RD_STORE_NCHAR(REGRIDDING_METHOD,'GRID_METH',I,L_TRUE,
	1            'Regridding method: TAU_SPACE, VEL_SPACE, FULL_R')
!
	CLOSE(UNIT=LUIN)
	CLOSE(UNIT=LUSCR)
!
! Decide here whether we will do an iteration or not.
!
	WRITE(T_OUT,'(A,I4)')'         Number of R revisons:',NO_R_REV
	WRITE(T_OUT,'(A,I4)')'                 Main counter:',MAIN_COUNTER
	WRITE(T_OUT,'(A,I4)')' Iteration to start revisions:',STRT_R_REV
	WRITE(T_OUT,'(A,I4)')'       Fequency of iterations:',FREQ_R_REV
        IF(NO_R_REV .EQ. 0 .OR. MAIN_COUNTER .LT. STRT_R_REV .OR.
	1  MOD( (MAIN_COUNTER-STRT_R_REV),FREQ_R_REV ) .NE. 0)THEN
	   WRITE(T_OUT,'(A,I4)')' No R revision required on this iteration'
	  CALL CLEAN_RD_STORE()
	  RETURN
	END IF
!
        CALL GET_LU(LU)
	OPEN(UNIT=LU,FILE='R_REGRIDDING_LOG',STATUS='UNKNOWN',ACTION='WRITE')
	CALL SET_LINE_BUFFERING(LU)		!Switches off line buffering
!
! Save existing grid, which will be used for the interplations.
!
	R_OLD(1:ND)=R(1:ND)
!
	ONLY_OB_DONE=.FALSE.
	IF(REGRIDDING_METHOD .EQ. 'TAU_SPACE')THEN
	  CALL DO_TAU_REGRID_V2(POPS,ESEC,DONE_R_REV,ND,NT,LU)
	ELSE IF(REGRIDDING_METHOD .EQ. 'VEL_SPACE')THEN
!
! VEL_SPACE interpolates R as a function of V (LIN_INTERP on the
! R(V) table in DO_VEL_REGRID_V2) and brackets the regridded
! interval by VHIGH/VLOW -- both assume V decreases monotonically.
! On a non-monotonic field (shock / velocity inversion) the R(V)
! map is not single-valued and the grid is corrupted silently.
! TAU_SPACE, FULL_R and OB_ONLY operate on monotone R/tau
! coordinates and are unaffected.  [DO_NONMONO_SOB,
! DO_NONMONO_KNOP]; default path unchanged.
!
	  IF(DO_NONMONO_SOB .OR. DO_NONMONO_KNOP)THEN
	    WRITE(T_OUT,'(A)')'Error in ADJUST_R_GRID_V4:'
	    WRITE(T_OUT,'(A)')'GRID_METH=VEL_SPACE is incompatible with'
	    WRITE(T_OUT,'(A)')'DO_NONMONO_SOB/DO_NONMONO_KNOP: the R(V)'
	    WRITE(T_OUT,'(A)')'interpolation in DO_VEL_REGRID_V2 assumes'
	    WRITE(T_OUT,'(A)')'a monotonic velocity field.'
	    WRITE(T_OUT,'(A)')'Use TAU_SPACE, FULL_R or OB_ONLY'
	    WRITE(T_OUT,'(A)')'(or set REV_RGRID=F).'
	    STOP
	  END IF
	  CALL DO_VEL_REGRID_V2(POPS,R,V,POP_ATOM,DONE_R_REV,ND,NT,LU)
	ELSE IF(REGRIDDING_METHOD .EQ. 'FULL_R')THEN
	  CALL DO_FULL_R_GRID_V1(R_OLD,DONE_R_REV,LU,ND)
	ELSE IF(REGRIDDING_METHOD .EQ. 'OB_ONLY')THEN
	  T1=0.5_LDP*(ROSS_MEAN(1)*CLUMP_FAC(1)+ROSS_MEAN(2)*CLUMP_FAC(2))
	  T2=0.5_LDP*(ROSS_MEAN(2)*CLUMP_FAC(2)+ROSS_MEAN(3)*CLUMP_FAC(3))
	  CALL RD_STORE_DBLE(DTAU2_ON_DTAU1,'D2OND1',L_FALSE,'~DTAU(2)/DTAU(1) at outer boudary')
	  R(2)=(T2*R(3)+DTAU2_ON_DTAU1*T1*R(1))/(T2+T1*DTAU2_ON_DTAU1)
	  IF(R(2) .GE. R(1) .OR. R(2) .LE. R(3))THEN
	    WRITE(T_OUT,'(A)')' Error -- invalid R(2) computed with OB_ONLY'
	    WRITE(T_OUT,'(A)')' Skipping outer boundary adjustment'
	  ELSE
	    ONLY_OB_DONE=.TRUE.
	  END IF
	ELSE
	  WRITE(T_OUT,'(A)')'Error- REGRIDDING_METHOD (GRID_METH option) not recognized in ADJUST_R_GRID_V4 '
	  WRITE(T_OUT,'(A)')'REGRIDDING_METHOD =',TRIM(REGRIDDING_METHOD)
	  CALL CLEAN_RD_STORE()
	  RETURN
	END IF
	CALL CLEAN_RD_STORE()
!
! Before refining the grid, we check to see if it is really necessary.
!
	T1=0.0_LDP
	DO I=2,ND-1
	  T1=MAX( T1,ABS(R(I)-R_OLD(I))/(R(I-1)-R(I+1)) )
	END DO
	T1=T1*2.0_LDP
	IF(T1 .LT. 1.0E-03_LDP .AND. .NOT. ONLY_OB_DONE)THEN
	  R=R_OLD
	  WRITE(T_OUT,*)'As grid is identical to 1 part in 1000, interpolation not necessary'
	  RETURN
	ELSE

! We now need to regrid all the populations. All interpolations (except
! sigma) are performed in the LOG-LOG plane. For SN this is ideal, since
! the density and velocity are power laws in r. For SIGMA, we do not take
! the log.
!
! We do not need to interpolate T, and ED directly, since these are part
! of POPS.
!
	  V_OLD(1:ND)=V(1:ND)
	  SIGMA_OLD(1:ND)=SIGMA(1:ND)
	  LOG_R=LOG(R)
	  LOG_R_OLD=LOG(R_OLD)
	  TA(1:ND)=LOG(V(1:ND))
!
	  CALL MON_INT_FUNS_V2(COEF,TA,LOG_R_OLD,ND)
	  J=1
	  ERROR=.FALSE.
	  DO I=2,ND-1
	    DO WHILE(R(I) .LT. R_OLD(J+1))
	      J=J+1
	    END DO
	    T1=LOG_R(I)-LOG_R_OLD(J)
	    V(I)=EXP( COEF(J,4)+T1*(COEF(J,3)+T1*(COEF(J,2)+T1*COEF(J,1))) )
	    SIGMA(I)=COEF(J,3)+T1*(2.0_LDP*COEF(J,2)+3.0_LDP*T1*COEF(J,1))
	    SIGMA(I)=SIGMA(I)-1.0_LDP
!
! edited by Devin, 13-Oct-2026: sigma < -1 on the interpolated grid is
! expected on non-monotonic velocity fields -- the reverse shock sits
! at V >> 0.2 km/s, so the upstream clamp/STOP below would abort a
! legitimate REV_RGRID revision.  Under [DO_NONMONO_SOB]/
! [DO_NONMONO_KNOP] the interpolated value is left as-is:
! RECOMPUTE_SIGMA_TRUE runs immediately after this routine returns
! (solve_for_pops.f) and reclamps SIGMA under SOB, while the KNOP
! formal solver handles sigma < -1 natively.  Default path unchanged.
!
	    IF(.NOT. (DO_NONMONO_SOB .OR. DO_NONMONO_KNOP))THEN
	      IF(SIGMA(I) .LT. -1.0_LDP .AND. V(I) .LT. 0.2E00_LDP)THEN
	        SIGMA(I)=-0.999_LDP
	      ELSE IF(SIGMA(I) .LT. -1.0_LDP)THEN
	        ERROR=.TRUE.
	      END IF
	    END IF
	  END DO
!
	  IF(ERROR)THEN
	    WRITE(6,*)'Error in ADJUST_R_GRID_V4 - SIGMA is not monotonic'
	    WRITE(6,*)'J, V(J) and SIGMA(J) follows'
	    DO J=1,I
	      WRITE(6,*)J,V(J),SIGMA(J)
	    END DO
	    STOP
	  END IF
!
	  WRITE(LU,'(A,8(7X,A))')'!Index','        R','     Rold','   Log(R)','Log(Rold)','        V',
	1                               '     Vold','    Sigma','Sigma_old'
	  WRITE(LU,'(A)')'!'
	  DO I=1,ND
	    WRITE(LU,'(I6,8ES16.5)')I,R(I),R_OLD(I),LOG_R(I),LOG_R_OLD(I),V(I),V_OLD(I),SIGMA(I),SIGMA_OLD(I)
	  END DO
	  CLOSE(LU)
!
! Now interpolate all populations, temperature etc.
!
	  DO I=1,NT
	    TA(1:ND)=LOG(POPS(I,1:ND))
	    CALL MON_INTERP(TB,ND,IONE,LOG_R,ND,TA,ND,LOG_R_OLD,ND)
	    POPS(I,1:ND)=EXP(TB(1:ND))
	  END DO
	END IF
!
	DONE_R_REV=.TRUE.
	NO_R_REV=NO_R_REV-1
        CALL UPDATE_KEYWORD(NO_R_REV,'[N_ITS]','ADJUST_R_DEFAULTS',L_TRUE,L_TRUE,LUIN)
	WRITE(T_OUT,*)'Adjusted R grid in ADJUST_R_GRID_V4'
!
	RETURN
	END
