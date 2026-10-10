!
! Subroutine to plot the Luminosity and auxiliary quanties from the file OBSFLUX.
! Routine is primarily designed to assist in diagnosing SN models.
!
! The final plot shows the "Corrected Luminosity" (not OBSERVED) which should be
! constant as a function of depth. Also shown is the difference between this curve,
! and the original "observed" luminosity (CMF).
!
! Unlike OBSFLUX, we intgerate inwards rather than outwards. This allows a more
! realistic examination of the errors relative to the observed flux.
!
! Program requires the following CMFGEN files:
!                                              OBSFLUX
!                                              RVTJ
!                                              MODEL
!
	PROGRAM PLT_CMF_LUM
	USE SET_KIND_MODULE
        USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 30-Aug-2023 -- Added LOC_SHOCK and cleaned up curve labeling.
!                            Zero curves now omitted from plots.
! Altered: 08-Ayg-2018 -- Fixed offsef in RAD_DECAY, ADI, DJDT, MECH -- single point.
!                            Values are zero at outer boundary.
! Altered: 19-Apr-2016 -- Added plot of dT/T and conserved L on same plot.
!                            Large dT/T (i.e., > 0.1) are a major reason
!                            why the "conserved" luminosity is not conserved.
! Altered: 19-Mar-2014 -- Inserted HJ_LUM to handle static relativistic models.
! Cleaned: 06-Nov-2011
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: POPS(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: XAXIS(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: LUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: HJ_LUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: MECH(:)
	REAL(KIND=LDP), ALLOCATABLE :: LOC_SHOCK(:)
	REAL(KIND=LDP), ALLOCATABLE :: ADI(:)
	REAL(KIND=LDP), ALLOCATABLE :: DJDT(:)
	REAL(KIND=LDP), ALLOCATABLE :: dT_ON_T(:)
	REAL(KIND=LDP), ALLOCATABLE :: RAD_DECAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: TOTAL(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHANGE(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_VEC(:)
!
	REAL(KIND=LDP) SPEED_OF_LIGHT
	REAL(KIND=LDP) T1
!
	EXTERNAL SPEED_OF_LIGHT
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=20)  XLABEL
!
	INTEGER, PARAMETER :: LU_RD=20
        INTEGER IREC,NITSF,RITE_N_TIMES,LAST_NG
	LOGICAL WRITE_RVSIG,NEWMOD
	LOGICAL PLT_AGAINST_DEPTH_INDX
	LOGICAL PLT_AGAINST_T
	LOGICAL NON_ZERO_HJ_LUM
	LOGICAL NON_ZERO_MECH
	LOGICAL NON_ZERO_ADI
	LOGICAL NON_ZERO_DJDt
	LOGICAL NON_ZERO_LOC_SHOCK
	LOGICAL NON_ZERO_RAD_DECAY
!
	INTEGER ND
	INTEGER NT
	INTEGER IOS
	INTEGER I
	LOGICAL FILE_OPEN
!
	WRITE(6,'(A)')' '
        OPEN(UNIT=LU_RD,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
          IF(IOS .EQ. 0)THEN
            DO WHILE(1 .EQ. 1)
              READ(LU_RD,'(A)',IOSTAT=IOS)STRING
              IF(IOS .NE. 0)EXIT
              IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
                READ(STRING,*)ND
                WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
              END IF
              IF(INDEX(STRING,'!Total number of variables') .NE. 0)THEN
                READ(STRING,*)NT
                WRITE(6,'(A,I4)')'    Number of variables in the model is:',NT
                EXIT
              END IF
            END DO
          END IF
          INQUIRE(UNIT=LU_RD,OPENED=FILE_OPEN)
        IF(FILE_OPEN)CLOSE(UNIT=LU_RD)
!
        IF(IOS .NE. 0)THEN
          WRITE(6,*)' Unable to open MODEL file to get # of depth points'
          CALL GEN_IN(ND,'Number of depth points')
          CALL GEN_IN(NT,'Number of variables (NT)')
        END IF
!
	IOS=0
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (XAXIS(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE (LUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (HJ_LUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (MECH(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LOC_SHOCK(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DJDT(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (dT_ON_T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ADI(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAD_DECAY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TOTAL(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CHANGE(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (WRK_VEC(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error -- unable to allocate vectors in PLT_CMF_LUM'
	  WRITE(6,*)'Error is ',IOS
	  STOP
	END IF
!
	LUM=0.0_LDP; HJ_LUM=0.0_LDP; MECH=0.0_LDP; DJDT=0.0_LDP; ADI=0.0_LDP; LOC_SHOCK=0.0_LDP
	RAD_DECAY=0.0_LDP; TOTAL=0.0_LDP; CHANGE=0.0_LDP; T=0.0_LDP; WRK_VEC=0.0_LDP
!
	OPEN(UNIT=LU_RD,FILE='RVTJ',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    STRING=' '
	    DO WHILE(INDEX(STRING,'Radius') .EQ. 0)
	      READ(LU_RD,'(A)')STRING
	    END DO
	    READ(LU_RD,*)(R(I),I=1,ND)
	    READ(LU_RD,'(A)')STRING
	    READ(LU_RD,*)(V(I),I=1,ND)
	    DO WHILE(INDEX(STRING,'Temperature') .EQ. 0)
	      READ(LU_RD,'(A)')STRING
	    END DO
	    READ(LU_RD,*)(T(I),I=1,ND)
	  END IF
	CLOSE(UNIT=LU_RD)
	IF(IOS .NE. 0)THEN
          ALLOCATE (POPS(NT,ND))
	  ALLOCATE (SIGMA(ND))
	  IREC=0                  ! Get last iteration
          CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LAST_NG,
	1                 WRITE_RVSIG,NT,ND,LU_RD,NEWMOD)
	  IF(NEWMOD)THEN
	    WRITE(6,*)'Unable to read R, V, etc from SCRTEMP'
	    STOP
	  END IF
	  T(1:ND)=POPS(NT,1:ND)
	END IF
!
	DO I=1,ND
	  R(I)=R(I)/R(ND)
	END DO
	V(1:ND)=0.001_LDP*V(1:ND)
	XAXIS(1:ND)=V(1:ND)
	XLABEL='V(Mm\u \ds\u-1\d)'
	PLT_AGAINST_DEPTH_INDX=.FALSE.
        CALL GEN_IN(PLT_AGAINST_DEPTH_INDX,'Plot against depth index instead of V?')
	IF(PLT_AGAINST_DEPTH_INDX)THEN
	  DO I=1,ND
	    XAXIS(I)=I
	  END DO
	  XLABEL='Depth index'
	END IF
!
	IF(T(1) .NE. 0.0_LDP)THEN
	  DO I=2,ND
	    dT_ON_T(I)=(T(I)-T(I-1))/T(I-1)
	  END DO
	  dT_ON_T(1)=dT_ON_T(2)
	END IF
!
	OPEN(UNIT=20,FILE='OBSFLUX',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Luminosity') .EQ. 0)
	    READ(20,'(A)')STRING
	  END DO
	  READ(20,*)(LUM(I),I=1,ND)
!
	  NON_ZERO_MECH=.FALSE.; NON_ZERO_LOC_SHOCK=.FALSE.
	  NON_ZERO_ADI=.FALSE.; NON_ZERO_DJDT=.FALSE.
	  NON_ZERO_HJ_LUM=.FALSE.; NON_ZERO_RAD_DECAY=.FALSE.
	  DO WHILE(1 .EQ. 1)
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ')
	      READ(20,'(A)')STRING
	    END DO
	    IF(INDEX(STRING,'Mechanical Luminosity') .NE. 0)THEN
	      READ(20,*)(MECH(I),I=1,ND)
	      IF(MAXVAL(ABS(MECH)) .GT. 0.0_LDP)NON_ZERO_MECH=.TRUE.
	    ELSE IF(INDEX(STRING,'locally by the shock') .NE. 0)THEN
	      READ(20,*)(LOC_SHOCK(I),I=1,ND)
	      IF(MAXVAL(ABS(LOC_SHOCK)) .GT. 0.0_LDP)NON_ZERO_LOC_SHOCK=.TRUE.
	    ELSE IF(INDEX(STRING,'Internal') .NE. 0)THEN
	      READ(20,*)(ADI(I),I=1,ND)
	      IF(MAXVAL(ABS(ADI)) .GT. 0.0_LDP)NON_ZERO_ADI=.TRUE.
	    ELSE IF(INDEX(STRING,'Luminosity [g.r^2.H + beta.g.r^2.J]') .NE. 0)THEN
	      READ(20,*)(HJ_LUM(I),I=1,ND)
	      IF(MAXVAL(ABS(HJ_LUM)) .GT. 0.0_LDP)NON_ZERO_HJ_LUM=.TRUE.
	    ELSE IF(INDEX(STRING,'Flux arrising') .NE. 0)THEN
	      READ(20,*)(DJDt(I),I=1,ND)
	      IF(MAXVAL(ABS(DJDt)) .GT. 0.0_LDP)NON_ZERO_DJDt=.TRUE.
	    ELSE IF (INDEX(STRING,'Energy deposited') .NE. 0)THEN
	      READ(20,*)(RAD_DECAY(I),I=1,ND)
	      IF(MAXVAL(ABS(RAD_DECAY)) .GT. 0.0_LDP)NON_ZERO_RAD_DECAY=.TRUE.
	    ELSE IF (INDEX(STRING,'Total Radiative') .NE. 0)THEN
	      EXIT
	    END IF
	    READ(20,'(A)')STRING
	  END DO
	  CLOSE(UNIT=20)
!
	WRITE(6,*)' '
	WRITE(6,'(3A)')PG_PEN(2),'          The CMF luminosity is plotted in red',DEF_PEN
	WRITE(6,'(3A)')PG_PEN(3),' The "conserved" luminosity is plotted in blue',DEF_PEN
	WRITE(6,*)' '
	T1=0.0_LDP
	TOTAL=LUM
	IF(NON_ZERO_HJ_LUM)TOTAL=HJ_LUM
	DO I=1,ND-1
	  T1=T1+LOC_SHOCK(I)+RAD_DECAY(I)-MECH(I)-ADI(I)-DJDT(I)
	  TOTAL(I+1)=TOTAL(I+1)+T1
	END DO
	CALL DP_CURVE_LAB(ND,XAXIS,LUM,'L(cmf)')
	CALL DP_CURVE_LAB(ND,XAXIS,TOTAL,'L(cons)')
	CALL GRAMON_PGPLOT(XLABEL,'Luminosity (L\d'//char(09)//'\u)',' ',' ')
!
	WRITE(6,*)' '
	WRITE(6,'(3A)')PG_PEN(2),'                        dT/T is plotted in red',DEF_PEN
	WRITE(6,'(3A)')PG_PEN(3),' The "conserved" luminosity is plotted in blue',DEF_PEN
	WRITE(6,'(3A)')PG_PEN(3),' The "conserved" luminosity has been scaled to 0.1 at the outer boundary',DEF_PEN
	WRITE(6,*)' '
	WRK_VEC=0.1_LDP*TOTAL/TOTAL(2)
	CALL DP_CURVE_LAB(ND,XAXIS,dT_ON_T,'dT/T')
	CALL DP_CURVE_LAB(ND,XAXIS,WRK_VEC,'L(cons)')
	CALL GRAMON_PGPLOT(XLABEL,'Scaled cons. L; dT/T',' ',' ')
	
	WRITE(6,*)' '
	WRITE(6,*)' Plotting the corrections at each depth'
	WRITE(6,*)' '
	I=1
	IF(NON_ZERO_RAD_DECAY)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,RAD_DECAY,'Rad. Dec.'); I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//' The "radioactive decay" contribution at each depth.',DEF_PEN
	END IF
	IF(NON_ZERO_LOC_SHOCK)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,LOC_SHOCK,'Shock');  I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//' The "local shock" contribution at each depth.',DEF_PEN
	END IF
	IF(NON_ZERO_MECH)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,MECH,'Mechanical');  I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//' The "mechanical" luminosity contribution at each depth.',DEF_PEN
	END IF
	IF(NON_ZERO_ADI)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,ADI,'Adibat.');    I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//' The "adiabatic" luminosity contribution at each depth.',DEF_PEN
	END IF
	IF(NON_ZERO_ADI)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,DJDT,'DJDt');   I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//' The "DJDt" luminosity contribution at each depth. ',DEF_PEN
	END IF
	WRITE(6,*)' '
	CALL GRAMON_PGPLOT(XLABEL,'\gDL(L\d'//char(09)//'\u)',' ',' ')

	LOC_SHOCK(2:ND)=LOC_SHOCK(1:ND-1); LOC_SHOCK(1)=0.0_LDP
	MECH(2:ND)=MECH(1:ND-1); MECH(1)=0.0_LDP
	DJDT(2:ND)=DJDT(1:ND-1); DJDT(1)=0.0_LDP
	ADI(2:ND)=ADI(1:ND-1); ADI(1)=0.0_LDP
	RAD_DECAY(2:ND)=RAD_DECAY(1:ND); RAD_DECAY(1)=0.0_LDP
	DO I=2,ND
	  LOC_SHOCK(I)=LOC_SHOCK(I)+LOC_SHOCK(I-1)
	  MECH(I)=MECH(I)+MECH(I-1)
	  DJDT(I)=DJDT(I)+DJDT(I-1)
	  ADI(I)=ADI(I)+ADI(I-1)
	  RAD_DECAY(I)=RAD_DECAY(I)+RAD_DECAY(I-1)
	END DO
	MECH=-MECH
	DJDT=-DJDT
	ADI=-ADI
!
	WRITE(6,*)
	WRITE(6,*)' '
	CALL WR_COL_STR(' |1*** |0Plotting cummalative contributions (intgerated inwards). |1***')
	CALL WR_COL_STR(' |1*** |0These represent the contributions to the "conserved luminosity)". |1***')
	WRITE(6,*)' '
!
	I=2
	WRITE(6,'(3A)')PG_PEN(2)//'            The "corrected luminosity" is in red',DEF_PEN
	CALL DP_CURVE_LAB(ND,XAXIS,TOTAL,'Cor. Lum.')
!
	IF(NON_ZERO_RAD_DECAY)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,RAD_DECAY,'Rad. decay'); I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//'    The term due to radioactive decay',DEF_PEN
	END IF
!
	IF(NON_ZERO_LOC_SHOCK)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,LOC_SHOCK,'Shock');  I=I+1
	WRITE(6,'(3A)')PG_PEN(I)//'  Local shock energy is in yellow',DEF_PEN
	END IF
!
	IF(NON_ZERO_DJDt)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,DJDT,'DJDt');   I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//'                       The DJDT term is in green',DEF_PEN
	END IF
!
	IF(NON_ZERO_MECH)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,MECH,'Mech.');    I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//'      The work on the gas term (MECH) is in purple',DEF_PEN
	END IF
!
	IF(NON_ZERO_ADI)THEN
	  CALL DP_CURVE_LAB(ND,XAXIS,ADI,'Adiab.');  I=I+1
	  WRITE(6,'(3A)')PG_PEN(I)//'  Adiabatic cooling/internal energy is in pink',DEF_PEN
	END IF
!
	I=I+1
	WRITE(6,'(3A)')PG_PEN(I)//' The total of the corrections terms in in yellow',DEF_PEN
	CHANGE=MECH+ADI+DJDT+RAD_DECAY+LOC_SHOCK
	CALL DP_CURVE_LAB(ND,XAXIS,CHANGE,'Total Cor.')
!
	WRITE(6,*)' '
	WRITE(6,*)'MECH',MECH(ND)
	WRITE(6,*)'ADI',ADI(ND)
	WRITE(6,*)'DJDT',DJDT(ND)
	WRITE(6,*)'RAD_DECAY',RAD_DECAY(ND)
	WRITE(6,*)'LOC_SHOCK',LOC_SHOCK(ND)
!
	CALL GRAMON_PGPLOT(XLABEL,'Luminosity (L\d\m9\u)',' ',' ')
!
	WRITE(6,*)' '
	WRITE(6,'(3A)')PG_PEN(2)//'Here we plot the total energy balance (normalized by the value at d=2).'
	WRITE(6,'(3A)')PG_PEN(2)//'I recommend subtracting 1 from plot 1, and multiply by 10 or 100.',DEF_PEN
	WRITE(6,'(3A)')PG_PEN(3)//'We also plot the normalized correction.',DEF_PEN
	WRITE(6,*)' '
!
	T1=TOTAL(2)
	TOTAL=TOTAL/T1
	CALL DP_CURVE_LAB(ND,XAXIS,TOTAL,'Cor. Lum.')
	CHANGE=CHANGE/T1
	CALL DP_CURVE_LAB(ND,XAXIS,CHANGE,'Total Cor.')
	CALL GRAMON_PGPLOT(XLABEL,' ',' ',' ')
!
	STOP
	END
