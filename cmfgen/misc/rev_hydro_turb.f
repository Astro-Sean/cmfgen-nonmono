!
! Auxilary program designed to modify the HYDRO file output from CMFGEN.
! Progam modifies the adopted stellar mass. Not that the percentage error
! os now defined so that it has a rang of pm 200%.
!
	PROGRAM REV_HYDRO_FILE
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
!
! Altered: 05-Oct-2026 - Validate V ordering in INT and log argument in TST.
! Altered: 20-Sep-2020 - Now use MODEL_SPEC to get ND.
! Altered: 05-Jun-2020 - Error now defined using ABS(GRAD)+ABS(G) rather than ABS(G_TOT) which can be
!                          small at sonic point.
!                          Output additional column: 100E/G_RAD
! Altered: 23-Jun-2003 - ND and STRING now initialized.
! Cleaned: 07-Nov-2000
!
	IMPLICIT NONE
!
	INTEGER I,J,K,IOS
	INTEGER NSTR
	INTEGER ND
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: IZERO=0
!
	REAL(KIND=LDP) R
	REAL(KIND=LDP) RSQ
	REAL(KIND=LDP) V
	REAL(KIND=LDP) E,E_ON_GRAD
	REAL(KIND=LDP) VdVdR
	REAL(KIND=LDP) dPdR
	REAL(KIND=LDP) g_tot
	REAL(KIND=LDP) g_rad
	REAL(KIND=LDP) g_elec
	REAL(KIND=LDP) Gamma
	REAL(KIND=LDP) MT
	REAL(KIND=LDP) RND
	REAL(KIND=LDP) MASS_NEW
	REAL(KIND=LDP) MASS_OLD
	REAL(KIND=LDP) GSUR_NEW
	REAL(KIND=LDP) GSUR_OLD
	REAL(KIND=LDP) GRAV_CON
	REAL(KIND=LDP) RPHOT
!
	REAL(KIND=LDP) dTPdR
        REAL(KIND=LDP) VTURB
	REAL(KIND=LDP) T1,T2
!
! These variables are used to estimate a new mass for the STAR so that the
! Hydrostatic equation is better satified (in a least squares sense) over
! some range of depths.
!
        REAL(KIND=LDP) SUM_ERROR,SUM_R,DENOM,DEL_M
        INTEGER LOW_LIM,HIGH_LIM
!
	REAL(KIND=LDP) GRAVITATIONAL_CONSTANT, MASS_SUN
	EXTERNAL GRAVITATIONAL_CONSTANT, MASS_SUN
!
	CHARACTER*132 STRING(500)
	CHARACTER*132 TMP_STRING
	CHARACTER*132 FMT
	CHARACTER*132 FILENAME,RVTJ_FILE_NAME
	CHARACTER(LEN=80) XLAB,YLAB,TIT
	CHARACTER(LEN=10) XOPT
!
	REAL(KIND=LDP) P_R(200)
	REAL(KIND=LDP) P_VEL(200)
	REAL(KIND=LDP) P_dVdR(200)
	REAL(KIND=LDP) P_dPdR(200)
	REAL(KIND=LDP) P_REQ(200)
	REAL(KIND=LDP) P_GRAD(200)
	REAL(KIND=LDP) P_GELEC(200)
	REAL(KIND=LDP) P_GTOT(200)
	REAL(KIND=LDP) P_GRAV(200)
	REAL(KIND=LDP) ED(200)
	REAL(KIND=LDP) TA(200)
	REAL(KIND=LDP) XVEC(200)
	REAL(KIND=LDP) YVEC(200)
	REAL(KIND=LDP) ZVEC(200)
!
	LOGICAL OLD_FORMAT
	LOGICAL PLANE_PARALLEL
!
	GRAV_CON=1.0E-20_LDP*GRAVITATIONAL_CONSTANT()*MASS_SUN()
!
	ND=0
	FILENAME='MODEL_SPEC'
	OPEN(UNIT=20,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)
        IF(IOS .EQ. 0)THEN
	  READ(20,'(A)',IOSTAT=IOS)STRING(1)
	 IF(INDEX(STRING(1),'[ND]') .NE. 0)THEN
	   READ(STRING(1),*)ND
	   WRITE(6,*)'Number of depth points is',ND
	 ELSE
	   WRITE(6,*)'Invalid error reading MODEL_SPEC'
	 END IF
	 CLOSE(UNIT=20)
	ELSE
	   WRITE(6,*)'Unable to open MODEL_SPEC. IOS=',IOS
	END IF
	IF(ND .EQ. 0)THEN
	  WRITE(6,'(A)',ADVANCE='NO')'Input number od depth points in model: '
	  READ(5,*)ND
	END IF
!
	STRING(:)=' '
	FILENAME='HYDRO'
	CALL GEN_IN(FILENAME,'Input hydro file')
	OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD')
	  NSTR=0
	  DO WHILE(1 .EQ. 1)
	    READ(10,'(A)',END=1000)STRING(NSTR+1)
	    NSTR=NSTR+1
	  END DO
1000	  CONTINUE
	CLOSE(UNIT=10)
!
        I=160
	FILENAME='REV_HYDRO'
	CALL GEN_IN(FILENAME,'Output hydro file')
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'UNKNOWN',' ',' ',I,IOS)
	WRITE(LU_OUT,'(5X,A,15X,A,6X,A, 6(4X,A,1X), 6X,A, 6X,A, 3X, A)')
	1        'R','V','% Error','    VdVdR',
	1                   ' dPdR/ROH',
	1                   'dTPdR/ROH',
	1                   '    g_TOT',
	1                   '    g_RAD',
	1                    '   g_ELEC','Gamma','M(t)','E/G_RAD'
!
	MASS_OLD=0.0_LDP
	DO I=ND+1,NSTR
	  IF(INDEX(STRING(I),'urface gravity is:') .NE. 0)THEN
	    J=INDEX(STRING(I),':')
	    READ(STRING(I)(J+1:),*)GSUR_OLD
	  ELSE IF(INDEX(STRING(I),'Photospheric radius is') .NE. 0)THEN
	    J=INDEX(STRING(I),':')
	    K=INDEX(STRING(I),'(')
	    IF(K .NE. 0)READ(STRING(I)(J+1:K-1),*)RPHOT
	    IF(K .EQ. 0)READ(STRING(I)(J+1:),*)RPHOT
	  ELSE IF(INDEX(STRING(I),'Stars mass is') .NE. 0)THEN
	    J=INDEX(STRING(I),':')
	    READ(STRING(I)(J+1:),*)MASS_OLD
	  END IF
	END DO
	IF(MASS_OLD .EQ. 0)THEN
	  WRITE(6,*)'Old mass not available'
	  CALL GEN_IN(MASS_OLD,'New mass in solar units')
	END IF
	WRITE(T_OUT,'(1X,A,F8.2)')'Old mass is ',MASS_OLD
	MASS_NEW=MASS_OLD
	CALL GEN_IN(MASS_NEW,'New mass in solar units')
	GSUR_NEW=GSUR_OLD*MASS_NEW/MASS_OLD
!
        VTURB=0.0_LDP; CALL GEN_IN(VTURB,'Turbulent velcity in km/s)')
        LOW_LIM=1; CALL GEN_IN(LOW_LIM,'Depth to begin revised mass estimate')
        HIGH_LIM=ND; CALL GEN_IN(HIGH_LIM,'Depth to end revised mass estimate')
!
	READ(STRING(2),*)T1
	PLANE_PARALLEL=.FALSE.
	IF(T1/RND .LT. 1.5_LDP)PLANE_PARALLEL=.TRUE.
	CALL GEN_IN(PLANE_PARALLEL,'Plane paraellel model?')
!
        SUM_ERROR=0
        SUM_R=0
	OLD_FORMAT=.TRUE.
	IF(INDEX(STRING(1),'dTPdR/ROH') .NE. 0)OLD_FORMAT=.FALSE.
	DO I=1,ND
	  IF(OLD_FORMAT)THEN
	    dTPdR=0.0_LDP
	    READ(STRING(I+1),*)R,V,E,VdVdR,dPdR,g_TOT,g_RAD,g_ELEC,Gamma
	  ELSE
	    READ(STRING(I+1),*)R,V,E,VdVdR,dPdR,dTPdR,g_TOT,g_RAD,g_ELEC,Gamma
	  END IF
	  RSQ=R*R
	  IF(PLANE_PARALLEL)RSQ=RND*RND
	  P_GRAV(I)=MASS_NEW*GRAV_CON/RSQ
	  g_TOT=g_RAD-P_GRAV(I)
	  Gamma=g_RAD/P_GRAV(I)
	  IF(VTURB .NE. 0.0_LDP)THEN
	    dTPdR=-0.5_LDP*VTURB*VTURB*(2.0_LDP/R+VdVdR/V/V)
	  END IF
          DENOM=ABS(VdVdR)+ ABS(dPdR)+ ABS(dTPdR)+ABS(P_GRAV(I))+ABS(G_RAD) !ABS(g_TOT))
          E=200.0_LDP*(VdVdR+dPdR+dTPdR-g_TOT)/DENOM
          E_ON_GRAD=100.0_LDP*(VdVdR+dPdR+dTPdR-g_TOT)/G_RAD
          MT=g_rad/g_elec-1.0_LDP
          IF(I .GE. LOW_LIM .AND. I .LE. HIGH_LIM)THEN
            SUM_ERROR=SUM_ERROR+0.005_LDP*E/R**2/DENOM
            SUM_R=SUM_R+GRAV_CON/R**4/DENOM**2
          END IF
!
	  P_R(I)=R
	  P_VEL(I)=V
	  P_dPdR(I)=dPdR
	  P_dVdR(I)=1.0E-05_LDP*VdVdR/V
	  P_REQ(I)=VdVdR+dPdR+dTPdR+P_GRAV(I)              !GSUR_NEW*(RND/R)**2
	  P_GRAD(I)=g_RAD
	  P_GELEC(I)=g_ELEC
	  P_GTOT(I)=g_TOT
!
	  IF(R .GT. 9.99E+04_LDP)THEN
	    FMT='(1X,ES12.6,ES13.4,F9.2,6(ES14.4),2F11.2)'
	  ELSE
	    FMT='(1X,F12.6,ES13.4,F9.2,6(ES14.4),3F11.3)'
	  END IF
	  WRITE(LU_OUT,FMT)R,V,E,VdVdR,dPdR,dTPdR,g_TOT,g_RAD,g_ELEC,Gamma,MT,E_ON_GRAD
	END DO
!
        DEL_M=SUM_ERROR/SUM_R
        WRITE(T_OUT,*)'Better fit to hydrostatic equation will be obtained with a mass ',MASS_NEW-DEL_M
!
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Momentum equation is:',
	1                    ' VdV/dr = - dPdR/ROH - g + g_RAD'
	WRITE(LU_OUT,'(1X,A,A)')'        or          :',
	1                    ' VdV/dr = - dPdR/ROH + g_tot'
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Note: The definition of the error has changed (Sep-2020) '
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,A)')'Error is 200.0D0*(VdVdR+dPdR_ON_ROH-g_TOT)/',
	1        '( ABS(VdVdR)+ ABS(dPdR_ON_ROH)+ g + g_RAD )'
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A)')'Gamma = g_rad/g [g=g_GRAV] '
C
	WRITE(LU_OUT,'(1X,A)')' '
	WRITE(LU_OUT,'(1X,A,ES14.4,A,F6.3,A)')
	1          'Surface gravity is: ',GSUR_NEW,' (i.e., log g=',LOG10(GSUR_NEW),')'
	WRITE(LU_OUT,'(1X,A,F8.2,A)')
	1          'Stars mass is: ',MASS_NEW,' Msun'
	CLOSE(LU_OUT)
!
	XVEC(1:ND)=P_VEL(1:ND)
	XLAB='V(km/s)'
	YLAB=' '
!
5000	CONTINUE
	XOPT='P'
	CALL GEN_IN(XOPT,'Option to set plot or to plot (P)')
	CALL SET_CASE_UP(XOPT,IZERO,IZERO)
	!
	IF(XOPT .EQ. 'P')THEN
	  IF(YLAB(1:1) .EQ. ';')YLAB(1:)=YLAB(2:)
	  CALL GRAMON_PGPLOT(XLAB,YLAB,' ',' ')
	  YLAB=' '
!
	ELSE IF(XOPT .EQ. 'H' .OR. XOPT(1:2) .EQ. 'HE' .OR. XOPT .EQ. '?')THEN
	   WRITE(6,*)RED_PEN
	   WRITE(6,*)'XdVdR   -- set X axis to dVdR'
	   WRITE(6,*)'XVEL    -- set X axis to V(km/s)'
	   WRITE(6,*)'XR      -- set X axis to R/R(ND)'
!
	   WRITE(6,*)'GRAD    -- plot g(rad)'
	   WRITE(6,*)'GELEC   -- plot g(elec)'
	   WRITE(6,*)'GRAV    -- plot g'
	   WRITE(6,*)'dPdR    -- plot (1/roh).dP/dr'
	   WRITE(6,*)'dVdR    -- plot dV/dR'
!
	   WRITE(6,*)'REQ     -- plot VdVdR+dPdR+dTPdR+g'
	   WRITE(6,*)'NGL     -- plot g_l/(g-g_e)'
	   WRITE(6,*)'NREQ    -- plot (VdVdR+dPdR+dTPdR+g)/g_e'
	   WRITE(6,*)'NGRAD   -- plot g_r/g_e'
!
	   WRITE(6,*)DEF_PEN
!
	ELSE IF(XOPT .EQ. 'XVEL')THEN
	  XVEC(1:ND)=P_VEL(1:ND)
	  XLAB='V(km/s)'
!
	ELSE IF(XOPT .EQ. 'XDVDR')THEN
	  XVEC(1:ND)=1000.0_LDP*P_dVdR(1:ND)
	  XLAB='dVdR(ks\u-1\d)'
!
	ELSE IF(XOPT .EQ. 'XT')THEN
	  IOS=1
	  RVTJ_FILE_NAME='RVTJ'
	  DO WHILE(IOS .NE. 0)
            CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
            OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
            IF(IOS .NE. 0)WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    RVTJ_FILE_NAME='../RVTJ'
          END DO
	  TMP_STRING=' '
	  DO WHILE(INDEX(TMP_STRING,'Electron density') .EQ. 0)
	    READ(LU_IN,'(A)')TMP_STRING
	  END DO
	  READ(LU_IN,*)(ED(I),I=1,ND)
          CLOSE(LU_IN)
	  XVEC(1:ND)=ED(1:ND)*6.65E-25_LDP*10.0E+05_LDP/P_dVdR(1:ND)
	  XLAB='t'
!
	ELSE IF(XOPT .EQ. 'XR')THEN
	  XVEC(1:ND)=P_R(1:ND)/P_R(ND)
	  WRITE(6,*)'R(ND)=',P_R(ND)
	  XLAB='R/R(ND)'
!
	ELSE IF(XOPT .EQ. 'VEL')THEN
	  YVEC(1:ND)=P_VEL(1:ND)
	  CALL DP_CURVE(ND,XVEC,YVEC)
	  YLAB='V(km/s)'
!
	ELSE IF(XOPT .EQ. 'INT')THEN
!
! INT walks outward to the unique V=30 km/s crossing.  A repeated or
! non-monotonic velocity has multiple crossings, so the depth selected
! by this diagnostic is not unique.
!
	  DO I=1,ND-1
	    IF(P_VEL(I) .LE. 0.0_LDP .OR.
	1       P_VEL(I+1) .GE. P_VEL(I))THEN
	      WRITE(T_OUT,*)'Error in REV_HYDRO_FILE: INT requires'
	      WRITE(T_OUT,*)'positive, strictly decreasing velocity'
	      WRITE(T_OUT,*)'V(I),V(I+1)=',P_VEL(I),P_VEL(I+1)
	      STOP
	    END IF
	  END DO
!
	  I=0
	  T1=0.0_LDP
	  DO WHILE(I .LT. ND)
	    IF(P_VEL(I+1) .LE. 30.0_LDP)EXIT
	    I=I+1
	    XVEC(I)=P_R(I)
	    YVEC(I)=(2.00_LDP*P_GRAD(I)+(1.0_LDP-T1)*P_GRAD(I))/P_VEL(I)/3.0_LDP
	    ZVEC(I)=(P_GRAV(I)-T1*P_GRAD(I)/3.0_LDP)/P_VEL(I)
	    WRITE(6,'(I5,4ES14.6)')I,P_GRAV(I),P_GRAD(I),YVEC(I),ZVEC(I)
	  END DO
	  IF(I .LT. 2)THEN
	    WRITE(T_OUT,*)'Error in REV_HYDRO_FILE: INT could not locate'
	    WRITE(T_OUT,*)'a V=30 km/s crossing with enough bracket points'
	    STOP
	  END IF
!
	  T1=0.0_LDP; T2=0.0_LDP
	  DO J=1,I-1
	    T1=T1+(XVEC(J)-XVEC(J+1))*(YVEC(J)+YVEC(J+1)) 	
	    T2=T2+(XVEC(J)-XVEC(J+1))*(ZVEC(J)+ZVEC(J+1)) 	
	  END DO
	  T1=0.5_LDP*T1; T2=0.5_LDP*T2
	  WRITE(6,*)T1,T2
	  T2=T2+(P_VEL(1)-P_VEL(I))
	  T1=T1/T2
	  WRITE(6,*)'Factor to revise mass loss rate us',T1
	  WRITE(6,*)'Vinf factor is',T2
	  DO J=1,I
	    YVEC(J)=2.3205_LDP*YVEC(J)*P_R(J)/P_VEL(1)
	    ZVEC(J)=2.3205_LDP*ZVEC(J)*P_R(J)/P_VEL(1)
	    XVEC(J)=LOG10(P_R(J)/P_R(I))
	  END DO
	  WRITE(6,*)'Scaling radius for X-axis is: ',P_VEL(I)
	  CALL DP_CURVE(I,XVEC,YVEC)
	  CALL DP_CURVE(I,XVEC,ZVEC)
	  XLAB='Log R'
!
	ELSE IF(XOPT .EQ. 'DVDR')THEN
	  CALL DP_CURVE(ND,P_VEL,P_dVdR)
	  YLAB=TRIM(YLAB)//';dV\dg\u/dr'
	ELSE IF(XOPT .EQ. 'DPDR')THEN
	  CALL DP_CURVE(ND,P_VEL,P_dPdR)
	  YLAB=TRIM(YLAB)//'; \gr\u-1\d dP\dg\u/dr'
	ELSE IF(XOPT .EQ. 'GRAD')THEN
	  CALL DP_CURVE(ND,P_VEL,P_GRAD)
	  YLAB=TRIM(YLAB)//'; g\dr\u'
	ELSE IF(XOPT .EQ. 'GELEC')THEN
	  CALL DP_CURVE(ND,P_VEL,P_GELEC)
	  YLAB=TRIM(YLAB)//'; g\de\u'
	ELSE IF(XOPT .EQ. 'GRAV')THEN
	  CALL DP_CURVE(ND,P_VEL,P_GELEC)
	  YLAB=TRIM(YLAB)//'; g'
	ELSE IF(XOPT .EQ. 'REQ')THEN
	  CALL DP_CURVE(ND,XVEC,P_REQ)
	  YLAB=TRIM(YLAB)//'; ( vdv/dr + \gr\u-1\d dP/dr + g )'
	ELSE IF(XOPT .EQ. 'TST')THEN
	  DO I=1,ND
	    T2=(P_GRAD(I)-P_GELEC(I))*P_VEL(I)*P_R(I)*P_R(I)*P_dVdR(I)
	    IF(T2 .LE. 0.0_LDP)THEN
	      WRITE(T_OUT,*)'Error in REV_HYDRO_FILE: TST log argument'
	      WRITE(T_OUT,*)'is non-positive at depth',I,' : ',T2
	      STOP
	    END IF
	    TA(I)=LOG(T2)
	  END DO
	  CALL DP_CURVE(ND,XVEC,TA)
	  YLAB=TRIM(YLAB)//' '
	ELSE IF(XOPT .EQ. 'NGL')THEN
	  DO I=1,ND
	    TA(I)=(P_GRAD(I)-P_GELEC(I))/(P_GRAV(I)-P_GELEC(I))
	  END DO
	  CALL DP_CURVE(ND,XVEC,TA)
	  YLAB=TRIM(YLAB)//'; g\dl\u/(g-g\de\u)'
	ELSE IF(XOPT .EQ. 'NREQ')THEN
	  TA(1:ND)=P_REQ(1:ND)/P_GELEC(1:ND)
	  CALL DP_CURVE(ND,XVEC,TA)
	  YLAB=TRIM(YLAB)//'; g\dh\u/g\de\u'
	ELSE IF(XOPT .EQ. 'NGRAD')THEN
	  TA(1:ND)=P_GRAD(1:ND)/P_GELEC(1:ND)
	  YLAB=TRIM(YLAB)//'; g\dr\u/g\de\u'
	  CALL DP_CURVE(ND,XVEC,TA)
	ELSE IF(XOPT .EQ. 'EX' .OR. XOPT(1:2) .EQ. 'ST')THEN
	  STOP
	ELSE
	  WRITE(6,*)'Unrecognized option'
	  GOTO 5000
	END IF
	GOTO 5000
!
	END
