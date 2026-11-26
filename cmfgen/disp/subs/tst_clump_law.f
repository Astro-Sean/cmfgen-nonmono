!
! Routine to be used with DISPGEN to test different clumping laws.
!
	SUBROUTINE TST_CLUMP_LAW(CLUMP_FAC,R,VEL,TAU,ND)
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Altered 05-Oct-2026 - Reject velocity clumping tests on non-monotonic V.
! Altered 20-May-2021 - Extra clump options added (updated from OSIRS
!                          21-May-2021.
! Created 01-Oct-2018
!
	INTEGER ND
	REAL(KIND=LDP) CLUMP_FAC(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) VEL(ND)
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP), ALLOCATABLE :: VNODE(:), FVAL(:)
!
	INTEGER, PARAMETER :: NPAR_MAX=6
	INTEGER, PARAMETER :: LUER=6
	REAL(KIND=LDP) CLUMP_PAR(NPAR_MAX)
	REAL(KIND=LDP) T1,T2
	INTEGER NPAR
	INTEGER I,J,K,IOS,LU
	CHARACTER(LEN=10) CLUMP_LAW
	CHARACTER(LEN=120) FILENAME
!
        CHARACTER(LEN=30) UC
        EXTERNAL UC
!
! Every test law is parameterized by V or spline nodes on V.  A repeated
! velocity gives the same factor at different radii, so the test is not
! meaningful without a single-valued V coordinate.
!
	DO I=1,ND
	  IF(VEL(I) .LE. 0.0_LDP .OR.
	1     (I .GT. 1 .AND. VEL(I) .GE. VEL(I-1)))THEN
	    WRITE(LUER,*)'Error in TST_CLUMP_LAW:'
	    WRITE(LUER,*)'VEL must be positive and strictly decreasing'
	    WRITE(LUER,*)'Depth, VEL(I-1), VEL(I)=',I,
	1          VEL(MAX(1,I-1)),VEL(I)
	    STOP
	  END IF
	END DO
!
	WRITE(LUER,'(A)')BLUE_PEN
	WRITE(LUER,'(A)')'EXPO:   F=C1 + (1-C1-C3)EXP(-V/C2) + C3.EXP(-V/C4) (4 params)'
	WRITE(LUER,'(A)')'EXPO:   F=C1 + (1-C1-C3)EXP(-A/C2) + C3.EXP(-A/C4), A=MAX[0,(V-C5)] C6 softens MAX function(6 params)'
	WRITE(LUER,'(A)')'MEXP:   F=(C1 + (1-C1-C3)EXP(-V/C2))EXP(-R/R(1))/(1.0+C2*EXP(T1); T1=(1-MAX(1-C4*R(1)/R)**C5'
	WRITE(LUER,'(A)')'REXP:   F=C1 + (1-C1)EXP(-V/C2)    + (1-C1)EXP((V-VINF)/C3)'
	WRITE(LUER,'(A)')'POW:    F=C1 + (1-C1)(-V/VINF)**C2'
	WRITE(LUER,'(A)')'RPOW:   F=1/[ (1/C1-1)*(-V/VINF)**C2 ]'
	WRITE(LUER,'(A)')DEF_PEN
!
	CALL USR_OPTION(CLUMP_LAW,'LAW','EXPO','Clumping law: EXPO, MEXP, REXP, POW, RPOW, SPLINE')
	CLUMP_LAW=UC(CLUMP_LAW)
	IF(CLUMP_LAW .NE. 'SPLINE')THEN
	  CALL USR_OPTION(NPAR,'NPAR','2','Number of clumping parameters [ <6 ]')
	  IF(NPAR .LE. 0 .OR. NPAR .GT. 6)THEN
	    WRITE(6,*)'Invalid  number of NPAR --- valid ranges is 1 to 6'
	    RETURN
	  END IF
	END IF
!
	IF(CLUMP_LAW(1:4) .EQ. 'EXPO')THEN
!
! CLUMP_PAR(1) is the clumping factor at infinity.
! CLUMP_PAR(2) is a velocity, and determines how fast the clumping factor
! approach CLUMP_PAR(1).
!
	  CLUMP_PAR=0.0_LDP
	  CLUMP_PAR(4)=1.0_LDP
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  IF(NPAR .GT. 2)THEN
	    CALL USR_OPTION(CLUMP_PAR(3),'CLP3','0.1','Second clumping factor')
	    CALL USR_OPTION(CLUMP_PAR(4),'CLP4','1.0','Second velocity scale')
	  END IF
	  IF(NPAR .EQ. 6)THEN
	    CALL USR_OPTION(CLUMP_PAR(5),'CLP5','0.0','Hard minimum velocity')
	    CALL USR_OPTION(CLUMP_PAR(6),'CLP6','0.2','Smoothing param')
	  END IF
	  T2=CLUMP_PAR(6)
	  DO K=1,ND
	    T1=VEL(K)
	    IF(NPAR .EQ. 6)T1=LOG(1.0_LDP+EXP(T2*(VEL(K)-CLUMP_PAR(5))))/T2
	    CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1)-CLUMP_PAR(3))*
	1                   EXP(-T1/CLUMP_PAR(2))+
	1                   CLUMP_PAR(3)*EXP(-T1/CLUMP_PAR(4))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:4) .EQ. 'MEXP')THEN
!
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity is CLP1/(1+CLP3)')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  CALL USR_OPTION(CLUMP_PAR(3),'CLP3','1.0','Second clumping factor')
	  CALL USR_OPTION(CLUMP_PAR(4),'CLP4','0.5','Radius scale (< 1.0')
	  CALL USR_OPTION(CLUMP_PAR(5),'CLP5','0.5','Exponent ')
	  DO K=1,ND
	    T1=1.0_LDP-(MAX(1.0_LDP,CLUMP_PAR(4)*R(1)/R(K)))**CLUMP_PAR(5)
	    CLUMP_FAC(K)=(CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*EXP(-VEL(K)/CLUMP_PAR(2)))/
	1                     (1.0_LDP+CLUMP_PAR(3)*EXP(T1))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:4) .EQ. 'REXP')THEN
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','500','Velocity scale factor')
	  CALL USR_OPTION(CLUMP_PAR(3),'CLP3','0.1','Second velocity factor')
	  CALL USR_OPTION(CLUMP_PAR(4),'CLP4','0.0','Hard minimum velocity')
	  CALL USR_OPTION(CLUMP_PAR(5),'CLP5','0.2','Smoothing param')
	  T2=CLUMP_PAR(5)
	  DO K=1,ND
	    T1=LOG(1.0_LDP+EXP(T2*(VEL(K)-CLUMP_PAR(4))))/T2
	    CLUMP_FAC(K)=CLUMP_PAR(1)+(1.0_LDP-CLUMP_PAR(1))*
	1         EXP(-T1/CLUMP_PAR(2))+
	1        (1.0_LDP-CLUMP_PAR(1))*EXP( (VEL(K)-VEL(1))/CLUMP_PAR(3))
	  END DO
!
	ELSE IF(CLUMP_LAW(1:3) .EQ. 'POW')THEN
	  IF(NPAR .NE. 2)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)' WRONG VALUE NPAR=',NPAR
	    WRITE(LUER,*)'NPAR assumed to be 2'
	  END IF
	  CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	  CALL USR_OPTION(CLUMP_PAR(2),'CLP2','1.0D0','Power law exponent (+ve)')
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0_LDP-(1.0_LDP-CLUMP_PAR(1))*(VEL(K)/VEL(1))**CLUMP_PAR(2)
	  END DO
        ELSE IF(CLUMP_LAW(1:4) .EQ. 'RPOW')THEN
	  IF(NPAR .EQ. 2 .OR. NPAR .EQ. 4)THEN
	    CALL USR_OPTION(CLUMP_PAR(1),'CLP1','0.1','Clumping factor at infinity')
	    CALL USR_OPTION(CLUMP_PAR(2),'CLP2','1.0D0','Power law exponent (+ve)')
	    CLUMP_PAR(3)=0.0_LDP; CLUMP_PAR(4)=1.0_LDP
	  END IF
	  IF(NPAR .EQ. 4)THEN
	    CALL USR_OPTION(CLUMP_PAR(3),'CLP3','0.0','Clumping factor at infinity')
	    CALL USR_OPTION(CLUMP_PAR(4),'CLP4','2.0D0','Power law exponent (+ve)')
	  END IF
	  IF(NPAR .NE. 2 .AND. NPAR .NE. 4)THEN
	    WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	    WRITE(LUER,*)' WRONG VALUE N_CLUMP_PAR=',NPAR
	    STOP
	  END IF
	  DO K=1,ND
	    T1=1.0_LDP/CLUMP_PAR(1)-1.0_LDP
	    CLUMP_FAC(K)=1.0_LDP/(1.0_LDP+T1*(VEL(K)/VEL(1))**CLUMP_PAR(2))/
	1         (1.0_LDP+CLUMP_PAR(3)*(VEL(K)/VEL(1))**CLUMP_PAR(4))
	  END DO
!
        ELSE IF(CLUMP_LAW(1:6) .EQ. 'SPLINE')THEN
	  CALL USR_OPTION(FILENAME,'FILE','CLUMP_NODES','File with clump claw (or TERM')
	  IF(UC(FILENAME(1:4)) .EQ. 'TERM')THEN
	    LU=5
	  ELSE
	    LU=11
	    OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS,ACTION='READ')
	  END IF
	  READ(LU,*)NPAR
	  NPAR=NPAR+2
	  ALLOCATE(VNODE(NPAR),FVAL(NPAR))
	  DO I=2,NPAR-1
	    READ(LU,*)VNODE(I),FVAL(I)
	  END DO
	  IF(LU .EQ. 11)CLOSE(UNIT=LU)
	  IF(VNODE(2) .LT. VNODE(3))THEN
	    DO I=1,NPAR/2
	      J=NPAR-I+1
	      T1=VNODE(I); VNODE(I)=VNODE(J); VNODE(J)=T1
	      T1=FVAL(I);  FVAL(I)=FVAL(J);   FVAL(J)=T1
	    END DO
	  END IF
	  VNODE(1)=VEL(1); FVAL(1)=FVAL(2)
	  VNODE(NPAR)=VEL(ND); FVAL(NPAR)=FVAL(NPAR-1)
	  I=1
	  CALL MON_INTERP(CLUMP_FAC,ND,I,VEL,ND,FVAL,NPAR,VNODE,NPAR)
	  DEALLOCATE(VNODE,FVAL)
!
	ELSE
	  WRITE(LUER,*)'Error in SET_ABUND_CLUMP'
	  WRITE(LUER,*)'Invalid law for computing clumping factor'
	  WRITE(LUER,*)'Set it to 1'
	  DO K=1,ND
	    CLUMP_FAC(K)=1.0_LDP
	  END DO
	END IF
!
	RETURN
	END
