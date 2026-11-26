!
! Subroutine to designed to create a NEW R grid given an old R grid, and optical depth scale on this
! grid. The routine places points places points logarithmically in R and TAU.
!
! Altered 05-Oct-2026: Validate positive, strictly decreasing OLD_V.
!
	SUBROUTINE ADJUST_ATM_R_GRID(R,OLD_R,OLD_V,OLD_TAU,
	1           V_CON,V_RAT_MAX,
	1           IB_RAT,OB_RAT,DTAU2_ON_DTAU1,N_IB_INS,N_OB_INS,ND,NS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NS			!For old grid
	INTEGER ND			!For final grid.
!
	REAL(KIND=LDP) R(ND)			!Returned
	REAL(KIND=LDP) OLD_R(NS)		!Input
	REAL(KIND=LDP) OLD_V(NS)		!Input
	REAL(KIND=LDP) OLD_TAU(NS)		!Input
!
	REAL(KIND=LDP) LOG_OLD_R(NS)
	REAL(KIND=LDP) LOG_OLD_V(NS)
	REAL(KIND=LDP) LOG_OLD_TAU(NS)
!
	REAL(KIND=LDP) LOG_R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) LOG_V(ND)
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) LOG_TAU(ND)
!
	REAL(KIND=LDP) V_CON
	REAL(KIND=LDP) V_RAT_MAX
	REAL(KIND=LDP) IB_RAT			!dTAU(I+1)/dTAU(I) near inner boudary.
	REAL(KIND=LDP) OB_RAT			!dTAU(I+1)/dRAU(I) near outer bunary
	REAL(KIND=LDP) DTAU2_ON_DTAU1		!DTAU(2)/DTAU(1) at outer boundary.
!
	INTEGER N_IB_INS
	INTEGER N_OB_INS
!
	REAL(KIND=LDP) LOG_V_CON
	REAL(KIND=LDP) R_CON, LOG_R_CON
	REAL(KIND=LDP) TAU_CON, LOG_TAU_CON
	INTEGER I_CON
!
	REAL(KIND=LDP) dTAU
	REAL(KIND=LDP) dLOG_TAU			!d(LOG(TAU))
	REAL(KIND=LDP) dLOG_TAU_SAVE
	REAL(KIND=LDP) dLOG_TAU1
	REAL(KIND=LDP) dLOG_TAU2
	REAL(KIND=LDP) T1,T2,T3
!
	INTEGER LU
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LUER=6
	INTEGER I,J,K,L
	INTEGER IBEG
	INTEGER ND1
!
	WRITE(6,'(A)')
	WRITE(6,'(A)')' Entering ADJUST_ATM_GRID to define R grid'
	WRITE(6,'(A)')' See R_GRID_SELECTION for computational information'
	WRITE(6,'(A)')
!
! The connection search and log interpolation use V as a coordinate, so a
! repeated velocity has no unique radius and a non-positive value is invalid.
!
	DO I=1,NS
	  IF(OLD_V(I) .LE. 0.0_LDP .OR.
	1     (I .GT. 1 .AND. OLD_V(I) .GE. OLD_V(I-1)))THEN
	    WRITE(6,*)'Error in ADJUST_ATM_R_GRID:'
	    WRITE(6,*)'OLD_V must be positive and strictly decreasing'
	    WRITE(6,*)'Depth, OLD_V(I-1), OLD_V(I)=',I,
	1          OLD_V(MAX(1,I-1)),OLD_V(I)
	    STOP
	  END IF
	END DO
!
	DO I=1,NS
	  T1=OLD_TAU(I)
	  OLD_TAU(I)=OLD_TAU(I)/(0.1_LDP+MIN(OLD_V(1)/2,OLD_V(I)))
	  WRITE(30,'(I4,3ES14.3)')I,OLD_V(I),T1,OLD_TAU(I)
	END DO
!
	LOG_OLD_TAU=LOG(OLD_TAU)
	LOG_OLD_R=LOG(OLD_R)
	LOG_OLD_V=LOG(OLD_V)
	IF(V_RAT_MAX .LT. 1.0_LDP .OR. V_RAT_MAX .GT. 2.0_LDP)THEN
	  WRITE(6,*)'Invalid V_RAT_MAX'
	END IF
!
	CALL GET_LU(LU,'ADJUST_ATM_R_GRID')
	OPEN(UNIT=LU,FILE='R_GRID_SELECTION',STATUS='UNKNOWN',ACTION='WRITE')
!
	WRITE(LU,'(A,2ES12.4)')' Outer boundary optical depth is:',OLD_TAU(1),LOG_OLD_TAU(1)
	WRITE(LU,'(A,2ES12.4)')' Inner boundary optical depth is:',OLD_TAU(NS),LOG_OLD_TAU(NS)
!
! Estimate average dTAU spacing.
!
	J=ND-1-N_OB_INS-N_IB_INS
	dLOG_TAU=(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1))/J
	dLOG_TAU_SAVE=dLOG_TAU
!
! Find connection point.
!
	I_CON=1
	DO I=1,NS
	  IF(OLD_V(I) .LT. V_CON)THEN
	    I_CON=I-1
	    EXIT
	  END IF
	END DO
!
! The new grid contains the connection point.
!
	LOG_V_CON=LOG(V_CON)
	I=I_CON
	T1=(LOG_V_CON-LOG_OLD_V(I+1))/(LOG_OLD_V(I)-LOG_OLD_V(I+1))
	LOG_R_CON=(1.0_LDP-T1)*LOG_OLD_R(I+1)+T1*LOG_OLD_R(I)
	LOG_TAU_CON=(1.0_LDP-T1)*LOG_OLD_TAU(I+1)+T1*LOG_OLD_TAU(I)
!
	J=ND-1-N_IB_INS-N_OB_INS
	T1=(LOG_TAU_CON-LOG_OLD_TAU(1))/(LOG_OLD_TAU(NS)-LOG_OLD_TAU(1))
	IF(T1 .LT. 0.5_LDP)T1=0.5_LDP
	IF(T1 .GT. 0.7_LDP)T1=0.7_LDP
	ND1=J*T1
!
	R_CON=EXP(LOG_R_CON)
	V_CON=EXP(LOG_V_CON)
	TAU_CON=EXP(LOG_TAU_CON)
	R(ND1)=R_CON
	V(ND1)=V_CON
	LOG_TAU(ND1)=LOG_TAU_CON
	TAU(ND1)=EXP(LOG_TAU_CON)
!
	WRITE(6,*)'I_CON=',I_CON,TAU_CON
	WRITE(6,*)OLD_R(I_CON),R_CON,OLD_R(I_CON+1)
	WRITE(6,*)'ND1=',ND1
!
! We define grid about transition point. We choose grid to satisfy V
! criterion, and then we make sure change in DTAU is not too large.
!
	K=ND1-1; L=I_CON
	dLOG_TAU1=dLOG_TAU
50      CONTINUE
	LOG_TAU(K)=LOG_TAU(ND1)-dLOG_TAU1
	DO I=I_CON,1,-1
	  IF( LOG_OLD_TAU(I) .LT. LOG_TAU(K))THEN
	    L=I
	    EXIT
	  END IF
	END DO
	T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	V(K)=EXP(LOG_V(K))
	IF(V(K)/V(K+1) .GT. V_RAT_MAX)THEN
	  dLOG_TAU1=dLOG_TAU1/1.1_LDP
	  WRITE(6,*)'dLOG_TAU1=',dLOG_TAU1
	  GOTO 50
	END IF
	LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	R(K)=EXP(LOG_R(K))
	TAU(K)=EXP(LOG_TAU(K))
	WRITE(6,*)'Done I_CON-1'
	WRITE(6,*)R(K),V(K),TAU(K)
	WRITE(6,*)R(K+1),V(K+1),TAU(K+1)
!
	WRITE(6,*)' V_RAT_MAX=', V_RAT_MAX
	dLOG_TAU2=dLOG_TAU
60      CONTINUE
	K=ND1+1; L=I_CON
	LOG_TAU(K)=LOG_TAU(ND1)+dLOG_TAU2
	DO I=I_CON,NS
	  IF( LOG_OLD_TAU(I+1) .GT. LOG_TAU(K))THEN
	    L=I
	    EXIT
	  END IF
	END DO
	T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	V(K)=EXP(LOG_V(K))
	WRITE(6,*)V(K-1),V(K)
	IF(V(K-1)/V(K) .GT. V_RAT_MAX)THEN
	  dLOG_TAU2=dLOG_TAU2/1.1_LDP
	  GOTO 60
	END IF
	LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	R(K)=EXP(LOG_R(K))
	TAU(K)=EXP(LOG_TAU(K))
!
	WRITE(6,*)R(K),V(K),TAU(K)
!
! Now check that the step sizes at the connection point are similar.
!
	K=ND1+1      !I_CON+1
	IF( TAU(K)-TAU(K-1) .GT. 1.25_LDP*(TAU(K-1)-TAU(K-2)) )THEN
          WRITE(6,*)TAU(K),TAU(K-1),TAU(K-2)
	  TAU(K)=TAU(K-1)+1.25_LDP*(TAU(K-1)-TAU(K-2))
	  LOG_TAU(K)=LOG(TAU(K))
	  DO I=I_CON,NS
	    IF( LOG_OLD_TAU(I+1) .GT. LOG_TAU(K))THEN
	      L=I
	      EXIT
	    END IF
	  END DO
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
!
	ELSE IF( 1.25_LDP*(TAU(K)-TAU(K-1)) .LT. (TAU(K-1)-TAU(K-2)) )THEN
	  K=ND1-1
          TAU(K)=TAU(K+1)-(TAU(K+2)-TAU(K+1))/1.25_LDP
	  LOG_TAU(K)=LOG(TAU(K))
	  DO I=I_CON,1,-1
	    IF( LOG_OLD_TAU(I) .LT. LOG_TAU(K))THEN
	      L=I
	      EXIT
	    END IF
	  END DO
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
	END IF
!
	K=ND1
	WRITE(6,*)'Done check around connection point'
	WRITE(6,*)R(K-1),V(K-1),TAU(K-1)
	WRITE(6,*)R(K),V(K),TAU(K)
	WRITE(6,*)R(K+1),V(K+1),TAU(K+1)
!
! We now do the spacing towards the inner boudary.
!
	IBEG=I_CON
	dLOG_TAU=(LOG_OLD_TAU(NS)-LOG_TAU(K))/(ND-K)
	WRITE(6,*)'dLOG_TAU for inner region is',dLOG_TAU
	DO WHILE(K .LE. ND-10)
	  dTAU=TAU(K)-TAU(K-1)
	  K=K+1
	  T2=LOG(TAU(K-1)+1.25_LDP*dTAU)
	  T3=LOG_TAU(K-1)+dLOG_TAU
	  LOG_TAU(K)=MIN(T2,T3)
          TAU(K)=EXP(LOG_TAU(K))
!
	  DO I=IBEG,NS-1
	    IF( LOG_OLD_TAU(I+1) .GT. LOG_TAU(K))THEN
	      L=I
	      EXIT
	    END IF
	  END DO
	  IBEG=L
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
	  WRITE(6,*)K,R(K),TAU(K-1),TAU(K)
!
	  IF(T3 .LT. T2)EXIT
	END DO
	WRITE(6,*)K,R(K),V(K),TAU(K)
!
! We now assumed the verlocity law etc is well behaved, and choose a
! uniform grid in log(tau).
!
	T1=0.0_LDP; T2=1.0_LDP
	DO I=1,N_IB_INS+1
	  T2=T2/IB_RAT
	  T1=T1+T2
	END DO
	dLOG_TAU=(LOG_OLD_TAU(NS)-LOG_TAU(K))/(ND-K-N_IB_INS+T1-1)
	DO J=K+1,ND-N_IB_INS-1
	  LOG_TAU(J)=LOG_TAU(J-1)+dLOG_TAU
	END DO
	T1=dLOG_TAU
	DO J=ND-N_IB_INS,ND-1
	  T1=T1/IB_RAT
	  LOG_TAU(J)=LOG_TAU(J-1)+T1
	END DO
!
	WRITE(6,*)'Doing final part of inner boudary'
	L=I_CON
	DO K=K,ND-1
	  DO WHILE(LOG_OLD_TAU(L+1) .LT. LOG_TAU(K))
	    L=L+1
	  END DO
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
	  TAU(K)=EXP(LOG_TAU(K))
	  WRITE(6,*)K,R(K),V(K),TAU(K)
	END DO
!
	R(ND)=OLD_R(NS); LOG_R(ND)=LOG(R(ND))
	TAU(ND)=OLD_TAU(NS); LOG_TAU(ND)=LOG(TAU(ND))
	V(ND)=OLD_V(NS); LOG_V(ND)=LOG(V(ND))
!
! We now set the grid towards the outer boudary.
!
	WRITE(6,*)' '
	WRITE(6,*)'Doing grid towards outer boundary'
	WRITE(6,*)' '
!
	IBEG=I_CON
	K=ND1-1
	dLOG_TAU=(LOG_TAU(K)-LOG_OLD_TAU(1))/(ND1-1-N_OB_INS)
	DO WHILE(K .GT. 15)
	  K=K-1
	  LOG_TAU(K)=LOG_TAU(K+1)-dLOG_TAU
          TAU(K)=EXP(LOG_TAU(K))
	  dTAU=TAU(K+1)-TAU(K)
	  IF( dTAU .GT. 1.25_LDP*(TAU(K+2)-TAU(K+1))) THEN
            TAU(K)=TAU(K+1)-1.25_LDP*(TAU(K+2)-TAU(K+1))
	    LOG_TAU(K)=LOG(TAU(K))
	    T3=LOG_TAU(K+1)-LOG_TAU(K)
	  END IF
!
	  DO I=IBEG,1,-1
	    IF( LOG_OLD_TAU(I) .LT. LOG_TAU(K) )THEN
	      L=I
	      EXIT
	    END IF
	  END DO
	  IBEG=L
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
!
	  WRITE(6,*)K,R(K),V(K),TAU(K)
	  IF(T3 .EQ. dLOG_TAU)EXIT
	END DO
	WRITE(6,*)'Done first part of out grid'
	IF(TAU(K) .LT. OLD_TAU(1))THEN
	  WRITE(6,*)'Error TAU(K) .LT. OLD_TAU(1)'
	  WRITE(6,*)TAU(K),OLD_TAU(1)
	  WRITE(6,*)'K=',K
	  STOP
	END IF
!
	IF(N_OB_INS .EQ. 1)THEN
	  dLOG_TAU=(LOG_TAU(K)-LOG_OLD_TAU(1))/(K-N_OB_INS-1)
	ELSE
	  T1=0.0_LDP; T2=1.0_LDP
	  DO I=1,N_OB_INS
	    T2=T2/OB_RAT
	    T1=T1+T2
	  END DO
	  dLOG_TAU=(LOG_TAU(K)-LOG_OLD_TAU(1))/(K-N_OB_INS-1)
	END IF
!
	DO J=K-1,N_OB_INS+1,-1
	  LOG_TAU(J)=LOG_TAU(J+1)-dLOG_TAU
	END DO
	IF(N_OB_INS .EQ. 1)THEN
	  LOG_TAU(2)=LOG_OLD_TAU(1)+dLOG_TAU/OB_RAT
	ELSE
	  T1=dLOG_TAU
	  DO J=N_OB_INS+1,2,-1
	    T1=T1/IB_RAT
	    LOG_TAU(J)=LOG_TAU(J+1)-T1
	  END DO
	END IF
!
	WRITE(6,*)'Finalizing outer part',LOG_OLD_TAU(1)
	L=I_CON
	DO K=K-1,2,-1
	  WRITE(6,*)K,LOG_TAU(K),L
	  DO WHILE( LOG_OLD_TAU(L) .GT. LOG_TAU(K))
	    L=L-1
	  END DO
	  T1=(LOG_TAU(K)-LOG_OLD_TAU(L))/(LOG_OLD_TAU(L+1)-LOG_OLD_TAU(L))
	  LOG_V(K)=T1*LOG_OLD_V(L+1)+(1.0_LDP-T1)*LOG_OLD_V(L)
	  V(K)=EXP(LOG_V(K))
	  LOG_R(K)=T1*LOG_OLD_R(L+1)+(1.0_LDP-T1)*LOG_OLD_R(L)
	  R(K)=EXP(LOG_R(K))
	  TAU(K)=EXP(LOG_TAU(K))
	END DO
!
	R(1)=OLD_R(1); LOG_R(1)=LOG(R(1))
	TAU(1)=OLD_TAU(1); LOG_TAU(1)=LOG(TAU(1))
	V(1)=OLD_V(1); LOG_V(1)=LOG(V(1))
!
! done up to here.

	WRITE(LU,'(/,A)')' '
	WRITE(LU,'(A,ES12.4)') ' Average spacing in Ln(tau)      is:',dTAU
	WRITE(LU,'(A,I3)')     ' Number of points inserted at IB is:',N_IB_INS
	WRITE(LU,'(A,I3)')     ' Number of points inserted at OB is:',N_OB_INS
!
	WRITE(LU,'(5X,A,17X,A,(4X,A))')
	1      'I','R','   V(km/s)','       Tau','      dTAU','    Log(Tau)',
	1      '  dLog(Tau)',' dTAU/dTAU','V(I)/V(I-1)'
	DO I=1,ND-1
	  IF(I .NE. 1)T1=(TAU(I+1)-TAU(I))/MAX(TAU(I)-TAU(I-1),1.0E-10_LDP)
	  IF(I .NE. 1)T2=V(I-1)/V(I)
	  WRITE(LU,'(I6,ES18.8,7ES14.4)')I,R(I),V(I),TAU(I),TAU(I+1)-TAU(I),
	1                     LOG10(TAU(I)),LOG10(TAU(I+1)/TAU(I)),T1,T2
	END DO
	WRITE(LU,'(I6,ES18.8,7ES14.4)')ND,R(ND),V(ND),TAU(ND),0.0D0,LOG_TAU(ND),0.0D0,0.0D0,V(ND-1)/V(ND)
	CLOSE(UNIT=LU)
!
! Make sure grid is monotonic. We do it here since we can check the actual grid
! computed in the output file.
!
	DO I=1,ND-1
	  IF(R(I) .LE. R(I+1))THEN
	    WRITE(6,*)' Error in SET_SN_R_GRID: R grid not monotonic'
	    WRITE(6,'(I4,3ES20.8)')I,R(MAX(1,I-1)),R(I),R(MIN(I+1,ND))
	    STOP
	  END IF
	END DO
!
	RETURN
	END
