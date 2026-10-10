!
! Designed to modify the density structure in DISPGEN.
! the WRRVSIG option can be then used to create a file as input for a CMFGEN model.
!
!Options:
!    SIN -- Clumping as specified in Hillier 1991 paper.
!
	SUBROUTINE MODIFY_DENSITY(MASS_DENSITY,R,V,TYPE_ATM,ND)
	USE SET_KIND_MODULE
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Created 22-Nov-2021
!
	INTEGER ND
	REAL(KIND=LDP) MASS_DENSITY(ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	CHARACTER(LEN=*) TYPE_ATM
!
! Local variables
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND)       !Temporary vectors
	REAL(KIND=LDP) NEW_DENSITY(ND)
	REAL(KIND=LDP) MASS(ND)
	REAL(KIND=LDP) TIME(ND)
	REAL(KIND=LDP) XV(ND)
	REAL(KIND=LDP) YV(ND)
!
	REAL(KIND=LDP) PI
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) SCALE_FAC
	REAL(KIND=LDP) SM_LEN
	INTEGER I,K
	INTEGER IST,IEND
	INTEGER NSHELL
	LOGICAL SMOOTH
	LOGICAL UP_STR_AVERAGING
	CHARACTER(LEN=10) DEN_OPTION
!
	DEN_OPTION='SIN'
	CALL USR_OPTION(SCALE_FAC,'SF','1.0D0','Factor to scale mass density by')
	CALL USR_OPTION(SMOOTH,'SMOOTH','T','Adjust density to prserve mdot')
	CALL USR_OPTION(SM_LEN,'SM_LEN','0.5','Smoothing length log space')
	CALL USR_OPTION(NSHELL,'NSHELL','20','Nuber of shells')
	CALL USR_OPTION(UP_STR_AVERAGING,'UP_AV','F','Only averae of 0.5*SM_LEN inwards')
!
	PI=ACOS(-1.0_LDP)
	XV=LOG10(R/R(ND))	
	IF(DEN_OPTION .EQ. 'SIN')THEN
          DO I=1,ND
            T1=MOD(NSHELL*LOG(R(I)/R(ND)),LOG(R(1)/R(ND)))/LOG(R(1)/R(ND))
            T2=V(I)*(SIN(2*PI*T1)-0.4_LDP)/V(1)
            TA(I)=10**T2
          END DO
	  NEW_DENSITY=MASS_DENSITY*TA
	END IF
!
	IF(SMOOTH)THEN
!
! We integrate to det the mass in each shell, and the time step, dt
!            dM=M(I+1)-M(I); dt=TIME(I+1)-TIME(I)
!
	  TA=NEW_DENSITY*R*R
	  CALL TORSCL(MASS,TA,R,TB,TC,ND,'ZERO',TYPE_ATM)
          TA=1.0_LDP/V
	  CALL TORSCL(TIME,TA,R,TB,TC,ND,'ZERO',TYPE_ATM)
!
	  IST=1
	  DO I=1,ND
	    DO K=IST,I
	      IF(LOG10(R(K)/R(I)) .LT. 0.5_LDP*SM_LEN)THEN
	         IST=K
	         EXIT
	      END IF
	    END DO
!
	    DO K=IST,ND
	      IF(LOG10(R(K)/R(I)) .LT. -0.5_LDP*SM_LEN)THEN
	        IEND=K
	        EXIT
	      END IF
	    END DO
	    IF(UP_STR_AVERAGING)IST=MIN(I,IEND-1)
	    YV(I)=(MASS(IEND)-MASS(IST))/(TIME(IEND)-TIME(IST))

	  END DO
	  YV=YV/YV(ND)                   !Normalize
	  CALL DP_CURVE(ND,XV,YV)
	END IF
	MASS_DENSITY=NEW_DENSITY/YV
!
	RETURN
	END
