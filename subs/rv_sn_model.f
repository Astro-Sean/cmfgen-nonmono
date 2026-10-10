!
! Routine to compute radius points to be used in the comoving frame
! integration. The radius points are chosen to be equally spaced in
! LOG(Tau) where Tau is assumed to be dominated by free-free
! processes and is consequently proportional to the integral of the
! density squared.
!
	SUBROUTINE RV_SN_MODEL(R,V,SIGMA,RMAX,RP,VCORE,BETA1,LU,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LU
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
!
	REAL(KIND=LDP) RMAX
	REAL(KIND=LDP) RP
	REAL(KIND=LDP) VCORE
	REAL(KIND=LDP) BETA1
!
! Local arrays.
!
	REAL(KIND=LDP) TA(ND),TB(ND),TC(ND)
	REAL(KIND=LDP) T1,DLNR
!
	INTEGER NBND_INS
	INTEGER I,J,MND
	INTEGER NOLD,NDOLD
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	LOGICAL RDINR
	CHARACTER*80 STRING
!
	NBND_INS=2
!
! Check whether the passed parameters are valid.
!
	IF(BETA1 .LT. 0.0_LDP)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in RV_SN_MODEL --- Invalid BETA'
	  STOP
	END IF
	IF(NBND_INS .LT. 1 .OR. NBND_INS .GT. 3)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V3 --- Invalid NBND_INS'
          WRITE(LUER,*)'NBND_INS should be 1, 2, or 3'
	END IF
!
	RDINR=.FALSE.
	IF(RDINR)THEN
	  OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR')
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file.
!
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LU,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU)
!
	  READ(LU,*)TA(1),TA(1),NOLD,NDOLD
! Check relative values.
	  IF(ND .NE. NDOLD)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error-NDOLD and ND are not equal in RDINR'
	    WRITE(LUER,*)'NDOLD=',NDOLD,' ND=',ND
	    STOP
	  END IF
! TA is used for everything but R which is all we want.
	  DO I=1,ND
	    READ(LU,*)R(I),TA(I),TA(I),TA(I)
	    READ(LU,*)(TA(J),J=1,NOLD)
	  END DO
	  R(1)=RMAX
!
! Compute Velocity and SIGMA
!
	  DO I=1,ND
	    V(I)=VCORE*(R(I)/R(ND))**BETA1
	    SIGMA(I)=BETA1-1.0_LDP
	  END DO
	  R(ND)=RP
	  CLOSE(UNIT=LU)
	  RETURN
	END IF
!
	MND=ND-2*NBND_INS
	T1=LOG(RMAX/RP)
	T1=EXP(T1/(MND-1))
	TA(MND)=RP
	DO I=MND-1,2,-1
	  TA(I)=RP*(T1**(MND-I))
	  WRITE(127,'(ES14.4)')TA(I)
	END DO
	TA(1)=RMAX
!
! Insert finer grid near both boundaries.
!
	DO I=2,MND-1
	  R(I+NBND_INS)=TA(I)
	END DO
	R(1)=TA(1)
	R(ND)=TA(MND)
	IF(NBND_INS .EQ. 1)THEN
	  R(2)=TA(1)-(TA(1)-TA(2))/20.0
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/20.0_LDP
	ELSE IF(NBND_INS .EQ. 2)THEN
!	  R(2)=TA(1)-(TA(1)-TA(2))/10.0D0
!	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/10.0D0
!	  R(3)=TA(1)-(TA(1)-TA(2))/3.0D0
!	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/3.0D0
	  R(2)=TA(1)-(TA(1)-TA(2))/50.0_LDP
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/10.0_LDP
	  R(3)=TA(1)-(TA(1)-TA(2))/5.0_LDP
	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/3.0_LDP
	ELSE IF(NBND_INS .EQ. 3)THEN
	  R(2)=TA(1)-(TA(1)-TA(2))/20.0_LDP
	  R(ND-1)=R(ND)+(TA(MND-1)-TA(MND))/20.0_LDP
	  R(3)=TA(1)-(TA(1)-TA(2))/8.0_LDP
	  R(ND-2)=R(ND)+(TA(MND-1)-TA(MND))/8.0_LDP
	  R(4)=TA(1)-(TA(1)-TA(2))/3.0_LDP
	  R(ND-3)=R(ND)+(TA(MND-1)-TA(MND))/3.0_LDP
	END IF
!
! Compute Velocity and SIGMA
!
	DO I=1,ND
	  V(I)=VCORE*(R(I)/R(ND))**BETA1
	  SIGMA(I)=BETA1-1.0_LDP
	END DO
!
	WRITE(127,*)RP,RMAX,VCORE,BETA1
	DO I=1,ND
	  WRITE(127,'(3ES12.4)')R(I),V(I),SIGMA(I)
	END DO
!
	RETURN
	END
