!
! Small program to read in BA_ASI_N_D? file. It can be then used to solve
! the set of simultaneous equations. Usefule for testing.
!
	PROGRAM SOLVE_BA_MAT
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER N
	LOGICAL USE_DC
	CHARACTER*80 FILENAME
!
	FILENAME='BA_ASCI_N_D1'
	CALL GEN_IN(FILENAME,'File with BA and STEQ data')
!
	N=100
	CALL GEN_IN(N,'Number of unknown populations (NT)')
!
! This option has been left in, but was for a special pure H model only.
! It will be hidden for most test cases.
!
	USE_DC=.FALSE.
	IF(N .LE. 30)THEN
	  CALL GEN_IN(USE_DC,'Use departue coefficients (pure H model only)?')
	END IF
!
	CALL SUB_SOLVE_BA(N,USE_DC,FILENAME)
!
	STOP
	END
!
	SUBROUTINE SUB_SOLVE_BA(N,USE_DC,FILENAME)
	USE SET_KIND_MODULE
!
	INTEGER N
	CHARACTER*80 FILENAME
	LOGICAL USE_DC
!
	INTEGER, PARAMETER :: NSNG=1
	CHARACTER*1, PARAMETER :: NO_TRANS='N'
!
	REAL(KIND=LDP) POPS(N)
	REAL(KIND=LDP) PLTE(N)
	REAL(KIND=LDP) STEQ(N)
	REAL(KIND=LDP) CMAT(N,N)
	REAL(KIND=LDP) STAT_WT(N)
	REAL(KIND=LDP) EDGE(N)
!
	REAL(KIND=LDP) SAV_CMAT(N,N)
	REAL(KIND=LDP) SAV_STEQ(N)
	REAL(KIND=LDP) RHS(N)
!
	REAL(KIND=LDP) NEW_POPS(N)
	REAL(KIND=LDP) NEW_LTE(N)
	REAL(KIND=LDP) NEW_SOL(N)
	REAL(KIND=LDP) ED_NEW,DI_NEW,T_NEW
!
	REAL(KIND=LDP) ROW_SF(N)
	REAL(KIND=LDP) COL_SF(N)
!
	REAL(KIND=LDP) ROW_CND
	REAL(KIND=LDP) COL_CND
	REAL(KIND=LDP) MAX_VAL
!
	REAL(KIND=LDP) HDKT
	REAL(KIND=LDP) GION,RGU,X,Y
	REAL(KIND=LDP) IONIZATION_ENERGY
	REAL(KIND=LDP) SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
!
	INTEGER IFAIL
	INTEGER IPIVOT(N)
	CHARACTER*132 STRING
	INTEGER I,J,K,L
	LOGICAL ITERATE
!
	HDKT=4.7994145_LDP
	ITERATE=.TRUE.
!
	OPEN(UNIT=10, FILE=FILENAME,STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'POP') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(10,'(A)')STRING
	  END DO
	  DO I=1,N
	    READ(STRING(10:),*)POPS(I)
	    READ(10,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read POPS'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'STEQ') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING
	  DO WHILE(STRING .EQ. ' ')
	    READ(10,'(A)')STRING
	  END DO
	  DO I=1,N
	    READ(STRING(10:),*)STEQ(I)
	    READ(10,'(A)')STRING
	  END DO
	  WRITE(6,*)'Successfully read STEQ'
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'C_MAT') .EQ. 0)		
	    READ(10,'(A)')STRING
	  END DO
!
	  DO K=1,N,5
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ')
	      READ(10,'(A)')STRING
	    END DO
	    DO I=1,N
	      IF(I .NE. 1)READ(10,'(A)')STRING
	      READ(STRING(10:),*)(CMAT(I,J),J=K,MIN(K+4,N))
	    END DO
	  END DO
	  WRITE(6,*)'Successfully read CMAT'
!
	CLOSE(UNIT=10)
	SAV_CMAT=CMAT
	SAV_STEQ=STEQ
!
! This section is for a pure hydrogen model only.
!
	IF(USE_DC)THEN
	  WRITE(6,*)'Using departure coefficients'
	  IONIZATION_ENERGY=109678.7640_LDP
	  DO I=1,10
	    EDGE(I)=IONIZATION_ENERGY*SPEED_OF_LIGHT()*1.0E-15_LDP/I/I
	    STAT_WT(I)=2.0_LDP*I*I
	  END DO
	  GION=1.0_LDP
!
	  RGU=LOG(2.07078E-22_LDP)
          X=HDKT/POPS(N)
	  Y=POPS(N-1)*POPS(N-2)*( POPS(N)**(-1.5_LDP) )/GION
          WRITE(6,*)RGU,X,Y
          DO I=1,N-3
            PLTE(I)=STAT_WT(I)*Y*EXP(EDGE(I)*X+RGU)
	    WRITE(6,*)EDGE(I),PLTE(I),POPS(I)
	  END DO
!
	  DO I=1,N
	    DO J=1,10
	      CMAT(I,N-2)=CMAT(I,N-2)+CMAT(I,J)
	      CMAT(I,N-1)=CMAT(I,N-1)+CMAT(I,J)
	      CMAT(I,N)=CMAT(I,N)-(1.5_LDP+HDKT*EDGE(J)/POPS(N))*CMAT(I,J)
	    END DO
	  END DO
	  SAV_CMAT=CMAT
	  CALL WR2D_MA(CMAT,N,N,'C_MAT_D61',96)
	END IF
!
! Perform the LU decomposition using DGETRF. We first equilibrize the matrix
! using DGEEQU sot the the maximum row and column values are approximately
! unity.
!
	WRITE(6,*)'Starting equilibriation'
	CALL DGEEQU(N,N,CMAT,N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  WRITE(6,*)'Error performing equilibriation onf CMAT'
	  WRITE(6,*)'IFAIL=',IFAIl
	  STOP
	END IF
	WRITE(6,*)'Ending equilibriation'
!
        DO J=1,N
          STEQ(J)=STEQ(J)*ROW_SF(J)
          DO I=1,N
            CMAT(I,J)=CMAT(I,J)*ROW_SF(I)*COL_SF(J)
          END DO
        END DO
!
	WRITE(6,*)'Beginning LU decomposition'
	CALL DGETRF(N,N,CMAT,N,IPIVOT,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error ins solution using DGETRF'
          WRITE(6,*)'IFAIL=',IFAIL
	  STOP
	END IF
	WRITE(6,*)'Finised LU decomposition'
!
! Now perform the solution.
!
	WRITE(6,*)'Beginning solution'
        CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,STEQ,N,IFAIL)
        IF(IFAIL .NE. 0)THEN
          WRITE(6,*)'Error ins solution using DGETRS'
          WRITE(6,*)'IFAIL=',IFAIL
	  STOP
	END IF
	WRITE(6,*)'End solution'

        DO J=1,N
          STEQ(J)=STEQ(J)*COL_SF(J)
        END DO
!
	RHS=0.0_LDP
	DO J=1,N
	  DO I=1,N
	    RHS(I)=RHS(I)+SAV_CMAT(I,J)*STEQ(J)
	  END DO
	END DO
!
	IF(ITERATE)THEN
	  RHS=(SAV_STEQ-RHS)*ROW_SF
	  WRITE(6,*)'Beginning solution iteration'
          CALL DGETRS(NO_TRANS,N,NSNG,CMAT,N,IPIVOT,RHS,N,IFAIL)
	  WRITE(6,*)'Finished solution iteration'
          IF(IFAIL .NE. 0)THEN
            WRITE(6,*)'Error ins solution using DGETRS'
            WRITE(6,*)'IFAIL=',IFAIL
	    STOP
	  END IF
	  RHS=RHS*COL_SF
	  NEW_SOL=STEQ+RHS
	  WRITE(6,*)'Writing solution and iterated solution to UNIT 49'
	  WRITE(49,'(5X,A,6X,A,15X,A,12X,A)')'I','% Change','Old solution','New Solution'
	  DO I=1,N
	    T1=0.0_LDP
	    IF(NEW_SOL(I) .NE. 0)T1=100.0_LDP*(NEW_SOL(I)-STEQ(I))/NEW_SOL(I)
	    WRITE(49,'(1X,I5,2X,ES12.4,3X,2ES24.14)')I,T1,STEQ(I),NEW_SOL(I)
	  END DO
	  STEQ=NEW_SOL
	END IF
!
	IF(USE_DC)THEN
	  DI_NEW=POPS(N-2)*(1.0_LDP-STEQ(N-2))
	  ED_NEW=POPS(N-1)*(1.0_LDP-STEQ(N-1))
	  T_NEW=POPS(N)*(1.0_LDP-STEQ(N))
	  RGU=LOG(2.07078E-22_LDP)
          X=HDKT/T_NEW
	  Y=ED_NEW*DI_NEW*( T_NEW**(-1.5_LDP) )/GION
          DO I=1,N-3
            NEW_LTE(I)=STAT_WT(I)*Y*EXP(EDGE(I)*X+RGU)
	  END DO
	  DO I=1,N-3
	    NEW_POPS(I)=(POPS(I)/PLTE(I))*(1.0_LDP-STEQ(I))*NEW_LTE(I)
	  END DO
	  NEW_POPS(N-2)=DI_NEW
	  NEW_POPS(N-1)=ED_NEW
	  NEW_POPS(N)=T_NEW
	  DO I=1,N
	    WRITE(6,'(1X,I3,4ES14.4)')I,STEQ(I),SAV_STEQ(I),RHS(I),
	1                          1.0D0-NEW_POPS(I)/POPS(I)
	  END DO
	
	ELSE
	  DO J=1,N
	    DO I=1,N
	      CMAT(I,J)=SAV_CMAT(I,J)*STEQ(J)
	    END DO
	  END DO
	  CALL WR2D_MA(CMAT,N,N,'C_MAT_D61',50)
	  WRITE(6,*)'CMAT x solution vector written to unit 50'
!
	  RHS=0.0_LDP
	  DO J=1,N
	    DO I=1,N
	      RHS(I)=RHS(I)+SAV_CMAT(I,J)*STEQ(J)
	    END DO
	  END DO
	  WRITE(51,'(4X,A,6X,A,11X,A,14X,A)')'I','Solution','Input RHS','Output RHS'
	  DO I=1,N
	    WRITE(51,'(1X,I4,ES14.4,3X,2ES24.15)')I,STEQ(I),SAV_STEQ(I),RHS(I)
	  END DO
	  WRITE(6,*)'Solution and RHS vectors written to unit 51'
	END IF
!
	STOP
	END
