	SUBROUTINE DO_dE_SPLIT(F_TO_S,INT_SEQ,ENERGY,NF,NS,LEVEL_NAMES,dE_OPTION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 23-Feb-2018: Desing to split the lower N LS sates into individual
!                      super levels.
!
	INTEGER NF
	INTEGER NS
!
	REAL(KIND=LDP) ENERGY(NF)
	REAL(KIND=LDP) dE_LIMIT
	INTEGER F_TO_S(NF)
	INTEGER INT_SEQ(NF)
	CHARACTER(LEN=*) LEVEL_NAMES(NF)
	CHARACTER(LEN=*) dE_OPTION
	LOGICAL, ALLOCATABLE :: DONE(:)
	INTEGER, ALLOCATABLE :: OLD_F_TO_S(:)
!
	REAL(KIND=LDP) LOC_dE_LIMIT
	INTEGER I,J,M
	INTEGER KL,KU
	INTEGER SPLIT_N
	INTEGER IOS,CNT
	INTEGER LAST_SL
	LOGICAL SAME_LS_STATE
	LOGICAL DO_LS_STATES
!
	READ(dE_OPTION,*,IOSTAT=IOS)dE_LIMIT
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error reading de_LIMIT in DO_dE_SPLIT'
	  STOP
	END IF
!
	ALLOCATE(DONE(NF))
	ALLOCATE(OLD_F_TO_S(NF))
	OLD_F_TO_S=F_TO_S
	LAST_SL=MAXVAL(F_TO_S)
	LOC_dE_LIMIT=ABS(dE_LIMIT)
	IF(dE_LIMIT .LT. 0)THEN
	  DO_LS_STATES=.TRUE.
	ELSE
	  DO_LS_STATES=.FALSE.
	END IF
!
	DO I=1,NF
	  KL=INDEX(LEVEL_NAMES(I),'[')
	  DO J=I+1,NF
	      KU=INDEX(LEVEL_NAMES(J),'[')
	      SAME_LS_STATE=.FALSE.
	      IF(KU*KL .GT. 0)THEN
	        IF(LEVEL_NAMES(J)(1:KU) .EQ.  LEVEL_NAMES(I)(1:KL))SAME_LS_STATE=.TRUE.
	      END IF
	      IF(F_TO_S(I) .NE. F_TO_S(J))THEN
	      ELSE IF(ENERGY(J)-ENERGY(I) .GT. LOC_dE_LIMIT .AND. DO_LS_STATES)THEN
	        LAST_SL=LAST_SL+1
	        F_TO_S(J)=LAST_SL
	        DO M=J+1,NF
	          IF(F_TO_S(M) .EQ. F_TO_S(I))F_TO_S(M)=LAST_SL
	        END DO
	        EXIT
!	      ELSE IF(ENERGY(J)-ENERGY(I) .GT. LOC_dE_LIMIT .AND. .NOT. SAME_LS_STATE)THEN
	      ELSE IF(ENERGY(J)-ENERGY(I) .GT. LOC_dE_LIMIT .AND. KU .EQ. 0)THEN
	        LAST_SL=LAST_SL+1
	        F_TO_S(J)=LAST_SL
	        DO M=J+1,NF
	          IF(F_TO_S(M) .EQ. F_TO_S(I))F_TO_S(M)=LAST_SL
	        END DO
	        EXIT
	      END IF
	   END DO
	END DO
!
! Make SL assignments sequential.
!
        CNT=0
        F_TO_S(1:NF)=-F_TO_S(1:NF)
        DO I=1,NF
          IF(F_TO_S(I) .LT. 0)THEN
            CNT=CNT+1
            DO J=I+1,NF
              IF(F_TO_S(J) .EQ. F_TO_S(I))F_TO_S(J)=CNT
            END DO
            F_TO_S(I)=CNT
	  END IF
	END DO
	NS=CNT
!
! Fix sequence number.
!
        DO I=1,NF
          IF(INT_SEQ(I) .NE. 0)THEN
            DO J=1,NF
              IF(INT_SEQ(I) .EQ. OLD_F_TO_S(J))THEN
                INT_SEQ(I)=F_TO_S(J)
                EXIT
              END IF
            END DO
          END IF
        END DO
!
	WRITE(145,'(A)')
	DO I=1,NF
	  WRITE(145,'(A,T40,3I6)')TRIM(LEVEL_NAMES(I)),I,F_TO_S(I),OLD_F_TO_S(I)
	END DO
	DEALLOCATE(DONE,OLD_F_TO_S)
!
	RETURN
	END
