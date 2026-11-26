	SUBROUTINE DO_SL_LS_SPLIT(F_TO_S,INT_SEQ,NF,NS,LEVEL_NAMES,SL_OPTION)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Created 23-Feb-2018: Desing to split the lower N LS sates into individual
!                      super levels.
!
	INTEGER NF
	INTEGER NS
	INTEGER F_TO_S(NF)
	INTEGER INT_SEQ(NF)
	CHARACTER(LEN=*) LEVEL_NAMES(NF)
	CHARACTER(LEN=*) SL_OPTION
	LOGICAL, ALLOCATABLE :: DONE(:)
	INTEGER, ALLOCATABLE :: OLD_F_TO_S(:)
!
	INTEGER I,J,M
	INTEGER KL,KU
	INTEGER SPLIT_N
	INTEGER IOS,CNT
	INTEGER LAST_SL
!
!	SPLIT_LS_N
!
	READ(SL_OPTION(10:),*,IOSTAT=IOS)SPLIT_N
	IF(IOS .NE. 0 .OR. SPLIT_N .LE. 0)THEN
	  WRITE(6,*)'Error reading SPLIT_N in DO_SL_LS_SPLIT'
	  STOP
	END IF
!
	ALLOCATE(DONE(NF))
	ALLOCATE(OLD_F_TO_S(NF))
	OLD_F_TO_S=F_TO_S
	LAST_SL=MAXVAL(F_TO_S)
!
	DONE=.FALSE.
	CNT=0
	DO I=1,NF
	  IF(CNT .EQ. SPLIT_N)EXIT
	  IF(.NOT. DONE(I))THEN
	    DONE(I)=.TRUE.
	    KL=INDEX(LEVEL_NAMES(I),'[')
	    IF(KL .NE. 0)THEN
	      DO J=I+1,NF
	        IF(DONE(J))THEN
	        ELSE IF(LEVEL_NAMES(I)(1:KL) .EQ. LEVEL_NAMES(J)(1:KL))THEN
	          DONE(J)=.TRUE.
	        ELSE IF(F_TO_S(I) .EQ. F_TO_S(J))THEN
	          LAST_SL=LAST_SL+1
	          F_TO_S(J)=LAST_SL
	          DONE(J)=.TRUE.
	          CNT=CNT+1
	          KU=INDEX(LEVEL_NAMES(J),'[')
	          IF(KU .NE.0)THEN
	            DO M=J+1,NF
	              IF(LEVEL_NAMES(J)(1:KU) .EQ. LEVEL_NAMES(M)(1:KU))THEN
	                DONE(M)=.TRUE.
	                F_TO_S(M)=LAST_SL
	              END IF
	            END DO
	          END IF
	        END IF
	      END DO
	      IF(CNT .GT. SPLIT_N)EXIT
	    END IF
	  END IF
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
