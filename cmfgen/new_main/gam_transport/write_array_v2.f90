        SUBROUTINE WRITE_ARRAY_V2(A,ND,N,B,FILENAME1,INDX,OPTION)
	USE SET_KIND_MODULE
        USE GAM_MU_MOD
        IMPLICIT NONE
!
        INTEGER :: ND
        INTEGER :: N
        INTEGER :: I,J,K
	INTEGER :: INDX
	INTEGER :: V
        REAL(KIND=LDP) :: A(N,ND)
        REAL(KIND=LDP) :: B(N)
        CHARACTER(LEN=40) :: FILENAME1
        CHARACTER(LEN=40) :: FILENAME2
        CHARACTER(LEN=15) :: STRING
        CHARACTER(LEN=*) :: OPTION
!
        DO I=1,ND
	  FILENAME2=FILENAME1
          FILENAME2=ADJUSTL(FILENAME2)
          J=INDEX(FILENAME2,' ')
          WRITE(STRING,'(I2.2)')INDX
          FILENAME2(J:J+1)=STRING(1:2)
          FILENAME2(J+2:J+2)='_'
          J=INDEX(FILENAME2,' ')
          WRITE(STRING,'(I3.3)')I
          FILENAME2(J:J+2)=STRING(1:3)
          FILENAME2(J+3:)='.dat'
          OPEN(UNIT=7,FILE=FILENAME2,STATUS='UNKNOWN',ACTION='WRITE')
          IF(OPTION .EQ. 'ANGLE')THEN
            V=R_MU(I)%MU_PTS
            DO K=1,V
              WRITE(7,'(F11.6,ES18.8)')R_MU(I)%MU_VECTOR(K),A(K,I)
            END DO
          ELSE IF(OPTION .EQ. 'FREQ')THEN
            DO K=1,N
              WRITE(7,'(ES18.8,ES18.8)')B(K)*4.135668D-18,A(K,I)
            END DO
	  ELSE
            WRITE(6,*)'ERROR IN WRITE_ARRAY'
            WRITE(6,*)'Invalid option for writing testing arrays'
            STOP
          END IF
          CLOSE(UNIT=7)
        END DO
	RETURN
        END SUBROUTINE WRITE_ARRAY_V2
