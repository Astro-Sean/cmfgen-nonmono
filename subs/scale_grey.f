!
! Routine reads in the ratio of T/TGREY for an old model. This
! ratio is then applied to the passed GREY temperature distributon.
!
	SUBROUTINE SCALE_GREY(TGREY,TAUROSS,IOS,LUIN,ND)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 21-Dec-2004: Call changed. Routine now returns IOS is successful.
! Altered 02-Oct-2004: GREY_SCL_FAC_IN is now default input file.
! Altered 24-Aug-2002: Bug fixed with handling of boundary values.
!
	INTEGER ND
	INTEGER LUIN
	INTEGER IOS
	REAL(KIND=LDP) TGREY(ND)
	REAL(KIND=LDP) TAUROSS(ND)
!
	INTEGER LU_ER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL(KIND=LDP) LOG_TAU(ND)
	REAL(KIND=LDP) SCALE_FAC(ND)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_TAU(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_SCALE_FAC(:)
	CHARACTER*80 STRING
!
	INTEGER I,J
	INTEGER ND_RD
!
	IOS=0
	OPEN(UNIT=LUIN,FILE='GREY_SCL_FAC_IN',IOSTAT=IOS,STATUS='OLD',ACTION='READ')
	IF(IOS .NE. 0)THEN
	   LU_ER=ERROR_LU()
	   WRITE(LU_ER,*)'Warning: unable to open file GREY_SCL_FAC_IN'
	   WRITE(LU_ER,*)'Will try to open GREY_SCL_FAC (older file name)'
	   OPEN(UNIT=LUIN,FILE='GREY_SCL_FAC',IOSTAT=IOS,STATUS='OLD',ACTION='READ')
	   IF(IOS .NE. 0)THEN
	     LU_ER=ERROR_LU()
	     WRITE(LU_ER,*)'Warning: unable to open file GREY_SCL_FAC'
	     WRITE(LU_ER,*)'No scaling of TGREY performed in SCALE_GREY'
	     RETURN
	  END IF
	  WRITE(LU_ER,*)'GREY_SCL_FAC successfully opened'
	END IF
	LOG_TAU(1:ND)=LOG10(TAUROSS(1:ND))
!
! Skip comments and blank lines
!
	STRING=' '
	DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error reading GREY_SCL_FAC_IN: IOSTAT=',IOS
	    RETURN
	  END IF
	END DO
!
! The first non-comment or blank line should be the number of depth points.
!
	READ(STRING,*)ND_RD
	ND_RD=ND_RD+2
	ALLOCATE (OLD_TAU(ND_RD))
	ALLOCATE (OLD_SCALE_FAC(ND_RD))
	DO I=2,ND_RD-1
	  READ(LUIN,*)OLD_TAU(I),OLD_SCALE_FAC(I)
	  IF(IOS .NE. 0)THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)'Error reading Tau & T from GREY_SCL_FAC_IN: IOSTAT=',IOS
	    RETURN
	  END IF
	END DO
!
	OLD_TAU(ND_RD)=MAX(LOG_TAU(ND),OLD_TAU(ND_RD-1)+0.1_LDP)
	OLD_SCALE_FAC(ND_RD)=1.0_LDP
!
	OLD_TAU(1)=MIN(LOG_TAU(1),OLD_TAU(2)-0.1_LDP)
	OLD_SCALE_FAC(1)=OLD_SCALE_FAC(2)
!
	DO I=2,ND_RD
	  IF(OLD_TAU(I-1) .GE. OLD_TAU(I))THEN
	    WRITE(6,*)'Error is SCALE_GREY -- old TAU scale is not monotonic'
	    WRITE(6,'(10ES14.6)')(OLD_TAU(J),J=1,ND_RD)
	    STOP
	  END IF
	END DO
	DO I=2,ND
	  IF(LOG_TAU(I-1) .GE. LOG_TAU(I) .OR. LOG_TAU(I) .NE. LOG_TAU(I))THEN
	    WRITE(6,*)'Error is SCALE_GREY -- new TAU scale is not monotonic'
	    WRITE(6,'(10ES14.6)')(LOG_TAU(J),J=1,ND)
	    STOP
	  END IF
	END DO
	WRITE(6,*)OLD_TAU(1),OLD_TAU(ND_RD),'In SCALE_GREY old then new'
	WRITE(6,*)LOG_TAU(1),LOG_TAU(ND_RD)
	  
	CALL LIN_INTERP(LOG_TAU,SCALE_FAC,ND,OLD_TAU,OLD_SCALE_FAC,ND_RD)
!
	DO I=1,ND
	  TGREY(I)=TGREY(I)*SCALE_FAC(I)
	END DO
!
! Clean up
!
	DEALLOCATE (OLD_TAU)
	DEALLOCATE (OLD_SCALE_FAC)
	CLOSE(LUIN)
!
	RETURN
	END
