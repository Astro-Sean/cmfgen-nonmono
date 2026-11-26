!
! Designed to be called by CMFGEN, and outputs:
!            ION_ID, Level name, SL, Level.
!
	SUBROUTINE WR_LEVEL_LINKS
	USE SET_KIND_MODULE
	USE MOD_CMFGEN
	IMPLICIT NONE
	INTEGER LU,ID,J,NL
!
! Altered: 24-Oct-2023 -- Cleaning; output file name changed.
! Created: 15-Oct-2023
!
	CALL GET_LU(LU,'In WR_LEVEL_LINKS')
	OPEN(UNIT=LU,FILE='LEVEL_SL_STEQ_LINKS',STATUS='UNKNOWN')
	WRITE(LU,'(X,A,3X,3(2X,A),2X,A,3X,A)')'Ion(ID)','   ID','Level','   SL','STEQ(lev)','Level name'
	DO ID=1,NUM_IONS
           IF(ATM(ID)%XzV_PRES)THEN
             DO J=1,ATM(ID)%NXzV_F
	       NL= ATM(ID)%F_TO_S_XzV(J)+ ATM(ID)%EQXzV-1
               WRITE(LU,'(1X,A10,3(1X,I6),4X,I6,4X,A)')ION_ID(ID),ID,
	1                   J,ATM(ID)%F_TO_S_XzV(J),NL,ATM(ID)%XzVLEVNAME_F(J)
             END DO
           END IF
        END DO
	WRITE(LU,'(A)')' '
	CLOSE(UNIT=LU)
!
	RETURN
	END
