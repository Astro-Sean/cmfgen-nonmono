!
	SUBROUTINE INSERT_PTS_FOR_3D(ND,OLD_R,NCLUMPS,D_SH_IND,SHELL_LOC,MAX_NEW_ND,
	1                            NEW_R,ADD_N_PTS,INSERT_PTS_IN_DIR_OF_SHIFT,
	1                            SHELL_LEN_SHIFT,NEW_ND, IERROR, DEBUG)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 01-May-2025 : Called changed (IERROR and DEBUG).
!
! INPUTS
!
	INTEGER ND             ! From original grid
	INTEGER MAX_NEW_ND         ! MAX_NEW_ND for NEW_R 
!	
	REAL(KIND=LDP)  OLD_R(ND)
!
	REAL(KIND=LDP)  SHELL_LEN_SHIFT     ! The shift of the shell in R grid 
!
	INTEGER NCLUMPS
	INTEGER D_SH_IND 
	INTEGER SHELL_LOC(NCLUMPS)
!
	INTEGER ADD_N_PTS
!
	LOGICAL INSERT_PTS_IN_DIR_OF_SHIFT
	LOGICAL DEBUG
!
! OUTPUTS
!
	INTEGER IERROR
	INTEGER NEW_ND          ! New ND that will be used in main code (not MAX_NEW_ND)
	REAL(KIND=LDP)  NEW_R(MAX_NEW_ND)       ! New Radial grid
!
! LOCAL VARIABLES
!
	REAL(KIND=LDP)  T1,T2
	INTEGER I,J,K,L
	INTEGER LUOUT
	INTEGER TMP_IND
	LOGICAL MONOTONIC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	  IF(.NOT. INSERT_PTS_IN_DIR_OF_SHIFT)THEN
	    IF(ADD_N_PTS .NE. 0)THEN
	       NEW_R(1) = OLD_R(1)
	       J = 1
	       K = 1
	       DO I = 2,ND-10
		  IF(K.GT.NCLUMPS)THEN

	          ELSE IF(I .GT. SHELL_LOC(K)-(D_SH_IND-4) .AND. 
	1            I .LE. SHELL_LOC(K)+(D_SH_IND-4)      )THEN

	          !ELSE IF(I .EQ. SHELL_LOC(K)+D_SH_IND)THEN
	          !  K=K+1
!
	          ELSE

		    IF(I .EQ. SHELL_LOC(K)+D_SH_IND+2) K=K+1
		    TMP_IND = 0
		    !IF( ( OLD_R(I-1)-OLD_R(I) ) .GE. 0.03)TMP_IND=3
		    !IF( I .EQ. SHELL_LOC(K)-D_SH_IND )TMP_IND=3
	            T2 = (OLD_R(I)-OLD_R(I-1))/(ADD_N_PTS + TMP_IND)
	            !T2 = (OLD_R(I+1)-OLD_R(I))/(ADD_N_PTS + TMP_IND)
	            DO L = 2,ADD_N_PTS+TMP_IND
	               J = J + 1
	               NEW_R(J)=OLD_R(I-1)+(L-1)*T2
	               !NEW_R(J)=OLD_R(I)+(L-1)*T2
	            END DO
	          END IF
!
	          J = J + 1
	          NEW_R(J) = OLD_R(I)
	       END DO
	    ELSE 
	       WRITE(6,*)'ADD_N_PTS can not be LE to 0 '
	       STOP
	    END IF
	  ELSE 
	    IF(ADD_N_PTS .NE. 0 .AND. SHELL_LEN_SHIFT .GE. 0.0)THEN
	      NEW_R(1) = OLD_R(1)
	      J = 1
	      K = 1
	      DO I = 2,ND-10
		IF(K.GT.NCLUMPS)THEN

	        ELSE IF(I .GT. SHELL_LOC(K)-(D_SH_IND+5) .AND. 
	1          I .LT. SHELL_LOC(K)-(D_SH_IND-4)      )THEN

		   TMP_IND = 0
		   IF( I .EQ. SHELL_LOC(K)-D_SH_IND )TMP_IND=3
	           T2 = (OLD_R(I)-OLD_R(I-1))/(ADD_N_PTS + TMP_IND)
	           !T2 = (OLD_R(I+1)-OLD_R(I))/(ADD_N_PTS + TMP_IND)
	           DO L = 2,ADD_N_PTS+TMP_IND !1,
	              J = J + 1
	              NEW_R(J)=OLD_R(I-1)+(L-1)*T2
	              !NEW_R(J)=OLD_R(I)+(L-1)*T2
	           END DO
!
		ELSE IF( I .GE. SHELL_LOC(K)-(D_SH_IND-4) .AND. 
	1                I .LT. SHELL_LOC(K)+D_SH_IND     .AND.
	1                MOD(I,2) .EQ. 0 )THEN
	!	   CYCLE

	        ELSE IF(I .EQ. SHELL_LOC(K)+D_SH_IND)THEN
	           K=K+1
	        END IF
!
	        J = J + 1
	        NEW_R(J) = OLD_R(I)
	      END DO
	    ELSE IF(ADD_N_PTS .NE. 0 .AND. SHELL_LEN_SHIFT .LT. 0.0)THEN
	      NEW_R(1) = OLD_R(1)
	      J = 1
	      K = 1
	      DO I = 2,ND-10
		IF( K .GT. NCLUMPS)THEN

	        ELSE IF(I .GT. SHELL_LOC(K)+(D_SH_IND -4) .AND. 
	1          I .LT. SHELL_LOC(K)+(D_SH_IND +2)      )THEN

		   TMP_IND = 0
		   IF( I .EQ. SHELL_LOC(K)+(D_SH_IND + 1) )TMP_IND=3
	           T2 = (OLD_R(I)-OLD_R(I-1))/(ADD_N_PTS + TMP_IND)
	           DO L = 2,ADD_N_PTS+TMP_IND
	              J = J + 1
	              NEW_R(J)=OLD_R(I-1)+(L-1)*T2
	           END DO
!
		ELSE IF( I .GE. SHELL_LOC(K)-D_SH_IND     .AND. 
	1                I .LE. SHELL_LOC(K)+(D_SH_IND-3) .AND.
	1                MOD(I,2) .EQ. 0 )THEN
	!	   CYCLE

	        ELSE IF(I .EQ. SHELL_LOC(K)+D_SH_IND+3)THEN
	           K=K+1
	        END IF
!
	        J = J + 1
	        NEW_R(J) = OLD_R(I)
	      END DO
	    END IF
	  END IF
!
	  NEW_R(J:J+10)=OLD_R(ND-10:ND)
	  J = J + 10
	  NEW_ND = J
!
! Check newly created R grid is monotonic.
!
	IERROR=0
	MONOTONIC=.TRUE.
	DO I=1,NEW_ND-1
	  IF(NEW_R(I+1) .EQ. NEW_R(I) .AND. I .NE. ND-1 .AND. I .NE. 1)THEN
	    WRITE(6,*)'Adjusting 2 equal R grid points in INSERT_PTS_FOR_3D'
	    NEW_R(I) = (NEW_R(I+1) + NEW_R(I-1))/2.0
	    WRITE(6,'(1X,3(I6,3X,ES18.10))')(J,NEW_R(J),J=I-1,I+1)
	  ELSE IF(NEW_R(I+1) .GT. NEW_R(I))THEN
	    WRITE(6,*)'New R grid is not monotonic -- check file NEW_R_GRID'
	    WRITE(6,*)'Location is :',I,I+1
	    MONOTONIC=.FALSE.
	    IERROR=1
	  END IF
	END DO
!
	IF(DEBUG .OR. .NOT. MONOTONIC)THEN
	  CALL GET_LU(LUOUT,'In INSERT_PTS_FOR_3D')
	  OPEN(UNIT=LUOUT,FILE='NEW_R_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	  DO J=1,NEW_ND-1
	    WRITE(LUOUT,'(I5,2ES16.6)')J,NEW_R(J),NEW_R(J)-NEW_R(J+1)
	  END DO
	  WRITE(LUOUT,'(I5,2ES16.6)')J,NEW_R(J)
	END IF
!
	RETURN
!
	END SUBROUTINE

