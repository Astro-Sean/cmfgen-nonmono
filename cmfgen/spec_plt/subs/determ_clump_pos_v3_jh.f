!
! This is a simple routine to find the shells position from its center.
! Thus, the shell's postion will be between the indices 
! CLUMP_POS(J) - D_HALF_SH and CLUMP_POS(J) + D_HALF_SH
!
	SUBROUTINE DETERM_CLUMP_POS_V3_JH(DENSITY,ND,
	1                                 NCLUMPS,D_SH_IND,CLUMP_POS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! From CMFGEN
!
	INTEGER ND
	REAL(KIND=LDP) DENSITY(ND)
!
! Values to be returned
!
        INTEGER CLUMP_POS(ND) ! nth Shell's index
	INTEGER NCLUMPS	      ! # of shells
        INTEGER D_SH_IND      ! # of points from center to edge
!
! Local variables
!
	REAL(KIND=LDP)  T1, T2

	INTEGER CNT, INNER_IND, OUTER_IND
        INTEGER I, J, K
!
	INTEGER ERROR_LU
	INTEGER LUER
 	LOGICAL EDGE_FOUND
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	LUER = ERROR_LU()
!
! First, count the number of shell that exist in the wind. 
!
	CLUMP_POS(:) = 0
	J = 0
	DO I= 3, ND-2 !ND-5
	  IF(DENSITY(I) .GT. DENSITY(I-1) .AND. DENSITY(I) .GT. DENSITY(I+1))THEN
	    J = J + 1
	    CLUMP_POS(J)=I
	  END IF
	END DO

	IF( J .EQ. 0 )THEN
	   DO I=1,ND
	     PRINT*,I,DLOG10(DENSITY(I))
	   END DO
	   WRITE(LUER,*)'ERROR IN DETERM_CLUMP_POS_V3:'
	   WRITE(LUER,*)'   SHELLS NOT FOUND IN WIND!!'
	   WRITE(LUER,*)'   CHECK IF USING SHELL MODEL'
	   STOP
	END IF

	NCLUMPS = J
!
! Next, we find the number of point from its center to its edge; this can be done 
! by looking at its inner side where a minimum exist. As such, D_SH_IND is the
! number of points from shell center to this minimum. Since all shells share the
! same profile, and have the same number points across, this should be fine. 
!

	D_SH_IND = 1
	DO I = 1,NCLUMPS
	   CNT = 0
 	   EDGE_FOUND = .FALSE.
	   J =  CLUMP_POS(I)
!
 	   DO WHILE(.NOT. EDGE_FOUND)
!	   DO WHILE(CNT .LE. ND  )
	      CNT=CNT+1
	      IF( DENSITY( J + CNT ) .LT. DENSITY( (J + 1) + CNT ) .AND. 
	1         DENSITY( J + CNT ) .LT. DENSITY( (J - 1) + CNT ) )THEN
!
	         EDGE_FOUND = .TRUE.
		 D_SH_IND = MAX(D_SH_IND,CNT)
		 EXIT
	      END IF
	      
!
	      IF(CNT .GT. ND)THEN
		 WRITE(LUER,*)' '
	 	 WRITE(LUER,*)'ERROR IN DETERM_CLUMP_POS_V3:'
		 WRITE(LUER,*)'   There are no edge on this shells in wind? '
		 WRITE(LUER,*)'   Shell position = ', I
	 	 WRITE(LUER,*)'   index position = ', J
		 STOP
	      END IF
	   END DO
	END DO
!
	D_SH_IND = D_SH_IND +2
	IF(D_SH_IND .EQ. 1)THEN
	   WRITE(LUER,*)' WARNING: D_SH_IND = 1 DETECTED!!'
	   WRITE(LUER,*)'          NCLUMPS might be equal to 0, causing D_SH_IND = 1 '
	END IF
!
	RETURN
	END SUBROUTINE
