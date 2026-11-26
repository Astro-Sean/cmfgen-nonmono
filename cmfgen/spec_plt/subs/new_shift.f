


	SUBROUTINE NEW_SHIFT(ND,ORIG_R,DEN,SHIFT_R,NCLUMPS,D_SH_IND,SHELL_LOC,
	1                    DIV_BY_SQ_DEN_SM,NEW_ND,NEW_R,NEW_DEN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! INPUTS
!
	INTEGER ND       
	REAL(KIND=LDP)  DEN(ND)      ! The Density  value being shifted
	REAL(KIND=LDP)  ORIG_R(ND)   ! Old R grid
	REAL(KIND=LDP)  SHIFT_R(ND)  ! Old R grid where shells are shifted

	INTEGER NCLUMPS
	INTEGER D_SH_IND
	INTEGER SHELL_LOC(NCLUMPS)

	INTEGER NEW_ND
	REAL(KIND=LDP)  NEW_R(NEW_ND)   ! New R grid where shifted shells are lin interp on

	LOGICAL DIV_BY_SQ_DEN_SM
!
! OUTPUTS
!
	REAL(KIND=LDP)  NEW_DEN(NEW_ND) ! New Density on NEW_R
!
! LOCAL VARIABLES 
!
	INTEGER I,J,K

	REAL(KIND=LDP)  SH_ONLY(ND)          ! Shell wind profile on old grid (ICM)
	REAL(KIND=LDP)  SM_WIND(ND)          ! Smooth wind profile on old grid

	REAL(KIND=LDP)  NEW_SH_ONLY(NEW_ND)  ! Shell wind profile on new grid
	REAL(KIND=LDP)  NEW_SM_WIND(NEW_ND)  ! Smooth wind profile on new grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	    CALL SM_WIND_PROF_JH(DEN,ORIG_R,ND,NCLUMPS,D_SH_IND,SHELL_LOC,SM_WIND)
!
	    IF(DIV_BY_SQ_DEN_SM)THEN
	      DO I=1,ND
	        SH_ONLY(I) = DEN(I)/(SM_WIND(I)*SM_WIND(I))
	      END DO
	    ELSE
	      DO I=1,ND
	        SH_ONLY(I) = DEN(I)/SM_WIND(I)
	      END DO
	    END IF
!
	    !CALL MON_INTERP(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
	    CALL MON_INTERP(NEW_SM_WIND,NEW_ND,1,NEW_R,NEW_ND,
	1                   SM_WIND,ND,SHIFT_R,ND)

!	    CALL LIN_INTERP(NEW_R,NEW_SM_WIND,NEW_ND,
!	1                   SHIFT_R,SM_WIND,ND)


	    CALL MON_INTERP(NEW_SH_ONLY,NEW_ND,1,NEW_R,NEW_ND,
	1                   SH_ONLY,ND,SHIFT_R,ND)

!            CALL LIN_INTERP(NEW_R,NEW_SH_ONLY,NEW_ND,
!	1                   SHIFT_R,SH_ONLY,ND)

	    IF(DIV_BY_SQ_DEN_SM)THEN
	      NEW_DEN(:) = NEW_SH_ONLY(:)*NEW_SM_WIND(:)*NEW_SM_WIND(:)
	    ELSE
	      NEW_DEN(:) = NEW_SH_ONLY(:)*NEW_SM_WIND(:)
	    END IF
!
	    RETURN
!
	END SUBROUTINE
