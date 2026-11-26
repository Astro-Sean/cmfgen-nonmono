!
! This is a simple routine to create a smooth wind profile from one that has
! shells in it. This is done by linearly interpolating between the shell's inner
! and outer edges - essentially drawing a straight line from neighboring ICMs. 
! Note that the most inner shell (at Z=0) we can not do this since the inner edge 
! will be above the ICM level.
!
	SUBROUTINE SM_WIND_PROF_JH(DENSITY,R,ND,NCLUMPS,
	1                          D_SH_IND,CLUMP_POS,NO_SH_DEN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! From CMFGEN
!
	INTEGER ND
!
	REAL(KIND=LDP) DENSITY(ND)
	REAL(KIND=LDP) R(ND)
!
! Values to be returned
!
	REAL(KIND=LDP)  NO_SH_DEN(ND)      ! Density proflie w/o shells
!
	INTEGER NCLUMPS		   ! # of shells
        INTEGER D_SH_IND           ! # of points of
!
        INTEGER CLUMP_POS(NCLUMPS) ! nth Shell's index
!
! Local variables
!
	REAL(KIND=LDP)  T1, T2
!
	INTEGER CNT, INNER_IND, OUTER_IND
        INTEGER I, J, K
!
	INTEGER ERROR_LU
	INTEGER LUER
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	LUER = ERROR_LU()
	NO_SH_DEN(1:ND)=DENSITY(1:ND)	
!
! Do normal interpolation between ICM (1) except for the most inner shell (2);
! we need to extrapolate on most inner shell
!
!  /\           /\
!_/  \         /  \
!     \       /    \  
!      \_____/      \___
!  (2)          (1)
!            ^       ^ OUTER
!            |
!           INNER
!
	DO I = 1,NCLUMPS
	   INNER_IND = CLUMP_POS(I)+D_SH_IND
	   OUTER_IND = CLUMP_POS(I)-D_SH_IND
!	  
	   DO K = OUTER_IND,INNER_IND
	      T2 = ( R(K) - R(OUTER_IND) )/( R(INNER_IND) - R(OUTER_IND) )
	      NO_SH_DEN(K) = T2*DENSITY(INNER_IND) + ( 1.0 - T2 )*DENSITY(OUTER_IND)	
	   END DO
	END DO
!
	RETURN
	END SUBROUTINE
