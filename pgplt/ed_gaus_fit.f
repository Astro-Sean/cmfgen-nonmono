!
! Subroutine designed to edit the modified Gaussian data. New Gaussians may also
! be added. After each edit, the full Gaussian list is output for inspection.
!
	SUBROUTINE ED_GAUS_FIT
	USE SET_KIND_MODULE
	USE GAUS_FIT_DATA
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 01-Feb-2022          !Now allow a maximum of 20 new Gaussians to be added.
! Created:  -Sep-2007		!Author: D. J. Hillier
!
	REAL(KIND=LDP), ALLOCATABLE :: TMP_PAR(:)
	INTEGER K			!Array index
	INTEGER IP			!# of Gaussian
	LOGICAL SCRAP
!
	IF( NUM_GAUS .NE. (NG_PAR-2)/4 )THEN
	  WRITE(6,*)'Error in ED_GAUS_FIT: NM_GAUS inconsistent with NG_PAR'
	  WRITE(6,*)'    NUM_GAUS=',NUM_GAUS
	  WRITE(6,*)'      NG_PAR=',NG_PAR
	  WRITE(6,*)'4*NUM_GAUS+2=',4*NUM_GAUS+2
	  STOP
	END IF
!
! We allow a maximum of 20 new Gaussians to be added.
!
	ALLOCATE(TMP_PAR(NG_PAR_MAX)); TMP_PAR(:)=0.0D0
	TMP_PAR(1:NG_PAR)=PAR(1:NG_PAR)
!
	CALL GEN_IN(TMP_PAR(1),'Y value at begining of fit range')
	CALL GEN_IN(TMP_PAR(2),'Baseline slope')
!
! Loop edit section until finished. NB: For historical reasons, the third parameter
! is the height. Howver we write it out second, since it is more important for
! seeing the importance of the Gaussian component.
!
	DO WHILE(1 .EQ. 1)
!
	  IP=0
	  WRITE(6,'(A,4(8X,A))')'Index','  Lambda','  Height','Sigma(a)','Exponent'
	  DO K=3,NG_PAR,4
	    WRITE(6,'(I5,4ES16.6)')1+(K-3)/4,TMP_PAR(K),TMP_PAR(K+2),TMP_PAR(K+1),TMP_PAR(K+3)
	  END DO
	  IF(NG_PAR .EQ. NG_PAR_MAX)THEN
	    WRITE(6,*)'You have reached the maximum number of Gaussians that you can add'
	    WRITE(6,*)'Fit the data, and renter edit routine (ED_GAUS_FIT)'
	  END IF
	  SCRAP=.FALSE.; CALL GEN_IN(SCRAP,'Discard all Gaussians?')
	  IF(SCRAP)THEN
	    NUM_GAUS=0
	    WRITE(6,*)'Enter 0 for next option to exit this routine'
	    WRITE(6,*)'You can still add a set of Gaussian parameters by hand'
	    IP=0
	  END IF
!
	  CALL GEN_IN(IP,'Gaussian data to edit: 0 to exit, -ve # to delete that Gaussian')
	  IF(IP .EQ. 0)EXIT
	  IF(IP .LT. 0 .AND. ABS(IP) .LE. NUM_GAUS)THEN
	    IP=-IP
	    K=3+4*(IP-1)
	    TMP_PAR(K:NG_PAR-4)=TMP_PAR(K+4:NG_PAR)
	    NUM_GAUS=NUM_GAUS-1
	    NG_PAR=2+NUM_GAUS*4
!
	    WRITE(6,'(A,4(8X,A))')'Index','  Lambda','  Height','Sigma(a)','Exponent'
	    DO K=3,NG_PAR,4
	      WRITE(6,'(I5,4ES16.6)')1+(K-3)/4,TMP_PAR(K),TMP_PAR(K+2),TMP_PAR(K+1),TMP_PAR(K+3)
	    END DO
!
! Can now edit the  Gaussian.
!
	  ELSE IF(IP .GT. 0)THEN
	    IP=MIN(NUM_GAUS+1,IP)
	    K=2+(IP-1)*4+1
	    CALL GEN_IN(TMP_PAR(K),'Central wavlength of Gaussian')
	    IF(TMP_PAR(K+2) .EQ. 0.0D0)TMP_PAR(K+2)=-0.2D0
	    CALL GEN_IN(TMP_PAR(K+2),'Offset from continuum (-ve for absorption)')
            IF(TMP_PAR(K+1) .EQ. 0.0D0)TMP_PAR(K+1)=0.2D0
	    CALL GEN_IN(TMP_PAR(K+1),'Sigma of Gassian')
	    IF(TMP_PAR(K+3) .EQ. 0.0D0)TMP_PAR(K+3)=2.0D0
	    CALL GEN_IN(TMP_PAR(K+3),'Guassian exponent')
	    NUM_GAUS=MAX(NUM_GAUS,IP)
	    NG_PAR=2+NUM_GAUS*4
	  END IF
	END DO
!
! Clean up and save the modified Gaussian parameters.
!
	NG_PAR=2+NUM_GAUS*4
	PAR(1:NG_PAR)=TMP_PAR(1:NG_PAR)
	DEALLOCATE (TMP_PAR)
!	
	RETURN
	END
