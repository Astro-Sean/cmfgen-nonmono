!
! Routine to read in analytical fits to partial photoionization
! cross-sections. Data is primarily need to handle ionizations by
! X-rays, although all fits are included.
!
! As a consequnce, need to be carefull that photoionization routes
! are not doble counted.
!
! Fits are from the Verner atomic data site.
!
! Ref:
!     Verner, D., Yakovlev, D.G., A\&A, 1995???
!     Analytic fits for partial photoinization cross-sections.
!
	MODULE XRAY_DATA_MOD
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 16-Jun-2025 -- ATNO_MAX_X updated from 28 to 60
! Altered: 28-Jul-2024 -- ATNO_MAX_X updated from 28 to 30
!
	INTEGER, PARAMETER :: ATNO_MAX_X=60
	INTEGER, PARAMETER :: PQN_MAX_X=3
	INTEGER, PARAMETER :: ANG_MAX_X=2
!
! We have added _X to make the names more unique. The notation is that
! of Verney and Yakovlev.
!
	REAL(KIND=LDP) E_THRESH_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL(KIND=LDP) E_0_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL(KIND=LDP) SIG_0_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL(KIND=LDP) y_a_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL(KIND=LDP) y_w_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
	REAL(KIND=LDP) P_X(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
!
! Gives the number of electrons ultimately ejected by the photoionization.
! At present we assume this is 1, or 2.
!
	INTEGER  N_ED_EJ(ATNO_MAX_X,ATNO_MAX_X,PQN_MAX_X,0:ANG_MAX_X)
!
	LOGICAL XRAY_PHOT_RD_IN
!
	END MODULE XRAY_DATA_MOD
!
!
!
! Subroutine read in all available data from the file XRAY_PHOT_FITS.
!
	SUBROUTINE RD_XRAY_FITS(LU)
	USE SET_KIND_MODULE
	USE XRAY_DATA_MOD
	IMPLICIT NONE
!
! Altered 28-July-2024 : Can now read in all 30 elements.
!
	INTEGER LU
	INTEGER I,J,K,L,IOS
	INTEGER IZ,INE,PQN,ANG
	INTEGER LAST_ATNO_RD
!
	EXTERNAL ERROR_LU
	INTEGER LU_ER,ERROR_LU
	LOGICAL NEW_FORMAT
	CHARACTER(LEN=100) STRING
!
	XRAY_PHOT_RD_IN=.FALSE.
	LU_ER=ERROR_LU()
	OPEN(UNIT=LU,FILE='XRAY_PHOT_FITS',STATUS='OLD',ACTION='READ',
	1             IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error opening XRAY_PHOT_FITS file'
	    WRITE(LU_ER,*)'IOSTAT=',IOS
	    STOP
	  END IF
!
	  NEW_FORMAT=.FALSE.
	  STRING=' '
	  DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)')STRING
	    IF(INDEX(STRING,'Format date:') .NE. 0)NEW_FORMAT=.TRUE.
	  END DO
!
	  SIG_0_X(:,:,:,:)=0.0_LDP
	  N_ED_EJ(:,:,:,:)=0
!
	  DO WHILE(1 .EQ. 1)	
	     READ(STRING,*,IOSTAT=IOS)IZ,INE,PQN,ANG
	     IF(IOS .NE. 0)THEN
	       WRITE(LU_ER,*)'Error reading IZ,INE... from XRAY_PHOT_FITS'
	       WRITE(LU_ER,*)'IOSTAT=',IOS
	       WRITE(LU_ER,'(A)')TRIM(STRING)
	       STOP
	     END IF
	     IF(PQN .GT. PQN_MAX_X)GOTO 50
	     IF(ANG .GT. ANG_MAX_X)GOTO 50
	     LAST_ATNO_RD=IZ
!
	     IF(NEW_FORMAT)THEN
	       READ(STRING,*,IOSTAT=IOS)I,J,K,L,E_THRESH_X(IZ,INE,PQN,ANG),
	1                          E_0_X(IZ,INE,PQN,ANG),
	1                          SIG_0_X(IZ,INE,PQN,ANG),
	1                          Y_A_X(IZ,INE,PQN,ANG),
	1                          P_X(IZ,INE,PQN,ANG),
	1                          Y_W_X(IZ,INE,PQN,ANG),
	1                          N_ED_EJ(IZ,INE,PQN,ANG)
	       IF(IOS .NE. 0)THEN
	         WRITE(LU_ER,*)'Error reading photoionization data from XRAY_PHOT_FITS'
	         WRITE(LU_ER,*)'IOSTAT=',IOS
	         STOP
	       END IF
	     ELSE
	       READ(STRING,*,IOSTAT=IOS)I,J,K,L,E_THRESH_X(IZ,INE,PQN,ANG),
	1                          E_0_X(IZ,INE,PQN,ANG),
	1                          SIG_0_X(IZ,INE,PQN,ANG),
	1                          Y_A_X(IZ,INE,PQN,ANG),
	1                          P_X(IZ,INE,PQN,ANG),
	1                          Y_W_X(IZ,INE,PQN,ANG)
	       IF(IOS .NE. 0)THEN
	         WRITE(LU_ER,*)'Error reading photoionization data from XRAY_PHOT_FITS'
	         WRITE(LU_ER,*)'IOSTAT=',IOS
	         STOP
	       END IF
	       N_ED_EJ(IZ,INE,PQN,ANG)=1
	       IF(PQN .EQ. 1 .AND. INE .GT. 3)N_ED_EJ(IZ,INE,PQN,ANG)=2
	       IF(IZ .LE. 12 .AND. PQN .EQ. 2 .AND. INE .GT. 11 .AND. ANG .EQ. 0)N_ED_EJ(IZ,INE,PQN,ANG)=2
	       IF(PQN .EQ. 2 .AND. INE .GE. 12)N_ED_EJ(IZ,INE,PQN,ANG)=2
	     END IF
!
50	    CONTINUE
	    READ(LU,'(A)',END=100)STRING
	  END DO
100	  CONTINUE
	  IF(LAST_ATNO_RD .GT. ATNO_MAX_X)THEN
	    WRITE(LU_ER,*)' '
	    WRITE(LU_ER,*)('*',I=1,80)
	    WRITE(LU_ER,*)('*',I=1,80)
	    WRITE(LU_ER,*)'Error reading photoionization data from XRAY_PHOT_FITS' 
	    WRITE(LU_ER,*)'Insuffient atomic data'
	    WRITE(LU_ER,*)'ATNO_MAX_X=',ATNO_MAX_X
	    WRITE(LU_ER,*)'LAST_ATNO_RD=',LAST_ATNO_RD
	    WRITE(LU_ER,*)('*',I=1,80)
	    WRITE(LU_ER,*)('*',I=1,80)
	    WRITE(LU_ER,*)' '
	  END IF
	  XRAY_PHOT_RD_IN=.TRUE.
!
	RETURN
	END
