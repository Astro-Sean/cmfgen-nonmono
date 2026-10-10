!
! Set of two routines to READ in desriptor populations.
! The Hydrogen and Helium populations aren now READ in from their
! own file.
!
	SUBROUTINE RD_RVTJ_PARAMS_V4(RMDOT,LUM,ABUNDH,TIME,NAME_CONV,
	1                  ND,NC,NP,FORMAT_DATE,FILNAME,LUIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 24-Mar-2023 : Simplified reading of variables after temperture.
!                         Changed to read MEAN_ABS opacity (20-Jan-2023).
!                         Should be backwards compatable.
!                         Not returned (done by separte read in DISPGEN).
! Altered 16-Aug-2019 : Changed to V4 (to maintain consistency with other routine in this file).
! Altered 15-Jun-2000 : Naming convention inserted.
!
	INTEGER ND,NC,NP,LUIN
	REAL(KIND=LDP) RMDOT,LUM,ABUNDH
	CHARACTER(LEN=*) TIME,FILNAME,NAME_CONV
	CHARACTER(LEN=*) FORMAT_DATE
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local Variables
!
	CHARACTER*11 LOC_FORMAT_DATE,PRODATE
	INTEGER NCF
	LOGICAL RD_FIX_T
!
	OPEN(UNIT=LUIN,FILE=FILNAME,STATUS='OLD',ACTION='READ')
!
	LOC_FORMAT_DATE='08-JAN-1996'
	READ(LUIN,'(T30,A11)')FORMAT_DATE
	CALL SET_CASE_UP(FORMAT_DATE,0,0)
	IF(   TRIM(FORMAT_DATE) .NE. '15-JUN-2000' .AND.
	1     TRIM(FORMAT_DATE) .NE. '10-NOV-2009' .AND.
	1     TRIM(FORMAT_DATE) .NE. '15-AUG-2019' .AND.
	1                      TRIM(FORMAT_DATE) .NE. LOC_FORMAT_DATE)THEN
	  WRITE(ERROR_LU(),*)'Wrong format date : RVTJ read failure'
	  WRITE(ERROR_LU(),*)'Subroutine called is RD_ASC_RVTJ_V4'
	  WRITE(ERROR_LU(),*)'Subroutine date 1 is: ',LOC_FORMAT_DATE
	  WRITE(ERROR_LU(),*)'Subroutine date 2 is: ','15-JUN-2000'
	  WRITE(ERROR_LU(),*)'Subroutine date 3 is: ','10-NOV-2009'
	  WRITE(ERROR_LU(),*)'Subroutine date 3 is: ','15-Aug-2019'
	  WRITE(ERROR_LU(),*)'File format date is: ',FORMAT_DATE
	  STOP
	END IF
	READ(LUIN,'(T30,A20)')TIME
	READ(LUIN,'(T30,A11)')PRODATE
	READ(LUIN,'(T30,BN,BZ,I5)')ND
	READ(LUIN,'(T30,BN,BZ,I5)')NC
	READ(LUIN,'(T30,BN,I5)')NP
	READ(LUIN,'(T30,BN,I5)')NCF
!
	READ(LUIN,'(T30,1PE12.5)')RMDOT
	READ(LUIN,'(T30,1PE12.5)')LUM
	READ(LUIN,'(T30,1PE12.5)')ABUNDH
	READ(LUIN,'(T30,L1)')RD_FIX_T
	IF(FORMAT_DATE .EQ. '15-JUN-2000'
	1   .OR. FORMAT_DATE .EQ. '10-NOV-2009'
	1   .OR. FORMAT_DATE .EQ. '15-AUG-2019')THEN
	   READ(LUIN,'(T30,A)')NAME_CONV
	ELSE
	   NAME_CONV='X_FOR_I'
	END IF
!
	RETURN
	END
!
! 
!
! This routine reads in the vecors, and returns the number of HI levels.
!
	SUBROUTINE RD_RVTJ_VEC_V4(R,V,SIGMA,ED,T,
	1       TGREY,dE_RAD_DECAY,
	1       ROSS_MEAN,FLUX_MEAN,PLANCK_MEAN,
	1       POP_ATOM,POP_ION,
	1       MASS_DENSITY,CLUMP_FAC,FORMAT_DATE,
	1       ND,LUIN)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
!Altered 16-Aug-2019 - changed to V4. PLANCK_MEAN added to call.
!
	INTEGER ND
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND),ED(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) TGREY(ND),dE_RAD_DECAY(ND)
	REAL(KIND=LDP) POP_ATOM(ND),POP_ION(ND)
	REAL(KIND=LDP) MASS_DENSITY(ND),CLUMP_FAC(ND)
	REAL(KIND=LDP) ROSS_MEAN(ND),FLUX_MEAN(ND),PLANCK_MEAN(ND)
	CHARACTER(LEN=*) FORMAT_DATE
!
	INTEGER LUIN
!
! Local variables
!
	REAL(KIND=LDP) T1
	CHARACTER*80 STRING
	INTEGER I
!
! At present, all vectors, are assumed to be output in same order.
! Could be changed by using descriptor string if required. We perform
! checks to confirm sequencing.
!
	CALL CHK_STRING(STRING,LUIN,'Radius','RD_RVTJ')
	READ(LUIN,*)(R(I),I=1,ND)
!
	CALL CHK_STRING(STRING,LUIN,'Velocity','RD_RVTJ')
	READ(LUIN,*)(V(I),I=1,ND)
!
	CALL CHK_STRING(STRING,LUIN,'dlnV/dlnr-1','RD_RVTJ')
	READ(LUIN,*)(SIGMA(I),I=1,ND)
!
	CALL CHK_STRING(STRING,LUIN,'Electron','RD_RVTJ')
	READ(LUIN,*)(ED(I),I=1,ND)
!
	CALL CHK_STRING(STRING,LUIN,'Temperature','RD_RVTJ')
	READ(LUIN,*)(T(I),I=1,ND)
!
! This format will allow new variables to be added after the temperature
! without affecting the reading of RVTJ.
!
	TGREY=0.0_LDP; dE_RAD_DECAY=0.0_LDP; ROSS_MEAN=0.0_LDP
	FLUX_MEAN=0.0_LDP; PLANCK_MEAN=0.0_LDP
	STRING=' '
	DO WHILE(INDEX(STRING,'Atom Density') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'Grey temperature') .NE. 0)THEN
	    READ(LUIN,*)(TGREY(I),I=1,ND)
	  ELSE IF(INDEX(STRING,'Heating: radioactive decay') .NE.  0)THEN
	    READ(LUIN,*)(dE_RAD_DECAY(I),I=1,ND)
	  ELSE IF(INDEX(STRING,'Rosseland Mean Opacity') .NE. 0)THEN
	    READ(LUIN,*)(ROSS_MEAN(I),I=1,ND)
	  ELSE IF(INDEX(STRING,'Flux Mean Opacity') .NE. 0)THEN
	    READ(LUIN,*)(FLUX_MEAN(I),I=1,ND)
	  ELSE IF(INDEX(STRING,'Planck Mean Opacity')  .NE. 0)THEN
	    READ(LUIN,*)(PLANCK_MEAN(I),I=1,ND)
	  ELSE IF(INDEX(STRING,'Absorption Mean Opacity')  .NE. 0)THEN
	    READ(LUIN,*)(T1,I=1,ND)
	  END IF
	END DO
!
	IF(SUM(TGREY) .EQ. 0.0D0)WRITE(6,*)'Grey temperature not set'
	IF(SUM(dE_RAD_DECAY) .EQ. 0.0D0)WRITE(6,*)'Radioactive decay heating not set'
	IF(SUM(ROSS_MEAN) .EQ. 0.0D0)WRITE(6,*)'Rosseland mean opacity not set'
	IF(SUM(FLUX_MEAN) .EQ. 0.0D0)WRITE(6,*)'Flux mean opacity not set'
	IF(SUM(PLANCK_MEAN) .EQ. 0.0D0)WRITE(6,*)'Planck mean opacity not set'
!
! Already read in header containg 'Atom Density'
!
	READ(LUIN,*)(POP_ATOM(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Ion Density','RD_RVTJ')
	READ(LUIN,*)(POP_ION(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Mass Density','RD_RVTJ')
	READ(LUIN,*)(MASS_DENSITY(I),I=1,ND)
	CALL CHK_STRING(STRING,LUIN,'Clumping Factor','RD_RVTJ')
	READ(LUIN,*)(CLUMP_FAC(I),I=1,ND)
!
	RETURN
	END
