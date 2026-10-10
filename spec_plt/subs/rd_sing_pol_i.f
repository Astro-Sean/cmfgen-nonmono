	SUBROUTINE RD_SING_POL_I(FLUX,OBSF,NOS,NOS_MAX,IREC,FILENAME,LU,IOS)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered: 17-Nov-2021: Change routine name from RD_SNG_POL_I to RD_SING_POL_I
!
	INTEGER NOS
	INTEGER NOS_MAX
	INTEGER IOS
	INTEGER IREC
	INTEGER LU
	REAL(KIND=LDP) FLUX(NOS)
	REAL(KIND=LDP) OBSF(NOS)
	CHARACTER(LEN=*) FILENAME
!
	REAL(KIND=LDP) T1
	INTEGER I,J
	INTEGER NBETA
	CHARACTER(LEN=132) STRING
!
	IOS=0
	OPEN(UNIT=LU,FILE=TRIM(FILENAME),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error opening polarization file: ',TRIM(FILENAME)
	  IOS=100
	  RETURN
	END IF
!
        DO WHILE( INDEX(STRING,'!Number of polar angles') .EQ. 0)
          READ(LU,'(A)')STRING
        END DO
        READ(STRING,*)NBETA
!
	DO WHILE( INDEX(STRING,'!Number of frequencies [NOS]') .EQ. 0)
         READ(LU,'(A)')STRING
        END DO
        READ(STRING,*)NOS
	IF(NOS .GT. NOS_MAX)THEN
	  WRITE(6,*)'Arrays are too small: (NOS, NOS_MAX)',NOS,NOS_MAX
	  CLOSE(UNIT=LU)
	  IOS=1000
	END IF
!
	DO WHILE( INDEX(STRING,'Polar Coordinates - [BETA]') .EQ. 0)
          READ(LU,'(A)')STRING
        END DO
        READ(LU,*)(OBSF(I),I=1,NBETA)
	WRITE(6,*)'Viewing angles are as folows:'
	T1=180.0_LDP/ACOS(-1.0_LDP)
	DO I=1,NBETA
	  WRITE(6,'(2X,I3,3X,ES12.4,3X,F8.2)')I,OBSF(I),OBSF(I)*T1
	END DO
!
        DO WHILE( INDEX(STRING,'Observers Frequencies') .EQ. 0)
          READ(LU,'(A)')STRING
        END DO
        READ(LU,*)(OBSF(I),I=1,NOS)
!
! Now read in the NBETA line profiles.
!
        DO WHILE( INDEX(STRING,'Normalized Line profiles') .EQ. 0 .AND.
	1         INDEX(STRING,'I spectrum') .EQ. 0)
          READ(LU,'(A)')STRING
        END DO
        READ(LU,*)((FLUX(I),I=1,NOS),J=1,IREC)
	CLOSE(UNIT=LU)
!
	RETURN
	END
