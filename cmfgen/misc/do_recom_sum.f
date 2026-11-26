	PROGRAM DO_RECOM_SUM
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered: 26-Jul-2024 -- Number of species increased to 28.
!
	INTEGER, PARAMETER :: NSPEC=28
	INTEGER, PARAMETER :: NION_MAX=21
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUOUT=10
!
	CHARACTER(LEN=12) SPECIES_ABR(NSPEC)
	CHARACTER(LEN=4) GEN_ION_ID(NION_MAX)
!
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: RECOM_RATE(:)
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) TMIN,TMAX
!
	INTEGER, PARAMETER :: NCHK=7
	REAL(KIND=LDP) T_REC_CHK(NCHK)
	REAL(KIND=LDP) RECOM_CHK(NCHK)
!
	INTEGER ND
	INTEGER I
	INTEGER ID
	INTEGER ISPEC
	INTEGER L,LST,LEND
	INTEGER K, KST, KEND
	INTEGER IOS
!
	LOGICAL SPEC_DONE
	LOGICAL FILE_EXISTS
	CHARACTER(LEN=100) FILE_NAME
	CHARACTER(LEN=100) STRING
!
	DATA T_REC_CHK/0.1_LDP,0.2_LDP,0.5_LDP,1.0_LDP,2.0_LDP,5.0_LDP,10.0_LDP/
!
! This is the same ID used in CMFGEN, and used to read files etc.
!
	DATA GEN_ION_ID /'0','I','2','III','IV','V',
	1                'SIX','SEV','VIII','IX','X','XI','XII',
	1                'XIII','XIV','XV','XSIX','XSEV','X8','X9','XX'/
!
	DATA SPECIES_ABR/'H','He','C','N','O','F','Ne',
	1         'Na','Mg','Al','Si','P','S','Cl','Ar',
	1         'K','Ca','Sc','Tk','V','Cr','Mn','Fe',
	1          'Co','Nk','Cu','Zn','Ba'/
!

	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')'This program requires the following files:'
	WRITE(6,'(A)')'    MODEL'
	WRITE(6,'(A)')'    RVTJ'
	WRITE(6,'(A)')'    XzVPRRR'
	WRITE(6,'(A)')'  '
	WRITE(6,'(A)')'Results will be output to RECOM_CHK'
	WRITE(6,'(A)')DEF_PEN
!
! Get number of depth points
!
	OPEN(UNIT=LUIN,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)EXIT
	    IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	      READ(STRING,*)ND
	      WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	      CLOSE(LUIN)
	      EXIT
	    END IF
	  END DO
	ELSE
	  WRITE(6,*)RED_PEN,'Unable to open file MODEL. IOS=',IOS
	  WRITE(6,*)DEF_PEN
	  STOP
	END IF
!
	ALLOCATE(T(ND))
	ALLOCATE(RECOM_RATE(ND))
	CALL RD_SING_VEC_RVTJ(T,ND,'Temperature','RVTJ',LUIN,IOS)
	IF(IOS .NE. 0)STOP
!
	TMIN=MINVAL(T);  TMAX=MAXVAL(T)
	WRITE(6,'(1X,A,1X,F7.3)')'T(min)=',TMIN
	WRITE(6,'(1X,A,1X,F7.3)')'T(max)=',TMAX
!
! Define the range of check temperatures covered by this model.
!
	KST=1; KEND=NCHK
	DO WHILE(TMIN .GT. T_REC_CHK(KST))
	  KST=KST+1
	  IF(KST .EQ. NCHK+1)THEN
	    WRITE(6,*)'Invalid temperature range for KST'
	    STOP
	  END IF
	END DO
	DO WHILE(TMAX .LT. T_REC_CHK(KEND))
	  KEND=KEND-1
	  IF(KEND .EQ. -1)THEN
	    WRITE(6,*)'Invalid temperature range for KEND'
	    STOP
	  END IF
	END DO
!
	OPEN(UNIT=LUOUT,FILE='RECOM_CHK',STATUS='UNKNOWN',IOSTAT=IOS,ACTION='WRITE')
	WRITE(LUOUT,'(A,T12,10F12.3)')'Species',(T_REC_CHK(K),K=KST,KEND)
	WRITE(LUOUT,'(A)')' '
!
! Loop over all species. Simple linear interpolation is used to get the
! recombination rate (which will include stimulated recombination).
!
	DO ISPEC=1,NSPEC
	  SPEC_DONE=.FALSE.
	  DO ID=1,NION_MAX
	    FILE_NAME=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(ID))//'PRRR'
	    INQUIRE(FILE=FILE_NAME,EXIST=FILE_EXISTS)
	    IF(FILE_EXISTS)THEN
	      SPEC_DONE=.TRUE.
	      OPEN(UNIT=LUIN,FILE=TRIM(FILE_NAME),STATUS='OLD',ACTION='READ')
	      DO L=1,ND,10
	        STRING=' '
	        DO WHILE(INDEX(STRING,'Radiative Recombination') .EQ. 0)
	          READ(LUIN,'(A)')STRING
	        END DO
	        LEND=MIN(L+9,ND)
	        READ(LUIN,*)(RECOM_RATE(I),I=L,LEND)
	      END DO
!
	     RECOM_CHK=0.0_LDP
	     DO K=KST,KEND
               DO I=1,ND-1
                 IF( (T_REC_CHK(K)-T(I))*(T(I+1)-T_REC_CHK(K)) .GE. 0 )THEN
                    T1=(T_REC_CHK(K)-T(I))/(T(I+1)-T(I))
                    RECOM_CHK(K)=(1.0_LDP-T1)*RECOM_RATE(I)+T1*RECOM_RATE(I+1)
                    EXIT
                 END IF
               END DO
	     END DO
	     WRITE(LUOUT,'(A,T12,10ES12.3)')TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(ID)),
	1                                   (RECOM_CHK(K),K=KST,KEND)
!
	    END IF
	  END DO
	  IF(SPEC_DONE)WRITE(LUOUT,'(A)')' '
	END DO
!
	STOP
	END
