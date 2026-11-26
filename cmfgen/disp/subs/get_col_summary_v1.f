!
! General routine to return collsion rates (in the form of collsion strengths
! OMEGA) for a single depth in the atmosphere. Two methods are used:
!	(1) Values are read in from a file. This file will not be read
!                in for every depth UNLESS it has been corrupted.
!                These values are used unless unavailable.
!                Interpolation is done in the Log-Log plane.
!	 (2) Use the approximate formulae which express the collsion
!                strengths interms of the oscilator values.
!
	SUBROUTINE GET_COL_SUMMARY_V1(OMEGA,EDGE,EIN_A,STAT_WT,
	1                       LEVELNAME,ZION,NLEV,TEMP,FILE_NAME)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER NLEV
	REAL(KIND=LDP) OMEGA(NLEV,NLEV)
!
	REAL(KIND=LDP) EIN_A(NLEV,NLEV)		!Einstein A coefficient
	REAL(KIND=LDP) EDGE(NLEV)		!Ionization frequency (10^15 Hz)
	REAL(KIND=LDP) STAT_WT(NLEV)		!Statistical weight.
	REAL(KIND=LDP) ZION
	CHARACTER*(*) LEVELNAME(NLEV)
!
	REAL(KIND=LDP) TEMP	       		!Temperature (10^4 K)
	CHARACTER*(*) FILE_NAME
!
	INTEGER ERROR_LU
	LOGICAL SAME_N
	EXTERNAL PHOT_SUB,ERROR_LU
!
	REAL(KIND=LDP), PARAMETER :: THREE=3.0_LDP
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE
	PARAMETER (MAX_TRANS=50000)
	PARAMETER (MAX_TVALS=40)
	PARAMETER (MAX_TAB_SIZE=MAX_TVALS*MAX_TRANS)
	INTEGER, SAVE :: ID_LOW(MAX_TRANS)
	INTEGER, SAVE :: ID_UP(MAX_TRANS)
	INTEGER, SAVE :: ID_INDX(MAX_TRANS)
	REAL(KIND=LDP), SAVE :: T_TABLE(MAX_TVALS)
	REAL(KIND=LDP), SAVE :: OMEGA_TABLE(MAX_TAB_SIZE)
	REAL(KIND=LDP), SAVE :: OMEGA_SET,OMEGA_SCALE
	INTEGER, SAVE :: NUM_TRANS,NUM_TVALS
!
! Local variables.
!
	INTEGER I		!Used for lower level.
	INTEGER J		!Used for upper level.
	INTEGER PHOT_ID
	INTEGER L,K,LUER,IFAIL
	INTEGER NL,NUP
!
	REAL(KIND=LDP) ALPHA
	REAL(KIND=LDP) T1
	REAL(KIND=LDP), PARAMETER :: ZERO=0.0_LDP
!
! Read in the atomic data. The data is not read in if the filename matches
! the filename of the previusly read data.
!
	CALL GEN_OMEGA_RD_V2(OMEGA_TABLE,T_TABLE,
	1                      ID_LOW,ID_UP,ID_INDX,
	1                      OMEGA_SCALE,OMEGA_SET,
	1                      LEVELNAME,STAT_WT,NLEV,FILE_NAME,
	1                      NUM_TRANS,NUM_TVALS,
	1                      MAX_TRANS,MAX_TVALS,MAX_TAB_SIZE)
!
! 
!
	LUER=ERROR_LU()
	OMEGA(:,:)=0.0_LDP			!NLEV*NLEV
!
	IF(NUM_TRANS .EQ. 0)GOTO 1000
	IF(TEMP .LE. T_TABLE(1))THEN
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K))
	  END DO
	ELSE IF(TEMP .GE. T_TABLE(NUM_TVALS))THEN
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K)+NUM_TVALS-1)
	  END DO
	ELSE
	  L=2
	  DO WHILE(TEMP .GT. T_TABLE(L))
	    L=L+1
	  END DO
	  T1=LOG(T_TABLE(L)/T_TABLE(L-1))
	  DO K=1,NUM_TRANS
	    I=ID_LOW(K)
	    J=ID_UP(K)
!
! NB: ID_INDX(K) refers to OMEGA for T_TABLE(1). Thus
! NB: ID_INDX(K)+L-1 refers to OMEGA for T_TABLE(L).
!
	    ALPHA=OMEGA_TABLE(ID_INDX(K)+L-1)/
	1                    OMEGA_TABLE(ID_INDX(K)+L-2)
	    ALPHA=LOG(ALPHA)/T1
	    OMEGA(I,J)=OMEGA_TABLE(ID_INDX(K)+L-2)*
	1                         (TEMP/T_TABLE(L-1))**ALPHA
	  END DO
	END IF
1000	CONTINUE
!
	WRITE(6,*)'Temp=',TEMP
	CALL OUT_COL_SUMMARY()
!
	CONTAINS
!
	SUBROUTINE OUT_COL_SUMMARY()
	USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER, SAVE :: LU_CS
	LOGICAL, SAVE :: FIRST_OUT=.TRUE.
!
	INTEGER K
	CHARACTER(LEN=100) FMT
	CHARACTER(LEN=5) TAB_LOC
!
	IF(FIRST_OUT)THEN
	  FIRST_OUT=.FALSE.
	  CALL GET_LU(LU_CS,'OUT_COL_SUMMARY -- gen_omega_rd_v2')
	  OPEN(UNIT=LU_CS,STATUS='UNKNOWN',ACTION='WRITE',FILE='NEW_COLLISION_SUMMARY')
	  WRITE(LU_CS,*)' '
	  WRITE(LU_CS,*)' NT_RD is the number of transitions read in from the collisonal data file.'
	  WRITE(LU_CS,*)' NT is the number of transitions computed'
	  WRITE(LU_CS,*)' It may be larger than NT_RD is the collision rates ar split'
	  WRITE(LU_CS,*)' '
	  WRITE(LU_CS,*)' Tranisitions to the ground term without a collision rate are listed.'
	  WRITE(LU_CS,*)' For permitted transiiotns, collison rates will be stimated using their gf value.'
	  WRITE(LU_CS,*)' Missing data for forbiddentransition is more omportant as their values cannot be'
	  WRITE(LU_CS,*)'     reliably estimated.'
	  WRITE(LU_CS,*)' '
	  WRITE(LU_CS,'(1X,A,T20,2X,A,5X,A,3(7X,A))')'Data file','NT_RD','NT',
	1               '  Temp','T(min)','T(max)'
	END IF
!
	WRITE(LU_CS,*)' '
	IF(NUM_TRANS .EQ. 0)THEN
	  WRITE(LU_CS,'(/,1X,2A)')'No collisional data read: ',TRIM(FILE_NAME)
	ELSE
	  WRITE(LU_CS,'(1X,A,T20,I6,3ES14.4)')TRIM(FILE_NAME),NUM_TRANS,TEMP,T_TABLE(1),T_TABLE(NUM_TVALS)
	END IF
	WRITE(LU_CS,*)' '
	K=MAXVAL(LEN_TRIM(LEVELNAME))
	WRITE(TAB_LOC,'(I3)')2*K+5
	FMT='(1X,A,T'//TRIM(ADJUSTL(TAB_LOC))//',8X,A,7X,A,5X,A)'
	WRITE(LU_CS,FMT)'Transition','Lam(A)','Ein A','Gamma'
!
	FMT='(1X,3A,T'//TRIM(ADJUSTL(TAB_LOC))//',ES14.4,ES12.2,F10.4)'
	DO I=1,MIN(5,NLEV)
	  DO J=I+1,MIN(10,NLEV)
	    T1=2.9979E+03_LDP/(EDGE(I)-EDGE(J))
	    WRITE(LU_CS,FMT)TRIM(LEVELNAME(I)),'-',
	1              TRIM(LEVELNAME(J)),T1,EIN_A(J,I),OMEGA(I,J)
	  END DO
	END DO
	FLUSH(LU_CS)
!
	RETURN
!
	END SUBROUTINE OUT_COL_SUMMARY
	END SUBROUTINE GET_COL_SUMMARY_V1
