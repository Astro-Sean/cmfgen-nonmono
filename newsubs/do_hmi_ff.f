!
! Subroutine to compute the free-free opacity associated with neutral H.
! NB: H- is the ground state; neutal hydrogen is the ion.
!
!
!
	MODULE HMI_FF_DATA
	IMPLICIT NONE
!
	INTEGER, SAVE :: NT
	INTEGER, SAVE :: NNU
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: CROSS(:,:)
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: LOG_T_TAB(:)
	REAL(KIND=LDP), SAVE, ALLOCATABLE :: LOG_NU_TAB(:)
!
	REAL(KIND=LDP) LOG_NU
	REAL(KIND=LDP) LOG_T
	INTEGER T_INDX
	INTEGER NU_INDX
!
	END MODULE HMI_FF_DATA
!
	SUBROUTINE DO_HMI_FF(ETA,CHI,ION_DEN,ED,TEMP,EMHNUKT,CONT_FREQ,LUIN,ND)
	USE SET_KIND_MODULE
	USE HMI_FF_DATA
	IMPLICIT NONE
!
! Created 26-Jun-2015 (Comments added 16-Aug-2015; cur_hmi)
!
	INTEGER ND
	INTEGER LUIN
!
	REAL(KIND=LDP) ETA(ND)
	REAL(KIND=LDP) CHI(ND)
	REAL(KIND=LDP) ION_DEN(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) TEMP(ND)
	REAL(KIND=LDP) EMHNUKT(ND)
	REAL(KIND=LDP) CONT_FREQ
!
	REAL(KIND=LDP) T1,T2,T3
	INTEGER I,J
	CHARACTER(LEN=200) STRING
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
	REAL(KIND=LDP) BOLTZMANN_CONSTANT
	EXTERNAL BOLTZMANN_CONSTANT
!
! Constants for opacity etc.
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	IF(FIRST_TIME)THEN
	  WRITE(6,*)'Opening H- free-free file'
	  OPEN(UNIT=LUIN,FILE='HMI_FF',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Format date') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'20-Jun-2015') .EQ. 0)THEN
	    WRITE(6,*)'Invalid format date when reading H- free-free data'
	    WRITE(6,*)TRIM(STRING)
	    STOP
	  END IF
	  READ(LUIN,'(A)')STRING
!
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'!NTHETA') .NE. 0)THEN
	    READ(STRING,*)NT
	  ELSE
	    WRITE(6,*)'Error -- NTHETA not found when reading H- free-free data'
	    WRITE(6,*)TRIM(STRING)
	    STOP
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'!NLAM') .NE. 0)THEN
	    READ(STRING,*)NNU
	  ELSE
	    WRITE(6,*)'Error -- NLAM not found when reading H- free-free data'
	    WRITE(6,*)TRIM(STRING)
	    STOP
	  END IF
	  WRITE(6,'(A,I3,4X,A,I3)')' NTHETA=',NT,'NLAM=',NNU
!
	  ALLOCATE (LOG_T_TAB(NT))
	  ALLOCATE (LOG_NU_TAB(NNU))
	  ALLOCATE (CROSS(NT,NNU))
!	
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LUIN,'(A)')STRING
	  END DO
	  WRITE(6,'(A)')STRING
!
! In the table, T is tabulated as 5040/T, is increasing with index.
! In the table, the wavelenth is in Angstroms.
! We change to T, NU(10^15Hz) for table axes -- both monotonically increase with index.
! Note: We read the T axis i backwards to ensure monotonically increasing.
! Recall: T is in units of 10^4 K, NU in units of 10^15 Hz.
!
! The cross-section tabulated is per hydrodegn atom per unit electron pressures (kT.Ne)
! and need to be multipled by 10^{-26}. There is a factor of 10^{10} to keep R.CHI unitless.
!
	  READ(STRING,*)(LOG_T_TAB(I),I=NT,1,-1)
	  LOG_T_TAB=LOG(0.504_LDP/LOG_T_TAB)
	  I=3; CALL WRITV_V2(LOG_T_TAB,NT,I,'LOG(T/10^4)_TAB',6)
	  DO J=1,NNU
	    READ(LUIN,*)LOG_NU_TAB(J),(CROSS(I,J),I=NT,1,-1)
	    LOG_NU_TAB(J)=LOG(2997.94_LDP/LOG_NU_TAB(J))	!Convert from Ang to 10^15Hz
	  END DO
	  T1=1.0E+04_LDP*1.0E-16_LDP*BOLTZMANN_CONSTANT()
	  CROSS=LOG(CROSS*T1)
	  FIRST_TIME=.FALSE.
	  I=3; CALL WRITV_V2(LOG_NU_TAB,NNU,I,'LOG(NU/10^15 Hz)_TAB',6)
	  WRITE(6,'(A)')' '
!
	END IF
!
	LOG_NU=LOG(CONT_FREQ)
	IF(LOG_NU .GT. LOG_NU_TAB(NNU))THEN
	  NU_INDX=NNU-1
	ELSE IF(LOG_NU .LT. LOG_NU_TAB(1))THEN
	  NU_INDX=1
	ELSE
	  NU_INDX=1
	  DO WHILE(1 .EQ. 1)
	    IF(LOG_NU .LE. LOG_NU_TAB(NU_INDX+1))EXIT
	    NU_INDX=NU_INDX+1
	  END DO
	END IF
!
	DO I=1,ND
	  LOG_T=LOG(TEMP(I))
	  IF(LOG_T .GE. LOG_T_TAB(NT))THEN
	    T_INDX=NT-1
	  ELSE IF(LOG_T .LE. LOG_T_TAB(1))THEN
	    T_INDX=1
	  ELSE
	    T_INDX=1
	    DO WHILE(1 .EQ. 1)
	      IF(LOG_T .LT. LOG_T_TAB(T_INDX+1))EXIT
	      T_INDX=T_INDX+1
	    END DO
	  END IF
!
! The tabulated cross-section already contains the free-free cross-section.  We need
! to multiply by TEMP(I) as its per unit electron pressurse. The constants were incorporated
! into the cross-sectins earlier.
!
	  T1=(LOG_NU-LOG_NU_TAB(NU_INDX))/(LOG_NU_TAB(NU_INDX+1)-LOG_NU_TAB(NU_INDX))
	  T2=(1.0_LDP-T1)*CROSS(T_INDX,NU_INDX)+T1*CROSS(T_INDX,NU_INDX+1)
	  T3=(1.0_LDP-T1)*CROSS(T_INDX+1,NU_INDX)+T1*CROSS(T_INDX+1,NU_INDX+1)
!
	  T1=(LOG_T-LOG_T_TAB(T_INDX))/(LOG_T_TAB(T_INDX+1)-LOG_T_TAB(T_INDX))
	  T1=TEMP(I)*EXP( (1.0_LDP-T1)*T2+T1*T3 )*ED(I)*ION_DEN(I)
!	  WRITE(6,'(2I4,7ES14.4)')T_INDX,NU_INDX,TEMP(I),ED(I),ION_DEN(I),CHI(I),T1,T2,T3
	  CHI(I)=CHI(I)+T1
	  ETA(I)=ETA(I)+T1*TWOHCSQ*(CONT_FREQ**3)*EMHNUKT(I)/(1.0_LDP-EMHNUKT(I))
!
	END DO
!
	RETURN
	END
!
	SUBROUTINE DO_HMI_FF_COOL(FF,ION_DEN,ED,TEMP,BPHOT_CR,JPHOT_CR,CONT_FREQ,ND)
	USE SET_KIND_MODULE
	USE HMI_FF_DATA
	IMPLICIT NONE
!
	INTEGER ND
	REAL(KIND=LDP) FF(ND)
	REAL(KIND=LDP) ION_DEN(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) TEMP(ND)
	REAL(KIND=LDP) BPHOT_CR(ND)
	REAL(KIND=LDP) JPHOT_CR(ND)
	REAL(KIND=LDP) CONT_FREQ
!
	REAL(KIND=LDP) T1,T2,T3
	INTEGER I,J
!
	LOG_NU=LOG(CONT_FREQ)
	IF(LOG_NU .GT. LOG_NU_TAB(NNU))THEN
	  NU_INDX=NNU-1
	ELSE IF(LOG_NU .LT. LOG_NU_TAB(1))THEN
	  NU_INDX=1
	ELSE
	  NU_INDX=1
	  DO WHILE(1 .EQ. 1)
	    IF(LOG_NU .LE. LOG_NU_TAB(NU_INDX+1))EXIT
	    NU_INDX=NU_INDX+1
	  END DO
	END IF
!
	DO I=1,ND
	  LOG_T=LOG(TEMP(I))
	  IF(LOG_T .GE. LOG_T_TAB(NT))THEN
	    T_INDX=NT-1
	  ELSE IF(LOG_T .LE. LOG_T_TAB(1))THEN
	    T_INDX=1
	  ELSE
	    T_INDX=1
	    DO WHILE(1 .EQ. 1)
	      IF(LOG_T .LT. LOG_T_TAB(T_INDX+1))EXIT
	      T_INDX=T_INDX+1
	    END DO
	  END IF
!
! The tabulated cross-section already contains the free-free cross-section.  We need to
! multiply by TEMP(I) as its per unit electron pressure. The constants were incorporated
! into the cross-sectins earlier.
!
	  T1=(LOG_NU-LOG_NU_TAB(NU_INDX))/(LOG_NU_TAB(NU_INDX+1)-LOG_NU_TAB(NU_INDX))
	  T2=(1.0_LDP-T1)*CROSS(T_INDX,NU_INDX)+T1*CROSS(T_INDX,NU_INDX+1)
	  T3=(1.0_LDP-T1)*CROSS(T_INDX+1,NU_INDX)+T1*CROSS(T_INDX+1,NU_INDX+1)
!
	  T1=(LOG_T-LOG_T_TAB(T_INDX))/(LOG_T_TAB(T_INDX+1)-LOG_T_TAB(T_INDX))
	  T1=TEMP(I)*EXP( (1.0_LDP-T1)*T2+T1*T3 )*ED(I)*ION_DEN(I)
!
! The constant in T2 is 4PI x 1.0E-10.
!
	  FF(I)=FF(I)+1.256637061E-09_LDP*T1*(BPHOT_CR(I)-JPHOT_CR(I))
!
	END DO
!
	RETURN
	END
