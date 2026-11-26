!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created
!
!     GENCOOL_SUM   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT Only the top rates are printed: Sorted using depths 1, 11, 21 etc.
!
	PROGRAM MOD_PRRR
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! ALtered 29-Jun-2024: Code now the ability to print out the recombination coefficient for each level.
! Altered 01-Mar-2023: PCUM and RCUM options installed.
! Altered 22-Sep-2022: Minor bug fix -- Recom. charge exchange was being output twice.
! Altered 10-Dec-2021: Improved output to ?_PRRR_SUM. Now much easier to read.
!                         Improved read in of header data.
! Altered 17-Nov-2021: Transfered from OSIRIS.
! Altered 22-Oct-2021: Major editing to improve optios (and earlier)
! Altered 27-May-2021: Small bug fixis (transferred from OSIRIS 21-Jul-20201).
! Altered 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered 08-Feb-2008: Extra terms (such as V term) sheck and output.
!
!
	REAL(KIND=LDP), ALLOCATABLE :: PHOT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RECOM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RECOM_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: PHOT_SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: COL_IR(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHG_IR(:)
	REAL(KIND=LDP), ALLOCATABLE :: COL_RR(:)
	REAL(KIND=LDP), ALLOCATABLE :: NT_IR_1E(:)
	REAL(KIND=LDP), ALLOCATABLE :: NT_IR_2E(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: XRAY_IR_2E(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RR_COEF(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHG_RR(:)
	REAL(KIND=LDP), ALLOCATABLE :: ADVEC_RR(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: TEMP(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: DI(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: XVEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: YVEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: TOT_SUM(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: YLEV(:)
	REAL(KIND=LDP), ALLOCATABLE :: XLEV(:)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,L
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	INTEGER IST,IEND,NLEV
	INTEGER NPRINT
	INTEGER NUM_HTDR_RATES
	INTEGER NREC
	INTEGER NHT
!
	REAL(KIND=LDP) T1,T2,T3,T4
	REAL(KIND=LDP) MIN_VAL
	LOGICAL FILE_OPEN
	LOGICAL NET_RECOM_PER_LEVEL
	LOGICAL DO_INDIV_RATES
	LOGICAL NORM
	LOGICAL RECOM_COEF 
	LOGICAL ABS_VALUE
	LOGICAL RD_NT_2E
	LOGICAL RD_XRAY_2E
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
        CHARACTER(LEN=2), PARAMETER :: FORMFEED=' '//CHAR(12)
!
	CHARACTER(LEN=132) TMP_STR
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) FILE_NAME
	CHARACTER(LEN=10) SPECIES
	CHARACTER(LEN=10) HTDR_SPECIES
	CHARACTER(LEN=5) OPTION
	CHARACTER(LEN=30) XLABEL
	CHARACTER(LEN=30) YLABEL
        CHARACTER(LEN=30) UC
        CHARACTER(LEN=30) FRMT
        EXTERNAL UC
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	ALLOCATE(PHOT_SUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(RECOM_SUM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(COL_IR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(CHG_IR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(COL_RR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(NT_IR_1E(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(NT_IR_2E(ND,2),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(XRAY_IR_2E(ND,2),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(RR_COEF(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(CHG_RR(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(ADVEC_RR(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE(R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(TEMP(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(DI(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE(TA(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(XVEC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(YVEC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(TOT_SUM(ND),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE(XVEC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(YVEC(ND),STAT=IOS)
!
	IOS=0
	OPEN(UNIT=20,FILE='RVTJ',STATUS='OLD',IOSTAT=IOS)
	  DO WHILE(INDEX(STRING,'Velocity (km/s)') .EQ. 0)
            READ(20,'(A)',IOSTAT=IOS)STRING
	  END DO
	  READ(20,*,IOSTAT=IOS)(V(I),I=1,ND)
	CLOSE(UNIT=20)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to get V(km/s) from RVTJ'
	  V(1:ND)=1.0
	END IF
!
	NET_RECOM_PER_LEVEL=.FALSE.
	DO_INDIV_RATES=.FALSE.
	RECOM_COEF=.FALSE.
	CALL GEN_IN(DO_INDIV_RATES,'Ouput recombination rate and photioization rate for each leve')
	CALL GEN_IN(RECOM_COEF,'Ouput recombination coef for for each leve')
	IF(.NOT. DO_INDIV_RATES)THEN
	   CALL GEN_IN(NET_RECOM_PER_LEVEL,'Ouput net recombination rate to each level?')
	END IF
!
	SPECIES='FeI'
1000	CONTINUE
!
	CHG_IR(1:ND)=0.0_LDP
	COL_IR(1:ND)=0.0_LDP
	CHG_RR(1:ND)=0.0_LDP
	COL_RR(1:ND)=0.0_LDP
	NT_IR_1E(1:ND)=0.0_LDP
	NT_IR_2E(1:ND,2)=0.0_LDP
	XRAY_IR_2E(1:ND,2)=0.0_LDP
	RECOM_SUM(1:ND)=0.0_LDP
	PHOT_SUM(1:ND)=0.0_LDP
	TOT_SUM(1:ND)=0.0_LDP
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(SPECIES,'File is assumed to be SPECIES//PRRR- EX to exit')
	  IF(UC(SPECIES) .EQ. 'EX')STOP
	  FILE_NAME=TRIM(SPECIES)//'PRRR'
	  OPEN(UNIT=20,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,'(2A,I5,A)')' Error -- unable to open file '//TRIM(FILE_NAME),' (IOS=',IOS,')'
	    WRITE(6,'(A)')' Check capatilzation of ION name'
	  END IF
	END DO
!
	STRING=' '
	DO WHILE(INDEX(STRING,'Photoionization Rate') .EQ. 0)
	  READ(20,'(A)')STRING
	END DO
	NLEV=-1
	DO WHILE(STRING .NE. ' ')
	   READ(20,'(A)')STRING
	   NLEV=NLEV+1
	END DO
	WRITE(6,'(A,I4)')' Number of levels in the model is',NLEV
	CLOSE(UNIT=20)
!
	IF(ALLOCATED(PHOT))DEALLOCATE(PHOT,RECOM)
	ALLOCATE (PHOT(ND,NLEV))
	ALLOCATE (RECOM(ND,NLEV))
	ALLOCATE(XLEV(NLEV),STAT=IOS)
	ALLOCATE(YlEV(NLEV),STAT=IOS)
	OPEN(UNIT=20,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
!
	IF(DO_INDIV_RATES .OR. NET_RECOM_PER_LEVEL)THEN
	  NPRINT=NLEV
	  CALL GEN_IN(NPRINT,'Number of levels to output for phot/red or net to each level?')
	END IF
	FILE_NAME=TRIM(FILE_NAME)//'_SUM'
	OPEN(UNIT=21,FILE=FILE_NAME,STATUS='UNKNOWN',ACTION='WRITE')
!
	STRING=' '
	FRMT='(1X,A,T21,A)'
	DO L=1,ND,10
	   RD_XRAY_2E=.FALSE.; RD_NT_2E=.FALSE.
	   IST=L; IEND=MIN(ND,IST+9)
	   WRITE(6,*)'Reading depths',IST,IEND
!
	   STRING=' '
	   WRITE(21,'(1X,A,T23,10(I11,'' ''))')'Depth',(I,I=IST,IEND)
	   DO WHILE(INDEX(STRING,'Photoionization Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Radius') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Radius(10^10cm)',TRIM(STRING)
	       READ(STRING,*)(R(I),I=IST,IEND)
	       WRITE(21,'(1X,A,T22,10ES12.4)')'V(km/s)',(V(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Temp(10^4K)',TRIM(STRING)
	       READ(STRING,*)(TEMP(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Electron Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'ED(cm^-3)',TRIM(STRING)
	       READ(STRING,*)(ED(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Ion Density') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Ion D.(cm^-3)',TRIM(STRING)
	       READ(STRING,*)(DI(I),I=IST,IEND)
	     END IF
	   END DO
	   WRITE(21,'(A)')' '
!
	   DO J=1,NLEV
	     READ(20,*)(PHOT(I,J),I=IST,IEND)
	     PHOT_SUM(IST:IEND)=PHOT_SUM(IST:IEND)+PHOT(IST:IEND,J)
	   END DO
	   WRITE(21,'(1X,A,T22,10ES12.4)')'Tot. Phot. IR',(PHOT_SUM(I),I=IST,IEND)
!
	   IF(DO_INDIV_RATES)THEN
	     DO J=1,NPRINT
	       WRITE(21,'(A,I4,A,T24,10ES12.4)')'  Phot(',J,')',(PHOT(I,J),I=IST,IEND)
	     END DO
	     IF(NPRINT .GT. 0 .AND. NPRINT .LT. NLEV)THEN
	       WRITE(TMP_STR,'(I3,A,I3,A)')NPRINT,'-',NLEV,')'
	       WRITE(21,'(A,A,T24,10ES12.4)')'  Phot(',TRIM(TMP_STR),(SUM(PHOT(I,NPRINT+1:NLEV)),I=IST,IEND)
	     END IF
	   END IF
!
	   DO WHILE(INDEX(STRING,'Recombination Rates') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Colisional Ionization Rate') .NE. 0)THEN
	        READ(20,'(A)')STRING
	        WRITE(21,FRMT)'Col. IR',TRIM(STRING)
	        READ(STRING,*)(COL_IR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Ionization Rate') .NE. 0)THEN
	        READ(20,'(A)')STRING
	        WRITE(21,FRMT)'Chg. TR. IR',TRIM(STRING)
	        READ(STRING,*)(CHG_IR(I),I=IST,IEND)
	     END IF
	   END DO
!
	   DO J=1,NLEV
	     READ(20,*)(RECOM(I,J),I=IST,IEND)
	     RECOM_SUM(IST:IEND)=RECOM_SUM(IST:IEND)+RECOM(IST:IEND,J)
	   END DO
!
	   WRITE(21,'(1X,A,T22,10ES12.4)')'Tot. Rad. Rec. R',(RECOM_SUM(I),I=IST,IEND)
	   IF(NET_RECOM_PER_LEVEL)THEN
	     DO J=1,NPRINT
	       WRITE(21,'(3X,A,I4,A,T24,10ES12.4)')'Net(',J,')',((RECOM(I,J)-PHOT(I,J)),I=IST,IEND)
	     END DO
	     IF(NPRINT .GT. 0 .AND. NPRINT .LT. NLEV)THEN
	       WRITE(TMP_STR,'(I3,A,I3,A)')NPRINT,'-',NLEV,')'
	       WRITE(21,'(3X,A,A,T24,10ES12.4)')'Phot(',TRIM(TMP_STR),
	1           ((SUM(RECOM(I,NPRINT+1:NLEV))-SUM(PHOT(I,NPRINT+1:NLEV))),I=IST,IEND)
	     END IF
	   ELSE IF(RECOM_COEF)THEN
	     DO J=1,NPRINT
	       WRITE(21,'(3X,A,I4,A,T22,10ES12.4)')'REC(',J,')',(RECOM(I,J)/ED(I)/DI(I),I=IST,IEND)
	     END DO
	     IF(NPRINT .GT. 0 .AND. NPRINT .LT. NLEV)THEN
	       WRITE(TMP_STR,'(I3,A,I3,A)')NPRINT,'-',NLEV,')'
	       WRITE(21,'(3X,A,A,T24,10ES12.4)')'Phot(',TRIM(TMP_STR),(SUM(RECOM(I,NPRINT+1:NLEV)/ED(I)/DI(I)),I=IST,IEND)
	     END IF
	   ELSE IF(DO_INDIV_RATES)THEN
	     DO J=1,NPRINT
	       WRITE(21,'(3X,A,I4,A,T22,10ES12.4)')'REC(',J,')',(RECOM(I,J),I=IST,IEND)
	     END DO
	     IF(NPRINT .GT. 0 .AND. NPRINT .LT. NLEV)THEN
	       WRITE(TMP_STR,'(I3,A,I3,A)')NPRINT,'-',NLEV,')'
	       WRITE(21,'(3X,A,A,T24,10ES12.4)')'Phot(',TRIM(TMP_STR),(SUM(RECOM(I,NPRINT+1:NLEV)),I=IST,IEND)
	     END IF
	   END IF
!
	   DO WHILE(INDEX(STRING,'Net Recombination Rate') .EQ. 0)
	     READ(20,'(A)')STRING
	     IF(INDEX(STRING,'Colisional Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Col. RR',TRIM(STRING)
	       READ(STRING,*)(COL_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Charge Transfer Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Chg TR RR',TRIM(STRING)
	       READ(STRING,*)(CHG_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Effective Advection Recombination Rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Eff. Adv. RR',TRIM(STRING)
	       READ(STRING,*)(ADVEC_RR(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'Non-Thermal Ionization') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       WRITE(21,FRMT)'Non therm. IR',TRIM(STRING)
	       READ(STRING,*)(NT_IR_1E(I),I=IST,IEND)
	     ELSE IF(INDEX(STRING,'NT ionization rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       J=1; IF(RD_NT_2E)J=2; RD_NT_2E=.TRUE.
	       WRITE(21,FRMT)'Net NT IR',TRIM(STRING)
	       READ(STRING,*)(NT_IR_2E(I,J),I=IST,IEND) 
!
! We also check the old X-ray format.
!
	     ELSE IF(INDEX(STRING,'Net X-ray recombination rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       J=1; IF(RD_XRAY_2E)J=2; RD_XRAY_2E=.TRUE.
	       READ(STRING,*)(XRAY_IR_2E(I,J),I=IST,IEND) 
	       XRAY_IR_2E(IST:IEND,J)=-XRAY_IR_2E(IST:IEND,J)
	       WRITE(21,'(A,T21,10ES12.4)')' Net Xray IR  ',XRAY_IR_2E(IST:IEND,J)
!
	     ELSE IF(INDEX(STRING,'Net X-ray ionization rate') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       J=1; IF(RD_XRAY_2E)J=2; RD_XRAY_2E=.TRUE.
	       WRITE(21,FRMT)'Net Xray IR',TRIM(STRING)
	       READ(STRING,*)(XRAY_IR_2E(I,J),I=IST,IEND)
!
	     END IF
	   END DO
!
	   READ(20,'(A)')STRING
	   WRITE(21,'(A)')' '
	   WRITE(21,FRMT)'Net dN/dt(%)',TRIM(STRING)
	   FLUSH(21)
!
	   STRING=' '
	   DO WHILE(STRING .EQ. ' ')
	     READ(20,'(A)')STRING
	   END DO
	   IF(INDEX(STRING,'Radiative Recombination Coefficient') .NE. 0)THEN
	     READ(20,'(A)')STRING
	     WRITE(21,FRMT)'Rad Rec alpha',TRIM(STRING)
	     READ(STRING,*)(RR_COEF(I),I=IST,IEND)
	   END IF
!
	   DO I=IST,IEND
	     TOT_SUM(I)=( PHOT_SUM(I)+RECOM_SUM(I) +
	1             COL_IR(I)+CHG_IR(I) +
	1             COL_RR(I)+CHG_RR(I) +
	1             ABS(XRAY_IR_2E(I,1)) + ABS(XRAY_IR_2E(I,2)) +
	1             ABS(NT_IR_1E(I)) +
	1             ABS(NT_IR_2E(I,1)) + ABS(NT_IR_2E(I,2)) +
	1             ABS(ADVEC_RR(I)) )/2.0_LDP
	   END DO
	   WRITE(21,'(A,T22,10ES12.4)')' Rate into/out of',(TOT_SUM(I),I=IST,IEND)
	   WRITE(21,'(/,/,A)')' '
!
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)',END=200)STRING
	     IF(STRING(1:1) .EQ. '1' .OR. STRING(1:1) .EQ. FORMFEED)THEN
	       EXIT
	       WRITE(21,'(A)')FORMFEED
	     END IF
	     WRITE(21,'(A)')TRIM(STRING)
	   END DO
	   FLUSH(21)
	END DO
200	CONTINUE
!
! Irrespective of the name, the terms now indicate the rate into the  ion.
! Thus -ve values indicate a loss.
!
!
	DO I=1,ND
	   PHOT_SUM(I)=-PHOT_SUM(I)
	   COL_IR(I)=-COL_IR(I)
	   CHG_IR(I)=-CHG_IR(I)
	   NT_IR_1E(I)=-NT_IR_1E(I)
	   NT_IR_2E(I,:)=-NT_IR_2E(I,:)
	   XRAY_IR_2E(I,:)=-XRAY_IR_2E(I,:)
!
	   RECOM_SUM(I)=RECOM_SUM(I)
	   COL_RR(I)=COL_RR(I)
	   CHG_RR(I)=CHG_RR(I)
	   ADVEC_RR(I)=ADVEC_RR(I)
	END DO
	YLABEL='Normalized rate'
	XVEC(1:ND)=V(1:ND)
	XLABEL='V(km/s)'
!
2000	CONTINUE
	OPTION='P'; WRITE(6,'(A)')' '
	CALL GEN_IN(OPTION,'Plot option XN, R, V, T, ED, EX(stop), NS (new species), NRR, RATES')
	OPTION=UC(OPTION)
	WRITE(6,*)'Option=',OPTION
	IF(OPTION(1:2) .EQ. 'XN')THEN
	  DO I=1,ND
	    XVEC(I)=I
	  END DO
	  XLABEL='Depth index'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'NS')THEN
	  GOTO 1000
	ELSE IF(OPTION .EQ. 'XR')THEN
	  XVEC(1:ND)=R(1:ND)/R(ND)
	  XLABEL='Radius(R/R\d*\u)         !(10\u10\d cm)'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'XV')THEN
	  XVEC(1:ND)=V(1:ND)
	  XLABEL='V(km/s)'
	  GOTO 2000
	ELSE IF(OPTION .EQ. 'XT')THEN
	  XVEC(1:ND)=TEMP(1:ND)
	  XLABEL='T(10\u4\dK)'
	ELSE IF(OPTION .EQ. 'XED')THEN
	  XVEC(1:ND)=ED(1:ND)
	  XLABEL='Ne(cm\u-3\d)'
	ELSE IF(OPTION(1:2) .EQ. 'EX')THEN
	  STOP
!
        ELSE IF(OPTION(1:2) .EQ. 'SP')THEN
	  CALL GEN_IN(K,'Number of levels to plot individually')
	  DO J=1,K
            DO I=1,ND
              YVEC(I)=(PHOT(I,J)-RECOM(I,J))/RECOM_SUM(I)
              YVEC(I)=PHOT(I,J)/RECOM_SUM(I)
            END DO
            CALL DP_CURVE(ND,XVEC,YVEC)
          END DO
!
	  CALL GEN_IN(K,'Final level to sum')
          DO J=K+1,L
            DO I=1,ND
!              YVEC(I)=YVEC(I)+(PHOT(I,J)-RECOM(I,J))/RECOM_SUM(I)
              YVEC(I)=YVEC(I)+PHOT(I,J)/RECOM_SUM(I)
            END DO
          END DO
          CALL DP_CURVE(ND,XVEC,YVEC)
!
	ELSE IF(OPTION(1:4) .EQ. 'PCUM')THEN
	  CALL GEN_IN(J,'Depth to examine')
	  YLEV=0.0_LDP; YLEV(1)=PHOT(J,1); XLEV(1)=1
	  DO I=2,NLEV
	    XLEV(I)=I
	    YLEV(I)=YLEV(I-1)+PHOT(J,I)
	  END DO
	  YLEV=-YLEV/PHOT_SUM(J)
	  YLABEL='\gs\di\un Phot. Rate.'
	  CALL DP_CURVE(NLEV,XLEV,YLEV)
!
	ELSE IF(OPTION(1:4) .EQ. 'RCUM')THEN
	  CALL GEN_IN(J,'Depth to examine')
	  YLEV=0; YLEV(1)=RECOM(J,1); XLEV(1)=1
	  DO I=2,NLEV
	    XLEV(I)=I
	    YLEV(I)=YLEV(I-1)+RECOM(J,I)
	  END DO
	  YLEV=YLEV/RECOM_SUM(J)
	  YLABEL='\gs\di\un Phot. Rate.'
	  CALL DP_CURVE(NLEV,XLEV,YLEV)
!
	
	
	ELSE IF(OPTION(1:4) .EQ. 'FRAC')THEN
	  DO I=1,ND
	    YVEC(I)=DI(I)/ED(I)
	  END DO
	  CALL DP_CURVE_LAB(ND,XVEC,YVEC,'DI/ED')
!
	ELSE IF(OPTION(1:5) .EQ. 'RHTDR')THEN
!
	   OPEN(UNIT=10,FILE='HTDR_RATES',STATUS='OLD',ACTION='READ')
	   STRING=' '
	   DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	     READ(10,'(A)')STRING
	   END DO
	   READ(STRING,*)NUM_HTDR_RATES
	   DO I=1,NUM_HTDR_RATES
	     READ(10,'(A)')HTDR_SPECIES
	     READ(10,*)NREC
	     READ(10,*)K,T1,T2
	     IF(HTDR_SPECIES .EQ. SPECIES)THEN
	        YVEC(1:ND)=T1*TEMP(1:ND)**(-T2)
	        CALL DP_CURVE_LAB(ND,XVEC,YVEC,'File recom. rate')
	     END IF
	     READ(10,*)NHT
	     READ(10,*)K,T1,T2,T3,T4
	     IF(HTDR_SPECIES .EQ. SPECIES)THEN
	       T1=T1/1.0E+06_LDP; T3=T3/1.0E+04_LDP; T4=T4/1.0E+04_LDP
	       DO K=1,ND
	         YVEC(K)=T1*(TEMP(K)**(-1.5_LDP))*EXP(-T3/TEMP(K))*(1.0_LDP+T2*EXP(-T4/TEMP(K)))
	       END DO
	       CALL DP_CURVE_LAB(ND,XVEC,YVEC,'HTDR recom. rate')
	     END IF
	     IF(HTDR_SPECIES .EQ. SPECIES)EXIT
	   END DO
	   CLOSE(UNIT=10)
!
	ELSE IF(OPTION(1:5) .EQ. 'NRR')THEN
	  DO I=1,ND
	    YVEC(I)=RECOM_SUM(I)/ED(I)/DI(I)
	  END DO
	  CALL DP_CURVE_LAB(ND,XVEC,YVEC,'Recom. Rate')
!
	ELSE IF(OPTION(1:5) .EQ. 'RATES')THEN
!
	  NORM=.TRUE.;  ABS_VALUE=.FALSE.
	  CALL GEN_IN(NORM,'Normalize rates')
	  CALL GEN_IN(ABS_VALUE,'Plot absolute values')
	  WRITE(6,*)' '
	  WRITE(6,*)'Calling DP_CURVE_LAB'
	  WRITE(6,*)'                 +ve means recombination'
	  WRITE(6,*)'                 -ve means ionizing'
	  WRITE(6,*)' '
!
	  YLABEL='Rates'; I=1
	  CALL PLT_RATE(XVEC,PHOT_SUM,TOT_SUM,'Photoionization rate','PR',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,RECOM_SUM,TOT_SUM,'Recombination rate','RR',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,NT_IR_1E,TOT_SUM,'Non-thermal ionization rat','NTI(1e)',NORM,ABS_VALUE,I,ND)
	  TA(1:ND)=NT_IR_2E(:,1); CALL PLT_RATE(XVEC,TA,TOT_SUM,'Non-thermal ionization rat (2e)','NTI(2e)',NORM,ABS_VALUE,I,ND)
	  TA(1:ND)=NT_IR_2E(:,2); CALL PLT_RATE(XVEC,TA,TOT_SUM,'Non-thermal ionization rat (2e)','NTI(2e)',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,ADVEC_RR,TOT_SUM,'Advection rate','AR',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,CHG_IR,TOT_SUM,'Charge exch. ionization rate','CEI',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,CHG_RR,TOT_SUM,'Charge exch. recombination rate','CER',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,COL_IR,TOT_SUM,'Collisional ionization rate','CI',NORM,ABS_VALUE,I,ND)
	  CALL PLT_RATE(XVEC,COL_RR,TOT_SUM,'Collisional recombination rate','CR',NORM,ABS_VALUE,I,ND)
	  TA(1:ND)=XRAY_IR_2E(:,1); CALL PLT_RATE(XVEC,TA,TOT_SUM,'Net X-ray ionization rate','XI(2e)',NORM,ABS_VALUE,I,ND)
	  TA(1:ND)=XRAY_IR_2E(:,2); CALL PLT_RATE(XVEC,TA,TOT_SUM,'Net X-ray ionization rate','XI(2e)',NORM,ABS_VALUE,I,ND)
!
	ELSE IF(OPTION .EQ. 'P')THEN
	  CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	ELSE IF(OPTION(1:5) .EQ. 'WEIRD')THEN
	  WRITE(6,*)ED
	  WRITE(6,*)DI
	  J=1
	  DO WHILE(J .NE. 0)
	    CALL GEN_IN(J,'Level - zero to cease')
	    IF(J .NE. 0)THEN
	      YVEC(1:ND)=1.0E+14_LDP*PHOT(1:ND,J)/ED/DI
	      CALL DP_CURVE(ND,XVEC,YVEC)
	      YVEC(1:ND)=1.0E+14_LDP*RECOM(1:ND,J)/ED/DI
	      CALL DP_CURVE(ND,XVEC,YVEC)
	    END IF
	  END DO
	  CALL GRAMON_PGPLOT(XLABEL,'1.0D+14\ga',' ',' ')
	  GOTO 2000
	ELSE
	  WRITE(6,*)'Unrecognized X-axis option'
	  GOTO 2000
	END IF
	GOTO 2000
!
	END
!
	SUBROUTINE PLT_RATE(XVEC,PHOT_SUM,TOT_SUM,DESC,CURVE_ID,NORM,ABS_VALUE,IP,ND)
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IP
!
	REAL(KIND=LDP) XVEC(ND)
	REAL(KIND=LDP) PHOT_SUM(ND)
	REAL(KIND=LDP) TOT_SUM(ND)
	REAL(KIND=LDP) YVEC(ND)
	CHARACTER(LEN=*) DESC
	CHARACTER(LEN=*) CURVE_ID
!
	INTEGER I
	LOGICAL NORM
	LOGICAL ABS_VALUE
!
	IF(MAXVAL(ABS(PHOT_SUM)) .EQ. 0.0_LDP)RETURN
!
	IF(NORM)THEN
	  YVEC=PHOT_SUM/TOT_SUM
	  CALL DP_CURVE_LAB(ND,XVEC,YVEC,CURVE_ID)
	ELSE IF(ABS_VALUE)THEN
	  YVEC=-40
	  DO I=1,ND
	    IF(PHOT_SUM(I) .GT. 0.0_LDP)YVEC(I)=LOG10(PHOT_SUM(I))
	  END DO
	  CALL DP_CURVE_LAB(ND,XVEC,YVEC,CURVE_ID)
	ELSE
	  YVEC=PHOT_SUM
	  CALL DP_CURVE_LAB(ND,XVEC,YVEC,CURVE_ID)
	END IF
	WRITE(6,'(A,I2,2A)')PG_PEN(IP+1)//'   Curve ',IP,TRIM(DESC),DEF_PEN
	IP=IP+1
!
	RETURN
	END
