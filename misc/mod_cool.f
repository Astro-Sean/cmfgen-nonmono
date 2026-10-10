!
! Program to put the GENCOOL file into a more user friendly format.
! Two files are created:
!
!     GENCOOL_SUM:   Same as GENCOOL but we have summed up over all bound-free rates.
!     GENSCOOL_SORT: Only the top rates are printed: Sorted using depths 1, 11, 21 etc.
!
	PROGRAM MOD_COOL
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 24-Sep-2023  Added SHOCK_POWER term (5-Sep-2023).
! Altered: 03-Aug-2019  Now create a GENCOOL_TYPE file whih summarizes cooling, at each depth,
!                             as a function of type (OSIRIS 21-Jul-2019).
! Altered: 18-Oct-2017  Fixed estimate of cooling time.
! Altered: 01-Mar-2016  Option to omit advection terms when not included in the model [25-Feb-2016].
! Altered: 17-Feb-2015  Improved estimate of cooling time by including ATOM_DENSITY.
! Altered:              Estimate of cooling time outout to GENCOOL_SORT
! Altered: 12-Mar-2014: Added read/ouput of non-thermal cooling.
! Altered: 17-Nov-2009: Now read in charge exchange cooling.
!                         Slight format change.
! Altered: 29-Jan-2009: ND is now read in from MODEL (if it exists).
! Altered: 08-Feb-2008: Extra terms (such as V term) sheck and output.
!
	INTEGER, PARAMETER :: MAX_RECS=1000
!
	CHARACTER*132 TMP_STR
	CHARACTER*132 STRING
	CHARACTER*132 STR_VEC(MAX_RECS)
	REAL(KIND=LDP) VALS(MAX_RECS,10)
	REAL(KIND=LDP) TA(MAX_RECS)
	INTEGER INDX(MAX_RECS)
!
	REAL(KIND=LDP), ALLOCATABLE :: BOUND(:)
	REAL(KIND=LDP), ALLOCATABLE :: SUM(:)
	REAL(KIND=LDP), ALLOCATABLE :: TOTAL_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: ATOM_DENSITY(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: NET_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: COOLING_TIME(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: RAD_DECAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: SHOCK_POWER(:)
	REAL(KIND=LDP), ALLOCATABLE :: COL_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: BF_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: FF_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: XKS_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: NT_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHG_RATE(:)
	REAL(KIND=LDP), ALLOCATABLE :: ART_HT(:)
	REAL(KIND=LDP), ALLOCATABLE :: AC_RV(:)
	REAL(KIND=LDP), ALLOCATABLE :: AC_RdT(:)
!
	INTEGER ND
	INTEGER NV
	INTEGER I,J,K,ID
	INTEGER N_INIT_RECS
	INTEGER NRECS
	INTEGER IOS
	REAL(KIND=LDP) T1
	LOGICAL FILE_OPEN
	LOGICAL ONLY_INCLUDED_TERMS
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	ONLY_INCLUDED_TERMS=.TRUE.
	CALL GEN_IN(ONLY_INCLUDED_TERMS,'Only output terms that are included?')
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of depth points')
	END IF
!
	ALLOCATE (SUM(ND)); SUM=0.0_LDP
	ALLOCATE (BOUND(ND));  BOUND=0.0_LDP
	ALLOCATE (TOTAL_RATE(ND)); TOTAL_RATE=0.0_LDP
	ALLOCATE (COOLING_TIME(ND)); COOLING_TIME=0.0_LDP
	ALLOCATE (ATOM_DENSITY(ND)); ATOM_DENSITY=0.0_LDP
	ALLOCATE (ED(ND)); ED=0.0_LDP
	ALLOCATE (T(ND)); T=0.0_LDP
	ALLOCATE (R(ND)); R=0.0_LDP
	ALLOCATE (V(ND)); V=0.0_LDP
	ALLOCATE (NET_RATE(ND)); NET_RATE=0.0_LDP
!
	ALLOCATE (RAD_DECAY(ND));    RAD_DECAY=0.0_LDP
	ALLOCATE (SHOCK_POWER(ND)); SHOCK_POWER=0.0_LDP
	ALLOCATE (COL_RATE(ND));       COL_RATE=0.0_LDP
	ALLOCATE (BF_RATE(ND));        BF_RATE=0.0_LDP
	ALLOCATE (FF_RATE(ND));        FF_RATE=0.0_LDP
	ALLOCATE (NT_RATE(ND));        NT_RATE=0.0_LDP
	ALLOCATE (XKS_RATE(ND));      XKS_RATE=0.0_LDP
	ALLOCATE (ART_HT(ND));          ART_HT=0.0_LDP
	ALLOCATE (CHG_RATE(ND));      CHG_RATE=0.0_LDP
	ALLOCATE (AC_RV(ND));            AC_RV=0.0_LDP
	ALLOCATE (AC_RdT(ND));          AC_RDT=0.0_LDP
!
!
	OPEN(UNIT=20,FILE='GENCOOL',STATUS='OLD',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='WRITE')
!
	DO I=1,1+(ND-1)/10
	   ID=1+(I-1)*10
!
	   READ(20,'(A)')STRING
	   WRITE(21,'(A)')TRIM(STRING)
	   DO WHILE(1 .EQ. 1)
	     READ(20,'(A)')TMP_STR
	     IF(INDEX(TMP_STR,'Radius') .NE. 0)THEN
	       WRITE(21,'(A,T13,10I12)')'Depth',(K,K=ID,MIN(ID+9,ND))
	       READ(20,'(A)')STRING
	       READ(STRING,*)(R(K),K=ID,MIN(ID+9,ND))
	       WRITE(21,'(A)')TMP_STR(4:9)//'     '//TRIM(STRING)
	     ELSE IF(INDEX(TMP_STR,'Velocity') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       READ(STRING,*)(V(K),K=ID,MIN(ID+9,ND))
	       WRITE(21,'(A)')TMP_STR(4:11)//'   '//TRIM(STRING)
	     ELSE IF(INDEX(TMP_STR,'Temperature') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       READ(STRING,*)(T(K),K=ID,MIN(ID+9,ND))
	       WRITE(21,'(A)')TMP_STR(4:11)//'.  '//TRIM(STRING)
	     ELSE IF(INDEX(TMP_STR,'Electron') .NE. 0)THEN
	       READ(20,'(A)')STRING
	       READ(STRING,*)(ED(K),K=ID,MIN(ID+9,ND))
	       WRITE(21,'(A)')TMP_STR(4:11)//'   '//TRIM(STRING)
	       EXIT
	     END IF
	   END DO
!
	   T1=0.0_LDP
	   DO WHILE(1 .EQ. 1)
	      READ(20,'(A)')STRING
	      IF( INDEX(STRING,'Coll') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,COL_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'COL ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Free-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,FF_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'FF ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Non-thermal') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,NT_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'NT ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'K-shell') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,XKS_RATE,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')TMP_STR(1:K)//'XKS ',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'V term') .NE. 0)THEN
	        IF(ONLY_INCLUDED_TERMS .AND. INDEX(STRING,'Not Incl') .NE. 0)THEN
	          READ(20,'(A)')STRING
	        ELSE
	          TMP_STR=STRING
	          TMP_STR=ADJUSTL(TMP_STR)
	          K=INDEX(TMP_STR,' ')
	          READ(20,'(A)')STRING
	          CALL SUM_RATES(TOTAL_RATE,AC_RV,STRING,ID,ND)
	          WRITE(21,'(A)')' '
	          WRITE(21,'(A,T12,A)')'AC.R(V).',TRIM(STRING)
	        END IF
	      ELSE IF( INDEX(STRING,'dTdR term') .NE. 0)THEN
	        IF(ONLY_INCLUDED_TERMS .AND. INDEX(STRING,'Not Incl') .NE. 0)THEN
	          READ(20,'(A)')STRING
	        ELSE
	          TMP_STR=STRING
	          TMP_STR=ADJUSTL(TMP_STR)
	          K=INDEX(TMP_STR,' ')
	          READ(20,'(A)')STRING
	          CALL SUM_RATES(TOTAL_RATE,AC_RdT,STRING,ID,ND)
	          WRITE(21,'(A,T12,A)')'AC.R(dT).',TRIM(STRING)
	        END IF
	      ELSE IF( INDEX(STRING,'decay') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,RAD_DECAY,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')'|R. decay|',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Shock Power Term') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        CALL SUM_RATES(TOTAL_RATE,SHOCK_POWER,STRING,ID,ND)
	        WRITE(21,'(A,T12,A)')'|Sh. Pow.|',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Artificial') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'|Art. HT|',TRIM(STRING)
	      ELSE IF( INDEX(STRING,'Rate') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'Net C.R.',TRIM(STRING)
	        READ(STRING,*)(NET_RATE(K),K=ID,MIN(ID+9,ND))
	     ELSE IF( INDEX(STRING,'Net') .NE. 0)THEN
	        TMP_STR=STRING
	        TMP_STR=ADJUSTL(TMP_STR)
	        K=INDEX(TMP_STR,' ')
	        READ(20,'(A)')STRING
	        WRITE(21,'(A,T12,A)')'% C.R.',TRIM(STRING)
	        IF(I .NE. 1+(ND-1)/10)THEN
	          READ(20,'(A)')STRING
	          WRITE(21,'(A)')STRING
	        END IF
	        EXIT
	      ELSE IF( INDEX(STRING,'Bound-Free') .NE. 0)THEN
	        TMP_STR=STRING
	        SUM=0.0_LDP
	        DO WHILE(1 .EQ. 1)
	          READ(20,'(A)')STRING
	          IF(STRING .EQ. ' ')THEN
	            TMP_STR=ADJUSTL(TMP_STR)
	            K=INDEX(TMP_STR,' ')
	            WRITE(21,'(A)')STRING
                    WRITE(21,'(A,T13,10ES12.4)')TMP_STR(1:K)//'BF ',(SUM(K),K=ID,MIN(ID+9,ND))
	            EXIT
	          ELSE
	            READ(STRING,*)(BOUND(K),K=ID,MIN(ID+9,ND))
	            CALL SUM_RATES(TOTAL_RATE,BF_RATE,STRING,ID,ND)
	            SUM=SUM+BOUND
	          END IF
	        END DO
	      ELSE IF( INDEX(STRING,'Charge exchange cooling rate') .NE. 0)THEN
	        TMP_STR=STRING
	        SUM=0.0_LDP
	        DO WHILE(1 .EQ. 1)
	          READ(20,'(A)')STRING
	          IF(STRING .EQ. ' ')THEN
	            TMP_STR=ADJUSTL(TMP_STR)
	            K=INDEX(TMP_STR,' ')
                    WRITE(21,'(A,T13,10ES12.4)')'Charge',(SUM(K),K=ID,MIN(ID+9,ND))
	            EXIT
	          ELSE
	            READ(STRING,*)(BOUND(K),K=ID,MIN(ID+9,ND))
	            CALL SUM_RATES(TOTAL_RATE,CHG_RATE,STRING,ID,ND)
	            SUM=SUM+BOUND
	          END IF
	        END DO
	      END IF
	   END DO
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' Cooling data has been written to GENCOOL_SUM'
	WRITE(6,'(A)')' Will now sort data to display most important terms'
	WRITE(6,'(A)')' 12 records is a reasonable number to output '
	WRITE(6,'(A)')' Cooling time is approximate and not valid at high densities'
!
! We have summed all rates without regard to sign. If in equilibrium, the cooling
! rate must be half this rate. This will not work at high densities where the is
! cancellation in heating/cooling rates for individual species/processes.
!
! For simplicty, we assume the energy per species is 1.5kT, and that the number
! of electrons is the same as the number of ions when the atom density is unavailable.
!
	TOTAL_RATE=TOTAL_RATE*0.5_LDP
	I=7
	CALL RD_SING_VEC_RVTJ(ATOM_DENSITY,ND,'Atom Density','RVTJ',I,IOS)
	IF(IOS .NE. 0)THEN
	  COOLING_TIME=3.0_LDP*1.3806E-12_LDP*T*ED/TOTAL_RATE
	ELSE
	  COOLING_TIME=1.5_LDP*1.3806E-12_LDP*T*(ED+ATOM_DENSITY)/TOTAL_RATE
	END IF
!
	OPEN(UNIT=20,FILE='GENCOOL_SUM',STATUS='UNKNOWN',ACTION='READ')
	OPEN(UNIT=21,FILE='GENCOOL_SORT',STATUS='UNKNOWN',ACTION='WRITE')
	NRECS=12
	CALL GEN_IN(NRECS,'Maximum number of records per depth to be output to sorted file')
!
	N_INIT_RECS=6
	VALS=0.0_LDP
	DO I=1,1+(ND-1)/10
	  ID=1+(I-1)*10
	  DO K=1,N_INIT_RECS
	    READ(20,'(A)')STR_VEC(K)
	  END DO
!
	  K=N_INIT_RECS
	  DO WHILE(1 .EQ. 1)
	    K=K+1
100	    READ(20,'(A)')STR_VEC(K)
	    IF(STR_VEC(K) .EQ. ' ')GOTO 100
	    READ(STR_VEC(K)(12:),*)(VALS(K,J),J=1,MIN(10,ND-(I-1)*10))
	    IF( INDEX(STR_VEC(K),'%') .NE. 0)THEN
	      IF(I .NE. 1+(ND-1)/10 )READ(20,'(A)')STRING			!Getting record with ^L
	      EXIT
	    END IF
	  END DO
!
	  NV=K-N_INIT_RECS-2
	  TA(:)=ABS(VALS(:,1))
	  CALL INDEXX(NV,TA(N_INIT_RECS+1),INDX,L_FALSE)
	  NV=NV+N_INIT_RECS+2
!
	  DO K=1,N_INIT_RECS
	    WRITE(21,'(A)')TRIM(STR_VEC(K))
	  END DO
	  WRITE(21,'(A)')' '
	  WRITE(21,'(A)')TRIM(STR_VEC(NV))
	  WRITE(21,'(A)')TRIM(STR_VEC(NV-1))
	  WRITE(21,'(A)')' '
	  WRITE(21,'(A,10ES12.4)')'Cool time(s)',(COOLING_TIME(K),K=ID,MIN(ID+9,ND))
	
	  WRITE(21,'(A)')' '
	  DO K=1,MIN(NRECS,NV-N_INIT_RECS-2)
	    WRITE(21,'(A)')TRIM(STR_VEC(INDX(K)+N_INIT_RECS))
	  END DO
!
	END DO
	CLOSE(UNIT=20)
	CLOSE(UNIT=21)
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Sorted cooling data has been written to GENCOOL_SORT'
	WRITE(6,'(A)')' '
!
	OPEN(UNIT=22,FILE='GENCOOL_TYPE',STATUS='UNKNOWN',ACTION='WRITE')
	TOTAL_RATE=COL_RATE+BF_RATE+FF_RATE+NT_RATE+XKS_RATE+AC_RV
	TOTAL_RATE=TOTAL_RATE+AC_RdT+CHG_RATE-RAD_DECAY-ART_HT-SHOCK_POWER
	DO ID=1,ND,10
	  WRITE(22,'(A,T12,10I12)')' Depth',(K,K=ID,MIN(ID+9,ND))
	  WRITE(22,'(A,T12,10ES12.4)')' R(10^10cm)',(R(K),K=ID,MIN(ID+9,ND))
	  IF(V(1) .GT. 0.0_LDP)THEN
	    WRITE(22,'(A,T12,10ES12.4)')' V(km/s)',(V(K),K=ID,MIN(ID+9,ND))
	  END IF
	  WRITE(22,'(A,T12,10ES12.4)')' T(10^4K)',(T(K),K=ID,MIN(ID+9,ND))
	  WRITE(22,'(A,T12,10ES12.4)')' ED', (ED(K),K=ID,MIN(ID+9,ND))
!
	  WRITE(22,'(A)')' '
	  WRITE(22,'(A,T12,10ES12.4)')' |R. decay|',(RAD_DECAY(K),K=ID,MIN(ID+9,ND))
	  WRITE(22,'(A,T12,10ES12.4)')' |Sh. Pow.|',(SHOCK_POWER(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' Art Ht.',(ART_HT(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' Col. cool',(COL_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' BF cool',(BF_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' FF cool',(FF_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' NT cool',(NT_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' K shell',(XKS_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' AC.R(V)',(AC_RV(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' AC.R(dT)',(AC_RdT(K),K=ID,MIN(ID+9,ND))
	  WRITE(22,'(A)')' '
          WRITE(22,'(A,T12,10ES12.4)')' Net col',(TOTAL_RATE(K),K=ID,MIN(ID+9,ND))
          WRITE(22,'(A,T12,10ES12.4)')' NET(file)',(NET_RATE(K),K=ID,MIN(ID+9,ND))
	  WRITE(22,'(A)')' '
	  WRITE(22,'(132A)')' ',('*',I=1,131)
	  WRITE(22,'(A)')' '
	END DO
!
	STOP
	END
!
!
!
	SUBROUTINE SUM_RATES(TOTAL_RATE,SPECIFIC_RATE,STRING,ID,ND)
	USE SET_KIND_MODULE
	INTEGER ID, ND
	REAL(KIND=LDP) TOTAL_RATE(ND)
	REAL(KIND=LDP) SPECIFIC_RATE(ND)
	REAL(KIND=LDP) TEMP_VEC(ND)
	CHARACTER(LEN=*) STRING
!
	INTEGER I
!
	READ(STRING,*)(TEMP_VEC(I),I=ID,MIN(ID+9,ND))
	DO I=ID,MIN(ID+9,ND)
	  TOTAL_RATE(I)=TOTAL_RATE(I)+ABS(TEMP_VEC(I))
	  SPECIFIC_RATE(I)=SPECIFIC_RATE(I)+TEMP_VEC(I)
	END DO
!
	RETURN
	END

