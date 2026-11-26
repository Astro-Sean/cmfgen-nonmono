	PROGRAM PLT_COOL_SORT
	USE SET_KIND_MODULE
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 17-Nov-2021 -- Copied from OSIRIS
! Altered: 16-Nov-2021 -- Cleaned and made option driven.
! Altered: 15-Jun-2021 --  Small bug fixis (transferred from OSIRIS 21-Jul-20201).
! Altered: 29-Mar-2020 -- Designed for nebular SN only (at present).
!
	INTEGER, PARAMETER :: ND_MAX=400
	INTEGER, PARAMETER :: NREC_MAX=200
	INTEGER, PARAMETER :: N_TIT=5
	INTEGER, PARAMETER :: LU_RD=7
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL(KIND=LDP) R(ND_MAX)
	REAL(KIND=LDP) V(ND_MAX)
	REAL(KIND=LDP) T(ND_MAX)
	REAL(KIND=LDP) ED(ND_MAX)
	REAL(KIND=LDP) PER_CR(ND_MAX)
	REAL(KIND=LDP) NET_CR(ND_MAX)
	REAL(KIND=LDP) COOL_TIME(ND_MAX)
	REAL(KIND=LDP) XV(ND_MAX)
	REAL(KIND=LDP) YV(ND_MAX)
	REAL(KIND=LDP) TOTAL(ND_MAX)
	REAL(KIND=LDP) TMP_VEC(ND_MAX)
	REAL(KIND=LDP) LUM(NREC_MAX)
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) LSTAR,LSUN
!
	REAL(KIND=LDP), ALLOCATABLE :: RATE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: VRATE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: dLUM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_RATE(:,:)
!
	CHARACTER(LEN=20)  LABEL(NREC_MAX)
	CHARACTER(LEN=80)  FILENAME
	CHARACTER(LEN=200) STRING
	CHARACTER(LEN=20)  TMP_STR
	CHARACTER(LEN=100) TITLE(N_TIT)
	CHARACTER(LEN=6)   COL
	CHARACTER(LEN=30)  XAXIS
	CHARACTER(LEN=30)  YAXIS
	CHARACTER(LEN=20) PLT_OPT
!
	INTEGER NREC
	INTEGER IOS
	INTEGER I,J,L,K,IBEG
	INTEGER IST,IEND
	INTEGER ND
	INTEGER LUM_STAR
	INTEGER KK(1)
	INTEGER INDX(2)
	INTEGER I_DECAY
	INTEGER POINT(NREC_MAX)
	INTEGER NPLT
	INTEGER IVEC(5)
!
	REAL(KIND=LDP) LUM_SUN
	INTEGER GET_INDX_DP
	EXTERNAL GET_INDX_DP,LUM_SUN
	LOGICAL USE_DECAY
	LOGICAL READ_RATE
	LOGICAL FIRST_TIME
	LOGICAL OMIT
	CHARACTER(LEN=3) XAX_VAR
!
	FIRST_TIME=.TRUE.
	LSUN=LUM_SUN()
!
	WRITE(6,'(A)')RED_PEN
	WRITE(6,'(A)')' This routine is still in devlopment and may need'
	WRITE(6,'(A)')' further changes.'
	WRITE(6,'(A)')' '
!
	WRITE(6,'(A)')DEF_PEN
        ND=0
	OPEN(UNIT=LU_RD,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LU_RD,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	      ELSE IF(INDEX(STRING,'[LSTAR]') .NE. 0)THEN
	        READ(STRING,*)LSTAR
	        EXIT
	      END IF
	    END DO
	    CLOSE(LU_RD)
	  END IF
	IF(ND .EQ. 0)CALL GEN_IN(ND,'Number of depth points in model')
        CALL RD_SING_VEC_RVTJ(V,ND,'Velocity','RVTJ',LU_RD,IOS)
!
! Initialize large arrays
!
	ALLOCATE (RATE(NREC_MAX,ND))
	ALLOCATE (VRATE(NREC_MAX,ND))
	ALLOCATE (dLUM(NREC_MAX,ND))
	RATE=0.0_LDP; VRATE=0.0_LDP; dLUM=0.0_LDP
!
 	FILENAME='GENCOOL_SORT'
	CALL GEN_IN(FILENAME,'GENCOOL_SORT file with data to be plotted')
	IF(FILENAME .EQ. ' ')STOP
!
	OPEN(UNIT=LU_RD,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	  NREC=0
	  RATE=0.0_LDP
	  DO IBEG=1,ND,10
	    DO WHILE(1 .EQ. 1)
	      READ(LU_RD,'(A)')STRING
	      K=INDEX(STRING,'  ')
	      IF(INDEX(STRING,'Radius') .NE. 0)THEN
	        READ(STRING(K:),*)(R(I),I=IBEG,MIN(IBEG+9,ND))
	      ELSE IF(INDEX(STRING,'Tempr') .NE. 0)THEN
	        READ(STRING(K:),*)(T(I),I=IBEG,MIN(IBEG+9,ND))
	      ELSE IF(INDEX(STRING,'Electr') .NE. 0)THEN
	        READ(STRING(K:),*)(ED(I),I=IBEG,MIN(IBEG+9,ND))
	      ELSE IF(INDEX(STRING,'% C.R.') .NE. 0)THEN
	        READ(STRING(K:),*)(PER_CR(I),I=IBEG,MIN(IBEG+9,ND))
	      ELSE IF(INDEX(STRING,'Net C.R.') .NE. 0)THEN
	        READ(STRING(K:),*)(NET_CR(I),I=IBEG,MIN(IBEG+9,ND))
	      ELSE IF(INDEX(STRING,'Cool time') .NE. 0)THEN
	        READ(STRING(K:),*)(COOL_TIME(I),I=IBEG,MIN(IBEG+9,ND))
	        EXIT
	      END IF
	    END DO
!
	    DO WHILE(1 .EQ. 1)
	      READ(LU_RD,'(A)',END=2000)STRING
	      READ_RATE=.FALSE.
	      IF(STRING .EQ. ' ')THEN
	      ELSE IF(INDEX(STRING,'Depth') .NE. 0)THEN
	        EXIT
	      ELSE
	        K=INDEX(STRING,'  ')
	        DO L=1,NREC
	          IF(STRING(1:K) .EQ. LABEL(L))THEN
	            READ(STRING(K:),*)(RATE(L,I),I=IBEG,MIN(IBEG+9,ND))
	            READ_RATE=.TRUE.
	            EXIT
	          END IF
	        END DO
	        IF(.NOT. READ_RATE)THEN
	          NREC=NREC+1
	          L=NREC
	          READ(STRING(K:),*)(RATE(L,I),I=IBEG,MIN(IBEG+9,ND))
	          LABEL(NREC)=STRING(1:K)
	        END IF
	      END IF
	    END DO
	  END DO
2000	  CONTINUE
	CLOSE(LU_RD)
	WRITE(6,*)'Read in data'
	WRITE(6,*)'ND=',ND
	WRITE(6,*)'NREC=',NREC
!
	I=0; J=0;
	DO L=1,NREC
	  IF(LABEL(L) .EQ. 'AC.R(V).')I=L
	  IF(LABEL(L) .EQ. 'AC.R(dT).')J=L
	END DO
	IF(I .NE. 0 .AND. J .NE. 0)THEN
	  RATE(I,:)=RATE(I,:)+RATE(J,:)
	  LABEL(I)='Ad. Cool.'
	  RATE(J,:)=0.0D0
	  WRITE(6,*)'Adiabatic terms combined'
	END IF
!
! Fix label names.
!
	DO L=1,NREC
	  IF(INDEX(LABEL(L),'2') .NE. 0)THEN
	    K=INDEX(LABEL(L),'2')
	    LABEL(L)(K:)='II'//LABEL(L)(K+1:)
	  ELSE IF(INDEX(LABEL(L),'SIX') .NE. 0)THEN
	    K=INDEX(LABEL(L),'SIX')
	    LABEL(L)(K:)='VI'//LABEL(L)(K+3:)
	  ELSE IF(INDEX(LABEL(L),'SEV') .NE. 0)THEN
	    K=INDEX(LABEL(L),'SEV')
	    LABEL(L)(K:K+2)='VII'
	  ELSE IF(INDEX(LABEL(L),'k') .NE. 0)THEN
	    K=INDEX(LABEL(L),'k')
	    LABEL(L)(K:)='i'//LABEL(L)(K+1:)
	  END IF
	END DO
!
! Determine the index for radiactive decay.
!
	I_DECAY=0
	DO L=1,NREC
	  IF( INDEX(LABEL(L),'decay') .NE. 0)THEN
	    I_DECAY=L
	    EXIT
	  END IF
	END DO
	IF(I_DECAY .NE. 0)WRITE(6,*)'I_DECAY=',I_DECAY
!
! Set default X axis
!
	XV(1:ND)=V(1:ND)
	XAXIS='V(km/s)'
	IF(V(1) .GT. 10000.0_LDP)THEN
	  XAXIS='V(Mm/s)'
	  XV(1:ND)=1.0E-03_LDP*V(1:ND)
	END IF
!
! Determine integrated luminosities
!
	dLUM=0.0_LDP
	T1=4.0_LDP*3.1459_LDP*1.0E+30_LDP
	T2=1.0_LDP
	IF(I_DECAY .EQ. 0)T2=LSUN
	DO L=1,NREC
	  YV(1:ND)=T1*RATE(L,1:ND)*R(1:ND)*R(1:ND)
	  CALL LUM_FROM_ETA(YV,R,ND)
	  dLUM(L,1:ND-1)=YV(1:ND-1)
	  LUM(L)=SUM(dLUM(L,1:ND-1))
	  IF(MOD(L,3) .EQ. 0)THEN
	    WRITE(6,'(A,3X,I5,ES16.4,6X)')LABEL(L)(1:12),L,LUM(L)/LSUN
	  ELSE
	    WRITE(6,'(A,3X,I5,ES16.4,6X)',ADVANCE='NO')LABEL(L)(1:12),L,LUM(L)/LSUN
	  END IF
	END DO
!
! Enter main plotting/examination loop.
!
	DO WHILE(1 .EQ. 1)
          WRITE(6,*)' '
          PLT_OPT='P'
          CALL GEN_IN(PLT_OPT,'Plot option')
	  CALL SET_CASE_UP(PLT_OPT,IZERO,IZERO)
!
	  IF(PLT_OPT .EQ. 'NRM' .OR. FIRST_TIME)THEN
	    USE_DECAY=.FALSE.
	    CALL GEN_IN(USE_DECAY,'Use radioactive decay for nomalization?')
	    IF(USE_DECAY)THEN
	      TOTAL(1:ND)=RATE(I_DECAY,1:ND)
	      YAXIS='Cooling rate/E[decay]'
	    ELSE
	      TOTAL=0.0_LDP
	      DO I=1,ND
	        DO L=1,NREC
	          TOTAL(I)=TOTAL(I)+ABS(RATE(L,I))
	        END DO
	      END DO
	      TOTAL(1:ND)=0.5_LDP*TOTAL(1:ND)                !Average of heating and cooling rate
	      YAXIS='Cooling rate/E[heat]'
	    END IF
	    WRITE(6,*)'Computed total cooling/heating rate'
!
! Normalize the data so can choose from all depths easily.
!
	    VRATE=RATE
	    DO L=1,NREC
	      VRATE(L,1:ND)=VRATE(L,1:ND)/TOTAL(1:ND)
	    END DO
	    IF(FIRST_TIME)THEN
	      WRITE(6,'(/,A)')' Need to do option again as needed to set some arrays'
	      FIRST_TIME=.FALSE.
	    END IF
!
	  ELSE IF(PLT_OPT .EQ. 'XVEL')THEN
	    XV(1:ND)=V(1:ND)
	    XAXIS='V(km/s)'
	    IF(V(1) .GT. 10000.0_LDP)THEN
	      XAXIS='V(Mm/s)'
	      XV(1:ND)=1.0E-03_LDP*V(1:ND)
	    END IF
	  ELSE IF(PLT_OPT .EQ. 'XR')THEN
	    XV(1:ND)=R(1:ND)/R(ND)
	    XAXIS='R/R(ND)'
	  ELSE IF(PLT_OPT .EQ. 'XN')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    XAXIS='Depth index'
!
	  ELSE IF(PLT_OPT .EQ. 'COOL')THEN
	    NPLT=10
	    CALL GEN_IN(NPLT,'Number of cooling curves to be illustrated (0 to read from file')
	    CALL GEN_IN(OMIT,'Omit depth from plotting?')
	    WRITE(6,'(A,A)')'X vector is ',TRIM(XAXIS)
	    IF(XV(1) .LT. XV(ND))THEN
	      T1=XV(1);  CALL GEN_IN(T1,'First depth for plotting')
	      T2=XV(ND); CALL GEN_IN(T2,'Last depth for plotting')
	    ELSE
	      T1=XV(ND); CALL GEN_IN(T1,'First depth for plotting')
	      T2=XV(1);  CALL GEN_IN(T2,'Last depth for plotting')
	    END IF
	    IF(XAXIS .EQ. 'Depth index')THEN
	      IST=NINT(T1); IEND=NINT(T2)
	    ELSE
	      IST=GET_INDX_DP(T1,XV,ND)
	      IEND=GET_INDX_DP(T2,XV,ND)
	      IF(IST .GT. IEND)THEN
	        I=IEND; IEND=IST; IST=I
	      END IF
	      IF(ABS(T2-XV(IEND)) .LT. ABS(T2-XV(IEND-1)))IEND=IEND+1
	    END IF
	    WRITE(6,*)'Depth indices are:',IST,IEND
	    IF(ALLOCATED(TMP_RATE))DEALLOCATE(TMP_RATE)
	    WRITE(6,*)'Deallocated TMP_RATE'
	    ALLOCATE(TMP_RATE(NREC,IEND-IST+1))
	    WRITE(6,*)'Allocated TMP_RATE'
	    TMP_RATE=0.0_LDP; TMP_RATE=VRATE(1:NREC,IST:IEND)
	    OPEN(UNIT=20,FILE='POINTER',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
	    WRITE(20,*)NPLT
	    DO L=1,NPLT
	      INDX=MAXLOC(ABS(TMP_RATE))
	      POINT(L)=INDX(1)
	      WRITE(6,'(A,20X,ES12.4)')TRIM(LABEL(POINT(L))),VRATE(INDX(1),INDX(2))
	      TMP_RATE(INDX(1),:)=0.0_LDP
	      WRITE(20,*)TRIM(LABEL(POINT(L)))
	    END DO
	    DO L=1,NPLT
	      IF(POINT(L) .NE. 0)THEN
	        IF(OMIT)THEN
	          YV(1:ND-1)=100.0_LDP*RATE(POINT(L),2:ND)/TOTAL(2:ND)
	          I=ND-1; CALL DP_CURVE_LAB(I,XV,YV,LABEL(POINT(L)))
	        ELSE
	          YV(1:ND)=100.0_LDP*RATE(POINT(L),1:ND)/TOTAL(1:ND)
	          CALL DP_CURVE_LAB(ND,XV,YV,LABEL(POINT(L)))
	        END IF
	      END IF
	    END DO
!
	
	    IVEC(1)=1; IVEC(2)=ND/4; IVEC(3)=ND/2; IVEC(4)=3*ND/4; IVEC(5)=ND
	    WRITE(6,*)IVEC
	    IF(IST .NE. 1)THEN
	      DO K=2,4
	        IF(IVEC(K) .GT. IST)IVEC(K)=IST
	      END DO
	    END IF
	    IF(IEND .NE. ND .AND. IST .NE. IEND)THEN
	      DO K=4,2,-1
	        IF(IEND .GT. IVEC(K))THEN
	          IVEC(K)=IEND
	          EXIT
	        END IF
	      END DO
	    END IF
	    WRITE(6,*)IVEC
!
	    TMP_STR='%Rate'
	    WRITE(6,'(A)')' '
	    WRITE(6,'(2X,A,T25,5(8X,I4))')'Depth index',(IVEC(K),K=1,5)
	    WRITE(6,'(2X,A,T25,5(ES12.4))')XAXIS,(XV(IVEC(K)),K=1,5)
	    WRITE(6,'(2X,A,3X,A,T25,5(7X,A))')'Index','Type',(TRIM(TMP_STR),K=1,5)
	    DO L=1,NPLT
	      IF(POINT(L) .NE. 0)THEN
	        WRITE(6,'(I7,3X,A,T25,5F12.3)')L,TRIM(LABEL(POINT(L))),
	1           (100.0D0*RATE(POINT(L),IVEC(K))/TOTAL(IVEC(K)), K=1,5)
	      END IF
	    END DO
!
	  ELSE IF(PLT_OPT .EQ. 'RDC')THEN
	    OPEN(UNIT=20,FILE='POINTER',STATUS='OLD',ACTION='READ')
	      READ(20,*)NPLT
	      DO L=1,NPLT
	        READ(20,'(A)')TMP_STR
	        DO K=1,NREC
	          IF(LABEL(K) .EQ. TMP_STR)THEN
	            POINT(L)=K
	            WRITE(6,*)TRIM(TMP_STR),POINT(L),K
	            EXIT
	          END IF
	        END DO
	      END DO
	    CLOSE(UNIT=20)
!
	  ELSE IF(PLT_OPT .EQ. 'LUM')THEN
	    T1=LSTAR*LSUN
	    IF(USE_DECAY)T1=LUM(I_DECAY)
	    DO L=1,NPLT
	      K=POINT(L)
	      YV(1:ND)=dLUM(K,1:ND)/T1
	      CALL DP_CURVE_LAB(ND,XV,YV,LABEL(L))
	    END DO
	    IF(USE_DECAY)THEN
	      YAXIS='Cooling rate/E[decay]'
	      WRITE(6,'(A)')' '
	      WRITE(6,'(A,ES14.4,A)')'Total radiactive decay energy is',LUM(1),'ergs'
	      WRITE(6,'(A,ES14.4,A)')'Total radiactive decay energy is',LUM(1),LUM(1)/3.826D+33,'Lsun'
	      WRITE(6,'(A)')' '
	    ELSE
	      YAXIS='Cooling rate (ergs/cm^3/s)/Lstar'
	    END IF
!
	  ELSE IF(PLT_OPT .EQ. 'H' .OR. PLT_OPT .EQ. 'HE' .OR. PLT_OPT .EQ. 'HELP')THEN
	     WRITE(6,'(A)')
	     WRITE(6,'(A)')'XVEL       Set X axis to V'
	     WRITE(6,'(A)')'XR         Set X axis to R/R*'
	     WRITE(6,'(A)')'XN         Set X axis to depth index'
	     WRITE(6,'(A)')'COOL       Plot curves showing cooling percentages'
	     WRITE(6,'(A)')
	     WRITE(6,'(A)')'P          Plot data stored on GRAMON buffer'
	     WRITE(6,'(A)')'EX(IT)     Exit program'
	     WRITE(6,'(A)')
!
	  ELSE IF(PLT_OPT .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XAXIS,YAXIS,' ',' ')
	  ELSE IF(PLT_OPT .EQ. 'EX' .OR. PLT_OPT .EQ. 'EXIT')THEN
	    STOP
	  ELSE IF(PLT_OPT .EQ. 'STOP')THEN
	    STOP
	  ELSE
	    WRITE(6,*)'Plot option not recognized'
	  END IF
	END DO
!
	END
