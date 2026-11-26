!
	PROGRAM PLT_ION_MF
	USE SET_KIND_MODULE
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Created 24-May-2020
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	REAL(KIND=LDP), ALLOCATABLE :: XV(:)
	REAL(KIND=LDP), ALLOCATABLE :: YV(:)
	REAL(KIND=LDP), ALLOCATABLE :: ZV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
!
	CHARACTER(LEN=10), ALLOCATABLE :: ION_ID(:)
	CHARACTER(LEN=10) ION_BEG,ION_END
	CHARACTER(LEN=30) LEVEL
!
	CHARACTER(LEN=20) FILE_DATE
!
	REAL(KIND=LDP), ALLOCATABLE :: ES_KAPPA(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: ION_FLUX_MEAN(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: DPTH_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_VEC(:)
	INTEGER, ALLOCATABLE :: INDX(:)
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) SUM_TO,SUM_FROM
	INTEGER ND
	INTEGER NIONS
!
	INTEGER I,J,K,L,IBEG,IDEPTH
	INTEGER ID,ID_BEG,ID_END
	INTEGER CNT
	INTEGER DPTH_INDX
	LOGICAL FILE_OPEN
	INTEGER IOS
	LOGICAL PERCENTAGE
!
	INTEGER, PARAMETER :: LU_RVTJ=10
	INTEGER, PARAMETER :: LUIN=11
	INTEGER, PARAMETER :: LUOUT=12	
	INTEGER IVEC(2)
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=80) TMP_FMT
	CHARACTER(LEN=120) DEFAULT
	CHARACTER(LEN=200) STRING
        CHARACTER(LEN=30) UC
	CHARACTER(LEN=20) PLT_OPT
	CHARACTER(LEN=20) XLABEL
	CHARACTER(LEN=20) YLABEL
        CHARACTER(LEN=200) TITLE(10)
	EXTERNAL UC
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot/examine ION_MF'
	WRITE(6,'(A)')' '
!
! Inialization.
!
	ION_BEG=' '; ION_END=' '
	YLABEL='M(t)'
!
100	CONTINUE
 	FILENAME='ION_FLUX_MEAN_OPAC'
	CALL GEN_IN(FILENAME,'File with data to be plotted/examined')
	IF(FILENAME .EQ. ' ')GOTO 100
!
! Get number of depth points
!
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	    READ(LUIN,*,IOSTAT=IOS)NIONS     !!Number of ions'
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading: ',TRIM(FILENAME)
	      STOP
	     END IF
	  ELSE
	    WRITE(6,*)'Error opening: ',TRIM(FILENAME)
	    STOP
	  END IF
!
	ALLOCATE (ES_KAPPA(ND),ROSS_MEAN(ND),FLUX_MEAN(ND))
	I=MAX(ND,NIONS)
	ALLOCATE (XV(I),YV(I),ZV(I),INDX(I))
	ALLOCATE (ION_FLUX_MEAN(ND,NIONS))
	ALLOCATE (ION_ID(NIONS))
!
	DO WHILE(INDEX(STRING,'Kappa electron scattering') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	READ(LUIN,*)ES_KAPPA
	READ(LUIN,'(A)')STRING
	READ(LUIN,*)ROSS_MEAN
	READ(LUIN,'(A)')STRING
	READ(LUIN,*)FLUX_MEAN
	DO ID=1,NIONS
	  READ(LUIN,'(A)')ION_ID(ID)
	  READ(LUIN,*)ION_FLUX_MEAN(:,ID)
	END DO
!
	FILENAME='RVTJ'
	ALLOCATE (R(ND),V(ND),ED(ND),SIGMA(ND))
	IF(IOS .EQ. 0)THEN
	   FILENAME='../RVTJ'
          CALL RD_SING_VEC_RVTJ(R,ND,'Radius',FILENAME,LU_RVTJ,IOS)
	END IF
        CALL RD_SING_VEC_RVTJ(V,ND,'Velocity',FILENAME,LU_RVTJ,IOS)
        CALL RD_SING_VEC_RVTJ(SIGMA,ND,'dlnV/dlnr-1',FILENAME,LU_RVTJ,IOS)
        CALL RD_SING_VEC_RVTJ(ED,ND,'Electron',FILENAME,LU_RVTJ,IOS)
!
5000	CONTINUE
!
! Set def axis.
!
	XV(1:ND)=V(1:ND)
	XLABEL='V(km/s)'
!
! Enter main plotting/examination loop.
!
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  PLT_OPT='P'
	  CALL GEN_IN(PLT_OPT,'Plot option')
!
	  IF(UC(PLT_OPT) .EQ. 'XLOGR')THEN
	    T1=R(ND)
	    DO I=1,ND
	      XV(I)=LOG10(R(I)/T1)
	    END DO
	    XLABEL='Log R/R(ND)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XV')THEN
	    XV(1:ND)=V(1:ND)
	    XLABEL='V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XLOGV')THEN
	    DO I=1,ND
	      IF(V(I) .LE. 0.0_LDP)THEN
	        WRITE(6,*)'Error in PLT_ION_MF: XLOGV requires positive V'
	        WRITE(6,*)'Depth, V=',I,V(I)
	        STOP
	      END IF
	    END DO
	    XV(1:ND)=LOG10(V(1:ND))
	    XLABEL='Log V(km/s)'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'XN')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    XLABEL='Depth index'
!
	  ELSE IF(UC(PLT_OPT(1:4)) .EQ. 'FLUX')THEN
	    CALL DP_CURVE(ND,XV,FLUX_MEAN)
!
	  ELSE IF(UC(PLT_OPT(1:4)) .EQ. 'ROSS')THEN
	    CALL DP_CURVE(ND,XV,ROSS_MEAN)
!
! Set the ion to be examined.
!
	  ELSE IF(UC(PLT_OPT(1:3)) .EQ. 'ION')THEN
	    CALL GEN_IN(ION_BEG,'ION that will be examined: e.g. C2, NIII')
	    YLABEL='M(t)'
	    PERCENTAGE=.FALSE.
	    CALL GEN_IN(PERCENTAGE,'Plot as a percentage (default is M(t)')
	    ION_BEG=UC(ION_BEG)
	    DO ID=1,NIONS
	      IF(UC(ION_ID(ID)) .EQ. ION_BEG)THEN
	        YV(1:ND)=ION_FLUX_MEAN(:,ID)
	        EXIT
	      END IF
	    END DO
	    IF(PERCENTAGE)YV(1:ND)=100.0_LDP*YV(1:ND)/FLUX_MEAN(1:ND)
	    CALL DP_CURVE(ND,XV,YV)
!
	  ELSE IF(UC(PLT_OPT(1:4)) .EQ. 'SPEC')THEN
	    CALL GEN_IN(ION_BEG,'Initial ION that will be examined: e.g. C2')
	    CALL GEN_IN(ION_END,'Final ION that will be examined: e.g. CIV')
	    YLABEL='M(t)'
	    PERCENTAGE=.FALSE.
	    CALL GEN_IN(PERCENTAGE,'Plot as a percentage (default is M(t)')
	    IF(PERCENTAGE)YLABEL='% contribution'
	    ION_BEG=UC(ION_BEG); ION_END=UC(ION_END)
	    ID_BEG=0; ID_END=0
	    DO ID=1,NIONS
	      IF(UC(ION_ID(ID)) .EQ. ION_BEG)THEN
	        ID_BEG=ID
	        EXIT
	      END IF
	    END DO
	    IF(ID_BEG .EQ. 0)THEN
	      WRITE(6,*)'Lower ionization stage not recognized: ',TRIM(ION_BEG)
	      GOTO 1000
	    END IF
	    DO ID=ID_BEG,NIONS
	      IF(UC(ION_ID(ID)) .EQ. ION_END)THEN
	        ID_END=ID
	        EXIT
	      END IF
	    END DO
	    IF(ID_END .EQ. 0)THEN
	      WRITE(6,*)'Upper ionization stage not recognized: ',TRIM(ION_END)
	      GOTO 1000
	    END IF
!
	    CNT=0; TITLE=' '
	    CALL GEN_IN(CNT,'Initial pen color (0 is no plot done befor hand)')
	    DO ID=ID_BEG,ID_END
	      YV(1:ND)=ION_FLUX_MEAN(:,ID)
	      IF(PERCENTAGE)YV(1:ND)=100.0_LDP*YV(1:ND)/FLUX_MEAN(1:ND)
              DEFAULT=TRIM(ION_ID(ID))
	      CALL DP_CURVE_LAB(ND,XV,YV,DEFAULT)
	    END DO
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'LIST')THEN
	    WRITE(6,'(8A)')(ION_ID(I),I=1,NIONS)
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'NORM')THEN
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EXD')THEN
	    DPTH_INDX=ND/2
	    CALL GEN_IN(DPTH_INDX,'Depth to be examined')
	    XV(1:NIONS)=ION_FLUX_MEAN(DPTH_INDX,1:NIONS)
	    CALL INDEXX(NIONS,XV,INDX,L_FALSE)
	    WRITE(6,'(A)')' '
	    WRITE(6,*)'     R/R*',R(DPTH_INDX)/R(ND),'R*',R(ND)
	    WRITE(6,*)'  V(km/s)',V(DPTH_INDX)
	    WRITE(6,*)'Ross Mean',ROSS_MEAN(DPTH_INDX)
	    WRITE(6,*)'Flux Mean',FLUX_MEAN(DPTH_INDX)
	    WRITE(6,'(A)')' '
	    DO I=1,12
	     WRITE(6,'(A,T15,F10.3,4X,F8.3)')ION_ID(INDX(I)),XV(INDX(I)),
	1          100.0D0*XV(INDX(I))/FLUX_MEAN(DPTH_INDX)
	    END DO
	    WRITE(6,'(A)')' '
	    T1=SUM(XV(1:NIONS))+1
	    WRITE(6,'(A,T16,F10.3,4X,F8.3)')'SUM(incl es)',T1,T1/FLUX_MEAN(DPTH_INDX)
!
	  ELSE IF(UC(PLT_OPT(1:1)) .EQ. 'H')THEN
!
	    WRITE(6,*)' '
	    WRITE(6,*)' XN:      Set X-axis to depth index'
	    WRITE(6,*)' XLOGR:   Set X-axis to Log R/R(ND)'
	    WRITE(6,*)' XLOGV:   Set X-axis to Log V'
	    WRITE(6,*)' XV:      Set X-axis to V'
	    WRITE(6,*)' '
	    WRITE(6,*)' SPEC:    Set the species and reads in SPCIES//PRRR file'
	    WRITE(6,*)' LEVEL:   Sets the level to be studied'
	    WRITE(6,*)' EXD:     Examine rates at a single depth'
	    WRITE(6,*)' NORM:    Plot rates as a function of depth'
	    WRITE(6,*)' E(X):    Exit routine'
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'P')THEN
	    CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
!
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
1000	  CONTINUE
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
