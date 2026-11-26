!
! Program reads in population estimates from a file and uses
! linear interpolation in the log plane to lay these estimates
! on the "new" radius grid.
!
	SUBROUTINE REGRID_LOG_DC_V1(DHEN,R,ED,T,DI,CLUMP_FAC,EDGE,F_TO_S,INT_SEQ,
	1               POPATOM,N,ND,LU_IN,INTERP_OPTION,FILE_NAME)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Aletred 18-Nov-2023 - Fixed issue with TR option.
! Altered 21-Aug-2023 - Fixed issue with TR option and rounding.
! Altered 20-Aug-2022 - Error/warning messages improved for TR option.
! Altered 31-May-2022 - Added RNS option: R grid not scaled.
! Altered 06-Dec-2021 - Removed write ED, OLD_E from ED option.
! Altered 18-Aug-2019 - Added TR option. Designed to interpolate DC's in T in the inner region, and
!                           in T/R in the outer region (where it it is assumed T has not changed).
!                           In the inner region, T will be monotonic. The aim is to get reduce changes
!                              in high ionization species in the outer regions
!                              where the DC's are a very strong function of T.
! Altered 18-Jan-2014 - Fixed minor bug -- needed ABS when checkikng equality of R(1) and OLD_R(1).
! Altered 16-Oct-2012 - Added BA_TX_CONV and NO_TX_CONV.
! Altered: 18-Nov-2011: Big fix. Not all DPOP was being set and this could cause a crash
!                            when taking LOGS (levels unused).
! Created: 18-Dec-2010: This routine replace REGRIDWSC_V3, REGRIDB_ON_NE, REGRID_TX_R.
!                       Interplation options are 'R', 'ED', 'SPH_TAU', and 'RTX'.
!                       Routine handled INPUT files with departure coefficients, or
!                       with LOG(DCs).
!
	INTEGER N,ND
	INTEGER LU_IN
	REAL(KIND=LDP) DHEN(N,ND)
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) ED(ND)
	REAL(KIND=LDP) DI(ND)
	REAL(KIND=LDP) EDGE(N)
	REAL(KIND=LDP) POPATOM(ND)
	REAL(KIND=LDP) CLUMP_FAC(ND)
	INTEGER F_TO_S(N)
	INTEGER INT_SEQ(N)
	CHARACTER(LEN=*) INTERP_OPTION
!
	REAL(KIND=LDP) TAU(ND)
	REAL(KIND=LDP) NEW_ED(ND)
	REAL(KIND=LDP) NEW_X(ND)
	REAL(KIND=LDP) NEW_T(ND)
	REAL(KIND=LDP) LOG_TEN
!
	REAL(KIND=LDP), ALLOCATABLE :: DPOP(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_CLUMP_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_DI(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_R(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_T(:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_TAU(:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: OLD_X(:)
!
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU,LUER,LUWARN,WARNING_LU
	EXTERNAL ERROR_LU,WARNING_LU
!
! Local Variables.
!
	REAL(KIND=LDP) T_EXCITE,FX,DELTA_T
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TMAX,TMIN
	REAL(KIND=LDP) RTAU1,RTAU1_OLD
	REAL(KIND=LDP), PARAMETER :: RONE=1.0_LDP
	REAL(KIND=LDP), PARAMETER :: SIG_TH=6.65E-15_LDP 	!10^10 times Thompson ross-section
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER JMIN,JINT
	INTEGER I,J,K
	INTEGER NZ,NOLD,NDOLD
	INTEGER NX,NX_ST,NX_END
	INTEGER COUNT,IOS
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL BAD_TX_CONV
	LOGICAL NO_TX_CONV(ND)
	LOGICAL TAKE_LOGS
	LOGICAL CHECK_DC
	LOGICAL CLUMP_PRES
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=*) FILE_NAME
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
	IF(FIRST)THEN
	  WRITE(6,'(/,A)')' Using REGRID_LOG_DC_V1 to read in initial population estimates'
	  WRITE(6,'(2A,/)')' Using option: ',TRIM(INTERP_OPTION)
	END IF
!
! Read in values from previous model.
!
	OPEN(UNIT=LU_IN,STATUS='OLD',FILE=FILE_NAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILE_NAME),' in REGRID_LOG_DC_V1'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. CLUMP_FAC is presently not uses in
! this routine, as we regrid in R.
!
! We now assume base 10 for the input -- there was a brief period with base e.
!
	I=0
	STRING=' '
	DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(LU_IN,'(A)')STRING
	END DO
	CHECK_DC=.FALSE.
	CLUMP_PRES=.TRUE.
	IF( INDEX(STRING,'!Format date') .EQ. 0)THEN
	  REWIND(LU_IN)
	  CHECK_DC=.TRUE.
	  CLUMP_PRES=.FALSE.
	END IF
	IF( INDEX(STRING,'!Format date') .NE. 0 .AND. INDEX(STRING,'18-Nov-2010') .NE. 0)THEN
	  TAKE_LOGS=.FALSE.
	  LOG_TEN=RONE
	ELSE IF( INDEX(STRING,'!Format date') .NE. 0 .AND. INDEX(STRING,'10-Dec-2010') .NE. 0)THEN
	  TAKE_LOGS=.FALSE.
	  LOG_TEN=LOG(10.0_LDP)
	ELSE
	  TAKE_LOGS=.TRUE.
	END IF
!
! T1=RP_OLD
!
	READ(LU_IN,*)T1,T2,NOLD,NDOLD
!
! Allocate required memory.
!
	I=MAX(NDOLD,NOLD,N,ND)
	ALLOCATE (TA(I))
	ALLOCATE (TB(I))
!
	ALLOCATE (DPOP(NOLD,NDOLD))
	ALLOCATE (OLD_R(NDOLD))
	ALLOCATE (OLD_T(NDOLD))
	ALLOCATE (OLD_ED(NDOLD))
	ALLOCATE (OLD_DI(NDOLD))
	ALLOCATE (OLD_X(NDOLD))
	ALLOCATE (OLD_TAU(NDOLD))
	ALLOCATE (OLD_CLUMP_FAC(NDOLD)); OLD_CLUMP_FAC=RONE
!
! NZ defines the number of atomic levels which can be found by
! direct interpolation
!
	NZ=N
	IF(N .GT. NOLD)NZ=NOLD
	DO I=1,NDOLD
	  IF(CLUMP_PRES)THEN
	    READ(LU_IN,*)OLD_R(I),OLD_DI(I),OLD_ED(I),OLD_T(I),T1,T2,OLD_CLUMP_FAC(I)
	  ELSE
	    READ(LU_IN,*)OLD_R(I),OLD_DI(I),OLD_ED(I),OLD_T(I)
	  END IF
	  READ(LU_IN,*)(TA(J),J=1,NOLD)
	  DO J=1,NOLD                         !We use NOLD for when we take Logs.
	    DPOP(J,I)=TA(J)
	  END DO
	END DO
	CLOSE(LU_IN)
!
! Decide if DPOP refers to b, b-1 or log b. Convert all to log b (base e).
!
	IF( ABS( TA(NOLD) ) .LT. 0.2_LDP .AND. CHECK_DC)DPOP=DPOP+RONE
	IF(TAKE_LOGS)THEN
	  DPOP=LOG(DPOP)
	ELSE IF(LOG_TEN .NE. RONE)THEN
	  DPOP=DPOP*LOG_TEN
	END IF
!
	NX_ST=1
	NX_END=ND
	IF(INTERP_OPTION .EQ. 'R' .OR. INTERP_OPTION .EQ. 'TR')THEN
	  IF(ABS(OLD_R(NDOLD)/R(ND)-RONE) .GT. 0.0001_LDP)THEN
	    IF(FIRST)THEN
	      WRITE(LUWARN,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	      WRITE(LUER,*)'Additional warnings will be output to WARNINGS'
	    ELSE
	      WRITE(LUWARN,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	    END IF
	    DO I=1,NDOLD
	      OLD_R(I)=R(ND)*( OLD_R(I)/OLD_R(NDOLD) )
	    END DO
	    OLD_R(NDOLD)=R(ND)
	  ELSE
	    OLD_R(NDOLD)=R(ND)
	  END IF
	  IF( ABS(RONE-OLD_R(1)/R(1)) .LE. 1.0E-10_LDP )OLD_R(1)=R(1)
	  IF(OLD_R(2) .GE. OLD_R(1))THEN
	    WRITE(LUER,*)'Reset OLD_R(1) in REGRID_T_ED but now OLD_R(2) .GE. OLD_R(1))'
	    STOP
	  END IF
	END IF
!
	IF(INTERP_OPTION .EQ. 'R')THEN
	  OLD_X=LOG(OLD_R)
	  NEW_X=LOG(R)
	  DO WHILE(R(NX_ST) .GT. OLD_R(1))
	    NX_ST=NX_ST+1
	  END DO
	  NX=ND-NX_ST+1
!
	ELSE IF(INTERP_OPTION .EQ. 'RNS')THEN
	  OLD_X=LOG(OLD_R)
	  NEW_X=LOG(R)
	  DO WHILE(R(NX_ST) .GT. OLD_R(1))
	    NX_ST=NX_ST+1
	  END DO
          NX_END=ND
          DO WHILE (NEW_X(NX_END) .LT. OLD_X(NDOLD))
            NX_END=NX_END-1
          END DO
          NX=NX_END-NX_ST+1
	  WRITE(6,*)'RNS option:',NX_ST,NX_END,NX
	  FLUSH(UNIT=6)
!
	ELSE IF(INTERP_OPTION .EQ. 'RSP')THEN
	  IF(ABS(OLD_R(NDOLD)/R(ND)-RONE) .GT. 0.0001_LDP)THEN
	    IF(FIRST)THEN
	      WRITE(LUER,*)'Using RSP option'
	      WRITE(LUWARN,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Warning - core radius not identical in REGRID_LOG_DC_V1'
	      WRITE(LUER,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	      WRITE(LUER,*)'Using RSP option'
	    ELSE
	      WRITE(LUWARN,*)'Rescaling to make Rcore identical --- ',TRIM(FILE_NAME)
	    END IF
	    DO I=1,NDOLD
	      OLD_R(I)=R(ND)*( OLD_R(I)/OLD_R(NDOLD) )
	    END DO
	    OLD_R(NDOLD)=R(ND)
	  ELSE
	    OLD_R(NDOLD)=R(ND)
	  END IF
	  IF( ABS(RONE-OLD_R(1)/R(1)) .LE. 1.0E-10_LDP )OLD_R(1)=R(1)
	  IF(OLD_R(2) .GE. OLD_R(1))THEN
	    WRITE(LUER,*)'Reset OLD_R(1) in REGRID_T_ED but now OLD_R(2) .GE. OLD_R(1))'
	    STOP
	  END IF
!
	  IF( ABS(R(1)-OLD_R(1)) .LT. 1.0E-10_LDP )THEN
	    OLD_R(1)=R(1)
	  ELSE
	    T1=LOG(R(1)/OLD_R(1))/(OLD_R(1)/OLD_R(NDOLD)-1)
	    DO I=2,NDOLD-1
	      OLD_R(I)=OLD_R(I)*EXP(T1*(OLD_R(I)/OLD_R(1)-RONE))
	    END DO
	    OLD_R(1)=R(1)
	  END IF
	  NX_ST=1; NX=ND
	  OLD_X=LOG(OLD_R)
	  NEW_X=LOG(R)
!
	ELSE IF(INTERP_OPTION .EQ. 'ED')THEN
	  OLD_X=LOG(OLD_ED*OLD_CLUMP_FAC)
	  NEW_X=LOG(ED)
!	  CALL WRITV(NEW_X,ND,'ED new',6)
!	  CALL WRITV(OLD_X,NDOLD,'ED old',6)
          DO WHILE (NEW_X(NX_ST) .LT. OLD_X(1))
            NX_ST=NX_ST+1
          END DO
          NX_END=ND
          DO WHILE (NEW_X(NX_END) .GT. OLD_X(NDOLD))
            NX_END=NX_END-1
          END DO
          NX=NX_END-NX_ST+1
!
	ELSE IF(INTERP_OPTION .EQ. 'TR')THEN
!
! This code may have issues if T is very nomonotonic.
!
	  TMIN=MINVAL(OLD_T)
	  TMAX=MAXVAL(OLD_T)
	  JMIN=MINLOC(OLD_T,IONE)
	  OLD_DI=LOG(OLD_DI)
	  OLD_X=OLD_ED*OLD_CLUMP_FAC
	  DO I=1,ND
	    IF(T(I) .LT. TMIN .AND. ABS(LOG(OLD_X(JMIN)/ED(I))) .LT. RONE)THEN
	       DHEN(1:NZ,I)=DPOP(1:NZ,JMIN)
	       DI(I)=EXP(OLD_DI(JMIN))
	    ELSE
	      JINT=0	
	      DO J=1,NDOLD-1
	        IF( (T(I)-OLD_T(J))*(OLD_T(J+1)-T(I)) .GE. -1.0E-10_LDP)THEN
	          IF(JINT .EQ. 0)THEN
	            JINT=J
	          ELSE
	            IF( ABS(R(I)/OLD_R(J)-RONE) .LT.  ABS(R(I)/OLD_R(JINT)-RONE))JINT=J
	          END IF
	        END IF
	      END DO
	      IF(JINT .EQ. 0 .AND. T(I) .GE. TMAX)THEN
	          IF(FIRST .AND. T(I) .GT. 1.0001_LDP*TMAX)THEN
	            WRITE(6,'(/,A)')'Warning: T at inner boundary outside old model range in REGRID_LOG_DC_V1'
	            WRITE(6,'(2A,3(10X,A))')'    I',' JINT','T(I)','TMIN','TMAX'
	            WRITE(6,'(2I5,3E16.6)')I,JINT,T(I),TMIN,TMAX
	            FLUSH(UNIT=6)
	          END IF
	          DHEN(1:NZ,I)=DPOP(1:NZ,NDOLD)
	          DI(I)=EXP(OLD_DI(NDOLD))
	      ELSE IF(JINT .EQ. 0 .AND. T(I) .LE. TMIN)THEN
	          IF(FIRST .AND. T(I) .LT. 0.9999_LDP*TMIN)THEN
	            WRITE(6,'(/,A)')'Warning: T outside old model range in REGRID_LOG_DC_V1'
	            WRITE(6,'(2A,3(10X,A))')'    I',' JINT','T(I)','TMIN','TMAX'
	            WRITE(6,'(2I5,3E16.6)')I,JINT,T(I),TMIN,TMAX
	            FLUSH(UNIT=6)
	          END IF
	          DHEN(1:NZ,I)=DPOP(1:NZ,JMIN)
	          DI(I)=EXP(OLD_DI(JMIN))
	      ELSE IF(JINT .EQ. 0)THEN
	          WRITE(6,'(/,A)')'Error in REGRID_LOG_DC_V1 - JINT=0'
	          WRITE(6,'(2A,3(10X,A))')'    I',' JINT','T(I)','TMIN','TMAX'
	          WRITE(6,'(2I5,3E16.6)')I,JINT,T(I),TMIN,TMAX
	          FLUSH(UNIT=6)
	          STOP
	      ELSE
	        T1= ABS(OLD_T(JINT+1)/OLD_T(JINT)-RONE)
	        IF(T1 .GT. 1.0E-08_LDP)THEN
	          T1=(T(I)-OLD_T(JINT))/(OLD_T(JINT+1)-OLD_T(JINT))
	        ELSE
	          T1=RONE
	        END IF
	        DO K=1,NZ
	           DHEN(K,I)=T1*DPOP(K,JINT+1)+(RONE-T1)*DPOP(K,JINT)
	        END DO
	        DI(I)=EXP(T1*OLD_DI(JINT+1)+(RONE-T1)*OLD_DI(JINT))
	      END IF
!	      WRITE(6,*)'Done ',I; FLUSH(UNIT=6)
	    END IF
	  END DO
!	  WRITE(6,*)'Done interp'; FLUSH(UNIT=6)
	  GOTO 1000
!
	ELSE IF(INTERP_OPTION .EQ. 'SPH_TAU')THEN
!
!	  WRITE(6,*)'Interp option is SPH_TAU'
!	  DO I=1,ND
!	    WRITE(6,'(I5,2ES14.4)')I,OLD_ED(I)*OLD_CLUMP_FAC(I),ED(I)
!	  END DO
!
! Determine radius at which optical depth is unity. We use NEW_ED as a
! temporary vector for ED*CLUMP_FAC.
!
          NEW_ED(1:ND)=ED(1:ND)*CLUMP_FAC(1:ND)
          TAU(1)=SIG_TH*NEW_ED(1)*R(1)
          K=1
	  DO I=2,ND
            TAU(I)=TAU(I-1)+SIG_TH*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5_LDP
            IF(TAU(I) .LE. RONE)K=I
          END DO
	  IF(K .EQ. 1 .OR. K .EQ. ND)THEN
	    LUER=ERROR_LU()
            WRITE(LUER,*)'Error computing RTAU1 in REGRID_LOG_DC_V1'
          END IF
	  T1=(RONE-TAU(K))/(TAU(K+1)-TAU(K))
          RTAU1=T1*R(K+1)+(RONE-T1)*R(K)
!
! Compute the spherical optical depth scale. Assume atmosphere has constant
! V at outer boundary.
!
          DO I=1,ND
            NEW_ED(I)=NEW_ED(I)*RTAU1*RTAU1/R(I)/R(I)
          END DO
          TAU(1)=SIG_TH*NEW_ED(1)*R(1)/3.0_LDP
          DO I=2,ND
            TAU(I)=TAU(I-1)+SIG_TH*(NEW_ED(I-1)+NEW_ED(I))*(R(I-1)-R(I))*0.5_LDP
          END DO
!
! Determine radius at which optical depth is unity in old model.
!
	  OLD_ED=OLD_ED*OLD_CLUMP_FAC
          OLD_TAU(1)=SIG_TH*OLD_ED(1)*OLD_R(1)
          K=1
	  DO I=2,NDOLD
            OLD_TAU(I)=OLD_TAU(I-1)+SIG_TH*(OLD_ED(I-1)+OLD_ED(I))*(OLD_R(I-1)-OLD_R(I))*0.5_LDP
            IF(OLD_TAU(I) .LE. RONE)K=I
          END DO
	  IF(K .EQ. 1 .OR. K .EQ. NDOLD)THEN
	    LUER=ERROR_LU()
            WRITE(LUER,*)'Error computing RTAU1_OLD in REGRID_LOG_DC_V1'
          END IF
          T1=(RONE-OLD_TAU(K))/(OLD_TAU(K+1)-OLD_TAU(K))
          RTAU1_OLD=T1*OLD_R(K+1)+(RONE-T1)*OLD_R(K)
!
! Compute the spherical optical depth scale.
!
          DO I=1,NDOLD
            OLD_ED(I)=OLD_ED(I)*(RTAU1_OLD/OLD_R(I))**2
          END DO
          OLD_TAU(1)=SIG_TH*OLD_ED(1)*OLD_R(1)/3.0_LDP
          DO I=2,NDOLD
            OLD_TAU(I)=OLD_TAU(I-1)+SIG_TH*(OLD_ED(I-1)+OLD_ED(I))*(OLD_R(I-1)-OLD_R(I))*0.5_LDP
          END DO
!
!	  WRITE(6,*)RTAU1,RTAU1_OLD
!	  DO I=1,NDOLD
!	    WRITE(6,'(I5,4ES14.4)')I,OLD_TAU(I),TAU(I),OLD_ED(I),NEW_ED(I)
!	  END DO
!
! Detrmine whether new mesh extends beyond oldmesh.
! NX is used to define the region over which we may use a linear
! interpolation. Beyond the old radius mesh we set the departure coefficients
! constant or equal to 1.
!
	  NX_ST=1
	  DO WHILE (TAU(NX_ST) .LT. OLD_TAU(1))
	    NX_ST=NX_ST+1
	  END DO
	  NX_END=ND
	  DO WHILE (TAU(NX_END) .GT. OLD_TAU(NDOLD))
	    NX_END=NX_END-1
	  END DO
	  NX=NX_END-NX_ST+1
!
	  OLD_X=LOG(OLD_TAU)
	  NEW_X=LOG(TAU)
	ELSE
	  WRITE(6,*)'Interp option not recognized:',TRIM(INTERP_OPTION)
	  STOP
	END IF
!
! Interpolate the ion density.
!
	DO I=1,NDOLD
	  OLD_DI(I)=LOG(OLD_DI(I))
	END DO
	CALL LINPOP(NEW_X(NX_ST),DI(NX_ST),NX,OLD_X,OLD_DI,NDOLD)
	DO I=NX_ST,NX_END
	  DI(I)=EXP(DI(I))
	END DO
!
! For locations outside grid, we assume ion density is fixed relative to the atom density.
!
        DO I=1,NX_ST-1
          DI(I)=POPATOM(I)*DI(NX_ST)/POPATOM(NX_ST)
        END DO
	DO I=NX_END+1,ND
	  DI(I)=POPATOM(I)*DI(NX_END)/POPATOM(NX_END)
	END DO
!
! Compute LOG(excitation temperatures), if INTERP_OPTION is RTX.
!
	IF(INTERP_OPTION .EQ. 'RTX')THEN
	  DO I=1,NZ
	    T_EXCITE=OLD_T(NDOLD)
	    DO J=NDOLD,1,-1
	      DELTA_T=100
	      DO WHILE(ABS(DELTA_T/T_EXCITE) .GT. 1.0E-08_LDP)
	        FX=DPOP(I,J)*EXP(HDKT*EDGE(I)*(RONE/OLD_T(J)-RONE/T_EXCITE))*
	1          (T_EXCITE/OLD_T(J))**1.5_LDP
	        DELTA_T=(FX-RONE)*T_EXCITE/FX/(1.5_LDP+HDKT*EDGE(I)/T_EXCITE)
	        IF(DELTA_T .GT.  0.8_LDP*T_EXCITE)DELTA_T=0.8_LDP*T_EXCITE
	        IF(DELTA_T .LT. -0.8_LDP*T_EXCITE)DELTA_T=-0.8_LDP*T_EXCITE
	        T_EXCITE=T_EXCITE-DELTA_T
	      END DO
	      DPOP(I,J)=LOG(T_EXCITE)
	    END DO
	  END DO
	END IF
!
	IF(INTERP_OPTION .EQ. 'RSP')THEN
	  CALL LINPOP(OLD_R,TA,NDOLD,R,T,ND)
	  DO J=1,NDOLD
	    IF(TA(J) .NE. TA(J))WRITE(6,*)'TA error',TA(J)
	    IF(T(J) .NE. T(J))WRITE(6,*)'T error',T(J)
	    DO I=1,NZ
	      IF(ABS(DPOP(I,J)) .GT. 0.3_LDP)THEN
	        DPOP(I,J)=DPOP(I,J)+1.5_LDP*(TA(J)/OLD_T(J))+HDKT*EDGE(I)*(RONE/OLD_T(J)-RONE/TA(J))
	      END IF
	    END DO
	  END DO
	END IF
!
! Interpolate the departure coefficients / excitation temperatures.
! As DPOP is in LOG plane, there is no need to take LOGS.
!
	DO J=1,NZ
	  DO I=1,NDOLD
	    TB(I)=DPOP(J,I)
	  END DO
	  CALL LINPOP(NEW_X(NX_ST),TA(NX_ST),NX,OLD_X,TB,NDOLD)
	  DO I=NX_ST,NX_END
	    DHEN(J,I)=TA(I)
	  END DO
	END DO
!
! Convert, if necessary, back from T_EXCITE to log(b).
!
	IF(INTERP_OPTION .EQ. 'RTX')THEN
	  DO I=1,NDOLD
	    TA(I)=LOG(OLD_T(I))
	  END DO
	  CALL LINPOP(NEW_X(NX_ST),NEW_T(NX_ST),NX,OLD_X,TA,NDOLD)
	  DO I=NX_ST,NX_END
	    NEW_T(I)=EXP(NEW_T(I))
	  END DO
          DO I=NX_ST,NX_END
            T1=T(I)+(EXP(DHEN(J,I))-NEW_T(I))
            DHEN(J,I)=HDKT*EDGE(J)*(RONE/T1-RONE/T(I)) + 1.5_LDP*LOG(T(I)/T1)
          END DO
	END IF
!
! When extending to larger radii, we assume that the ionization
! state of the gas does not change. This supersedes the earlier assumption
! of a fixed departure coefficient. Excited populations are adjusted
! logarithmically. Code assumes ED is proportional to POP_ATOM, and
! hence may break down for SN models (with variable composition).
!
	DO I=1,NX_ST-1
	  DO J=1,NZ
	    T1=POPATOM(NX_ST)/POPATOM(I)
	    T2=MIN(RONE, ABS(DHEN(J,NX_ST))/ABS(DHEN(1,NX_ST)) )
	    DHEN(J,I)=DHEN(J,NX_ST)+T2*LOG(T1)
	  END DO
	END DO
	DO I=NX_END+1,ND
	  DO J=1,NZ
	    DHEN(J,I)=DHEN(J,NX_END)
	  END DO
	END DO
!
1000	CONTINUE
!
! Compute departure coefficients for N>NZ. These levels are set to have these
! same excitation temperature as the highest level.
!
	NO_TX_CONV=.FALSE.
	BAD_TX_CONV=.FALSE.
	IF(N .GT. NZ)THEN
	  T_EXCITE=T(ND)
	  DO I=ND,1,-1
!
! We first compute the excitation temperature on level NZ.
!
	    DELTA_T=100.0_LDP
	    COUNT=0
!	    DO WHILE(ABS(DELTA_T) .GT. 1.0D-08 .AND. COUNT .LT. 100)
	    DO WHILE(ABS(DELTA_T/T_EXCITE) .GT. 1.0E-04_LDP .AND. COUNT .LT. 100)
	      FX=EXP(DHEN(NZ,I)+HDKT*EDGE(NZ)*(RONE/T(I)-RONE/T_EXCITE))*(T_EXCITE/T(I))**1.5_LDP
	      DELTA_T=(FX-RONE)*T_EXCITE/FX/(1.5_LDP+HDKT*EDGE(NZ)/T_EXCITE)
	      IF(DELTA_T .GT.  0.8_LDP*T_EXCITE)DELTA_T=0.8_LDP*T_EXCITE
	      IF(DELTA_T .LT. -0.8_LDP*T_EXCITE)DELTA_T=-0.8_LDP*T_EXCITE
	      T_EXCITE=T_EXCITE-DELTA_T
	      COUNT=COUNT+1
	      IF(COUNT .GE. 60)THEN
	        WRITE(6,*)I,COUNT,T(I),T_EXCITE,DELTA_T
	        FLUSH(UNIT=6)
	      END IF
	    END DO
!
! We can now compute the Departure coefficients.
!
	    IF(COUNT .GE. 100)THEN
	      NO_TX_CONV(I)=.TRUE.
	      BAD_TX_CONV=.TRUE.
	      DO J=NZ+1,N
	        DHEN(J,I)=DHEN(NZ,I)
	      END DO
	    ELSE
	      DO J=NZ+1,N
	        DHEN(J,I)=HDKT*EDGE(J)*(RONE/T_EXCITE-RONE/T(I))+1.5_LDP*LOG(T(I)/T_EXCITE)
	      END DO
	    END IF
	  END DO
	  IF(BAD_TX_CONV)THEN
	    WRITE(LUER,*)'Warning in REGRID_LOG_DC_V1 - TX did not converge in 100 iterations'
	    WRITE(LUER,*)'This might affect convergence. File is ',TRIM(FILE_NAME)
	    WRITE(LUER,*)'It occurred at the following depths:'
	    DO I=1,ND
	      IF(NO_TX_CONV(I))WRITE(LUER,'(I5)',ADVANCE='NO')I
	      IF(MOD(I,15) .EQ. 0)WRITE(LUER,'(A)')' '
	    END DO
          END IF
	END IF
!
	CLOSE(UNIT=8)
!
! Ensure that all levels belonging to the same level have the
! same DC. We take the DC of the lowest state.
!
	IF(SUM(INT_SEQ) .EQ. 0)THEN
	  DO I=1,ND
	    TA(:)=0.0_LDP
	    DO J=1,N
	      IF(TA(F_TO_S(J)) .EQ. 0.0_LDP)TA(F_TO_S(J))=DHEN(J,I)
	    END DO
	    DO J=1,N
	      DHEN(J,I)=TA(F_TO_S(J))
	    END DO
	  END DO
	END IF
!
	WRITE(170,'(A)')FILE_NAME
	WRITE(170,'(5ES16.6)')(DHEN(1:N,1))
!
! Free up memory.
!
	DEALLOCATE (DPOP)
	DEALLOCATE (OLD_CLUMP_FAC)
	DEALLOCATE (OLD_ED)
	DEALLOCATE (OLD_DI)
	DEALLOCATE (OLD_R)
	DEALLOCATE (OLD_T)
	DEALLOCATE (OLD_TAU)
	DEALLOCATE (OLD_X)
	DEALLOCATE (TA)
	DEALLOCATE (TB)
!
	FIRST=.FALSE.
	RETURN
	END
