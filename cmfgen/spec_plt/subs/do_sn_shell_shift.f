!
! Subroutine designed to offset SN shells in V as a function of Beta. The
! shells are (generally) shifted radomly in V by an amount dV or by a
! constant scale factor dV/V. Two main output files (plus some diagnostic 
! files) are created:
!
! As this routines use V, it only works (as written) for SN models
!    with homologous expansion.
!
!    3D_DATA -- Large file contaning R, V, SIGMA, T, ESEC, CHI, ETA
!                   R is the same for all beta
!                   VR, SIGMA, T, and ESEC are functions of R, BETA
!                   ETA and CHI -- function of R, BETA and frequency
!                   ETA -- no scattering contribution.
!
!    LINE_MOM_DATA -- Primarily used for the 5 polarization moments.
!                     J_L and J_R are assumed to be 0.5J
!                     K= J/3. Other moments are zero. Cannot be used for
!                       polarization calculations.
! 
! Requires RVTJ, ETA_DATA, CHI_DATA and ES_J_CONV data if both output
! files are being output. 3D_DATA and LINE_MOM_DATA do not need to be
! created.
! 
	MODULE SN_SHIFT_MODULE
	USE SET_KIND_MODULE
!
! Cleaned: 01-May-2024
!
	INTEGER ND_MULTI
	INTEGER NB_MULTI
	INTEGER NCF_MULTI
!
	REAL(KIND=LDP), ALLOCATABLE :: R_MULTI(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: BETA_MULTI(:)
	REAL(KIND=LDP), ALLOCATABLE :: FREQ_MULTI(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: ETA_MULTI(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_MULTI(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: J_MULTI(:,:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: VR_MULTI(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_MULTI(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TEMP_MULTI(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_MULTI(:,:)
!
        REAL(KIND=LDP), ALLOCATABLE :: dV_BETA(:)
        INTEGER, ALLOCATABLE :: NS_VEC(:)
        INTEGER, ALLOCATABLE :: IST_VEC(:)
        INTEGER, ALLOCATABLE :: IEND_VEC(:)
	LOGICAL LOG_V_SHIFT
!
	PUBLIC DO_SN_SHELL_SHIFT
	PRIVATE CREATE_SHIFTED_SN_GRID, COMPUTE_BETA, DO_SHIFT
	PRIVATE WRITE_MULTI,  WRITE_JH_POL
	CONTAINS
!
	SUBROUTINE DO_SN_SHELL_SHIFT(R,T,V,SIGMA,ESEC,CLUMP_FAC,ND,NC,
	1           ETA, CHI, OPAC_NU, NCF,
	1           ES_J, ES_NU, ES_V, ES_ND, ES_NCF)
	USE MOD_USR_LIST_OPT
	IMPLICIT NONE
!
! Altered: 05-Oct-2026: Reject non-monotonic V; interpolation uses V as a coordinate.
! Altered: 29-May-2024: CLUMP_FAC added to call.
!
	INTEGER ND
	INTEGER NC
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) CLUMP_FAC(ND)
!
	INTEGER NCF
	REAL(KIND=LDP) ETA(ND,NCF)
	REAL(KIND=LDP) CHI(ND,NCF)
	REAL(KIND=LDP) OPAC_NU(NCF)
!
	INTEGER ES_ND
	INTEGER ES_NCF
	REAL(KIND=LDP) ES_J(ES_ND,ES_NCF)
	REAL(KIND=LDP) ES_NU(ES_ND,ES_NCF)
	REAL(KIND=LDP) ES_V(ES_ND)
!
! Local variables
!
	REAL(KIND=LDP) dV
	REAL(KIND=LDP) OP_SCL_EXPONENT
	REAL(KIND=LDP) LAM_START
	REAL(KIND=LDP) LAM_END
	REAL(KIND=LDP) ES_RES_KMS
!
	INTEGER NB
	INTEGER NPHI
	CHARACTER(LEN=20) BETA_LAW
	LOGICAL TOP_BOT_SYM
	LOGICAL DO_ETA_SHIFT
	LOGICAL DO_J_SHIFT
	LOGICAL USE_OLD_SEED 
!
! This utility shifts the model in velocity space and interpolates T,
! ESEC, SIGMA, ETA, and CHI on V.  R(V) must therefore be single-valued;
! reject a supplied non-monotonic field instead of silently choosing one
! radial branch.
!
	IF(ANY(V(1:ND) .LE. 0.0_LDP) .OR.
	1      ANY(V(1:ND-1) .LE. V(2:ND)))THEN
	  WRITE(6,'(A)')' Error in DO_SN_SHELL_SHIFT'
	  WRITE(6,'(A)')' V must be positive and strictly decreasing with depth'
	  WRITE(6,'(A)')' Velocity-space shell shifts require a homologous/monotonic field'
	  STOP
	END IF
!
! We use USE_LIST_OPT as input/output is clearer when we have a large number of options.
!
	CALL USR_LIST_OPT(dV,'dV','200.0','Velocity shift (km/s) (-ve for scaled shift [e.g., 1.1])')
	CALL USR_LIST_OPT(OP_SCL_EXPONENT,'OP_EXP','3.0','Exponent to scale op. when shifting in V')
	CALL USR_LIST_OPT(NB,'NB','11','Number of angles in polar grid')
	CALL USR_LIST_OPT(BETA_LAW,'BETA_LAW','UNIFORM_BETA','Beta distribution - UNIFORM_BETA, UNIFORM_COSB, POW, BMIN')
	CALL USR_LIST_OPT(TOP_BOT_SYM,'TOP_BOT_SYM','T','Top-bottom symmetric?')
	CALL USR_LIST_OPT(LAM_START,'LAM_ST','400.0', 'Minimum wavelength of 2D grid (Ang)')
	CALL USR_LIST_OPT(LAM_END,'LAM_END','7000.0','Maximum wavelength of 2D grid (Ang)')
	CALL USR_LIST_OPT(ES_RES_KMS,'ES_RES_KMS','200.0','Resolution (km/s) to sample elec. scat. J')
	CALL USR_LIST_OPT(NPHI,'NPHI','11','Number of azimuthal angles (3,7, 11, 15 etc')
	CALL USR_LIST_OPT(DO_ETA_SHIFT,'DO_ETA_SHIFT','T','Shift ETA and output 3D_DATA?')
	CALL USR_LIST_OPT(DO_J_SHIFT,'DO_J_SHIFT','T','Shift H and output LINE_MOM_DATA?')
	CALL USR_LIST_OPT(USE_OLD_SEED,'OLD_SEED','F','Use old see for RAN generator')
!
        CALL CREATE_SHIFTED_SN_R_GRID(R,V,dV,USE_OLD_SEED,ND,NB)
        CALL COMPUTE_BETA(BETA_LAW,TOP_BOT_SYM,NB)
        CALL DO_SHIFT(R,V,SIGMA,T,ESEC,CLUMP_FAC,ETA,CHI,OPAC_NU,LAM_START,LAM_END,dV,
	1                 OP_SCL_EXPONENT,DO_ETA_SHIFT,ND,NB,NCF)
        IF(DO_ETA_SHIFT)CALL WRITE_MULTI()
!
	IF(DO_J_SHIFT)THEN
           CALL WRITE_JH_POL(ES_J,ES_NU,ES_V,LAM_START,LAM_END,ES_RES_KMS,ES_NCF,NC,ES_ND,NPHI)
        END IF
!
	RETURN
	END SUBROUTINE DO_SN_SHELL_SHIFT
!
CONTAINS
	SUBROUTINE COMPUTE_BETA(BETA_LAW,TOP_BOT_SYM,NB)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Revised: 15-Apr-2024 _ BETA_LAW,TOP_BOT_SYM included in call.
!
! Input values
!
	INTEGER NB
	LOGICAL TOP_BOT_SYM
	CHARACTER(LEN=*) BETA_LAW 
!
	INTEGER, PARAMETER :: NPAR_MAX=10
	INTEGER NPAR
	REAL(KIND=LDP) PARAMS(NPAR_MAX)
!
! Local variables
!
	REAL(KIND=LDP) PI,T1,T2
	REAL(KIND=LDP) dBETA_MIN, dBETA_MAX
	INTEGER I,K,LU_OUT
!
	NB_MULTI=NB
	IF(MOD(NB_MULTI,2) .NE. 1)THEN
	  WRITE(6,*)'Error in COMPUTE_BETA_CL_V2'
	  WRITE(6,*)'NBETA_MULTI must be ODD'
	  STOP
	END IF
!
!
! Two options are present to specify the azimuth angle Beta.
! (i) BETA uniformly spaced from 0 to PI/2 (0 to PI if not TOP_BOT_SYM).
! (i) COS_BETA uniformly spaced from 1 to 0 (-1 to 1 if not TOP_BOT_SYM).
!
! Interpolation section in LONG_CHAR may need to be changed, if BETA selction
! is changed.
!
! NB. Make sure this is a double precision calculation.
!
	IF(ALLOCATED(BETA_MULTI))DEALLOCATE(BETA_MULTI)
	ALLOCATE(BETA_MULTI(NB))
	PI=4.0*ATAN(1.0D0)
	IF(BETA_LAW .EQ. 'UNIFORM_BETA')THEN
          IF(TOP_BOT_SYM)THEN
            T1=0.5D0*PI/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI
              BETA_MULTI(I)=(I-1)*T1
            END DO
          ELSE
            T1=PI/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI
              BETA_MULTI(I)=(I-1)*T1
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'UNIFORM_COSB')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI
              T2=1.0D0-(I-1)*T1
              BETA_MULTI(I)=ACOS(T2)
            END DO
          ELSE
            T1=2.0D0/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI
              T2=1.0D0-(I-1)*T1
              BETA_MULTI(I)=ACOS( T2 )
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'POWER')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI
              T2=((I-1)*T1)**PARAMS(1)
              BETA_MULTI(I)=0.5*PI*T2
            END DO
          ELSE
            T1=2.0D0/(NB_MULTI-1.0D0)
            DO I=1,NB_MULTI/2+1
              T2=((I-1)*T1)**PARAMS(1)
              BETA_MULTI(I)=0.5*PI*T2
	      BETA_MULTI(NB_MULTI-I+1)=PI-BETA_MULTI(I)
            END DO
          END IF
        ELSE IF(BETA_LAW .EQ. 'BMIN')THEN
          IF(TOP_BOT_SYM)THEN
            T1=1.0D0/(NB_MULTI-1.0D0)
	    K=NINT(PARAMS(1))
            DO I=1,K
              BETA_MULTI(I)=(I-1)*PARAMS(2)
	    END DO
	    T2=(0.5*PI-BETA_MULTI(K))/(NB_MULTI-K)
            DO I=K+1,NB_MULTI
              T1=(I-K)*T2
              BETA_MULTI(I)=BETA_MULTI(K)+T1
            END DO
          ELSE
	  END IF
        ELSE 
	  WRITE(6,*)'Unrecognized law for choosing Beta in COMPUTE_BETA'
	  WRITE(6,*)'BETA_LAW=',BETA_LAW
          STOP
	END IF
!
	CALL GET_LU(LU_OUT,'COMPUTE_BETA')
	OPEN(UNIT=LU_OUT,FILE='BETA_GRID',STATUS='UNKNOWN')
	dBETA_MIN=100.0D0
	dBETA_MAX=-100.0D0
	WRITE(6,'(A7,2A14)')'Index','Beta','dBETA'
	DO I=1,NB_MULTI-1
	  WRITE(LU_OUT,'(I7,2ES14.4)')I,BETA_MULTI(I),BETA_MULTI(I+1)-BETA_MULTI(I)
	  dBETA_MIN=MIN(dBETA_MIN,BETA_MULTI(I+1)-BETA_MULTI(I))
	  dBETA_MAX=MAX(dBETA_MAX,BETA_MULTI(I+1)-BETA_MULTI(I))
	END DO
	WRITE(LU_OUT,'(I7,3ES14.4)')NB_MULTI,BETA_MULTI(NB_MULTI)
	CLOSE(UNIT=LU_OUT)
!
	WRITE(6,*)' '
	WRITE(6,*)'Grid written to BETA_GRID'
	WRITE(6,*)'dBETA (min and max) are', dBETA_MIN, dBETA_MAX
	WRITE(6,*)'BETA(1), BETA(NMAX) are', BETA_MULTI(1), BETA_MULTI(NB_MULTI)
!
	RETURN
	END SUBROUTINE COMPUTE_BETA
!
	SUBROUTINE CREATE_SHIFTED_SN_R_GRID(R,V,dV,USE_OLD_SEED,ND,NB)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NB
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) dV
	LOGICAL USE_OLD_SEED
!
! Local variables
!
	REAL(KIND=LDP) VSCL_FAC
	REAL(KIND=LDP) VMAX
	REAL(KIND=LDP) VMIN
	REAL(KIND=LDP) RAND_NUM
!
	REAL(KIND=LDP) dR
	REAL(KIND=LDP) dR_OUT 
	REAL(KIND=LDP) dR_IN
	REAL(KIND=LDP) T1
!
	REAL(KIND=LDP), ALLOCATABLE :: VR(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_V(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_WRK(:)
	INTEGER, ALLOCATABLE ::  ILOC(:)
!
	INTEGER NSEED
	INTEGER, ALLOCATABLE :: ISEED(:)

	INTEGER NW
	INTEGER NS
	INTEGER ID,IB
	INTEGER IST,IEND
	INTEGER K,L,LUOUT
	INTEGER NDM2
!
	INTEGER, PARAMETER :: IONE=1
!
	NB_MULTI=NB
	IF(ALLOCATED(NS_VEC))THEN
	  DEALLOCATE(NS_VEC,IST_VEC,IEND_VEC,dV_BETA)
	END IF
	ALLOCATE(NS_VEC(NB))
	ALLOCATE(IST_VEC(NB))
	ALLOCATE(IEND_VEC(NB))
	ALLOCATE(dV_BETA(NB))
!
	CALL GET_LU(LUOUT,'CREATE_SHIFTED_R_GRID')
	CALL RANDOM_SEED(SIZE = NSEED)
	IF(ALLOCATED(ISEED))DEALLOCATE(ISEED)
	ALLOCATE(ISEED(NSEED))
	IF(USE_OLD_SEED)THEN
	  OPEN(UNIT=LUOUT,FILE='OLD_RANDOM_SEED',STATUS='OLD',ACTION='READ')
	    READ(LUOUT,*)ISEED
	  CLOSE(LUOUT)
	ELSE
	  CALL RANDOM_SEED()
	  CALL RANDOM_SEED(GET=ISEED)
	  OPEN(UNIT=LUOUT,FILE='OLD_RANDOM_SEED',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(LUOUT,*)ISEED
	  CLOSE(LUOUT)
	END IF
!
! We omit depths 2 and ND-2 from the V grid as these are usually have very small
! spacings.
!
	NDM2=ND-2
	IF(ALLOCATED(NEW_V))DEALLOCATE(NEW_V,ILOC,V_WRK,VR)
	ALLOCATE(NEW_V(ND*NB))
	ALLOCATE(ILOC(NB))
	ALLOCATE(V_WRK(1:NDM2))
	ALLOCATE(VR(NDM2,NB))
	V_WRK(1)=V(1); V_WRK(2:NDM2-1)=V(3:NDM2); V_WRK(NDM2)=V(ND)
!
	IF(dV .GE. 0.0)THEN
	  DO IB=1,NB
	    CALL RANDOM_NUMBER(RAND_NUM)
	    dV_BETA(IB)=dV*(2*RAND_NUM-1.0_LDP)
	    IF(IB .EQ. 1)dV_BETA(IB)=0.0_LDP
	    DO ID=1,NDM2
	      VR(ID,IB)=V_WRK(ID)+dV_BETA(IB)
	    END DO
	  END DO
	  LOG_V_SHIFT=.FALSE.
	ELSE
	  T1=DLOG(-dV)
	  DO IB=1,NB
	    CALL RANDOM_NUMBER(RAND_NUM)
	    dV_BETA(IB)=EXP(T1*(2*RAND_NUM-1.0_LDP))
	    IF(IB .EQ. 1)dV_BETA(IB)=1.0_LDP
	    DO ID=1,NDM2
	      VR(ID,IB)=V_WRK(ID)*dV_BETA(IB)
	    END DO
	  END DO
	  LOG_V_SHIFT=.TRUE.
	END IF
!
	ID=1
	ILOC(1:NB)=1
	VMIN=MINVAL(VR)
	VMAX=MAXVAL(VR)
	NEW_V(1)=VMAX
!
	WRITE(6,*)' '
	WRITE(6,*)'VMIN,VMAX:',MINVAL(V),MAXVAL(V),' (original V grid)'
	WRITE(6,*)'VMIN,VMAX:',VMIN,VMAX,' (shifted V grid)'
!
! Define the new velcity grid, using the old velocity grid to set the
! step sizes in V.
!
	DO WHILE(NEW_V(ID) .GT. VMIN)
	  dV=VMAX
	  DO IB=1,NB
	    L=ILOC(IB)
	    IF(L .LE. NDM2)THEN
	      DO WHILE(VR(L,IB) .GT. NEW_V(ID) .AND. L .LT. NDM2)
	        L=L+1
	        ILOC(IB)=L
	        IF(L .GE. NDM2)EXIT
	      END DO
	      T1=VMAX
	      IF(L .NE. 1)T1=VR(L-1,IB)-VR(L,IB)
	      dV=MIN(dV,T1)
	      IF(L .NE. NDM2)T1=VR(L,IB)-VR(L+1,IB)
	      dV=MIN(dV,T1)
	    END IF
	  END DO
	  ID=ID+1
	  NEW_V(ID)=NEW_V(ID-1)-dV
	  IF(NEW_V(ID) .LT. VMIN)THEN
	    NEW_V(ID)=VMIN
	    EXIT
	  END IF
	END DO
	NW=ID
!
	WRITE(6,'(/,A,/)')' New velocity grid written to NEW_V_GRID'
	OPEN(LUOUT,FILE='NEW_V_GRID',STATUS='UNKNOWN')
	WRITE(LUOUT,'(A,2A15)')' Depth','V(km/s)','dV(km/s)'
	DO ID=1,NW-1
	  WRITE(LUOUT,'(I6,2ES15.6)')ID,NEW_V(ID),NEW_V(ID)-NEW_V(ID+1)
	END DO
	WRITE(LUOUT,'(I6,2ES15.6)')NW,NEW_V(NW)
	CLOSE(LUOUT)
!
	IF(ALLOCATED(VR_MULTI))DEALLOCATE(VR_MULTI,R_MULTI)
	ALLOCATE(VR_MULTI(NW,NB),R_MULTI(NW))
	DO IB=1,NB_MULTI
	  VR_MULTI(1:NW,IB)=NEW_V(1:NW)
	END DO
        IST=1
!
	T1=1.0_LDP-1.0E-10_LDP
        DO WHILE(NEW_V(IST) .GT. T1*V(1))
           IST=IST+1
        END DO
!
        IEND=NW-1
	T1=1.0_LDP+1.0E-10_LDP
        DO WHILE(NEW_V(IEND) .LT. T1*V(ND))
           IEND=IEND-1
        END DO
!
	DO ID=1,IST-1
	  R_MULTI(ID)=R(1)*NEW_V(ID)/V(1)
	END DO
	DO ID=IEND+1,NW
	  R_MULTI(ID)=R(ND)*NEW_V(ID)/V(ND)
	END DO
	NS=IEND-IST+1
	CALL MON_INTERP_FAST(R_MULTI(IST),NS,IONE,NEW_V(IST),NS,R,ND,V,ND,'CREATE_SHIFTED_R_GRID')	
	ND_MULTI=NW
!
! Clean up allocated arrays
!
	DEALLOCATE(V_WRK,VR,NEW_V,ILOC)
!
! We can only use interpolation inside the original V boundaries. Below we determne what 
! range this corresponds to in the new angular dependent grid. Thes will be used in
! routines such as DO_SHIFT to do the shifting.
!
	ALLOCATE(V_WRK(ND))
	WRITE(6,'(3A5,4A14)')'IB','IST','IEND','VR_MULTI(IST)','V_WRK(1)','V(1)','dV_BETA(IB)'
	DO IB=1,NB
	  IF(LOG_V_SHIFT)THEN
	    V_WRK(1:ND)=V(1:ND)*dV_BETA(IB)
	  ELSE
	    V_WRK(1:ND)=V(1:ND)+dV_BETA(IB)
	  END IF
!
	  IST=1; T1=1.0_LDP-1.0E-10_LDP
	  DO WHILE(VR_MULTI(IST,IB) .GT. T1*V_WRK(1))
	    IST=IST+1
	  END DO
!
	  IEND=NW; T1=1.0_LDP+1.0E-10_LDP
	  DO WHILE(VR_MULTI(IEND,IB) .LE. T1*V_WRK(ND))
	   IEND=IEND-1
	  END DO 
	  NS_VEC(IB)=IEND-IST+1
	  IST_VEC(IB)=IST
	  IEND_VEC(IB)=IEND
	  WRITE(6,'(3I5,4ES14.4)')IB,IST,IEND,VR_MULTI(IST,IB),V_WRK(1),V(1),dV_BETA(IB)
	END DO
	WRITE(6,*)' '
!
! Clean up allocated arrays
!
	DEALLOCATE(V_WRK)
!
	WRITE(6,*)'Exiting CREATE_SHIFTED_SN_R_GRID'
!
	RETURN
	END SUBROUTINE CREATE_SHIFTED_SN_R_GRID


	SUBROUTINE DO_SHIFT(R,V,SIGMA,T,ESEC,CLUMP_FAC,ETA,CHI,NU,LAM_ST,LAM_END,dV,
	1            OP_SCL_EXPONENT,DO_ETA_SHIFT,ND,NB,NCF)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NB
	INTEGER NCF
!
	REAL(KIND=LDP) R(ND)
	REAL(KIND=LDP) V(ND)
	REAL(KIND=LDP) T(ND)
	REAL(KIND=LDP) SIGMA(ND)
	REAL(KIND=LDP) ESEC(ND)
	REAL(KIND=LDP) CLUMP_FAC(ND)
!
	REAL(KIND=LDP) NU(NCF)
	REAL(KIND=LDP) ETA(ND,NCF)
	REAL(KIND=LDP) CHI(ND,NCF)
!
	REAL(KIND=LDP) LAM_ST
	REAL(KIND=LDP) LAM_END
	REAL(KIND=LDP) dV
	REAL(KIND=LDP) OP_SCL_EXPONENT
	LOGICAL DO_ETA_SHIFT
!
! Local data
!
	REAL(KIND=LDP) MIN_FREQ
	REAL(KIND=LDP) MAX_FREQ
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) ALPHA
!
	REAL(KIND=LDP) RAND_NUM
!
	REAL(KIND=LDP), ALLOCATABLE :: OLD_ETA(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: OLD_CHI(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_SCL(:)
!
	INTEGER NS,NW
	INTEGER I, J, ML
	INTEGER ML_ST,ML_END
	INTEGER IST,IEND
	INTEGER IB
!
	INTEGER, PARAMETER ::IONE=1
	CHARACTER(LEN=80) TMP_STR,CURVE_LAB
!
	NW=ND_MULTI
	NB_MULTI=NB
	MIN_FREQ=2.99792E+03_LDP/MAX(LAM_ST,LAM_END)
	MAX_FREQ=2.99792E+03_LDP/MIN(LAM_ST,LAM_END)
!
	ML_ST=1
	DO ML=1,NCF
	  IF(NU(ML) .LT. MAX_FREQ)EXIT
	  ML_ST=ML
	END DO
	ML_END=NCF
	DO ML=NCF,1,-1
	  IF(NU(ML) .GT. MIN_FREQ)EXIT
	  ML_END=ML
	END DO
	NCF_MULTI=ML_END-ML_ST+1
!
	IF(ALLOCATED(FREQ_MULTI))THEN
	  DEALLOCATE(FREQ_MULTI,SIGMA_MULTI,TEMP_MULTI,ESEC_MULTI)
	END IF
!
	IF(ALLOCATED(ETA_MULTI))THEN
	  DEALLOCATE(ETA_MULTI,CHI_MULTI)
	END IF
!
	ALLOCATE(FREQ_MULTI(NCF_MULTI)); FREQ_MULTI(1:NCF_MULTI)=NU(ML_ST:ML_END)
	ALLOCATE(SIGMA_MULTI(NW,NB))
	ALLOCATE(TEMP_MULTI(NW,NB))
	ALLOCATE(ESEC_MULTI(NW,NB))
	IF(DO_ETA_SHIFT)THEN
	  ALLOCATE(ETA_MULTI(NW,NB,NCF_MULTI))
	  ALLOCATE(CHI_MULTI(NW,NB,NCF_MULTI))
	END IF
!
! In the following, ALPHA is defined to be positive when the variable decreases with increasing R.
! 
	DO IB=1,NB
	  NS=NS_VEC(IB); IST=IST_VEC(IB); IEND=IEND_VEC(IB)
	  WRITE(6,*)'Doing ETA and CHI shifts for IB, IST, IEND=',IB, IST, IEND
	  IF(ALLOCATED(V_VEC))DEALLOCATE(V_VEC,WRK_VEC,WRK_MAT,V_SCL)
	  ALLOCATE(V_VEC(NS),WRK_VEC(NS))
	  ALLOCATE(WRK_MAT(NS,NCF_MULTI),V_SCL(NS))
	  IF(LOG_V_SHIFT)THEN
	    V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)/dV_BETA(IB)
	    V_SCL=(V_VEC/(V_VEC*dV_BETA(IB)))**OP_SCL_EXPONENT	
	  ELSE
	    V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)-dV_BETA(IB)
	    V_SCL=(V_VEC/(V_VEC+dV_BETA(IB)))**OP_SCL_EXPONENT
	  END IF
!
	  CALL MON_INTERP_FAST(WRK_VEC,NS,IONE,V_VEC,NS,T,ND,V,ND,'T_MULTI')
	  ALPHA=LOG(T(1)/T(5))/LOG(R(5)/R(1)); ALPHA=MAX(2.0_LDP,ALPHA)
	  TEMP_MULTI(IST:IEND,IB)=WRK_VEC(1:NS)
	  IF(IST  .NE. 1) TEMP_MULTI(1:IST-1,IB)  =TEMP_MULTI(IST,IB)*(R_MULTI(1:IST-1)/R(1))**ALPHA
	  ALPHA=LOG(T(ND)/T(ND-5))/LOG(R(ND-5)/R(ND))
	  IF(IEND .NE. NW)TEMP_MULTI(IEND+1:NW,IB)=TEMP_MULTI(IEND,IB)/(R_MULTI(IEND)/R_MULTI(IEND+1:NW))**ALPHA
!
	  CALL MON_INTERP_FAST(WRK_VEC,NS,IONE,V_VEC,NS,ESEC,ND,V,ND,'ESEC_MULTI')
	  ESEC_MULTI(IST:IEND,IB)=WRK_VEC(1:NS)*V_SCL(1:NS)
	  ALPHA=LOG(ESEC(1)/ESEC(5))/LOG(R(5)/R(1)); ALPHA=MAX(2.0_LDP,ALPHA)
	  IF(IST  .NE. 1) ESEC_MULTI(1:IST-1,IB)  =ESEC_MULTI(IST,IB)*(R_MULTI(IST)/R_MULTI(1:IST-1))**ALPHA
	  ALPHA=LOG(ESEC(ND)/ESEC(ND-5))/LOG(R(ND-5)/R(ND))
	  IF(IEND .NE. NW)ESEC_MULTI(IEND+1:NW,IB)=ESEC_MULTI(IEND,IB)/(R_MULTI(IEND)/R_MULTI(IEND+1:NW))**ALPHA
!
	  CALL MON_INTERP_FAST(WRK_VEC,NS,IONE,V_VEC,NS,SIGMA,ND,V,ND,'SIGMA_MULTI')
	  IF(IST  .NE. 1)SIGMA_MULTI(1:IST-1,IB)=SIGMA(1)
	  IF(IEND .NE. NW)SIGMA_MULTI(IEND+1:NW,IB)=SIGMA(ND)
	  SIGMA_MULTI(IST:IEND,IB)=WRK_VEC(1:NS)
	END DO
!
!Plot ESEC as a function of Beta. Usefule to see what the shited data lools like.
!
	CALL DP_CURVE_LAB(ND,V,ESEC,'Orig')
	T1=45.0_LDP/ATAN(1.0_LDP)               !Radians to degrees
	DO IB=1,NB
	  IF(IB .LE. 3 .OR. MOD(IB,5) .EQ. 0)THEN
	    WRITE(TMP_STR,'(F5.1)')BETA_MULTI(IB)*T1; TMP_STR=ADJUSTL(TMP_STR)
	    WRITE(CURVE_LAB,'(I3,A)')IB,')'; CURVE_LAB='('//ADJUSTL(CURVE_LAB)
	    CURVE_LAB=TRIM(TMP_STR)//CURVE_LAB
	    CALL DP_CURVE_LAB(ND_MULTI,VR_MULTI,ESEC_MULTI(1,IB),CURVE_LAB)
	  END IF
	END DO
	CALL GRAMON_PGPLOT('V(km/s)','Esec',' ',' ')
!
	IF(DO_ETA_SHIFT)THEN
	  ALLOCATE(OLD_CHI(ND,NCF_MULTI));  OLD_CHI(:,1:NCF_MULTI)=CHI(:,ML_ST:ML_END)
	  ALLOCATE(OLD_ETA(ND,NCF_MULTI));  OLD_ETA(:,1:NCF_MULTI)=ETA(:,ML_ST:ML_END)
	  DO ML=1,NCF_MULTI
	     OLD_CHI(:,ML)=OLD_CHI(:,ML)*CLUMP_FAC(:)
	     OLD_ETA(:,ML)=OLD_ETA(:,ML)*CLUMP_FAC(:)
	  END DO
	  DO ML=1,NCF_MULTI
	    OLD_CHI(:,ML)=OLD_CHI(:,ML)-ESEC(:)
	  END DO
!
	  DO IB=1,NB
	    NS=NS_VEC(IB); IST=IST_VEC(IB); IEND=IEND_VEC(IB)
	    WRITE(6,*)'Doing ETA and CHI shifts for IB, IST, IEND=',IB, IST, IEND
	    IF(ALLOCATED(V_VEC))DEALLOCATE(V_VEC,WRK_VEC,WRK_MAT,V_SCL)
	    ALLOCATE(V_VEC(NS),WRK_VEC(NS))
	    ALLOCATE(WRK_MAT(NS,NCF_MULTI),V_SCL(NS))
	    IF(LOG_V_SHIFT)THEN
	      V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)/dV_BETA(IB)
	      V_SCL=(V_VEC/(V_VEC*dV_BETA(IB)))**OP_SCL_EXPONENT
	    ELSE
	      V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)-dV_BETA(IB)
	      V_SCL=(V_VEC/(V_VEC+dV_BETA(IB)))**OP_SCL_EXPONENT
	    END IF
!
	    CALL MON_INTERP_FAST(WRK_MAT,NS,NCF_MULTI,V_VEC,NS,OLD_ETA,ND,V,ND,'ETA')
	    DO ML=1,NCF_MULTI
	      ETA_MULTI(IST:IEND,IB,ML) =WRK_MAT(1:NS,ML)*V_SCL(1:NS)
	      ALPHA=LOG(OLD_ETA(ND,ML)/OLD_ETA(ND-5,ML))/LOG(R(ND-5)/R(ND))
	      DO I=IEND+1,NW
	        ETA_MULTI(I,IB,ML)=ETA_MULTI(IEND,IB,ML)*(R_MULTI(IEND)/R_MULTI(I))**ALPHA
	      END DO
	      ALPHA=LOG(OLD_ETA(5,ML)/OLD_ETA(1,ML))/LOG(R(1)/R(5)); ALPHA=MAX(4.0_LDP,ALPHA)
	      DO I=1,IST-1
              ETA_MULTI(I,IB,ML)  =ETA_MULTI(IST,IB,ML)*(R_MULTI(IST)/R_MULTI(I))**ALPHA
	      END DO
	    END DO
!
	    CALL MON_INTERP_FAST(WRK_MAT,NS,NCF_MULTI,V_VEC,NS,OLD_CHI,ND,V,ND,'CHI')
	    DO ML=1,NCF_MULTI
	      CHI_MULTI(IST:IEND,IB,ML) =WRK_MAT(1:NS,ML)*V_SCL(1:NS)
	      T1=OLD_CHI(ND,ML)/OLD_CHI(ND-5,ML)
	      IF(T1 .LE. 0)THEN
	        ALPHA=0
	      ELSE
	        ALPHA=LOG(OLD_CHI(ND,ML)/OLD_CHI(ND-5,ML))/LOG(R(ND-5)/R_MULTI(ND))
	      END IF
	      ALPHA=MAX(0.0_LDP,ALPHA)
	      DO I=IEND+1,NW
	        CHI_MULTI(I,IB,ML)=CHI_MULTI(IEND,IB,ML)*(R_MULTI(IEND)/R_MULTI(I))**ALPHA
	      END DO
	      T1=OLD_CHI(5,ML)/OLD_CHI(1,ML)
	      IF(T1 .LE. 0)THEN
	        ALPHA=2
	      ELSE
	        ALPHA=LOG(OLD_CHI(5,ML)/OLD_CHI(1,ML))/LOG(R(1)/R_MULTI(5))
	      END IF
	      ALPHA=MAX(2.0_LDP,ALPHA)
	      DO I=1,IST-1
              CHI_MULTI(I,IB,ML)=CHI_MULTI(IST,IB,ML)*(R_MULTI(IST)/R_MULTI(I))**ALPHA
	      END DO
	    END DO
	    DO ML=1,NCF_MULTI
	      CHI_MULTI(:,IB,ML)=CHI_MULTI(:,IB,ML)+ESEC_MULTI(:,IB)
	    END DO
!
	  END DO
	  DEALLOCATE(OLD_CHI,OLD_ETA)
	  WRITE(6,'(/,A,2ES14.6,A)')' MIN/MAX ETA',MINVAL(ETA_MULTI),MAXVAL(ETA_MULTI),' (small frequency grid)'
	  WRITE(6,'(A,2ES14.6)')    ' MIN/MAX CHI',MINVAL(CHI_MULTI),MAXVAL(CHI_MULTI)
	END IF
	DEALLOCATE(V_VEC, WRK_VEC, V_SCL, WRK_MAT)
!
	RETURN
	END SUBROUTINE DO_SHIFT
!
	SUBROUTINE WRITE_MULTI
        USE SET_KIND_MODULE
	IMPLICIT NONE
	INTEGER LU,J
!
	WRITE(6,'(/,1X,A,/)')'Opening 3D_DATA for output'
        WRITE(6,*)'       Number of depth points is',SHAPE(R_MULTI)
        WRITE(6,*)'     Shape of ETA [ND,NB,NCF] is',SHAPE(ETA_MULTI)
        WRITE(6,*)'Shape of VR_MULTI [ND,NB,NCF] is',SHAPE(ETA_MULTI)
!
	CALL GET_LU(LU,'In WRITE_MULTI')
	OPEN(UNIT=LU,FILE='3D_DATA',FORM='UNFORMATTED',ACCESS='SEQUENTIAL',ACTION='WRITE',STATUS='UNKNOWN')
          WRITE(LU)ND_MULTI,NB_MULTI,NCF_MULTI
	  WRITE(LU)(R_MULTI,J=1,NB_MULTI)
	  WRITE(LU)BETA_MULTI
	  WRITE(LU)FREQ_MULTI
	  WRITE(LU)TEMP_MULTI
	  WRITE(LU)VR_MULTI
	  WRITE(LU)SIGMA_MULTI
          WRITE(LU)ESEC_MULTI
          WRITE(LU)CHI_MULTI
          WRITE(LU)ETA_MULTI
        CLOSE(LU)
!
	RETURN
	END SUBROUTINE WRITE_MULTI
!
	SUBROUTINE WRITE_JH_POL(ES_J,FINE_FREQ,V,LAM_ST,LAM_END,ES_RES_KMS,NCF,NC,ND,NPHI)
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Input values
!
	INTEGER ND
	INTEGER NC
	INTEGER NCF
	INTEGER NPHI
!
	REAL(KIND=LDP) ES_J(ND,NCF)
	REAL(KIND=LDP) FINE_FREQ(NCF)
	REAL(KIND=LDP) V(1:ND)
!
	REAL(KIND=LDP) ES_RES_KMS
	REAL(KIND=LDP) LAM_ST, MIN_FREQ
	REAL(KIND=LDP) LAM_END, MAX_FREQ
!
! Local variables
!
	REAL(KIND=LDP), ALLOCATABLE :: J_L_ALT(:,:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CRSE_ES_J(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ORIG_J(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_ES_J(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: CRSE_J_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: WRK_MULTI(:)
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: CRSE_FREQ(:)
!
	REAL(KIND=LDP)  PHI_VEC(NPHI)
	REAL(KIND=LDP)  P(NCF+ND)
!
	REAL(KIND=LDP)  DEL_NU
	INTEGER NPNT_CRSE
	INTEGER NF_CRSE
	INTEGER NC_CRSE
	INTEGER NP_CRSE
	INTEGER NP_INS
	INTEGER NR_MAX
	INTEGER, PARAMETER :: IONE=1
!
	REAL(KIND=LDP)  PI
	REAL(KIND=LDP)  T1, T2, T3
	INTEGER ID, IB	
	INTEGER NS, NW, IST, IEND
	INTEGER I, J, K, ML, NP
!
	CHARACTER FORMFEED*2
	CHARACTER TIME*20,FMT*132
	CHARACTER*12 PRODATE
	PARAMETER (PRODATE='13-Nov-2018') !Must be changed after alterations
!
	INTEGER, PARAMETER :: LU   =300
	INTEGER, PARAMETER :: LUOUT=9
        INTEGER, PARAMETER :: LUMOM=10
        INTEGER, PARAMETER :: LUMOD=20
        INTEGER, PARAMETER :: LUOBS=23
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	  WRITE(6,'(/,A)')' Entered WRITE_JH_POL'
	  I=12
          FORMFEED=' '//CHAR(I)
!
! Compute CRSE frequency grid for upon which we store the electron scattering
! source function (/ESEC)
!
	  MIN_FREQ=2.99792E+03_LDP/MAX(LAM_ST,LAM_END)
          MAX_FREQ=2.99792E+03_LDP/MIN(LAM_ST,LAM_END)
          DEL_NU=0.5D0*(MAX_FREQ+MIN_FREQ)*ES_RES_KMS/3.0E+05
          NF_CRSE=(MAX_FREQ-MIN_FREQ)/DEL_NU+1
	  DEL_NU=(MAX_FREQ-MIN_FREQ)/(NF_CRSE-1)
	  ALLOCATE(CRSE_FREQ(NF_CRSE))
	  DO ML=1,NF_CRSE
            CRSE_FREQ(ML)=MAX_FREQ-(ML-1)*DEL_NU
          END DO
          CRSE_FREQ(NF_CRSE)=MIN_FREQ
	  WRITE(6,*)'Set coarse freq. grid, NF_CRSE=',NF_CRSE
!
	  PI=ACOS(-1.0D0)
	  DO I=1,NPHI
            PHI_VEC(I)=(I-1)*PI/(NPHI-1)
          END DO
!
! Two options are present to specify the azimuth angle Beta.
! Create the file LINE_MOM_DATA here
!
	  CALL DATE_TIME(TIME)
          OPEN(UNIT=LUMOD,FILE='LINE_MOM_DATA',STATUS='UNKNOWN')
          WRITE(LUMOD,'(/,'' Model Started on:'',15X,(A))')TIME
          WRITE(LUMOD,'('' Main program last changed on:'',3X,(A))')PRODATE
          WRITE(LUMOD,'('' Main Format Date:'',15X,(A))')'05-Feb-1994'
!
! We will use a blank line in the outputfile to indicate the end of a section.
!
! *****************************
!
! Spatial Co-ordinates.
!
          WRITE(LUMOD,'()')
          FMT='(7X,I8,3X,''!Number of depth points [ND]'')'
          WRITE(LUMOD,FMT)ND_MULTI
          FMT='(7X,I8,3X,''!Number of polar angles [NBETA]'')'
          WRITE(LUMOD,FMT)NB_MULTI
!
! Frequency dependance.
!
          WRITE(LUMOD,'()')
          FMT='(7X,I8,3X,''!Number of CMF frequencies [NF_CRSE]'')'
          WRITE(LUMOD,FMT)NF_CRSE
!
! Co-ordiantes for angular description of radiation field (at a given
! spatial depth).
!
	  NC_CRSE = NC
	  NP_CRSE = NC+ND_MULTI
	  NP_INS  = 0
	  NPNT_CRSE = 0
	  NR_MAX  = ND_MULTI
          FMT='(7X,I8,3X,''!Number of azimuthal angles [NPHI]'')'
          WRITE(LUMOD,FMT)NPHI
          FMT='(7X,I8,3X,''!Number of point source rays [NPNT_CRSE]'')'
          WRITE(LUMOD,FMT)NPNT_CRSE
          FMT='(7X,I8,3X,''!Number of core rays [NC_CRSE]'')'
          WRITE(LUMOD,FMT)NC_CRSE
          FMT='(7X,I8,3X,''!Total number of rays [NP_CRSE]'')'
          WRITE(LUMOD,FMT)NP_CRSE
          FMT='(7X,I8,3X,''!Number of rays inserted [NP_INS]'')'
          WRITE(LUMOD,FMT)NP_INS
          FMT='(7X,I8,3X,''!Number of ray points [NR_MAX]'')'
          WRITE(LUMOD,FMT)NR_MAX
!
! *****************************
!
! Compute coarse impact grid.
!
	  NW=ND_MULTI;  NP=NC+NW
	  ALLOCATE (TA(NW)); TA(1:NW)=R_MULTI(1:NW)
	  CALL IMPAR(P,TA,TA(NW),NC,NW,NP)
!
          CALL RITE_2DJ_ASC_JH(R_MULTI,NW,1,'Radius Coordinates - [R]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(BETA_MULTI,NB_MULTI,1,'Polar Coordinates - [BETA]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(CRSE_FREQ,NF_CRSE,1,'Coarse CMF frequency grid - [FREQ]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(PHI_VEC,NPHI,1,'Azimuthal Coordinates - [PHI]:',LUMOD)
          CALL RITE_2DJ_ASC_JH(P,NP,1,'Impact Parameters - [P]:',LUMOD)
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')FORMFEED
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Boundary Parameters Format Date:'',14X,(A))'
          WRITE(LUMOD,FMT)'06-Dec-1993'
          FMT='(14X,L1,3X,''!Optically thin core?'')'
          WRITE(LUMOD,FMT).FALSE.   !THIN_CORE
          FMT='(14X,L1,3X,''!Schuster boundary condition at core?'')'
          WRITE(LUMOD,FMT).FALSE.   !SCHUS_BC
          FMT='(X,1PE14.6,3X,''!Schuster intensity'')'
          WRITE(LUMOD,FMT)1.0E+00   !I_SCHUS
!
          FMT='(14X,L1,3X,''!Diffusion boundary condition at core?'')'
          WRITE(LUMOD,FMT).TRUE.    !DIFF
          FMT='(X,1PE14.6,3X,''!Diffusion flux'')'
          WRITE(LUMOD,FMT)0.0E+00   !DBB
!
          FMT='(14X,L1,3X,''!Point source?'')'
          WRITE(LUMOD,FMT).FALSE.   !PNT_SRCE
          FMT='(X,1PE14.6,3X,''!Radius of point source'')'
          WRITE(LUMOD,FMT)0.0E+00   !R_PNT_SRCE
          FMT='(X,1PE14.6,3X,''!B Point source'')'
          WRITE(LUMOD,FMT)0.0E+00   !B_PNT_SRCE
!
          FMT='(14X,L1,3X,''!Optically thick at outer boundary?'')'
          WRITE(LUMOD,FMT).FALSE.   !THK
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Angular Parameters Format Date:'',15X,(A))'
          WRITE(LUMOD,FMT)'30-Dec-1993'
          WRITE(LUMOD,'('' Angular Law:'',15X,(A))')'POW_E3' 	!DEN_FN_NAME
          FMT='(7X,I8,3X,''!Number of density parameters'')'
          WRITE(LUMOD,FMT) 3               			!N_FN_PRMS
          DO I=1,3                         			!N_FN_PRMS
            FMT='(X,1PE14.6,3X,''!A'',I1)'
            WRITE(LUMOD,FMT)0.0E+0,I       			!DEN_FN_PRMS(I),I
          END DO
!
          FMT='(X,1PE14.6,3X,''!TH_ANG_EXP'')'
          WRITE(LUMOD,FMT)0.0E+00          			!TH_ANG_EXP
          FMT='(X,1PE14.6,3X,''!ETA_ANG_EXP'')'
          WRITE(LUMOD,FMT)0.0E+00          			!ETA_ANG_EXP
          FMT='(X,1PE14.6,3X,''!ESEC_ANG_EXP'')'
          WRITE(LUMOD,FMT)0.0E+00          			!ESEC_ANG_EXP
          FMT='(X,1PE14.6,3X,''!ESEC_ANG_EXP'')'
          WRITE(LUMOD,FMT)0.0E+00          			!ETAL_ANG_EXP
          FMT='(X,1PE14.6,3X,''!ESEC_ANG_EXP'')'
          WRITE(LUMOD,FMT)0.0E+00     				!CHIL_ANG_EXP
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Opacity Parameters Format Date:'',15X,(A))'
          WRITE(LUMOD,FMT)'06-Dec-1993'
          WRITE(LUMOD,'('' Opacity Law:'',15X,(A))')'File'
          FMT='(X,1PE14.6,3X,''!Constant for thermal Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MCHI(1)
          FMT='(X,1PE14.6,3X,''!Exponent for thermal Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MCHI(2)
          FMT='(X,1PE14.6,3X,''!Constant for thermal Emissivity'')'
          WRITE(LUMOD,FMT)0.0E+00   !META(1)
!
          FMT='(X,1PE14.6,3X,''!Exponent for thermal Emissivity'')'
          WRITE(LUMOD,FMT)0.0E+00   !META(2)
          FMT='(X,1PE14.6,3X,''!Constant for electron scattering'//
	1       ' Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MESEC(1)
          FMT='(X,1PE14.6,3X,''!Exponent for electron scattering'//
	1       ' Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MESEC(2)
          FMT='(X,1PE14.6,3X,''!Constant for line Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MCHIL(1)
          FMT='(X,1PE14.6,3X,''!Exponent for line Opacity'')'
          WRITE(LUMOD,FMT)0.0E+00   !MCHIL(2)
          FMT='(X,1PE14.6,3X,''!Constant for line Emissivity'')'
          WRITE(LUMOD,FMT)0.0E+00   !METAL(1)
          FMT='(X,1PE14.6,3X,''!Exponent for line Emissivity'')'
          WRITE(LUMOD,FMT)0.0E+00   !METAL(2)
!
! *****************************
!
          !CALL RITE_2DJ_ASC_JH(CHI,ND,1,'Thermal Opacity.',LUMOD)
          !CALL RITE_2DJ_ASC_JH(ETA,ND,1,'Thermal Emissivity.',LUMOD)
!
	  TA(1:NW) = 0.0E+00
          CALL RITE_2DJ_ASC_JH(TA,NW,1,'Thermal Opacity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(TA,NW,1,'Thermal Emissivity.',LUMOD)
          TA(1:NW)=ESEC_MULTI(1:NW,1)
	  CALL RITE_2DJ_ASC_JH(TA,NW,1,'Electron Scattering Opacity.',LUMOD)
!
	  TA(1:NW) = 0.0E+00
          CALL RITE_2DJ_ASC_JH(TA,NW,1,'Line Opacity.',LUMOD)
          CALL RITE_2DJ_ASC_JH(TA,NW,1,'Line Emissivity.',LUMOD)
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')FORMFEED
!
! *****************************
!
          WRITE(LUMOD,'()')
          FMT='('' Structure Parameters Format Date:'',15X,(A))'
          WRITE(LUMOD,FMT)'06-Dec-1993'
!
          TA(1:NW)=TEMP_MULTI(1:NW,1);   CALL RITE_2DJ_ASC_JH(TA,NW,1,'Temperature (10^4K)',LUMOD)
          TA(1:NW)=VR_MULTI(1:NW,1);     CALL RITE_2DJ_ASC_JH(TA,NW,1,'Velocity (km/s)',LUMOD)
          TA(1:NW)=SIGMA_MULTI(1:NW,1);  CALL RITE_2DJ_ASC_JH(TA,NW,1,'Velocity gradient [dlnV/Dlnr-1]',LUMOD)
!
! Place J on a coarser frequency grid. We can do this if dominated by electron scattering.
!
	  WRITE(6,'(1X,A,2ES16.6,A)')  'Min/Max coarse NU:',MINVAL(CRSE_FREQ),MAXVAL(CRSE_FREQ)
	  WRITE(6,'(1X,A,2ES16.6,A)')  '  Min/Max fine NU:',MINVAL(FINE_FREQ),MAXVAL(FINE_FREQ)
	  WRITE(6,'(1X,A,2ES16.6,A,/)')'     Min/Max ES_J:',MINVAL(ES_J),MAXVAL(ES_J),' (full frequency range)'
	  ALLOCATE(TMP_ES_J(NF_CRSE,ND),CRSE_ES_J(ND,NF_CRSE),ORIG_J(NCF,ND))
	  ORIG_J=TRANSPOSE(ES_J)
	  WRITE(6,*)'Done transpose of J.'
          CALL MON_INTERP_FAST_V2(TMP_ES_J,NF_CRSE,ND,CRSE_FREQ,NF_CRSE,ORIG_J,NCF,FINE_FREQ,NCF,'J_ND')
	  WRITE(6,*)'Done interpolation of J onto coarse frequency grid.'
	  CRSE_ES_J=TRANSPOSE(TMP_ES_J)
	  DEALLOCATE(TMP_ES_J,ORIG_J)
	  WRITE(6,'(/,1X,A,2ES14.6,A)')'Min/max CRSE_ES_J is',
	1              MINVAL(CRSE_ES_J),MAXVAL(CRSE_ES_J),' (before regridding in depth)'
!
! Now need to put J on the fine spatial grid.
!
	  ALLOCATE(J_L_ALT(NF_CRSE,NB_MULTI,ND_MULTI))
	  ALLOCATE(V_VEC(NW),WRK_VEC(NW),WRK_MAT(NF_CRSE,NW),CRSE_J_VEC(ND))
!
	  DO IB=1,NB_MULTI
	    NS=NS_VEC(IB); IST=IST_VEC(IB); IEND=IEND_VEC(IB)
            IF(LOG_V_SHIFT)THEN
	      V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)/dV_BETA(IB)
	    ELSE
              V_VEC(1:NS)=VR_MULTI(IST:IEND,IB)-dV_BETA(IB)
	    END IF
	    DO ML=1,NF_CRSE
	      CRSE_J_VEC(1:ND)=CRSE_ES_J(1:ND,ML)
              CALL MON_INTERP_FAST_V2(WRK_VEC(IST),NS,IONE,V_VEC,NS,CRSE_J_VEC,ND,V,ND,'J_MULTI')
	      IF(IST  .NE. 1)WRK_VEC(1:IST-1)=CRSE_J_VEC(IST)
              IF(IEND .NE. NW)WRK_VEC(IEND+1:NW)=CRSE_J_VEC(ND)
	      WRK_MAT(ML,1:NW)=WRK_VEC(1:NW)
	    END DO
	    DO ID=1,NW
	      J_L_ALT(1:NF_CRSE,IB,ID)=WRK_MAT(1:NF_CRSE,ID)
	    END DO
	  END DO
	  DEALLOCATE(V_VEC,WRK_VEC,WRK_MAT,CRSE_J_VEC)
!
	  WRITE(6,'(1X,A,2ES14.6,A)')'  Min/max J_L_ALT is',
	1       MINVAL(J_L_ALT),MAXVAL(J_L_ALT),' (after regridding J in depth)'
!
! J=J_L + J_R
!
	  J_L_ALT=0.5_LDP*J_L_ALT
          WRITE(LUMOD,'()')
          WRITE(LUMOD,'(A)')' Moments of Radiation field:'
          FMT='(5X,I8,5X,''!Number of radiation moments'')'
          WRITE(LUMOD,FMT)5
          I=NB_MULTI*NF_CRSE       !Temporary Variable
          CALL RITE_2DJ_ASC_JH(J_L_ALT,I,NW,'J_L moment.',LUMOD)
	  WRITE(6,'(/,1X,A)')'Written J_L moment'; FLUSH(UNIT=6)
          CALL RITE_2DJ_ASC_JH(J_L_ALT,I,NW,'J_R moment.',LUMOD)
	  WRITE(6,*)'Written J_R moment'; FLUSH(UNIT=6)
!
! We now use J_L_ALT asstorage for the other variables.
!
	  J_L_ALT=J_L_ALT/3.0_LDP
          I=NB_MULTI*NF_CRSE       !Temporary Variable
          CALL RITE_2DJ_ASC_JH(J_L_ALT,I,NW,'K20_L moment.',LUMOD)
	  WRITE(6,'(1X,A)')'Written K20 moment'; FLUSH(UNIT=6)
!
	  J_L_ALT=0.0_LDP
	  CALL RITE_2DJ_ASC_JH(J_L_ALT,I,NW,'GT moment.',LUMOD)
          CALL RITE_2DJ_ASC_JH(J_L_ALT,I,NW,'KT moment.',LUMOD)
	  WRITE(6,*)'Written GT and KT  moments'; FLUSH(UNIT=6)
	  DEALLOCATE(J_L_ALT)
!
	  RETURN
!
	END SUBROUTINE WRITE_JH_POL
	END MODULE SN_SHIFT_MODULE
