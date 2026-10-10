!
! Calling routine to compute the Observer's fluxes. The fluxes are computed
! using an OBSERVER'S FRAME formulation. J, ETA, and CHI, computed in the
! comoving frame, must be supplied.
!
	PROGRAM OBS_FRAME_V2
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! Altered 05-Oct-2026 : Use ABS(V)/MAXVAL(ABS(V)) for VTURB (Devin).
! Altered 30-Apr-2020 : Now include RD_RAY option
!                         New routine implemented to READ CHI_DATA etc.
!                         Input file changed to OBS_FRAME_INPUT
! Altered 20-Aug-2019 : Updated RD_RVTJ to V4 routines.
!                         TGREY,dE_RAD_DECAY,PLANCK_MEAN, and FORMAT_DATE added.
! Altered 13-Apr-2017 --- Updated one call to WRITE_DIRECT_INFO_V3 (RJ call, was orig.)
! Altered 06-Jan-1999 --- Now call OBS_FRAME_SUB_V2.
!                         TAU_MAX, ES_DTAU, and INT_METHOD included.
! Altered 03-Jan-1999 --- MAX_DEL_V_RES_ZONE was being incorrectly computed.
! Altered 11-Feb-1998 --- RJ_VEC was being set outside valid range.
! Created 9-Dec-1998 :
!
	INTEGER ND,NP,NC

	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: IFOUR=4
	INTEGER, PARAMETER :: ISIX=6
	INTEGER, PARAMETER :: ITEN=10
	INTEGER, PARAMETER :: T_OUT=6		!Terminal IO
!
	INTEGER, PARAMETER :: LUIN=8
	INTEGER, PARAMETER :: LUMOD=9
	INTEGER, PARAMETER :: LU_FLUX=20
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	REAL(KIND=LDP), ALLOCATABLE :: R(:)
	REAL(KIND=LDP), ALLOCATABLE :: V(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA(:)
	REAL(KIND=LDP), ALLOCATABLE :: T(:)
	REAL(KIND=LDP), ALLOCATABLE :: ED(:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC(:)
!
! These vectors are neeed so that we can use the standard RVTJ read routine.
! Only CLUMP_FAC is utilized.
!
	REAL(KIND=LDP), ALLOCATABLE :: TGREY(:)
	REAL(KIND=LDP), ALLOCATABLE :: dE_RAD_DECAY(:)
	REAL(KIND=LDP), ALLOCATABLE :: ROSS_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: FLUX_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: PLANCK_MEAN(:)
	REAL(KIND=LDP), ALLOCATABLE :: POP_ATOM(:)
	REAL(KIND=LDP), ALLOCATABLE :: MASS_DENSITY(:)
	REAL(KIND=LDP), ALLOCATABLE :: CLUMP_FAC(:)
	REAL(KIND=LDP), ALLOCATABLE :: POPION(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: VTURB(:)
	REAL(KIND=LDP), ALLOCATABLE :: MAX_DEL_V_RES_ZONE(:)
!
! P is the impact parameter. PQW is R^2 p dp
!
	REAL(KIND=LDP), ALLOCATABLE :: P(:)
	REAL(KIND=LDP), ALLOCATABLE :: MU_AT_RMAX(:)
	REAL(KIND=LDP), ALLOCATABLE :: HQW_AT_RMAX(:)
!
	REAL(KIND=LDP) RMDOT
	REAL(KIND=LDP) RLUM
	REAL(KIND=LDP) ABUND_HYD
	CHARACTER(LEN=20) FILE_DATE
	CHARACTER(LEN=20) TIME
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER(LEN=10) NAME_CONVENTION
	CHARACTER(LEN=80) STRING
!
! 
	INTEGER NLF
	REAL(KIND=LDP), ALLOCATABLE :: FREQ_CMF(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_CMF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_CMF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RAY_CMF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_RJ_CMF(:,:)
!
	REAL(KIND=LDP), ALLOCATABLE :: FREQ_RJ(:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_CMF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: RJ_VEC(:)
	REAL(KIND=LDP), ALLOCATABLE :: NEW_RJ_VEC(:)
!
	INTEGER NOS
	REAL(KIND=LDP), ALLOCATABLE :: OBS_FREQ(:)
	REAL(KIND=LDP), ALLOCATABLE :: OBS_FLUX(:)
! 
!
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
	INTEGER ERROR_LU,LUER
	REAL(KIND=LDP) SPEED_OF_LIGHT,C_KMS
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
!
	REAL(KIND=LDP) TAU_MAX
	REAL(KIND=LDP) ES_DTAU
	REAL(KIND=LDP) FRAC_DOP
	REAL(KIND=LDP) DEL_V_OBS
	REAL(KIND=LDP) VTURB_MIN
	REAL(KIND=LDP) VTURB_MAX
	REAL(KIND=LDP) MIN_FREQ
	REAL(KIND=LDP) MAX_FREQ
	REAL(KIND=LDP) MIN_WAVE
	REAL(KIND=LDP) MAX_WAVE
	REAL(KIND=LDP) T1
	REAL(KIND=LDP) OBS_TAU_MAX
	REAL(KIND=LDP) OBS_ES_DTAU
	REAL(KIND=LDP) TAU_REF
!
	CHARACTER(LEN=10) OBS_INT_METHOD
!
	INTEGER NCF
	INTEGER NCF_RJ
	INTEGER ND_RD
	INTEGER N_INS_OBS
	INTEGER INIT_REC
	INTEGER IOS
	INTEGER I,K,LS,ML
	INTEGER CONT_REC
	INTEGER REC_LENGTH
!
	LOGICAL RD_RAYLEIGH_DATA
	LOGICAL INTERP_RJ_NEC
	LOGICAL CONVOLVE_J
	LOGICAL DO_CLUMP
	LOGICAL WRITE_IP
	LOGICAL WRITE_RTAU
	LOGICAL WRITE_dFR
	LOGICAL DO_REL_IN_OBSFRAME
	LOGICAL TMP_LOG
! 
!
! Set constants.
!
	CHIBF=2.815E-06_LDP
	CHIFF=3.69E-29_LDP
	HDKT=4.7994145_LDP
	TWOHCSQ=0.0147452575_LDP
	OPLIN=2.6540081E+08_LDP
	EMLIN=5.27296E-03_LDP
	OPLIN=2.6540081E+08_LDP
	EMLIN=5.27296E-03_LDP
!
	CONT_REC=3
	LUER=ERROR_LU()
	C_KMS=1.0E-05_LDP*SPEED_OF_LIGHT()
!
	CALL GEN_ASCI_OPEN(LUMOD,'OUT_FRAME','UNKNOWN',' ','WRITE',IZERO,IOS)
	CALL GEN_ASCI_OPEN(LUIN,'OBS_FRAME_INPUT','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening OBS_FRAME_INPUT in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUMOD)
	CLOSE(LUIN)
!
	CALL RD_STORE_DBLE(MIN_WAVE,'MIN_WAVE',L_TRUE,'Minimum observers wavelenghth (Ang)')
	CALL RD_STORE_DBLE(MAX_WAVE,'MAX_WAVE',L_TRUE,'Maximum observers wavelenghth (Ang)')
	CALL RD_STORE_DBLE(DEL_V_OBS,'DEL_V_OBS',L_TRUE,
	1       'Spacing (km/s) for observers frame frequencies ')
!
! NB: VTURB_MIN and VTURB_MAX do NOT directly effect the spectrum. They
! only effect the accuracy of its computation.
!
	CALL RD_STORE_DBLE(VTURB_MIN,'VTURB_MIN',L_TRUE,
	1       'Minimum turbulent velocity (km/s)')
	CALL RD_STORE_DBLE(VTURB_MAX,'VTURB_MAX',L_TRUE,
	1       'Maximum turbulent velocity (units of Vinf)')
	CALL RD_STORE_DBLE(FRAC_DOP,'FRAC_DOP',L_TRUE,
	1       'Fractional spacing (Doppler widths) across resonance one')
!
	CALL RD_STORE_LOG(CONVOLVE_J,'CONV_ES_J',L_TRUE,
	1       'Convolve J with e.s. redistribution function?')
!
	CALL RD_STORE_DBLE(TAU_MAX,'TAU_MAX',L_TRUE,
	1       'Maximum TAU along ray at which integration ceases.')
	CALL RD_STORE_DBLE(ES_DTAU,'ES_DTAU',L_TRUE,
	1       'Maximum DTAU on Electron scattering optical depth scale.')
	CALL RD_STORE_NCHAR(OBS_INT_METHOD,'INT_METH',ITEN,L_TRUE,
	1            'Integration method for computing I along ray')
!
	CALL RD_STORE_DBLE(OBS_TAU_MAX,'TAU_MAX',L_TRUE,
	1    'Optical depth at which observers frame integration is terminated')
	CALL RD_STORE_DBLE(OBS_ES_DTAU,'ES_DTAU',L_TRUE,
	1      'Maximum increments in e.s. optical depth scale')
	CALL RD_STORE_INT(N_INS_OBS,'N_INS_OBS',L_TRUE,
	1          'Mininum number of points to be inserted in '//
	1          ' observers frame grid (>= 0)')
!
	CALL RD_STORE_LOG(WRITE_IP,'WR_IP',L_TRUE,
	1        'Output I as a functio of p and frequency?')
	WRITE_RTAU=.FALSE.
	CALL RD_STORE_LOG(WRITE_RTAU,'WR_RTAU',L_FALSE,
	1        'Output R(Tau=Tau_ref) as a function of p and frequency?')
	WRITE_dFR=.FALSE.
	CALL RD_STORE_LOG(WRITE_dFR,'WR_dFR',L_FALSE,
	1        'Output dFR as a functioin of R and frequency?')
	CALL RD_STORE_DBLE(TAU_REF,'TAU_REF',WRITE_RTAU,
	1        'Reference tau for WR_TAU')
	CALL RD_STORE_LOG(RD_RAYLEIGH_DATA,'RD_RAY',L_TRUE,'Reference tau for WR_TAU')
!
	DO_REL_IN_OBSFRAME=.TRUE.	
	CALL RD_STORE_LOG(DO_REL_IN_OBSFRAME,'DO_RELO',L_FALSE,
	1        'Use all relativistic terms in Observer''s frame calculation.')

!
	CALL CLEAN_RD_STORE()
!
!
! Read in atmosphere parameters.
!
	CALL RD_RVTJ_PARAMS_V4(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1                             ND,NC,NP,FORMAT_DATE,'RVTJ',LUIN)
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
	ALLOCATE (T(ND))
	ALLOCATE (ED(ND))
	ALLOCATE (TGREY(ND))
	ALLOCATE (dE_RAD_DECAY(ND))
	ALLOCATE (ROSS_MEAN(ND))
	ALLOCATE (FLUX_MEAN(ND))
	ALLOCATE (PLANCK_MEAN(ND))
	ALLOCATE (POP_ATOM(ND))
	ALLOCATE (MASS_DENSITY(ND))
	ALLOCATE (POPION(ND))
	ALLOCATE (CLUMP_FAC(ND))
	CALL RD_RVTJ_VEC_V4(R,V,SIGMA,ED,T,TGREY,dE_RAD_DECAY,
	1       ROSS_MEAN,FLUX_MEAN,PLANCK_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,FORMAT_DATE,ND,LUIN)
	CLOSE(LUIN)
!
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,'ETA_DATA',LUIN,IOS)
	IF(IOS .NE. 0)STOP
	  OPEN(UNIT=LUIN,FILE='ETA_DATA',STATUS='OLD',RECL=REC_LENGTH,
	1        ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	  READ(LUIN,REC=CONT_REC)INIT_REC,NCF,ND_RD
	CLOSE(LUIN)
!
	ALLOCATE (FREQ_CMF(NCF))
	ALLOCATE (ETA_CMF(ND,NCF))
	ALLOCATE (CHI_CMF(ND,NCF))
	ALLOCATE (RAY_CMF(ND,NCF))
	CALL READ_CHI_DATA_FILE('ETA_DATA',ETA_CMF,FREQ_CMF,L_TRUE,NCF,ND,LUIN,CONT_REC)
	CALL READ_CHI_DATA_FILE('CHI_DATA',CHI_CMF,FREQ_CMF,L_FALSE,NCF,ND,LUIN,CONT_REC)
	IF(RD_RAYLEIGH_DATA)THEN
	  CALL READ_CHI_DATA_FILE('RAY_DATA',RAY_CMF,FREQ_CMF,L_FALSE,NCF,ND,LUIN,CONT_REC)
	ELSE
	  RAY_CMF=0.0_LDP
	END IF
!
! 
!
! ***************************************************************************
!
! Read in mean intensity as a function of depth and frequency. If requested,
! J will be convolved with the electron scattering redistribution function.
! If the frequency grid for J differs from that for ETA and CHI, it is
! interpolated onto the ETA (and CHI) frequency grid.
!
! ***************************************************************************
!
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,'RJ_DATA',LUIN,IOS)
	OPEN(UNIT=LUIN,FILE='RJ_DATA',STATUS='OLD',RECL=REC_LENGTH,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	READ(LUIN,REC=CONT_REC)INIT_REC,NCF_RJ,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in OBS_FRAME'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	ALLOCATE (RJ_CMF(ND,NCF_RJ))
	ALLOCATE (FREQ_RJ(NCF_RJ))
	DO ML=1,NCF_RJ
	  READ(LUIN,REC=INIT_REC+ML-1)(RJ_CMF(K,ML),K=1,ND),FREQ_RJ(ML)
	END DO
!
! Detemine if RJ_CMF is on same grid as ETA and CHI. If not put on same grid
! via interpolation. We then combine RJ with ETA.
!
! At present the ETA_CMF frequency range must be included with RJ_CMF frequency range.
! In practice this could be changed since we later restrinct the frequency range
! for which we are computing the spectrum.
!
	IF(NCF .EQ. NCF_RJ)THEN
	  INTERP_RJ_NEC=.FALSE.
	  DO ML=1,NCF
	    IF(FREQ_CMF(ML) .NE. FREQ_RJ(ML))INTERP_RJ_NEC=.TRUE.
	  END DO
	ELSE
	  INTERP_RJ_NEC=.TRUE.
	END IF
	IF(INTERP_RJ_NEC)THEN
	  ALLOCATE (RJ_VEC(NCF_RJ))
	  ALLOCATE (NEW_RJ_VEC(NCF))
	  ALLOCATE (NEW_RJ_CMF(ND,NCF))
	END IF
!
	IF(RD_RAYLEIGH_DATA)THEN
	  IF(.NOT. CONVOLVE_J)THEN
	    WRITE(6,*)'As Rayleigh scattering is coherent, you should read'
	    WRITE(6,*)'in data that is not already convolved with the e.s. '
	    WRITE(6,*)'redistribution function'
	    STOP
	  END IF
	  IF(INTERP_RJ_NEC)THEN
	    DO K=1,ND
	      RJ_VEC(1:NCF_RJ)=RJ_CMF(K,1:NCF_RJ)
	      CALL MON_INTERP(NEW_RJ_VEC,NCF,IONE,FREQ_CMF,NCF,
	1                             RJ_VEC,NCF_RJ,FREQ_RJ,NCF_RJ)
	      NEW_RJ_CMF(K,1:NCF)=NEW_RJ_VEC(1:NCF)
	    END DO
	    DO ML=1,NCF
	       DO I=1,ND
	         ETA_CMF(I,ML)=ETA_CMF(I,ML)+NEW_RJ_CMF(I,ML)*RAY_CMF(I,ML)
	      END DO
	    END DO
	  ELSE
	    DO ML=1,NCF
	      DO I=1,ND
	        ETA_CMF(I,ML)=ETA_CMF(I,ML)+NEW_RJ_CMF(I,ML)*RAY_CMF(I,ML)
	      END DO
	    END DO
	  END IF
	END IF
!
! If requested, convolve J with the electron-scattering redistribution
! funtion. K is used for LUIN, and LUOUT but is not accessed.
!
	IF(CONVOLVE_J)THEN
	  I=ND*NCF_RJ
	  CALL COMP_J_CONV_V2(RJ_CMF,I,FREQ_RJ,T,ND,NCF_RJ,
	1         K,'J PASSED VIA CALL',CONT_REC,L_FALSE,L_FALSE,
	1         K,'RETURN J VIA CALL')
	END IF
!
! We do not need to multiply ESEC (and RAY_CMF) by the climping factor since
! they are used to update ETA, and ETA is multiplied.
!
	ALLOCATE (ESEC(ND))
	ESEC(1:ND)=6.65E-15_LDP*ED(1:ND)
	IF(INTERP_RJ_NEC)THEN
	  DO K=1,ND
	    RJ_VEC(1:NCF_RJ)=RJ_CMF(K,1:NCF_RJ)
	    CALL MON_INTERP(NEW_RJ_VEC,NCF,IONE,FREQ_CMF,NCF,
	1                             RJ_VEC,NCF_RJ,FREQ_RJ,NCF_RJ)
	    NEW_RJ_CMF(K,1:NCF)=NEW_RJ_VEC(1:NCF)
	  END DO
	  DO ML=1,NCF
	    DO I=1,ND
	      ETA_CMF(I,ML)=ETA_CMF(I,ML)+NEW_RJ_CMF(I,ML)*ESEC(I)
	    END DO
	  END DO
	ELSE
	  DO ML=1,NCF
	    DO I=1,ND
	      ETA_CMF(I,ML)=ETA_CMF(I,ML)+RJ_CMF(I,ML)*ESEC(I)
	    END DO
	  END DO
	END IF
!
	DEALLOCATE (RJ_CMF)
	IF(INTERP_RJ_NEC)THEN
	  DEALLOCATE (NEW_RJ_CMF)
	  DEALLOCATE (RJ_VEC)
	  DEALLOCATE (NEW_RJ_VEC)
	END IF
!
! 
! ***************************************************************************
! ***************************************************************************
!
! Determine the turbulent velocity as a funtion of depth. This is only
! used to define the points inserted along each ray via MAX_DEL_V_RES_ZONE.
!
	ALLOCATE (VTURB(ND))
	ALLOCATE (MAX_DEL_V_RES_ZONE(ND))
!
! Bug fix, not flag-gated: this is a speed law, so normalize by the
! full velocity magnitude range.  The old V(I)/V(1) form can divide by
! zero or produce a negative turbulent width for a non-monotonic field.
!
	T1=MAXVAL(ABS(V(1:ND)))
	IF(T1 .GT. 1.0E-12_LDP)THEN
	  DO I=1,ND
	    VTURB(I)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*
	1        MIN(1.0_LDP,ABS(V(I))/T1)
	  END DO
	ELSE
	  VTURB(1:ND)=VTURB_MIN
	END IF
!
	DO I=1,ND
	  MAX_DEL_V_RES_ZONE(I)=VTURB(I)*FRAC_DOP
	END DO
!
! If CLUMPING is important, need to correct emissivities and opacities.
! NB: ETA and CHI read from the files do NOT contain the clumping factor.
!
	DO_CLUMP=.FALSE.
	DO I=1,ND
	  IF( ABS(CLUMP_FAC(I)-1.0_LDP) .GT. 1.0E-06_LDP)DO_CLUMP=.TRUE.
	END DO
!
	IF(DO_CLUMP)THEN
	  DO ML=1,NCF
	    ETA_CMF(:,ML)=ETA_CMF(:,ML)*CLUMP_FAC(:)
	    CHI_CMF(:,ML)=CHI_CMF(:,ML)*CLUMP_FAC(:)
	  END DO
	END IF
!
! 
! ***************************************************************************
!
! Define the impact parameters and appropriate quarature weights for obtaing
! the observed flux. At present the non-core rays should be defined by the
! radius grid.
!
! ***************************************************************************
!
	NC=5+ND/4
	NP=ND+NC
	ALLOCATE (P(NP))
	ALLOCATE (MU_AT_RMAX(NP))
	ALLOCATE (HQW_AT_RMAX(NP))
!
! Compute impact parameter valiues.
!
	CALL IMPAR(P,R,R(ND),NC,ND,NP)
	T1=R(ND)/(NC-1)
	DO I=2,NC-1
	  P(I)=P(I-1)+T1
	END DO
	P(NC)=P(NC+1)-0.001_LDP*(P(NC+1)-P(NC))
!
! By definition, p * dp equals R**2 * mu * dmu. Integration over mu is
! more stable, and is to be preferred. To get better accuracy with the
! integration, NORDULUND weights will be used (Changed 11-Dec-1986).
!
	T1=R(1)*R(1)
	DO LS=1,NP
	  MU_AT_RMAX(LS)=SQRT(T1-P(LS)*P(LS))/R(1)
	END DO
	CALL HWEIGHT(MU_AT_RMAX,HQW_AT_RMAX,NP)
!
! ***************************************************************************
! ***************************************************************************
!
! Define the observer's frame frequency grid. This grid is in units of
! 10^15 Hz and must be ordered from highest to lowest (consistent with
! the comoving frame frequencies).
!
! ***************************************************************************
! ***************************************************************************
!
! 0.01D0=(1.0D+3/1.0D-10)/1.0D+15
!
	MAX_FREQ=0.01_LDP*C_KMS/MIN_WAVE
	MIN_FREQ=0.01_LDP*C_KMS/MAX_WAVE
	NOS=LOG(MAX_FREQ/MIN_FREQ)/LOG(1.0_LDP+DEL_V_OBS/C_KMS)-1
!
	ALLOCATE (OBS_FLUX(NOS))
	ALLOCATE (OBS_FREQ(NOS))
!
	OBS_FREQ(1)=MAX_FREQ
	DO I=2,NOS
	  OBS_FREQ(I)=MAX_FREQ/(1.0_LDP+DEL_V_OBS/C_KMS)**(I-1)
	END DO
!
	WRITE(LUER,*)' '
	WRITE(LUER,'(A,I3)')'                 Number of depth points is: ',ND
	WRITE(LUER,'(A,I3)')'                    Number of core rays is: ',NC
	WRITE(LUER,'(A,I3)')'            Number of impact parameters is: ',NP
	WRITE(LUER,'(A,I6)')'         Number of CMF frequency points is: ',NCF
	WRITE(LUER,'(A,I6)')'  Number of observer''s frequency points is: ',NOS
	WRITE(LUER,*)' '
!
! ***************************************************************************
! ***************************************************************************
!
! Compute the observer's frame fluxes. The fluxes are returned in Janskies.
!
	TMP_LOG=.FALSE.
	IF(NP .LT. NC+4)TMP_LOG=.TRUE.                  !Indicates plane-parallel
	CALL OBS_FRAME_SUB_V9(ETA_CMF,CHI_CMF,FREQ_CMF,
	1            R,V,T,ED,ND,NCF,
	1            P,MU_AT_RMAX,HQW_AT_RMAX,NC,NP,
	1            OBS_FREQ,OBS_FLUX,NOS,
	1            MAX_DEL_V_RES_ZONE,OBS_TAU_MAX,OBS_ES_DTAU,
	1            N_INS_OBS,OBS_INT_METHOD,
	1            TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_dFR,
	1            DO_REL_IN_OBSFRAME,TMP_LOG)
	CALL TUNE(3,' ')
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFRAME','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(STRING,'(I10)')NOS
	  STRING=ADJUSTL(STRING)
	  STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	  CALL WRITV_V2(OBS_FREQ,NOS,ISIX,STRING,LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,NOS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
!
	STOP
	END
