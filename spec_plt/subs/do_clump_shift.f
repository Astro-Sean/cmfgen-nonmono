!
! Subroutine designed to offset shells in R as a function of Beta. The
! shells are (generally) shifted radomly in R -- the distance moved being
! related to the shell width. Two main output files (plus some diagnoistic 
! files) are created:
!
!    3D_DATA -- Large file contaning R, V, SIGMA, T, ESEC, CHI, ETA
!                   R is the same for all beta
!                   VR, SIGMA, T, and ESES are functions of R, BETA
!                   ETA and CHI -- function of R, BETA and frequency
!                   ETA -- no scattering cintribution
!
!    LINE_MOM_DATA -- Primarily used for the 5 polarization moments.
!                     J_L and J_R are assumed to be 0.5J
!                     K= J/3. Other moments are zero. Cannot be used for
!                       polarization calculations.
! 
! Requires RVTJ, ETA_DATA, CHI_DATA and ES_J_CONV data if both output
! files are being output. 3D_DATA and LINE_MOM_DATA do not beed to be
! created.
! 
        MODULE CLUMP_SHIFT_MODULE
        USE SET_KIND_MODULE
!
! Crated 01-May-2025   -   Based on routines developed by Brin Flores.
!                          Cleaned and simplified.
! 
        INTEGER ND_MULTI
        INTEGER NB_MULTI
        INTEGER NCF_MULTI
        INTEGER NCF_CRSE
        INTEGER NPHI_MULTI
!       
        REAL(KIND=LDP), ALLOCATABLE :: R_MULTI(:)
!       
        REAL(KIND=LDP), ALLOCATABLE :: BETA_MULTI(:)
        REAL(KIND=LDP), ALLOCATABLE :: FREQ_MULTI(:)
        REAL(KIND=LDP), ALLOCATABLE :: FREQ_CRSE(:)
!         
        REAL(KIND=LDP), ALLOCATABLE :: ETA_MULTI(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: CHI_MULTI(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: J_MULTI(:,:,:)
        REAL(KIND=LDP), ALLOCATABLE :: ESEC_MULTI(:,:)
!
        REAL(KIND=LDP), ALLOCATABLE :: VR_MULTI(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: SIGMA_MULTI(:,:)
        REAL(KIND=LDP), ALLOCATABLE :: TEMP_MULTI(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: R_SHIFT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SHELL_SHIFT(:)
	REAL(KIND=LDP), ALLOCATABLE :: HALF_SH_LEN(:)
	INTEGER, ALLOCATABLE :: SHELL_LOC(:)
!
	INTEGER NCLUMPS
	INTEGER D_SH_INDX

	END MODULE CLUMP_SHIFT_MODULE

	SUBROUTINE DO_CLUMP_SHIFT (R,T,V,SIGMA,ESEC,MASS_DENSITY,NC,ND,
	1           ETA, CHI, NU, NCF, RJ_ES, R_ES, ND_ES)
	USE SET_KIND_MODULE
	USE CLUMP_SHIFT_MODULE
	USE MOD_USR_LIST_OPT
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
	INTEGER, INTENT(IN) :: NC
	INTEGER, INTENT(IN) :: ND
	REAL(KIND=LDP), INTENT(IN) :: R(ND)
	REAL(KIND=LDP), INTENT(IN) :: T(ND)
	REAL(KIND=LDP), INTENT(IN) :: V(ND)
	REAL(KIND=LDP), INTENT(IN) :: SIGMA(ND)
	REAL(KIND=LDP), INTENT(IN) :: ESEC(ND)
	REAL(KIND=LDP), INTENT(IN) :: MASS_DENSITY(ND)
!
	INTEGER, INTENT(IN) :: NCF
	REAL(KIND=LDP), INTENT(IN) :: ETA(ND,NCF)
	REAL(KIND=LDP), INTENT(IN) :: CHI(ND,NCF)
	REAL(KIND=LDP), INTENT(IN) :: NU(NCF)
!
	INTEGER, INTENT(IN) :: ND_ES
	REAL(KIND=LDP), INTENT(IN) :: RJ_ES(ND_ES,NCF)
	REAL(KIND=LDP), INTENT(IN) :: R_ES(ND_ES)
!
! Local data
!
	REAL(KIND=LDP), ALLOCATABLE :: YV(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: J_CRSE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_MAT(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: TMP_VEC(:)
!
	REAL(KIND=LDP) LAM_START
	REAL(KIND=LDP) LAM_END
	REAL(KIND=LDP) HALF_SH_LEN_SCL_FAC
	REAL(KIND=LDP) DIV_BY_DEN_SQ
	REAL(KIND=LDP) BETA_PARAMS(10)
!
	REAL(KIND=LDP) MIN_FREQ,MAX_FREQ,DEL_NU
	REAL(KIND=LDP) T1,T2,T3
	REAL(KIND=LDP) ES_RES_KMS
!
	INTEGER N_BETA_PARAMS
	INTEGER ADD_N_PTS
	INTEGER INIT_SEED
	INTEGER IERROR
!
	LOGICAL TOP_BOT_SYM
	LOGICAL RAND_IN_BETA
	LOGICAL DIV_BY_SQ_DEN_SM
	LOGICAL INSERT_PTS_IN_DIR_OF_SHIFT 
	LOGICAL DO_OPAC_SHIFT
	LOGICAL DO_J_SHIFT
	LOGICAL DO_ETA_SHIFT
	LOGICAL RD_IN_SHIFTS
	LOGICAL RD_R_MULTI_GRID
	LOGICAL OUTPUT_3D_CHI
	LOGICAL OUTPUT_3D_J
	LOGICAL DEBUG 
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER MIN_INDX_NU
	INTEGER MAX_INDX_NU
	INTEGER I,J,K,ML
	INTEGER LU_IN
	INTEGER GET_INDX_DP
	EXTERNAL GET_INDX_DP
!
	CHARACTER(LEN=20) BETA_LAW
	CHARACTER(LEN=20) XAXIS,YAXIS,PLT_TITLE
	CHARACTER(LEN=80) TMP_STR,CURVE_LAB
!
	WRITE(6,*)' '
	CALL USR_LIST_OPT(LAM_START,'LAMS','4000.0','Start vacuum wavelength in Ang')
	CALL USR_LIST_OPT(LAM_END,  'LAME','7000.0','End vacuum wavelength in Ang')
	CALL USR_LIST_OPT(TOP_BOT_SYM,'TOP_BOT_SYM','T','BETA space only from 0 to PI/2')
	CALL USR_LIST_OPT(RAND_IN_BETA,'RAND_IN_BETA','T','Set variables to be randomly shifted in beta?')
	IF(RAND_IN_BETA)THEN
	  CALL USR_LIST_OPT(INIT_SEED,'INIT_SEED','245','Seed number for RNG')
	END IF
	CALL USR_LIST_OPT(HALF_SH_LEN_SCL_FAC,'HALF_SH_LEN_SCL_FAC','1.0','Scale shift by ? ( x Half shell length )')
	CALL USR_LIST_OPT( DIV_BY_SQ_DEN_SM,'DIV_BY_SQ_DEN_SM','F','Divide wind by den_sm^2?')
	CALL USR_LIST_OPT(ADD_N_PTS,'ADD_N_PTS','5','Add n pts between pts in ICM')
	CALL USR_LIST_OPT(INSERT_PTS_IN_DIR_OF_SHIFT,'INS_DIR_SHIFT','F','Insert points in direction of shift')
	CALL USR_LIST_OPT(ES_RES_KMS,'ES_RES_KMS','200','Resolution for electron scattering J in km/s')
	CALL USR_LIST_OPT(NPHI_MULTI,'NPHI','25','Number of pts in Phi')
	NB_MULTI =2 
	DO WHILE(MOD(NB_MULTI,2).EQ.0)
	  CALL USR_LIST_OPT(NB_MULTI,'NB_MULTI','11','Number of pts in BETA -- must be odd')
	  IF(MOD(NB_MULTI,2).EQ.0)WRITE(6,*)'NB_MULTI cannot be even'
	  N_BETA_PARAMS=1
	END DO
	CALL USR_LIST_OPT(BETA_LAW,'BETA_LAW','UNIFORM_BETA','Law for choosing BETA (UNIFORM_BETA, POW, BMIN)')
	IF(BETA_LAW .EQ. 'POWER')THEN
	  N_BETA_PARAMS=1
	  CALL USR_LIST_OPT(BETA_PARAMS(1),'BP_1','1.5,','Parameter 1 for BETA_LAW=POWER')
	ELSE IF(BETA_LAW .EQ. 'BMIN')THEN
	  N_BETA_PARAMS=2
	  CALL USR_LIST_OPT(BETA_PARAMS(1),'BP_1','3',   'Parameter 1 for BETA_LAW=BMIN')
	  CALL USR_LIST_OPT(BETA_PARAMS(1),'BP_2','1.5,','Parameter 2 for BETA_LAW=BMIN')
	END IF

	CALL USR_LIST_OPT(RD_R_MULTI_GRID,'RD_R_GRID','F','Read in new R grid from R_MULTI_GRID?')
	CALL USR_LIST_OPT(RD_IN_SHIFTS, 'RD_SHIFTS','F','Read in shifts for each BETA ?')
!
	CALL USR_LIST_OPT(OUTPUT_3D_J,'OUT_3D_J','F','Output LINE_MOM_DATA?')!
	IF(OUTPUT_3D_J)THEN
	  DO_J_SHIFT=.TRUE.
	ELSE
	  CALL USR_LIST_OPT(DO_J_SHIFT,'DO_J_SHIFT','F','Apply shell shft to J')
	END IF
!
	CALL USR_LIST_OPT(OUTPUT_3D_CHI,'OUT_3D_CHI','F','Output 3D_DATA?')
	IF(OUTPUT_3D_CHI)THEN
	  DO_ETA_SHIFT=.TRUE.
	ELSE
	  CALL USR_LIST_OPT(DO_ETA_SHIFT,'DO_ETA_SHIFT','F','Apply shift to CHI and ETA data?')
	END IF
	CALL USR_LIST_OPT(DEBUG,'DEBUG','F','Verbose output for debugging')
	WRITE(6,*)' '
!
! Find new wavelength range to reduce NU grid.
!
        MIN_FREQ=2.99792E+03_LDP/MAX(LAM_START,LAM_END)
        MAX_FREQ=2.99792E+03_LDP/MIN(LAM_START,LAM_END)
!
        MIN_INDX_NU=GET_INDX_DP(MAX_FREQ,NU,NCF)-1;   MIN_INDX_NU=MAX(MIN_INDX_NU,1)
        MAX_INDX_NU=GET_INDX_DP(MIN_FREQ,NU,NCF)+1;   MAX_INDX_NU=MIN(NCF,MAX_INDX_NU)
	MAX_FREQ=NU(MIN_INDX_NU); MIN_FREQ=NU(MAX_INDX_NU)
	NCF_MULTI=MAX_INDX_NU-MIN_INDX_NU+1
!
	IF(ALLOCATED(FREQ_MULTI))DEALLOCATE(FREQ_MULTI)
	ALLOCATE(FREQ_MULTI(NCF_MULTI))
	FREQ_MULTI(1:NCF_MULTI)=NU(MIN_INDX_NU:MAX_INDX_NU)
!
! Get grid on which electron scattered J is stored. ES_RES_KMS is the the minimum
! resolution. Because electron scattering has a broad scattering kernal, the frequency
! resolution can be much lower than for ETA and CHI. Grid is assumed to be equally
! spaced in frequency.
!
	DEL_NU=MIN_FREQ*ES_RES_KMS/3.0E+05_LDP
	NCF_CRSE=(MAX_FREQ-MIN_FREQ)/DEL_NU+1
	DEL_NU=(MAX_FREQ-MIN_FREQ)/(NCF_CRSE-1)
	IF(ALLOCATED(FREQ_CRSE))DEALLOCATE(FREQ_CRSE)
	ALLOCATE(FREQ_CRSE(NCF_CRSE))
	DO ML=1,NCF_CRSE
	  FREQ_CRSE(ML)=MAX_FREQ-(ML-1)*DEL_NU
	END DO
!
	WRITE(6,*)' '
	WRITE(6,*)'Freq (10^15 Hz) range is',NU(MIN_INDX_NU),NU(MAX_INDX_NU)
	WRITE(6,*)'New NCF is',NCF_MULTI
	WRITE(6,*)'New ESE NCF is',NCF_CRSE
!
	WRITE(6,'(/,1X,A)')'Computing Beta grid'
	IF(ALLOCATED(BETA_MULTI))DEALLOCATE(BETA_MULTI)
	ALLOCATE(BETA_MULTI(NB_MULTI))
	CALL COMPUTE_BETA_CL_V2(BETA_MULTI,NB_MULTI,TOP_BOT_SYM,BETA_LAW,BETA_PARAMS,N_BETA_PARAMS)
	WRITE(6,*)'Computed Beta grid'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Do some prep work for shifting shells (i.e. set up arrays, find where shells 
! are located, and measure half shell length)
!
! We alloccate ND storage locatios, since the number of shells is certainly less than ND.
!
	IF(ALLOCATED(SHELL_LOC))DEALLOCATE(SHELL_LOC)
	ALLOCATE (SHELL_LOC(ND))
	CALL DETERM_CLUMP_POS_V3_JH(MASS_DENSITY,ND,NCLUMPS,D_SH_INDX,SHELL_LOC,DEBUG)
!
! Here we insert grid points  in the ICM in order to shift the shells. if
! USE_SCALE_SIGN_TO_ADD_PTS (JNK_LOG) = T, then I insert grid pts in the direction of
! which the shell is shifted
!
	IF(RD_R_MULTI_GRID)THEN
	  WRITE(6,*)'Size of R_MULTI is',SIZE(R_MULTI)
	  OPEN(UNIT=LU_IN,FILE='R_MULTI_GRID',STATUS='OLD',ACTION='READ')
	    READ(LU_IN,*)ND_MULTI
	    IF(ALLOCATED(R_MULTI))DEALLOCATE(R_MULTI)
	    ALLOCATE(R_MULTI(ND_MULTI))
	    DO I=1,ND_MULTI
	      READ(LU_IN,*)R_MULTI(I)
	    END DO
	  CLOSE(UNIT=LU_IN)
	ELSE
	  T3 = 1.0_LDP
	  IF( .NOT. RAND_IN_BETA ) T3 = HALF_SH_LEN_SCL_FAC
	  IF(ALLOCATED(R_MULTI))DEALLOCATE(R_MULTI)
	  ALLOCATE(R_MULTI(50*ND))
	  CALL INSERT_PTS_FOR_3D(ND,R,NCLUMPS,D_SH_INDX,SHELL_LOC,
	1                        SIZE(R_MULTI),R_MULTI,ADD_N_PTS,
	1                        INSERT_PTS_IN_DIR_OF_SHIFT,T3,ND_MULTI, IERROR, DEBUG)
	  IF(IERROR .NE. 0)THEN
	    WRITE(6,*)'Exiting shift routine due to major error'
	    RETURN
	  END IF
	  WRITE(6,*)'ND_MULTI = ',ND_MULTI
	END IF
!
! Here we measure the half width of each shell and store it.
! Note: DH_SH_INDX is same for each shell since all shells share the same
!         profile.
!
	IF(ALLOCATED(HALF_SH_LEN))DEALLOCATE(HALF_SH_LEN)
	ALLOCATE (HALF_SH_LEN(NCLUMPS))
!
	HALF_SH_LEN(1:NCLUMPS) = 0.0
	DO I = 1,NCLUMPS
          K = SHELL_LOC(I)
          J = SHELL_LOC(I)-(D_SH_INDX-4)
          HALF_SH_LEN(I) = ABS(R(K) - R(J))
	END DO
!
	IF(ALLOCATED(ESEC_MULTI))THEN
	  DEALLOCATE(ESEC_MULTI,J_MULTI,VR_MULTI,SIGMA_MULTI)
	  DEALLOCATE(TEMP_MULTI)
	END IF
!
	ALLOCATE (ESEC_MULTI(ND_MULTI,NB_MULTI))
	ALLOCATE (VR_MULTI(ND_MULTI,NB_MULTI))
	ALLOCATE (TEMP_MULTI(ND_MULTI,NB_MULTI))
	ALLOCATE (SIGMA_MULTI(ND_MULTI,NB_MULTI))
	ALLOCATE (TMP_VEC(ND_MULTI))
	ALLOCATE (R_SHIFT(ND,NCLUMPS))
!
! For each angle, we shift the shift using RNG 
!
	WRITE(6,*)'Seting random seed, NBETA=',NB_MULTI
	CALL RANDOM_SEED(INIT_SEED)
	CALL CREATE_R_SHIFT(R,R_SHIFT,SHELL_LOC,HALF_SH_LEN,HALF_SH_LEN_SCL_FAC,
	1                   RD_IN_SHIFTS, D_SH_INDX, NB_MULTI, NCLUMPS, ND, DEBUG)
	WRITE(6,*)'Calculated R shift vectoris'
!
! We do the shifting on shells for ED, to calculated ESEC_MULTI
!
	T1=45.0_LDP/ATAN(1.0_LDP)		!Radians to degrees.
	DO J=1,NB_MULTI
	  CALL NEW_SHIFT(ND,R,ESEC,R_SHIFT(1,J),
	1                    NCLUMPS,D_SH_INDX,SHELL_LOC(1:NCLUMPS),
	1                    DIV_BY_SQ_DEN_SM,ND_MULTI,
	1                    R_MULTI,ESEC_MULTI(1,J)) 
!
	  IF(J .LE. 3 .OR. MOD(J,5) .EQ. 0)THEN
	    WRITE(TMP_STR,'(F5.1)')BETA_MULTI(J)*T1; TMP_STR=ADJUSTL(TMP_STR)
	    WRITE(CURVE_LAB,'(I3,A)')J,')'; CURVE_LAB='('//ADJUSTL(CURVE_LAB)
	    CURVE_LAB=TRIM(TMP_STR)//CURVE_LAB
	    CALL DP_CURVE_LAB(ND_MULTI,R_MULTI,ESEC_MULTI(1,J),CURVE_LAB)
	  END IF
	END DO
	WRITE(6,'(1X,A,/)')'Calulated shifted electron scattering curves - will need to take logs for plotting'
	CALL GRAMON_PGPLOT('R_MULTI','ESEC',' ',' ')
!
! Now lin interp V, SIGMA, and T onto new grid. 
! 
! NOTE: T uses OLD_R_SHIT to interp onto new grid since T is connected to ED, 
!       V and SIGMA are not. 
!
	DO J=1,NB_MULTI
          CALL LIN_INTERP(R_MULTI,TEMP_MULTI(1,J),ND_MULTI,R_SHIFT(1,J),T,ND)
        END DO
	CALL LIN_INTERP(R_MULTI,VR_MULTI,ND_MULTI,R,V,ND)  
	DO J=1,NB_MULTI; VR_MULTI(:,J)=VR_MULTI(:,1); END DO
        CALL LIN_INTERP(R_MULTI,SIGMA_MULTI,ND_MULTI,R,SIGMA,ND)
	DO J=1,NB_MULTI; SIGMA_MULTI(:,J)=SIGMA_MULTI(:,1); END DO
!
! Now do the shifting on shells for each frequency index for J_MULTI
!
	 IF(DO_J_SHIFT)THEN
!
! We first reduce the size of J in frequency space.
!
	  WRITE(6,*)'Creating smaller J set'
	  ALLOCATE (J_CRSE(NCF_CRSE,ND_ES))
	  ALLOCATE (TMP_MAT(NCF_MULTI,ND_ES))
	  DO I=1,ND_ES
	    DO ML=MIN_INDX_NU,MAX_INDX_NU
	      TMP_MAT(ML-MIN_INDX_NU+1,I)=RJ_ES(I,ML)
	    END DO
	  END DO
	  WRITE(6,*)MAXVAL(TMP_MAT),MINVAL(TMP_MAT)
	  CALL MON_INTERP_FAST_V2(J_CRSE,NCF_CRSE,ND_ES,FREQ_CRSE,NCF_CRSE,
	1                 TMP_MAT,NCF_MULTI,FREQ_MULTI,NCF_MULTI,'Interp. J in nu')
	  WRITE(6,*)'Created smaller J set'
	  DEALLOCATE(TMP_MAT); ALLOCATE(TMP_MAT(ND_ES,NCF_CRSE))
	  TMP_MAT=TRANSPOSE(J_CRSE); DEALLOCATE(J_CRSE)	
	  WRITE(6,*)MAXVAL(TMP_MAT),MINVAL(TMP_MAT)
!  	
	  WRITE(6,*)'Doing J shift'
	   ALLOCATE (J_MULTI(NCF_CRSE,NB_MULTI,ND_MULTI))
	   ALLOCATE (YV(MAX(ND,ND_ES)))
	   DO J=1,NB_MULTI
	     DO ML=1,NCF_CRSE
	       DO I=1,ND
	          YV(I)=TMP_MAT(I,ML)
	       END DO
	       CALL NEW_SHIFT(ND_ES,R_ES,YV,R_SHIFT(1,J),
	1                      NCLUMPS,D_SH_INDX,SHELL_LOC(1:NCLUMPS),
	1                      DIV_BY_SQ_DEN_SM,ND_MULTI,
	1                      R_MULTI,TMP_VEC)
	       J_MULTI(ML,J,:) = TMP_VEC(:)
	     END DO
	   END DO
	   DEALLOCATE(TMP_MAT)
	 END IF
!
	IF(OUTPUT_3D_J)THEN
	  K=ND_MULTI+NC
	  WRITE(6,*)'Outputting shifted J'
	  CALL WRITE_LINE_DATA_V3(R_MULTI,VR_MULTI,TEMP_MULTI,SIGMA_MULTI,
	1                            ESEC_MULTI,BETA_MULTI,
	1                            J_MULTI,FREQ_CRSE,NCF_CRSE,
	1                            ND_MULTI,NC,K,NPHI_MULTI,NB_MULTI)
	END IF
!
	IF(DO_ETA_SHIFT)THEN
	  WRITE(6,*)'Doing ETA shift'
	  IF(ALLOCATED(ETA_MULTI))DEALLOCATE(ETA_MULTI)
	  ALLOCATE(ETA_MULTI(ND_MULTI,NB_MULTI,NCF_MULTI))
	  DO J=1,NB_MULTI
	    DO ML=MIN_INDX_NU,MAX_INDX_NU
	      DO I=1,ND
	        YV(I)=ETA(I,ML)
	      END DO
	      CALL NEW_SHIFT(ND,R,YV,R_SHIFT(1,J),
	1               NCLUMPS,D_SH_INDX,SHELL_LOC,
	1               DIV_BY_SQ_DEN_SM,ND_MULTI,
	1               R_MULTI,TMP_VEC) 
	      K=ML-MIN_INDX_NU+1
	      ETA_MULTI(1:ND_MULTI,J,K) = TMP_VEC(1:ND_MULTI)
	    END DO
	  END DO
!
	  WRITE(6,*)'Doing CHI shift'
	  IF(ALLOCATED(CHI_MULTI))DEALLOCATE(CHI_MULTI)
	  ALLOCATE(CHI_MULTI(ND_MULTI,NB_MULTI,NCF_MULTI))
	  DO J=1,NB_MULTI
	    DO ML=MIN_INDX_NU,MAX_INDX_NU
	      DO I=1,ND
	        YV(I)=CHI(I,ML)
	      END DO
	      CALL NEW_SHIFT(ND,R,YV,R_SHIFT(1,J),
	1                     NCLUMPS,D_SH_INDX,SHELL_LOC,
	1                     DIV_BY_SQ_DEN_SM,ND_MULTI,
	1                     R_MULTI,TMP_VEC)
	      K=ML-MIN_INDX_NU+1
	      CHI_MULTI(1:ND_MULTI,J,K) = TMP_VEC(1:ND_MULTI)
	    END DO
	  END DO
	  DEALLOCATE(TMP_VEC)
	END IF
!
	IF(OUTPUT_3D_CHI)THEN
	  CALL TUNE(1,'3D_CHI_OUT')
	  CALL WRITE_MULTI_FULL_CHI_JH(ND_MULTI,NB_MULTI,NCF_MULTI,R_MULTI,
	1                              FREQ_MULTI,BETA_MULTI,VR_MULTI,
	1                              SIGMA_MULTI,TEMP_MULTI,ESEC_MULTI,
	1                              ETA_MULTI,CHI_MULTI)
	  CALL TUNE(2,'3D_CHI_OUT')
	END IF
	CALL TUNE(3,' ')
!
! Clean up sotrage locations.
!
	IF(ALLOCATED(YV))DEALLOCATE(YV)
!
! We deallocate all data variables defined in CLUMP_SHIFT_MODULE. These
! could be deleted if we wish to have the variables accessed elsewhere.
!
        IF(ALLOCATED(R_MULTI))DEALLOCATE(R_MULTI)
!       
        IF(ALLOCATED(BETA_MULTI))DEALLOCATE(BETA_MULTI)
        IF(ALLOCATED(FREQ_MULTI))DEALLOCATE(FREQ_MULTI)
        IF(ALLOCATED(FREQ_CRSE))DEALLOCATE(FREQ_CRSE)
!         
        IF(ALLOCATED(ETA_MULTI))DEALLOCATE(ETA_MULTI)
        IF(ALLOCATED(CHI_MULTI))DEALLOCATE(CHI_MULTI)
        IF(ALLOCATED(J_MULTI))DEALLOCATE(J_MULTI)
        IF(ALLOCATED(ESEC_MULTI))DEALLOCATE(ESEC_MULTI)
!
        IF(ALLOCATED(VR_MULTI))DEALLOCATE(VR_MULTI)
        IF(ALLOCATED(SIGMA_MULTI))DEALLOCATE(SIGMA_MULTI)
        IF(ALLOCATED(TEMP_MULTI))DEALLOCATE(TEMP_MULTI)
	IF(ALLOCATED(R_SHIFT))DEALLOCATE(R_SHIFT)
	IF(ALLOCATED(SHELL_SHIFT))DEALLOCATE(SHELL_SHIFT)
	IF(ALLOCATED(HALF_SH_LEN))DEALLOCATE(HALF_SH_LEN)
	IF(ALLOCATED(SHELL_LOC))DEALLOCATE(SHELL_LOC)
!
	RETURN
	END
