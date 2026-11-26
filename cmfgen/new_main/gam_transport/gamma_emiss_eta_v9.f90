! Subroutine that will calculate the intrinsic emissivity as a function of
! depth and frequency. Depth is the first index and frequency is the first.
!
!
! Created July 11, 2014
! Edited 15 May 2017:	Now scales line energies to total energy per decay (accounting for
!			Kin. En. from leptons). Correct units should be compared.
! Edited 25 May 2017:   Doing testing to make sure no error with energy being created
!
! Edited 20 June 2017:	Added normalized gamma-ray lines as a parameter
! Edited 20 May  2023:  ETT_EMISS set to 1.0E-120 when zero.
!------------------------------------------------------------------------------------------------
!
	SUBROUTINE GAMMA_INT_EMISS_V9(ETA_EMISS,NU_GRID_VEC,NF_GRID_PTS,&
			R,ND,V_GAUSS,DECAY_KIN_E,NORM_GAM_LINES,&
			NORM_LINES_ONLY,INST_DECAY)
	USE SET_KIND_MODULE
	USE GAMMA_NUC_DECAY_V2
	USE NUC_ISO_MOD
!	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
	INTEGER :: I,J,K,L,ML,LS
	INTEGER :: NF_GRID_PTS,ND
	INTEGER :: JK
!
	LOGICAL :: NORM_GAM_LINES
	LOGICAL :: NORM_LINES_ONLY
	LOGICAL :: INST_DECAY
!
	REAL(KIND=LDP) :: NU_GRID_VEC(NF_GRID_PTS)
	REAL(KIND=LDP) :: ETA_EMISS(ND,NF_GRID_PTS)
	REAL(KIND=LDP) :: R(ND)
	REAL(KIND=LDP) :: DECAY_KIN_E(ND)
	REAL(KIND=LDP) :: GAUSS_INT
	REAL(KIND=LDP) :: V_GAUSS
	REAL(KIND=LDP) :: T1,T2,T3,T4,T5,T6,T7
	REAL(KIND=LDP) :: TOT_EGAM,DE
	REAL(KIND=LDP) :: DELTA_T
	REAL(KIND=LDP) :: PI
	REAL(KIND=LDP) :: FOURPI
	REAL(KIND=LDP) :: NORM_FAC
	REAL(KIND=LDP) :: DEL_T
	REAL(KIND=LDP), PARAMETER :: PLANCK = 4.135668E-021_LDP ! UNITS OF MeV*s SINCE PHOTON ENERGIES IN MeV
	REAL(KIND=LDP), PARAMETER :: MeV_ERG = 1.602177E-06_LDP
	REAL(KIND=LDP), PARAMETER :: ERGS_TO_MEV = 624150.9
	EXTERNAL GAUSS_INT
!
! Testing variables
!
	REAL(KIND=LDP) :: WORK1_EN(ND)
	REAL(KIND=LDP) :: DECAY_EN
	REAL(KIND=LDP) :: WORK2_EN(ND)
	REAL(KIND=LDP) :: DECAY_ETA_EN
	REAL(KIND=LDP), ALLOCATABLE :: TA(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB(:)
	REAL(KIND=LDP), ALLOCATABLE :: TC(:)
!
!
! I will need to modify this to only use the lines with the species present which I can use with an IF
! statement later on
!
	ETA_EMISS=0.0_LDP
	DECAY_KIN_E=0.0_LDP
	PI=ACOS(-1.0_LDP)
	FOURPI=4.0_LDP*ACOS(-1.0_LDP)
	WRITE(6,'(/,A)')" Using gamma_int_emiss_v9.f90 to compute the gamma-ray emission"
!
! DEL_T is saved from the subroutine DO_SPECIES_DECAYS
!
	DEL_T=NUC_ISO_MOD_DELTA_T
	WRITE(6,*)"DELTA_T:",DEL_T
!
! Looping through frequency then elements read in by RD_NUC_DECAY_DATA_GAM and
! READ_NUC_DECAY_DATA_V3
!
	WRITE(6,'(A)')' Calculating and normalizing energies of G-ray lines...'
!
!
! There is a problem with line energies and Kin. En. not summing
! to the total energy per decay. We account for the difference
! below by scaling the line energies by the fractional difference
! between the total energy as summed from lines and Kin. En.
! and the read in ENERGY_PER_DECAY (same units must be compared)
!
	IF(NORM_GAM_LINES)THEN
20	  FORMAT(1X,A7,1X,A10,2X,4(A10,6X),A10)
21	  FORMAT(1X,A7,1X,F12.4,4(F10.5,6X),F10.5)
	  WRITE(6,20)'Species','Mass','Kin. E','TOT_EGAM','E/DECAY','%diff','NORM_FAC'
	  DO K=1,N_SPECIES
	    DO L=1,NUM_DECAY_PATHS
	      IF(GAM_ISO(K)%SPECIES .EQ. NUC(L)%SPECIES .AND. &
			GAM_ISO(K)%ATOMIC_MASS .EQ. NUC(L)%MASS)THEN
	        LS=NUC(L)%LNK_TO_ISO
	        TOT_EGAM=0.0_LDP
	        NORM_FAC=1.0E+0_LDP
	        DO ML=1,GAM_ISO(K)%NUM_GAMMA
	          T1=GAM_ISO(K)%E_GAMMA(ML)
		  T2=GAM_ISO(K)%PROB(ML)
		  TOT_EGAM=TOT_EGAM+T1*T2
	        END DO
!
	        T1=NUC(L)%ENERGY_PER_DECAY*ERGS_TO_MEV
!
! Normalize just the gamma-ray lines or both the lines and the kinetic
! energy from positrons
!
		IF(NORM_LINES_ONLY)THEN
		  T1=T1-GAM_ISO(K)%KIN_ENERGY
	          IF(T1 .GT. 0.0_LDP)THEN
		    DE=ABS(TOT_EGAM-T1)/T1
	          ELSE
		    DE=0.0_LDP
	          END IF
	          IF(DE .GT. 1.0E-05_LDP .AND. TOT_EGAM .NE. 0.0_LDP)THEN
		    NORM_FAC=T1/TOT_EGAM
	          END IF
	          T2=GAM_ISO(K)%ATOMIC_MASS
	          T3=GAM_ISO(K)%KIN_ENERGY
	          WRITE(6,21)GAM_ISO(K)%SPECIES,T2,T3,&
		  	TOT_EGAM,T1,DE,NORM_FAC
	          GAM_ISO(K)%PROB=GAM_ISO(K)%PROB*NORM_FAC
		ELSE
		  TOT_EGAM=(TOT_EGAM+GAM_ISO(K)%KIN_ENERGY)
		  IF(T1 .GT. 0.0_LDP)THEN
		    DE=ABS(TOT_EGAM-T1)/T1
		  ELSE
		    DE=0.0_LDP
	          END IF
	          IF(DE .GT. 1.0E-05_LDP .AND. TOT_EGAM .NE. 0.0_LDP)THEN
		    NORM_FAC=T1/TOT_EGAM
	          END IF
	          T2=GAM_ISO(K)%ATOMIC_MASS
	          T3=GAM_ISO(K)%KIN_ENERGY
	          WRITE(6,21)GAM_ISO(K)%SPECIES,T2,T3,&
		  	TOT_EGAM,T1,DE,NORM_FAC
	          GAM_ISO(K)%PROB=GAM_ISO(K)%PROB*NORM_FAC
		  GAM_ISO(K)%KIN_ENERGY=GAM_ISO(K)%KIN_ENERGY*NORM_FAC
		END IF
	      END IF
	    END DO !decay paths loop
	  END DO !species loop
	END IF
!
! End of normalizing gamma-ray lines
!
!
	DO J=1,NF_GRID_PTS
	  T2=NU_GRID_VEC(J)
	  DO I=1,ND
	    DO K=1,N_SPECIES
	      DO L=1,NUM_DECAY_PATHS
		IF(GAM_ISO(K)%SPECIES .EQ. NUC(L)%SPECIES .AND. &
			GAM_ISO(K)%ATOMIC_MASS .EQ. NUC(L)%MASS)THEN
		  LS=NUC(L)%LNK_TO_ISO
		  DO ML=1,GAM_ISO(K)%NUM_GAMMA
		    T1=GAM_ISO(K)%E_GAMMA(ML)/PLANCK
		    T3=GAUSS_INT(T2,T1,V_GAUSS)
		    T4=GAM_ISO(K)%PROB(ML)*GAM_ISO(K)%E_GAMMA(ML)*T3
		    IF(INST_DECAY)THEN
		      ETA_EMISS(I,J)=ETA_EMISS(I,J)+T4*ISO(LS)%DECAY_LUM(I)
		    ELSE
		      ETA_EMISS(I,J)=ETA_EMISS(I,J)+T4*ISO(LS)%NUM_DECAYS(I)
		    END IF
		  END DO
		END IF
	      END DO !decay paths loop
	    END DO !species loop
	    IF(ETA_EMISS(I,J) .EQ. 0.0_LDP)ETA_EMISS(I,J)=1.0E-120_LDP
	  END DO !depth loop
	END DO !frequency loop
!
! Next piece is to compare the used isotopes and store the decay kinetic energy read in from
! GAM_NUC_DECAY_DATA_SUB_V3
!
	DO I=1,ND
	  DO K=1,N_SPECIES
	    DO L=1,NUM_DECAY_PATHS
	      IF(GAM_ISO(K)%SPECIES .EQ. NUC(L)%SPECIES .AND. &
			GAM_ISO(K)%ATOMIC_MASS .EQ. NUC(L)%MASS)THEN
		LS=NUC(L)%LNK_TO_ISO
		IF(INST_DECAY)THEN
		  DECAY_KIN_E(I)=DECAY_KIN_E(I)+ISO(LS)%DECAY_LUM(I)*GAM_ISO(K)%KIN_ENERGY
		ELSE
		  DECAY_KIN_E(I)=DECAY_KIN_E(I)+ISO(LS)%NUM_DECAYS(I)*GAM_ISO(K)%KIN_ENERGY
		END IF
	      END IF
	    END DO !decay paths loop
	  END DO !species loop
	END DO !depth loop
!
!==========================================================================================
!
! Following part is used to test/diagnose if energy is being conserved
!
	WORK1_EN=0.0_LDP
	DECAY_EN=0.0_LDP
	DO I=1,ND
	  DO J=1,NUM_DECAY_PATHS
	    LS=NUC(J)%LNK_TO_ISO
	    IF(INST_DECAY)THEN
	      T1=ISO(LS)%DECAY_LUM(I)
	    ELSE
	      T1=ISO(LS)%NUM_DECAYS(I)
	    END IF
	    T2=NUC(J)%ENERGY_PER_DECAY
	    WORK1_EN(I)=WORK1_EN(I)+T1*T2
	  END DO
	END DO
!	DO I=1,ND-1
!	  T1=R(I)-R(I+1)
!	  T2=R(I)*R(I)
!	  T3=R(I+1)*R(I+1)
!	  DECAY_EN=DECAY_EN+0.5D0*T1*(WORK1_EN(I)*T2+WORK1_EN(I+1)*T3)
!	END DO
	ALLOCATE(TA(ND))
	TA=0.0_LDP
	DO I=1,ND
	  TA(I)=WORK1_EN(I)*R(I)*R(I)
	END DO
	CALL LUM_FROM_ETA(TA,R,ND)
	DECAY_EN=SUM(TA)
	IF(INST_DECAY)THEN
	  DECAY_EN=DECAY_EN*FOURPI*1.0E+30_LDP
	ELSE
	  DECAY_EN=DECAY_EN*FOURPI/DEL_T*1.0E+30_LDP
	END IF
	ALLOCATE(TB(NF_GRID_PTS))
	TB=0.0_LDP
	WORK2_EN=0.0_LDP
	DO I=1,ND
!	  DO J=1,NF_GRID_PTS-1
!	    T1=ABS(NU_GRID_VEC(J)-NU_GRID_VEC(J+1))
!	    T2=ETA_EMISS(I,J)+ETA_EMISS(I,J+1)
!	    WORK2_EN(I)=WORK2_EN(I)+0.5D0*T1*T2
!	  END DO
	  DO J=1,NF_GRID_PTS
	    TB(J)=ETA_EMISS(I,J)
	  END DO
	  CALL LUM_FROM_ETA(TB,NU_GRID_VEC,NF_GRID_PTS)
	  WORK2_EN(I)=SUM(TB)
	END DO
	WORK2_EN=(WORK2_EN+DECAY_KIN_E)*MeV_ERG
	DECAY_ETA_EN=0.0_LDP
!	DO I=1,ND-1
!	  T1=R(I)-R(I+1)
!	  T2=R(I)*R(I)
!	  T3=R(I+1)*R(I+1)
!	  DECAY_ETA_EN=DECAY_ETA_EN+0.5D0*T1*(WORK2_EN(I)*T2+WORK2_EN(I+1)*T3)
!	END DO
	TA=0.0_LDP
	DO I=1,ND
	  TA(I)=WORK2_EN(I)*R(I)*R(I)
	END DO
	DECAY_ETA_EN=SUM(TA)
	IF(INST_DECAY)THEN
	  DECAY_ETA_EN=DECAY_ETA_EN*FOURPI*1.0E+30_LDP
	ELSE
	  DECAY_ETA_EN=DECAY_ETA_EN*FOURPI/DEL_T*1.0E+30_LDP
	END IF
	OPEN(UNIT=7,FILE='check_decays_energy_compare.dat',ACTION='WRITE',&
	     STATUS='UNKNOWN')
	WRITE(7,'(A)')'! Comparing the direct energy obtained from number of decays'
	WRITE(7,'(A)')'! and the energy per decay to the energy obtained by'
	WRITE(7,'(A)')'! integrating the isotropic gamma-ray emissision before we'
	WRITE(7,'(A)')'! put it in the correct units (i.e. /s/4pi)'
	WRITE(7,'(A)')'! Energies should be in ergs/cm^3'
	WRITE(7,'(A,ES16.6)')'! Total nuclear Luminosity counted from decays in ergs/s:',&
				DECAY_EN
	WRITE(7,'(A,ES16.6)')'! Total nuclear Luminosity integrated from iso emiss in ergs/s:',&
				DECAY_ETA_EN
	WRITE(7,'(A1,A16,1X,A16,1X,A10,3X,A10)')'!','From decays','From em. integ.',&
		'% diff','frac. diff'
	DO I=1,ND
	  T1=ABS(WORK1_EN(I)-WORK2_EN(I))/WORK1_EN(I)
	  WRITE(7,'(1X,ES16.6,1X,ES16.6,1X,F10.4,3X,F10.4)')WORK1_EN(I),WORK2_EN(I),&
		T1*1.E+02_LDP,WORK1_EN(I)/WORK2_EN(I)
	END DO
	CLOSE(UNIT=7)
!
!==========================================================================================
!
	IF(INST_DECAY)THEN
	  ETA_EMISS=ETA_EMISS*MeV_ERG/4.0_LDP/PI
	  DECAY_KIN_E=DECAY_KIN_E*MeV_ERG
	ELSE
	  ETA_EMISS=ETA_EMISS*MeV_ERG/DEL_T/4.0_LDP/PI
	  DECAY_KIN_E=DECAY_KIN_E*MeV_ERG/DEL_T
	END IF
!
	DEALLOCATE(TA,TB)
!
	RETURN
	END SUBROUTINE GAMMA_INT_EMISS_V9
!
! Function is to ensure the emission preserves gaussian profiles with given
! spread based on velocity
!
	FUNCTION GAUSS_INT(X,XO,VAR)
	USE SET_KIND_MODULE
	IMPLICIT NONE
	REAL(KIND=LDP) ::  GAUSS_INT
	REAL(KIND=LDP) :: MU
	REAL(KIND=LDP) :: X,XO
	REAL(KIND=LDP), PARAMETER :: C=299792 !km/s
	REAL(KIND=LDP) :: SIGMA,VAR ! VAR in km/s
	REAL(KIND=LDP), PARAMETER :: ONE=1.0_LDP
	REAL(KIND=LDP), PARAMETER :: TWO=2.0_LDP
	REAL(KIND=LDP) :: PI
!
	PI=ACOS(-ONE)
	SIGMA=VAR*XO/C
!
	GAUSS_INT=EXP( -(X-XO)*(X-XO)/TWO/SIGMA/SIGMA )/SIGMA/SQRT(TWO*PI)
	RETURN
	END FUNCTION
