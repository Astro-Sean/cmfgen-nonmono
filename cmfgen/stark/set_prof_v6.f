! Routine to compute:
!
!          (1) Doppler line profiles for any species
!          (2) Stark and Voigt profiles (of varrying accuracy) for individual lines.
!          (3) Approximate STARK profiles for any Hydrogenic species of charge
!                Z. The theory is excellent for high Balmer lines, but only
!                approximate for Ha. Uses the GRIEM theory as modified by Auer
!                and Mihalas AP J S 24 1972.
!
!
! Output:
!        PROF - Profile as a function of depth. Either:
!                    (a) Doppler profile or
!                    (b) STARK profile convolved with a Doppler profile.
!
!                    The Doppler profile can have a turbulent contribution
!                    given by VTURB.
! Input:
!       NU         - Freqency (in units of 10^15 Hz)
!       ED_IN      - Electron density (/cm^3) (Vector, length ND)
!       POP_PROTON - Proton density (/cm^3) (Vector, length ND)
!	POP_HEPLUS - He+ density (/cm^3) (Vector, length ND)
!       TEMP_IN    - Temperature (10^4 K)     (Vector, length ND)
!       Z_POP      - Level charge (not charge on ion) [=0 for H levels].
!       VTURB_IN   - Turbulent velocity in km/s (function of depth).
!       ND         - Nuber of (Ne,T) values profile is to be computed for.
!
!	VTURB		!Turbulent velocity (km/s): same at all depths
!	TDOP            !Temp for Doppler profile (if same at all depths).
!	AMASS_DOP	!Same for all species when used.
!	MAX_PROF_ED	!Limits witdt of Voigt/Stark profiles
!
	SUBROUTINE SET_PROF_V6(LINE_PROF_SIM,AMASS_SIM,
	1             RESONANCE_ZONE,END_RESONANCE_ZONE,
	1             SIM_LINE_POINTER,MAX_SIM,
	1             Z_POP,ML_CUR,
	1             ED_IN,POP_PROTON,POP_HEPLUS,
	1             TEMP_IN,VTURB_IN,ND,NT,
	1             TDOP,AMASS_DOP,VTURB,MAX_PROF_ED,
	1             NORM_PROFILE,LU_STK)
	USE SET_KIND_MODULE
!
! Altered 05-Jul-2022: Revamped to help improve speed and reduce call overhead
!                         OMP used for GRIEM profile computation.
!                         Originally OMP NOT used for VOIGT section as it slowed the
!                         code. However we now store the full Voigt profile, which
!                         allows better paraellization.LOOP of lines is now contained in
!                         this routine.
!                         Based on SET_PROF_V5.
!
! Altered 18-May-2015: LOC_GAM_COL is set to C4 unless specifically set in STRK_LIST.
! Altered 20-May-2014: VERBOSE option introduced, and sone diagnostic output was modified.
! Altered 14-May-2014: MAX_PROF_ED now applies to all profiles
! Altered 31-Jan-2014: Added MAX_PROF_ED to call (changed to V5).
! Altered 06-Jan-2014: VTURB_FIX replaced by TDOP,AMASS_DOP,VTURB in call (changed to V4).
!                          (taken from cur_cm_25jun13 development version).
! Altered 03-Jan-2001: Check for unrecognized profile type.
! Altered 07-Jan-1999: ML_CUR now passed in call.
!                      Profile is now recomputed whenever ML_CUR=ML_ST,
!                      or when profile is unavailable.
!
	USE PROF_MOD
	USE MOD_FREQ_OBS
	USE MOD_STRK_LIST
!
	IMPLICIT NONE
	INTEGER ML_CUR
	INTEGER MAX_SIM
	INTEGER ND
	INTEGER NT
	INTEGER LU_STK
!
	REAL(KIND=LDP) LINE_PROF_SIM(ND,MAX_SIM)
	REAL(KIND=LDP) AMASS_SIM(MAX_SIM)
	INTEGER SIM_LINE_POINTER(MAX_SIM)
	LOGICAL RESONANCE_ZONE(MAX_SIM)
	LOGICAL END_RESONANCE_ZONE(MAX_SIM)
!
	REAL(KIND=LDP) ED_IN(ND)
	REAL(KIND=LDP) POP_PROTON(ND)
	REAL(KIND=LDP) POP_HEPLUS(ND)
	REAL(KIND=LDP) TEMP_IN(ND)
	REAL(KIND=LDP) VTURB_IN(ND)
	REAL(KIND=LDP) Z_POP(NT)
!
	REAL(KIND=LDP) VTURB		!Turbulent velocity (km/s): same at all depths
	REAL(KIND=LDP) TDOP
	REAL(KIND=LDP) AMASS_DOP	!Same for all species when used.
	REAL(KIND=LDP) MAX_PROF_ED	!Limits witdt of Voigt/Stark profiles
!
	LOGICAL NORM_PROFILE
!
	LOGICAL, PARAMETER :: RET_LOG=.FALSE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Local variables
!
	INTEGER I,IS,ML,ID,ITR
	INTEGER LOC_INDX		!Indicates which storage
	INTEGER NF
	INTEGER NL
	INTEGER NUP
	INTEGER ML_ST
	INTEGER ML_END
	INTEGER LN_PNT
	INTEGER SIM_INDX
!
	REAL(KIND=LDP), ALLOCATABLE ::  XNU(:)
	REAL(KIND=LDP) VTURB_SQ(ND)
	REAL(KIND=LDP) ED_MOD(ND)
!
	REAL(KIND=LDP) AMASS_IN
	REAL(KIND=LDP) Z_IN
	REAL(KIND=LDP) NU_ZERO
	REAL(KIND=LDP) GAM_RAD
	REAL(KIND=LDP) C4_INTER
!
	REAL(KIND=LDP) T1,T2
	REAL(KIND=LDP) TMP_ED,NU_DOP
	REAL(KIND=LDP) A_VOIGT
	REAL(KIND=LDP) V_VOIGT
	REAL(KIND=LDP) LOC_GAM_RAD
	REAL(KIND=LDP) LOC_GAM_COL
!
	LOGICAL PROF_COMPUTED
	LOGICAL VERBOSE
!
! External functions
!
	REAL(KIND=LDP) VOIGT
!
        VTURB_SQ(1:ND)=( VTURB_IN(1:ND)/12.85_LDP  )**2
!
	DO SIM_INDX=1,MAX_SIM
	  IF(RESONANCE_ZONE(SIM_INDX))THEN
!
	    AMASS_IN=AMASS_SIM(SIM_INDX)
	    LN_PNT=SIM_LINE_POINTER(SIM_INDX)
!
	    Z_IN=Z_POP(VEC_NL(LN_PNT))+1.0_LDP
	    NU_ZERO=VEC_FREQ(LN_PNT)
	    GAM_RAD=VEC_ARAD(LN_PNT)
	    C4_INTER=VEC_C4(LN_PNT)
	    NL=VEC_MNL_F(LN_PNT)
	    NUP=VEC_MNUP_F(LN_PNT)
	    ML_ST=LINE_ST_INDX_IN_NU(LN_PNT)
	    ML_END=LINE_END_INDX_IN_NU(LN_PNT)
	    PROF_COMPUTED=.FALSE.
!
! Doppler profile is the same for all species, and is the same at all depths.
!
	    IF(PROF_TYPE(LN_PNT) .EQ. 'DOP_FIX')THEN
              T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
	      T2=12.85_LDP*SQRT( TDOP/AMASS_DOP + (VTURB/12.85_LDP)**2 )
              NU_DOP=NU_ZERO*T2/C_KMS
              DO I=1,ND
                LINE_PROF_SIM(I,SIM_INDX)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	      END DO
	      PROF_COMPUTED=.TRUE.
!
! Doppler profile is the same at all depth, but not the same for all species.
! A reasonable choice for T_DOP is Teff.
!
	    ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'DOP_SPEC')THEN
              T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
	      T2=12.85_LDP*SQRT( TDOP/AMASS_IN + (VTURB/12.85_LDP)**2 )
              NU_DOP=NU_ZERO*T2/C_KMS
              DO I=1,ND
                LINE_PROF_SIM(I,SIM_INDX)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	      END DO
	      PROF_COMPUTED=.TRUE.
!
! Doppler profile varies with species and depth.
!
	    ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'DOPPLER')THEN
              T1=1.0E-15_LDP/1.77245385095516_LDP         !1.0D-15/SQRT(PI)
              T2=NU_ZERO*12.85_LDP/C_KMS
              DO I=1,ND
                NU_DOP=T2*SQRT( TEMP_IN(I)/AMASS_IN + VTURB_SQ(I) )
                LINE_PROF_SIM(I,SIM_INDX)=EXP( -( (NU(ML_CUR)-NU_ZERO)/NU_DOP )**2 )*T1/NU_DOP
	      END DO
	      PROF_COMPUTED=.TRUE.
!
	    ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'VOIGT')THEN
!
	      IF(VGT_POINTER(SIM_INDX) .EQ. 0)THEN
	        DO IS=1,NVGT_ST
	          IF(VGT(IS)%NF .EQ. 0)THEN
	            VGT(IS)%NF=ML_END-ML_ST+1
	            ALLOCATE(VGT(IS)%VGT_PROF(ND,VGT(IS)%NF))
	            VGT_POINTER(SIM_INDX)=IS
!	            WRITE(6,*)IS,SIM_INDX,ML_END-ML_ST+1; FLUSH(UNIT=6)
!
	            T2=NU_ZERO*12.85_LDP/C_KMS
	            LOC_GAM_RAD=GAM_RAD
	            LOC_GAM_COL=1.55E+04_LDP*(ABS(C4_INTER)**0.667_LDP)
	            ID=PROF_LIST_LOCATION(LN_PNT)
	            IF(ID .NE. 0)THEN
	              LOC_GAM_RAD=LST_GAM_RAD(ID)
	              LOC_GAM_COL=LST_GAM_COL(ID)
	            END IF
!	            WRITE(6,*)'Starting parallel loop'; FLUSH(UNIT=6)
!
!$OMP PARALLEL DO PRIVATE(I,ML,TMP_ED,NU_DOP,A_VOIGT,V_VOIGT)
	            DO ML=ML_ST,ML_END
	              DO I=1,ND
	                TMP_ED=MIN(ED_IN(I),MAX_PROF_ED)
                        NU_DOP=T2*SQRT( TEMP_IN(I)/AMASS_SIM(SIM_INDX) + VTURB_SQ(I) )
	                A_VOIGT=0.25E-15_LDP*(LOC_GAM_RAD+LOC_GAM_COL*TMP_ED)/PI/NU_DOP
	                V_VOIGT=(NU(ML)-NU_ZERO)/NU_DOP
	                VGT(IS)%VGT_PROF(I,ML-ML_ST+1)=1.0E-15_LDP*VOIGT(A_VOIGT,V_VOIGT)/NU_DOP
	              END DO
	            END DO
!$OMP END PARALLEL DO
!
	            EXIT
	          END IF
	        END DO
	      END IF

!	      WRITE(6,*)'Setting profile'; FLUSH(UNIT=6)
	      ML=ML_CUR-ML_ST+1
	      IS=VGT_POINTER(SIM_INDX)
!	      WRITE(6,*)IS,ML,SIM_INDX; FLUSH(UNIT=6)
	      DO I=1,ND
	         LINE_PROF_SIM(I,SIM_INDX)=VGT(IS)%VGT_PROF(I,ML)
	      END DO
	      IF(END_RESONANCE_ZONE(SIM_INDX))THEN
	        VGT(IS)%NF=0
	        DEALLOCATE (VGT(IS)%VGT_PROF)
	        VGT_POINTER(SIM_INDX)=0
	      END IF
	      PROF_COMPUTED=.TRUE.
!
!
!
! If we reach here we have two choices:
!
! (1) Profile has already been computed, hence we can use
!       the tabulated values.
! (2) We need to compute full profile, using a variety
!       of different methods.
!`
	    ELSE IF(ML_CUR .NE. ML_ST)THEN
	      DO LOC_INDX=1,NSTORE_MAX
	        IF( .NOT. STORE_AVAIL(LOC_INDX) )THEN
	          IF( NU_ZERO .EQ. NU_ZERO_STORE(LOC_INDX)    .AND.
	1             NL .EQ. NL_STORE(LOC_INDX)              .AND.
	1             NUP .EQ. NUP_STORE(LOC_INDX)             )THEN
!
!
! Check store correct, then set profile data
!
	            LST_FREQ_LOC(LOC_INDX)=LST_FREQ_LOC(LOC_INDX)+1
	            IF(NU_STORE(LST_FREQ_LOC(LOC_INDX), LOC_INDX) .NE. NU(ML_CUR))THEN
	              WRITE(LUER,*)'Error in SET_PROF'
	              WRITE(LUER,*)'Frequencies don''t match'
	              WRITE(LUER,*)NU(ML_CUR),NU_STORE(LST_FREQ_LOC(LOC_INDX),LOC_INDX)
	              STOP
	            END IF
	            LINE_PROF_SIM(1:ND,SIM_INDX)=PROF_STORE(1:ND,LST_FREQ_LOC(LOC_INDX),LOC_INDX)
	            PROF_COMPUTED=.TRUE.
!
! If possible, free up STORE:
!
	            IF( END_RESONANCE_ZONE(SIM_INDX) )THEN
	              STORE_AVAIL(LOC_INDX)=.TRUE.
	            END IF
	            EXIT
	          END IF
	        END IF
	      END DO
!
	    ELSE
!
! We always compute new profile if ML_CUR .EQ. ML_ST
!
	      DO LOC_INDX=1,NSTORE_MAX
	        IF( .NOT. STORE_AVAIL(LOC_INDX) )THEN
	          IF( NU_ZERO .EQ. NU_ZERO_STORE(LOC_INDX)    .AND.
	1             NL .EQ. NL_STORE(LOC_INDX)              .AND.
	1             NUP .EQ. NUP_STORE(LOC_INDX)             )THEN
	             STORE_AVAIL(LOC_INDX)=.TRUE.
	          END IF
	          EXIT
	        END IF
	      END DO
!
! 
!
! In this section of the routine, we compute general STARK profiles
! over the entire line profile, since it is generally computationally
! prohibitive to compute the STARK profile for each frequency. Rather
! we compute the STARK profile [f(Ne,T,Vturb) ] once, and save the data
! for subsequent calls.
!
! Determine the location where the STARK (i.e. Line) profile can be stored.
!
	      CALL GET_VERBOSE_INFO(VERBOSE)
	      DO I=1,ND
	        ED_MOD(I)=MIN(ED_IN(I),MAX_PROF_ED)
	      END DO
	      LOC_INDX=0
	      DO ML=1,NSTORE_MAX
	        IF( STORE_AVAIL(ML) )THEN
	          LOC_INDX=ML
	          EXIT
	        END IF
	      END DO
	      IF(LOC_INDX .EQ. 0)THEN
	        WRITE(LUER,*)'Error in SET_PROF --- insufficient storage'
	        WRITE(LUER,*)'Maximum number of storage locations is:',NSTORE_MAX
	        DO ML=1,NSTORE_MAX
	          WRITE(LUER,'(I4,ES14.5,2I8)')ML,NU_ZERO_STORE(ML),NL_STORE(ML),NUP_STORE(ML)
	        END DO
	        STOP
	      END IF
	      STORE_AVAIL(LOC_INDX)=.FALSE.
	      PROF_STORE(:,:,LOC_INDX)=0.0_LDP
!
! Determine how many frequencies we need to compute profile at, and check
! we have sufficent storage.
!
	      IF(ML_END-ML_ST+1 .GT. NFREQ_MAX)THEN
	        WRITE(LUER,*)'Error in SET_PROF --- NFREQ_MAX too small'
	        WRITE(LUER,*)'NFREQ_MAX=',NFREQ_MAX
	        WRITE(LUER,*)'Required NFREQ_MAX is',ML_END-ML_ST+1
	        STOP
	      END IF
!
! We can now compute the individual line profiles. New methods can
! be added to this section.
!
	      NF=ML_END-ML_ST+1
	      IF(NF .GT. NFREQ_MAX)THEN
	        WRITE(LUER,*)'Error in SET_PROF_V6 --- insufficient storage'
	        WRITE(LUER,*)'Max # of frequency storage locations is:',NFREQ_MAX
	        WRITE(LUER,*)'# of storage locatins requires is:',NF
	        STOP
	      END IF
	      ID=PROF_LIST_LOCATION(LN_PNT)		!Location in arrays
	      I=ML_END-ML_ST+1
	      IF(ALLOCATED(XNU))THEN
	        IF(SIZE(XNU) .LT. I)THEN
	          DEALLOCATE(XNU)
	          ALLOCATE(XNU(I))
	        END IF
	      ELSE
	         ALLOCATE(XNU(I))
	      END IF
	       XNU(1:NF)=NU(ML_ST:ML_END)
900	       FORMAT(X,A,T20,A,T26,F20.10,2I6,3I10,F6.1)
	       IF(PROF_TYPE(LN_PNT) .EQ. 'BS_HHE' .OR. PROF_TYPE(LN_PNT) .EQ. 'LEMKE_HI')THEN
	         IF(PROF_TYPE(LN_PNT) .EQ. 'BS_HHE')THEN
                   IF(VERBOSE)WRITE(LUER,900)'Using BS_HHE: ',LST_SPECIES(ID),LST_WAVE(ID),
	1                LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
                 ELSE
                   IF(VERBOSE)WRITE(LUER,900)'Using LEMKE:  ',LST_SPECIES(ID),LST_WAVE(ID),
	1                LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	         END IF
	         ITR=LST_PID(ID)		!Location in file.
	         CALL RD_BS_STRK_TAB(LST_SPECIES(ID),NU_ZERO,
	1               LST_NL(ID),LST_NUP(ID),PROF_TYPE(LN_PNT),ITR,LU_STK)
	         CALL STRK_BS_HHE(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1               ED_MOD,TEMP_IN,VTURB_IN,ND,
	1               NU_ZERO,AMASS_IN,L_FALSE)
!
	       ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'HeI_IR_STRK')THEN
                 IF(VERBOSE)WRITE(LUER,900)'Using HeI_IR_STRK: ',LST_SPECIES(ID),LST_WAVE(ID),
	1                   LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	         CALL STRK_HEI_IR(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1                   ED_MOD,TEMP_IN,VTURB_IN,
	1                   POP_PROTON,POP_HEPLUS,ND,NU_ZERO,
	1                   LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                   AMASS_IN,PROF_TYPE(LN_PNT),LU_STK)
!
	       ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'DS_STRK')THEN
                 IF(VERBOSE)WRITE(LUER,900)'Using DS_STRK: ',LST_SPECIES(ID),LST_WAVE(ID),
	1              LST_NL(ID),LST_NUP(ID),ML_CUR,ML_ST,ML_END
	         CALL STRK_DS(PROF_STORE(1,1,LOC_INDX),XNU,NF,
	1	         ED_MOD,TEMP_IN,VTURB_IN,
	1                ND,NU_ZERO,
	1                LST_SPECIES(ID),LST_NL(ID),LST_NUP(ID),
	1                AMASS_IN,PROF_TYPE(LN_PNT),LU_STK)

	       ELSE IF(PROF_TYPE(LN_PNT) .EQ. 'HZ_STARK')THEN
!
! Check validity of passed parameters for HI and HeII lines.
!
	         IF(NUP .LE. NL)THEN
	           WRITE(LUER,*)'Error in SET_PROF: NUP < NL'
	           STOP
	         END IF
                 IF(VERBOSE)WRITE(LUER,900)'Using HZ_STARK: ',' ',0.01_LDP*C_KMS/NU_ZERO,
	1              NL,NUP,ML_CUR,ML_ST,ML_END,Z_IN
!
! Convert from Frequency to Angstrom space, measured from line center.
!
	         DO ML=ML_ST,ML_END
	           DWS_GRIEM(ML-ML_ST+1)=0.01_LDP*C_KMS*(1.0_LDP/NU(ML)-1.0_LDP/NU_ZERO)
	         END DO
!
! Now compute, and store the Doppler broadened Stark profile for each
! depth. We assume that the profile is zero outside the computational range.
!
	         PROF_STORE(1:ND,1:NF,LOC_INDX)=0.0_LDP
!$OMP PARALLEL DO PRIVATE(I,PR_GRIEM,TMP_ED)
	         DO I=1,ND
	           TMP_ED=MIN(ED_IN(I),MAX_PROF_ED)
                   CALL GRIEM_V2(PR_GRIEM,DWS_GRIEM,NF,
	1               TMP_ED,TEMP_IN(I),VTURB_IN(I),
	1               NL,NUP,Z_IN,AMASS_IN,RET_LOG)
	            PROF_STORE(I,1:NF,LOC_INDX)=PR_GRIEM(1:NF)
	         END DO
!$OMP END PARALLEL DO
!
	       ELSE IF(.NOT. PROF_COMPUTED)THEN
	         WRITE(LUER,*)'Error in SET_PROF_V3'
	         WRITE(LUER,*)'Unrecognized profile type'
	         WRITE(LUER,*)'PROF_TYPE=',PROF_TYPE(LN_PNT)
	         WRITE(LUER,*)NU_ZERO,NL,NUP,AMASS_IN
	         STOP
	      END IF			!Profile type
!
! Normalize profile to an integral of unity assuming a trapazoidal integral.
! Since NU is units of 10^15 Hz, we need to multiply the integral by
! 10^15 before scaling.
!
	      IF(NORM_PROFILE)THEN
	        NF=ML_END-ML_ST+1
!$OMP PARALLEL DO PRIVATE (I,T1,ML)
	        DO I=1,ND
	          PROF_STORE(I,1,LOC_INDX)=0.0_LDP
	          PROF_STORE(I,NF,LOC_INDX)=0.0_LDP
	          T1=0.0_LDP
	          DO ML=1,NF-1
	              T1=T1+0.5_LDP*(NU(ML_ST+ML-1)-NU(ML_ST+ML))*
	1                (PROF_STORE(I,ML,LOC_INDX)+PROF_STORE(I,ML+1,LOC_INDX))
	          END DO
	          T1=ABS(T1)*1.0E+15_LDP
	          PROF_STORE(I,1:NF,LOC_INDX)=PROF_STORE(I,1:NF,LOC_INDX)/T1
	        END DO
!$OMP END PARALLEL DO
	      END IF
!
! Remember to return profile.
!
	      LINE_PROF_SIM(1:ND,SIM_INDX)=PROF_STORE(1:ND,1,LOC_INDX)
!
! Save information to describe profile so that we don't have to compute it
! all over again.
!
	      NU_ZERO_STORE(LOC_INDX)=NU_ZERO
	      AMASS_STORE(LOC_INDX)=AMASS_IN
	      NL_STORE(LOC_INDX)=NL
	      NUP_STORE(LOC_INDX)=NUP
	      LST_FREQ_LOC(LOC_INDX)=1
	      NF_STORE(LOC_INDX)=NF
	      NU_STORE(1:NF,LOC_INDX)=NU(ML_ST:ML_END)
!
	    END IF
	  ELSE
	    LINE_PROF_SIM(1:ND,SIM_INDX)=0.0_LDP
	  END IF
	END DO
!
!	I=ND-7
!	WRITE(99,*)AMASS_IN,NL,NUP
!	WRITE(99,*)ED_MOD(I),TEMP_IN(I),VTURB_IN(I)
!	WRITE(99,'(1P,5E15.5)')0.01D0*C_KMS/NU_STORE(1:NF,LOC_INDX)
!	WRITE(99,'(1P,5E15.5)')PROF_STORE(I,1:NF,LOC_INDX)
!
	RETURN
	END
