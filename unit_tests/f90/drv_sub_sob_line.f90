!
! Driver for SUB_SOB_LINE_V3 (new_main/mod_subs/sub_sob_line_v3.f).
! Manufactures a minimal 8-depth atmosphere with two simultaneous
! Sobolev lines, sets the module state the routine reads, and prints
! machine-parseable outputs for the Python harness.
!
! The driver always takes the LAMBDA_ITERATION branch (COMPUTE_BA=T,
! LAMBDA_ITERATION=T, VAR_SOB_JC=F, IMPURITY_CODE=F), which exercises the
! altered code without needing the DO_VAR_CONT machinery:
!   - the ETAL_MAT>0/NEG_OPACITY VB_SIM/VC_SIM guard (05-Oct-2026),
!   - the STEQ residual update carrying -DRC_SIM,
!   - the SRAT zero-denominator guard and the expanded dRATE_dUP/dT
!     forms that avoid 0*Inf,
!   - the SIGMA vs SIGMA_TRUE dispatch under DO_NONMONO_SOB /
!     DO_SOB_CLAMPED_LINES,
!   - the SOB_RH_COUPLE invocation under DO_NONLOC_SOB,
!   - the LST_ITERATION SOBEW dispatch.
!
! argv(1) selects the case (see below); outputs print as
!   VB  <i> <vb_sim(i,1)> <vc_sim(i,1)>
!   V2  <i> <vb_2(i)> <vc_2(i)>
!   ZN  <i> <znet_sim(i,1)> <jbar(i)>
!   DR  <i> <drc_sim(i,1)>
!   ST  <k> <steq(1,k)> <steq(2,k)> <steq_t(k)>
!   BAU <k> <ba(2,4,k,k)> <ba(2,3,k,k)> <ba(1,4,k,k)> <ba(2,4,k,k)>
!   BAT <k> <ba_t(4,k,k)> <ba_t(3,k,k)> <ba_t(nt,k,k)>
!   NAN <count of non-finite entries>
!   EW  <ew> <cont_int>                        (LST_ITERATION cases)
!
! Cases:
!   base       healthy inputs, all flags off
!   sig_true   DO_NONMONO_SOB=T, clamped F, SIGMA_TRUE /= SIGMA
!   sig_clmp   DO_NONMONO_SOB=T, clamped T (must equal base bitwise)
!   neg_etal   ETAL_MAT(3,1)<0: VB_SIM/VC_SIM zeroed at that depth
!   neg_opac   NEG_OPACITY(4)=T: VB_SIM zeroed at that depth
!   zeropop    POPS(NUP,K)=0: expanded dRATE_dUP stays finite
!   zerorat    U_STAR_RATIO(:,1)=0: expanded dRATE_dT stays finite
!   nonloc     DO_NONLOC_SOB=T + velocity dip: DRC_SIM /= 0 somewhere
!   ew_true    LST_ITERATION=T, DO_NONMONO_SOB=T, clamped F
!   ew_clmp    LST_ITERATION=T, DO_NONMONO_SOB=T, clamped T
!
      PROGRAM DRV_SUB_SOB
      USE SET_KIND_MODULE
      USE ANG_QW_MOD
      USE CONTROL_VARIABLE_MOD
      USE CMF_SOB_MOD
      USE LINE_MOD
      USE LINE_VEC_MOD
      USE MOD_CMFGEN
      USE OPAC_MOD
      USE RADIATION_MOD
      USE STEQ_DATA_MOD
      USE VAR_RAD_MOD
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: ND=8, NC=2, NP=4, NT=4, NIV=4
      INTEGER, PARAMETER :: NSE=4, MAX_SIM=2, NBND=ND
      INTEGER, PARAMETER :: IDIAG=1, NLF=1
!
      REAL(KIND=LDP) POPS(NT,ND), CHIL(ND), ETAL(ND)
      REAL(KIND=LDP) EW, CONT_INT, FL, CONT_FREQ, AMASS
      LOGICAL NEG_OPACITY(ND)
      CHARACTER*32 CASE
!
      INTEGER I, K, NNAN
      INTEGER ARG_NL, ARG_NUP
      REAL(KIND=LDP) X
      LOGICAL LST_IT
!
      COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
      COMMON/LINE/ OPLIN,EMLIN
      REAL(KIND=LDP) CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
      CASE='base'
      CALL GET_COMMAND_ARGUMENT(1,CASE)
!
! ---- allocate module state ------------------------------------------
!
      ALLOCATE(R(ND)); ALLOCATE(V(ND))
      ALLOCATE(SIGMA(ND)); ALLOCATE(SIGMA_TRUE(ND))
      ALLOCATE(SOURCE(ND)); ALLOCATE(CHI(ND)); ALLOCATE(CHI_SCAT(ND))
      ALLOCATE(RJ(ND)); ALLOCATE(TA(ND))
      ALLOCATE(P(NP)); ALLOCATE(AQW(ND,NP)); ALLOCATE(HQW(ND,NP))
      ALLOCATE(JBAR(ND))
      ALLOCATE(EINA(MAX_SIM)); ALLOCATE(OSCIL(MAX_SIM))
      ALLOCATE(GLDGU(MAX_SIM)); ALLOCATE(FL_SIM(MAX_SIM))
      ALLOCATE(AMASS_SIM(MAX_SIM))
      ALLOCATE(SIM_NL(MAX_SIM)); ALLOCATE(SIM_NUP(MAX_SIM))
      ALLOCATE(SIM_LINE_POINTER(MAX_SIM))
      ALLOCATE(CHIL_MAT(ND,MAX_SIM)); ALLOCATE(ETAL_MAT(ND,MAX_SIM))
      ALLOCATE(BB_COR(ND,MAX_SIM))
      ALLOCATE(VB_SIM(ND,MAX_SIM)); ALLOCATE(VC_SIM(ND,MAX_SIM))
      ALLOCATE(VB_2(ND)); ALLOCATE(VC_2(ND))
      ALLOCATE(BETA(ND)); ALLOCATE(BETAC(ND))
      ALLOCATE(BETAC_SIM(ND,MAX_SIM))
      ALLOCATE(ZNET_SIM(ND,MAX_SIM)); ALLOCATE(DRC_SIM(ND,MAX_SIM))
      ALLOCATE(L_STAR_RATIO(ND,MAX_SIM))
      ALLOCATE(U_STAR_RATIO(ND,MAX_SIM))
      ALLOCATE(dL_RAT_dT(ND,MAX_SIM))
      ALLOCATE(dU_RAT_dT(ND,MAX_SIM))
      ALLOCATE(AVE_ENERGY(NT))
      ALLOCATE(VEC_ID(MAX_SIM)); ALLOCATE(VEC_MNL_F(MAX_SIM))
      ALLOCATE(VEC_MNUP_F(MAX_SIM))
      ALLOCATE(VJ(NT,NBND,ND))
      ALLOCATE(STEQ_T(ND)); ALLOCATE(BA_T(NT,NBND,ND))
      ALLOCATE(ATM(1)%F_TO_S_XzV(4))
      ALLOCATE(SE(1)%STEQ(NSE,ND))
      ALLOCATE(SE(1)%BA(NSE,NIV,NBND,ND))
      ALLOCATE(SE(1)%LNK_TO_F(NIV)); ALLOCATE(SE(1)%LNK_TO_IV(NSE))
!
! ---- healthy baseline state ------------------------------------------
!
      R=(/ 9.0E4_LDP,7.0E4_LDP,5.0E4_LDP,3.5E4_LDP, &
           2.4E4_LDP,1.6E4_LDP,1.1E4_LDP,7.0E3_LDP /)
      V=(/ 8.0E3_LDP,6.5E3_LDP,5.2E3_LDP,4.1E3_LDP, &
           3.2E3_LDP,2.4E3_LDP,1.7E3_LDP,1.0E3_LDP /)
      SIGMA=(/ 0.20_LDP,0.30_LDP,0.40_LDP,0.55_LDP, &
               0.70_LDP,0.85_LDP,1.00_LDP,1.20_LDP /)
      SIGMA_TRUE=SIGMA
      CHI=(/ 1.0E-3_LDP,1.2E-3_LDP,1.5E-3_LDP,1.9E-3_LDP, &
           2.5E-3_LDP,3.4E-3_LDP,4.8E-3_LDP,7.0E-3_LDP /)
      CHI_SCAT=0.3_LDP*CHI
      SOURCE=(/ 2.0E-4_LDP,1.7E-4_LDP,1.4E-4_LDP,1.1E-4_LDP, &
                8.0E-5_LDP,5.0E-5_LDP,2.5E-5_LDP,1.0E-5_LDP /)
      RJ=1.0E-4_LDP
      TA=0.0_LDP
! Core rays (LS<=NC) need P < R(ND); non-core ray LS is tangent at shell
! NI=ND-(LS-NC-1), i.e. P(LS)=R(NI) so Z(NI)=sqrt(R(NI)^2-P^2)=0.
! R(ND)=7E3, R(ND-1)=1.1E4.
      P=(/ 2.0E3_LDP,5.0E3_LDP,7.0E3_LDP,1.1E4_LDP /)
      AQW=0.25_LDP
      HQW=0.25_LDP
!
      POPS=1.0E-2_LDP
      POPS(1,:)=(/ 8.0E-1_LDP,7.5E-1_LDP,6.8E-1_LDP,6.0E-1_LDP, &
                  5.0E-1_LDP,3.8E-1_LDP,2.5E-1_LDP,1.0E-1_LDP /)
      POPS(2,:)=(/ 4.0E-1_LDP,3.7E-1_LDP,3.3E-1_LDP,2.8E-1_LDP, &
                  2.2E-1_LDP,1.5E-1_LDP,9.0E-2_LDP,4.0E-2_LDP /)
      POPS(3,:)=0.7_LDP*POPS(2,:)
      POPS(4,:)=0.2_LDP*POPS(2,:)
      CHIL=(/ 3.0E4_LDP,4.5E4_LDP,6.5E4_LDP,9.5E4_LDP, &
              1.4E5_LDP,2.1E5_LDP,3.2E5_LDP,5.0E5_LDP /)
      ETAL=1.3_LDP*CHIL
      NEG_OPACITY=.FALSE.
!
! Two simultaneous lines.
!
      NUM_SIM_LINES=2
      EINA=(/ 1.0E8_LDP,6.0E7_LDP /)
      OSCIL=(/ 1.0E-1_LDP,5.0E-2_LDP /)
      GLDGU=(/ 1.5_LDP,1.0_LDP /)
      FL_SIM=(/ 4.0E-5_LDP,5.5E-5_LDP /)
      AMASS_SIM=(/ 1.0_LDP,12.0_LDP /)
      SIM_NL=(/ 3,1 /); SIM_NUP=(/ 4,2 /)
      SIM_LINE_POINTER=(/ 1,2 /)
      VEC_ID=(/ 1,1 /); VEC_MNL_F=(/ 1,1 /); VEC_MNUP_F=(/ 2,2 /)
      ATM(1)%F_TO_S_XzV=(/ 1,2,3,4 /)
!
      CHIL_MAT(:,1)=0.6_LDP*CHIL
      CHIL_MAT(:,2)=0.4_LDP*CHIL
      ETAL_MAT(:,1)=0.7_LDP*ETAL
      ETAL_MAT(:,2)=0.5_LDP*ETAL
      BB_COR(:,1)=1.10_LDP
      BB_COR(:,2)=0.90_LDP
!
      L_STAR_RATIO=0.8_LDP
      U_STAR_RATIO=1.2_LDP
      dL_RAT_dT=1.0E-2_LDP
      dU_RAT_dT=2.0E-2_LDP
      AVE_ENERGY=(/ 1.0E-5_LDP,1.5E-5_LDP,3.0E-5_LDP,4.0E-5_LDP /)
!
      VB_SIM=0.0_LDP; VC_SIM=0.0_LDP
      VB_2=0.0_LDP; VC_2=0.0_LDP
      BETA=0.0_LDP; BETAC=0.0_LDP; BETAC_SIM=0.0_LDP
      ZNET_SIM=0.0_LDP; DRC_SIM=0.0_LDP
      JBAR=0.0_LDP
      VJ=0.0_LDP
      STEQ_T=0.0_LDP; BA_T=0.0_LDP
      SE(1)%STEQ=0.0_LDP; SE(1)%BA=0.0_LDP
      SE(1)%LNK_TO_F=(/ 1,2,3,4 /)
      SE(1)%LNK_TO_IV=(/ 1,2,3,4 /)
      SE(1)%N_IV=NIV; SE(1)%N_SE=NSE
!
      CHIBF=1.0_LDP; CHIFF=1.0_LDP; HDKT=4.8E-6_LDP
      TWOHCSQ=0.01474_LDP; OPLIN=3.0E4_LDP; EMLIN=1.0E-5_LDP
!
      COMPUTE_BA=.TRUE.
      LAMBDA_ITERATION=.TRUE.
      DO_NONMONO_SOB=.FALSE.
      DO_NONLOC_SOB=.FALSE.
      DO_SOB_CLAMPED_LINES=.FALSE.
      SCL_LINE_COOL_RATES=.FALSE.
      DIF=.FALSE.
      THK_CONT=.TRUE.
      IC=0.0_LDP
      DBB=0.0_LDP
      METHOD='LOGMON'
      DO_SOB_RH_EXTINCT=.TRUE.
!
! ---- case mutations ---------------------------------------------------
!
      LST_IT=.FALSE.
      SELECT CASE(TRIM(CASE))
      CASE('base')
      CASE('sig_true')
        DO_NONMONO_SOB=.TRUE.
        SIGMA_TRUE(4)=-0.6_LDP
        SIGMA_TRUE(5)=-0.3_LDP
      CASE('sig_clmp')
        DO_NONMONO_SOB=.TRUE.
        DO_SOB_CLAMPED_LINES=.TRUE.
        SIGMA_TRUE(4)=-0.6_LDP
        SIGMA_TRUE(5)=-0.3_LDP
      CASE('neg_etal')
        ETAL_MAT(3,1)=-2.0E4_LDP
      CASE('neg_opac')
        NEG_OPACITY(4)=.TRUE.
      CASE('zeropop')
        POPS(4,:)=0.0_LDP
      CASE('zerorat')
        U_STAR_RATIO(:,1)=0.0_LDP
      CASE('nonloc')
        DO_NONLOC_SOB=.TRUE.
        DO_NONMONO_SOB=.TRUE.
        V(4)=1.2E3_LDP; V(5)=3.4E3_LDP
        SIGMA_TRUE=SIGMA
! SOBEW needs an escaping continuum: with the thick baseline CHI,
! NOES*R exceeds 50 at the outer shell and every T3/T4 factor is zero,
! collapsing EW to 0 regardless of SIGMA.
      CASE('ew_true')
        DO_NONMONO_SOB=.TRUE.
        SIGMA_TRUE(4)=-0.6_LDP
        SIGMA_TRUE(5)=-0.3_LDP
        CHI=CHI*1.0E-4_LDP; CHI_SCAT=CHI_SCAT*1.0E-4_LDP
        LST_IT=.TRUE.
      CASE('ew_clmp')
        DO_NONMONO_SOB=.TRUE.
        DO_SOB_CLAMPED_LINES=.TRUE.
        SIGMA_TRUE(4)=-0.6_LDP
        SIGMA_TRUE(5)=-0.3_LDP
        CHI=CHI*1.0E-4_LDP; CHI_SCAT=CHI_SCAT*1.0E-4_LDP
        LST_IT=.TRUE.
      CASE DEFAULT
        STOP 'bad case'
      END SELECT
!
      FL=4.0E-5_LDP; CONT_FREQ=1.0_LDP; AMASS=1.0_LDP
      EW=0.0_LDP; CONT_INT=0.0_LDP
! NL and NUP are intent(inout) style dummies (the routine overwrites
! them from SIM_NL/SIM_NUP), so pass writable variables not literals.
      ARG_NL=3; ARG_NUP=4
!
      CALL SUB_SOB_LINE_V3('BA',POPS,CHIL,ETAL,NEG_OPACITY, &
          FL,CONT_FREQ,AMASS, &
          .FALSE.,.FALSE.,.FALSE.,LST_IT, &
          EW,CONT_INT, &
          ARG_NL,ARG_NUP,NT,ND,NC,NP, ND,NC,NP, NLF,IDIAG,NBND)
!
! ---- report -----------------------------------------------------------
!
      DO I=1,ND
        WRITE(*,'(A,I3,2ES20.12)')'VB ',I,VB_SIM(I,1),VC_SIM(I,1)
      END DO
      DO I=1,ND
        WRITE(*,'(A,I3,2ES20.12)')'V2 ',I,VB_2(I),VC_2(I)
      END DO
      DO I=1,ND
        WRITE(*,'(A,I3,2ES20.12)')'ZN ',I,ZNET_SIM(I,1),JBAR(I)
      END DO
      DO I=1,ND
        WRITE(*,'(A,I3,ES20.12)')'DR ',I,DRC_SIM(I,1)
      END DO
      DO K=1,ND
        WRITE(*,'(A,I3,3ES20.12)')'ST ',K, &
            SE(1)%STEQ(1,K),SE(1)%STEQ(2,K),STEQ_T(K)
      END DO
      DO K=1,ND
        WRITE(*,'(A,I3,4ES20.12)')'BAU ',K, &
            SE(1)%BA(2,4,K,K),SE(1)%BA(2,3,K,K), &
            SE(1)%BA(1,4,K,K),SE(1)%BA(2,4,K,K)
      END DO
      DO K=1,ND
        WRITE(*,'(A,I3,3ES20.12)')'BAT ',K, &
            BA_T(4,K,K),BA_T(3,K,K),BA_T(NT,K,K)
      END DO
!
      NNAN=0
      DO I=1,ND
        X=SE(1)%STEQ(1,I)+SE(1)%STEQ(2,I)+STEQ_T(I)+VB_SIM(I,1)+ &
          VC_SIM(I,1)+VB_2(I)+VC_2(I)+ZNET_SIM(I,1)+JBAR(I)+ &
          DRC_SIM(I,1)
        IF(.NOT. ABS(X) .LT. 1.0E300_LDP)NNAN=NNAN+1
      END DO
      DO I=1,NSE
        DO K=1,ND
          X=SUM(ABS(SE(1)%BA(I,1:NIV,K,1:ND)))+ &
            SUM(ABS(BA_T(1:NT,K,1:ND)))
          IF(.NOT. ABS(X) .LT. 1.0E300_LDP)NNAN=NNAN+1
        END DO
      END DO
      WRITE(*,'(A,I5)')'NAN ',NNAN
!
      IF(LST_IT)THEN
        WRITE(*,'(A,2ES20.12)')'EW  ',EW,CONT_INT
      END IF
!
      END PROGRAM DRV_SUB_SOB
