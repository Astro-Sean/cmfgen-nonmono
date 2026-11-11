! Driver for DO_VAR_CONT (new_main/mod_subs/do_var_cont.f).
!
! Alteration under test: the nonmonotonic fixed-J branch (the
! ELSE IF(DO_NONMONO_KNOP .OR. DO_NONMONO_SOB) arm).  Instead of
! VAR_MOM_J_CMF_V12 -- which linearizes the monotonic V11 moment
! system and is invalid on reversal profiles -- the branch copies
! the diagonal approximate-Lambda derivatives stored by whichever
! forward solver ran last:
!   KNOP + USE_LAM_ES=T -> GET_FG_ALO_DIAG   (FG_J_CMF_NONMONO_V1)
!   KNOP + USE_LAM_ES=F -> GET_MOM_ALO_DIAG  (MOM_J_CMF_NONMONO_V1)
!   SOB  + USE_LAM_ES=T -> GET_FG13_ALO_DIAG under DO_FG13_ALO,
!                          GET_FG_ALO_DIAG   under DO_SOB_FG_KNOP
!   SOB  + USE_LAM_ES=F -> GET_MOM_ALO_DIAG
! The stored TA=dJ/dCHI, TB=dJ/dETA land on TX(:,:,1:2) after the
! optional DO_ES_ALI_BA Jacobi amplification 1/(1-TB*CHI_SCAT_CLUMP)
! floored by ES_ALI_BA_MIN, DO_FG13_ALO_TRI fills the K+-1 bands
! from GET_FG13_ALO_TRI (needs NUM_BNDS>=3), and the moment
! quarantine mask (GET_MOM_CELL_QUAR) zeroes frozen rows.  The
! shared accumulator then moves 1:2 into TX(:,:,3:5) and VJ.
!
! The driver injects the stored state through the same SAVE
! variables the solvers write, then reads TX/VJ back out.
!
! Usage: drv CASE
!   nostore : KNOP, no store      -> ALO_OK=F, everything stays 0
!   fg      : KNOP + LAM_ES + FG store   -> TX diag = stored vals
!   mom     : KNOP + mom store           -> moment diagonal used
!   momquar : mom store + QUAR mask      -> masked rows zeroed
!   fg13    : SOB + LAM_ES + FG13 store  -> V13 diagonal used
!   sobfgk  : SOB + LAM_ES + SOB_FG_KNOP -> FG store used
!   sobmom  : SOB + LAM_ES=F + mom store -> moment diagonal used
!   esali   : FG store + DO_ES_ALI_BA    -> Jacobi amplification
!   esalif  : esali with DJDETA straddling the ES_ALI_BA_MIN floor
!             -> DENB capped at 1/MIN, floor-hit counter = 2
!   tri     : SOB + DO_FG13_ALO + TRI store -> off-diag fill
!
      PROGRAM DRV_DVC
      USE SET_KIND_MODULE
      USE MOD_CMFGEN
      USE OPAC_MOD
      USE LINE_MOD
      USE RADIATION_MOD
      USE VAR_RAD_MOD
      USE CONTROL_VARIABLE_MOD
      USE ANG_QW_MOD
      USE FG_J_CMF_NONMONO_MOD_V1, ONLY: &
     &    ALO_DJDCHI,ALO_DJDETA,ALO_VALID,ALO_ND
      USE MOD_MOM_J_NONMONO_V1, ONLY: &
     &    MOM_ALO_DCHI,MOM_ALO_DETA,MOM_ALO_VALID,MOM_ALO_ND, &
     &    MOM_CELL_QUAR
      USE FG_J_CMF_MOD_V13, ONLY: &
     &    V13_DJDCHI => ALO_DJDCHI, V13_DJDETA => ALO_DJDETA, &
     &    V13_DJDCHI_M => ALO_DJDCHI_M, &
     &    V13_DJDETA_M => ALO_DJDETA_M, &
     &    V13_DJDCHI_P => ALO_DJDCHI_P, &
     &    V13_DJDETA_P => ALO_DJDETA_P, &
     &    V13_VALID => ALO_VALID, V13_ND => ALO_ND
      IMPLICIT NONE
      INTEGER, PARAMETER :: ND=4,NT=2,NB=4,DI=2,NM=6
      REAL(KIND=LDP) POPS(NT,ND)
      INTEGER I
      CHARACTER(LEN=16) CASE
      CALL GETARG(1,CASE)
!
      ALLOCATE(R(ND),V(ND),SIGMA(ND),T(ND),ED(ND),ESEC(ND))
      ALLOCATE(TA(ND),TB(ND),TC(ND),RJ(ND),RJ_ES(ND))
      ALLOCATE(CHI_CLUMP(ND),ETA_CLUMP(ND),CHI_SCAT_CLUMP(ND))
      ALLOCATE(ESEC_CLUMP(ND),ES_COH_VEC(ND),ETA_CONT(ND))
      ALLOCATE(CHI_NOSCAT_PREV(ND),ETA_PREV(ND),ETA_C_EVAL(ND))
      ALLOCATE(VCHI(NT,ND),VETA(NT,ND),VCHI_ALL(NT,ND))
      ALLOCATE(VETA_ALL(NT,ND),VCHI_SAV(NT,ND),VETA_SAV(NT,ND))
      ALLOCATE(VCHI_ALL_SAV(NT,ND),VETA_ALL_SAV(NT,ND))
      ALLOCATE(TX(ND,ND,NM),TVX(ND-1,ND,NM),VJ(NT,NB,ND))
      ALLOCATE(dJ_LOC(NM,NB,ND),DIFFW(NT),FA(ND))
      ALLOCATE(dJ_DIF_d_T(ND),dJ_DIF_d_dTdR(ND))
      ALLOCATE(dRSQH_DIF_d_T(ND),dRSQH_DIF_d_dTdR(ND))
      ALLOCATE(DO_THIS_TX_MATRIX(NM),VAR_IN_USE_CNT(NM))
      ALLOCATE(IMP_TRANS_VEC(NM),VAR_LEV_ID(NM))
!
      DO I=1,ND
        R(I)=30.0_LDP-5.0_LDP*I
        V(I)=800.0_LDP-100.0_LDP*I
        SIGMA(I)=-0.5_LDP
        T(I)=3.0_LDP
        ED(I)=1.0_LDP
        ESEC(I)=0.5_LDP
        TA(I)=0.0_LDP; TB(I)=0.0_LDP; TC(I)=0.0_LDP
        RJ(I)=1.0_LDP; RJ_ES(I)=0.8_LDP
        CHI_CLUMP(I)=1.0_LDP; ETA_CLUMP(I)=0.2_LDP
        CHI_SCAT_CLUMP(I)=0.3_LDP; ESEC_CLUMP(I)=0.3_LDP
        ES_COH_VEC(I)=0.0_LDP; ETA_CONT(I)=0.2_LDP
        CHI_NOSCAT_PREV(I)=0.7_LDP; ETA_PREV(I)=0.2_LDP
        ETA_C_EVAL(I)=0.1_LDP
      END DO
      VCHI=0.05_LDP; VETA=0.02_LDP
      VCHI_ALL=0.05_LDP; VETA_ALL=0.02_LDP
      VCHI_SAV=0.05_LDP; VETA_SAV=0.02_LDP
      VCHI_ALL_SAV=0.05_LDP; VETA_ALL_SAV=0.02_LDP
      DIFFW(1:NT)=0.0_LDP; FA(1:ND)=0.0_LDP
      DO_THIS_TX_MATRIX(:)=.TRUE.
      VAR_IN_USE_CNT(:)=0
      VAR_IN_USE_CNT(3:6)=1
      IMP_TRANS_VEC(:)=.FALSE.
      VAR_LEV_ID(:)=1
      SPECIES_PRES(1)=.FALSE.
      POPS=1.0_LDP
!
      CONT_VEL=.TRUE.; ACCURATE=.FALSE.; THIS_FREQ_EXT=.FALSE.
      COHERENT_ES=.FALSE.; MIXED_ES_VAR=.FALSE.
      PLANE_PARALLEL_NO_V=.FALSE.; PLANE_PARALLEL=.FALSE.
      USE_J_REL=.FALSE.; USE_DJDT_RTE=.FALSE.; DIF=.FALSE.
      DO_CLUMP_MODEL=.FALSE.; DO_SRCE_VAR_ONLY=.FALSE.
      INCLUDE_dSLdT=.FALSE.; COMPUTE_NEW_CROSS=.FALSE.
      DO_ES_ALI_BA=(CASE .EQ. 'esali' .OR. CASE .EQ. 'esalif')
      ES_ALI_BA_MIN=0.05_LDP
      ES_ALI_BA_NFLOOR=0
      DO_FG13_ALO_TRI=(CASE .EQ. 'tri')
      DO_FG13_ALO=(CASE .EQ. 'fg13' .OR. CASE .EQ. 'tri')
      DO_SOB_FG_KNOP=(CASE .EQ. 'sobfgk')
!
      DO_NONMONO_SOB=(CASE .EQ. 'fg13' .OR. CASE .EQ. 'sobfgk' &
     &     .OR. CASE .EQ. 'sobmom' .OR. CASE .EQ. 'tri')
      DO_NONMONO_KNOP=.NOT. DO_NONMONO_SOB
      USE_LAM_ES=(CASE .EQ. 'fg' .OR. CASE .EQ. 'fg13' &
     &     .OR. CASE .EQ. 'sobfgk' .OR. CASE .EQ. 'tri' &
     &     .OR. CASE .EQ. 'esali' .OR. CASE .EQ. 'esalif')
!
      IF(CASE .EQ. 'fg' .OR. CASE .EQ. 'esali' &
     &     .OR. CASE .EQ. 'sobfgk' .OR. CASE .EQ. 'esalif')THEN
        ALLOCATE(ALO_DJDCHI(ND),ALO_DJDETA(ND))
        ALO_DJDCHI=(/ -1.0E-02_LDP,-2.0E-02_LDP, &
     &                -3.0E-02_LDP,-4.0E-02_LDP /)
        IF(CASE .EQ. 'esalif')THEN
          ALO_DJDETA=(/ 0.4_LDP,3.1666667_LDP,4.0_LDP,0.7_LDP /)
        ELSE
          ALO_DJDETA=(/ 0.4_LDP,0.5_LDP,0.6_LDP,0.7_LDP /)
        END IF
        ALO_VALID=.TRUE.; ALO_ND=ND
      END IF
      IF(CASE .EQ. 'mom' .OR. CASE .EQ. 'momquar' &
     &     .OR. CASE .EQ. 'sobmom')THEN
        ALLOCATE(MOM_ALO_DCHI(ND),MOM_ALO_DETA(ND))
        MOM_ALO_DCHI=(/ -5.0E-02_LDP,-6.0E-02_LDP, &
     &                  -7.0E-02_LDP,-8.0E-02_LDP /)
        MOM_ALO_DETA=(/ 0.9_LDP,0.8_LDP,0.7_LDP,0.6_LDP /)
        MOM_ALO_VALID=.TRUE.; MOM_ALO_ND=ND
      END IF
      IF(CASE .EQ. 'momquar')THEN
        ALLOCATE(MOM_CELL_QUAR(ND))
        MOM_CELL_QUAR=(/ .FALSE.,.TRUE.,.FALSE.,.TRUE. /)
      END IF
      IF(CASE .EQ. 'fg13' .OR. CASE .EQ. 'tri')THEN
        ALLOCATE(V13_DJDCHI(ND),V13_DJDETA(ND))
        ALLOCATE(V13_DJDCHI_M(ND),V13_DJDETA_M(ND))
        ALLOCATE(V13_DJDCHI_P(ND),V13_DJDETA_P(ND))
        V13_DJDCHI=(/ -1.1E-02_LDP,-2.2E-02_LDP, &
     &                -3.3E-02_LDP,-4.4E-02_LDP /)
        V13_DJDETA=(/ 0.41_LDP,0.51_LDP,0.61_LDP,0.71_LDP /)
        V13_DJDCHI_M=(/ -1.0E-03_LDP,-1.1E-03_LDP, &
     &                  -1.2E-03_LDP,-1.3E-03_LDP /)
        V13_DJDETA_M=(/ 0.11_LDP,0.12_LDP,0.13_LDP,0.14_LDP /)
        V13_DJDCHI_P=(/ -2.0E-03_LDP,-2.1E-03_LDP, &
     &                  -2.2E-03_LDP,-2.3E-03_LDP /)
        V13_DJDETA_P=(/ 0.21_LDP,0.22_LDP,0.23_LDP,0.24_LDP /)
        V13_VALID=.TRUE.; V13_ND=ND
      END IF
!
      CALL DO_VAR_CONT(POPS,'CONTINUUM',.FALSE., &
     &     1.0_LDP,1.0_LDP,1,.TRUE.,2, &
     &     ND,ND,ND,NB,DI,NT,NM, &
     &     ND,ND,ND,0,NM)
!
      WRITE(*,'(A,4ES14.6)')'TX3 ',(TX(I,I,3),I=1,ND)
      WRITE(*,'(A,4ES14.6)')'TX4 ',(TX(I,I,4),I=1,ND)
      WRITE(*,'(A,4ES14.6)')'TX5 ',(TX(I,I,5),I=1,ND)
      WRITE(*,'(A,3ES14.6)')'TX3M ',(TX(I+1,I,3),I=1,ND-1)
      WRITE(*,'(A,3ES14.6)')'TX3P ',(TX(I,I+1,3),I=1,ND-1)
      WRITE(*,'(A,4ES14.6)')'TX4M ',(TX(I+1,I,4),I=1,ND-1),0.0_LDP
      WRITE(*,'(A,4ES14.6)')'TX4P ',0.0_LDP,(TX(I,I+1,4),I=1,ND-1)
      WRITE(*,'(A,4ES14.6)')'VJD ',(VJ(1,I,I),I=1,ND)
      WRITE(*,'(A,I10)')'NFLOOR ',ES_ALI_BA_NFLOOR
      WRITE(*,'(A)')'DONE_'
      END
