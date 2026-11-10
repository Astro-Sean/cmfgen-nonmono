! Driver for FG_J_CMF_V13 (plane/fg_j_cmf_v13.f).
!
! Alteration under test: the diagonal/tridiagonal approximate-Lambda
! store (ALO_DJDETA, ALO_DJDCHI, LAMSTAR_DIAG and the _M/_P off-
! diagonals) accumulated during the INTEGRAL solve - the same
! construction FG_J_CMF_NONMONO_V1 uses for the Knop solver.  It is
! fork-only code (no upstream counterpart), so the meaningful check is
! the stored response against a finite-difference measurement of the
! real solver:
!
!   'base'  - one solve, print JNU/FEDD plus all six ALO arrays.
!   'etaK'  - ETA(K) scaled by 1.01 before the solve.
!   'chiK'  - CHI(K) scaled by 1.01 before the solve.
!
! The test compares (J_pert-J_base)/delta against the stored
! diagonal at the perturbed depth and the _P band at depth K-1.
!
! Atmosphere identical to drv_fg12 (12 depths, 6 core + 12 tangent
! rays, production quadrature weights, METHOD='ZERO', INTEGRAL,
! THK=F, INIT=T, NEW_FREQ=T, INNER_BND_METH='DIFFUSION', DBB=0).
!
      PROGRAM DRV_FG13
      USE SET_KIND_MODULE
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: ND=12, NC=6, NP=ND+NC
      REAL(KIND=LDP) :: ETA(ND),CHI(ND),ESEC(ND),V(ND),SIGMA(ND),R(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: JNU(ND),FEDD(ND),IPLUS_P(NP)
      REAL(KIND=LDP) :: JQW(ND,NP),HQW(ND,NP),KQW(ND,NP),NQW(ND,NP)
      REAL(KIND=LDP) :: HMIDQW(ND,NP),NMIDQW(ND,NP)
      REAL(KIND=LDP) :: VDOP_VEC(ND)
      REAL(KIND=LDP) :: DJDETA(ND),DJDCHI(ND)
      REAL(KIND=LDP) :: DJDETA_M(ND),DJDCHI_M(ND)
      REAL(KIND=LDP) :: DJDETA_P(ND),DJDCHI_P(ND)
      REAL(KIND=LDP) :: RP,DELMU,DBB,IC,FREQ,dLOG_NU
      REAL(KIND=LDP) :: REXT_FAC,VDOP_FRAC
      REAL(KIND=LDP) :: RETURNED_IN_HBC,RETURNED_OUT_HBC
      CHARACTER(LEN=32) :: CASE
      INTEGER :: I,K
      LOGICAL :: MID,ALO_OK,TRI_OK
      EXTERNAL JWEIGHT_V2,HWEIGHT_V2,KWEIGHT_V2,NWEIGHT_V2
!
      CALL GETARG(1,CASE)
!
      DO I=1,ND
        R(I)=20.0_LDP-18.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        V(I)=2000.0_LDP-1800.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        SIGMA(I)=-0.5_LDP
        ETA(I)=1.0E-6_LDP
        CHI(I)=1.0E-8_LDP
        ESEC(I)=1.0E-10_LDP
        VDOP_VEC(I)=30.0_LDP
        JNU(I)=1.0E-8_LDP
        FEDD(I)=1.0_LDP/3.0_LDP
      END DO
      IPLUS_P(:)=0.0_LDP
!
      RP=R(ND)
      P(1)=0.0_LDP
      DELMU=1.0_LDP/NC
      DO I=2,NC
        P(I)=RP*SQRT(1.0_LDP-(DELMU*(NC-I+1))**2)
      END DO
      DO I=NC+1,NP
        P(I)=R(ND-I+NC+1)
      END DO
!
      CALL NORDANGQW_V2(JQW,R,P,NC,ND,NP,JWEIGHT_V2)
      CALL NORDANGQW_V2(HQW,R,P,NC,ND,NP,HWEIGHT_V2)
      CALL NORDANGQW_V2(KQW,R,P,NC,ND,NP,KWEIGHT_V2)
      CALL NORDANGQW_V2(NQW,R,P,NC,ND,NP,NWEIGHT_V2)
      MID=.TRUE.
      CALL GENANGQW_V2(HMIDQW,R,P,NC,ND,NP,HWEIGHT_V2,MID)
      CALL GENANGQW_V2(NMIDQW,R,P,NC,ND,NP,NWEIGHT_V2,MID)
!
      IC=1.0_LDP
      FREQ=1.0_LDP
      dLOG_NU=1.0E-6_LDP
      REXT_FAC=1.0_LDP
      VDOP_FRAC=0.5_LDP
      DBB=0.0_LDP
!
! Perturb one array at one depth: 'eta5' or 'chi9'.  The step is +1%
! so the FD ratio lands next to the stored derivative.
!
      K=0
      IF(CASE(1:3) .EQ. 'eta' .OR. CASE(1:3) .EQ. 'chi') &
         READ(CASE(4:),*)K
      IF(K .GE. 1 .AND. K .LE. ND)THEN
        IF(CASE(1:3) .EQ. 'eta')ETA(K)=ETA(K)*1.10_LDP
        IF(CASE(1:3) .EQ. 'chi')CHI(K)=CHI(K)*1.10_LDP
      END IF
!
      CALL FG_J_CMF_V13(ETA,CHI,ESEC,V,SIGMA,R,P, &
         JNU,FEDD,JQW,HQW,KQW,NQW,HMIDQW,NMIDQW, &
         RETURNED_IN_HBC,RETURNED_OUT_HBC,IPLUS_P, &
         FREQ,dLOG_NU,'DIFFUSION',DBB,IC, &
         VDOP_VEC,VDOP_FRAC,REXT_FAC, &
         'ZERO','INTEGRAL',.FALSE.,.TRUE.,.TRUE.,NC,NP,ND)
!
      DJDETA=0.0_LDP
      DJDCHI=0.0_LDP
      DJDETA_M=0.0_LDP
      DJDCHI_M=0.0_LDP
      DJDETA_P=0.0_LDP
      DJDCHI_P=0.0_LDP
      CALL GET_FG13_ALO_DIAG(DJDCHI,DJDETA,ND,ALO_OK)
      CALL GET_FG13_ALO_TRI(DJDCHI_M,DJDETA_M,DJDCHI_P,DJDETA_P, &
                            ND,TRI_OK)
!
      WRITE(6,'(A,*(ES20.12))')'JNU=',JNU
      WRITE(6,'(A,*(ES20.12))')'FEDD=',FEDD
      WRITE(6,*)'ALO_VALID=',ALO_OK,' TRI_VALID=',TRI_OK
      WRITE(6,'(A,*(ES20.12))')'DJDETA=',DJDETA
      WRITE(6,'(A,*(ES20.12))')'DJDCHI=',DJDCHI
      WRITE(6,'(A,*(ES20.12))')'DJDETA_M=',DJDETA_M
      WRITE(6,'(A,*(ES20.12))')'DJDCHI_M=',DJDCHI_M
      WRITE(6,'(A,*(ES20.12))')'DJDETA_P=',DJDETA_P
      WRITE(6,'(A,*(ES20.12))')'DJDCHI_P=',DJDCHI_P
      WRITE(6,*)'DONE_FG13'
!
      END
