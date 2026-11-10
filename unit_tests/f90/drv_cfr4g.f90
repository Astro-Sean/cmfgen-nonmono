! Driver for CMF_FORMAL_REL_V4_GAM (plane/cmf_formal_rel_v4_gam.f).
!
! Same coastal guards as CMF_FORMAL_REL_V4 under test:
!   1. SQRT(MAX(0,(R-P)*(R+P))) mu diagnostics at both boundaries.
!   2. N_STORE=MAX(1,NINT(2*V(ND)/5)) under INNER_BND_METH='HOLLOW':
!      a coasting V(ND) truncates the upstream value to 0 and the
!      frequency wrap then dies in MOD(CUR_LOC+1,N_STORE).
!
! GAM extras: GAM_ETA(ND,NA) is the per-ray emissivity (column per
! angle/ray), GAM_INT(ND,NA) the returned intensity, ANG_MULT the
! angular refinement multiplier.  VDOP_VEC is set large so the
! velocity-field insertion never fires: SOLVE_CMF_FORMAL_V3_GAM
! hard-stops unless each ray keeps exactly its small-grid point
! count (NRAY == ND or NP-IP+1).
!
! Cases match drv_cfr4:
!   ctl   - HOLLOW boundary, normal expanding V; identical path.
!   coast - V(ND)=0.1 so NINT(2*V(ND)/5)=0.
!   nudge - outermost tangent ray's P pushed 1e-7 above R(1).
!
      PROGRAM DRV_CFR4G
      USE SET_KIND_MODULE
      USE MOD_RD_GAMRAY_CNTRL_VARIABLES, ONLY: VERBOSE_GAMMA
      IMPLICIT NONE
!
! NA must cover the mirrored index J=2*(NP+1-ID)-IP the solver uses
! to pair each ray with its oppositely directed twin - worst case
! J=2*NP-1 at ID=1,IP=1.
!
      INTEGER, PARAMETER :: ND=12, NC=6, NP=ND+NC, NA=2*NP
      REAL(KIND=LDP) :: ETA(ND),CHI(ND),ESEC(ND),V(ND),SIGMA(ND),R(ND)
      REAL(KIND=LDP) :: P(NP)
      REAL(KIND=LDP) :: JNU(ND),FEDD(ND),IPLUS(NP)
      REAL(KIND=LDP) :: VDOP_VEC(ND)
      REAL(KIND=LDP) :: GAM_ETA(ND,NA),GAM_INT(ND,NA)
      REAL(KIND=LDP) :: RP,DELMU,DBB,IC,B_PLANCK,FREQ,dLOG_NU
      REAL(KIND=LDP) :: REXT_FAC,VDOP_FRAC
      REAL(KIND=LDP) :: RET_HNU_AT_IB,RET_HNU_AT_OB
      CHARACTER(LEN=32) :: CNAME
      INTEGER :: I,NU_INDX,ANG_MULT
!
      CALL GETARG(1,CNAME)
!
! define_grid_v2_gam fills R_MU unconditionally but only allocates it
! under VERBOSE_GAMMA - left F the mu-grid loop touches unallocated
! memory.  It also writes ./data/GAM_MU_GRID, so the test creates a
! data/ subdir per case.
!
      VERBOSE_GAMMA=.TRUE.
!
      DO I=1,ND
        R(I)=20.0_LDP-18.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        V(I)=2000.0_LDP-1800.0_LDP*(I-1)/REAL(ND-1,KIND=LDP)
        SIGMA(I)=-0.5_LDP
        ETA(I)=1.0E-6_LDP
        CHI(I)=1.0E-8_LDP
        ESEC(I)=1.0E-10_LDP
        VDOP_VEC(I)=2000.0_LDP
        JNU(I)=1.0E-8_LDP
        FEDD(I)=1.0_LDP/3.0_LDP
      END DO
      IPLUS(:)=0.0_LDP
      DO I=1,NA
        GAM_ETA(:,I)=ETA(:)
      END DO
      GAM_INT(:,:)=0.0_LDP
      NU_INDX=1
      ANG_MULT=1
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
      IC=1.0_LDP
      FREQ=1.0_LDP
      dLOG_NU=1.0E-6_LDP
      REXT_FAC=1.0_LDP
      VDOP_FRAC=0.5_LDP
      DBB=0.0_LDP
      B_PLANCK=1.0E-4_LDP
!
      IF(CNAME .EQ. 'coast')THEN
        V(ND)=0.1_LDP
      ELSE IF(CNAME .EQ. 'nudge')THEN
        P(NP)=R(1)*1.0000001_LDP
      END IF
!
      CALL CMF_FORMAL_REL_V4_GAM(ETA,CHI,ESEC,V,SIGMA,R,P, &
         JNU,FEDD,RET_HNU_AT_IB,RET_HNU_AT_OB,IPLUS, &
         FREQ,dLOG_NU,B_PLANCK,DBB, &
         'HOLLOW',.FALSE., &
         VDOP_VEC,VDOP_FRAC,REXT_FAC, &
         'ZERO',.TRUE.,.TRUE.,NC,NP,ND, &
         NA,NU_INDX,GAM_ETA,GAM_INT,ANG_MULT)
!
      WRITE(6,'(A,*(ES14.5))')'JNU=',JNU
      WRITE(6,'(A,*(ES14.5))')'FEDD=',FEDD
      WRITE(6,'(A,2ES14.5)')'HNU=',RET_HNU_AT_IB,RET_HNU_AT_OB
      WRITE(6,*)'DONE_CFR4G'
!
      END
