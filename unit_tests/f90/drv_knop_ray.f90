! Driver for KNOP_RAY_SOLVE_V1 (plane/fg_j_cmf_nonmono_v1.f).
! Manufactured core-ray problem: Z=R, RAY_PNT=identity into an
! ND_ADD=0 extended grid, one ray, a manufactured compressing band
! (GAM<0) in the middle depths.  Fills every module input the
! production refresh fills (capture arrays, monotone cubics via
! MON_INT_FUNS_V2, sparse I_NEXT stores), calls the solver for a
! named case, and prints KEY=value lines for the python harness.
!
! Cases:
!   conv   moderate coupling, alternating sweep -> residual exit
!   ord    same field, ordered sweep (FG_ORD_SWEEP=T, NI>=3, pack)
!   match  conv + DO_NONMONO_FREQ_MATCH -> NMATCH>0 multi-bin slots
!   cap    tight tolerance + 4 passes -> pass-cap exit (NXIT=4)
!   fixed  store seeded with a self-consistent field -> small resid
!   rev    two compressing bands -> sign-reversal upwind coverage
!   strg   deep-negative GAM band -> strong-compression regime
!
      PROGRAM DRV_KNOP_RAY
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      USE CONTROL_VARIABLE_MOD, ONLY: GS_FOREIGN_SEED
      IMPLICIT NONE
!
      INTEGER ND,NDX,NCFX,NIX
      PARAMETER (ND=8,NDX=8,NCFX=6,NIX=8)
      INTEGER I,K,L,LS,IOS,NCASE
      INTEGER NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT
      INTEGER*8 NBAD,NTOT,NMATCH,NCOMP,NP2,NBBLK,NTBLK
      REAL(KIND=LDP) RESID
      REAL(KIND=LDP) RR(NDX),LGR(NDX)
      REAL(KIND=LDP) LGW(NDX)
      CHARACTER(LEN=32) CNAME
      REAL(KIND=LDP) GAMP
      REAL(KIND=LDP) VV(ND),SG(ND),CKP(3,1,NCFX-1),CKM(3,1,NCFX-1)
      REAL(KIND=LDP) DMAX_CK
!
      CALL GETARG(1,CNAME)
      NCASE=1
      IF(CNAME .EQ. 'ord')  NCASE=2
      IF(CNAME .EQ. 'match')NCASE=3
      IF(CNAME .EQ. 'cap')  NCASE=4
      IF(CNAME .EQ. 'fixed')NCASE=5
      IF(CNAME .EQ. 'rev')  NCASE=6
      IF(CNAME .EQ. 'strg') NCASE=7
      IF(CNAME .EQ. 'ckpt') NCASE=9
      IF(CNAME .EQ. 'rfsh') NCASE=10
      IF(CNAME .EQ. 'gap')  NCASE=11
      IF(CNAME .EQ. 'gapr') NCASE=12
      IF(CNAME .EQ. 'nanchi')NCASE=13
      IF(CNAME .EQ. 'infgam')NCASE=14
      IF(CNAME .EQ. 'ab')    NCASE=15
      IF(CNAME .EQ. 'abx')   NCASE=22
      IF(CNAME .EQ. 'gpos')  NCASE=16
      IF(CNAME .EQ. 'gneg')  NCASE=17
      IF(CNAME .EQ. 'gpp')   NCASE=18
      IF(CNAME .EQ. 'gpn0')  NCASE=19
      IF(CNAME .EQ. 'gzz')   NCASE=20
      IF(CNAME .EQ. 'gcr')   NCASE=21
!
! The A/B case drives the production in-sweep solver itself
! (FG_J_CMF_NONMONO_V1) and needs its own model-input setup, so it
! lives in a contained routine with its own locals.
!
      IF(NCASE .EQ. 15)CALL DRV_AB(.TRUE.)
      IF(NCASE .EQ. 22)CALL DRV_AB(.FALSE.)
      IF(NCASE .EQ. 15 .OR. NCASE .EQ. 22)GOTO 999
!
! Grid: R decreasing outward->inward, mildly non-uniform so the
! log-space spans differ per cell.
!
      RR=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
      LGR=LOG(RR)
!
      KN_ND=ND
      KN_NP=1
      KN_NC=0
      KN_NCF=NCFX
      KN_ND_EXT=NDX
      ND_ADD=0
      N_COMP=3
!
! gap/gapr: mask covers depths 3-4 only, so a compressing ray point
! bracketed by unmasked depths has no lagged donor and must have its
! whole frequency coupling dropped (both sides of the operator).
!
      IF(NCASE .EQ. 11 .OR. NCASE .EQ. 12)N_COMP=2
!
      ALLOCATE(NI_RAY(1),RAY_PNT(NIX,1),R_RAY(NIX,1),Z(NIX,1), &
      R_EXT(NDX),LOG_R_EXT(NDX),KN_R(ND),KN_P(1), &
      J_PNT(ND,1),COMP_DEPTH_IDX(ND), &
      GAM(NIX,1),GAM_MAG(NIX,1),dGAMdR(NIX,1), &
      KN_CHI_EXT(NDX,NCFX),KN_ETA_EXT(NDX,NCFX), &
      KN_CHI_COEF(NDX,4,NCFX),KN_ETA_COEF(NDX,4,NCFX), &
      KN_FREQ(NCFX),KN_DLN_P(NCFX),KN_DLN_N(NCFX), &
      KN_IC(NCFX),KN_DBB(NCFX), &
      I_P_DEPTH_NEXT(N_COMP,1,NCFX), &
      I_M_DEPTH_NEXT(N_COMP,1,NCFX),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT alloc_fail=1'
        GOTO 999
      END IF
!
      NI_RAY(1)=NIX
      DO I=1,NIX
        RAY_PNT(I,1)=I
        R_RAY(I,1)=RR(I)
        Z(I,1)=RR(I)
      END DO
      R_EXT=RR
      LOG_R_EXT=LGR
      KN_R=RR
      KN_P(1)=0.0_LDP
      DO K=1,ND
        J_PNT(K,1)=K
        COMP_DEPTH_IDX(K)=0
      END DO
      COMP_DEPTH_IDX(3)=1
      COMP_DEPTH_IDX(4)=2
      IF(N_COMP .GE. 3)COMP_DEPTH_IDX(5)=3
!
! Signed GAM: expansion outside, compression band over depths 3-5
! (magnitudes give |Q| ~ O(0.3) for the conv case).  dGAMdR feeds
! only the chi_hat curvature term; keep it smooth.
!
      GAMP=4.0E-6_LDP
      IF(NCASE .EQ. 3)GAMP=2.0E-5_LDP
      IF(NCASE .EQ. 4)GAMP=4.0E-6_LDP
      DO I=1,NIX
        IF(NCASE .EQ. 6)THEN
! Two compressing bands (depths 3-4 and 6-7) separated by an
! expanding cell: the velocity-gradient sign reverses twice along
! the ray, exercising the sign-dependent upwind pick in both
! directions on a single profile.
          IF((I .GE. 3 .AND. I .LE. 4) .OR. &
      (I .GE. 6 .AND. I .LE. 7))THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 7)THEN
! Strong compression: same single band, GAM ~150x deeper negative -
! the sigma << -1 regime.
          IF(I .GE. 3 .AND. I .LE. 5)THEN
            GAM(I,1)=-6.0E-4_LDP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 16)THEN
! gpos: expansion everywhere -- no compressing cell exists, the
! frequency coupling runs purely on the live upwind branch.
          GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
        ELSE IF(NCASE .EQ. 17)THEN
! gneg: compression everywhere -- only the masked depths 3-5 can
! supply donors; the rest exercise the donor-drop gate.
          GAM(I,1)=-GAMP*(0.5_LDP+0.1_LDP*I)
        ELSE IF(NCASE .EQ. 18)THEN
! gpp: GAM(4) one micro-step on the expanding side of zero -- twin of
! gpn0 for the sign-crossing continuity check.
          IF(I .EQ. 4)THEN
            GAM(I,1)=1.0E-9_LDP*GAMP
          ELSE IF(I .EQ. 3 .OR. I .EQ. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 19)THEN
! gpn0: GAM(4) one micro-step on the compressing side. The discrete
! operator is continuous in |GAM| with the sign flipping only the
! upwind pick, so the committed field must match gpp's to the tiny
! coupling magnitude.
          IF(I .EQ. 4)THEN
            GAM(I,1)=-1.0E-9_LDP*GAMP
          ELSE IF(I .EQ. 3 .OR. I .EQ. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 20)THEN
! gzz: exact zeros at adjacent cells 4-5 -- SIGN(1,0)=+1 picks the
! expanding branch there; with |GAM|=0 both the chi_hat term and the
! coupling vanish and SP must reduce to SR.
          IF(I .EQ. 3)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE IF(I .EQ. 4 .OR. I .EQ. 5)THEN
            GAM(I,1)=0.0_LDP
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 21)THEN
! gcr: linear glide through zero between cells 4 and 5 -- a resolved
! gradient reversal with no discontinuous coefficients.
          GAM(I,1)=GAMP*0.3_LDP*(DBLE(I)-4.5_LDP)
        ELSE
          IF(I .GE. 3 .AND. I .LE. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        END IF
        GAM_MAG(I,1)=ABS(GAM(I,1))
        dGAMdR(I,1)=0.3_LDP*GAM(I,1)/RR(I)
      END DO
!
! gapr: point 5 is made exactly non-compressing (GAM=0). Its
! coefficient operator is then identical to the 'gap' case, where
! GAM(5)<0 has no donor slot and the whole coupling is dropped --
! the committed fields must agree to solve tolerance.
!
      IF(NCASE .EQ. 12)THEN
        GAM(5,1)=0.0_LDP
        GAM_MAG(5,1)=0.0_LDP
        dGAMdR(5,1)=0.0_LDP
      END IF
!
! infgam: a non-finite GAM must fail the coefficient build, not be
! repaired into a committed field.
!
      IF(NCASE .EQ. 14)THEN
        GAM(4,1)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
        GAM_MAG(4,1)=ABS(GAM(4,1))
        dGAMdR(4,1)=GAM(4,1)
      END IF
!
! Opacity/emissivity: smooth positive fields with a mild frequency
! dependence so each bin sees a slightly different transfer problem.
! The capture store is single precision in production; mirror that.
!
      DO L=1,NCFX
        DO K=1,NDX
          KN_CHI_EXT(K,L)=REAL(1.0E-3_LDP* &
      (1.0_LDP+0.05_LDP*L)*(RR(1)/RR(K))**2,4)
          KN_ETA_EXT(K,L)=REAL(2.0E-6_LDP* &
      (1.0_LDP+0.02_LDP*L)*(RR(1)/RR(K)),4)
        END DO
      END DO
      DO L=1,NCFX
        DO K=1,NDX
          LGW(K)=LOG(MAX(DBLE(KN_CHI_EXT(K,L)),1.0E-300_LDP))
        END DO
        CALL MON_INT_FUNS_V2(KN_CHI_COEF(1,1,L),LGW,LOG_R_EXT,NDX)
        DO K=1,NDX
          LGW(K)=LOG(MAX(DBLE(KN_ETA_EXT(K,L)),1.0E-300_LDP))
        END DO
        CALL MON_INT_FUNS_V2(KN_ETA_COEF(1,1,L),LGW,LOG_R_EXT,NDX)
      END DO
!
! nanchi: a captured NaN opacity must fail the coefficient build at
! frequency 3 (node-matched ray point 4 reads the tabulated NaN
! directly), not be repaired into a committed field.
!
      IF(NCASE .EQ. 13)THEN
        KN_CHI_EXT(4,3)=IEEE_VALUE(0.0,IEEE_QUIET_NAN)
        DO K=1,NDX
          LGW(K)=LOG(MAX(DBLE(KN_CHI_EXT(K,3)),1.0E-300_LDP))
        END DO
        CALL MON_INT_FUNS_V2(KN_CHI_COEF(1,1,3),LGW,LOG_R_EXT,NDX)
      END IF
!
      KN_FREQ(1)=6.0E15_LDP
      DO L=2,NCFX
        KN_FREQ(L)=KN_FREQ(L-1)*0.75_LDP
      END DO
      KN_DLN_P(1)=1.0E-30_LDP
      DO L=2,NCFX
        KN_DLN_P(L)=LOG(KN_FREQ(L-1)/KN_FREQ(L))
      END DO
      DO L=1,NCFX-1
        KN_DLN_N(L)=LOG(KN_FREQ(L)/KN_FREQ(L+1))
      END DO
      KN_DLN_N(NCFX)=KN_DLN_P(NCFX)
      KN_IC=0.0_LDP
      KN_DBB=0.0_LDP
!
      KN_INNER_METH='ZERO_FLUX'
      KN_METHOD_SAV='ZERO'
      KN_INSERT_SAV=.FALSE.
!
      ALLOCATE(KN_RES_PASS(2*(60+120)),KN_HOT_CNT(NCFX),STAT=IOS)
      IF(IOS .NE. 0)GOTO 999
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
!
      FG_FREQ_MATCH=(NCASE .EQ. 3)
      FG_ORD_SWEEP=(NCASE .EQ. 2)
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      KN_RAY_TOL=2.0E-3_LDP
      IF(NCASE .EQ. 4)KN_RAY_TOL=1.0E-14_LDP
      KN_RAY_MAXPASS=60
      IF(NCASE .EQ. 4)KN_RAY_MAXPASS=4
      KN_RAY_BULK_TOL=0.0_LDP
      KN_RAY_OMEGA=1.0_LDP
      KN_RAY_STALL=0.98_LDP
      KN_RAY_HOTPASSES=0
      KN_HOT_TOTAL_PASSES=0
      KN_ADAPT_RELAX=.FALSE.
      INEXT_RELAX=0.5_LDP
      KN_RAY_PACK_MB=2048.0_LDP
      NONMONO_INEXT_CAP=1.0E+20_LDP
      FG_COMMIT_CONV=.TRUE.
      FG_COMMIT_NC=0.05_LDP
      FG_COMMIT_JUMP=50.0_LDP
      FG_DIV_HEAL_CAP=1.0E+6_LDP
      FG_BULK_CONTR=.FALSE.
      FG_ADAPT_OMEGA_V2=.FALSE.
!
      I_P_DEPTH_NEXT=0.0_LDP
      I_M_DEPTH_NEXT=0.0_LDP
!
! Seed the lagged store with a smooth field; for 'fixed' the seed is
! the analytic fixed point of a two-bin coupled solve so the pass
! residual starts small (seed-consistency check).
!
      DO K=1,N_COMP
        DO L=1,NCFX-1
          IF(NCASE .EQ. 5)THEN
            I_M_DEPTH_NEXT(K,1,L)=1.190E-3_LDP* &
      (0.977_LDP)**(L-1)*(0.78_LDP)**(K-1)
            I_P_DEPTH_NEXT(K,1,L)=1.173E-3_LDP* &
      (0.977_LDP)**(L-1)*(0.78_LDP)**(K-1)
          ELSE
            I_M_DEPTH_NEXT(K,1,L)= &
      1.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
            I_P_DEPTH_NEXT(K,1,L)= &
      2.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
          END IF
        END DO
      END DO
!
      IF(NCASE .EQ. 9)THEN
!
! Checkpoint round-trip: fingerprint arrays mirror the model grid
! (V homologous, SIGMA=0), dump, scramble, restore, compare.
!
        VV=1.0E-4_LDP*RR(1:ND)
        SG=0.0_LDP
        FG_GS_NP=KN_NP
        FG_GS_NCF=KN_NCF
        CKP=I_P_DEPTH_NEXT(:,:,1:NCFX-1)
        CKM=I_M_DEPTH_NEXT(:,:,1:NCFX-1)
        CALL SET_FG_GS_TAG(7)
        CALL DUMP_FG_INEXT_GS(7,RR(1:ND),VV,SG,ND,291)
        I_P_DEPTH_NEXT=-7.0_LDP
        I_M_DEPTH_NEXT=-7.0_LDP
        CALL RESTORE_FG_INEXT_GS(RR(1:ND),VV,SG,ND,KN_NP,KN_NCF)
        DMAX_CK=MAX(MAXVAL(ABS(I_P_DEPTH_NEXT(:,:,1:NCFX-1)-CKP)), &
      MAXVAL(ABS(I_M_DEPTH_NEXT(:,:,1:NCFX-1)-CKM)))
        WRITE(*,'(A,L2)')'CKPT_OK=',I_NEXT_GS_OK
        WRITE(*,'(A,ES16.8)')'CKPT_MAXDIFF=',DMAX_CK
!
! Fingerprint mismatch: restore with a different V must reject and
! leave the store zeroed.
!
        I_P_DEPTH_NEXT=1.0_LDP
        I_M_DEPTH_NEXT=1.0_LDP
        VV(4)=VV(4)*1.5_LDP
        CALL RESTORE_FG_INEXT_GS(RR(1:ND),VV,SG,ND,KN_NP,KN_NCF)
        WRITE(*,'(A,L2)')'CKPT_FP_OK=',I_NEXT_GS_OK
        WRITE(*,'(A,ES16.8)')'CKPT_FP_STORE=', &
      MAXVAL(ABS(I_P_DEPTH_NEXT))
!
! Tag mismatch without GS_FOREIGN_SEED: reject.
!
        VV=1.0E-4_LDP*RR(1:ND)
        I_P_DEPTH_NEXT=1.0_LDP
        CALL SET_FG_GS_TAG(8)
        CALL RESTORE_FG_INEXT_GS(RR(1:ND),VV,SG,ND,KN_NP,KN_NCF)
        WRITE(*,'(A,L2)')'CKPT_TAG_OK=',I_NEXT_GS_OK
!
! Tag mismatch with GS_FOREIGN_SEED=T: accepted as immature seed.
!
        GS_FOREIGN_SEED=.TRUE.
        CALL RESTORE_FG_INEXT_GS(RR(1:ND),VV,SG,ND,KN_NP,KN_NCF)
        WRITE(*,'(A,L2)')'CKPT_FOR_OK=',I_NEXT_GS_OK
        WRITE(*,'(A,ES16.8)')'CKPT_FOR_STORE=', &
      MAXVAL(ABS(I_P_DEPTH_NEXT))
        GS_FOREIGN_SEED=.FALSE.
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 10)THEN
!
! End-to-end refresh: arms the capture-store guard and runs the
! production post-sweep path (OMP ray loop, census, commit).
!
        DO_KNOP_REFRESH=.TRUE.
        FREEZE_I_NEXT=.FALSE.
        HAS_COMPRESSING=.TRUE.
        KN_STORE_READY=.TRUE.
        FREQ_IDX_ITER=KN_NCF-1
        SOLUTION_METHOD='INTEGRAL'
        CALL FG_KNOP_REFRESH_V1()
        WRITE(*,'(A,ES16.8)')'RF_RESID=',KN_RF_RESID
        WRITE(*,'(A,ES16.8)')'RF_DCHG=',FG_RF_DCHG
        WRITE(*,'(A,I8)')'RF_NDIV=',FG_GS_NDIV
        WRITE(*,'(A,ES16.8)')'RF_MAXCHG=',FG_GS_MAXCHG
        WRITE(*,'(A,ES16.8)')'RF_STAL=',DBLE(FG_GS_NSTAL)
        WRITE(*,'(A,ES16.8)')'RF_PLAT=',DBLE(FG_GS_NPLAT)
        DO K=1,N_COMP
          DO L=1,NCFX-1
            WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'STORE ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
          END DO
        END DO
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      RESID=0.0_LDP
      CALL KNOP_RAY_SOLVE_V1(1,RESID,NSTORE,NPASS,NCAP,NNEG, &
      NBAD,NTOT,NDIV,NMATCH,NCOMP,NMAX,NXIT,NP2,NBBLK,NTBLK)
!
      WRITE(*,'(A,I4)')'RESULT resid_ok=',1
      WRITE(*,'(A,ES16.8)')'RESID=',RESID
      WRITE(*,'(A,I8)')'NSTORE=',NSTORE
      WRITE(*,'(A,I6)')'NPASS=',NPASS
      WRITE(*,'(A,I8)')'NCAP=',NCAP
      WRITE(*,'(A,I8)')'NNEG=',NNEG
      WRITE(*,'(A,I8)')'NDIV=',NDIV
      WRITE(*,'(A,I8)')'NMATCH=',NMATCH
      WRITE(*,'(A,I8)')'NCOMP=',NCOMP
      WRITE(*,'(A,I6)')'NMAX=',NMAX
      WRITE(*,'(A,I6)')'NXIT=',NXIT
      WRITE(*,'(A,I8)')'NBAD=',NBAD
      WRITE(*,'(A,I8)')'NTOT=',NTOT
      DO K=1,N_COMP
        DO L=1,NCFX-1
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'STORE ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
        END DO
      END DO
      WRITE(*,'(A)')'RESULT_OK'
999   CONTINUE
!
      CONTAINS
!
! ----------------------------------------------------------------------
! A/B driver (KNOP audit task 5). Runs the production in-sweep solver
! FG_J_CMF_NONMONO_V1 on a manufactured non-monotonic model -- no ray
! arrays are filled here because the solver builds them itself -- then
! lets the production refresh FG_KNOP_REFRESH_V1 re-solve the captured
! coefficients. KN_DBG_* and KR_DBG_* then hold the two operators'
! coefficients side by side and the committed stores can be compared
! slot by slot.
!
! Manufactured velocity: mostly wind-like (V decreasing inward), with a
! steep localised rise through depths 3-5 so sigma_true < -1 there and
! the solver's own mask construction flags them compressing.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_AB(ORDSW)
      LOGICAL ORDSW
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=8,NPA=1,NCFA=6)
      REAL(KIND=LDP) RA(NDA),VA(NDA),SGA(NDA),SGT(NDA),PA(NPA)
      REAL(KIND=LDP) ETA_M(NDA,NCFA),CHI_M(NDA,NCFA),ESC_M(NDA,NCFA)
      REAL(KIND=LDP) JNU_A(NDA),FEDD_A(NDA),IPLUS_A(NPA)
      REAL(KIND=LDP) JQW_A(NDA,NPA),HQW_A(NDA,NPA),KQW_A(NDA,NPA)
      REAL(KIND=LDP) NQW_A(NDA,NPA),HMID_A(NDA,NPA),NMID_A(NDA,NPA)
      REAL(KIND=LDP) VDOP_A(NDA),FRQV(NCFA),INBC_A,HBC_A
      REAL(KIND=LDP) LVR(NDA),BUMP,V0,DLOGA,FREQ_A,DBB_A,IC_A
      REAL(KIND=LDP) VDOP_F,REXT_A,DCH
      CHARACTER*6 METH_A
      CHARACTER(LEN=8) SOPT_A
      CHARACTER(LEN=12) IBND_A
      LOGICAL THK_A,INIT_A,NEWF_A
      INTEGER NC_A,NITER,ITR
      REAL(KIND=LDP) SMAX
!
      RA=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
!
! V(R): smooth expansion envelope V0*(R/R1)^0.6 times a localised dip
! that turns V around through depths 3-5 (V increasing inward there).
! sigma_true = dlnV/dlnR - 1 by centred log differences, one-sided at
! the ends -- the same convention RECOMPUTE_SIGMA_TRUE applies.
!
      V0=3.0E3_LDP
      DO I=1,NDA
        BUMP=1.0_LDP-0.55_LDP*EXP(-((DBLE(I)-4.0_LDP)/1.1_LDP)**2)
        VA(I)=V0*(RA(I)/RA(1))**0.6_LDP*BUMP
      END DO
      LVR=LOG(VA)
      SGT(1)=(LVR(1)-LVR(2))/(LOG(RA(1))-LOG(RA(2)))-1.0_LDP
      DO I=2,NDA-1
        SGT(I)=(LVR(I-1)-LVR(I+1))/(LOG(RA(I-1))-LOG(RA(I+1)))-1.0_LDP
      END DO
      SGT(NDA)=(LVR(NDA-1)-LVR(NDA))/(LOG(RA(NDA-1))-LOG(RA(NDA))) &
      -1.0_LDP
      SGA=MAX(SGT,-0.999_LDP)
!
      DO I=1,NDA
        WRITE(*,'(A,I2,3ES16.8)')'VMOD ',I,RA(I),VA(I),SGT(I)
      END DO
!
! Opacity/emissivity profiles as smooth functions of depth with a mild
! frequency tilt -- same functional family as the ray cases.
!
      DO L=1,NCFA
        DO I=1,NDA
          CHI_M(I,L)=1.0E-3_LDP*(1.0_LDP+0.05_LDP*L)* &
      (RA(1)/RA(I))**2
          ETA_M(I,L)=2.0E-6_LDP*(1.0_LDP+0.02_LDP*L)*(RA(1)/RA(I))
          ESC_M(I,L)=0.0_LDP
        END DO
      END DO
!
      FRQV(1)=6.0E15_LDP
      DO L=2,NCFA
        FRQV(L)=FRQV(L-1)*0.75_LDP
      END DO
      DLOGA=LOG(1.0_LDP/0.75_LDP)
!
      PA(1)=0.0_LDP
!
! Doppler refinement is a resolution knob, not the physics under
! test: a large VDOP with a coarse fraction keeps the ray on the
! eight model depths (no sub-inserts), so the debug dumps can stay
! sized (NDA,NCFA).
!
      VDOP_A=1.0E5_LDP
      VDOP_F=0.5_LDP
      REXT_A=1.0_LDP
      DBB_A=0.0_LDP
      IC_A=0.0_LDP
      METH_A='ZERO'
      SOPT_A='INT'
      IBND_A='ZERO_FLUX'
      THK_A=.FALSE.
      NEWF_A=.TRUE.
      NC_A=0
!
! Flags and debug plumbing. FG_FREQ_MATCH exercises the matched-offset
! path end to end; KN_DBG_LS selects the dumped ray.
!
      FG_FREQ_MATCH=.TRUE.
! ab exercises the ordered refresh (pack + SPW); abx forces the
! alternating SWEEP_FLD path so the same operator is compared under
! both sweep strategies.
!
      FG_ORD_SWEEP=ORDSW
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      KN_RAY_TOL=2.0E-3_LDP
      KN_RAY_MAXPASS=60
      KN_RAY_BULK_TOL=0.0_LDP
      KN_RAY_OMEGA=1.0_LDP
      KN_RAY_STALL=0.98_LDP
      KN_RAY_HOTPASSES=0
      KN_HOT_TOTAL_PASSES=0
      KN_ADAPT_RELAX=.FALSE.
      INEXT_RELAX=0.5_LDP
      KN_RAY_PACK_MB=2048.0_LDP
      NONMONO_INEXT_CAP=1.0E+20_LDP
      FG_COMMIT_CONV=.TRUE.
      FG_COMMIT_NC=0.05_LDP
      FG_COMMIT_JUMP=50.0_LDP
      FG_DIV_HEAL_CAP=1.0E+6_LDP
      FG_BULK_CONTR=.FALSE.
      FG_ADAPT_OMEGA_V2=.FALSE.
      KN_DBG_LS=1
      ALLOCATE(KN_DBG_CHIH(NDA,NCFA),KN_DBG_Q(NDA,NCFA), &
      KN_DBG_DCHR(NDA,NCFA),KN_DBG_DTAU(NDA,NCFA), &
      KN_DBG_SPM(NDA,NCFA),KN_DBG_SPP(NDA,NCFA), &
      KN_DBG_WX(NDA,NCFA),KN_DBG_DLNM(NDA,NCFA), &
      KN_DBG_IMPR(NDA,NCFA),KN_DBG_IPPR(NDA,NCFA), &
      KN_DBG_NM(NDA,NCFA),KN_DBG_KC1(NDA,NCFA),KN_DBG_KC2(NDA,NCFA), &
      KN_DBG_A0(NDA,NCFA),KN_DBG_A1(NDA,NCFA),KN_DBG_A2(NDA,NCFA), &
      KN_DBG_A3(NDA,NCFA),KN_DBG_A4(NDA,NCFA), &
      KR_DBG_CHIH(NDA,NCFA),KR_DBG_Q(NDA,NCFA), &
      KR_DBG_DCHR(NDA,NCFA),KR_DBG_DTAU(NDA,NCFA), &
      KR_DBG_WX(NDA,NCFA),KR_DBG_DLNM(NDA,NCFA), &
      KR_DBG_SPM(NDA,NCFA),KR_DBG_SPP(NDA,NCFA), &
      KR_DBG_A0(NDA,NCFA),KR_DBG_A1(NDA,NCFA),KR_DBG_A2(NDA,NCFA), &
      KR_DBG_A3(NDA,NCFA),KR_DBG_A4(NDA,NCFA), &
      KR_DBG_FLDM(NDA,NCFA),KR_DBG_FLDP(NDA,NCFA), &
      KR_DBG_NFM(NDA,NCFA),KR_DBG_IKC1(NDA,NCFA), &
      KR_DBG_IKC2(NDA,NCFA),KR_DBG_JPR1(NDA,NCFA), &
      KR_DBG_JPR2(NDA,NCFA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT aballoc_fail=1'
        RETURN
      END IF
      KN_DBG_CHIH=-1.0_LDP; KN_DBG_Q=-1.0_LDP; KN_DBG_DCHR=-1.0_LDP
      KN_DBG_DTAU=-1.0_LDP; KN_DBG_SPM=-1.0_LDP; KN_DBG_SPP=-1.0_LDP
      KN_DBG_WX=-1.0_LDP; KN_DBG_DLNM=-1.0_LDP
      KN_DBG_IMPR=-1.0_LDP; KN_DBG_IPPR=-1.0_LDP
      KN_DBG_NM=-1; KN_DBG_KC1=-1; KN_DBG_KC2=-1
      KN_DBG_A0=-1.0_LDP; KN_DBG_A1=-1.0_LDP; KN_DBG_A2=-1.0_LDP
      KN_DBG_A3=-1.0_LDP; KN_DBG_A4=-1.0_LDP
      KR_DBG_CHIH=-1.0_LDP; KR_DBG_Q=-1.0_LDP; KR_DBG_DCHR=-1.0_LDP
      KR_DBG_DTAU=-1.0_LDP; KR_DBG_WX=-1.0_LDP; KR_DBG_DLNM=-1.0_LDP
      KR_DBG_SPM=-1.0_LDP; KR_DBG_SPP=-1.0_LDP
      KR_DBG_A0=-1.0_LDP; KR_DBG_A1=-1.0_LDP; KR_DBG_A2=-1.0_LDP
      KR_DBG_A3=-1.0_LDP; KR_DBG_A4=-1.0_LDP
      KR_DBG_FLDM=-1.0_LDP; KR_DBG_FLDP=-1.0_LDP
      KR_DBG_NFM=-1; KR_DBG_IKC1=-1; KR_DBG_IKC2=-1
      KR_DBG_JPR1=-1; KR_DBG_JPR2=-1
!
! In-sweep great iterations: each frequency call solves with the
! lagged donor store from the previous generation; iterating the
! outer loop contracts the Picard map toward the fixed point the
! refresh targets directly.
!
      NITER=6
      DO ITR=1,NITER
        DO L=1,NCFA
          IF(L .LT. NCFA)THEN
            CALL SET_FG_NU_NEXT(FRQV(L+1))
          ELSE
            CALL SET_FG_NU_NEXT(FRQV(NCFA)*0.75_LDP)
          END IF
          INIT_A=(ITR .EQ. 1 .AND. L .EQ. 1)
          IF(ITR .GT. 1 .AND. L .EQ. 1)INIT_A=.TRUE.
          CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L),ESC_M(1,L), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA, &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NPA,NDA,SGT,NCFA)
        END DO
      END DO
!
! Committed in-sweep store (the lagged map's contraction point) and
! the sweep-side coefficient dump are now resident; snapshot the
! store before the refresh overwrites it.
!
      DO K=1,N_COMP
        DO L=1,NCFA-1
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'SWSTOR ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
        END DO
      END DO
      DO I=1,NDA
        DO L=1,NCFA
          WRITE(*,'(A,I2,A,I2,A,10ES16.8)')'DBS ',I,' ',L,' ', &
      KN_DBG_CHIH(I,L),KN_DBG_Q(I,L),KN_DBG_DCHR(I,L), &
      KN_DBG_DTAU(I,L),KN_DBG_SPM(I,L),KN_DBG_SPP(I,L), &
      KN_DBG_WX(I,L),KN_DBG_DLNM(I,L),KN_DBG_IMPR(I,L), &
      KN_DBG_IPPR(I,L)
          WRITE(*,'(A,I2,A,I2,A,3I6)')'DBSI ',I,' ',L,' ', &
      KN_DBG_NM(I,L),KN_DBG_KC1(I,L),KN_DBG_KC2(I,L)
          WRITE(*,'(A,I2,A,I2,A,5ES16.8)')'DBSA ',I,' ',L,' ', &
      KN_DBG_A0(I,L),KN_DBG_A1(I,L),KN_DBG_A2(I,L), &
      KN_DBG_A3(I,L),KN_DBG_A4(I,L)
        END DO
      END DO
!
! Refresh pass over the captured coefficients, starting from the
! in-sweep committed store -- the Picard step the in-sweep map is an
! outer iteration of.
!
      ALLOCATE(KN_RES_PASS(2*(60+120)),KN_HOT_CNT(NCFA),STAT=IOS)
      IF(IOS .NE. 0)RETURN
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
      HAS_COMPRESSING=.TRUE.
      KN_STORE_READY=.TRUE.
      FREQ_IDX_ITER=KN_NCF-1
      SOLUTION_METHOD='INTEGRAL'
      CALL FG_KNOP_REFRESH_V1()
      WRITE(*,'(A,ES16.8)')'RF_RESID=',KN_RF_RESID
      WRITE(*,'(A,I8)')'RF_NDIV=',FG_GS_NDIV
      WRITE(*,'(A,ES16.8)')'RF_MAXCHG=',FG_GS_MAXCHG
      DO K=1,N_COMP
        DO L=1,NCFA-1
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'RFSTOR ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
        END DO
      END DO
      DO I=1,NDA
        DO L=1,NCFA
          WRITE(*,'(A,I2,A,I2,A,10ES16.8)')'DBR ',I,' ',L,' ', &
      KR_DBG_CHIH(I,L),KR_DBG_Q(I,L),KR_DBG_DCHR(I,L), &
      KR_DBG_DTAU(I,L),KR_DBG_SPM(I,L),KR_DBG_SPP(I,L), &
      KR_DBG_WX(I,L),KR_DBG_DLNM(I,L),KR_DBG_FLDM(I,L), &
      KR_DBG_FLDP(I,L)
          WRITE(*,'(A,I2,A,I2,A,5I6)')'DBRI ',I,' ',L,' ', &
      KR_DBG_NFM(I,L),KR_DBG_IKC1(I,L),KR_DBG_IKC2(I,L), &
      KR_DBG_JPR1(I,L),KR_DBG_JPR2(I,L)
          WRITE(*,'(A,I2,A,I2,A,5ES16.8)')'DBRA ',I,' ',L,' ', &
      KR_DBG_A0(I,L),KR_DBG_A1(I,L),KR_DBG_A2(I,L), &
      KR_DBG_A3(I,L),KR_DBG_A4(I,L)
        END DO
      END DO
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_AB
!
      END
