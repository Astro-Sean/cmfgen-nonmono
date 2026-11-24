! Driver for the marginal-coupling experiment (K1 audit, Stage-2
! block-GS scoping). Manufactures a ray whose compressing cells sit in
! the |GAM|/DLN >> chi regime, so the lagged-donor term dominates the
! operator the way the m15 cold-start model does (refresh residual
! pinned ~4-6e-2, rho ~ 0.9, bulk-exit masking). One ray, Z=R,
! RAY_PNT=identity, same poking pattern as drv_knop_ray.f90.
!
! Cases:
!   blk   production controller (bulk floor 6e-2, stall on, hot
!         passes armed) -- must reproduce the masking: early bulk
!         exit, hot budget nearly unspent.
!   blkh  same physics, bulk floor off, stall on, hot passes armed --
!         the ray must stall into TRY_HOT; measures whether windowed
!         refinement contracts the hot band to tolerance.
!   blkn  bulk floor off, stall/plateau off, hot off, 400 passes --
!         clean asymptotic contraction measurement.
!   blkx  blkn at 1e-14 tolerance -- fixed-point depth.
!   blko  blk physics under the ordered sweep.
!   blkb  marginality confined to a frequency band (chi steep in L)
!         so TRY_HOT exercises a window, not the whole grid.
!   loop  40 production refreshes on the same ray -- the store
!         fixed-point map with frozen coefficients (the level where
!         the m15 residual oscillates 2.4e-2 <-> 3.2e-2).
!   loopn loop with no early exits -- the map polished to tol each
!         refresh; separates map non-contraction from exit masking.
!
      PROGRAM DRV_KNOP_BLOCK
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
      IMPLICIT NONE
!
      INTEGER ND,NDX,NCFX,NIX
      PARAMETER (ND=8,NDX=8,NCFX=64,NIX=8)
      INTEGER I,K,L,LS,IOS,KNP
      INTEGER NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT,NCASE
      INTEGER NFIRST,NLAST,NRUN,NLEN,NLMAX,NRF,IRF
      INTEGER*8 NBAD,NTOT,NMATCH,NCOMP,NP2,NBBLK,NTBLK
      INTEGER*8 NCAP_LST,NNG_LST
      REAL(KIND=LDP) RESID,RHO_M,RHO_P
      REAL(KIND=LDP) RR(NDX),LGR(NDX)
      REAL(KIND=LDP) LGW(NDX)
      CHARACTER(LEN=32) CNAME
      REAL(KIND=LDP) GAMP,GNEG,T1
!
      CALL GETARG(1,CNAME)
      NCASE=0
      IF(CNAME .EQ. 'blk')  NCASE=1
      IF(CNAME .EQ. 'blkh') NCASE=2
      IF(CNAME .EQ. 'blkn') NCASE=3
      IF(CNAME .EQ. 'blkx') NCASE=4
      IF(CNAME .EQ. 'blko') NCASE=5
      IF(CNAME .EQ. 'blkb') NCASE=6
      IF(CNAME .EQ. 'loop') NCASE=7
      IF(CNAME .EQ. 'loopn')NCASE=8
      IF(CNAME .EQ. 'loopf')NCASE=9
      IF(CNAME .EQ. 'loopg')NCASE=10
      IF(NCASE .EQ. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL unknown_case'
        GOTO 999
      END IF
!
      RR=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
      LGR=LOG(RR)
!
      KN_ND=ND
      KNP=1
      KN_NP=KNP
      KN_NC=0
      KN_NCF=NCFX
      KN_ND_EXT=NDX
      ND_ADD=0
      N_COMP=6
!
      ALLOCATE(NI_RAY(KNP),RAY_PNT(NIX,KNP),R_RAY(NIX,KNP), &
      Z(NIX,KNP), &
      R_EXT(NDX),LOG_R_EXT(NDX),KN_R(ND),KN_P(KNP), &
      J_PNT(ND,KNP),COMP_DEPTH_IDX(ND), &
      GAM(NIX,KNP),GAM_MAG(NIX,KNP),dGAMdR(NIX,KNP), &
      KN_CHI_EXT(NDX,NCFX),KN_ETA_EXT(NDX,NCFX), &
      KN_CHI_COEF(NDX,4,NCFX),KN_ETA_COEF(NDX,4,NCFX), &
      KN_FREQ(NCFX),KN_DLN_P(NCFX),KN_DLN_N(NCFX), &
      KN_IC(NCFX),KN_DBB(NCFX), &
      I_P_DEPTH_NEXT(N_COMP,KNP,NCFX), &
      I_M_DEPTH_NEXT(N_COMP,KNP,NCFX), &
      FG_RAY_SEEDED(KNP),KN_CS_LS(KNP),KN_CS_REA(KNP),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT alloc_fail=1'
        GOTO 999
      END IF
!
      NI_RAY=NIX
      DO LS=1,KNP
        DO I=1,NIX
          RAY_PNT(I,LS)=I
          R_RAY(I,LS)=RR(I)
          Z(I,LS)=RR(I)
        END DO
      END DO
      R_EXT=RR
      LOG_R_EXT=LGR
      KN_R=RR
      KN_P=0.0_LDP
      DO LS=1,KNP
        DO K=1,ND
          J_PNT(K,LS)=K
        END DO
      END DO
      COMP_DEPTH_IDX=0
      DO K=3,8
        COMP_DEPTH_IDX(K)=K-2
      END DO
!
! Signed GAM: mild expansion outside, a deep compressing band over
! depths 3-5. GNEG picks the coupling strength: -3e-2 gives
! |GAM|/DLN ~ 0.1-0.5 versus chi ~ 4e-3-7e-3 on the band, so the
! lagged-donor coefficient Q dominates chi_hat -- the production
! sigma_true ~ -155 regime. Expansion cells keep the small GAMP
! profile so the fixed point is not a trivial wall.
!
! Signed GAM. The marginal regime needs BOTH |GAM|/DLN >> chi
! (lagged-donor weight x/(1+x) -> 1 in the compressing cells) AND a
! donor offset NFM ~ |GAM|*dZ/DLN of a few bins (a real Jacobi chain,
! not a clipped edge read). With dZ ~ 1e4 and DLN = 0.288, GNEG ~
! -3e-4 gives NFM ~ 10-18 inside the 64-bin grid; chi ~ 1e-5 on the
! band then keeps x ~ 30-100. The compression band covers depths
! 3-8 (N_COMP=6) so the donor cells of compressing cells are
! themselves compressing -- a closed marginal chain, the shape the
! m15 model takes (sigma_true < -1 at all depths).
!
      GAMP=4.0E-6_LDP
      GNEG=-3.0E-4_LDP
      DO I=1,NIX
        IF(I .GE. 3)THEN
          GAM(I,1)=GNEG*(1.0_LDP+0.3_LDP*I)
        ELSE
          GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
        END IF
        GAM_MAG(I,1)=ABS(GAM(I,1))
        dGAMdR(I,1)=0.3_LDP*GAM(I,1)/RR(I)
      END DO
!
! Opacity/emissivity. 'blkb' steepens chi in frequency so the
! marginal band confines to low-L bins (chi overtakes |GAM|/DLN
! around L~25-30); other cases use the mild slope that leaves the
! whole grid marginal.
!
      DO L=1,NCFX
        DO K=1,NDX
          IF(NCASE .EQ. 6)THEN
            KN_CHI_EXT(K,L)=REAL(1.0E-5_LDP* &
      (1.0_LDP+2.0_LDP*L)*(RR(1)/RR(K))**2,4)
          ELSE
            KN_CHI_EXT(K,L)=REAL(1.0E-5_LDP* &
      (1.0_LDP+0.05_LDP*L)*(RR(1)/RR(K))**2,4)
          END IF
          KN_ETA_EXT(K,L)=REAL(2.0E-8_LDP* &
      (1.0_LDP+0.02_LDP*L)*(RR(K)/RR(1))**3,4)
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
      ALLOCATE(KN_RES_PASS(2*(400+120)),KN_HOT_CNT(NCFX),STAT=IOS)
      IF(IOS .NE. 0)GOTO 999
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
!
! Controller defaults. FG_FREQ_MATCH=T reproduces the production
! donor-offset path (DO_NONMONO_FREQ_MATCH=T is a KNOP structural
! override).
!
      FG_FREQ_MATCH=.TRUE.
      FG_ORD_SWEEP=(NCASE .EQ. 5)
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      KN_RAY_TOL=2.0E-3_LDP
      KN_RAY_MAXPASS=60
      KN_RAY_BULK_TOL=6.0E-2_LDP
      KN_RAY_OMEGA=1.0_LDP
      KN_RAY_STALL=0.98_LDP
      KN_RAY_HOTPASSES=120
      KN_HOT_TOTAL_PASSES=480
      KN_ADAPT_RELAX=.FALSE.
      INEXT_RELAX=0.5_LDP
      KN_RAY_PACK_MB=2048.0_LDP
      NONMONO_INEXT_CAP=1.0E+20_LDP
      FG_COMMIT_CONV=.TRUE.
      FG_COMMIT_NC=0.05_LDP
      CALL SET_FG_CERT_CONV(.TRUE.)
      FG_COMMIT_JUMP=50.0_LDP
      FG_DIV_HEAL_CAP=1.0E+6_LDP
      FG_BULK_CONTR=.FALSE.
      FG_ADAPT_OMEGA_V2=.FALSE.
      FG_PLATEAU_DET=.TRUE.
!
! blkh: no bulk floor so the ray must grind to the stall exit where
! TRY_HOT fires. blkn: every early exit off, 400 passes for the
! clean contraction measure. blkx adds the deep tolerance. blkb
! shares blkh's exits so the confined band reaches TRY_HOT.
!
      IF(NCASE .EQ. 2 .OR. NCASE .EQ. 6 .OR. NCASE .EQ. 8)THEN
        KN_RAY_BULK_TOL=0.0_LDP
      END IF
      IF(NCASE .EQ. 3 .OR. NCASE .EQ. 4 .OR. NCASE .EQ. 8)THEN
        KN_RAY_STALL=1.1_LDP
        FG_PLATEAU_DET=.FALSE.
        KN_RAY_MAXPASS=400
        KN_RAY_HOTPASSES=0
        KN_HOT_TOTAL_PASSES=0
      END IF
      IF(NCASE .EQ. 4)KN_RAY_TOL=1.0E-14_LDP
      IF(NCASE .EQ. 6)KN_RAY_MAXPASS=200
!
      I_P_DEPTH_NEXT=0.0_LDP
      I_M_DEPTH_NEXT=0.0_LDP
      FG_RAY_SEEDED=1
      DO LS=1,KNP
        DO K=1,N_COMP
          DO L=1,NCFX-1
            I_M_DEPTH_NEXT(K,LS,L)= &
      1.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
            I_P_DEPTH_NEXT(K,LS,L)= &
      2.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
          END DO
        END DO
      END DO
!
! loopf/loopg seed the store far off the solved field: the smooth
! seed above sits within one pass of the fixed point, so the inner
! solve always reaches tol before the bulk exit can engage. The
! production runs arrive at each refresh carrying a store that is
! O(1) off -- the bulk floor can then fire first, which is exactly
! the plateau mechanism under test.
!
      IF(NCASE .EQ. 9)THEN
        I_M_DEPTH_NEXT=I_M_DEPTH_NEXT*30.0_LDP
        I_P_DEPTH_NEXT=I_P_DEPTH_NEXT*30.0_LDP
      END IF
      IF(NCASE .EQ. 10)THEN
        DO K=1,N_COMP
          DO L=1,NCFX-1
            T1=1.0E-1_LDP*(1.0_LDP+0.05_LDP*L)* &
      (1.0_LDP-2.0_LDP*MOD(K+L,2))
            I_M_DEPTH_NEXT(K,1,L)=T1
            I_P_DEPTH_NEXT(K,1,L)=T1
          END DO
        END DO
      END IF
!
      IF(ALLOCATED(JMP_SIB_M))DEALLOCATE(JMP_SIB_M)
      IF(ALLOCATED(JMP_SIB_P))DEALLOCATE(JMP_SIB_P)
      ALLOCATE(JMP_SIB_M(KN_NCF),JMP_SIB_P(KN_NCF),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL sibref_alloc'
        GOTO 999
      END IF
      CALL KN_JMP_SIBREF_V1(I_M_DEPTH_NEXT,JMP_SIB_M,IOS)
      CALL KN_JMP_SIBREF_V1(I_P_DEPTH_NEXT,JMP_SIB_P,IOS)
!
! loop/loopn: iterate the production refresh on the same ray with
! frozen coefficients -- the pure store fixed-point map. KN_RF_RESID
! is the refresh residual the OUTGEN table reports; the per-refresh
! trace shows whether the map contracts, plateaus, or oscillates.
! KN_Q_RESID arms the independent fixed-point residual read on the
! committed field (the qres diagnostic), so the true outer-map
! residual prints alongside the solve-reported one.
!
      IF(NCASE .GE. 7)THEN
        DO_KNOP_REFRESH=.TRUE.
        FREEZE_I_NEXT=.FALSE.
        HAS_COMPRESSING=.TRUE.
        KN_STORE_READY=.TRUE.
        FREQ_IDX_ITER=KN_NCF-1
        FG_GS_NP=KN_NP
        FG_GS_NCF=KN_NCF
        SOLUTION_METHOD='INTEGRAL'
        KN_Q_RESID=.TRUE.
        NRF=40
        DO IRF=1,NRF
          CALL FG_KNOP_REFRESH_V1()
          WRITE(*,'(A,I3,6ES16.8)')'RF ',IRF,KN_RF_RESID, &
      FG_RF_DCHG,FG_GS_MAXCHG,DBLE(FG_GS_NSTAL), &
      DBLE(FG_GS_NPLAT),INEXT_RELAX
          WRITE(*,'(A,I3,A,I10,A,I10)')'RF2 ',IRF,' ', &
      FG_GS_NBAD,' ',FG_GS_NTOT
          IF(ALLOCATED(FG_RAY_VC))THEN
            WRITE(*,'(A,I3,A,I3,A,ES16.8)')'VC ',IRF,' ', &
      FG_RAY_VC(1),' ',FG_RAY_DCH(1)
          END IF
        END DO
        WRITE(*,'(A,I8)')'HOT_SPENT=',INT(KN_HOT_SPENT,4)
        DO K=1,N_COMP
          WRITE(*,'(A,I2,A,4ES16.8)')'STORE_K ',K,' ', &
      I_P_DEPTH_NEXT(K,1,1),I_P_DEPTH_NEXT(K,1,2), &
      I_P_DEPTH_NEXT(K,1,NCFX/2),I_P_DEPTH_NEXT(K,1,NCFX-1)
        END DO
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      CALL KNOP_RAY_SOLVE_V1(1,RESID,NSTORE,NPASS,NCAP,NNEG, &
      NBAD,NTOT,NDIV,NMATCH,NCOMP,NMAX,NXIT,NP2,NBBLK,NTBLK, &
      RHO_M,RHO_P,NCAP_LST,NNG_LST)
!
      WRITE(*,'(A,ES16.8)')'RESID=',RESID
      WRITE(*,'(A,ES16.8)')'RHO_M=',RHO_M
      WRITE(*,'(A,ES16.8)')'RHO_P=',RHO_P
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
      WRITE(*,'(A,I10)')'NCAP_LST=',NCAP_LST
      WRITE(*,'(A,I10)')'NNG_LST=',NNG_LST
      WRITE(*,'(A,I8)')'HOT_SPENT=',INT(KN_HOT_SPENT,4)
!
! Per-pass residual history: the contraction curve the stall and
! plateau detectors read. Zero-padded tail trimmed on print.
!
      DO I=1,2*(400+120)
        IF(KN_RES_PASS(I) .LE. 0.0_LDP)EXIT
        WRITE(*,'(A,I4,A,ES16.8)')'PASSRES ',I,' ',KN_RES_PASS(I)
      END DO
!
! Hot-run census: contiguous stretches of KN_HOT_CNT>0 are the
! production candidate blocks. Report run count, longest run, and
! the first/last hot index so 'blkb' confinement is readable.
!
      IF(ALLOCATED(KN_HOT_CNT))THEN
        L=COUNT(KN_HOT_CNT(1:KN_NCF) .GT. 0)
        WRITE(*,'(A,I8)')'HOT_CELLFREQ=',L
        NFIRST=0
        NLAST=0
        DO I=1,KN_NCF
          IF(KN_HOT_CNT(I) .GT. 0)THEN
            IF(NFIRST .EQ. 0)NFIRST=I
            NLAST=I
          END IF
        END DO
        WRITE(*,'(A,I6,A,I6)')'HOT_RANGE=',NFIRST,' ',NLAST
        NRUN=0
        NLEN=0
        NLMAX=0
        DO I=1,KN_NCF
          IF(KN_HOT_CNT(I) .GT. 0)THEN
            IF(NLEN .EQ. 0)NRUN=NRUN+1
            NLEN=NLEN+1
            NLMAX=MAX(NLMAX,NLEN)
          ELSE
            NLEN=0
          END IF
        END DO
        WRITE(*,'(A,I6,A,I6)')'HOT_RUNS=',NRUN,' ',NLMAX
      END IF
      IF(ALLOCATED(KN_CS_LS))THEN
        WRITE(*,'(A,I3,A,I3)')'CLS_R=',KN_CS_LS(1), &
      ' ',KN_CS_REA(1)
      END IF
      DO K=1,N_COMP
        WRITE(*,'(A,I2,A,4ES16.8)')'STORE_K ',K,' ', &
      I_P_DEPTH_NEXT(K,1,1),I_P_DEPTH_NEXT(K,1,2), &
      I_P_DEPTH_NEXT(K,1,NCFX/2),I_P_DEPTH_NEXT(K,1,NCFX-1)
      END DO
      WRITE(*,'(A)')'RESULT_OK'
999   CONTINUE
      END PROGRAM DRV_KNOP_BLOCK
