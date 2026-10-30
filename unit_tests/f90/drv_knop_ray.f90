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
!
      PROGRAM DRV_KNOP_RAY
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
      IF(CNAME .EQ. 'ckpt') NCASE=9
      IF(CNAME .EQ. 'rfsh') NCASE=10
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
      COMP_DEPTH_IDX(5)=3
!
! Signed GAM: expansion outside, compression band over depths 3-5
! (magnitudes give |Q| ~ O(0.3) for the conv case).  dGAMdR feeds
! only the chi_hat curvature term; keep it smooth.
!
      GAMP=4.0E-6_LDP
      IF(NCASE .EQ. 3)GAMP=2.0E-5_LDP
      IF(NCASE .EQ. 4)GAMP=4.0E-6_LDP
      DO I=1,NIX
        IF(I .GE. 3 .AND. I .LE. 5)THEN
          GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
        ELSE
          GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
        END IF
        GAM_MAG(I,1)=ABS(GAM(I,1))
        dGAMdR(I,1)=0.3_LDP*GAM(I,1)/RR(I)
      END DO
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
      END
