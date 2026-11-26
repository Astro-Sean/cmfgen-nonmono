! Driver for the ordering-stall reproduction (Knop remediation WP-0).
!
! Manufactures a ray whose compressing-cell donor weight sits on the
! deeper, unmarched row -- the stale-donor geometry the m15 stall
! diagnosis assigns to SWEEP_FLD_ORD. Every interior row I carries
! RAY_PNT(I)=I and R_RAY between depth nodes I and I+1 at WINX=0.95,
! so the ordered march's compressing donor
!     (1-WINX)*FLD(JPR1,L+NFM) + WINX*FLD(JPR2,L+NFM)
! reads JPR1=J (own row, fresh on the row's descending L pass) and
! JPR2=J+1 (deeper row, a full pass stale on the ascending J march).
! With the near-unity |Q| share parked on the unmarched row the
! refresh solve contracts at ~WINX*|Q| per pass and hits the pass
! cap far above KN_RAY_TOL -- the production signature (stalled rays
! exit IXIT=3/4/6 at residual ~3e-2). SWEEP_FLD computes the whole
! donor column inside one pass, so the same operator is block lower
! triangular for this uniform-sign ray and converges exactly.
!
! Grid: NI=ND=140, J_PNT(K)=K, RAY_PNT(J)=J, KN_R(K)=141-K,
! R_RAY(J)=KN_R(J)+0.95*(KN_R(J+1)-KN_R(J)) for J<NI (uniform Z step
! 1.0, uniform WINX 0.95); the last row self-donates (K+1 off-grid)
! and anchors the chain. Depths 2..140 carry the compressing sparse
! slots so SEED_FLD seeds rows 2..NI from the store.
!
! Cases:
!   ord   ordered sweep (FG_ORD_SWEEP=T), production exits armed,
!         matched spans (NFM=3): stall signature.
!   uno   same physics, FG_ORD_SWEEP=F: clean exit in a few passes.
!   ordx/unox: all early exits off (stall/plateau/bulk/hot), the
!         clean 60-pass contraction measure.
!   ordb/unob: unmatched spans (FG_FREQ_MATCH=F) with drift
!         |GAM|*dz/DLN=30 -- the thin-cell, drift-dominated algebra
!         (chi_hat*dz~30, |Q|/(chi+|Q|)=0.999) with the offset-1
!         donor; NFM=3 under matching is incompatible with a thick
!         cell because matching widens DLNM to NFM*DLN.
!   ordm/unom: mixed-sign ray; arg2 = expanding fraction in tenths
!         (2/5/8). Expanding cells read the same-row L-1 donor, so
!         they stay fresh under either march; compressing runs keep
!         the stale-J+1 defect. unom must converge <=6 passes --
!         the alternating pass direction is what handles both signs.
!   m2u   WP-2 mode 2 (FG_ORD_MODE=2), uniform compressing: the
!         first pass descends so the L+NFM donor column is fresh
!         from the start -- the mode-1 pass 1 is a full stale read.
!         KN_FRESH_DIAG arms and FRESH lines print per pass.
!   m2m   WP-2 mode 2 on the mixed ray (arg2 as ordm/unom): the
!         dominant-sign rule keeps the mode-1 ascending start, so
!         the case must reproduce unom bitwise.
!   capok WP-3 [KN_GATE_COMMIT_RES] on a converging ray: the
!         certified commit stores the fixed-point field, so the
!         post-commit verify reads ~tol and the census is clean.
!   nccm  flag on, ordered march: the stalled exit is not
!         certified, so FG_COMMIT_NC writes a bounded provisional
!         step from the 1e-6 seed -- the store parks far off the
!         fixed point and the committed verify reads >>tol while
!         the solve's own bad census can read clean.
!   nclg  same provisional commit with the flag off: no verify,
!         no census -- the legacy gate sees only the solve's
!         reported metrics.
!
! FG_FRESH_DIAG is armed on ord/uno/m2u/m2m; KN_FRESH_FR(p,d,r)
! prints as FRESH lines after the solve (dir 1=I_M, 2=I_P, -1
! marks a pass that never ran).
!
      PROGRAM DRV_KNOP_ORDERING
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
      IMPLICIT NONE
!
      INTEGER ND,NDX,NCFX,NIX,NCOMPX
      PARAMETER (ND=140,NDX=140,NCFX=200,NIX=140,NCOMPX=139)
      INTEGER I,K,L,LS,IOS,KNP
      INTEGER NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT,NCASE
      INTEGER NFIRST,NLAST,NRUN,NLEN,NLMAX,NEXP,NPAT
      INTEGER*8 NBAD,NTOT,NMATCH,NCOMP,NP2,NBBLK,NTBLK
      INTEGER*8 NCAP_LST,NNG_LST
      REAL(KIND=LDP) RESID,RHO_M,RHO_P
      REAL(KIND=LDP) LGW(NDX)
      REAL(KIND=LDP) GMAG,CHI0,ETA0,DLN0,T1
      CHARACTER(LEN=32) CNAME
      CHARACTER(LEN=32) CARG
!
      CALL GETARG(1,CNAME)
      CALL GETARG(2,CARG)
      NCASE=0
      IF(CNAME .EQ. 'ord')  NCASE=1
      IF(CNAME .EQ. 'uno')  NCASE=2
      IF(CNAME .EQ. 'ordx') NCASE=3
      IF(CNAME .EQ. 'unox') NCASE=4
      IF(CNAME .EQ. 'ordb') NCASE=5
      IF(CNAME .EQ. 'unob') NCASE=6
      IF(CNAME .EQ. 'ordm') NCASE=7
      IF(CNAME .EQ. 'unom') NCASE=8
      IF(CNAME .EQ. 'ordmb')NCASE=9
      IF(CNAME .EQ. 'unomb')NCASE=10
      IF(CNAME .EQ. 'm2u')  NCASE=11
      IF(CNAME .EQ. 'm2m')  NCASE=12
      IF(CNAME .EQ. 'capok')NCASE=13
      IF(CNAME .EQ. 'nccm') NCASE=14
      IF(CNAME .EQ. 'nclg') NCASE=15
      IF(NCASE .EQ. 0)THEN
        WRITE(*,'(A)')'RESULT_FAIL unknown_case'
        GOTO 999
      END IF
      NEXP=0
      IF(NCASE .GE. 7 .AND. NCASE .NE. 11 .AND. NCASE .LT. 13)THEN
        READ(CARG,*,IOSTAT=IOS)NEXP
        IF(IOS .NE. 0 .OR. (NEXP .NE. 2 .AND. NEXP .NE. 5 .AND. &
          NEXP .NE. 8))THEN
          WRITE(*,'(A)')'RESULT_FAIL bad_expanding_fraction'
          GOTO 999
        END IF
      END IF
!
      KN_ND=ND
      KNP=1
      KN_NP=KNP
      KN_NC=0
      KN_NCF=NCFX
      KN_ND_EXT=NDX
      ND_ADD=0
      N_COMP=NCOMPX
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
! Geometry. KN_R descends by 1 per node; every interior row sits 95%
! of the way to the deeper node so WINX=0.95 uniformly and JPR2=J+1
! carries the stale share. Row NI self-donates through IKC1 alone
! (K+1=NI+1 leaves the grid), so it is the converged anchor at the
! bottom of the chain. Z mirrors R_RAY: uniform step 1.0 everywhere.
!
      NI_RAY=NIX
      DO K=1,ND
        KN_R(K)=DBLE(ND+1-K)
      END DO
      DO LS=1,KNP
        DO I=1,NIX
          RAY_PNT(I,LS)=I
          IF(I .LT. NIX)THEN
            R_RAY(I,LS)=KN_R(I)+0.95_LDP*(KN_R(I+1)-KN_R(I))
          ELSE
            R_RAY(I,LS)=R_RAY(I-1,LS)-1.0_LDP
          END IF
          Z(I,LS)=R_RAY(I,LS)
        END DO
      END DO
      R_EXT=KN_R
      LOG_R_EXT=LOG(KN_R)
      KN_P=0.0_LDP
      DO LS=1,KNP
        DO K=1,ND
          J_PNT(K,LS)=K
        END DO
      END DO
      COMP_DEPTH_IDX=0
      DO K=2,ND
        COMP_DEPTH_IDX(K)=K-1
      END DO
!
! Signed GAM. Compressing everywhere (GAM<0) for the uniform cases;
! ordm/unom flip the sign on an evenly spaced subset of interior
! rows (expanding cells always read the same-row L-1 donor and stay
! fresh, so they anchor the compressing runs between them). The
! magnitude sets the drift: |GAM|*dz/DLN_N=2.5 -> matched NFM=3
! (DLNM=0.3, g=0.833) for the matched cases, 30 for the unmatched
! drift-dominated pair (DLNM=DLN_N, g=30). chi=g/999 keeps
! |Q|/(chi+|Q|)=0.999 in both; ETA=chi+g keeps SR=1.
!
      DLN0=0.1_LDP
      IF(NCASE .EQ. 5 .OR. NCASE .EQ. 6 .OR. NCASE .EQ. 9 .OR. &
        NCASE .EQ. 10)THEN
        GMAG=3.0_LDP
      ELSE
        GMAG=0.25_LDP
      END IF
      DO LS=1,KNP
        DO I=1,NIX
          T1=-1.0_LDP
          IF(NCASE .GE. 7 .AND. NCASE .NE. 11 .AND. &
          NCASE .LT. 13 .AND. &
          I .GE. 2 .AND. I .LE. NIX-1)THEN
            NPAT=10/NEXP
            IF(NEXP .EQ. 8)THEN
              IF(MOD(I-2,5) .NE. 0)T1=1.0_LDP
            ELSE
              IF(MOD(I-2,NPAT) .EQ. 0)T1=1.0_LDP
            END IF
          END IF
          GAM(I,LS)=T1*GMAG
          GAM_MAG(I,LS)=ABS(GAM(I,LS))
          dGAMdR(I,LS)=0.0_LDP
        END DO
      END DO
!
! Opacity/emissivity: constant cubic (chi under matching is
! chi_hat - |GAM|/DLNM = g/999 for |Q|=0.999; ETA=chi_hat -> SR=1).
!
      IF(NCASE .EQ. 5 .OR. NCASE .EQ. 6 .OR. NCASE .EQ. 9 .OR. &
        NCASE .EQ. 10)THEN
        CHI0=GMAG/DLN0/999.0_LDP
        ETA0=CHI0+GMAG/DLN0
      ELSE
        CHI0=GMAG/(3.0_LDP*DLN0)/999.0_LDP
        ETA0=CHI0+GMAG/(3.0_LDP*DLN0)
      END IF
      DO L=1,NCFX
        DO K=1,NDX
          KN_CHI_EXT(K,L)=REAL(CHI0,4)
          KN_ETA_EXT(K,L)=REAL(ETA0,4)
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
        KN_FREQ(L)=KN_FREQ(L-1)*EXP(-DLN0)
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
      SOLUTION_METHOD='INTEGRAL'
!
      ALLOCATE(KN_RES_PASS(2*(400+120)),KN_HOT_CNT(NCFX),STAT=IOS)
      IF(IOS .NE. 0)GOTO 999
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
!
! Controller defaults: production exits armed (stall, plateau, hot
! passes, bulk floor). ordx/unox disarm them all for the clean
! 60-pass contraction measure. ordb/unob drop the matched span so
! the drift-dominated thin cell is honest input.
!
      FG_FREQ_MATCH=(NCASE .NE. 5 .AND. NCASE .NE. 6 .AND. &
        NCASE .NE. 9 .AND. NCASE .NE. 10)
      FG_ORD_SWEEP=(NCASE .EQ. 1 .OR. NCASE .EQ. 3 .OR. &
        NCASE .EQ. 5 .OR. NCASE .EQ. 7 .OR. &
        NCASE .EQ. 14 .OR. NCASE .EQ. 15)
! WP-2: cases 11/12 push the concrete order mode; every other
! case keeps the -1 sentinel so the solve resolves it from
! FG_ORD_SWEEP exactly as before. m2m on a mixed ray resolves
! to the mode-1 ascending start, so it must match unom bitwise.
      IF(NCASE .GE. 11 .AND. NCASE .LT. 13)FG_ORD_MODE=2
!
! WP-2 fresh census: arm it on the four diagnostic cases and
! allocate the (pass,dir,ray) store the solve fills. -1 marks a
! pass that never ran.
!
      FG_FRESH_DIAG=(NCASE .EQ. 1 .OR. NCASE .EQ. 2 .OR. &
        NCASE .GE. 11)
      IF(FG_FRESH_DIAG)THEN
        KN_FRESH_MAXP=2*(400+120)
        IF(ALLOCATED(KN_FRESH_FR))DEALLOCATE(KN_FRESH_FR)
        ALLOCATE(KN_FRESH_FR(KN_FRESH_MAXP,2,KNP),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        KN_FRESH_FR=-1.0_LDP
      END IF
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=.FALSE.
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
! WP-3: nclg disarms the committed-state verify entirely -- the
! census slots then stay at -1 and the gate falls back to the
! pre-commit metrics.
      IF(NCASE .EQ. 15)CALL SET_FG_GATE_COMMIT_RES(.FALSE.)
      FG_COMMIT_CONV=.TRUE.
      FG_COMMIT_NC=0.05_LDP
      CALL SET_FG_CERT_CONV(.TRUE.)
      FG_COMMIT_JUMP=50.0_LDP
      FG_DIV_HEAL_CAP=1.0E+6_LDP
      FG_BULK_CONTR=.FALSE.
      FG_ADAPT_OMEGA_V2=.FALSE.
      FG_PLATEAU_DET=.TRUE.
!
      IF(NCASE .EQ. 3 .OR. NCASE .EQ. 4)THEN
        KN_RAY_BULK_TOL=0.0_LDP
        KN_RAY_STALL=1.1_LDP
        FG_PLATEAU_DET=.FALSE.
        KN_RAY_HOTPASSES=0
        KN_HOT_TOTAL_PASSES=0
      END IF
!
! Store seed 1e-6: the solved field is O(1), so the ray starts far
! from its fixed point the way the production model does.
!
      I_P_DEPTH_NEXT=1.0E-6_LDP
      I_M_DEPTH_NEXT=1.0E-6_LDP
! WP-3: the commit-gate cases drop the seeded mark so the
! first-commit waiver lets the provisional or certified write
! through without the jump veto reading the 1e-6 seed delta as
! a runaway.
      FG_RAY_SEEDED=1
      IF(NCASE .GE. 13)FG_RAY_SEEDED=0
!
! WP-3: on a bare solve call the refresh-level allocs never ran,
! so the driver arms the committed-census slots and the residual
! scalar itself. -1 marks a verify that never ran.
!
      IF(NCASE .GE. 13)THEN
        IF(ALLOCATED(FG_NBAD_CM_RAY))DEALLOCATE(FG_NBAD_CM_RAY)
        IF(ALLOCATED(FG_NTOT_CM_RAY))DEALLOCATE(FG_NTOT_CM_RAY)
        ALLOCATE(FG_NBAD_CM_RAY(KNP),FG_NTOT_CM_RAY(KNP),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        FG_NBAD_CM_RAY=-1_8
        FG_NTOT_CM_RAY=-1_8
        FG_RF_COMMIT_RES=-1.0_LDP
        IF(ALLOCATED(FG_RAY_VC))DEALLOCATE(FG_RAY_VC)
        ALLOCATE(FG_RAY_VC(KNP),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        FG_RAY_VC=-1
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
      DO I=1,2*(400+120)
        IF(KN_RES_PASS(I) .LE. 0.0_LDP)EXIT
        WRITE(*,'(A,I4,A,ES16.8)')'PASSRES ',I,' ',KN_RES_PASS(I)
      END DO
!
      IF(FG_FRESH_DIAG .AND. ALLOCATED(KN_FRESH_FR))THEN
        DO I=1,KN_FRESH_MAXP
          IF(KN_FRESH_FR(I,1,1) .LT. -0.5_LDP)EXIT
          WRITE(*,'(A,I4,A,ES16.8)')'FRESH_M ',I,' ', &
        KN_FRESH_FR(I,1,1)
        END DO
        DO I=1,KN_FRESH_MAXP
          IF(KN_FRESH_FR(I,2,1) .LT. -0.5_LDP)EXIT
          WRITE(*,'(A,I4,A,ES16.8)')'FRESH_P ',I,' ', &
        KN_FRESH_FR(I,2,1)
        END DO
      END IF
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
!
! WP-3 committed-state exports: the post-commit verify max-norm,
! the per-ray bad census, and the commit verdict the chain
! reached. -1 marks a verify that never ran (flag off, veto, or
! a ray with no compressing cells).
!
      IF(NCASE .GE. 13)THEN
        WRITE(*,'(A,ES16.8)')'CRES=',FG_RF_COMMIT_RES
        WRITE(*,'(A,I14)')'CM_NB=',FG_NBAD_CM_RAY(1)
        WRITE(*,'(A,I14)')'CM_NT=',FG_NTOT_CM_RAY(1)
        IF(ALLOCATED(FG_RAY_VC))THEN
          WRITE(*,'(A,I4)')'VC=',FG_RAY_VC(1)
        END IF
      END IF
      DO K=1,N_COMP
        IF(K .GT. 4 .AND. K .LT. N_COMP-1)CYCLE
        WRITE(*,'(A,I3,A,4ES16.8)')'STORE_K ',K,' ', &
      I_P_DEPTH_NEXT(K,1,1),I_P_DEPTH_NEXT(K,1,2), &
      I_P_DEPTH_NEXT(K,1,NCFX/2),I_P_DEPTH_NEXT(K,1,NCFX-1)
      END DO
      WRITE(*,'(A)')'RESULT_OK'
999   CONTINUE
      END PROGRAM DRV_KNOP_ORDERING
