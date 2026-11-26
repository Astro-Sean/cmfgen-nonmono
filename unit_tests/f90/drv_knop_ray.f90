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
!   sibref three-ray poisoned store -> KN_JMP_SIBREF_V1 reference
!   rf3    three identical rays -> per-ray refresh census
!   jmpx   floored ray + Inf sibling -> sibling-ref jump veto
!   prb    ray 3 armed as diagnostic probe -> MODE 1 isolation
!   prbv   probe with 1E6 store -> probe veto stays off counters
!   prbs   probe + KN_SEED_ARM=1 -> arm-independent production
!   punit  MODE 1 helper unit checks (resid pool, omega EXCL)
!
      PROGRAM DRV_KNOP_RAY
      USE, INTRINSIC :: IEEE_ARITHMETIC
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      USE CONTROL_VARIABLE_MOD, ONLY: GS_FOREIGN_SEED,CONV_DIAG_LVL, &
      JFNK_IN_EVAL
      IMPLICIT NONE
!
      INTEGER ND,NDX,NCFX,NIX
      PARAMETER (ND=8,NDX=8,NCFX=6,NIX=8)
      INTEGER I,K,L,LS,IOS,I2,NCASE,KNP,NRF,IRF
      INTEGER NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT
      INTEGER*8 NBAD,NTOT,NMATCH,NCOMP,NP2,NBBLK,NTBLK
      INTEGER*8 NCAP_LST,NNG_LST
      INTEGER*8 NBDQ(1)
      LOGICAL MSKQ(1)
      INTEGER*8 BNC,BNP,BINT,BNEW,BGON,BNCC(2),BNIC(2),BNUN
      REAL(KIND=LDP) BJAC
      REAL(KIND=LDP) RESID,RHO_M,RHO_P
      REAL(KIND=LDP) RR(NDX),LGR(NDX)
      REAL(KIND=LDP) LGW(NDX)
      CHARACTER(LEN=32) CNAME
      REAL(KIND=LDP) GAMP,T1
      REAL(KIND=LDP) VV(ND),SG(ND),CKP(3,1,NCFX-1),CKM(3,1,NCFX-1)
      REAL(KIND=LDP) CHI_D(8),DCH_D(8),Z_D(8),R_D(8),DT_D(8)
      REAL(KIND=LDP) DMAX_CK
      LOGICAL LEX
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
      IF(CNAME .EQ. 'ab2')   NCASE=32
      IF(CNAME .EQ. 'abz')   NCASE=33
      IF(CNAME .EQ. 'gzp1')  NCASE=23
      IF(CNAME .EQ. 'gzn1')  NCASE=24
      IF(CNAME .EQ. 'gzp2')  NCASE=25
      IF(CNAME .EQ. 'gzn2')  NCASE=26
      IF(CNAME .EQ. 'gzp3')  NCASE=27
      IF(CNAME .EQ. 'gzn3')  NCASE=28
      IF(CNAME .EQ. 'gzp4')  NCASE=29
      IF(CNAME .EQ. 'gzn4')  NCASE=30
      IF(CNAME .EQ. 'gzb')   NCASE=31
      IF(CNAME .EQ. 'gpos')  NCASE=16
      IF(CNAME .EQ. 'gposo') NCASE=16
      IF(CNAME .EQ. 'gnego') NCASE=17
      IF(CNAME .EQ. 'gzzo')  NCASE=16
      IF(CNAME .EQ. 'gzzf')  NCASE=16
      IF(CNAME .EQ. 'gneg')  NCASE=17
      IF(CNAME .EQ. 'gpp')   NCASE=18
      IF(CNAME .EQ. 'gpn0')  NCASE=19
      IF(CNAME .EQ. 'gzz')   NCASE=20
      IF(CNAME .EQ. 'gcr')   NCASE=21
      IF(CNAME .EQ. 'slw')   NCASE=34
      IF(CNAME .EQ. 'slwo')  NCASE=35
      IF(CNAME .EQ. 'rfa')   NCASE=36
      IF(CNAME .EQ. 'capx')  NCASE=37
      IF(CNAME .EQ. 'sibref')NCASE=38
      IF(CNAME .EQ. 'rf3')   NCASE=39
      IF(CNAME .EQ. 'jmpx')  NCASE=40
      IF(CNAME .EQ. 'ompos') NCASE=41
      IF(CNAME .EQ. 'omneg') NCASE=42
      IF(CNAME .EQ. 'omone') NCASE=43
      IF(CNAME .EQ. 'omfew') NCASE=44
      IF(CNAME .EQ. 'omout') NCASE=45
      IF(CNAME .EQ. 'jmp2')  NCASE=46
      IF(CNAME .EQ. 'ni0')   NCASE=47
      IF(CNAME .EQ. 'ni1')   NCASE=48
      IF(CNAME .EQ. 'ni1f')  NCASE=97
      IF(CNAME .EQ. 'ni1s')  NCASE=98
      IF(CNAME .EQ. 'heal')  NCASE=49
      IF(CNAME .EQ. 'esav')  NCASE=50
      IF(CNAME .EQ. 'pmon')  NCASE=51
      IF(CNAME .EQ. 'pnom')  NCASE=52
      IF(CNAME .EQ. 'sfail') NCASE=53
      IF(CNAME .EQ. 'rfail') NCASE=54
      IF(CNAME .EQ. 'scold') NCASE=55
      IF(CNAME .EQ. 'dtau')  NCASE=56
      IF(CNAME .EQ. 'p1t1')  NCASE=57
      IF(CNAME .EQ. 'p1t1b') NCASE=58
      IF(CNAME .EQ. 'p1t1bf')NCASE=58
      IF(CNAME .EQ. 'p1t1n') NCASE=58
      IF(CNAME .EQ. 'repl')  NCASE=59
      IF(CNAME .EQ. 'replo') NCASE=60
      IF(CNAME .EQ. 'replc') NCASE=61
      IF(CNAME .EQ. 'replh') NCASE=62
      IF(CNAME .EQ. 'topo')  NCASE=63
      IF(CNAME .EQ. 'prb')   NCASE=64
      IF(CNAME .EQ. 'prbv')  NCASE=65
      IF(CNAME .EQ. 'prbs')  NCASE=66
      IF(CNAME .EQ. 'punit') NCASE=67
      IF(CNAME .EQ. 'p1t1s') NCASE=68
      IF(CNAME .EQ. 'p1t1r') NCASE=72
      IF(CNAME .EQ. 'p1t1j') NCASE=73
      IF(CNAME .EQ. 'p1t1h') NCASE=74
      IF(CNAME .EQ. 'p1t1d') NCASE=75
      IF(CNAME .EQ. 'p1t1c') NCASE=76
      IF(CNAME .EQ. 'p1t1x') NCASE=77
      IF(CNAME .EQ. 'prbz')  NCASE=69
      IF(CNAME .EQ. 'prbl')  NCASE=70
      IF(CNAME .EQ. 'prbr')  NCASE=71
      IF(CNAME .EQ. 'csut')  NCASE=78
      IF(CNAME .EQ. 'cerb')  NCASE=79
      IF(CNAME .EQ. 'cer0')  NCASE=80
      IF(CNAME .EQ. 'cern')  NCASE=81
      IF(CNAME .EQ. 'cer2')  NCASE=82
      IF(CNAME .EQ. 'cerp')  NCASE=83
      IF(CNAME .EQ. 'rfmx')  NCASE=84
      IF(CNAME .EQ. 'rfab')  NCASE=85
      IF(CNAME .EQ. 'invdln')NCASE=86
      IF(CNAME .EQ. 'inval') NCASE=87
      IF(CNAME .EQ. 'badid') NCASE=88
      IF(CNAME .EQ. 'bido')  NCASE=89
      IF(CNAME .EQ. 'lagc')  NCASE=90
      IF(CNAME .EQ. 'lagw')  NCASE=91
      IF(CNAME .EQ. 'lagr')  NCASE=92
      IF(CNAME .EQ. 'lags')  NCASE=93
      IF(CNAME .EQ. 'sdab')  NCASE=94
      IF(CNAME .EQ. 'sdp')   NCASE=95
      IF(CNAME .EQ. 'qres')  NCASE=96
!
! The A/B case drives the production in-sweep solver itself
! (FG_J_CMF_NONMONO_V1) and needs its own model-input setup, so it
! lives in a contained routine with its own locals.
!
      IF(NCASE .EQ. 15)CALL DRV_AB(.TRUE.,1,1)
      IF(NCASE .EQ. 22)CALL DRV_AB(.FALSE.,1,1)
      IF(NCASE .EQ. 32)CALL DRV_AB(.TRUE.,2,2)
      IF(NCASE .EQ. 33)CALL DRV_AB(.TRUE.,1,1)
      IF(NCASE .EQ. 15 .OR. NCASE .EQ. 22 .OR. &
      NCASE .EQ. 32 .OR. NCASE .EQ. 33)GOTO 999
      IF(NCASE .EQ. 57)THEN
        CALL DRV_P1T1
        GOTO 999
      END IF
      IF(NCASE .EQ. 58)THEN
        CALL DRV_P1T1B
        GOTO 999
      END IF
!
! The omega-pool cases exercise the adaptive-relaxation estimator
! directly on synthetic sample sets; no ray grid is needed.
!
      IF(NCASE .GE. 41 .AND. NCASE .LE. 45)THEN
        CALL DRV_OMEGA(NCASE)
        GOTO 999
      END IF
      IF(NCASE .EQ. 67)THEN
        CALL DRV_M1UNIT()
        GOTO 999
      END IF
      IF(NCASE .EQ. 68)THEN
        CALL DRV_P1T1S()
        GOTO 999
      END IF
      IF(NCASE .EQ. 72)THEN
        CALL DRV_P1T1R('all')
        GOTO 999
      END IF
      IF(NCASE .EQ. 73)THEN
        CALL DRV_P1T1R('jac')
        GOTO 999
      END IF
      IF(NCASE .EQ. 74)THEN
        CALL DRV_P1T1R('hist')
        GOTO 999
      END IF
      IF(NCASE .EQ. 75)THEN
        CALL DRV_P1T1R('testd')
        GOTO 999
      END IF
      IF(NCASE .EQ. 76)THEN
        CALL DRV_P1T1R('ctrl')
        GOTO 999
      END IF
      IF(NCASE .EQ. 77)THEN
        CALL DRV_P1T1R('cert')
        GOTO 999
      END IF
      IF(NCASE .EQ. 78)THEN
        CALL DRV_CERT()
        GOTO 999
      END IF
!
! Grid: R decreasing outward->inward, mildly non-uniform so the
! log-space spans differ per cell.
!
      RR=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
      LGR=LOG(RR)
!
      KN_ND=ND
! sibref/rf3/jmpx run a three-ray set so the sibling-reference
! snapshot and the jump veto operate across real sibling columns.
! rfa needs >=KN_RHO_MIN_SAM=10 pooled ratio samples, so it runs
! six rays (12 samples) or the omega estimator takes the small-
! sample retain path it was added to guard.
      KNP=1
      IF(NCASE .EQ. 36)KNP=6
      IF(NCASE .GE. 38 .AND. NCASE .LE. 40)KNP=3
      IF(NCASE .GE. 46 .AND. NCASE .LE. 48)KNP=3
      IF(NCASE .EQ. 97)KNP=3
      IF(NCASE .GE. 59 .AND. NCASE .LE. 62)KNP=3
      IF(NCASE .EQ. 63)KNP=3
      IF(NCASE .GE. 64 .AND. NCASE .LE. 66)KNP=3
      IF(NCASE .GE. 79 .AND. NCASE .LE. 87)KNP=3
      KN_NP=KNP
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
!
! ni0/ni1: ray 1 is a degenerate ray. NI=0 is a construction failure
! (the refresh must fail closed on it); NI=1 is a legitimate
! single-point geometry with no marchable transfer domain and stays
! in the plain unsolved census.
!
      IF(NCASE .EQ. 47)NI_RAY(1)=0
      IF(NCASE .EQ. 48 .OR. NCASE .EQ. 97 .OR. NCASE .EQ. 98) &
        NI_RAY(1)=1
! inval (87): an absurd NI_RAY forces the ray's working-array
! ALLOCATE to fail inside the solve -- the ray must classify
! INVALID_STATE (allocation failure), not divergence.
      IF(NCASE .EQ. 87)NI_RAY(1)=2000000000
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
      IF(NCASE .EQ. 88 .OR. NCASE .EQ. 89)GAMP=4.0E-6_LDP
      DO I=1,NIX
        IF(NCASE .EQ. 6 .OR. NCASE .EQ. 92)THEN
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
        ELSE IF(NCASE .EQ. 7 .OR. NCASE .EQ. 93)THEN
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
        ELSE IF(NCASE .GE. 23 .AND. NCASE .LE. 30)THEN
! gz[pn][1-4]: resolved +/- perturbations of GAM(4) across zero at
! magnitudes 1e-2..1e-5 of GAMP (odd=positive, even=negative). The
! |GAM|/DLN contribution at 1e-2*GAMP is ~1e-7 vs CHI ~ 1e-3 -- a
! measurable operator change at 8-digit resolution, converging to
! the zero-coupling operator as the magnitude falls.
          IF(I .EQ. 4)THEN
            T1=1.0E-2_LDP*10.0_LDP**(-DBLE((NCASE-23)/2))
            IF(MOD(NCASE,2) .EQ. 0)T1=-T1
            GAM(I,1)=T1*GAMP
          ELSE IF(I .EQ. 3 .OR. I .EQ. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 33)THEN
! abz: the ab profile with GAM(4)=0 exactly -- runs the production
! sweep plus refresh under DRV_AB so KN_DBG_NDON can prove the
! zero cell takes no expanding-donor read on either path.
          IF(I .EQ. 4)THEN
            GAM(I,1)=0.0_LDP
          ELSE IF(I .EQ. 3 .OR. I .EQ. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE IF(NCASE .EQ. 31)THEN
! gzb: exact zero at cell 4 between opposite-sign neighbours
! (expanding at 3, compressing at 5) -- the crossing sits inside the
! compressing bracket, so both bracket ends must keep their own
! upwind rules while cell 4 carries no coupling.
          IF(I .EQ. 4)THEN
            GAM(I,1)=0.0_LDP
          ELSE IF(I .EQ. 3)THEN
            GAM(I,1)=GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE IF(I .EQ. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        ELSE
          IF(I .GE. 3 .AND. I .LE. 5)THEN
            GAM(I,1)=-GAMP*(1.0_LDP+0.3_LDP*I)
          ELSE
            GAM(I,1)=GAMP*(0.5_LDP+0.1_LDP*I)
          END IF
        END IF
        IF(CNAME .EQ. 'gzzo' .OR. CNAME .EQ. 'gzzf') &
          GAM(I,1)=0.0_LDP
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
! Multi-ray cases give every ray the same atmosphere so a genuine
! solve lands on the same field the sibling reference reports.
!
      DO LS=2,KNP
        GAM(:,LS)=GAM(:,1)
        GAM_MAG(:,LS)=GAM_MAG(:,1)
        dGAMdR(:,LS)=dGAMdR(:,1)
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
! invdln (86): a degenerate (zero) frequency span is an invalid
! interval -- the BUILD_FREQ DLN guard must flag INVALID_STATE,
! not divergence.
      IF(NCASE .EQ. 86)THEN
        KN_DLN_P(2)=0.0_LDP
        KN_DLN_N(2)=0.0_LDP
      END IF
      KN_IC=0.0_LDP
      KN_DBB=0.0_LDP
!
      KN_INNER_METH='ZERO_FLUX'
      KN_METHOD_SAV='ZERO'
      KN_INSERT_SAV=.FALSE.
!
! pmon/pnom: non-ZERO captured method, MON-class vs plain. The
! refresh BUILD_FREQ clamps DTAU unconditionally now (INTEGRAL-
! path parity), so the two must run the identical operator.
!
      IF(NCASE .EQ. 51)KN_METHOD_SAV='XYZMON'
      IF(NCASE .EQ. 52)KN_METHOD_SAV='XYZ   '
!
      ALLOCATE(KN_RES_PASS(2*(60+120)),KN_HOT_CNT(NCFX),STAT=IOS)
      IF(IOS .NE. 0)GOTO 999
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
!
      FG_FREQ_MATCH=(NCASE .EQ. 3)
      FG_ORD_SWEEP=(NCASE .EQ. 2 .OR. NCASE .EQ. 35 .OR. &
        NCASE .EQ. 89 .OR. CNAME .EQ. 'gposo' .OR. &
        CNAME .EQ. 'gnego' .OR. CNAME .EQ. 'gzzo')
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
      KN_RAY_TOL=2.0E-3_LDP
      IF(NCASE .EQ. 4 .OR. NCASE .EQ. 88 .OR. NCASE .EQ. 89) &
        KN_RAY_TOL=1.0E-14_LDP
      KN_RAY_MAXPASS=60
      IF(NCASE .EQ. 4 .OR. NCASE .EQ. 88 .OR. NCASE .EQ. 89) &
        KN_RAY_MAXPASS=4
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
!
! The KR_DBG dump blocks hold the production sweep's last pass;
! the certification verify pass is an extra scratch sweep that
! would overwrite them, so the ab* operator-parity cases keep
! certification off. Production leaves it on (default).
!
      CALL SET_FG_CERT_CONV(NCASE .NE. 15 .AND. NCASE .NE. 22 &
         .AND. NCASE .NE. 32 .AND. NCASE .NE. 33)
      FG_COMMIT_JUMP=50.0_LDP
      FG_DIV_HEAL_CAP=1.0E+6_LDP
      FG_BULK_CONTR=.FALSE.
      FG_ADAPT_OMEGA_V2=.FALSE.
!
! slw/slwo/rfa (34-36): a non-converging ray that grinds the full
! pass budget so the per-solve contraction-history export (audit A1)
! accumulates post-transient samples; stall/plateau exits are
! disabled or the run stops before the estimator sees anything.
! 'rfa' also turns the refresh-level adaptive relaxation on so the
! pooled P90 estimate and the omega update run end to end.
!
      IF(NCASE .GE. 34)THEN
        KN_RAY_TOL=1.0E-14_LDP
        KN_RAY_MAXPASS=14
        KN_RAY_STALL=1.1_LDP
        FG_PLATEAU_DET=.FALSE.
      END IF
      IF(NCASE .EQ. 36)THEN
        KN_ADAPT_RELAX=.TRUE.
        FG_ADAPT_OMEGA_V2=.TRUE.
        CONV_DIAG_LVL=2
      END IF
!
! capx (37, audit C1): a low NONMONO_INEXT_CAP makes the march cap
! cells every pass, so the all-pass NCAP and the same-census
! NCAP_LST must separate by the pass multiplicity.
!
      IF(NCASE .EQ. 37)NONMONO_INEXT_CAP=1.0E-4_LDP
!
      I_P_DEPTH_NEXT=0.0_LDP
      I_M_DEPTH_NEXT=0.0_LDP
!
! FG_RAY_SEEDED=1 reproduces the pre-certification harness state
! (the array was unallocated, so the seed-fill waiver never armed
! and the jump veto always evaluated). The certification cases
! (79-83) model fresh rays instead so the waiver and the SEED
! marks are observable.
!
      FG_RAY_SEEDED=1
      IF(NCASE .EQ. 79 .OR. NCASE .EQ. 80 .OR. &
         NCASE .EQ. 81 .OR. NCASE .EQ. 82 .OR. NCASE .EQ. 83) &
         FG_RAY_SEEDED=0
!
! Seed the lagged store with a smooth field; for 'fixed' the seed is
! the analytic fixed point of a two-bin coupled solve so the pass
! residual starts small (seed-consistency check).
!
      DO LS=1,KNP
        DO K=1,N_COMP
          DO L=1,NCFX-1
            IF(NCASE .EQ. 5)THEN
              I_M_DEPTH_NEXT(K,LS,L)=1.190E-3_LDP* &
      (0.977_LDP)**(L-1)*(0.78_LDP)**(K-1)
              I_P_DEPTH_NEXT(K,LS,L)=1.173E-3_LDP* &
      (0.977_LDP)**(L-1)*(0.78_LDP)**(K-1)
            ELSE
              I_M_DEPTH_NEXT(K,LS,L)= &
      1.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
              I_P_DEPTH_NEXT(K,LS,L)= &
      2.0E-4_LDP*(1.0_LDP+0.1_LDP*L)*(1.0_LDP+0.1_LDP*K)
            END IF
          END DO
        END DO
      END DO
!
! lagc/lagr/lags: cold start -- the store carries no generation at
! all, the worst-converged legal seed state. FG_RAY_SEEDED=0 marks
! the ray uncommitted so the cold-path waiver and SEED census fire
! the way a production cold start does.
!
      IF(NCASE .EQ. 90 .OR. NCASE .EQ. 92 .OR. NCASE .EQ. 93)THEN
        I_M_DEPTH_NEXT=0.0_LDP
        I_P_DEPTH_NEXT=0.0_LDP
        FG_RAY_SEEDED=0
      END IF
!
! sdp: the retained store arrives poisoned -- a NaN slot, an Inf
! slot, negatives, and overflow-magnitude values on every
! compressing depth. A refresh must either commit only finite
! values or refuse the commit; poison must never pass through.
!
      IF(NCASE .EQ. 95)THEN
        FG_RAY_SEEDED=0
        DO K=1,N_COMP
          I_M_DEPTH_NEXT(K,1,2)=IEEE_VALUE(0.0_LDP, &
      IEEE_POSITIVE_INF)
          I_P_DEPTH_NEXT(K,1,2)=IEEE_VALUE(0.0_LDP, &
      IEEE_POSITIVE_INF)
          I_M_DEPTH_NEXT(K,1,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_P_DEPTH_NEXT(K,1,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_M_DEPTH_NEXT(K,1,4)=-5.0_LDP
          I_P_DEPTH_NEXT(K,1,4)=-7.0_LDP
          I_M_DEPTH_NEXT(K,1,5)=1.0E200_LDP
          I_P_DEPTH_NEXT(K,1,5)=1.0E200_LDP
        END DO
      END IF
!
      IF(NCASE .EQ. 38)THEN
!
! sibref: probe KN_JMP_SIBREF_V1 on a three-ray store with a
! poisoned middle ray -- +Inf (freq 2), NaN (freq 3), finite-huge
! (freq 4) and an all-floored column (freq 5). The reference must
! stay finite, equal the median of the sane ray peaks, and return
! -1 where no sane sibling exists.
!
        DO K=1,N_COMP
          I_M_DEPTH_NEXT(K,2,2)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
          I_P_DEPTH_NEXT(K,2,2)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
          I_M_DEPTH_NEXT(K,2,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_P_DEPTH_NEXT(K,2,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_M_DEPTH_NEXT(K,2,4)=1.0E200_LDP
          I_P_DEPTH_NEXT(K,2,4)=1.0E200_LDP
        END DO
        DO LS=1,KNP
          DO K=1,N_COMP
            I_M_DEPTH_NEXT(K,LS,5)=0.0_LDP
            I_P_DEPTH_NEXT(K,LS,5)=0.0_LDP
          END DO
        END DO
        ALLOCATE(JMP_SIB_M(KN_NCF),JMP_SIB_P(KN_NCF),STAT=IOS)
        CALL KN_JMP_SIBREF_V1(I_M_DEPTH_NEXT,JMP_SIB_M,IOS)
        CALL KN_JMP_SIBREF_V1(I_P_DEPTH_NEXT,JMP_SIB_P,IOS)
        WRITE(*,'(A,I3)')'SSTAT=',IOS
        DO L=1,KN_NCF
          WRITE(*,'(A,I2,A,2ES16.8)')'SREF ',L,' ', &
      JMP_SIB_M(L),JMP_SIB_P(L)
        END DO
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 50)THEN
!
! esav: JFNK eval-bracket fidelity. Sentinel-fill a representative
! spread of the pushed fields (scalars, census tallies, flags,
! counters), push, overwrite with eval-side junk, pop -- every
! printed field must come back to the production sentinel.
!
        FG_GS_MAXCHG=1.25_LDP
        FG_SW_MAXCHG=0.75_LDP
        FG_RF_DCHG=3.5_LDP
        FG_RF_PEND=2.5_LDP
        FG_GS_NCONV=42_8
        FG_GS_NBAD=7_8
        FG_GS_NSTAL=3_8
        FG_GS_NDIV=2
        FG_DIV_HEAL=9
        FG_COMMIT_REJ=4
        KN_HOT_SPENT=11_8
        KN_STORE_READY=.TRUE.
        FREQ_IDX_ITER=4
        CALL FG_EVAL_STATE_PUSH_V1()
        FG_GS_MAXCHG=999.0_LDP
        FG_SW_MAXCHG=888.0_LDP
        FG_RF_DCHG=0.0_LDP
        FG_RF_PEND=0.0_LDP
        FG_GS_NCONV=0_8
        FG_GS_NBAD=0_8
        FG_GS_NSTAL=0_8
        FG_GS_NDIV=0
        FG_DIV_HEAL=0
        FG_COMMIT_REJ=0
        KN_HOT_SPENT=0_8
        KN_STORE_READY=.FALSE.
        FREQ_IDX_ITER=0
        CALL FG_EVAL_STATE_POP_V1()
        WRITE(*,'(A,ES16.8)')'EV_MAXCHG=',FG_GS_MAXCHG
        WRITE(*,'(A,ES16.8)')'EV_SWMX=',FG_SW_MAXCHG
        WRITE(*,'(A,ES16.8)')'EV_RFD=',FG_RF_DCHG
        WRITE(*,'(A,ES16.8)')'EV_RFP=',FG_RF_PEND
        WRITE(*,'(A,ES16.8)')'EV_NCONV=',DBLE(FG_GS_NCONV)
        WRITE(*,'(A,ES16.8)')'EV_NBAD=',DBLE(FG_GS_NBAD)
        WRITE(*,'(A,ES16.8)')'EV_NSTAL=',DBLE(FG_GS_NSTAL)
        WRITE(*,'(A,ES16.8)')'EV_NDIV=',DBLE(FG_GS_NDIV)
        WRITE(*,'(A,ES16.8)')'EV_HEAL=',DBLE(FG_DIV_HEAL)
        WRITE(*,'(A,ES16.8)')'EV_JREJ=',DBLE(FG_COMMIT_REJ)
        WRITE(*,'(A,ES16.8)')'EV_HOT=',DBLE(KN_HOT_SPENT)
        WRITE(*,'(A,L2)')'EV_SRDY=',KN_STORE_READY
        WRITE(*,'(A,I8)')'EV_FRQ=',FREQ_IDX_ITER
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 56)THEN
!
! dtau: exercise the shared non-ZERO DTAU operator
! (KN_DTAU_INT_V1) directly -- the same routine the INTEGRAL
! sweep block and the refresh BUILD_FREQ both call. Profile A is
! monotone increasing CHI with a steep kink (the clamp must bind
! where the cubic overshoots), profile B is non-monotone.
!
        DO I=1,8
          R_D(I)=1.0_LDP+0.5_LDP*(I-1)
          Z_D(I)=DBLE(9-I)
        END DO
        CHI_D=(/1.0_LDP,2.0_LDP,4.0_LDP,16.0_LDP,32.0_LDP, &
                64.0_LDP,128.0_LDP,256.0_LDP/)
        DCH_D=0.0_LDP
        DCH_D(1:7)=CHI_D(2:8)-CHI_D(1:7)
        DCH_D(8)=DCH_D(7)
        CALL KN_DTAU_INT_V1(8,CHI_D,DCH_D,Z_D,R_D,DT_D)
        DO I=1,7
          WRITE(*,'(A,I2,A,3ES16.8)')'DTA ',I,' ',DT_D(I), &
      CHI_D(I)*(Z_D(I)-Z_D(I+1)),CHI_D(I+1)*(Z_D(I)-Z_D(I+1))
        END DO
        CHI_D=(/1.0_LDP,8.0_LDP,2.0_LDP,64.0_LDP,4.0_LDP, &
                32.0_LDP,16.0_LDP,128.0_LDP/)
        DCH_D=0.0_LDP
        DCH_D(1:7)=CHI_D(2:8)-CHI_D(1:7)
        DCH_D(8)=DCH_D(7)
        CALL KN_DTAU_INT_V1(8,CHI_D,DCH_D,Z_D,R_D,DT_D)
        DO I=1,7
          WRITE(*,'(A,I2,A,3ES16.8)')'DTB ',I,' ',DT_D(I), &
      CHI_D(I)*(Z_D(I)-Z_D(I+1)),CHI_D(I+1)*(Z_D(I)-Z_D(I+1))
        END DO
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 53 .OR. NCASE .EQ. 55)THEN
!
! sfail/scold: probe the sibling-reference status contract. sfail
! arms the injectable alloc-failure hook -- ISTAT must come back 1
! (never conflated with a cold start). scold floors every sibling
! store -- a genuine cold start: ISTAT=0, every REF entry -1.
!
        IF(NCASE .EQ. 53)FG_JMP_FAIL_INJ=.TRUE.
        IF(NCASE .EQ. 55)THEN
          I_M_DEPTH_NEXT=0.0_LDP
          I_P_DEPTH_NEXT=0.0_LDP
        END IF
        ALLOCATE(JMP_SIB_M(KN_NCF),JMP_SIB_P(KN_NCF),STAT=IOS)
        CALL KN_JMP_SIBREF_V1(I_M_DEPTH_NEXT,JMP_SIB_M,IOS)
        CALL KN_JMP_SIBREF_V1(I_P_DEPTH_NEXT,JMP_SIB_P,I2)
        WRITE(*,'(A,2I3)')'SSTAT=',IOS,I2
        DO L=1,KN_NCF
          WRITE(*,'(A,I2,A,2ES16.8)')'SREF ',L,' ', &
      JMP_SIB_M(L),JMP_SIB_P(L)
        END DO
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 40 .OR. NCASE .EQ. 46)THEN
!
! jmpx/jmp2: ray 1 has a floored frequency-1 column, ray 2 carries
! +Inf there, ray 3 holds a small sane peak. The refresh must
! exclude the Inf sibling, take ray 3 as the reference, and
! jump-reject ray 1's ~1E-4 fill (>> 50*1E-9) while rays 2 and 3
! commit normally. jmp2 runs the refresh twice for the repeated-
! rejection maturity check.
!
        DO K=1,N_COMP
          I_M_DEPTH_NEXT(K,1,1)=1.0E-30_LDP
          I_P_DEPTH_NEXT(K,1,1)=1.0E-30_LDP
          I_M_DEPTH_NEXT(K,3,1)=1.0E-9_LDP
          I_P_DEPTH_NEXT(K,3,1)=1.0E-9_LDP
        END DO
        I_M_DEPTH_NEXT(1,2,1)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
        I_P_DEPTH_NEXT(1,2,1)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
        I_M_DEPTH_NEXT(2,2,1)=1.0E-9_LDP
        I_M_DEPTH_NEXT(3,2,1)=1.0E-9_LDP
        I_P_DEPTH_NEXT(2,2,1)=1.0E-9_LDP
        I_P_DEPTH_NEXT(3,2,1)=1.0E-9_LDP
      END IF
!
      IF(NCASE .EQ. 49)THEN
!
! heal: a diverged ray carries a mixed poisoned store. I_M rows hold
! {1e-9, 1e-9, 0, -5, +Inf}: the zero slot is already healed and must
! not count; the negative and Inf slots do (2 per row). I_P rows
! hold {1e-9, 0, NaN, -7, 1e200}: NaN, negative and above-cap slots
! count (3 per row), the zero does not. GAM Inf forces the
! divergence so the heal path runs. Expect FG_DIV_HEAL = 3*(2+3)=15.
!
        GAM(4,1)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
        DO K=1,N_COMP
          I_M_DEPTH_NEXT(K,1,1)=1.0E-9_LDP
          I_M_DEPTH_NEXT(K,1,2)=1.0E-9_LDP
          I_M_DEPTH_NEXT(K,1,3)=0.0_LDP
          I_M_DEPTH_NEXT(K,1,4)=-5.0_LDP
          I_M_DEPTH_NEXT(K,1,5)=IEEE_VALUE(0.0_LDP,IEEE_POSITIVE_INF)
          I_P_DEPTH_NEXT(K,1,1)=1.0E-9_LDP
          I_P_DEPTH_NEXT(K,1,2)=0.0_LDP
          I_P_DEPTH_NEXT(K,1,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_P_DEPTH_NEXT(K,1,4)=-7.0_LDP
          I_P_DEPTH_NEXT(K,1,5)=1.0E200_LDP
        END DO
      END IF
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
! P1-T1g repeatability probe: snapshot the refresh input state, run
! the production refresh twice with a restore in between, and compare
! every committed store element. 'repl'/'replo' force the hot-pass
! budget and a non-converging tolerance so the shared KN_HOT_SPENT
! counter and the commit vetoes are exercised under OMP; 'replo' adds
! the ordered sweep. 'replc' keeps the converging baseline. Cases
! 59-61 are >=39 so they must exit here before the rfsh block.
!
      IF(NCASE .GE. 59 .AND. NCASE .LE. 62)THEN
        IF(NCASE .EQ. 60)FG_ORD_SWEEP=.TRUE.
        IF(NCASE .LE. 60)THEN
          KN_RAY_TOL=1.0E-14_LDP
          KN_RAY_MAXPASS=14
          KN_RAY_STALL=1.1_LDP
          KN_RAY_HOTPASSES=120
          KN_HOT_TOTAL_PASSES=480
          FG_PLATEAU_DET=.FALSE.
        END IF
! 'replh' (62) stresses the shared OMP hot-pass budget: a near-zero
! per-frequency tolerance keeps every ray's hot window open and a tiny
! refresh-wide KN_HOT_TOTAL_PASSES forces the rays to compete for a
! handful of hot passes, so SCHEDULE(DYNAMIC) ordering decides which
! rays get them. If the two DRV_REPL passes disagree, the shared
! budget is a genuine cross-ray, schedule-dependent edge.
!
        IF(NCASE .EQ. 62)THEN
          KN_RAY_TOL=1.0E-30_LDP
          KN_RAY_MAXPASS=4
          KN_RAY_STALL=0.0_LDP
          KN_RAY_HOTPASSES=60
          KN_HOT_TOTAL_PASSES=2
          KN_ADAPT_RELAX=.FALSE.
          FG_PLATEAU_DET=.FALSE.
        END IF
        CALL DRV_REPL()
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      IF(NCASE .EQ. 10 .OR. NCASE .EQ. 36 .OR. NCASE .EQ. 63 .OR. &
     &    NCASE .GE. 39)THEN
!
! End-to-end refresh: arms the capture-store guard and runs the
! production post-sweep path (OMP ray loop, census, commit).
!
        IF(NCASE .EQ. 54)FG_JMP_FAIL_INJ=.TRUE.
! 'topo' (63) arms the P1-T1h compact topology checkpoint: the
! refresh writes KNOP_TOPO at the serial point before the ray loop.
        IF(NCASE .EQ. 63)KNOP_DUMP_TOPO=.TRUE.
        DO_KNOP_REFRESH=.TRUE.
        FREEZE_I_NEXT=.FALSE.
        HAS_COMPRESSING=.TRUE.
        KN_STORE_READY=.TRUE.
        FREQ_IDX_ITER=KN_NCF-1
        FG_GS_NP=KN_NP
        FG_GS_NCF=KN_NCF
        SOLUTION_METHOD='INTEGRAL'
        NRF=1
        IF(NCASE .EQ. 46 .OR. NCASE .EQ. 71)NRF=2
! lag*/sdab/qres: multi-refresh lagged-map probes. Three refreshes
! track the outer lagged-field iteration on a fixed captured
! problem; sdab re-seeds garbage between refresh 1 and 2.
        IF(NCASE .GE. 90 .AND. NCASE .LE. 93)NRF=3
        IF(NCASE .EQ. 94)NRF=2
! qres: arm the independent fixed-point residual diagnostic so the
! committed field's true map residual is printed per ray.
        IF(NCASE .EQ. 96)KN_Q_RESID=.TRUE.
!
! prb/prbv/prbs (64-66): MODE 1 probe isolation end-to-end. Ray 3
! is armed as the diagnostic probe and its retained store is
! scaled so its solve seeds, residuals and telemetry differ from
! the two identical production rays. The sibling medians behind
! JMP_SIB are unaffected (median of {s,s,k*s} is s), so rays 1-2
! see the same inputs as 'rf3' and every production export must
! match the rf3 baseline exactly. 'prbv' scales the probe store
! down by 1E-4 so its converged candidate lands ~1E4 times the
! retained slice -- a real jump veto (relative-to-store) whose
! rejection must not reach the production counters; 'prbs'
! additionally arms KN_SEED_ARM=1 (sibling-median probe seed) to
! show the production controller values are seed-arm independent.
!
        IF(NCASE .GE. 64 .AND. NCASE .LE. 66)THEN
          CALL SET_FG_NOCP_PROBE(.TRUE.,3)
          IF(NCASE .EQ. 65)THEN
            I_M_DEPTH_NEXT(:,3,:)=I_M_DEPTH_NEXT(:,3,:)*1.0E-4_LDP
            I_P_DEPTH_NEXT(:,3,:)=I_P_DEPTH_NEXT(:,3,:)*1.0E-4_LDP
          ELSE
            I_M_DEPTH_NEXT(:,3,:)=I_M_DEPTH_NEXT(:,3,:)*2.5_LDP
            I_P_DEPTH_NEXT(:,3,:)=I_P_DEPTH_NEXT(:,3,:)*2.5_LDP
          END IF
          IF(NCASE .EQ. 66)CALL SET_FG_SEED_ARM(1)
        END IF
!
! prbz/prbl/prbr (69-71): MODE 1 empty-production omega guard.
! The single ray is armed as the diagnostic probe, leaving zero
! production rays: KN_RF_RESID takes the 1.0 fail-closed sentinel
! and the controller must hold INEXT_RELAX and KN_RF_RESID_PRV at
! their seeded values. 'prbz' runs the V2 rho-pool path; 'prbl'
! the legacy cross-refresh ratio where an unguarded sentinel
! would move omega; 'prbr' un-arms the probe on pass 2 so a real
! production refresh must resume without the sentinel as its
! residual-ratio denominator.
!
        IF(NCASE .EQ. 69 .OR. NCASE .EQ. 70)THEN
          CALL SET_FG_NOCP_PROBE(.TRUE.,1)
          KN_ADAPT_RELAX=.TRUE.
          KN_RF_RESID_PRV=0.4_LDP
          INEXT_RELAX=0.62_LDP
          FG_ADAPT_OMEGA_V2=(NCASE .EQ. 69)
        END IF
        IF(NCASE .EQ. 71)THEN
          KN_ADAPT_RELAX=.TRUE.
          KN_RF_RESID_PRV=0.7_LDP
          INEXT_RELAX=0.62_LDP
          FG_ADAPT_OMEGA_V2=.FALSE.
        END IF
!
! cerb/cer0/cern/cer2/cerp (79-83): Phase 2.15 commit gate
! end-to-end. A negative KN_RAY_TOL disables the residual
! exit so the rays march to the NPO>=6 bulk test, where the
! huge positive KN_RAY_BULK_TOL floor accepts them -- and the
! residual still reads above the (negative) tolerance, so the
! classifier lands on BULK_ONLY: capped provisional writes
! (VC=3) under FG_COMMIT_NC>0, retained store (VC=2) under
! NC=0, never a converged commit or seeded mark. 'cern'
! keeps the strict BTOL=0 path where the same rays certify
! FULL_CONVERGENCE; 'cer2' runs two refreshes to pin repeated
! provisional behaviour; 'cerp' arms ray 3 as the probe so
! its suppressed verdict stays out of the cert census.
!
        IF(NCASE .EQ. 79 .OR. NCASE .EQ. 80 .OR. &
           NCASE .EQ. 82 .OR. NCASE .EQ. 83)THEN
          KN_RAY_BULK_TOL=1.0E30_LDP
          KN_RAY_TOL=-1.0E-30_LDP
        END IF
! rfmx (84): cap every ray at one pass so all exits are
! iteration-limited -- the all-ray flag must read F.
! rfab (85): two refreshes; the first converges (flag T), the
! second is pass-capped (must reset to F, not stay stale T).
        IF(NCASE .EQ. 84)KN_RAY_MAXPASS=1
        IF(NCASE .EQ. 85)NRF=2
        IF(NCASE .EQ. 80)FG_COMMIT_NC=0.0_LDP
        IF(NCASE .EQ. 82)NRF=2
        IF(NCASE .EQ. 83)THEN
          CALL SET_FG_NOCP_PROBE(.TRUE.,3)
          I_M_DEPTH_NEXT(:,3,:)=I_M_DEPTH_NEXT(:,3,:)*2.5_LDP
          I_P_DEPTH_NEXT(:,3,:)=I_P_DEPTH_NEXT(:,3,:)*2.5_LDP
        END IF
!
! badid/bido (88/89): bad-cell identity tracking end-to-end. One
! pass-capped refresh leaves a nonzero last-pass census; the
! production merge inside FG_KNOP_REFRESH_V1 must reproduce
! FG_GS_NBAD exactly from the packed per-(ray,channel) lists.
! bido runs the same config through SWEEP_FLD_ORD.
!
        IF(NCASE .EQ. 88 .OR. NCASE .EQ. 89)THEN
          KN_BADID_ON=.TRUE.
          KN_RAY_MAXPASS=1
          IF(ALLOCATED(KN_BADRAY))DEALLOCATE(KN_BADRAY)
          ALLOCATE(KN_BADRAY(KN_NP,2),STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(*,'(A)')'RESULT_FAIL badid_alloc'
            GOTO 999
          END IF
          KN_BADRAY%N=0_8
        END IF
        DO IRF=1,NRF
        IF(NCASE .EQ. 71)CALL SET_FG_NOCP_PROBE(IRF .EQ. 1,1)
        IF(NCASE .EQ. 85 .AND. IRF .EQ. 2)KN_RAY_MAXPASS=1
!
! sdab (seed A/B): refresh 2 starts from a garbage retained field
! -- finite but ~1E3x oversized with an alternating sign sawtooth
! plus a NaN and an Inf slot, the shape of a badly-converged or
! corrupted generation. If the per-ray solve converges the
! committed field is seed-independent; the python side compares
! the refresh-1 and refresh-2 committed stores.
!
        IF(NCASE .EQ. 94 .AND. IRF .EQ. 2)THEN
          DO K=1,N_COMP
            DO L=1,NCFX-1
              T1=1.0E-1_LDP*(1.0_LDP+0.05_LDP*L)* &
      (1.0_LDP-2.0_LDP*MOD(K+L,2))
              I_M_DEPTH_NEXT(K,1,L)=T1
              I_P_DEPTH_NEXT(K,1,L)=T1
            END DO
          END DO
          I_M_DEPTH_NEXT(1,1,3)=IEEE_VALUE(0.0_LDP,IEEE_QUIET_NAN)
          I_P_DEPTH_NEXT(2,1,4)=IEEE_VALUE(0.0_LDP, &
      IEEE_POSITIVE_INF)
        END IF
        IF(NRF .GT. 1)WRITE(*,'(A,I2)')'RF_PASS ',IRF
        CALL FG_KNOP_REFRESH_V1()
        WRITE(*,'(A,ES16.8)')'RF_RESID=',KN_RF_RESID
        WRITE(*,'(A,ES16.8)')'RF_RESPRV=',KN_RF_RESID_PRV
        WRITE(*,'(A,ES16.8)')'RF_DCHG=',FG_RF_DCHG
        WRITE(*,'(A,ES16.8)')'RF_PEND=',FG_RF_PEND
        WRITE(*,'(A,I8)')'RF_NDIV=',FG_GS_NDIV
        WRITE(*,'(A,ES16.8)')'RF_MAXCHG=',FG_GS_MAXCHG
        WRITE(*,'(A,ES16.8)')'RF_STAL=',DBLE(FG_GS_NSTAL)
        WRITE(*,'(A,ES16.8)')'RF_PLAT=',DBLE(FG_GS_NPLAT)
        WRITE(*,'(A,ES16.8)')'RF_INEXT_RELAX=',INEXT_RELAX
        WRITE(*,'(A,I8)')'RF_HEAL=',FG_DIV_HEAL
        WRITE(*,'(A,L2)')'RF_FAIL=',FG_RF_FAILED
        WRITE(*,'(A,L2)')'RF_ALLCONV=',FG_RF_ALL_CONV
        WRITE(*,'(A,I8)')'RF_JREJ=',FG_COMMIT_REJ
        WRITE(*,'(A,I8)')'RF_NCREJ=',FG_NC_REJ
        WRITE(*,'(A,I8)')'RF_NCCAP=',FG_NC_CAP
        WRITE(*,'(A,I14)')'RF_NBAD=',FG_GS_NBAD
        WRITE(*,'(A,I14)')'RF_NTOT=',FG_GS_NTOT
        WRITE(*,'(A,I14)')'RF_NBB=',FG_GS_NBAD_BLK
        WRITE(*,'(A,I14)')'RF_NTB=',FG_GS_NTOT_BLK
        WRITE(*,'(A,I14)')'RF_CAPPED=',FG_RF_CAPPED
        WRITE(*,'(A,I14)')'RF_ICAP=',FG_INEXT_CAPPED
        IF(ALLOCATED(FG_RAY_VC))THEN
          DO LS=1,KN_NP
            WRITE(*,'(A,I2,A,I3,A,ES16.8)')'VC ',LS,' ', &
      FG_RAY_VC(LS),' ',FG_RAY_DCH(LS)
          END DO
        END IF
        IF(ALLOCATED(KN_CS_LS))THEN
          DO LS=1,KN_NP
            WRITE(*,'(A,I2,A,I3,A,I3)')'CLS ',LS,' ', &
      KN_CS_LS(LS),' ',KN_CS_REA(LS)
          END DO
        END IF
        IF(ALLOCATED(FG_RAY_SEEDED))THEN
          DO LS=1,KN_NP
            WRITE(*,'(A,I2,A,I2)')'SEED ',LS,' ',FG_RAY_SEEDED(LS)
          END DO
        END IF
        DO K=1,N_COMP
          DO L=1,NCFX-1
            WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'STORE ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
          END DO
        END DO
        DO LS=1,KN_NP
          DO K=1,N_COMP
            DO L=1,NCFX-1
              WRITE(*,'(A,I2,A,I2,A,I2,A,2ES16.8)') &
      'RSTORE ',LS,' ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,LS,L),I_M_DEPTH_NEXT(K,LS,L)
            END DO
          END DO
        END DO
        CALL DRV_MATURITY()
        END DO
        IF(NCASE .EQ. 88 .OR. NCASE .EQ. 89)THEN
          WRITE(*,'(A,I14)')'BID_N1=',KN_BADRAY(1,1)%N
          WRITE(*,'(A,I14)')'BID_N2=',KN_BADRAY(1,2)%N
          WRITE(*,'(A,I14)')'BID_SUM=',KN_BADRAY(1,1)%N+ &
      KN_BADRAY(1,2)%N
          CALL GET_FG_BADID(BNC,BNP,BINT,BNEW,BGON,BJAC,BNCC, &
      BNIC,BNUN)
          WRITE(*,'(A,I14)')'BID_NC=',BNC
          WRITE(*,'(A,I14)')'BID_INT=',BINT
          WRITE(*,'(A,I14)')'BID_NP=',BNP
          WRITE(*,'(A,I14)')'BID_UNTR=',BNUN
          WRITE(*,'(A,ES16.8)')'BID_JAC=',BJAC
        END IF
        IF(NCASE .EQ. 63)THEN
          INQUIRE(FILE='KNOP_TOPO',EXIST=LEX)
          WRITE(*,'(A,L2)')'TOPO_FILE=',LEX
        END IF
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
      RESID=0.0_LDP
!
! GAM-crossing cases also dump the refresh-side operator: the
! coefficient-level assertions in test_knop_ray (audit follow-up B)
! read CHI_HAT/Q/SP/SR cell by cell, which the committed store alone
! cannot express.
!
      IF(NCASE .GE. 16)THEN
        KN_DBG_LS=1
        ALLOCATE(KR_DBG_CHIH(ND,NCFX),KR_DBG_Q(ND,NCFX), &
        KR_DBG_DCHR(ND,NCFX),KR_DBG_DTAU(ND,NCFX), &
        KR_DBG_WX(ND,NCFX),KR_DBG_DLNM(ND,NCFX), &
        KR_DBG_SPM(ND,NCFX),KR_DBG_SPP(ND,NCFX), &
        KR_DBG_NFM(ND,NCFX),KR_DBG_IKC1(ND,NCFX), &
        KR_DBG_IKC2(ND,NCFX),KR_DBG_JPR1(ND,NCFX), &
        KR_DBG_JPR2(ND,NCFX),KR_DBG_SR(ND,NCFX), &
        KR_DBG_A0(ND,NCFX),KR_DBG_A1(ND,NCFX), &
        KR_DBG_A2(ND,NCFX),KR_DBG_A3(ND,NCFX),KR_DBG_A4(ND,NCFX), &
        KR_DBG_FLDM(ND,NCFX),KR_DBG_FLDP(ND,NCFX), &
        KN_DBG_NDON(ND,NCFX),STAT=IOS)
        KR_DBG_CHIH=-1.0_LDP; KR_DBG_Q=-1.0_LDP
        KR_DBG_DCHR=-1.0_LDP; KR_DBG_DTAU=-1.0_LDP
        KR_DBG_WX=-1.0_LDP; KR_DBG_DLNM=-1.0_LDP
        KR_DBG_SPM=-1.0_LDP; KR_DBG_SPP=-1.0_LDP
        KR_DBG_NFM=-1; KR_DBG_IKC1=-1; KR_DBG_IKC2=-1
        KR_DBG_JPR1=-1; KR_DBG_JPR2=-1; KR_DBG_SR=-1.0_LDP
        KR_DBG_A0=-1.0_LDP; KR_DBG_A1=-1.0_LDP; KR_DBG_A2=-1.0_LDP
        KR_DBG_A3=-1.0_LDP; KR_DBG_A4=-1.0_LDP
        KR_DBG_FLDM=-1.0_LDP; KR_DBG_FLDP=-1.0_LDP
        KN_DBG_NDON=0
      END IF
!
! FG_JUMP_MAX now reads its sibling reference from the JMP_SIB
! snapshot arrays (built by FG_KNOP_REFRESH_V1 in production); a
! direct ray-solve call must fill them the same way.
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
      WRITE(*,'(A,I4)')'RESULT resid_ok=',1
      WRITE(*,'(A,I0,A,I0)')'DISPATCH ordered=',KN_CENS_ORD_RAY, &
      ' frequency=',KN_CENS_UNO_RAY
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
      IF(ALLOCATED(KN_CS_LS))THEN
        WRITE(*,'(A,I3,A,I3)')'CLS_R=',KN_CS_LS(1), &
      ' ',KN_CS_REA(1)
      END IF
      IF(NCASE .GE. 16 .AND. ALLOCATED(KR_DBG_CHIH))THEN
        DO L=1,NCFX
          DO I=1,ND
            WRITE(*,'(A,I2,A,I2,A,11ES16.8)')'DBR ',I,' ',L,' ', &
      KR_DBG_CHIH(I,L),KR_DBG_Q(I,L),KR_DBG_DCHR(I,L), &
      KR_DBG_DTAU(I,L),KR_DBG_SPM(I,L),KR_DBG_SPP(I,L), &
      KR_DBG_WX(I,L),KR_DBG_DLNM(I,L),KR_DBG_FLDM(I,L), &
      KR_DBG_FLDP(I,L),KR_DBG_SR(I,L)
            WRITE(*,'(A,I2,A,I2,A,5I6)')'DBRI ',I,' ',L,' ', &
      KR_DBG_NFM(I,L),KR_DBG_IKC1(I,L),KR_DBG_IKC2(I,L), &
      KR_DBG_JPR1(I,L),KR_DBG_JPR2(I,L)
!
! DBRC: the march recurrence coefficients A0..A4 plus DTAU, so the
! operator identities A0=exp(-DTAU) and A0+A1+A2=1 (constant-source
! equilibrium) are checkable cell by cell.
!
            WRITE(*,'(A,I2,A,I2,A,6ES16.8)')'DBRC ',I,' ',L,' ', &
      KR_DBG_A0(I,L),KR_DBG_A1(I,L),KR_DBG_A2(I,L), &
      KR_DBG_A3(I,L),KR_DBG_A4(I,L),KR_DBG_DTAU(I,L)
            WRITE(*,'(A,I2,A,I2,A,I6)')'DBRD ',I,' ',L,' ', &
      KN_DBG_NDON(I,L)
          END DO
        END DO
        DO I=1,ND
          WRITE(*,'(A,I2,A,ES16.8)')'GAMV ',I,' ',GAM(I,1)
        END DO
      END IF
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
      SUBROUTINE DRV_AB(ORDSW,NRAYS,LSDBG)
      LOGICAL ORDSW
      INTEGER NRAYS,LSDBG
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=8,NPA=2,NCFA=6)
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
      PA(2)=1.2E4_LDP
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
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
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
      KN_DBG_LS=LSDBG
      ALLOCATE(KN_DBG_CHIH(NDA,NCFA),KN_DBG_Q(NDA,NCFA), &
      KN_DBG_DCHR(NDA,NCFA),KN_DBG_DTAU(NDA,NCFA), &
      KN_DBG_SPM(NDA,NCFA),KN_DBG_SPP(NDA,NCFA), &
      KN_DBG_WX(NDA,NCFA),KN_DBG_DLNM(NDA,NCFA), &
      KN_DBG_IMPR(NDA,NCFA),KN_DBG_IPPR(NDA,NCFA), &
      KN_DBG_NM(NDA,NCFA),KN_DBG_KC1(NDA,NCFA),KN_DBG_KC2(NDA,NCFA), &
      KN_DBG_A0(NDA,NCFA),KN_DBG_A1(NDA,NCFA),KN_DBG_A2(NDA,NCFA), &
      KN_DBG_A3(NDA,NCFA),KN_DBG_A4(NDA,NCFA),KN_DBG_SR(NDA,NCFA), &
      KR_DBG_CHIH(NDA,NCFA),KR_DBG_Q(NDA,NCFA), &
      KR_DBG_DCHR(NDA,NCFA),KR_DBG_DTAU(NDA,NCFA), &
      KR_DBG_WX(NDA,NCFA),KR_DBG_DLNM(NDA,NCFA), &
      KR_DBG_SPM(NDA,NCFA),KR_DBG_SPP(NDA,NCFA), &
      KR_DBG_A0(NDA,NCFA),KR_DBG_A1(NDA,NCFA),KR_DBG_A2(NDA,NCFA), &
      KR_DBG_A3(NDA,NCFA),KR_DBG_A4(NDA,NCFA),KR_DBG_SR(NDA,NCFA), &
      KR_DBG_FLDM(NDA,NCFA),KR_DBG_FLDP(NDA,NCFA), &
      KR_DBG_NFM(NDA,NCFA),KR_DBG_IKC1(NDA,NCFA), &
      KR_DBG_IKC2(NDA,NCFA),KR_DBG_JPR1(NDA,NCFA), &
      KR_DBG_JPR2(NDA,NCFA),KN_DBG_NDON(NDA,NCFA),STAT=IOS)
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
      KN_DBG_A3=-1.0_LDP; KN_DBG_A4=-1.0_LDP; KN_DBG_SR=-1.0_LDP
      KR_DBG_CHIH=-1.0_LDP; KR_DBG_Q=-1.0_LDP; KR_DBG_DCHR=-1.0_LDP
      KR_DBG_DTAU=-1.0_LDP; KR_DBG_WX=-1.0_LDP; KR_DBG_DLNM=-1.0_LDP
      KR_DBG_SPM=-1.0_LDP; KR_DBG_SPP=-1.0_LDP
      KR_DBG_A0=-1.0_LDP; KR_DBG_A1=-1.0_LDP; KR_DBG_A2=-1.0_LDP
      KR_DBG_A3=-1.0_LDP; KR_DBG_A4=-1.0_LDP
      KR_DBG_FLDM=-1.0_LDP; KR_DBG_FLDP=-1.0_LDP
      KR_DBG_SR=-1.0_LDP
      KR_DBG_NFM=-1; KR_DBG_IKC1=-1; KR_DBG_IKC2=-1
      KR_DBG_JPR1=-1; KR_DBG_JPR2=-1
      KN_DBG_NDON=0
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
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
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
          WRITE(*,'(A,I2,A,I2,A,11ES16.8)')'DBS ',I,' ',L,' ', &
      KN_DBG_CHIH(I,L),KN_DBG_Q(I,L),KN_DBG_DCHR(I,L), &
      KN_DBG_DTAU(I,L),KN_DBG_SPM(I,L),KN_DBG_SPP(I,L), &
      KN_DBG_WX(I,L),KN_DBG_DLNM(I,L),KN_DBG_IMPR(I,L), &
      KN_DBG_IPPR(I,L),KN_DBG_SR(I,L)
          WRITE(*,'(A,I2,A,I2,A,3I6)')'DBSI ',I,' ',L,' ', &
      KN_DBG_NM(I,L),KN_DBG_KC1(I,L),KN_DBG_KC2(I,L)
          WRITE(*,'(A,I2,A,I2,A,5ES16.8)')'DBSA ',I,' ',L,' ', &
      KN_DBG_A0(I,L),KN_DBG_A1(I,L),KN_DBG_A2(I,L), &
      KN_DBG_A3(I,L),KN_DBG_A4(I,L)
          WRITE(*,'(A,I2,A,I2,A,I6)')'DBSD ',I,' ',L,' ', &
      KN_DBG_NDON(I,L)
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
          WRITE(*,'(A,I2,A,I2,A,11ES16.8)')'DBR ',I,' ',L,' ', &
      KR_DBG_CHIH(I,L),KR_DBG_Q(I,L),KR_DBG_DCHR(I,L), &
      KR_DBG_DTAU(I,L),KR_DBG_SPM(I,L),KR_DBG_SPP(I,L), &
      KR_DBG_WX(I,L),KR_DBG_DLNM(I,L),KR_DBG_FLDM(I,L), &
      KR_DBG_FLDP(I,L),KR_DBG_SR(I,L)
          WRITE(*,'(A,I2,A,I2,A,5I6)')'DBRI ',I,' ',L,' ', &
      KR_DBG_NFM(I,L),KR_DBG_IKC1(I,L),KR_DBG_IKC2(I,L), &
      KR_DBG_JPR1(I,L),KR_DBG_JPR2(I,L)
          WRITE(*,'(A,I2,A,I2,A,5ES16.8)')'DBRA ',I,' ',L,' ', &
      KR_DBG_A0(I,L),KR_DBG_A1(I,L),KR_DBG_A2(I,L), &
      KR_DBG_A3(I,L),KR_DBG_A4(I,L)
          WRITE(*,'(A,I2,A,I2,A,I6)')'DBRD ',I,' ',L,' ', &
      KN_DBG_NDON(I,L)
        END DO
      END DO
      DO K=1,NDA
        WRITE(*,'(A,I2,A,2I6)')'JMAP ',K,' ',J_PNT(K,LSDBG), &
      COMP_DEPTH_IDX(K)
      END DO
      DO I=1,NI_RAY(LSDBG)
        WRITE(*,'(A,I2,A,ES16.8)')'GAMV ',I,' ',GAM(I,LSDBG)
      END DO
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_AB
!
! ----------------------------------------------------------------------
! P1-T1 field-sensitivity probe. Drives the production in-sweep solver
! (FG_J_CMF_NONMONO_V1) on a manufactured non-monotonic model whose
! compressing cells reach |Q| -> 1 (opacity depressed in the reversal
! zone so |GAM|/DLN >> CHI). Dumps, for the selected ray and every
! frequency, the in-sweep operator coefficients (KN_DBG_*) plus the
! production I_M/I_P intensities so a python re-march can be validated
! cell-by-cell before finite-differencing the lagged I_NEXT channel.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_P1T1
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=8,NPA=2,NCFA=6)
      REAL(KIND=LDP) RA(NDA),VA(NDA),SGA(NDA),SGT(NDA),PA(NPA)
      REAL(KIND=LDP) ETA_M(NDA,NCFA),CHI_M(NDA,NCFA),ESC_M(NDA,NCFA)
      REAL(KIND=LDP) JNU_A(NDA),FEDD_A(NDA),IPLUS_A(NPA)
      REAL(KIND=LDP) JQW_A(NDA,NPA),HQW_A(NDA,NPA),KQW_A(NDA,NPA)
      REAL(KIND=LDP) NQW_A(NDA,NPA),HMID_A(NDA,NPA),NMID_A(NDA,NPA)
      REAL(KIND=LDP) VDOP_A(NDA),FRQV(NCFA),INBC_A,HBC_A
      REAL(KIND=LDP) LVR(NDA),BUMP,V0,DLOGA,DBB_A,IC_A
      REAL(KIND=LDP) VDOP_F,REXT_A,DEPF
      CHARACTER*6 METH_A
      CHARACTER(LEN=8) SOPT_A
      CHARACTER(LEN=12) IBND_A
      LOGICAL THK_A,INIT_A,NEWF_A
      INTEGER NC_A,NITER,ITR,L,I,LSDBG,IOS,NRAYS
      REAL(KIND=LDP) SMAX
      REAL(KIND=LDP), ALLOCATABLE :: IM0(:,:),IP0(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: SP0(:,:),SM0(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: KQ0(:),KIM0(:),KIP0(:)
      INTEGER, ALLOCATABLE :: KKC10(:),KKC20(:),KNM0(:)
      REAL(KIND=LDP) EPSV,OLDV,DIMV,DIPV,GAIN,AMPF,RMAX
      INTEGER KC,SL,FS,I2,IA,IMX
!
      LSDBG=1
      NRAYS=1
      RA=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
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
! Opacity: base profile times a depth-dependent depression centred on
! the reversal zone (depths 3-5) so the compressing cells reach
! |GAM|/DLN >> CHI and |Q| -> 1. DEPF is the multiplicative factor.
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
      DEPF=1.0E-1_LDP
      CHI_M(2,:)=CHI_M(2,:)*DEPF
      CHI_M(3,:)=CHI_M(3,:)*1.0E-3_LDP
      CHI_M(4,:)=CHI_M(4,:)*1.0E-5_LDP
      CHI_M(5,:)=CHI_M(5,:)*1.0E-7_LDP
!
      FRQV(1)=6.0E15_LDP
      DO L=2,NCFA
        FRQV(L)=FRQV(L-1)*0.75_LDP
      END DO
      DLOGA=LOG(1.0_LDP/0.75_LDP)
!
      PA(1)=0.0_LDP
      PA(2)=1.2E4_LDP
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
      FG_FREQ_MATCH=.TRUE.
      FG_ORD_SWEEP=.TRUE.
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
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
      KN_DBG_LS=LSDBG
      ALLOCATE(KN_DBG_CHIH(NDA,NCFA),KN_DBG_Q(NDA,NCFA), &
      KN_DBG_DCHR(NDA,NCFA),KN_DBG_DTAU(NDA,NCFA), &
      KN_DBG_SPM(NDA,NCFA),KN_DBG_SPP(NDA,NCFA), &
      KN_DBG_WX(NDA,NCFA),KN_DBG_DLNM(NDA,NCFA), &
      KN_DBG_IMPR(NDA,NCFA),KN_DBG_IPPR(NDA,NCFA), &
      KN_DBG_NM(NDA,NCFA),KN_DBG_KC1(NDA,NCFA),KN_DBG_KC2(NDA,NCFA), &
      KN_DBG_A0(NDA,NCFA),KN_DBG_A1(NDA,NCFA),KN_DBG_A2(NDA,NCFA), &
      KN_DBG_A3(NDA,NCFA),KN_DBG_A4(NDA,NCFA),KN_DBG_SR(NDA,NCFA), &
      KN_DBG_INEXTM(NDA,NCFA),KN_DBG_INEXTP(NDA,NCFA), &
      KN_DBG_NDON(NDA,NCFA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1alloc_fail=1'
        RETURN
      END IF
      KN_DBG_CHIH=-1.0_LDP; KN_DBG_Q=-1.0_LDP; KN_DBG_DCHR=-1.0_LDP
      KN_DBG_DTAU=-1.0_LDP; KN_DBG_SPM=-1.0_LDP; KN_DBG_SPP=-1.0_LDP
      KN_DBG_WX=-1.0_LDP; KN_DBG_DLNM=-1.0_LDP
      KN_DBG_IMPR=-1.0_LDP; KN_DBG_IPPR=-1.0_LDP
      KN_DBG_NM=-1; KN_DBG_KC1=-1; KN_DBG_KC2=-1
      KN_DBG_A0=-1.0_LDP; KN_DBG_A1=-1.0_LDP; KN_DBG_A2=-1.0_LDP
      KN_DBG_A3=-1.0_LDP; KN_DBG_A4=-1.0_LDP; KN_DBG_SR=-1.0_LDP
      KN_DBG_INEXTM=-1.0_LDP; KN_DBG_INEXTP=-1.0_LDP
      KN_DBG_NDON=0
!
! In-sweep great iterations. On the last iteration, dump the
! production I_M/I_P intensities for the selected ray after each
! frequency call so the python re-march can be validated against the
! same solve state that produced the KN_DBG coefficients.
!
      NITER=8
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
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
          IF(ITR .EQ. NITER)THEN
            DO I=1,NI_RAY(LSDBG)
              WRITE(*,'(A,I3,A,I3,A,2ES16.8)')'IMP ',L,' ',I,' ', &
      I_M(I,LSDBG),I_P(I,LSDBG)
            END DO
          END IF
        END DO
      END DO
!
      DO K=1,N_COMP
        DO L=1,NCFA-1
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'SWSTOR ',K,' ',L,' ', &
      I_P_DEPTH_NEXT(K,1,L),I_M_DEPTH_NEXT(K,1,L)
        END DO
      END DO
      DO I=1,NDA
        DO L=1,NCFA
          WRITE(*,'(A,I2,A,I2,A,11ES16.8)')'DBS ',I,' ',L,' ', &
      KN_DBG_CHIH(I,L),KN_DBG_Q(I,L),KN_DBG_DCHR(I,L), &
      KN_DBG_DTAU(I,L),KN_DBG_SPM(I,L),KN_DBG_SPP(I,L), &
      KN_DBG_WX(I,L),KN_DBG_DLNM(I,L),KN_DBG_IMPR(I,L), &
      KN_DBG_IPPR(I,L),KN_DBG_SR(I,L)
          WRITE(*,'(A,I2,A,I2,A,3I6)')'DBSI ',I,' ',L,' ', &
      KN_DBG_NM(I,L),KN_DBG_KC1(I,L),KN_DBG_KC2(I,L)
          WRITE(*,'(A,I2,A,I2,A,5ES16.8)')'DBSA ',I,' ',L,' ', &
      KN_DBG_A0(I,L),KN_DBG_A1(I,L),KN_DBG_A2(I,L), &
      KN_DBG_A3(I,L),KN_DBG_A4(I,L)
          WRITE(*,'(A,I2,A,I2,A,I6)')'DBSD ',I,' ',L,' ', &
      KN_DBG_NDON(I,L)
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'DBSN ',I,' ',L,' ', &
      KN_DBG_INEXTM(I,L),KN_DBG_INEXTP(I,L)
        END DO
      END DO
      DO K=1,NDA
        WRITE(*,'(A,I2,A,2I6)')'JMAP ',K,' ',J_PNT(K,LSDBG), &
      COMP_DEPTH_IDX(K)
      END DO
      DO I=1,NI_RAY(LSDBG)
        WRITE(*,'(A,I2,A,ES16.8)')'GAMV ',I,' ',GAM(I,LSDBG)
      END DO
!
! ---- P1-T1 finite-difference probe (diagnostic only) ----
! Re-run the frequency sweep with FREEZE_I_NEXT=T so the lagged I_NEXT
! store stays fixed at its post-warmup value. For each compressing
! cell, perturb the store elements it reads and re-solve the same
! frequency with NEW_FREQ=F to measure dI_formal/dI_NEXT. JFNK_IN_EVAL
! keeps the L=1 INIT from disarming the freeze. I_*_STORE (the
! previous-frequency upwind donor) is restored after each frequency so
! the next frequency's I_*_PREV stays clean.
!
      ALLOCATE(IM0(SIZE(I_M,1),SIZE(I_M,2)),IP0(SIZE(I_P,1),SIZE(I_P,2)),&
      SP0(SIZE(I_P_STORE,1),SIZE(I_P_STORE,2)),SM0(SIZE(I_M_STORE,1),&
      SIZE(I_M_STORE,2)),KQ0(NDA),KIM0(NDA),KIP0(NDA),KKC10(NDA),&
      KKC20(NDA),KNM0(NDA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1fd_alloc_fail=1'
        RETURN
      END IF
      JFNK_IN_EVAL=.TRUE.
      FREEZE_I_NEXT=.TRUE.
      DO L=1,NCFA
        IF(L .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)*0.75_LDP)
        END IF
        INIT_A=(L .EQ. 1)
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L),ESC_M(1,L), &
    VA,SGA,RA,PA, &
    JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
    INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA, &
    IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
    METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
        IM0=I_M; IP0=I_P; SP0=I_P_STORE; SM0=I_M_STORE
!
! Snapshot the KN_DBG diagnostics for this frequency before any
! perturbation: each re-solve overwrites KN_DBG(:,L) for every cell,
! so the baseline values must be captured up front.
!
        KQ0=KN_DBG_Q(:,L); KIM0=KN_DBG_INEXTM(:,L)
        KIP0=KN_DBG_INEXTP(:,L); KKC10=KN_DBG_KC1(:,L)
        KKC20=KN_DBG_KC2(:,L); KNM0=KN_DBG_NM(:,L)
!
! Probe compressing cells (Q<0 selects the I_NEXT branch). For each,
! perturb the lagged store coherently and re-solve at fixed frequency.
!
        DO I=1,NI_RAY(LSDBG)
          IF(KQ0(I) .LT. -0.3_LDP .AND. KKC10(I) .GE. 1)THEN
            FS=FREQ_IDX_ITER+KNM0(I)
            IF(FS .LT. 1 .OR. FS .GT. UBOUND(I_P_DEPTH_NEXT,3))THEN
              WRITE(*,'(A,I3,A,I3,A,I3,A,I3)')'FDSKIP ',L,' ',I,' ', &
          FS,' ',UBOUND(I_P_DEPTH_NEXT,3)
              CYCLE
            END IF
            KC=KKC10(I); SL=KKC20(I)
            IF(KC .LT. 1 .OR. KC .GT. N_COMP)THEN
              WRITE(*,'(A,I3,A,I3,A,I3)')'FDSKIP ',L,' ',I,' ',KC
              CYCLE
            END IF
            IF(SL .LT. 1 .OR. SL .GT. N_COMP)SL=0
            DO IA=1,2
              AMPF=10.0_LDP**(-(IA+3))
!
! Direction P: perturb I_P_DEPTH_NEXT, measure I_P response.
!
              EPSV=AMPF*MAX(ABS(IP0(I,LSDBG)),1.0E-30_LDP)
              OLDV=I_P_DEPTH_NEXT(KC,LSDBG,FS)
              I_P_DEPTH_NEXT(KC,LSDBG,FS)=OLDV+EPSV
              IF(SL .GE. 1 .AND. SL .NE. KC) &
          I_P_DEPTH_NEXT(SL,LSDBG,FS)=I_P_DEPTH_NEXT(SL,LSDBG,FS)+EPSV
              CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L), &
          ESC_M(1,L),VA,SGA,RA,PA,JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A, &
          HMID_A,NMID_A,INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA,IBND_A,DBB_A, &
          IC_A,VDOP_A,VDOP_F,REXT_A,METH_A,SOPT_A,THK_A,.FALSE., &
          .FALSE.,NC_A,NRAYS,NDA,SGT,NCFA)
              DIPV=I_P(I,LSDBG)-IP0(I,LSDBG)
              IMX=MAXLOC(ABS(I_P(1:NI_RAY(LSDBG),LSDBG)- &
          IP0(1:NI_RAY(LSDBG),LSDBG)),DIM=1)
              RMAX=MAXVAL(ABS(I_P(1:NI_RAY(LSDBG),LSDBG)- &
          IP0(1:NI_RAY(LSDBG),LSDBG)))
              I_P_DEPTH_NEXT(KC,LSDBG,FS)=OLDV
              IF(SL .GE. 1 .AND. SL .NE. KC) &
          I_P_DEPTH_NEXT(SL,LSDBG,FS)=I_P_DEPTH_NEXT(SL,LSDBG,FS)-EPSV
              GAIN=DIPV/EPSV
              WRITE(*,*)'FDRP',L,I,FS,EPSV,KNM0(I),KQ0(I),KIP0(I), &
          IP0(I,LSDBG),DIPV,GAIN,IMX,RMAX
!
! Direction M: perturb I_M_DEPTH_NEXT, measure I_M response.
!
              EPSV=AMPF*MAX(ABS(IM0(I,LSDBG)),1.0E-30_LDP)
              OLDV=I_M_DEPTH_NEXT(KC,LSDBG,FS)
              I_M_DEPTH_NEXT(KC,LSDBG,FS)=OLDV+EPSV
              IF(SL .GE. 1 .AND. SL .NE. KC) &
          I_M_DEPTH_NEXT(SL,LSDBG,FS)=I_M_DEPTH_NEXT(SL,LSDBG,FS)+EPSV
              CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L), &
          ESC_M(1,L),VA,SGA,RA,PA,JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A, &
          HMID_A,NMID_A,INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA,IBND_A,DBB_A, &
          IC_A,VDOP_A,VDOP_F,REXT_A,METH_A,SOPT_A,THK_A,.FALSE., &
          .FALSE.,NC_A,NRAYS,NDA,SGT,NCFA)
              DIMV=I_M(I,LSDBG)-IM0(I,LSDBG)
              IMX=MAXLOC(ABS(I_M(1:NI_RAY(LSDBG),LSDBG)- &
          IM0(1:NI_RAY(LSDBG),LSDBG)),DIM=1)
              RMAX=MAXVAL(ABS(I_M(1:NI_RAY(LSDBG),LSDBG)- &
          IM0(1:NI_RAY(LSDBG),LSDBG)))
              I_M_DEPTH_NEXT(KC,LSDBG,FS)=OLDV
              IF(SL .GE. 1 .AND. SL .NE. KC) &
          I_M_DEPTH_NEXT(SL,LSDBG,FS)=I_M_DEPTH_NEXT(SL,LSDBG,FS)-EPSV
              GAIN=DIMV/EPSV
              WRITE(*,*)'FDRM',L,I,FS,EPSV,KNM0(I),KQ0(I),KIM0(I), &
          IM0(I,LSDBG),DIMV,GAIN,IMX,RMAX
            END DO
          END IF
        END DO
        I_P_STORE=SP0; I_M_STORE=SM0
      END DO
      FREEZE_I_NEXT=.FALSE.
      JFNK_IN_EVAL=.FALSE.
      DEALLOCATE(IM0,IP0,SP0,SM0,KQ0,KIM0,KIP0,KKC10,KKC20,KNM0)
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_P1T1
!
!
! ----------------------------------------------------------------------
! Omega-pool probe (audit Issues 2 + 5, tests T25-T28). Feeds
! KN_OMEGA_POOL_V1 synthetic per-ray contraction samples and prints
! the returned omega, sample count and pooled statistics so the
! python side can check the small-sample retain path, the P90
! outlier rejection and the damped-vs-recovered regimes.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_OMEGA(NCS)
      INTEGER NCS
      INTEGER NRAY,I2,NDV_T(40)
      REAL(KIND=LDP) RHOM_T(40),RHOP_T(40)
      REAL(KIND=LDP) OM0,OM1,RHO_P90,RHO_MED,RHO_MX
      INTEGER NSAM
!
      NRAY=0
      NDV_T=0
      RHOM_T=0.0_LDP; RHOP_T=0.0_LDP
      OM0=0.75_LDP
      IF(NCS .EQ. 41)THEN
! ompos: slowly contracting positive-mode pool (rho < 1). The
! heuristic should recover omega toward 1/(1+rho_p90), lightly
! damped against the seed.
        NRAY=20
        DO I2=1,NRAY
          RHOM_T(I2)=0.60_LDP
          RHOP_T(I2)=0.62_LDP
        END DO
      ELSE IF(NCS .EQ. 42)THEN
! omneg: near-unity and expansive magnitudes -- the signature a
! sign-flipping or diverging mode leaves in an unsigned ratio.
! omega must be damped below the seed.
        NRAY=20
        DO I2=1,NRAY
          RHOM_T(I2)=1.05_LDP+0.005_LDP*I2
          RHOP_T(I2)=1.20_LDP
        END DO
      ELSE IF(NCS .EQ. 43)THEN
! omone: a single valid sample must not update omega.
        NRAY=1
        RHOM_T(1)=0.7_LDP
      ELSE IF(NCS .EQ. 44)THEN
! omfew: 8 samples (4 rays, both directions) is still below the
! P90 minimum -- retain the seed.
        NRAY=4
        DO I2=1,NRAY
          RHOM_T(I2)=0.4_LDP
          RHOP_T(I2)=0.9_LDP
        END DO
      ELSE IF(NCS .EQ. 45)THEN
! omout: 40 sane samples + one 100x outlier + 2 diverged rays.
! The P90 index must ignore the outlier and the diverged rays.
        NRAY=23
        DO I2=1,20
          RHOM_T(I2)=0.55_LDP+0.01_LDP*I2
          RHOP_T(I2)=0.55_LDP+0.01_LDP*I2
        END DO
        RHOM_T(21)=100.0_LDP
        RHOP_T(21)=95.0_LDP
        NDV_T(22)=1; NDV_T(23)=1
        RHOM_T(22)=500.0_LDP; RHOP_T(22)=500.0_LDP
        RHOM_T(23)=400.0_LDP; RHOP_T(23)=400.0_LDP
      END IF
!
      CALL KN_OMEGA_POOL_V1(RHOM_T,RHOP_T,NDV_T,OM0,OM1, &
      NSAM,RHO_P90,RHO_MED,RHO_MX)
      WRITE(*,'(A,ES16.8)')'OM_PREV=',OM0
      WRITE(*,'(A,ES16.8)')'OM_NEW=',OM1
      WRITE(*,'(A,I8)')'OM_NSAM=',NSAM
      WRITE(*,'(A,ES16.8)')'OM_P90=',RHO_P90
      WRITE(*,'(A,ES16.8)')'OM_MED=',RHO_MED
      WRITE(*,'(A,ES16.8)')'OM_MAX=',RHO_MX
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_OMEGA
!
! ----------------------------------------------------------------------
! DRV_M1UNIT ('punit', 67): unit checks for the MODE 1 isolation
! helpers on synthetic arrays -- KN_PROD_RESID_V1 fallback and
! exclusion behaviour, KN_OMEGA_POOL_V1 EXCL pooling, and the
! KN_NOCP_PROBE_RAY predicate. No ray grid needed.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_M1UNIT()
      INTEGER NDV_U(40)
      REAL(KIND=LDP) RES_U(40),RHOM_U(40),RHOP_U(40)
      REAL(KIND=LDP) OM0,OM1,OMB,RHO_P90,RHO_MED,RHO_MX
      REAL(KIND=LDP) OM1B,P90B,MXB,MEDB
      INTEGER NSAM,NSB,I2,ISTU
      LOGICAL L1,L2,L3
!
! --- residual pool ---
!
      NDV_U=0
      RES_U(1)=0.1_LDP; RES_U(2)=0.3_LDP; RES_U(3)=5.0_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_PX=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),3)
      WRITE(*,'(A,ES16.8)')'M1U_RES_P0=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),0)
      RES_U(3)=0.2_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_LOW=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),3)
      NDV_U(1)=1; NDV_U(2)=1
      RES_U(1)=1.0E30_LDP; RES_U(2)=1.0E30_LDP; RES_U(3)=5.0_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_ALLDV=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),3)
      NDV_U(3)=1
      RES_U(1)=7.0_LDP; RES_U(2)=8.0_LDP; RES_U(3)=9.0_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_ALLDVP=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),3)
      NDV_U(1:3)=0
      RES_U(1)=9.0_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_SOLO=', &
        KN_PROD_RESID_V1(RES_U(1:1),NDV_U(1:1),1)
      WRITE(*,'(A,I4)')'M1U_IST_SOLO=', &
        KN_PROD_POOL_STAT(NDV_U(1:1),1)
      RES_U(1)=1.0E-9_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_SOLO2=', &
        KN_PROD_RESID_V1(RES_U(1:1),NDV_U(1:1),1)
      WRITE(*,'(A,I4)')'M1U_IST_SOLO2=', &
        KN_PROD_POOL_STAT(NDV_U(1:1),1)
      RES_U(1)=1.0E30_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_SOLO3=', &
        KN_PROD_RESID_V1(RES_U(1:1),NDV_U(1:1),1)
      NDV_U(1)=1; NDV_U(2)=1
      RES_U(1)=1.0E30_LDP; RES_U(2)=1.0E30_LDP; RES_U(3)=5.0_LDP
      WRITE(*,'(A,ES16.8)')'M1U_RES_DVIST=', &
        KN_PROD_RESID_V1(RES_U(1:3),NDV_U(1:3),3)
      WRITE(*,'(A,I4)')'M1U_IST_ALLDV=', &
        KN_PROD_POOL_STAT(NDV_U(1:3),3)
!
! --- omega pool exclusion ---
!
      NDV_U=0
      DO I2=1,20
        RHOM_U(I2)=0.60_LDP
        RHOP_U(I2)=0.62_LDP
      END DO
      RHOM_U(21)=1.4_LDP; RHOP_U(21)=1.5_LDP
      OM0=0.75_LDP
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:20),RHOP_U(1:20), &
        NDV_U(1:20),OM0,OMB,NSB,P90B,MEDB,MXB)
      WRITE(*,'(A,ES16.8)')'M1U_OM_BASE=',OMB
      WRITE(*,'(A,I8)')'M1U_NSAM_BASE=',NSB
      WRITE(*,'(A,ES16.8)')'M1U_P90_BASE=',P90B
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:21),RHOP_U(1:21), &
        NDV_U(1:21),OM0,OM1,NSAM,RHO_P90,RHO_MED,RHO_MX,21)
      WRITE(*,'(A,ES16.8)')'M1U_OM_EXHI=',OM1
      WRITE(*,'(A,I8)')'M1U_NSAM_EXHI=',NSAM
      WRITE(*,'(A,ES16.8)')'M1U_P90_EXHI=',RHO_P90
      WRITE(*,'(A,ES16.8)')'M1U_RHOP_DIAG=',RHOP_U(21)
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:21),RHOP_U(1:21), &
        NDV_U(1:21),OM0,OM1B,NSB,P90B,MEDB,MXB)
      WRITE(*,'(A,ES16.8)')'M1U_OM_INCL=',OM1B
      WRITE(*,'(A,I8)')'M1U_NSAM_INCL=',NSB
      WRITE(*,'(A,ES16.8)')'M1U_P90_INCL=',P90B
      RHOM_U(21)=0.05_LDP; RHOP_U(21)=0.06_LDP
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:21),RHOP_U(1:21), &
        NDV_U(1:21),OM0,OM1,NSAM,RHO_P90,RHO_MED,RHO_MX,21)
      WRITE(*,'(A,ES16.8)')'M1U_OM_EXLO=',OM1
      WRITE(*,'(A,I8)')'M1U_NSAM_EXLO=',NSAM
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:5),RHOP_U(1:5),NDV_U(1:5), &
        OM0,OM1B,NSB,P90B,MEDB,MXB,5)
      WRITE(*,'(A,ES16.8)')'M1U_OM_FEW=',OM1B
      WRITE(*,'(A,I8)')'M1U_NSAM_FEW=',NSB
      RHOM_U(1)=0.05_LDP; RHOP_U(1)=0.06_LDP
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:1),RHOP_U(1:1),NDV_U(1:1), &
        OM0,OM1B,NSB,P90B,MEDB,MXB,1)
      WRITE(*,'(A,ES16.8)')'M1U_OM_NOPROD=',OM1B
      WRITE(*,'(A,I8)')'M1U_NSAM_NOPROD=',NSB
      RHOM_U(1)=99.0_LDP; RHOP_U(1)=99.0_LDP
      CALL KN_OMEGA_POOL_V1(RHOM_U(1:1),RHOP_U(1:1),NDV_U(1:1), &
        OM0,OM1B,NSB,P90B,MEDB,MXB,1)
      WRITE(*,'(A,ES16.8)')'M1U_OM_NOPROD2=',OM1B
!
! --- omega gate helper: status 0/1 admit, status 2 blocks ---
!
      WRITE(*,'(A,L2)')'M1U_OMOK0=',KN_PROD_OMEGA_OK(0)
      WRITE(*,'(A,L2)')'M1U_OMOK1=',KN_PROD_OMEGA_OK(1)
      WRITE(*,'(A,L2)')'M1U_OMOK2=',KN_PROD_OMEGA_OK(2)
      WRITE(*,'(A,I4)')'M1U_IST_NORM=', &
        KN_PROD_POOL_STAT(NDV_U,3)
!
! --- probe predicate ---
!
      CALL SET_FG_NOCP_PROBE(.TRUE.,3)
      L1=KN_NOCP_PROBE_RAY(3)
      L2=KN_NOCP_PROBE_RAY(1)
      CALL SET_FG_NOCP_PROBE(.FALSE.,3)
      L3=KN_NOCP_PROBE_RAY(3)
      WRITE(*,'(A,3L2)')'M1U_PRED=',L1,L2,L3
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_M1UNIT
!
! ----------------------------------------------------------------------
! Phase 2.15 certification classifier: deterministic table over
! KN_CERTIFY covering every class, every reason code, both
! severity folds and the LCERT=F legacy mode. The commit gate in
! KNOP_RAY_SOLVE_V1 reads only ICLS==KN_CS_FULL, so these rows
! pin the exact class the gate sees; 'cerb'/'cer0'/'cern'/'cer2'
! run the same classes end-to-end through the refresh.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_CERT()
      INTEGER I_C,ICLS_T,IREA_T
      INTEGER IX_T(2),ND_T,INVR_T
      REAL(KIND=LDP) RE_T(2),QQ_T(2)
      REAL(KIND=LDP) RQNAN,RQINF,RQZ
      LOGICAL LC_T
      CHARACTER(LEN=8) CN_T
!
      KN_RAY_TOL=2.0E-3_LDP
      KN_RAY_BULK_TOL=6.0E-2_LDP
      RQZ=0.0_LDP
      RQNAN=RQZ/RQZ
      RQINF=1.0_LDP/RQZ
      ND_T=0
!
! c01: both channels strict residual + verified -> FULL.
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=1.9E-4_LDP; RE_T(2)=1.9E-4_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=1.0E-4_LDP; LC_T=.TRUE.
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c01 ',ICLS_T,IREA_T
!
! c02: residual passes, verify residual above tol -> STAG.
      QQ_T(1)=5.0E-3_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c02 ',ICLS_T,IREA_T
!
! c03: bulk floored exits, TOL < RES <= BTOL, valid verify -> BULK.
! (Phase 2.16: under LCERT a bulk class also needs an evaluated,
! finite, non-negative Q -- see c23/c24 for the fail-closed rows.)
      IX_T(1)=2; IX_T(2)=2
      RE_T(1)=8.3E-3_LDP; RE_T(2)=3.0E-2_LDP
      QQ_T(1)=3.0E-2_LDP; QQ_T(2)=4.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c03 ',ICLS_T,IREA_T
!
! c04: same inputs with BTOL=0 -> bulk rule disabled, an IXIT=2
! arriving anyway has no contract -> INV (fail closed).
      KN_RAY_BULK_TOL=0.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c04 ',ICLS_T,IREA_T
      KN_RAY_BULK_TOL=6.0E-2_LDP
!
! c05: pass-cap exit -> ITLIM.
      IX_T(1)=4; IX_T(2)=4
      RE_T(1)=1.0E-2_LDP; RE_T(2)=1.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c05 ',ICLS_T,IREA_T
!
! c06: same-direction stall exit -> STAG.
      IX_T(1)=3; IX_T(2)=3
      RE_T(1)=5.0E-2_LDP; RE_T(2)=5.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c06 ',ICLS_T,IREA_T
!
! c07: NDIV set -> DIV regardless of channel exits.
      IX_T(1)=1; IX_T(2)=1; ND_T=1
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c07 ',ICLS_T,IREA_T
      ND_T=0
!
! c08: non-finite (NaN) M residual -> INV.
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=RQNAN; RE_T(2)=1.0E-4_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c08 ',ICLS_T,IREA_T
!
! c09: non-finite (Inf) P residual -> INV. The finite M
! channel gets a valid verify so the worst-channel reason is
! the P channel's non-finite state, not the missing verify.
      RE_T(1)=1.0E-4_LDP; RE_T(2)=RQINF
      QQ_T(1)=0.0_LDP; QQ_T(2)=0.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c09 ',ICLS_T,IREA_T
!
! c10: nominal exits with the verify never run -> INV.
      RE_T(2)=1.0E-4_LDP
      QQ_T(1)=-1.0_LDP; QQ_T(2)=-1.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c10 ',ICLS_T,IREA_T
!
! c11: verify ran but returned a non-finite value -> INV.
      QQ_T(1)=RQNAN
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c11 ',ICLS_T,IREA_T
!
! c12: severity fold -- M channel FULL, P channel BULK -> BULK.
! (Phase 2.16: the BULK channel's verify must be evaluated.)
      IX_T(1)=1; IX_T(2)=2
      RE_T(1)=1.0E-4_LDP; RE_T(2)=8.3E-3_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=3.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c12 ',ICLS_T,IREA_T
!
! c13: LCERT=F -> legacy contract, no verify required.
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=1.0E-4_LDP; RE_T(2)=1.0E-4_LDP
      QQ_T(2)=-1.0_LDP; LC_T=.FALSE.
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c13 ',ICLS_T,IREA_T
!
! c14: LCERT=F does not rescue a bulk exit either -> BULK.
      IX_T(1)=2; IX_T(2)=2
      RE_T(1)=8.3E-3_LDP; RE_T(2)=8.3E-3_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c14 ',ICLS_T,IREA_T
      LC_T=.TRUE.
!
! c15: bulk code arriving with strict residual met -> FULL when
! the verify also passes (the class follows the residual
! contract, not the exit code alone).
      RE_T(1)=1.0E-4_LDP; RE_T(2)=1.0E-4_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=0.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c15 ',ICLS_T,IREA_T
!
! c16: historical IXIT=0 exact exit still certifies.
      IX_T(1)=0; IX_T(2)=0
      RE_T(1)=0.0_LDP; RE_T(2)=0.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c16 ',ICLS_T,IREA_T
!
! c17: unset exit code -> INV.
      IX_T(1)=-1; IX_T(2)=-1
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c17 ',ICLS_T,IREA_T
!
! c18: plateau exit (IXIT=6) -> STAG.
      IX_T(1)=6; IX_T(2)=6
      RE_T(1)=5.0E-2_LDP; RE_T(2)=5.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c18 ',ICLS_T,IREA_T
!
! c19: severity fold -- BULK channel vs ITLIM channel -> ITLIM.
      IX_T(1)=2; IX_T(2)=4
      RE_T(1)=8.3E-3_LDP; RE_T(2)=1.0E-2_LDP
      QQ_T(1)=3.0E-2_LDP; QQ_T(2)=-1.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c19 ',ICLS_T,IREA_T
!
! c20: explicit channel divergence code -> DIV.
      IX_T(1)=5; IX_T(2)=5
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c20 ',ICLS_T,IREA_T
!
! Phase 2.16 rows: invalid-state separation, negative residuals,
! complete verify validation under LCERT, and the INVR override.
!
      INVR_T=0
!
! c21: negative production residual -> INV (a max norm cannot be
! negative; a negative RES is a broken residual, reason 3).
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=-1.0E-4_LDP; RE_T(2)=1.0E-4_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=0.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c21 ',ICLS_T,IREA_T
!
! c22: negative verify residual on the strict path -> INV(2).
      RE_T(1)=1.0E-4_LDP
      QQ_T(1)=-5.0E-3_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c22 ',ICLS_T,IREA_T
!
! c23: bulk exits but the verify never ran (Q=-1) -> INV, not
! BULK -- the missing-verification sentinel cannot hide behind
! the provisional class (reason 2).
      IX_T(1)=2; IX_T(2)=2
      RE_T(1)=8.3E-3_LDP; RE_T(2)=3.0E-2_LDP
      QQ_T(1)=-1.0_LDP; QQ_T(2)=-1.0_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c23 ',ICLS_T,IREA_T
!
! c24: bulk exits with a non-finite verify -> INV, not BULK.
      QQ_T(1)=RQNAN; QQ_T(2)=3.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c24 ',ICLS_T,IREA_T
!
! c25: bulk exits with an evaluated but above-tol verify -> BULK.
! BULK is already a provisional class; only missing/invalid Q
! fails it closed.
      QQ_T(1)=5.0E-2_LDP; QQ_T(2)=4.0E-2_LDP
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c25 ',ICLS_T,IREA_T
!
! c26: invalid coefficient (INVR=7) plus finite divergence
! (NDIV=1) -> INV wins over DIV.
      IX_T(1)=5; IX_T(2)=5
      RE_T(1)=1.0E-2_LDP; RE_T(2)=1.0E-2_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=0.0_LDP
      ND_T=1; INVR_T=7
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c26 ',ICLS_T,IREA_T
      ND_T=0
!
! c27: invalid frequency span (INVR=8) on nominal exits -> INV.
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=1.0E-4_LDP; RE_T(2)=1.0E-4_LDP
      INVR_T=8
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c27 ',ICLS_T,IREA_T
!
! c28: allocation failure (INVR=9) -> INV.
      INVR_T=9
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T,INVR_T)
      WRITE(*,'(A,2I3)')'CERT c28 ',ICLS_T,IREA_T
      INVR_T=0
!
! c29: NDIV divergence reported without the optional INVR argument
! -> DIV. Pairs with c01 (no INVR, FULL) so the absent-argument
! path is exercised on both sides of the NDIV gate.
      IX_T(1)=1; IX_T(2)=1
      RE_T(1)=1.0E-4_LDP; RE_T(2)=1.0E-4_LDP
      QQ_T(1)=0.0_LDP; QQ_T(2)=0.0_LDP
      ND_T=1
      CALL KN_CERTIFY(IX_T(1),IX_T(2),RE_T(1),RE_T(2), &
      QQ_T(1),QQ_T(2),ND_T,LC_T,ICLS_T,IREA_T)
      WRITE(*,'(A,2I3)')'CERT c29 ',ICLS_T,IREA_T
      ND_T=0
!
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
      END SUBROUTINE DRV_CERT
!
! ----------------------------------------------------------------------
! Maturity read-back (audit Issue 3): calls the production
! GET_FG_GS_MATURITY so the test sees the RF_DCHG the release gate
! in SOLVE_FOR_POPS is actually fed -- committed-or-pending after
! the fold, not the committed-only store delta.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_MATURITY()
      LOGICAL GS_OK_T,GS_ALLC_T
      REAL(KIND=LDP) GS_CHG_T,GS_MAX_T,GS_SWM_T,GS_RFD_T
      INTEGER GS_CAPF_T,GS_CAPP_T,GS_NDIV_T,GS_NRAY_T
      INTEGER*8 GS_NBAD_T,GS_NTOT_T,GS_NCONV_T,GS_NBLK_T,GS_SWN_T
      INTEGER*8 GS_NSTL_T,GS_NPLT_T
      INTEGER*8 GS_SWC_T,GS_SWT_T,GS_RFC_T
!
      CALL GET_FG_GS_MATURITY(GS_OK_T,GS_CHG_T,GS_CAPF_T, &
      GS_MAX_T,GS_CAPP_T,GS_NBAD_T,GS_NTOT_T,GS_NCONV_T, &
      GS_NDIV_T,GS_NRAY_T,GS_NBLK_T, &
      GS_SWM_T,GS_SWN_T,GS_RFD_T, &
      GS_NSTL_T,GS_NPLT_T, &
      GS_SWC_T,GS_SWT_T,GS_RFC_T,GS_ALLC_T)
      WRITE(*,'(A,L2)')'GS_OK=',GS_OK_T
      WRITE(*,'(A,ES16.8)')'GS_RFD=',GS_RFD_T
      WRITE(*,'(A,L2)')'GS_ALLCONV=',GS_ALLC_T
      RETURN
      END SUBROUTINE DRV_MATURITY
!
! ----------------------------------------------------------------------
! P1-T1g repeatability probe (diagnostic only). Snapshots the complete
! refresh input state (lagged stores + per-ray seed flags), runs the
! production FG_KNOP_REFRESH_V1, restores the snapshot, re-runs, and
! compares every committed store element plus the per-ray commit
! verdicts. A bit-identical second pass proves the refresh is a
! deterministic function of its captured inputs on this configuration;
! drift identifies state a store-only snapshot misses. KN_HOT_SPENT is
! recorded after each pass: it is the shared OMP hot-pass budget and
! the only refresh input that can couple rays through scheduling.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_REPL()
      IMPLICIT NONE
!
      INTEGER LS,K,L,IOSR,NX_DIF
      INTEGER, ALLOCATABLE :: SEE_IN(:),VC1(:),VC2(:)
      REAL(KIND=LDP), ALLOCATABLE :: SP_IN(:,:,:),SM_IN(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: SP_O1(:,:,:),SM_O1(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: SP_O2(:,:,:),SM_O2(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RPMX(:)
      REAL(KIND=LDP) DMAXA,DMAXR,RMS,FRQ0,FRQ1,V1,V2,FD,DN,RES1
      INTEGER NMIS_P,NMIS_M
      INTEGER*8 HSP1
!
      INTEGER OMP_GET_MAX_THREADS
      EXTERNAL OMP_GET_MAX_THREADS
!
! Refresh entry gates identical to the production post-sweep path.
!
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
      HAS_COMPRESSING=.TRUE.
      KN_STORE_READY=.TRUE.
      FREQ_IDX_ITER=KN_NCF-1
      SOLUTION_METHOD='INTEGRAL'
!
! Arm the per-ray seed-fill state the production path maintains; the
! jump veto engages only on rays already seeded, so seed them all.
!
      IF(ALLOCATED(FG_RAY_SEEDED))DEALLOCATE(FG_RAY_SEEDED)
      ALLOCATE(FG_RAY_SEEDED(KN_NP),STAT=IOSR)
      IF(IOSR .EQ. 0)FG_RAY_SEEDED=1
!
      ALLOCATE(SP_IN(N_COMP,KN_NP,KN_NCF),SM_IN(N_COMP,KN_NP,KN_NCF), &
               SP_O1(N_COMP,KN_NP,KN_NCF),SM_O1(N_COMP,KN_NP,KN_NCF), &
               SP_O2(N_COMP,KN_NP,KN_NCF),SM_O2(N_COMP,KN_NP,KN_NCF), &
               SEE_IN(KN_NP),VC1(KN_NP),VC2(KN_NP),RPMX(KN_NP), &
               STAT=IOSR)
      IF(IOSR .NE. 0)THEN
        WRITE(*,'(A)')'REPL_FAIL alloc'
        RETURN
      END IF
      VC1=-99
      VC2=-99
      SP_IN=I_P_DEPTH_NEXT
      SM_IN=I_M_DEPTH_NEXT
      SEE_IN=FG_RAY_SEEDED
!
! Input fingerprint over the frozen KN_* opacity, the frequency grid
! and the ray geometry. Recomputed before pass 2 to prove the kernel
! left its inputs untouched; the GAM_MAG/dGAMdR sanitisation in
! BUILD_FREQ is the only in-place geometry write and must be a no-op
! on finite manufactured data.
!
      FRQ0=0.0_LDP
      DO L=1,KN_NCF
        FRQ0=FRQ0+KN_FREQ(L)+KN_DLN_P(L)+KN_DLN_N(L)+KN_DBB(L)+KN_IC(L)
      END DO
      DO K=1,SIZE(KN_CHI_EXT,1)
        DO L=1,KN_NCF
          FRQ0=FRQ0+DBLE(KN_CHI_EXT(K,L))+DBLE(KN_ETA_EXT(K,L))
        END DO
      END DO
      DO LS=1,KN_NP
        DO K=1,SIZE(GAM_MAG,1)
          FRQ0=FRQ0+GAM_MAG(K,LS)+dGAMdR(K,LS)+GAM(K,LS)
        END DO
      END DO
!
! ---- refresh pass 1 ----
!
      CALL FG_KNOP_REFRESH_V1()
      SP_O1=I_P_DEPTH_NEXT
      SM_O1=I_M_DEPTH_NEXT
      RES1=KN_RF_RESID
      HSP1=KN_HOT_SPENT
      IF(ALLOCATED(FG_RAY_VC))VC1=FG_RAY_VC
      WRITE(*,'(A,I4)')'REPL_OMP=',OMP_GET_MAX_THREADS()
      WRITE(*,'(A,I8)')'REPL_HOT1=',HSP1
!
! ---- restore and refresh pass 2 ----
!
      I_P_DEPTH_NEXT=SP_IN
      I_M_DEPTH_NEXT=SM_IN
      FG_RAY_SEEDED=SEE_IN
!
      FRQ1=0.0_LDP
      DO L=1,KN_NCF
        FRQ1=FRQ1+KN_FREQ(L)+KN_DLN_P(L)+KN_DLN_N(L)+KN_DBB(L)+KN_IC(L)
      END DO
      DO K=1,SIZE(KN_CHI_EXT,1)
        DO L=1,KN_NCF
          FRQ1=FRQ1+DBLE(KN_CHI_EXT(K,L))+DBLE(KN_ETA_EXT(K,L))
        END DO
      END DO
      DO LS=1,KN_NP
        DO K=1,SIZE(GAM_MAG,1)
          FRQ1=FRQ1+GAM_MAG(K,LS)+dGAMdR(K,LS)+GAM(K,LS)
        END DO
      END DO
      WRITE(*,'(A,ES16.8)')'REPL_FRQFP=',ABS(FRQ1-FRQ0)
!
      CALL FG_KNOP_REFRESH_V1()
      SP_O2=I_P_DEPTH_NEXT
      SM_O2=I_M_DEPTH_NEXT
      IF(ALLOCATED(FG_RAY_VC))VC2=FG_RAY_VC
!
! ---- element-wise compare over the writable store slots ----
! Store slot s holds frequency s+1; slot KN_NCF is never written and
! stays at its seed, so it is excluded as a constant boundary input.
!
      DMAXA=0.0_LDP
      DMAXR=0.0_LDP
      RMS=0.0_LDP
      DN=0.0_LDP
      NMIS_P=0
      NMIS_M=0
      RPMX=0.0_LDP
      DO LS=1,KN_NP
        DO K=1,N_COMP
          DO L=1,KN_NCF-1
            V1=SP_O1(K,LS,L)
            V2=SP_O2(K,LS,L)
            FD=ABS(V1-V2)
            DMAXA=MAX(DMAXA,FD)
            IF(ABS(V1) .GT. 1.0E-30_LDP) &
               DMAXR=MAX(DMAXR,FD/MAX(ABS(V1),1.0E-30_LDP))
            RMS=RMS+FD*FD
            DN=DN+1.0_LDP
            RPMX(LS)=MAX(RPMX(LS),FD)
            IF(FD .GT. 1.0E-12_LDP*MAX(1.0_LDP,ABS(V1)))NMIS_P=NMIS_P+1
            V1=SM_O1(K,LS,L)
            V2=SM_O2(K,LS,L)
            FD=ABS(V1-V2)
            DMAXA=MAX(DMAXA,FD)
            IF(ABS(V1) .GT. 1.0E-30_LDP) &
               DMAXR=MAX(DMAXR,FD/MAX(ABS(V1),1.0E-30_LDP))
            RMS=RMS+FD*FD
            DN=DN+1.0_LDP
            RPMX(LS)=MAX(RPMX(LS),FD)
            IF(FD .GT. 1.0E-12_LDP*MAX(1.0_LDP,ABS(V1)))NMIS_M=NMIS_M+1
          END DO
        END DO
      END DO
      IF(DN .GT. 0.0_LDP)RMS=SQRT(RMS/DN)
!
      NX_DIF=0
      IF(ALLOCATED(VC1) .AND. ALLOCATED(FG_RAY_VC) .AND. &
         SIZE(VC1) .EQ. SIZE(FG_RAY_VC))NX_DIF=COUNT(VC1 .NE. VC2)
!
      WRITE(*,'(A,ES16.8)')'REPL_MAXABS=',DMAXA
      WRITE(*,'(A,ES16.8)')'REPL_MAXREL=',DMAXR
      WRITE(*,'(A,ES16.8)')'REPL_RMS=',RMS
      WRITE(*,'(A,I8)')'REPL_NMIS_P=',NMIS_P
      WRITE(*,'(A,I8)')'REPL_NMIS_M=',NMIS_M
      WRITE(*,'(A,I8)')'REPL_VCDIF=',NX_DIF
      WRITE(*,'(A,I8)')'REPL_HOTSPENT=',KN_HOT_SPENT
      WRITE(*,'(A,ES16.8)')'REPL_RES1=',RES1
      WRITE(*,'(A,ES16.8)')'REPL_RES2=',KN_RF_RESID
      DO LS=1,KN_NP
        WRITE(*,'(A,I3,A,ES16.8)')'REPL_RAYMX ',LS,' ',RPMX(LS)
      END DO
!
      DEALLOCATE(SP_IN,SM_IN,SP_O1,SM_O1,SP_O2,SM_O2, &
                 SEE_IN,VC1,VC2,RPMX)
      RETURN
      END SUBROUTINE DRV_REPL
!
      END

! ----------------------------------------------------------------------
! P1-T1b global radiation-map Jacobian (diagnostic only). One map
! application F(X) is the complete production ordered radiation
! iteration on the same manufactured |Q|->1 model as DRV_P1T1: the
! frequency sweep (INIT at L=1, NEW_FREQ=T per frequency, the normal
! in-sweep relaxed I_NEXT store update with FREEZE_I_NEXT=F) followed
! by the production post-sweep FG_KNOP_REFRESH_V1. The state vector X
! is the set of lagged-store elements the map actually writes:
! I_P/I_M_DEPTH_NEXT at compressing (KC,ray,slot) for slot=1..NCF-1.
! I_*_PREV/I_*_STORE are slaved, not independent state: INIT zeroes
! PREV, each NEW_FREQ copies STORE into PREV, and STORE is rewritten
! at the end of every frequency call before the next NEW_FREQ reads
! it, so no previous-sweep value of either array can reach the next
! sweep. Store slot NCF is never written by the sweep or the refresh
! and stays identically zero, so it is excluded as a constant
! boundary input rather than dynamic state.
!
! Jacobian variants: IV=1 full map (sweep+refresh, omega=0.5),
! IV=2 sweep only (omega=0.5), IV=3 sweep only (omega=1.0) for the
! (1-w)I + w J_T algebraic cross-check. Dense central finite
! differences at two amplitudes; LAPACK dgeev for the spectrum;
! direct nonlinear iteration for the transient-growth check.
! ----------------------------------------------------------------------
      SUBROUTINE DRV_P1T1B
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=8,NPA=2,NCFA=6)
      REAL(KIND=LDP) RA(NDA),VA(NDA),SGA(NDA),SGT(NDA),PA(NPA)
      REAL(KIND=LDP) ETA_M(NDA,NCFA),CHI_M(NDA,NCFA),ESC_M(NDA,NCFA)
      REAL(KIND=LDP) JNU_A(NDA),FEDD_A(NDA),IPLUS_A(NPA)
      REAL(KIND=LDP) JQW_A(NDA,NPA),HQW_A(NDA,NPA),KQW_A(NDA,NPA)
      REAL(KIND=LDP) NQW_A(NDA,NPA),HMID_A(NDA,NPA),NMID_A(NDA,NPA)
      REAL(KIND=LDP) VDOP_A(NDA),FRQV(NCFA),INBC_A,HBC_A
      REAL(KIND=LDP) LVR(NDA),BUMP,V0,DLOGA,DBB_A,IC_A
      REAL(KIND=LDP) VDOP_F,REXT_A,DEPF
      CHARACTER*6 METH_A
      CHARACTER(LEN=8) SOPT_A
      CHARACTER(LEN=12) IBND_A
      LOGICAL THK_A,INIT_A,NEWF_A
      INTEGER NC_A,NITER,ITR,L,I,LSDBG,IOS,NRAYS
      REAL(KIND=LDP) SMAX
!
! P1-T1b map/Jacobian work arrays. POSP/POSM map (KC,slot) to the X
! index or 0; WKCH lists the writable compressing slots. XPR/XFP/JW
! are FD scratch; the SP*/SM* arrays snapshot the slaved
! I_*_STORE/I_*_PREV donors so every MAP_F call starts from identical
! bookkeeping state.
!
      INTEGER NCOMP,NSNP,NSLOT,NW,NX,LWK8
      INTEGER, ALLOCATABLE :: POSP(:,:),POSM(:,:),WKCH(:),FDIF(:)
      REAL(KIND=LDP), ALLOCATABLE :: XO(:),XNEW(:),X0V(:),FX0(:)
      REAL(KIND=LDP), ALLOCATABLE :: XPR(:),XFP(:),JW(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: STO_P(:,:,:),STO_M(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: SPSP(:,:),SMSP(:,:),SMPP(:,:), &
      SMMP(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: J1(:,:),J2(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: XST(:),DXP(:),JN(:,:),ANW(:,:)
      INTEGER, ALLOCATABLE :: IPI8(:)
      INTEGER NWIT
      REAL(KIND=LDP) RDX
      REAL(KIND=LDP), ALLOCATABLE :: RESV(:)
      REAL(KIND=LDP), ALLOCATABLE :: A8(:,:),VR8(:,:),WK8(:),WR8(:), &
      WI8(:),VL8(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: VDOM(:),VPT(:)
      REAL(KIND=LDP) EPSV,AMPF,XFLOOR,RN,RD,TNRM,EMOD,RNPRV
      INTEGER KC,SL,J,IV,IA,INFO8,IMX,JP,K
      INTEGER*8 NEVAL
      INTEGER IRFI(3)
      REAL(KIND=LDP) OMV(3)
      CHARACTER(LEN=8) VNAME(3)
      CHARACTER(LEN=32) CARG
!
! CNAME is not host-associated into this routine, so read the case
! name locally; 'p1t1bf' selects the FREEZE_I_PREV ablation below.
      CALL GETARG(1,CARG)
!
      LSDBG=1
      NRAYS=1
      RA=(/8.0E4_LDP,6.4E4_LDP,5.0E4_LDP,3.9E4_LDP, &
      3.0E4_LDP,2.2E4_LDP,1.5E4_LDP,1.0E4_LDP/)
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
      DO L=1,NCFA
        DO I=1,NDA
          CHI_M(I,L)=1.0E-3_LDP*(1.0_LDP+0.05_LDP*L)* &
      (RA(1)/RA(I))**2
          ETA_M(I,L)=2.0E-6_LDP*(1.0_LDP+0.02_LDP*L)*(RA(1)/RA(I))
          ESC_M(I,L)=0.0_LDP
        END DO
      END DO
!
      DEPF=1.0E-1_LDP
      CHI_M(2,:)=CHI_M(2,:)*DEPF
      CHI_M(3,:)=CHI_M(3,:)*1.0E-3_LDP
      CHI_M(4,:)=CHI_M(4,:)*1.0E-5_LDP
      CHI_M(5,:)=CHI_M(5,:)*1.0E-7_LDP
!
      FRQV(1)=6.0E15_LDP
      DO L=2,NCFA
        FRQV(L)=FRQV(L-1)*0.75_LDP
      END DO
      DLOGA=LOG(1.0_LDP/0.75_LDP)
!
      PA(1)=0.0_LDP
      PA(2)=1.2E4_LDP
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
      FG_FREQ_MATCH=.TRUE.
      FG_ORD_SWEEP=.TRUE.
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
!
! P1-T1b: suppress the INIT-time compression-set hysteresis rebuild.
! SIGMA_TRUE_IN/V are identical on every MAP_F call, so the compressing
! set cannot physically change; without this the check misfires at the
! set boundary and NEW_R_GRID deallocates+zeroes I_*_DEPTH_NEXT each
! sweep, leaving the lagged coupling dead (F independent of X).
!
      KNOP_DIAG_NO_REGRID=.TRUE.
! P1-T1e: 'p1t1bf' reruns the identical Jacobian with FREEZE_I_PREV=T
! so the expanding I_PREV frequency-feedback channel is ablated
! (I_PREV stays at its INIT-zeroed value). Comparing the same-slot
! entries against the normal 'p1t1b' run isolates whether those edges
! are I_PREV-mediated. Diagnostic only; production default is .FALSE.
      FREEZE_I_PREV=(CARG .EQ. 'p1t1bf')
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
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
      KN_DBG_LS=LSDBG
      ALLOCATE(KN_DBG_CHIH(NDA,NCFA),KN_DBG_Q(NDA,NCFA), &
      KN_DBG_DCHR(NDA,NCFA),KN_DBG_DTAU(NDA,NCFA), &
      KN_DBG_SPM(NDA,NCFA),KN_DBG_SPP(NDA,NCFA), &
      KN_DBG_WX(NDA,NCFA),KN_DBG_DLNM(NDA,NCFA), &
      KN_DBG_IMPR(NDA,NCFA),KN_DBG_IPPR(NDA,NCFA), &
      KN_DBG_NM(NDA,NCFA),KN_DBG_KC1(NDA,NCFA),KN_DBG_KC2(NDA,NCFA), &
      KN_DBG_A0(NDA,NCFA),KN_DBG_A1(NDA,NCFA),KN_DBG_A2(NDA,NCFA), &
      KN_DBG_A3(NDA,NCFA),KN_DBG_A4(NDA,NCFA),KN_DBG_SR(NDA,NCFA), &
      KN_DBG_INEXTM(NDA,NCFA),KN_DBG_INEXTP(NDA,NCFA), &
      KN_DBG_NDON(NDA,NCFA), &
      KR_DBG_CHIH(NDA,NCFA),KR_DBG_Q(NDA,NCFA), &
      KR_DBG_DCHR(NDA,NCFA),KR_DBG_DTAU(NDA,NCFA), &
      KR_DBG_WX(NDA,NCFA),KR_DBG_DLNM(NDA,NCFA), &
      KR_DBG_SPM(NDA,NCFA),KR_DBG_SPP(NDA,NCFA), &
      KR_DBG_A0(NDA,NCFA),KR_DBG_A1(NDA,NCFA), &
      KR_DBG_A2(NDA,NCFA),KR_DBG_A3(NDA,NCFA),KR_DBG_A4(NDA,NCFA), &
      KR_DBG_FLDM(NDA,NCFA),KR_DBG_FLDP(NDA,NCFA), &
      KR_DBG_SR(NDA,NCFA), &
      KR_DBG_NFM(NDA,NCFA),KR_DBG_IKC1(NDA,NCFA), &
      KR_DBG_IKC2(NDA,NCFA),KR_DBG_JPR1(NDA,NCFA), &
      KR_DBG_JPR2(NDA,NCFA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1balloc_fail=1'
        RETURN
      END IF
      KN_DBG_CHIH=-1.0_LDP; KN_DBG_Q=-1.0_LDP; KN_DBG_DCHR=-1.0_LDP
      KN_DBG_DTAU=-1.0_LDP; KN_DBG_SPM=-1.0_LDP; KN_DBG_SPP=-1.0_LDP
      KN_DBG_WX=-1.0_LDP; KN_DBG_DLNM=-1.0_LDP
      KN_DBG_IMPR=-1.0_LDP; KN_DBG_IPPR=-1.0_LDP
      KN_DBG_NM=-1; KN_DBG_KC1=-1; KN_DBG_KC2=-1
      KN_DBG_A0=-1.0_LDP; KN_DBG_A1=-1.0_LDP; KN_DBG_A2=-1.0_LDP
      KN_DBG_A3=-1.0_LDP; KN_DBG_A4=-1.0_LDP; KN_DBG_SR=-1.0_LDP
      KN_DBG_INEXTM=-1.0_LDP; KN_DBG_INEXTP=-1.0_LDP
      KN_DBG_NDON=0
      KR_DBG_NFM=-1; KR_DBG_IKC1=-1; KR_DBG_IKC2=-1
      KR_DBG_JPR1=-1; KR_DBG_JPR2=-1
      KR_DBG_CHIH=-1.0_LDP; KR_DBG_Q=-1.0_LDP; KR_DBG_DCHR=-1.0_LDP
      KR_DBG_DTAU=-1.0_LDP; KR_DBG_WX=-1.0_LDP; KR_DBG_DLNM=-1.0_LDP
      KR_DBG_SPM=-1.0_LDP; KR_DBG_SPP=-1.0_LDP
      KR_DBG_A0=-1.0_LDP; KR_DBG_A1=-1.0_LDP; KR_DBG_A2=-1.0_LDP
      KR_DBG_A3=-1.0_LDP; KR_DBG_A4=-1.0_LDP; KR_DBG_SR=-1.0_LDP
      KR_DBG_FLDM=-1.0_LDP; KR_DBG_FLDP=-1.0_LDP
!
! Deterministic solver inputs. The QW arrays are read-only inputs
! that feed only the diagnostic J-moment outputs (JNU_STORE -> FEDD,
! JPLUS_*, LAMSTAR_DIAG), never the intensity/store path, so unit
! weights are safe and keep JNU_STORE nonzero for the FEDD quotient.
! The JNU/FEDD/IPLUS/HBC dummies are solver outputs re-initialised
! so no previous evaluation's values can leak into the next call.
!
      JQW_A=1.0_LDP; HQW_A=1.0_LDP; KQW_A=1.0_LDP
      NQW_A=1.0_LDP; HMID_A=1.0_LDP; NMID_A=1.0_LDP
      JNU_A=0.0_LDP; FEDD_A=0.0_LDP; IPLUS_A=0.0_LDP
      INBC_A=0.0_LDP; HBC_A=0.0_LDP
!
! Cold start: one full production map application from the zeroed
! store so the coupling arrays are allocated and populated. This is
! iteration 1 of the fixed-point sequence X^(n+1)=F(X^(n)).
!
      DO L=1,NCFA
        IF(L .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)*0.75_LDP)
        END IF
        INIT_A=(L .EQ. 1)
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L),ESC_M(1,L), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA, &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
      END DO
      CALL FG_KNOP_REFRESH_V1()
!
! Snapshot the post-cold store as the inactive-element baseline and
! build the active-state index. The store is (N_COMP,NP,NCF_IN); this
! driver runs one ray (NP=1), so only LS=1 is indexed.
!
      IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT))THEN
        WRITE(*,'(A)')'RESULT p1t1bstore_missing=1'
        RETURN
      END IF
      NCOMP=SIZE(I_P_DEPTH_NEXT,1)
      NSNP=SIZE(I_P_DEPTH_NEXT,2)
      NSLOT=SIZE(I_P_DEPTH_NEXT,3)
      IF(NSNP .NE. 1)THEN
        WRITE(*,'(A)')'RESULT p1t1bnp_unsupported=1'
        RETURN
      END IF
      ALLOCATE(STO_P(NCOMP,NSNP,NSLOT),STO_M(NCOMP,NSNP,NSLOT), &
      POSP(NCOMP,NSLOT),POSM(NCOMP,NSLOT),WKCH(NCOMP),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1bidx_alloc_fail=1'
        RETURN
      END IF
      STO_P=I_P_DEPTH_NEXT
      STO_M=I_M_DEPTH_NEXT
!
! Writable (KC,slot) set: the sweep writes slot FIDX for FIDX=1..NCF-1
! at every compressing physical depth K whose J_PNT(K,1) lands on a
! valid ray point -- exactly the write gate in the store-update block.
! A compressing depth failing that gate is never written, so its slot
! is a constant input, not state; including it would inject spurious
! unit eigenvalues into J.
!
      NW=0
      POSP=0; POSM=0
      DO K=1,NDA
        KC=COMP_DEPTH_IDX(K)
        IF(KC .GE. 1)THEN
          I=J_PNT(K,1)
          IF(I .GE. 1 .AND. I .LE. NI_RAY(1))THEN
            NW=NW+1
            WKCH(NW)=KC
          END IF
        END IF
      END DO
!
! Active slots are 1..NSLOT-1: slot NSLOT is never written by the
! sweep or the refresh and stays at its zero initial value.
!
      NX=2*NW*(NSLOT-1)
      J=0
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          POSP(WKCH(I),SL)=J
        END DO
      END DO
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          POSM(WKCH(I),SL)=J
        END DO
      END DO
      IF(J .NE. NX)THEN
        WRITE(*,'(A)')'RESULT p1t1bidx_mismatch=1'
        RETURN
      END IF
      ALLOCATE(XO(NX),XNEW(NX),X0V(NX),FX0(NX),XPR(NX),XFP(NX), &
      RESV(64),J1(NX,NX),J2(NX,NX),JW(NX,NX),VDOM(NX),VPT(NX),FDIF(NX), &
      SPSP(SIZE(I_P_STORE,1),SIZE(I_P_STORE,2)), &
      SMSP(SIZE(I_M_STORE,1),SIZE(I_M_STORE,2)), &
      SMPP(SIZE(I_P_PREV,1),SIZE(I_P_PREV,2)), &
      SMMP(SIZE(I_M_PREV,1),SIZE(I_M_PREV,2)),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1bx_alloc_fail=1'
        RETURN
      END IF
      SPSP=I_P_STORE
      SMSP=I_M_STORE
      SMPP=I_P_PREV
      SMMP=I_M_PREV
!
! Dump the state indexing for the python side.
!
      OPEN(91,FILE='p1t1b_idx.txt',STATUS='REPLACE')
      WRITE(91,'(A,4I8)')'P1BIDX ',NX,NW,NSLOT,NCOMP
      J=0
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I8,1X,A1,1X,2I6)')'ELEM',J,'P',WKCH(I),SL
        END DO
      END DO
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I8,1X,A1,1X,2I6)')'ELEM',J,'M',WKCH(I),SL
        END DO
      END DO
      CLOSE(91)
!
! Warm-up: iterate the full production map from the cold-sweep state.
! Record ||F(X)-X|| per application; X0 is the final iterate.
!
      NWARM=30
      CALL PACK_X(X0V)
      NEVAL=0
      DO ITR=2,NWARM
        CALL MAP_F(X0V,XNEW,1)
        NEVAL=NEVAL+1
        RN=0.0_LDP
        DO J=1,NX
          RN=RN+(XNEW(J)-X0V(J))**2
        END DO
        RN=SQRT(RN)
        IF(ITR .LE. 64)RESV(ITR)=RN
        X0V=XNEW
      END DO
!
! Baseline repeatability: two evaluations at X0 must agree.
!
      CALL MAP_F(X0V,XO,1)
      NEVAL=NEVAL+1
      CALL MAP_F(X0V,XNEW,1)
      NEVAL=NEVAL+1
      RD=0.0_LDP; TNRM=0.0_LDP
      DO J=1,NX
        RD=MAX(RD,ABS(XNEW(J)-XO(J)))
        TNRM=MAX(TNRM,ABS(XO(J)))
      END DO
      WRITE(*,'(A,ES16.8)')'P1B_REP_MAXDIFF=',RD
      WRITE(*,'(A,ES16.8)')'P1B_REP_XNORM=',TNRM
      WRITE(*,'(A,ES16.8)')'P1B_REP_REL=',RD/MAX(TNRM,1.0E-30_LDP)
      FX0=XO
!
      OPEN(92,FILE='p1t1b_resid.txt',STATUS='REPLACE')
      WRITE(92,'(A,I4,A,ES16.8)')'WARMUP_ITER NITER_EVALS',NWARM-1
      DO ITR=2,MIN(NWARM,64)
        WRITE(92,'(A,I4,A,ES16.8)')'RES ',ITR,' ',RESV(ITR)
      END DO
      WRITE(92,'(A,ES16.8)')'X0_NORM ',TNRM
      WRITE(92,'(A,ES16.8)')'FX0_NORM ',MAXVAL(ABS(FX0))
      WRITE(92,'(A,ES16.8)')'REP_MAXDIFF ',RD
      WRITE(92,'(A,ES16.8)')'REP_REL ',RD/MAX(TNRM,1.0E-30_LDP)
      CLOSE(92)
      OPEN(93,FILE='p1t1b_X0.txt',STATUS='REPLACE')
      DO J=1,NX
        WRITE(93,'(I8,ES25.16)')J,X0V(J)
      END DO
      CLOSE(93)
!
! Dump the warm-store KN_DBG/KR_DBG coefficient state for the
! boundary-incidence analysis (geometry is store-independent, so the
! post-warm-up capture is representative).
!
      DO I=1,NI_RAY(LSDBG)
        DO L=1,NCFA
          WRITE(*,'(A,I2,A,I2,A,11ES16.8)')'DBS ',I,' ',L,' ', &
      KN_DBG_CHIH(I,L),KN_DBG_Q(I,L),KN_DBG_DCHR(I,L), &
      KN_DBG_DTAU(I,L),KN_DBG_SPM(I,L),KN_DBG_SPP(I,L), &
      KN_DBG_WX(I,L),KN_DBG_DLNM(I,L),KN_DBG_IMPR(I,L), &
      KN_DBG_IPPR(I,L),KN_DBG_SR(I,L)
          WRITE(*,'(A,I2,A,I2,A,3I6)')'DBSI ',I,' ',L,' ', &
      KN_DBG_NM(I,L),KN_DBG_KC1(I,L),KN_DBG_KC2(I,L)
          WRITE(*,'(A,I2,A,I2,A,2ES16.8)')'DBSN ',I,' ',L,' ', &
      KN_DBG_INEXTM(I,L),KN_DBG_INEXTP(I,L)
          WRITE(*,'(A,I2,A,I2,A,5I6)')'DBRI ',I,' ',L,' ', &
      KR_DBG_NFM(I,L),KR_DBG_IKC1(I,L),KR_DBG_IKC2(I,L), &
      KR_DBG_JPR1(I,L),KR_DBG_JPR2(I,L)
        END DO
      END DO
      DO K=1,NDA
        WRITE(*,'(A,I2,A,2I6)')'JMAP ',K,' ',J_PNT(K,LSDBG), &
      COMP_DEPTH_IDX(K)
      END DO
!
! Dense Jacobian for the three variants. Two amplitudes per column;
! central difference where X0_j-eps stays positive, forward where the
! positivity floor forbids the minus side.
!
      IRFI(1)=1; IRFI(2)=0; IRFI(3)=0
      OMV(1)=0.5_LDP; OMV(2)=0.5_LDP; OMV(3)=1.0_LDP
      VNAME(1)='fullw05'; VNAME(2)='swp_w05'; VNAME(3)='swp_w10'
      XFLOOR=1.0E-3_LDP
      LWK8=8*NX
!
      DO IV=1,3
        IF(CARG .EQ. 'p1t1bf' .AND. IV .NE. 2) CYCLE
        INEXT_RELAX=OMV(IV)
        CALL MAP_F(X0V,FX0,IRFI(IV))
        NEVAL=NEVAL+1
!
! Two amplitudes: 1E-4 and 1E-5 times the per-element scale floor.
! JW holds the current amplitude's Jacobian; J1/J2 keep both for the
! consistency check and the dgeev call below.
!
        DO IA=1,2
          AMPF=10.0_LDP**(-(IA+3))
          JW=0.0_LDP
          DO J=1,NX
            EPSV=AMPF*MAX(ABS(X0V(J)),XFLOOR)
            XPR=X0V
            XPR(J)=XPR(J)+EPSV
            CALL MAP_F(XPR,XFP,IRFI(IV))
            NEVAL=NEVAL+1
            IF(X0V(J)-EPSV .GE. 0.0_LDP)THEN
              FDIF(J)=0
              XPR=X0V
              XPR(J)=XPR(J)-EPSV
              CALL MAP_F(XPR,XO,IRFI(IV))
              NEVAL=NEVAL+1
              DO JP=1,NX
                JW(JP,J)=(XFP(JP)-XO(JP))/(2.0_LDP*EPSV)
              END DO
            ELSE
!
! Forward difference where the minus perturbation would take the
! store element negative; FDIF flags the column for the python side.
!
              FDIF(J)=1
              DO JP=1,NX
                JW(JP,J)=(XFP(JP)-FX0(JP))/EPSV
              END DO
            END IF
          END DO
          IF(IA .EQ. 1)THEN
            J1=JW
          ELSE
            J2=JW
          END IF
!
! Dump this variant/amp Jacobian (list-directed; python parses the
! whitespace-separated rows).
!
          OPEN(94,FILE='p1t1b_J'//TRIM(VNAME(IV))//'_a'// &
          ACHAR(48+IA)//'.txt',STATUS='REPLACE')
          DO JP=1,NX
            WRITE(94,*) (JW(JP,J),J=1,NX)
          END DO
          CLOSE(94)
        END DO
!
! Column consistency between the two amplitudes for this variant:
! max elementwise |J1-J2|, max elementwise |J|, their ratio, the
! column max of J2, and the forward-difference flag.
!
        OPEN(95,FILE='p1t1b_cons_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        DO J=1,NX
          RN=0.0_LDP; RD=0.0_LDP
          DO JP=1,NX
            RN=MAX(RN,ABS(J1(JP,J)-J2(JP,J)))
            RD=MAX(RD,MAX(ABS(J1(JP,J)),ABS(J2(JP,J))))
          END DO
          WRITE(95,'(I8,4ES16.8,1X,I2)')J,RN,RD, &
      RN/MAX(RD,1.0E-30_LDP),MAXVAL(ABS(J2(:,J))),FDIF(J)
        END DO
        CLOSE(95)
!
! LAPACK dgeev on the finer-amplitude Jacobian (LDP is real*8 kind).
! Eigenvalues dump in natural order; the python side sorts by
! modulus and forms residuals against the dumped J.
!
        ALLOCATE(A8(NX,NX),VR8(NX,NX),WK8(LWK8),WR8(NX),WI8(NX), &
        VL8(1,NX),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(*,'(A)')'RESULT p1t1bdgeev_alloc=1'
          RETURN
        END IF
        A8=J2
        CALL DGEEV('N','V',NX,A8,NX,WR8,WI8,VL8,1,VR8,NX,WK8,LWK8, &
        INFO8)
        WRITE(*,'(A,I2,A,I8)')'P1B_DGEEV IV INFO NEL',IV,INFO8,NX
        OPEN(96,FILE='p1t1b_eig_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        DO J=1,NX
          EMOD=SQRT(WR8(J)**2+WI8(J)**2)
          WRITE(96,'(I8,3ES25.16)')J,WR8(J),WI8(J),EMOD
        END DO
        CLOSE(96)
!
! Dominant right eigenvector of the full-map Jacobian feeds the
! nonlinear-iteration perturbation. For a complex pair the stored VR
! column is the real part of the pair; documented approximation.
!
        IF(IV .EQ. 1)THEN
          IMX=1
          EMOD=-1.0_LDP
          DO J=1,NX
            RD=SQRT(WR8(J)**2+WI8(J)**2)
            IF(RD .GT. EMOD)THEN
              EMOD=RD; IMX=J
            END IF
          END DO
          TNRM=0.0_LDP
          DO J=1,NX
            VDOM(J)=VR8(J,IMX)
            TNRM=TNRM+VDOM(J)**2
          END DO
          TNRM=SQRT(TNRM)
          IF(TNRM .GT. 0.0_LDP)VDOM=VDOM/TNRM
          WRITE(*,'(A,I4,A,3ES16.8)')'P1B_DOMEIG K WR WI MOD',IMX, &
      WR8(IMX),WI8(IMX),SQRT(WR8(IMX)**2+WI8(IMX)**2)
        END IF
        DEALLOCATE(A8,VR8,WK8,WR8,WI8,VL8)
!
! p1t1n: Newton solve of the implemented fixed-point equation
! F(X)=X on the packed compressing store, run per variant inside
! the IV loop so each map is solved with its own J2 and its own
! refresh/relaxation setting (IRFI(IV), OMV(IV)). The three
! variants isolate the physics: fullw05 packs the post-refresh
! store, so its map already contains the inner fixed-point solve
! and comes out locally constant; swp_w05/swp_w10 pack the
! post-sweep store and are the raw lagged inner map whose spectral
! radius the refresh exists to overcome. For an affine map one
! Newton step is the exact fixed point; ||F(X*)-X*|| then reports
! whether a stable fixed point exists and whether the map is
! affine to FD accuracy. Step 0 reuses J2 (Jacobian at X0V);
! later steps rebuild the FD Jacobian at the current iterate.
! Diagnostic only - the production solver is untouched and only
! this driver's own store is perturbed.
!
        IF(CARG .EQ. 'p1t1n')THEN
          ALLOCATE(XST(NX),DXP(NX),IPI8(NX),JN(NX,NX),ANW(NX,NX), &
          STAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(*,'(A)')'RESULT p1t1n_alloc=1'
            RETURN
          END IF
          XST=X0V
          RDX=0.0_LDP
          OPEN(98,FILE='p1t1b_newton_'//TRIM(VNAME(IV))//'.txt', &
          STATUS='REPLACE')
          DO NWIT=0,4
            CALL MAP_F(XST,XNEW,IRFI(IV))
            NEVAL=NEVAL+1
            RN=0.0_LDP; RD=0.0_LDP
            DO J=1,NX
              RN=RN+(XNEW(J)-XST(J))**2
              RD=RD+XST(J)**2
            END DO
            RN=SQRT(RN); RD=SQRT(RD)
            WRITE(98,'(I4,3ES16.8)')NWIT,RN,RDX, &
      RN/MAX(RD,1.0E-30_LDP)
            WRITE(*,'(A,I2,A,I2,A,3ES16.8)')'P1B_NWT IV',IV, &
      ' IT',NWIT,' RN DX REL',RN,RDX, &
      RN/MAX(RD,1.0E-30_LDP)
            IF(NWIT .EQ. 4)EXIT
            IF(RN .EQ. 0.0_LDP)EXIT
            IF(NWIT .EQ. 0)THEN
              JN=J2
            ELSE
!
! FD Jacobian at the current iterate, same construction and
! amplitude as J2 above.
!
              AMPF=1.0E-5_LDP
              DO J=1,NX
                EPSV=AMPF*MAX(ABS(XST(J)),XFLOOR)
                XPR=XST
                XPR(J)=XPR(J)+EPSV
                CALL MAP_F(XPR,XFP,IRFI(IV))
                NEVAL=NEVAL+1
                IF(XST(J)-EPSV .GE. 0.0_LDP)THEN
                  XPR=XST
                  XPR(J)=XPR(J)-EPSV
                  CALL MAP_F(XPR,XO,IRFI(IV))
                  NEVAL=NEVAL+1
                  DO JP=1,NX
                    JN(JP,J)=(XFP(JP)-XO(JP))/(2.0_LDP*EPSV)
                  END DO
                ELSE
                  DO JP=1,NX
                    JN(JP,J)=(XFP(JP)-XNEW(JP))/EPSV
                  END DO
                END IF
              END DO
            END IF
!
! Solve (I - J) * DXP = F(XST) - XST; DXP is overwritten with the
! step. For an affine map XST+DXP is the exact fixed point.
!
            ANW=-JN
            DO J=1,NX
              ANW(J,J)=ANW(J,J)+1.0_LDP
            END DO
            DXP=XNEW-XST
            CALL DGESV(NX,1,ANW,NX,IPI8,DXP,NX,INFO8)
            WRITE(*,'(A,I2,A,I2,A,I4)')'P1B_NWT IV',IV, &
      ' IT',NWIT,' DGESV_INFO',INFO8
            IF(INFO8 .NE. 0)EXIT
            RDX=0.0_LDP
            DO J=1,NX
              RDX=RDX+DXP(J)**2
            END DO
            RDX=SQRT(RDX)
            XST=XST+DXP
          END DO
          CLOSE(98)
          OPEN(99,FILE='p1t1b_Xstar_'//TRIM(VNAME(IV))//'.txt', &
          STATUS='REPLACE')
          DO J=1,NX
            WRITE(99,'(I8,ES25.16)')J,XST(J)
          END DO
          CLOSE(99)
          WRITE(*,'(A,I2,A,I8)')'P1B_NWT IV',IV,' NEVAL',NEVAL
          DEALLOCATE(XST,DXP,IPI8,JN,ANW)
        END IF
      END DO
      INEXT_RELAX=0.5_LDP
      WRITE(*,'(A,I8)')'P1B_NEVAL ',NEVAL
!
! Direct nonlinear iteration on the full production map. Three
! perturbation types around X0: dominant-eigenvector aligned, broad
! (all-ones), and the largest stored outward-channel element. delta
! is 1E-3*||X0||_2. Per step we record ||X_{n+1}-X0||,
! ||X_{n+1}-X_n|| and the step-to-step contraction ratio; the first
! ratio of each type is meaningless (RNPRV=0) and is ignored on the
! python side. p1t1n skips this block: its product is the per-variant
! Newton residual history, not the relaxed-iteration probes.
!
      IF(CARG .NE. 'p1t1n')THEN
      TNRM=0.0_LDP
      DO J=1,NX
        TNRM=TNRM+X0V(J)**2
      END DO
      TNRM=SQRT(TNRM)
!
! Localised perturbation: the largest stored outward-channel element,
! i.e. the strongest compressing cell's lagged field. P1-T1 showed
! the outward channel carries the near-unity gain.
!
      IMX=1
      EMOD=-1.0_LDP
      DO SL=1,NSLOT-1
        DO I=1,NW
          JP=POSP(WKCH(I),SL)
          IF(JP .GT. 0 .AND. ABS(X0V(JP)) .GT. EMOD)THEN
            EMOD=ABS(X0V(JP)); IMX=JP
          END IF
        END DO
      END DO
!
      OPEN(97,FILE='p1t1b_nl.txt',STATUS='REPLACE')
      DO IRF=0,1
      DO IV=1,3
        IF(CARG .EQ. 'p1t1bf' .AND. IV .NE. 2) CYCLE
!
! Build the perturbation vector for this type.
!
        VPT=0.0_LDP
        IF(IV .EQ. 1)THEN
          VPT=VDOM
        ELSE IF(IV .EQ. 2)THEN
          VPT=1.0_LDP
        ELSE
          VPT(IMX)=1.0_LDP
        END IF
        RN=0.0_LDP
        DO J=1,NX
          RN=RN+VPT(J)**2
        END DO
        RN=SQRT(RN)
        IF(RN .GT. 0.0_LDP)VPT=VPT/RN
        XNEW=X0V+VPT*(1.0E-3_LDP*TNRM)
        RNPRV=0.0_LDP
!
        DO J=1,25
          CALL MAP_F(XNEW,XO,IRF)
          NEVAL=NEVAL+1
          RD=0.0_LDP; RN=0.0_LDP
          DO JP=1,NX
            RD=RD+(XO(JP)-X0V(JP))**2
            RN=RN+(XO(JP)-XNEW(JP))**2
          END DO
          RD=SQRT(RD); RN=SQRT(RN)
          WRITE(97,'(I2,I2,I4,3ES16.8)')IRF,IV,J,RD,RN, &
      RN/MAX(RNPRV,1.0E-30_LDP)
          RNPRV=RN
          XNEW=XO
        END DO
      END DO
      END DO
      CLOSE(97)
      END IF
!
      DEALLOCATE(XO,XNEW,X0V,FX0,XPR,XFP,STO_P,STO_M,J1,J2,JW, &
      RESV,VDOM,VPT,WKCH,POSP,POSM,FDIF,SPSP,SMSP,SMPP,SMMP)
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
!
! ------------------------------------------------------------------
! MAP_F: evaluate the complete ordered radiation map. Restores the
! slaved bookkeeping arrays to their post-cold-sweep snapshot,
! installs XV into the lagged store over the inactive-element
! baseline, runs the production frequency sweep with the normal
! in-sweep relaxed store update, optionally runs the production
! post-sweep refresh, then packs the resulting active store
! elements into XO. Loop variables are local so the host's DO-loop
! counters are not clobbered.
! IRF=1 includes FG_KNOP_REFRESH_V1 (full production map);
! IRF=0 stops after the sweep (bare lagged map).
! ------------------------------------------------------------------
      CONTAINS
!
      SUBROUTINE MAP_F(XV,XO,IRF)
      REAL(KIND=LDP) XV(:),XO(:)
      INTEGER IRF
      INTEGER L2,I2,SL2,J2,JP2
      LOGICAL LI,LF
!
! I_*_STORE/I_*_PREV are slaved: PREV is zeroed at INIT and copied
! from STORE at each NEW_FREQ, STORE is rewritten at the end of every
! frequency call before the next NEW_FREQ reads it. Restoring them
! to the post-cold snapshot only guarantees identical starting
! bookkeeping; no previous-sweep value can reach the next sweep.
!
      I_P_STORE=SPSP
      I_M_STORE=SMSP
      I_P_PREV=SMPP
      I_M_PREV=SMMP
!
! Same deterministic inputs as the cold sweep: the QW arrays feed
! only diagnostic J-moment outputs, and the output dummies are
! re-initialised so no previous evaluation leaks into this one.
!
      JQW_A=1.0_LDP; HQW_A=1.0_LDP; KQW_A=1.0_LDP
      NQW_A=1.0_LDP; HMID_A=1.0_LDP; NMID_A=1.0_LDP
      JNU_A=0.0_LDP; FEDD_A=0.0_LDP; IPLUS_A=0.0_LDP
      INBC_A=0.0_LDP; HBC_A=0.0_LDP
      I_P_DEPTH_NEXT=STO_P
      I_M_DEPTH_NEXT=STO_M
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          JP2=POSP(WKCH(I2),SL2)
          IF(JP2 .GT. 0) &
      I_P_DEPTH_NEXT(WKCH(I2),1,SL2)=XV(JP2)
          JP2=POSM(WKCH(I2),SL2)
          IF(JP2 .GT. 0) &
      I_M_DEPTH_NEXT(WKCH(I2),1,SL2)=XV(JP2)
        END DO
      END DO
      DO L2=1,NCFA
        IF(L2 .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L2+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)*0.75_LDP)
        END IF
        LI=(L2 .EQ. 1)
        LF=.TRUE.
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L2),CHI_M(1,L2),ESC_M(1,L2), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L2),DLOGA, &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,LI,LF,NC_A,NRAYS,NDA,SGT,NCFA)
      END DO
      IF(IRF .EQ. 1)CALL FG_KNOP_REFRESH_V1()
!
! Pack the post-map store into XO (slot-major order, matching POSP).
!
      J2=0
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_P_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_M_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE MAP_F
!
! PACK_X: extract the current store's active elements into XO
! (no map evaluation; local loop variables only).
!
      SUBROUTINE PACK_X(XO)
      REAL(KIND=LDP) XO(:)
      INTEGER L2,I2,SL2,J2
      J2=0
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_P_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_M_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE PACK_X
!
      END SUBROUTINE DRV_P1T1B
!
! ------------------------------------------------------------------
! DRV_P1T1S: production-coefficient packed-map experiment (Phase 2.9).
! Identical machinery to DRV_P1T1B but the physical input deck is read
! from 'p1t1s_in.txt' (extracted caseK2 J_RESCUE_DUMP window plus the
! nearest OUT_FLUX ESEC/CHI_SCAT table) instead of built inline. The
! FD dense Jacobian used at NX=30 is infeasible here (NX ~ 1.5E3), so
! the Jacobian action is probed by map matvecs: spectral radius by
! power iteration, Newton direction by GMRES on (I-J). The donor dump
! (DBS/DBSI/JMAP lines) lets the python side build the exact packed
! graph and test acyclicity without an FD pass.
! Diagnostic only; no production path changes.
! ------------------------------------------------------------------
      SUBROUTINE DRV_P1T1S
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=70,NPA=2,NCFA=13)
      REAL(KIND=LDP) RA(NDA),VA(NDA),SGA(NDA),SGT(NDA),PA(NPA)
      REAL(KIND=LDP) ETA_M(NDA,NCFA),CHI_M(NDA,NCFA),ESC_M(NDA,NCFA)
      REAL(KIND=LDP) JNU_A(NDA),FEDD_A(NDA),IPLUS_A(NPA)
      REAL(KIND=LDP) JQW_A(NDA,NPA),HQW_A(NDA,NPA),KQW_A(NDA,NPA)
      REAL(KIND=LDP) NQW_A(NDA,NPA),HMID_A(NDA,NPA),NMID_A(NDA,NPA)
      REAL(KIND=LDP) VDOP_A(NDA),FRQV(NCFA),DLNV(NCFA),TARR(NDA)
      REAL(KIND=LDP) DBB_A,IC_A,VDOP_F,REXT_A,INBC_A,HBC_A,DLOGA
      CHARACTER*6 METH_A
      CHARACTER(LEN=8) SOPT_A
      CHARACTER(LEN=12) IBND_A
      CHARACTER(LEN=64) WLN
      LOGICAL THK_A,INIT_A,NEWF_A
      INTEGER NC_A,NITER,ITR,L,I,LSDBG,IOS,NRAYS
      INTEGER NCOMP,NSNP,NSLOT,NW,NX
      INTEGER, ALLOCATABLE :: POSP(:,:),POSM(:,:),WKCH(:)
      REAL(KIND=LDP), ALLOCATABLE :: XO(:),XNEW(:),X0V(:),FX0(:)
      REAL(KIND=LDP), ALLOCATABLE :: STO_P(:,:,:),STO_M(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: SPSP(:,:),SMSP(:,:),SMPP(:,:), &
      SMMP(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: VDOM(:),VPT(:),RESV(:)
      REAL(KIND=LDP), ALLOCATABLE :: XSEED(:,:),XFIX(:,:)
      REAL(KIND=LDP) RN,RD,TNRM,EPSV,RLAM,RPREV,RLO,RHI
      INTEGER KC,SL,J,IV,IA,JP,K,NWARM,ISEED,KGM,NGM,NWT
      INTEGER*8 NEVAL
      INTEGER IRFI(3)
      REAL(KIND=LDP) OMV(3)
      CHARACTER(LEN=8) VNAME(3)
!
! Deck reader: token/section layout produced by the extractor script.
!
      OPEN(89,FILE='p1t1s_in.txt',STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1s_deck_missing=1'
        RETURN
      END IF
      READ(89,'(A)')WLN
      READ(89,'(A)')WLN
      READ(89,*)I,L
      IF(I .NE. NDA .OR. L .NE. NCFA)THEN
        WRITE(*,'(A,2I6)')'RESULT p1t1s_dim_mismatch=',I,L
        RETURN
      END IF
      CALL RDSECT(89,'R',RA)
      CALL RDSECT(89,'V',VA)
      CALL RDSECT(89,'SIG',SGA)
      CALL RDSECT(89,'SIGT',SGT)
      CALL RDSECT(89,'ESEC',ESC_M(1:NDA,1))
      CALL RDSECT(89,'CHI_SCAT',TARR)
      CALL RDSECT(89,'FL',FRQV)
      DO L=1,NCFA
        CALL RDSECT(89,'CHI',CHI_M(1:NDA,L))
      END DO
      DO L=1,NCFA
        CALL RDSECT(89,'ETA',ETA_M(1:NDA,L))
      END DO
      CLOSE(89)
      ESC_M=SPREAD(ESC_M(1:NDA,1),2,NCFA)
      FRQV=FRQV*1.0E15_LDP
      DO L=1,NCFA-1
        DLNV(L)=LOG(FRQV(L)/FRQV(L+1))
      END DO
      DLNV(NCFA)=DLNV(NCFA-1)
!
      PA(1)=0.0_LDP
      PA(2)=1.2E4_LDP
      VDOP_A=0.01_LDP*MAXVAL(ABS(VA))
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
      LSDBG=1
      NRAYS=1
!
      FG_FREQ_MATCH=.TRUE.
      FG_ORD_SWEEP=.TRUE.
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
      KNOP_DIAG_NO_REGRID=.TRUE.
      FREEZE_I_PREV=.FALSE.
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
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
      KN_DBG_LS=LSDBG
      ALLOCATE(KN_DBG_CHIH(NDA,NCFA),KN_DBG_Q(NDA,NCFA), &
      KN_DBG_DCHR(NDA,NCFA),KN_DBG_DTAU(NDA,NCFA), &
      KN_DBG_SPM(NDA,NCFA),KN_DBG_SPP(NDA,NCFA), &
      KN_DBG_WX(NDA,NCFA),KN_DBG_DLNM(NDA,NCFA), &
      KN_DBG_IMPR(NDA,NCFA),KN_DBG_IPPR(NDA,NCFA), &
      KN_DBG_NM(NDA,NCFA),KN_DBG_KC1(NDA,NCFA),KN_DBG_KC2(NDA,NCFA), &
      KN_DBG_A0(NDA,NCFA),KN_DBG_A1(NDA,NCFA),KN_DBG_A2(NDA,NCFA), &
      KN_DBG_A3(NDA,NCFA),KN_DBG_A4(NDA,NCFA),KN_DBG_SR(NDA,NCFA), &
      KN_DBG_INEXTM(NDA,NCFA),KN_DBG_INEXTP(NDA,NCFA), &
      KN_DBG_NDON(NDA,NCFA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1salloc_fail=1'
        RETURN
      END IF
      KN_DBG_CHIH=-1.0_LDP; KN_DBG_Q=-1.0_LDP; KN_DBG_DCHR=-1.0_LDP
      KN_DBG_DTAU=-1.0_LDP; KN_DBG_SPM=-1.0_LDP; KN_DBG_SPP=-1.0_LDP
      KN_DBG_WX=-1.0_LDP; KN_DBG_DLNM=-1.0_LDP
      KN_DBG_IMPR=-1.0_LDP; KN_DBG_IPPR=-1.0_LDP
      KN_DBG_NM=-1; KN_DBG_KC1=-1; KN_DBG_KC2=-1
      KN_DBG_A0=-1.0_LDP; KN_DBG_A1=-1.0_LDP; KN_DBG_A2=-1.0_LDP
      KN_DBG_A3=-1.0_LDP; KN_DBG_A4=-1.0_LDP; KN_DBG_SR=-1.0_LDP
      KN_DBG_INEXTM=-1.0_LDP; KN_DBG_INEXTP=-1.0_LDP
      KN_DBG_NDON=0
!
      JQW_A=1.0_LDP; HQW_A=1.0_LDP; KQW_A=1.0_LDP
      NQW_A=1.0_LDP; HMID_A=1.0_LDP; NMID_A=1.0_LDP
      JNU_A=0.0_LDP; FEDD_A=0.0_LDP; IPLUS_A=0.0_LDP
      INBC_A=0.0_LDP; HBC_A=0.0_LDP
!
! Cold sweep: iteration 1 of the packed fixed-point sequence.
!
      DO L=1,NCFA
        IF(L .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)* &
      EXP(-DLNV(NCFA)))
        END IF
        INIT_A=(L .EQ. 1)
        DLOGA=DLNV(L)
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L),ESC_M(1,L), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA, &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
      END DO
      CALL FG_KNOP_REFRESH_V1()
!
      IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT))THEN
        WRITE(*,'(A)')'RESULT p1t1sstore_missing=1'
        RETURN
      END IF
      NCOMP=SIZE(I_P_DEPTH_NEXT,1)
      NSNP=SIZE(I_P_DEPTH_NEXT,2)
      NSLOT=SIZE(I_P_DEPTH_NEXT,3)
      IF(NSNP .NE. 1)THEN
        WRITE(*,'(A)')'RESULT p1t1snp_unsupported=1'
        RETURN
      END IF
      ALLOCATE(STO_P(NCOMP,NSNP,NSLOT),STO_M(NCOMP,NSNP,NSLOT), &
      POSP(NCOMP,NSLOT),POSM(NCOMP,NSLOT),WKCH(NCOMP),STAT=IOS)
      STO_P=I_P_DEPTH_NEXT
      STO_M=I_M_DEPTH_NEXT
      NW=0
      POSP=0; POSM=0
      DO K=1,NDA
        KC=COMP_DEPTH_IDX(K)
        IF(KC .GE. 1)THEN
          I=J_PNT(K,1)
          IF(I .GE. 1 .AND. I .LE. NI_RAY(1))THEN
            NW=NW+1
            WKCH(NW)=KC
          END IF
        END IF
      END DO
      NX=2*NW*(NSLOT-1)
      J=0
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          POSP(WKCH(I),SL)=J
        END DO
      END DO
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          POSM(WKCH(I),SL)=J
        END DO
      END DO
      ALLOCATE(XO(NX),XNEW(NX),X0V(NX),FX0(NX),VDOM(NX),VPT(NX), &
      RESV(300),XSEED(NX,3),XFIX(NX,3), &
      SPSP(SIZE(I_P_STORE,1),SIZE(I_P_STORE,2)), &
      SMSP(SIZE(I_M_STORE,1),SIZE(I_M_STORE,2)), &
      SMPP(SIZE(I_P_PREV,1),SIZE(I_P_PREV,2)), &
      SMMP(SIZE(I_M_PREV,1),SIZE(I_M_PREV,2)),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1sx_alloc_fail=1'
        RETURN
      END IF
      SPSP=I_P_STORE
      SMSP=I_M_STORE
      SMPP=I_P_PREV
      SMMP=I_M_PREV
      WRITE(*,'(A,5I8)')'P1S_DIMS ',NX,NW,NSLOT,NCOMP,NI_RAY(1)
!
      OPEN(91,FILE='p1t1s_idx.txt',STATUS='REPLACE')
      WRITE(91,'(A,4I8)')'P1SIDX ',NX,NW,NSLOT,NCOMP
      J=0
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I8,1X,A1,1X,2I6)')'ELEM',J,'P',WKCH(I),SL
        END DO
      END DO
      DO SL=1,NSLOT-1
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I8,1X,A1,1X,2I6)')'ELEM',J,'M',WKCH(I),SL
        END DO
      END DO
      CLOSE(91)
!
! Warm up the packed map so X0V is a realistic iterate.
!
      NWARM=30
      CALL PACK_X(X0V)
      NEVAL=0
      DO ITR=2,NWARM
        CALL MAP_F(X0V,XNEW,1)
        NEVAL=NEVAL+1
        RN=SQRT(SUM((XNEW-X0V)**2))
        RESV(ITR)=RN
        X0V=XNEW
      END DO
!
! Donor/topology dump: one record set per (ray point, freq) so the
! python side can build the exact packed-graph edge list
! (write (KC_w,L) <- read (KC1|KC2,L+NM)).
!
      OPEN(96,FILE='p1t1s_dbs.txt',STATUS='REPLACE')
      DO I=1,NI_RAY(LSDBG)
        DO L=1,NCFA
          WRITE(96,'(A,I3,A,I3,A,11ES16.8)')'DBS ',I,' ',L,' ', &
      KN_DBG_CHIH(I,L),KN_DBG_Q(I,L),KN_DBG_DCHR(I,L), &
      KN_DBG_DTAU(I,L),KN_DBG_SPM(I,L),KN_DBG_SPP(I,L), &
      KN_DBG_WX(I,L),KN_DBG_DLNM(I,L),KN_DBG_IMPR(I,L), &
      KN_DBG_IPPR(I,L),KN_DBG_SR(I,L)
          WRITE(96,'(A,I3,A,I3,A,3I6)')'DBSI ',I,' ',L,' ', &
      KN_DBG_NM(I,L),KN_DBG_KC1(I,L),KN_DBG_KC2(I,L)
          WRITE(96,'(A,I3,A,I3,A,2ES16.8)')'DBSN ',I,' ',L,' ', &
      KN_DBG_INEXTM(I,L),KN_DBG_INEXTP(I,L)
        END DO
      END DO
      DO K=1,NDA
        WRITE(96,'(A,I3,A,2I6)')'JMAP ',K,' ',J_PNT(K,LSDBG), &
      COMP_DEPTH_IDX(K)
      END DO
      CLOSE(96)
!
      IRFI(1)=1; IRFI(2)=0; IRFI(3)=0
      OMV(1)=0.5_LDP; OMV(2)=0.5_LDP; OMV(3)=1.0_LDP
      VNAME(1)='fullw05'; VNAME(2)='swp_w05'; VNAME(3)='swp_w10'
!
      DO IV=1,3
        INEXT_RELAX=OMV(IV)
        CALL MAP_F(X0V,FX0,IRFI(IV))
        NEVAL=NEVAL+1
        TNRM=SQRT(SUM(X0V**2))
!
! Power iteration on the packed Jacobian via map matvecs. For a
! nilpotent J the iterate collapses in <= nilpotency-degree steps;
! that collapse is itself the test, not a failure of the estimator.
!
        OPEN(97,FILE='p1t1s_pow_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        DO J=1,NX
          VPT(J)=SIN(0.37_LDP*DBLE(J))
        END DO
        VPT=VPT/SQRT(SUM(VPT**2))
        RPREV=1.0_LDP
        RLO=1.0E30_LDP; RHI=0.0_LDP
        DO K=1,60
          CALL MATVEC(X0V,FX0,VPT,VDOM,IRFI(IV))
          NEVAL=NEVAL+1
          RLAM=SQRT(SUM(VDOM**2))
          WRITE(97,'(I4,3ES20.12)')K,RLAM,RLAM/RPREV
          IF(K .GT. 45)THEN
            RLO=MIN(RLO,RLAM); RHI=MAX(RHI,RLAM)
          END IF
          RPREV=MAX(RLAM,1.0E-300_LDP)
          IF(RLAM .GT. 0.0_LDP)VPT=VDOM/RLAM
        END DO
        WRITE(97,'(A,3ES20.12)')'RHO_EST ',RLAM,RLO,RHI
        CLOSE(97)
!
! Nilpotence chain: F^d(X0) becomes constant when d exceeds the
! longest slot-directed dependency path. Residual collapse to the
! map noise floor is the structural certificate.
!
        IF(IRFI(IV) .EQ. 0)THEN
          XNEW=X0V
          DO K=1,NSLOT+2
            CALL MAP_F(XNEW,XO,0)
            NEVAL=NEVAL+1
            RN=SQRT(SUM((XO-XNEW)**2))
            WRITE(*,'(A,I2,A,I4,A,ES16.8)')'P1S_NILP iv',IV, &
      ' k',K,' res',RN
            XNEW=XO
          END DO
        END IF
!
! JFNK-Newton from three seeds. Inner solve: GMRES on (I-J) via
! matvecs; no explicit matrix, no damping, no shift. Linear-solve
! residual and true map residual are both recorded.
!
        OPEN(98,FILE='p1t1s_newton_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        XSEED(1:NX,1)=X0V
        XSEED(1:NX,2)=0.5_LDP*X0V
        DO J=1,NX
          XSEED(J,2)=XSEED(J,2)+1.0E-3_LDP*MAX(ABS(X0V(J)),1.0E-3_LDP)* &
      SIN(1.7_LDP*DBLE(J))
        END DO
        XSEED(1:NX,3)=0.0_LDP
        DO ISEED=1,3
          XNEW=XSEED(1:NX,ISEED)
          DO NWT=1,6
            CALL MAP_F(XNEW,XO,IRFI(IV))
            NEVAL=NEVAL+1
            RN=SQRT(SUM((XO-XNEW)**2))
            IF(RN .LT. 1.0E-13_LDP*MAX(1.0_LDP,TNRM))EXIT
            VPT=XO-XNEW
            CALL GMRESV(XNEW,XO,VPT,VDOM,IRFI(IV),NGM,RD,NEVAL)
            XNEW=XNEW+VDOM
            CALL MAP_F(XNEW,XO,IRFI(IV))
            NEVAL=NEVAL+1
            RD=SQRT(SUM((XO-XNEW)**2))
            WRITE(98,'(A,I2,A,I2,A,I5,A,2ES16.8)')'NWT seed',ISEED, &
      ' it',NWT,' kgm',NGM,' res',RN,RD
          END DO
          CALL MAP_F(XNEW,XO,IRFI(IV))
          NEVAL=NEVAL+1
          XFIX(1:NX,ISEED)=XNEW
          RN=SQRT(SUM((XO-XNEW)**2))
          WRITE(98,'(A,I2,A,2ES16.8)')'SEED ',ISEED, &
      ' final',RN,SQRT(SUM(XNEW**2))
        END DO
        DO ISEED=2,3
          WRITE(98,'(A,I2,2ES16.8)')'XFIX_DIFF seed',ISEED, &
      SQRT(SUM((XFIX(1:NX,ISEED)-XFIX(1:NX,1))**2)), &
      MAXVAL(ABS(XFIX(1:NX,ISEED)-XFIX(1:NX,1)))
        END DO
        CLOSE(98)
      END DO
      WRITE(*,'(A,I10)')'P1S_NEVAL=',NEVAL
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
!
      CONTAINS
!
      SUBROUTINE RDSECT(U,KEY,ARR)
      INTEGER U
      CHARACTER(LEN=*) KEY
      REAL(KIND=LDP) ARR(:)
      CHARACTER(LEN=64) HDR
      READ(U,*)HDR
      IF(TRIM(HDR) .NE. TRIM(KEY))THEN
        WRITE(*,'(2A)')'RDSECT key mismatch: ',TRIM(HDR)
        STOP 9
      END IF
      READ(U,*)ARR
      RETURN
      END SUBROUTINE RDSECT
!
      SUBROUTINE MAP_F(XV,XO,IRF)
      REAL(KIND=LDP) XV(:),XO(:)
      INTEGER IRF
      INTEGER L2,I2,SL2,J2,JP2
      LOGICAL LI,LF
      I_P_STORE=SPSP
      I_M_STORE=SMSP
      I_P_PREV=SMPP
      I_M_PREV=SMMP
      JQW_A=1.0_LDP; HQW_A=1.0_LDP; KQW_A=1.0_LDP
      NQW_A=1.0_LDP; HMID_A=1.0_LDP; NMID_A=1.0_LDP
      JNU_A=0.0_LDP; FEDD_A=0.0_LDP; IPLUS_A=0.0_LDP
      INBC_A=0.0_LDP; HBC_A=0.0_LDP
      I_P_DEPTH_NEXT=STO_P
      I_M_DEPTH_NEXT=STO_M
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          JP2=POSP(WKCH(I2),SL2)
          IF(JP2 .GT. 0) &
      I_P_DEPTH_NEXT(WKCH(I2),1,SL2)=XV(JP2)
          JP2=POSM(WKCH(I2),SL2)
          IF(JP2 .GT. 0) &
      I_M_DEPTH_NEXT(WKCH(I2),1,SL2)=XV(JP2)
        END DO
      END DO
      DO L2=1,NCFA
        IF(L2 .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L2+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)*EXP(-DLNV(NCFA)))
        END IF
        LI=(L2 .EQ. 1)
        LF=.TRUE.
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L2),CHI_M(1,L2),ESC_M(1,L2), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L2),DLNV(L2), &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,LI,LF,NC_A,NRAYS,NDA,SGT,NCFA)
      END DO
      IF(IRF .EQ. 1)CALL FG_KNOP_REFRESH_V1()
      J2=0
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_P_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_M_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE MAP_F
!
      SUBROUTINE PACK_X(XO)
      REAL(KIND=LDP) XO(:)
      INTEGER L2,I2,SL2,J2
      J2=0
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_P_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      DO SL2=1,NSLOT-1
        DO I2=1,NW
          J2=J2+1
          XO(J2)=I_M_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE PACK_X
!
! MATVEC: J*v via one forward-difference map probe. The amplitude is
! scaled to the seed-point magnitude so vanishing elements still see
! a resolvable perturbation.
!
      SUBROUTINE MATVEC(XB,FB,V,W,IRF)
      REAL(KIND=LDP) XB(:),FB(:),V(:),W(:)
      INTEGER IRF
      REAL(KIND=LDP) EPSM,VS
      REAL(KIND=LDP), ALLOCATABLE :: XP(:),FP(:)
      ALLOCATE(XP(NX),FP(NX))
      VS=SQRT(SUM(V**2))
      EPSM=1.0E-5_LDP*MAX(SQRT(SUM(XB**2)),1.0E-3_LDP)/MAX(VS,1.0E-30_LDP)
      XP=XB+EPSM*V
      CALL MAP_F(XP,FP,IRF)
      W=(FP-FB)/EPSM
      DEALLOCATE(XP,FP)
      RETURN
      END SUBROUTINE MATVEC
!
! GMRESV: solve (I-J)*dx = rhs by unrestarted-then-restarted GMRES
! using MATVEC for the J action. Returns the Krylov count and the
! achieved linear residual; no shift or damping is applied.
!
      SUBROUTINE GMRESV(XB,FB,RHS,DX,IRF,NGM,RRES,NEV)
      REAL(KIND=LDP) XB(:),FB(:),RHS(:),DX(:)
      INTEGER IRF,NGM
      REAL(KIND=LDP) RRES
      INTEGER*8 NEV
      INTEGER KM,IU,JU,NRES
      PARAMETER (KM=80)
      REAL(KIND=LDP), ALLOCATABLE :: VV(:,:),HH(:,:),CS(:),SN(:), &
      GG(:),YY(:),W2(:)
      REAL(KIND=LDP) BETA,HT,T1,T2,RTOL
      ALLOCATE(VV(NX,KM+1),HH(KM+1,KM),CS(KM),SN(KM),GG(KM+1), &
      YY(KM),W2(NX))
      RTOL=1.0E-9_LDP
      DX=0.0_LDP
      NGM=0
      DO NRES=1,3
        CALL MATVEC(XB,FB,DX,W2,IRF)
        NEV=NEV+1
        W2=RHS-(DX-W2)
        BETA=SQRT(SUM(W2**2))
        IF(BETA .LT. RTOL*MAX(SQRT(SUM(RHS**2)),1.0E-30_LDP))EXIT
        VV(1:NX,1)=W2/BETA
        GG=0.0_LDP; GG(1)=BETA
        DO IU=1,KM
          CALL MATVEC(XB,FB,VV(1:NX,IU),W2,IRF)
          NEV=NEV+1
          W2=VV(1:NX,IU)-W2
          DO JU=1,IU
            HH(JU,IU)=SUM(W2*VV(1:NX,JU))
            W2=W2-HH(JU,IU)*VV(1:NX,JU)
          END DO
          HH(IU+1,IU)=SQRT(SUM(W2**2))
          IF(HH(IU+1,IU) .GT. 0.0_LDP) &
      VV(1:NX,IU+1)=W2/HH(IU+1,IU)
          DO JU=1,IU-1
            T1=CS(JU)*HH(JU,IU)+SN(JU)*HH(JU+1,IU)
            T2=-SN(JU)*HH(JU,IU)+CS(JU)*HH(JU+1,IU)
            HH(JU,IU)=T1; HH(JU+1,IU)=T2
          END DO
          HT=SQRT(HH(IU,IU)**2+HH(IU+1,IU)**2)
          IF(HT .EQ. 0.0_LDP)THEN
            CS(IU)=1.0_LDP; SN(IU)=0.0_LDP
          ELSE
            CS(IU)=HH(IU,IU)/HT; SN(IU)=HH(IU+1,IU)/HT
          END IF
          HH(IU,IU)=HT
          HH(IU+1,IU)=0.0_LDP
          T1=CS(IU)*GG(IU)
          GG(IU+1)=-SN(IU)*GG(IU)
          GG(IU)=T1
          NGM=NGM+1
          IF(ABS(GG(IU+1)) .LT. &
      RTOL*MAX(BETA,1.0E-30_LDP))EXIT
        END DO
        DO JU=IU,1,-1
          T1=GG(JU)
          DO K=JU+1,IU
            T1=T1-HH(JU,K)*YY(K)
          END DO
          IF(HH(JU,JU) .NE. 0.0_LDP)THEN
            YY(JU)=T1/HH(JU,JU)
          ELSE
            YY(JU)=0.0_LDP
          END IF
        END DO
        DO JU=1,IU
          DX=DX+YY(JU)*VV(1:NX,JU)
        END DO
        IF(ABS(GG(IU+1)) .LT. &
      RTOL*MAX(BETA,1.0E-30_LDP))EXIT
      END DO
      CALL MATVEC(XB,FB,DX,W2,IRF)
      NEV=NEV+1
      W2=RHS-(DX-W2)
      RRES=SQRT(SUM(W2**2))
      DEALLOCATE(VV,HH,CS,SN,GG,YY,W2)
      RETURN
      END SUBROUTINE GMRESV
!
      END SUBROUTINE DRV_P1T1S
!
! ------------------------------------------------------------------
! DRV_P1T1R: refresh inner-operator experiment (Phase 2.10).
! The iterated state is the per-ray field FLD(NI,NCF), seeded from
! the lagged store by SEED_FLD and iterated in place by the ordered
! march. With KN_RAY_MAXPASS=1 a KNOP_RAY_SOLVE_V1 call is exactly
! one production pass; KR_DBG_FLDM/KR_DBG_FLDP export the resulting
! field, so the diagnostic map
!     Y = seeded FLD at (depth-point, slot)  ->  G(Y) = post-pass FLD
! reproduces the original in-place ordering with no code change.
! Free state per channel: NWD depth points x NCFA-1 slots (slot 1 is
! never seeded and slots map L=s+1), packed P-then-M like p1t1s.
! OMEGA_RAY is frozen per experiment through KN_RAY_OMEGA; the
! pass-end adaptation is outside the one-pass operator.
! Diagnostic only; no production path changes.
! ------------------------------------------------------------------
      SUBROUTINE DRV_P1T1R(MODE)
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      CHARACTER(LEN=*) MODE
      INTEGER NDA,NPA,NCFA
      PARAMETER (NDA=70,NPA=2,NCFA=13)
      REAL(KIND=LDP) RA(NDA),VA(NDA),SGA(NDA),SGT(NDA),PA(NPA)
      REAL(KIND=LDP) ETA_M(NDA,NCFA),CHI_M(NDA,NCFA),ESC_M(NDA,NCFA)
      REAL(KIND=LDP) JNU_A(NDA),FEDD_A(NDA),IPLUS_A(NPA)
      REAL(KIND=LDP) JQW_A(NDA,NPA),HQW_A(NDA,NPA),KQW_A(NDA,NPA)
      REAL(KIND=LDP) NQW_A(NDA,NPA),HMID_A(NDA,NPA),NMID_A(NDA,NPA)
      REAL(KIND=LDP) VDOP_A(NDA),FRQV(NCFA),DLNV(NCFA),TARR(NDA)
      REAL(KIND=LDP) DBB_A,IC_A,VDOP_F,REXT_A,INBC_A,HBC_A,DLOGA
      CHARACTER*6 METH_A
      CHARACTER(LEN=8) SOPT_A
      CHARACTER(LEN=12) IBND_A
      CHARACTER(LEN=64) WLN
      LOGICAL THK_A,INIT_A,NEWF_A
      INTEGER NC_A,L,I,LSDBG,IOS,NRAYS
      INTEGER NI,NW,NX
      INTEGER, ALLOCATABLE :: WKCH(:),DPNT(:)
      REAL(KIND=LDP), ALLOCATABLE :: Y0(:),YO(:),YB(:),YV(:),YW(:)
      REAL(KIND=LDP), ALLOCATABLE :: YS(:,:),YF(:,:),HV(:,:),HH(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: STO_P(:,:,:),STO_M(:,:,:)
      REAL(KIND=LDP), ALLOCATABLE :: RESV(:)
      REAL(KIND=LDP) RN,RD,RLAM,RPREV,RLO,RHI,TNRM,AFFN,AFFD
      INTEGER KC,SL,J,IV,IA,JP,K,ISEED,KGM,NPASSSV
      INTEGER NCAPSV,NNGSV,NSTSV,NCPSV,NNLSV
      INTEGER*8 NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NP28,NBB8,NTB8
      INTEGER*8 NCLS8,NNLS8
      INTEGER*8 NEVAL
      INTEGER NXITSV,NMXSV
      REAL(KIND=LDP) RESIDSV,RHO_MSV,RHO_PSV
      INTEGER OMI(4)
      REAL(KIND=LDP) OMV(4),OMH(5),EPJ,QRES,T1,T2,T3,T4,RR(5)
      CHARACTER(LEN=8) VNAME(4),VHNAME(5)
      INTEGER NFLOOR,NCAPN
      EXTERNAL KNOP_RAY_SOLVE_V1
!
      OPEN(78,FILE='p1t1r_deck.txt',STATUS='OLD',IOSTAT=IOS)
      IF(IOS .EQ. 0)THEN
        READ(78,'(A)')WLN
        CLOSE(78)
      ELSE
        WLN='p1t1s_in.txt'
      END IF
      OPEN(89,FILE=TRIM(WLN),STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1r_deck_missing=1'
        RETURN
      END IF
      READ(89,'(A)')WLN
      READ(89,'(A)')WLN
      READ(89,*)I,L
      IF(I .NE. NDA .OR. L .NE. NCFA)THEN
        WRITE(*,'(A,2I6)')'RESULT p1t1r_dim_mismatch=',I,L
        RETURN
      END IF
      CALL RSEC(89,'R',RA)
      CALL RSEC(89,'V',VA)
      CALL RSEC(89,'SIG',SGA)
      CALL RSEC(89,'SIGT',SGT)
      CALL RSEC(89,'ESEC',ESC_M(1:NDA,1))
      CALL RSEC(89,'CHI_SCAT',TARR)
      CALL RSEC(89,'FL',FRQV)
      DO L=1,NCFA
        CALL RSEC(89,'CHI',CHI_M(1:NDA,L))
      END DO
      DO L=1,NCFA
        CALL RSEC(89,'ETA',ETA_M(1:NDA,L))
      END DO
      CLOSE(89)
      ESC_M=SPREAD(ESC_M(1:NDA,1),2,NCFA)
      FRQV=FRQV*1.0E15_LDP
      DO L=1,NCFA-1
        DLNV(L)=LOG(FRQV(L)/FRQV(L+1))
      END DO
      DLNV(NCFA)=DLNV(NCFA-1)
!
      PA(1)=0.0_LDP
      PA(2)=1.2E4_LDP
      VDOP_A=0.01_LDP*MAXVAL(ABS(VA))
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
      LSDBG=1
      NRAYS=1
!
      FG_FREQ_MATCH=.TRUE.
      FG_ORD_SWEEP=.TRUE.
      DO_KNOP_REFRESH=.TRUE.
      FREEZE_I_NEXT=.FALSE.
      KNOP_DIAG_NO_REGRID=.TRUE.
      FREEZE_I_PREV=.FALSE.
      FG_ORD_AITKEN=.FALSE.
      FG_ORD_AA=.FALSE.
      FG_ORD_P2DAMP=.FALSE.
      FG_NI1_FULL=(NCASE .EQ. 97 .OR. NCASE .EQ. 98)
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
      KN_DBG_LS=LSDBG
!
      JQW_A=1.0_LDP; HQW_A=1.0_LDP; KQW_A=1.0_LDP
      NQW_A=1.0_LDP; HMID_A=1.0_LDP; NMID_A=1.0_LDP
      JNU_A=0.0_LDP; FEDD_A=0.0_LDP; IPLUS_A=0.0_LDP
      INBC_A=0.0_LDP; HBC_A=0.0_LDP
!
      DO L=1,NCFA
        IF(L .LT. NCFA)THEN
          CALL SET_FG_NU_NEXT(FRQV(L+1))
        ELSE
          CALL SET_FG_NU_NEXT(FRQV(NCFA)*EXP(-DLNV(NCFA)))
        END IF
        INIT_A=(L .EQ. 1)
        DLOGA=DLNV(L)
        CALL FG_J_CMF_NONMONO_V1(ETA_M(1,L),CHI_M(1,L),ESC_M(1,L), &
      VA,SGA,RA,PA, &
      JNU_A,FEDD_A,JQW_A,HQW_A,KQW_A,NQW_A,HMID_A,NMID_A, &
      INBC_A,HBC_A,IPLUS_A,FRQV(L),DLOGA, &
      IBND_A,DBB_A,IC_A,VDOP_A,VDOP_F,REXT_A, &
      METH_A,SOPT_A,THK_A,INIT_A,NEWF_A,NC_A,NRAYS,NDA,SGT,NCFA)
      END DO
      CALL FG_KNOP_REFRESH_V1()
!
      IF(.NOT. ALLOCATED(I_P_DEPTH_NEXT))THEN
        WRITE(*,'(A)')'RESULT p1t1r_store_missing=1'
        RETURN
      END IF
      NI=NI_RAY(1)
      ALLOCATE(KR_DBG_FLDM(NI,NCFA),KR_DBG_FLDP(NI,NCFA),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1r_dbg_alloc=1'
        RETURN
      END IF
      ALLOCATE(STO_P(SIZE(I_P_DEPTH_NEXT,1),1, &
      SIZE(I_P_DEPTH_NEXT,3)), &
      STO_M(SIZE(I_M_DEPTH_NEXT,1),1,SIZE(I_M_DEPTH_NEXT,3)), &
      WKCH(NDA),DPNT(NDA),STAT=IOS)
      STO_P=I_P_DEPTH_NEXT
      STO_M=I_M_DEPTH_NEXT
!
! Writable set: compressing depths whose J_PNT lands on a valid ray
! point - the SEED_FLD source cells. Y packs slot 2..NCF of the
! field at those points, P channel first then M.
!
      NW=0
      DO K=1,NDA
        KC=COMP_DEPTH_IDX(K)
        IF(KC .GE. 1)THEN
          I=J_PNT(K,1)
          IF(I .GE. 1 .AND. I .LE. NI)THEN
            NW=NW+1
            WKCH(NW)=KC
            DPNT(NW)=I
          END IF
        END IF
      END DO
      NX=2*NW*(NCFA-1)
      ALLOCATE(Y0(NX),YO(NX),YB(NX),YV(NX),YW(NX),YS(NX,3),YF(NX,3), &
      HV(NX,42),HH(42,41),RESV(300),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT p1t1r_alloc_fail=1'
        RETURN
      END IF
      WRITE(*,'(A,4I8)')'P1R_DIMS ',NX,NW,NCFA,NI
!
      OPEN(91,FILE='p1t1r_idx.txt',STATUS='REPLACE')
      WRITE(91,'(A,3I8)')'P1RIDX ',NX,NW,NCFA
      J=0
      DO SL=2,NCFA
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I6,1X,A1,1X,3I6)')'ELEM',J,'P',WKCH(I),DPNT(I),SL
        END DO
      END DO
      DO SL=2,NCFA
        DO I=1,NW
          J=J+1
          WRITE(91,'(A,I6,1X,A1,1X,3I6)')'ELEM',J,'M',WKCH(I),DPNT(I),SL
        END DO
      END DO
      CLOSE(91)
!
! Baseline seed: the post-refresh committed store -> seed field.
!
      CALL YSEED(Y0)
      CALL GEVAL(Y0,YB,1.0_LDP)
      NEVAL=0
      TNRM=SQRT(SUM(Y0**2))
      WRITE(*,'(A,2ES16.8)')'P1R_Y0 ',TNRM,SQRT(SUM((YB-Y0)**2))
!
      OMV(1)=1.0_LDP; OMV(2)=0.5_LDP; OMV(3)=0.507_LDP; OMV(4)=0.507_LDP
      VNAME(1)='om100'; VNAME(2)='om050'; VNAME(3)='om507'
      VNAME(4)='om507r'
      OMH(1)=1.0_LDP; OMH(2)=0.75_LDP; OMH(3)=0.507_LDP
      OMH(4)=0.5_LDP; OMH(5)=0.25_LDP
      VHNAME(1)='om100'; VHNAME(2)='om075'; VHNAME(3)='om507'
      VHNAME(4)='om050'; VHNAME(5)='om025'
!
! jac mode: dense columnwise finite-difference Jacobian at om100
! and om050 for the full spectrum / SVD / non-normality analysis.
!
      IF(MODE .EQ. 'jac')THEN
        DO IV=1,2
          CALL GEVAL(Y0,YB,OMV(IV))
          NEVAL=NEVAL+1
          OPEN(96,FILE='p1t1r_jac_'//TRIM(VNAME(IV))//'.unf', &
      STATUS='REPLACE',FORM='UNFORMATTED')
          WRITE(96)NX
          DO J=1,NX
            EPJ=1.0E-5_LDP*MAX(ABS(Y0(J)),1.0E-3_LDP)
            YV=Y0
            YV(J)=YV(J)+EPJ
            CALL GEVAL(YV,YO,OMV(IV))
            NEVAL=NEVAL+1
            YW=(YO-YB)/EPJ
            WRITE(96)YW
          END DO
          CLOSE(96)
        END DO
        WRITE(*,'(A,I10)')'P1R_JAC_NEVAL=',NEVAL
        GOTO 995
      END IF
!
! hist mode: iterate G at fixed omega from a cold start and record
! the step-residual ratio sequence. This is the production pass
! history with the pass-end OMEGA_RAY adaptation disabled, which is
! exactly the fixed-omega experiment of Part 7.
!
      IF(MODE .EQ. 'hist')THEN
        DO IV=1,5
          OPEN(97,FILE='p1t1r_hist_'//TRIM(VHNAME(IV))//'.txt', &
      STATUS='REPLACE')
          YV=0.0_LDP
          RPREV=-1.0_LDP
          DO K=1,50
            CALL GEVAL(YV,YO,OMH(IV))
            NEVAL=NEVAL+1
            RN=SQRT(SUM((YO-YV)**2))/MAX(SQRT(SUM(YO**2)), &
      1.0E-30_LDP)
            WRITE(97,'(I4,3ES20.12)')K,RN, &
      MERGE(RN/RPREV,0.0_LDP,RPREV .GT. 0.0_LDP), &
      SQRT(SUM(YO**2))
            RPREV=RN
            YV=YO
          END DO
          CLOSE(97)
        END DO
        WRITE(*,'(A,I10)')'P1R_HIST_NEVAL=',NEVAL
        GOTO 995
      END IF
!
! testd mode: clean seeded run. The seed vector comes from a file
! (or 'NORMAL' for the un-seeded case); every invocation is a fresh
! process so no module state carries over. KN_Q_RESID is armed so
! the production verify residual lands in KN_Q_LS(1).
!
      IF(MODE .EQ. 'testd')THEN
        ALLOCATE(KN_Q_LS(1),KN_QS_STAT(1),KN_QA_Q(1,2), &
      KN_QA_J(1,2),KN_QA_L(1,2),KN_QA_VOLD(1,2),KN_QA_VN(1,2), &
      KN_QA_VM(1,2),STAT=IOS)
        KN_Q_RESID=.TRUE.
        OPEN(77,FILE='p1t1d_seed.txt',STATUS='OLD',IOSTAT=IOS)
        IF(IOS .EQ. 0)THEN
          READ(77,'(A)')WLN
          IF(WLN(1:1) .EQ. 'N' .OR. WLN(1:1) .EQ. 'n')THEN
            YV=Y0
          ELSE
            REWIND(77)
            READ(77,*)YV
          END IF
          CLOSE(77)
        ELSE
          YV=Y0
        END IF
        CALL GEVAL(YV,YO,1.0_LDP)
        NEVAL=NEVAL+1
        WRITE(*,'(A,3ES16.8)')'P1TD_GRES ', &
      SQRT(SUM((YO-YV)**2)), &
      SQRT(SUM((YO-YV)**2))/MAX(SQRT(SUM(YO**2)),1.0E-30_LDP), &
      MAXVAL(ABS(YO-YV))
        CALL YINJECT(YV)
        KN_RAY_OMEGA=1.0_LDP
        CALL KNOP_RAY_SOLVE_V1(1,RESIDSV,NSTSV,NPASSSV,NCAPSV,NNGSV, &
      NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NMXSV,NXITSV,NP28,NBB8,NTB8, &
      RHO_MSV,RHO_PSV,NCLS8,NNLS8)
        QRES=-1.0_LDP
        IF(ALLOCATED(KN_Q_LS))QRES=KN_Q_LS(1)
        WRITE(*,'(A,ES16.8,A,I4,A,ES16.8)')'P1TD_RESID=',RESIDSV, &
      ' npass=',NPASSSV,' qres=',QRES
        GOTO 995
      END IF
!
! ctrl mode: structure-aware controller strategies, implemented as
! an external per-pass omega policy over the same one-pass G. Each
! pass calls the production solve exactly once, so the in-place
! update order and branch logic are unchanged; only the omega
! sequence differs. Strategies:
!   1 = fixed omega=1 (nilpotent reference)
!   2 = transient window W=13 then production-style damping
!   3 = median-of-5 ratio controller
!   4 = structure-aware: omega=1 while the iterate contracts or
!       stays inside the transient window, else fall back to the
!       median controller
!
      IF(MODE .EQ. 'ctrl')THEN
        DO IV=1,4
          OPEN(97,FILE='p1t1r_ctrl_s'//CHAR(ICHAR('0')+IV)//'.txt', &
      STATUS='REPLACE')
          YV=0.0_LDP
          T1=1.0_LDP
          RPREV=-1.0_LDP
          RR=0.0_LDP
          DO K=1,40
            CALL GEVAL(YV,YO,T1)
            NEVAL=NEVAL+1
            RN=SQRT(SUM((YO-YV)**2))/MAX(SQRT(SUM(YO**2)), &
      1.0E-30_LDP)
            NFLOOR=COUNT(YO .LE. 1.0E-29_LDP)
            T2=MERGE(RN/RPREV,0.0_LDP,RPREV .GT. 0.0_LDP)
            WRITE(97,'(I4,4ES20.12,I6)')K,RN,T2,T1, &
      SQRT(SUM(YO**2)),NFLOOR
            RR(1:4)=RR(2:5)
            RR(5)=T2
            RPREV=RN
            YV=YO
            IF(RN .LT. 1.0E-13_LDP)EXIT
            IF(IV .EQ. 2)THEN
              IF(K .GT. 13 .AND. T2 .GT. 1.05_LDP) &
      T1=MAX(0.45_LDP,T1*0.7_LDP)
            ELSE IF(IV .EQ. 3)THEN
              IF(K .GE. 5)THEN
                CALL MED5(RR,T4)
                IF(T4 .GT. 1.05_LDP)T1=MAX(0.45_LDP,T1*0.7_LDP)
                IF(T4 .LT. 0.5_LDP)T1=MIN(1.0_LDP,T1*1.15_LDP)
              END IF
            ELSE IF(IV .EQ. 4)THEN
              IF(K .GT. 13)THEN
                CALL MED5(RR,T4)
                IF(T4 .GT. 1.05_LDP)T1=MAX(0.45_LDP,T1*0.7_LDP)
              END IF
            END IF
          END DO
          CLOSE(97)
        END DO
        WRITE(*,'(A,I10)')'P1R_CTRL_NEVAL=',NEVAL
        GOTO 995
      END IF
!
! cert mode: production solve with KN_Q_RESID armed, then classify
! the outcome under the Phase 2.13 convergence contract using the
! returned exit class and residual quantities only. Optional
! p1t1x_opts.txt supplies a KN_RAY_BULK_TOL override (strategy D).
!
      IF(MODE .EQ. 'cert')THEN
        OPEN(76,FILE='p1t1x_opts.txt',STATUS='OLD',IOSTAT=IOS)
        IF(IOS .EQ. 0)THEN
          READ(76,*,IOSTAT=IOS)T3
          IF(IOS .EQ. 0)KN_RAY_BULK_TOL=T3
          READ(76,*,IOSTAT=IOS)T3
          IF(IOS .EQ. 0)KN_RAY_MAXPASS=NINT(T3)
          CLOSE(76)
        END IF
        ALLOCATE(KN_Q_LS(1),KN_QS_STAT(1),KN_QA_Q(1,2), &
      KN_QA_J(1,2),KN_QA_L(1,2),KN_QA_VOLD(1,2),KN_QA_VN(1,2), &
      KN_QA_VM(1,2),STAT=IOS)
        KN_Q_RESID=.TRUE.
        OPEN(77,FILE='p1t1d_seed.txt',STATUS='OLD',IOSTAT=IOS)
        IF(IOS .EQ. 0)THEN
          READ(77,'(A)')WLN
          IF(WLN(1:1) .EQ. 'N' .OR. WLN(1:1) .EQ. 'n')THEN
            YV=Y0
          ELSE
            REWIND(77)
            READ(77,*)YV
          END IF
          CLOSE(77)
        ELSE
          YV=Y0
        END IF
        CALL YINJECT(YV)
        KN_RAY_OMEGA=1.0_LDP
        CALL KNOP_RAY_SOLVE_V1(1,RESIDSV,NSTSV,NPASSSV,NCAPSV,NNGSV, &
      NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NMXSV,NXITSV,NP28,NBB8,NTB8, &
      RHO_MSV,RHO_PSV,NCLS8,NNLS8)
        QRES=-1.0_LDP
        IF(ALLOCATED(KN_Q_LS))QRES=KN_Q_LS(1)
        WRITE(*,'(A,I3,A,ES16.8,A,I4,A,I10,A,I10,A,ES16.8)') &
      'P1X_RAW ixit=',NXITSV,' resid=',RESIDSV,' npass=',NPASSSV, &
      ' nbad=',NBAD8,' ntot=',NTOT8,' qres=',QRES
        IF(NDIV8 .GT. 0 .OR. RESIDSV .NE. RESIDSV)THEN
          WLN='INVALID_STATE'
        ELSE IF(NXITSV .EQ. 5)THEN
          WLN='DIVERGENCE'
        ELSE IF(NXITSV .EQ. 4)THEN
          WLN='ITERATION_LIMIT'
        ELSE IF(NXITSV .EQ. 3 .OR. NXITSV .EQ. 6)THEN
          WLN='STAGNATION'
        ELSE IF(NXITSV .EQ. 1 .OR. (NXITSV .EQ. 2 .AND. &
      RESIDSV .LE. KN_RAY_TOL .AND. QRES .LE. KN_RAY_TOL))THEN
          WLN='FULL_CONVERGENCE'
        ELSE IF(NXITSV .EQ. 2)THEN
          WLN='BULK_ONLY'
        ELSE
          WLN='INVALID_STATE'
        END IF
        WRITE(*,'(A,A12,A,ES16.8)')'P1X_STATUS ',WLN, &
      ' btol=',KN_RAY_BULK_TOL
        GOTO 995
      END IF
!
      DO IV=1,3
!
! Power iteration on J_G at the warm state.
!
        OPEN(97,FILE='p1t1r_pow_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        DO J=1,NX
          YV(J)=SIN(0.37_LDP*DBLE(J))
        END DO
        YV=YV/SQRT(SUM(YV**2))
        CALL GEVAL(Y0,YB,OMV(IV))
        NEVAL=NEVAL+1
        RPREV=1.0_LDP
        RLO=1.0E30_LDP; RHI=0.0_LDP
        DO K=1,60
          CALL GMATVEC(Y0,YB,YV,YW,OMV(IV))
          NEVAL=NEVAL+1
          RLAM=SQRT(SUM(YW**2))
          WRITE(97,'(I4,3ES20.12)')K,RLAM,RLAM/RPREV
          IF(K .GT. 45)THEN
            RLO=MIN(RLO,RLAM); RHI=MAX(RHI,RLAM)
          END IF
          RPREV=MAX(RLAM,1.0E-300_LDP)
          IF(RLAM .GT. 0.0_LDP)YV=YW/RLAM
        END DO
        WRITE(97,'(A,3ES20.12)')'RHO_EST ',RLAM,RLO,RHI
        CLOSE(97)
!
! Arnoldi Ritz values: 40-step Krylov basis on J_G; the small H
! eigenproblem goes to DGEEV externally by the python side; here we
! only dump the Arnoldi vectors' Rayleigh quotients.
!
        OPEN(97,FILE='p1t1r_arnoldi_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        DO J=1,NX
          HV(J,1)=SIN(0.71_LDP*DBLE(J))+0.3_LDP*COS(0.19_LDP*DBLE(J))
        END DO
        HV(1:NX,1)=HV(1:NX,1)/SQRT(SUM(HV(1:NX,1)**2))
        DO K=1,40
          CALL GMATVEC(Y0,YB,HV(1:NX,K),YW,OMV(IV))
          NEVAL=NEVAL+1
          DO IA=1,K
            HH(IA,K)=SUM(YW*HV(1:NX,IA))
            YW=YW-HH(IA,K)*HV(1:NX,IA)
          END DO
          HH(K+1,K)=SQRT(SUM(YW**2))
          IF(HH(K+1,K) .GT. 0.0_LDP .AND. K .LT. 41) &
      HV(1:NX,K+1)=YW/HH(K+1,K)
        END DO
        DO K=1,41
          WRITE(97,*) (HH(K,IA),IA=1,MIN(40,K))
        END DO
        CLOSE(97)
!
! JFNK-Newton on G(Y)=Y from three seeds.
!
        OPEN(98,FILE='p1t1r_newton_'//TRIM(VNAME(IV))//'.txt', &
        STATUS='REPLACE')
        YS(1:NX,1)=Y0
        YS(1:NX,2)=0.5_LDP*Y0
        DO J=1,NX
          YS(J,2)=YS(J,2)+1.0E-3_LDP*MAX(ABS(Y0(J)),1.0E-3_LDP)* &
      SIN(1.7_LDP*DBLE(J))
        END DO
        YS(1:NX,3)=0.0_LDP
        DO ISEED=1,3
          YV=YS(1:NX,ISEED)
          DO K=1,6
            CALL GEVAL(YV,YO,OMV(IV))
            NEVAL=NEVAL+1
            RN=SQRT(SUM((YO-YV)**2))
            IF(RN .LT. 1.0E-13_LDP*MAX(1.0_LDP,TNRM))EXIT
            YW=YO-YV
            CALL GMRESG(YV,YO,YW,YB,OMV(IV),KGM,RD,NEVAL)
            YV=YV+YB
            CALL GEVAL(YV,YO,OMV(IV))
            NEVAL=NEVAL+1
            RD=SQRT(SUM((YO-YV)**2))
            WRITE(98,'(A,I2,A,I2,A,I5,A,2ES16.8)')'NWT s',ISEED, &
      ' i',K,' kg',KGM,' r',RN,RD
          END DO
          YF(1:NX,ISEED)=YV
          WRITE(98,'(A,I2,A,2ES16.8)')'SEED ',ISEED,' fin',RD, &
      SQRT(SUM(YV**2))
        END DO
        DO ISEED=2,3
          WRITE(98,'(A,I2,2ES16.8)')'YF_DIFF s',ISEED, &
      SQRT(SUM((YF(1:NX,ISEED)-YF(1:NX,1))**2)), &
      MAXVAL(ABS(YF(1:NX,ISEED)-YF(1:NX,1)))
        END DO
        CLOSE(98)
      END DO
!
! Affinity test at omega=1: G(Y1+Y2)-G(0) vs G(Y1)-G(0)+G(Y2)-G(0).
!
      OPEN(99,FILE='p1t1r_affine.txt',STATUS='REPLACE')
      CALL GEVAL(YS(1:NX,3),YB,1.0_LDP)
      NEVAL=NEVAL+1
      YS(1:NX,3)=YB
      CALL GEVAL(Y0,YO,1.0_LDP)
      NEVAL=NEVAL+1
      YS(1:NX,1)=YO-YB
      DO J=1,NX
        YV(J)=0.3_LDP*Y0(J)*COS(0.9_LDP*DBLE(J))+ &
      1.0E-4_LDP*SIN(2.1_LDP*DBLE(J))
      END DO
      CALL GEVAL(YV,YO,1.0_LDP)
      NEVAL=NEVAL+1
      YS(1:NX,2)=YO-YB
      YW=Y0+YV
      CALL GEVAL(YW,YO,1.0_LDP)
      NEVAL=NEVAL+1
      YW=YO-YB
      AFFD=SQRT(SUM((YW-(YS(1:NX,1)+YS(1:NX,2)))**2))
      AFFN=SQRT(SUM(YS(1:NX,1)**2))+SQRT(SUM(YS(1:NX,2)**2))
      WRITE(99,'(A,3ES16.8)')'AFFINE_ERR ',AFFD,AFFN,AFFD/MAX(AFFN, &
      1.0E-30_LDP)
      CLOSE(99)
!
! Branch census: count G(Y0) outputs sitting on the floor or cap,
! and dump the omega=1 fixed points so testd runs can seed from
! them in fresh processes.
!
      NFLOOR=COUNT(YB .LE. 1.0E-29_LDP)
      NCAPN=COUNT(YB .GE. 0.5_LDP*NONMONO_INEXT_CAP)
      WRITE(*,'(A,2I8,A,I8)')'P1R_CENSUS ',NFLOOR,NCAPN,' of ',NX
      OPEN(88,FILE='p1t1r_ystar_om100.txt',STATUS='REPLACE')
      WRITE(88,'(ES24.16)')YF(1:NX,1)
      CLOSE(88)
!
! TEST D: seed the store with the omega=1 fixed point and run the
! original refresh path (full pass budget) - does KN_Q_RESID/RESID
! collapse relative to the un-seeded original iteration?
!
      IF(NW .GT. 0)THEN
        ALLOCATE(KN_Q_LS(1),KN_QS_STAT(1),KN_QA_Q(1,2), &
      KN_QA_J(1,2),KN_QA_L(1,2),KN_QA_VOLD(1,2),KN_QA_VN(1,2), &
      KN_QA_VM(1,2),STAT=IOS)
        KN_Q_RESID=.TRUE.
        YW=YF(1:NX,1)
        CALL YINJECT(YW)
        KN_RAY_OMEGA=1.0_LDP
        CALL KNOP_RAY_SOLVE_V1(1,RESIDSV,NSTSV,NPASSSV,NCAPSV,NNGSV, &
      NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NMXSV,NXITSV,NP28,NBB8,NTB8, &
      RHO_MSV,RHO_PSV,NCLS8,NNLS8)
        QRES=-1.0_LDP
        IF(ALLOCATED(KN_Q_LS))QRES=KN_Q_LS(1)
        WRITE(*,'(A,ES16.8,A,I4,A,ES16.8)')'P1R_SEEDED_RESID=', &
      RESIDSV,' npass=',NPASSSV,' qres=',QRES
        CALL YINJECT(Y0)
        CALL KNOP_RAY_SOLVE_V1(1,RESIDSV,NSTSV,NPASSSV,NCAPSV,NNGSV, &
      NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NMXSV,NXITSV,NP28,NBB8,NTB8, &
      RHO_MSV,RHO_PSV,NCLS8,NNLS8)
        QRES=-1.0_LDP
        IF(ALLOCATED(KN_Q_LS))QRES=KN_Q_LS(1)
        WRITE(*,'(A,ES16.8,A,I4,A,ES16.8)')'P1R_UNS_RESID=',RESIDSV, &
      ' npass=',NPASSSV,' qres=',QRES
      END IF
  995 CONTINUE
      WRITE(*,'(A,I10)')'P1R_NEVAL=',NEVAL
      WRITE(*,'(A)')'RESULT_OK'
      RETURN
!
      CONTAINS
!
      SUBROUTINE RSEC(U,KEY,ARR)
      INTEGER U
      CHARACTER(LEN=*) KEY
      REAL(KIND=LDP) ARR(:)
      CHARACTER(LEN=64) HDR
      READ(U,*)HDR
      IF(TRIM(HDR) .NE. TRIM(KEY))THEN
        WRITE(*,'(2A)')'RSEC key mismatch: ',TRIM(HDR)
        STOP 9
      END IF
      READ(U,*)ARR
      RETURN
      END SUBROUTINE RSEC
!
! YSEED: pack the seed field implied by the current store: SEED_FLD
! puts STORE(KC,LS,s) into FLD(J_PNT,s+1), so Y element (KC,s) is the
! store slot s for s=1..NCFA-1.
!
      SUBROUTINE YSEED(Y)
      REAL(KIND=LDP) Y(:)
      INTEGER I2,SL2,J2
      J2=0
      DO SL2=1,NCFA-1
        DO I2=1,NW
          J2=J2+1
          Y(J2)=I_P_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      DO SL2=1,NCFA-1
        DO I2=1,NW
          J2=J2+1
          Y(J2)=I_M_DEPTH_NEXT(WKCH(I2),1,SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE YSEED
!
! YINJECT: write Y into the store slots the seed reads.
!
      SUBROUTINE YINJECT(Y)
      REAL(KIND=LDP) Y(:)
      INTEGER I2,SL2,J2
      J2=0
      DO SL2=1,NCFA-1
        DO I2=1,NW
          J2=J2+1
          I_P_DEPTH_NEXT(WKCH(I2),1,SL2)=Y(J2)
        END DO
      END DO
      DO SL2=1,NCFA-1
        DO I2=1,NW
          J2=J2+1
          I_M_DEPTH_NEXT(WKCH(I2),1,SL2)=Y(J2)
        END DO
      END DO
      RETURN
      END SUBROUTINE YINJECT
!
! GEVAL: one production pass. Restores the post-refresh store
! snapshot, injects the seed vector, runs KNOP_RAY_SOLVE_V1 with a
! single-pass budget at the given OMEGA_RAY, then repacks the
! exported fields at the seeded ray points.
!
      SUBROUTINE GEVAL(YI,YO,OMG)
      REAL(KIND=LDP) YI(:),YO(:),OMG
      INTEGER I2,SL2,J2,MP_SV
      I_P_DEPTH_NEXT=STO_P
      I_M_DEPTH_NEXT=STO_M
      CALL YINJECT(YI)
      MP_SV=KN_RAY_MAXPASS
      KN_RAY_MAXPASS=1
      KN_RAY_OMEGA=OMG
      CALL KNOP_RAY_SOLVE_V1(1,RESIDSV,NSTSV,NPASSSV,NCAPSV,NNGSV, &
      NBAD8,NTOT8,NDIV8,NMAT8,NCMP8,NMXSV,NXITSV,NP28,NBB8,NTB8, &
      RHO_MSV,RHO_PSV,NCLS8,NNLS8)
      KN_RAY_MAXPASS=MP_SV
      J2=0
      DO SL2=2,NCFA
        DO I2=1,NW
          J2=J2+1
          YO(J2)=KR_DBG_FLDP(DPNT(I2),SL2)
        END DO
      END DO
      DO SL2=2,NCFA
        DO I2=1,NW
          J2=J2+1
          YO(J2)=KR_DBG_FLDM(DPNT(I2),SL2)
        END DO
      END DO
      RETURN
      END SUBROUTINE GEVAL
!
      SUBROUTINE MED5(R,MED)
      REAL(KIND=LDP) R(5),MED,S(5)
      INTEGER I1,J1
      S=R
      DO I1=1,4
        DO J1=I1+1,5
          IF(S(J1) .LT. S(I1))THEN
            MED=S(I1); S(I1)=S(J1); S(J1)=MED
          END IF
        END DO
      END DO
      MED=S(3)
      RETURN
      END SUBROUTINE MED5
!
      SUBROUTINE GMATVEC(YB,FB,V,W,OMG)
      REAL(KIND=LDP) YB(:),FB(:),V(:),W(:),OMG
      REAL(KIND=LDP) EPSM,VS
      REAL(KIND=LDP), ALLOCATABLE :: XP(:),FP(:)
      ALLOCATE(XP(NX),FP(NX))
      VS=SQRT(SUM(V**2))
      EPSM=1.0E-5_LDP*MAX(SQRT(SUM(YB**2)),1.0E-3_LDP)/MAX(VS,1.0E-30_LDP)
      XP=YB+EPSM*V
      CALL GEVAL(XP,FP,OMG)
      W=(FP-FB)/EPSM
      DEALLOCATE(XP,FP)
      RETURN
      END SUBROUTINE GMATVEC
!
      SUBROUTINE GMRESG(YB,FB,RHS,DX,OMG,NGM,RRES,NEV)
      REAL(KIND=LDP) YB(:),FB(:),RHS(:),DX(:),OMG
      INTEGER NGM
      REAL(KIND=LDP) RRES
      INTEGER*8 NEV
      INTEGER KM,IU,JU,NRES
      PARAMETER (KM=80)
      REAL(KIND=LDP), ALLOCATABLE :: VV(:,:),HHG(:,:),CS(:),SN(:), &
      GG(:),YY(:),W2(:)
      REAL(KIND=LDP) BETA,HT,T1,T2,RTOL
      ALLOCATE(VV(NX,KM+1),HHG(KM+1,KM),CS(KM),SN(KM),GG(KM+1), &
      YY(KM),W2(NX))
      RTOL=1.0E-9_LDP
      DX=0.0_LDP
      NGM=0
      DO NRES=1,3
        CALL GMATVEC(YB,FB,DX,W2,OMG)
        NEV=NEV+1
        W2=RHS-(DX-W2)
        BETA=SQRT(SUM(W2**2))
        IF(BETA .LT. RTOL*MAX(SQRT(SUM(RHS**2)),1.0E-30_LDP))EXIT
        VV(1:NX,1)=W2/BETA
        GG=0.0_LDP; GG(1)=BETA
        DO IU=1,KM
          CALL GMATVEC(YB,FB,VV(1:NX,IU),W2,OMG)
          NEV=NEV+1
          W2=VV(1:NX,IU)-W2
          DO JU=1,IU
            HHG(JU,IU)=SUM(W2*VV(1:NX,JU))
            W2=W2-HHG(JU,IU)*VV(1:NX,JU)
          END DO
          HHG(IU+1,IU)=SQRT(SUM(W2**2))
          IF(HHG(IU+1,IU) .GT. 0.0_LDP) &
      VV(1:NX,IU+1)=W2/HHG(IU+1,IU)
          DO JU=1,IU-1
            T1=CS(JU)*HHG(JU,IU)+SN(JU)*HHG(JU+1,IU)
            T2=-SN(JU)*HHG(JU,IU)+CS(JU)*HHG(JU+1,IU)
            HHG(JU,IU)=T1; HHG(JU+1,IU)=T2
          END DO
          HT=SQRT(HHG(IU,IU)**2+HHG(IU+1,IU)**2)
          IF(HT .EQ. 0.0_LDP)THEN
            CS(IU)=1.0_LDP; SN(IU)=0.0_LDP
          ELSE
            CS(IU)=HHG(IU,IU)/HT; SN(IU)=HHG(IU+1,IU)/HT
          END IF
          HHG(IU,IU)=HT
          HHG(IU+1,IU)=0.0_LDP
          T1=CS(IU)*GG(IU)
          GG(IU+1)=-SN(IU)*GG(IU)
          GG(IU)=T1
          NGM=NGM+1
          IF(ABS(GG(IU+1)) .LT. &
      RTOL*MAX(BETA,1.0E-30_LDP))EXIT
        END DO
        DO JU=IU,1,-1
          T1=GG(JU)
          DO K=JU+1,IU
            T1=T1-HHG(JU,K)*YY(K)
          END DO
          IF(HHG(JU,JU) .NE. 0.0_LDP)THEN
            YY(JU)=T1/HHG(JU,JU)
          ELSE
            YY(JU)=0.0_LDP
          END IF
        END DO
        DO JU=1,IU
          DX=DX+YY(JU)*VV(1:NX,JU)
        END DO
        IF(ABS(GG(IU+1)) .LT. &
      RTOL*MAX(BETA,1.0E-30_LDP))EXIT
      END DO
      CALL GMATVEC(YB,FB,DX,W2,OMG)
      NEV=NEV+1
      W2=RHS-(DX-W2)
      RRES=SQRT(SUM(W2**2))
      DEALLOCATE(VV,HHG,CS,SN,GG,YY,W2)
      RETURN
      END SUBROUTINE GMRESG
!
      END SUBROUTINE DRV_P1T1R
