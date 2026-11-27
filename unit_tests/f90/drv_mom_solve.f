! Driver for MOM_J_CMF_NONMONO_V1 (plane/mom_j_cmf_nonmono_v1.f).
! Manufactured model atmosphere: V profile integrated from a prescribed
! sigma_true(r) so compression bands of controlled depth/width land on
! the expanded grid.  Ray-moment store (MOD_RAY_MOM_STORE) is seeded
! with self-consistent Eddington factors, then the production solver is
! driven over NCF frequencies for NIT "great iterations" (sweeps).
! Dumps the internal coefficient arrays for independent verification.
!
! Usage: drv_mom_solve <mode> <ntype> <niter> <outbc> <inbm> <hchk>
!                      <match> <tag>
!   mode : homo | beta | dec | comp1 | comp2 | strong | jump | zero |
!          sigx | vconst | compw | recall (comp1 + NEW_FREQ=F recall
!          per frequency -- the comp_j_blank INACCURATE pattern; the
!         sweep census must not exceed GS_NTOT) | nurf (comp1 +
!         NUSOLV-off vs NUSOLV-on refresh A/B -- fixed-point
!         invariance of the implicit frequency-chain pass) | rfst
!         (compw + one plain ALT refresh on the thick scattering
!         atmosphere -- the compressing-branch stall fixture)
!   ntype: nonj | gonly
!   niter: number of sweeps over all NCF frequencies
!   outbc: 1 (Eddington closure) or 2 (negated interior row)
!   inbm : diff | zfx | pnt
!   hchk : none | av | max
!   match: 0/1 -> MOM_FREQ_MATCH
!   tag  : output file prefix (written to cwd)
!   Optional: arg9 = ESEC/CHI albedo (default 0.6), arg10 = CHI
!   multiplier (advection-dominated thin limit), arg11 = log opacity
!   modulation across frequency (blanketing analogue), arg12 =
!   per-sweep ray-store drift amplitude (moving-fixed-point
!   emulation of the production coefficient update).
!
      PROGRAM DRV_MOM_SOLVE
      USE SET_KIND_MODULE
      USE MOD_MOM_J_NONMONO_V1
      USE MOD_RAY_MOM_STORE
      USE CONTROL_VARIABLE_MOD, ONLY: CONV_DIAG_LVL
      IMPLICIT NONE
!
      INTEGER ND,NCF
#ifndef ND_IN
#define ND_IN 24
#endif
#ifndef NCF_IN
#define NCF_IN 8
#endif
      PARAMETER (ND=ND_IN,NCF=NCF_IN)
      INTEGER I,L,IT,IOS,NIT,NCASE,OUTBC
      REAL(KIND=LDP) R(ND),V(ND),SIG(ND),SIGT(ND)
      REAL(KIND=LDP) CHI(ND),ETA(ND),ESEC(ND),VDOP(ND),CHI_L(ND)
      REAL(KIND=LDP) JNU(ND),RSQH(ND)
      REAL(KIND=LDP) FREQ,DLN,DLNC,V0,RC,BETA,T1,T2,S0,DLNR
      REAL(KIND=LDP) ALB,CHISCL,CHIWAV,PERT,PFAC,BALB
      REAL(KIND=LDP) DLNA,V0A,RFBLK,RFMAXP,RFRLX,RFAA,RFDIR
      REAL(KIND=LDP) GF_ARG,RFGS,RFRLX2
      REAL(KIND=LDP) VDOP_FRAC,DBB,IC,IBSTAB,FJ,HF,GF
      REAL(KIND=LDP) FRQ(NCF),NU_NXT
      CHARACTER(LEN=32) MODE,NTYPE,INBM,HCHK,OUTTAG,CARG
      CHARACTER(LEN=6) NTYP6,METH,HCHK6
      CHARACTER(LEN=12) INBM6
      LOGICAL COHR,INIT,NEW_FREQ,DOMATCH
      INTEGER FIDX
      INTEGER IB1,IB2
      REAL(KIND=LDP), ALLOCATABLE :: CKJ(:,:),CKH(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: JRF0(:,:),HRF0(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: JDSC(:,:),JASC(:,:),JNUV(:,:)
      REAL(KIND=LDP) JV(ND),HSC(ND),V11MAX,V11REL
      INTEGER(KIND=2), ALLOCATABLE :: CKQ(:,:)
      REAL(KIND=LDP) DMAX
!
! Array dump helper (writes one column of doubles to tag_name.dat)
!
      CALL GETARG(1,MODE)
      CALL GETARG(2,NTYPE)
      CALL GETARG(3,CARG); READ(CARG,*)NIT
      CALL GETARG(4,CARG); READ(CARG,*)OUTBC
      CALL GETARG(5,INBM)
      CALL GETARG(6,HCHK)
      CALL GETARG(7,CARG); READ(CARG,*)I
      DOMATCH=(I .NE. 0)
      CALL GETARG(8,OUTTAG)
      IF(LEN_TRIM(OUTTAG) .EQ. 0)OUTTAG='mom'
!
! Optional arg 9: scattering albedo ESEC/CHI (default 0.6). The
! production m15 wind is electron-scattering dominated (albedo -> 1),
! where the lagged fixed-point map goes marginal; the sweep lets the
! driver measure the contraction rate against (1-albedo) instead of
! hardcoding the 0.4-absorption regime.
!
      ALB=0.6_LDP
      CALL GETARG(9,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)ALB
!
! Optional arg 10: CHI multiplier. Scaling opacity down moves the
! compressing rows toward the advection-dominated limit, where
! chi_hat loses its chi floor, Q -> GAM_FAC/(GAM_FAC sign) ~ +/-1,
! and the lagged map degenerates toward a frequency shift. The
! production worst cells sit in this thin compressing regime.
!
      CHISCL=1.0_LDP
      CALL GETARG(10,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)CHISCL
!
! Optional arg 11: log-space opacity-modulation amplitude across
! frequency (blanketing analogue): CHI_L *= EXP(CHIWAV*SIN(2 pi L/97)).
! Production line blanketing swings opacity by decades over tens of
! bins, so adjacent-frequency J fields differ steeply and a matched
! lagged read spans that gradient.
!
      CHIWAV=0.0_LDP
      CALL GETARG(11,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)CHIWAV
!
! Optional arg 12: per-sweep Eddington-drift amplitude. Production
! rebuilds the ray store from the formal solve every great
! iteration, so the lagged map's coefficients -- and its fixed
! point -- move between sweeps. PERT scales JNU/KNU/HNU/NNU by
! 1+PERT*sin(7.31*IT+depth-index phase) each sweep before the
! frequency loop: a deterministic moving target for measuring the
! residual floor a contracting map sustains under coefficient
! drift.
!
      PERT=0.0_LDP
      CALL GETARG(12,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)PERT
!
! Optional arg 13: rfst band albedo (default 1.0 = pure scattering
! CDS). Values below 1 restore absorption inside the band, which
! bounds the compressing solve and moves the fixture from the
! runaway-cap stall toward the smooth marginal-map stall.
!
      BALB=1.0_LDP
      CALL GETARG(13,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)BALB
!
! Optional arg 14: log-frequency grid spacing (default 5E-2). The
! production compressing rows sit at tight bins (line-edge doublets
! reach 2E-11); every advection coefficient carries 1/DLN, so a
! small DLN moves the fixture toward the production marginal-map
! regime without needing ALB=1 exactly.
!
      DLNA=5.0E-2_LDP
      CALL GETARG(14,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)DLNA
!
! Optional arg 15: outer velocity V0 (default 4000). CON_GAM~V/R
! scales every advection coefficient; production V in the CDS is
! ~2400-6400, but a larger V0 lifts |GAM| into the |Q|~1 regime at
! moderate opacity.
!
      V0A=4.0E3_LDP
      CALL GETARG(15,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)V0A
!
! Optional args 16-18: rfst/rfj refresh-schedule overrides. RFBLK<0
! leaves the production bulk-exit default alone; >=0 wires
! MOM_RF_BULK_EXIT via SET_MOM_RF_BULK (0 keeps the march running
! past the bulk tail so the marginal compressing rows get the full
! KN_M_MAXPASS window). RFMAXP/RFRLX>0 route through
! SET_MOM_KN_OPTS; unset knobs keep the module values.
!
      RFBLK=-1.0_LDP
      CALL GETARG(16,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFBLK
      RFMAXP=0.0_LDP
      CALL GETARG(17,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFMAXP
      RFRLX=0.0_LDP
      CALL GETARG(18,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFRLX
!
! Optional arg 19: rfst/rfj AA toggle. RFAA<0 leaves MOM_RF_AA at
! its module default (TRUE); 0 switches the pass-boundary Anderson
! acceleration off so the marginal-map drift keeps its geometric
! signature for the stall detector.
!
      RFAA=-1.0_LDP
      CALL GETARG(19,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFAA
!
! Optional arg 20: rfst/rfj pass-order selector. 0 keeps the
! production failing-arm schedule 'ALT'; 1 forces all-descending
! (the live-GS compressing read under MOM_DOM_CMP) and 2 forces
! all-ascending. Used to separate the Jacobi snapshot lag from the
! compressing map's own contractivity.
!
      RFDIR=0.0_LDP
      CALL GETARG(20,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFDIR
!
! Optional arg 21: midpoint N/H closure factor GF (default 0.5).
! On the USE_LAM_ES=F route GF comes from V13 ray moments, which
! are expansion-biased at sigma<-1 (comp_j_blank.f:891, F10): a
! sweep of GF probes whether the stall signature is sensitive to
! a mispriced compressing-branch closure. GET_MOMS_NON_REL fills
! the degenerate |H|<1E-3|N| case with SIGN(1E10,N), so GF=1E10
! also covers that branch.
!
      GF_ARG=0.5_LDP
      CALL GETARG(21,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)GF_ARG
!
! Optional arg 22: rfst/rfj live-GS toggle. RFGS>=0 wires
! MOM_RF_GS_DESC via SET_MOM_RF_OPTS (KEEPNEG untouched): 1 makes
! compressing rows read the live J_NEXT_FRQ on descending passes
! instead of the frozen JNEXT_RF_READ snapshot - the WP-M2
! candidate fix for the frozen-Jacobi non-contractive map.
!
      RFGS=-1.0_LDP
      CALL GETARG(22,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFGS
!
! Optional arg 23: rfst/rfj refresh relax factor. RFRLX2>=0 wires
! MOM_RF_RLX via SET_MOM_RF_RLX (refresh damp schedule only, the
! sweep write keeps JNEXT_RELAX): the WP-M5 fix pairs
! MOM_RF_GS_DESC with a full-step refresh.
!
      RFRLX2=-1.0_LDP
      CALL GETARG(23,CARG)
      IF(LEN_TRIM(CARG) .GT. 0)READ(CARG,*)RFRLX2
!
      NTYP6='N_ON_J'
      IF(NTYPE .EQ. 'gonly')NTYP6='G_ONLY'
      IF(NTYPE .EQ. 'mixed')NTYP6='MIXED'
      INBM6='DIFFUSION'
      IF(INBM .EQ. 'zfx')INBM6='ZERO_FLUX'
      IF(INBM .EQ. 'pnt')INBM6='PNT_SRCE'
      IF(INBM .EQ. 'schu')INBM6='SCHUSTER'
      HCHK6='NONE'
      IF(HCHK .EQ. 'av')HCHK6='AV_VAL'
      IF(HCHK .EQ. 'max')HCHK6='MAX_VAL'
      METH='LOGMON'
      COHR=.TRUE.
      VDOP_FRAC=1.0_LDP
      DBB=1.0E-3_LDP
      IC=2.0E-4_LDP
      IBSTAB=0.5_LDP
      DLN=DLNA
!
! Model R grid: log-spaced, R(1) outer.
!
      DO I=1,ND
        R(I)=4.0_LDP*(1.0_LDP/4.0_LDP)**(DBLE(I-1)/DBLE(ND-1))
      END DO
      DLNR=LOG(R(1)/R(ND))/DBLE(ND-1)
!
! sigma_true profile and integrated V: dlnV/dlnR = 1 + sigma.
! V(R(1)) = V0; integrate inward (lnR decreasing).
!
      NCASE=0
      SIGT=0.0_LDP
!
! Compression-band edges as grid fractions: cells 10-13 at ND=24.
! Defined in fractional depth so the seed case probes the same
! physical configuration under grid refinement (test E).
!
      IB1=NINT(9.0_LDP*DBLE(ND-1)/23.0_LDP)+1
      IB2=NINT(12.0_LDP*DBLE(ND-1)/23.0_LDP)+1
      V11MAX=0.0_LDP
      IF(MODE .EQ. 'homo' .OR. MODE .EQ. 'v11')THEN
        NCASE=1
      ELSE IF(MODE .EQ. 'beta')THEN
        NCASE=2
      ELSE IF(MODE .EQ. 'dec')THEN
        NCASE=3; SIGT=-2.0_LDP
      ELSE IF(MODE .EQ. 'vconst')THEN
        NCASE=4; SIGT=-1.0_LDP
      ELSE IF(MODE .EQ. 'comp1' .OR. MODE .EQ. 'ckpt' .OR.
     1        MODE .EQ. 'rfd' .OR. MODE .EQ. 'nurf' .OR.
     2        MODE .EQ. 'recall')THEN
        NCASE=5; SIGT(10:13)=-4.0_LDP
      ELSE IF(MODE .EQ. 'comp2')THEN
        NCASE=6; SIGT(7:9)=-4.0_LDP; SIGT(15:17)=-6.0_LDP
      ELSE IF(MODE .EQ. 'strong')THEN
        NCASE=7; SIGT(10:14)=-12.0_LDP
      ELSE IF(MODE .EQ. 'compw' .OR. MODE .EQ. 'nuw' .OR.
     1        MODE .EQ. 'rfst')THEN
        NCASE=8; SIGT(8:16)=-4.0_LDP
      ELSE IF(MODE .EQ. 'jump' .OR. MODE .EQ. 'rfj')THEN
        NCASE=9; SIGT(12:12)=-60.0_LDP
      ELSE IF(MODE .EQ. 'zero')THEN
        NCASE=10; SIGT(10:14)=-3.0_LDP
      ELSE IF(MODE .EQ. 'sigx')THEN
        NCASE=11
        DO I=8,16
          SIGT(I)=-4.0_LDP*MOD(I,2)
        END DO
      ELSE IF(MODE .EQ. 'dbl')THEN
        NCASE=12; SIGT(10:14)=-4.0_LDP
      ELSE IF(MODE .EQ. 'zsrc')THEN
        NCASE=13; SIGT(10:13)=-4.0_LDP
      ELSE IF(MODE .EQ. 'zedg')THEN
        NCASE=14; SIGT(10:13)=-4.0_LDP
      ELSE IF(MODE .EQ. 'seed')THEN
        NCASE=15; SIGT(IB1:IB2)=-4.0_LDP
      ELSE
        WRITE(*,'(A)')'RESULT unknown_mode=1'
        GOTO 999
      END IF
!
      V0=V0A
      V(1)=V0
      DO I=1,ND-1
        T1=0.5_LDP*(SIGT(I)+SIGT(I+1))
        V(I+1)=V(I)*EXP((1.0_LDP+T1)*LOG(R(I+1)/R(I)))
      END DO
      IF(MODE .EQ. 'beta')THEN
!
! Monotonic beta law: V = V0*(1-(RC/R))^BETA, RC=0.8, BETA=1.
! sigma = (R/V) dV/dR - 1 = RC*BETA/(R-RC) - 1 ... evaluated
! directly as (R/V)*dV/dR-1: dV/dR = V0*BETA*RC/R^2*(1-RC/R)^(BETA-1)
!
        RC=0.8_LDP; BETA=1.0_LDP
        DO I=1,ND
          V(I)=V0*(1.0_LDP-(RC/R(I)))**BETA
          T1=V0*BETA*(RC/R(I)/R(I))*
     1       (1.0_LDP-(RC/R(I)))**(BETA-1.0_LDP)
          SIGT(I)=R(I)*T1/V(I)-1.0_LDP
        END DO
      END IF
      SIG(1:ND)=MAX(SIGT(1:ND),-0.999_LDP)
!
! Opacity/emissivity: mildly radial, optically thick-ish inside.
! SOURCE=ETA/CHI ~ S0; ESEC/CHI ~ 0.6 (coherent scattering).
!
      S0=1.0E-4_LDP
      DO I=1,ND
        CHI(I)=CHISCL*2.0_LDP*(R(1)/R(I))**1.5_LDP
        ETA(I)=CHI(I)*S0*(1.0_LDP+0.2_LDP*(R(1)-R(I))/R(1))
        ESEC(I)=ALB*CHI(I)
        VDOP(I)=50.0_LDP
      END DO
!
! rfst: production depth structure -- the compressing band is the
! thick electron-scattering CDS (arg10 scales CHI inside the band,
! ESEC/CHI -> 1 there); the expanding rows outside keep the arg9
! albedo. Compressing rows then carry the marginal PSI/TB ~ 0.9
! lagged map while the expanding rows stay absorption-anchored and
! contractive.
!
      IF(MODE .EQ. 'rfst' .OR. MODE .EQ. 'rfj')THEN
        DO I=1,ND
          IF(SIGT(I) .LT. -1.0_LDP)THEN
            CHI(I)=CHI(I)*CHISCL
            ETA(I)=CHI(I)*S0*(1.0_LDP+0.2_LDP*(R(1)-R(I))/R(1))
            ESEC(I)=BALB*CHI(I)
          END IF
        END DO
      END IF
!
! zsrc: zero-emissivity compression band. The production iter-57 dump
! (MOM_J_NONMONO_ERRORS) shows the first negative-J cells sit at
! compressing depths with ETA ~ 1E-133; the lagged H/EPS debit then
! drives the assembled RHS negative. This mode reproduces that state:
! comp1 geometry with ETA~0 inside the band only.
!
      IF(MODE .EQ. 'zsrc')ETA(10:13)=1.0E-30_LDP
!
! zedg: same compression band and zero in-band source as zsrc, plus a
! frequency-dependent opacity edge: at even frequencies the band is
! opaque (CHI*EDGEFAC) so its J/H field is small inside the band and
! large outside. At an odd (transparent) frequency, compressing rows
! read the upwind redder store at that steeply-dropping field -- the
! production configuration behind the negative assembled RHS in
! MOM_J_NONMONO_ERRORS. ETA for the band stays zero at all
! frequencies.
!
      IF(MODE .EQ. 'zedg')ETA(10:13)=1.0E-30_LDP
!
! seed: comp1 geometry + zero in-band source. After one warm sweep the
! lagged redder store is overwritten with a crafted state matching the
! production failure record (MOM_J_NONMONO_ERRORS): J at the upwind
! bin ~0 inside the band while the lagged H carries a strong
! inward-growing negative flux there. A single transparent-frequency
! solve then shows whether the assembled RHS goes negative at the
! band, which is the verified first-failure mechanism of the it57
! model (negative XM_SAVE before the solve, not a DGTSV defect).
!
      IF(MODE .EQ. 'seed')ETA(IB1:IB2)=1.0E-30_LDP
!
! Frequency grid (blue -> red).
!
      FRQ(1)=1.0E15_LDP
      DO L=2,NCF
        FRQ(L)=FRQ(L-1)*EXP(-DLN)
      END DO
!
! dbl: bound-free-edge doublet probe. FRQ(4) sits 1E-11 below FRQ(3)
! in log frequency, so the L=3 solve's compressing rows get
! MOM_DLN ~ 1E-11 (floored at MOM_DLN_MIN) and expanding rows at
! L=4 get the same near-null backward spacing.
!
      IF(MODE .EQ. 'dbl')FRQ(4)=FRQ(3)*EXP(-1.0E-11_LDP)
!
! Seed MOD_RAY_MOM_STORE on the model grid (ND_STORE=ND -> identity
! interpolation branch in GET_MOMS_NON_REL, HN_DEF_ON_NODES=F -> the
! midpoint H,N are read straight from the stores).
!
      FJ=1.0_LDP/3.0_LDP
      HF=0.30_LDP
      GF=GF_ARG
      ND_STORE=ND
      ALLOCATE(R_STORE(ND),RMID_STORE(ND),EXT_RMID_STORE(ND+2),
     1 JNU_STORE(ND),HNU_STORE(ND),NNU_STORE(ND),KNU_STORE(ND),
     2 GAM_REL_STORE(ND),STAT=IOS)
      IF(IOS .NE. 0)GOTO 999
      R_STORE(1:ND)=R(1:ND)
      DO I=1,ND-1
        RMID_STORE(I)=0.5_LDP*(R(I)+R(I+1))
      END DO
      RMID_STORE(ND)=0.9_LDP*R(ND)
      DO I=1,ND
        JNU_STORE(I)=S0*(R(1)/R(I))**1.2_LDP
        KNU_STORE(I)=FJ*JNU_STORE(I)
      END DO
      DO I=1,ND
        T1=0.5_LDP*(JNU_STORE(MAX(1,I-1))+JNU_STORE(MIN(ND,I)))
        HNU_STORE(I)=HF*T1
        NNU_STORE(I)=GF*HNU_STORE(I)
      END DO
      GAM_REL_STORE=0.0_LDP
      EXT_RMID_STORE(1)=R(1)+0.1_LDP
      EXT_RMID_STORE(2:ND+1)=RMID_STORE(1:ND)
      EXT_RMID_STORE(ND+2)=0.8_LDP*R(ND)
      HN_DEF_ON_NODES=.FALSE.
!
      HNU_AT_OB=HNU_STORE(1)
      NNU_AT_OB=NNU_STORE(1)
      HNU_AT_IB=HF*JNU_STORE(ND)
      NNU_AT_IB=GF*HNU_AT_IB
      HBC=HNU_AT_OB/JNU_STORE(1)
      NBC=NNU_AT_OB/JNU_STORE(1)
      IN_HBC=0.4_LDP
      JPLUS_IB=JNU_STORE(ND)*0.9_LDP
      JMIN_IB=JNU_STORE(ND)*0.1_LDP
!
      MOM_FREQ_MATCH=DOMATCH
      CONV_DIAG_LVL=2
!
      IF(MODE .EQ. 'dbl')
     1   OPEN(UNIT=305,FILE=TRIM(OUTTAG)//'_jprof.txt',
     2        STATUS='REPLACE')
!
      WRITE(*,'(A,I4,A,I4)')'GRID ND=',ND,' NCF=',NCF
      WRITE(*,'(A,ES14.6,A,ES14.6)')'V(1)=',V(1),' V(ND)=',V(ND)
      WRITE(*,'(A,ES14.6)')'SIGT_MIN=',MINVAL(SIGT)
!
! Sweeps: each call solves one frequency; the lagged J_NEXT/H_NEXT
! stores converge across NIT great iterations.
!
      IF(MODE .EQ. 'ckpt')THEN
!
! Checkpoint round-trip: run 2 sweeps to populate J_NEXT/H_NEXT,
! dump, scramble the stores, restore with the matching tag, compare.
!
        NIT=2
      END IF
      DO IT=1,NIT
        IF(PERT .NE. 0.0_LDP .AND. IT .GT. 1)THEN
          DO I=1,ND
            PFAC=1.0_LDP+PERT*SIN(7.31_LDP*IT+0.37_LDP*I)
            JNU_STORE(I)=JNU_STORE(I)*PFAC
            KNU_STORE(I)=FJ*JNU_STORE(I)
          END DO
          DO I=1,ND
            T1=0.5_LDP*(JNU_STORE(MAX(1,I-1))+JNU_STORE(MIN(ND,I)))
            HNU_STORE(I)=HF*T1
            NNU_STORE(I)=GF*HNU_STORE(I)
          END DO
          HNU_AT_OB=HNU_STORE(1)
          NNU_AT_OB=NNU_STORE(1)
          HNU_AT_IB=HF*JNU_STORE(ND)
          NNU_AT_IB=GF*HNU_AT_IB
          HBC=HNU_AT_OB/JNU_STORE(1)
          NBC=NNU_AT_OB/JNU_STORE(1)
          JPLUS_IB=JNU_STORE(ND)*0.9_LDP
          JMIN_IB=JNU_STORE(ND)*0.1_LDP
        END IF
        DO L=1,NCF
          FREQ=FRQ(L)
          FIDX=L
          IF(L .LT. NCF)THEN
            NU_NXT=FRQ(L+1)
          ELSE
            NU_NXT=-1.0_LDP
          END IF
          CALL SET_MOM_NU_NEXT(NU_NXT)
          INIT=(IT .EQ. 1 .AND. L .EQ. 1)
          NEW_FREQ=.TRUE.
! dbl: pass the true backward spacing LOG(F_{L-1}/F_L) so expanding
! rows at the doublet see the same near-null spacing the
! compressing rows derive from NU_NEXT.
          IF(MODE .EQ. 'dbl' .AND. L .GT. 1)THEN
            DLNC=LOG(FRQ(L-1)/FRQ(L))
          ELSE
            DLNC=DLN
          END IF
!
!         zedg: opaque compression band at even frequencies only.
!
          CHI_L(1:ND)=CHI(1:ND)
          IF(CHIWAV .NE. 0.0_LDP)
     1       CHI_L(1:ND)=CHI_L(1:ND)*EXP(CHIWAV*
     2         SIN(6.283185307_LDP*DBLE(L)/97.0_LDP))
          IF(MODE .EQ. 'zedg' .AND. MOD(L,2) .EQ. 0)
     1       CHI_L(10:13)=CHI(10:13)*1.0E+4_LDP
          CALL MOM_J_CMF_NONMONO_V1(ETA,CHI_L,ESEC,
     1         V,SIG,SIGT,R,
     2         JNU,RSQH,
     3         VDOP,VDOP_FRAC,
     4         FREQ,DLNC,
     5         INBM6,DBB,IC,IBSTAB,
     6         NTYP6,HCHK6,METH,COHR,OUTBC,
     7         INIT,NEW_FREQ,0,1,ND,FIDX,NCF)
!
!         recall: mimic the comp_j_blank INACCURATE loop, which
!         re-calls the solver at the same frequency with NEW_FREQ=F.
!         The sweep census (GS_NBAD5 etc.) must stay within the cell
!         count -- a per-call census would recount compressing rows
!         (frozen store -> identical drift) and inflate past NTOT.
!
          IF(MODE .EQ. 'recall')THEN
            CALL MOM_J_CMF_NONMONO_V1(ETA,CHI_L,ESEC,
     1           V,SIG,SIGT,R,
     2           JNU,RSQH,
     3           VDOP,VDOP_FRAC,
     4           FREQ,DLNC,
     5           INBM6,DBB,IC,IBSTAB,
     6           NTYP6,HCHK6,METH,COHR,OUTBC,
     7           .FALSE.,.FALSE.,0,1,ND,FIDX,NCF)
          END IF
!
!         v11: A/B against the stock upstream solver on the same
!         monotonic input -- the non-mono discretisation must reduce
!         to V11 where sigma_true >= -1 everywhere.
!
          IF(MODE .EQ. 'v11')THEN
            CALL MOM_J_CMF_V11(ETA,CHI_L,ESEC,
     1         V,SIG,R,
     2         JV,HSC,
     3         VDOP,VDOP_FRAC,
     4         FREQ,DLNC,
     5         INBM6,DBB,IC,IBSTAB,
     6         NTYP6,HCHK6,METH,COHR,OUTBC,
     7         INIT,NEW_FREQ,0,1,ND)
            T1=MAXVAL(ABS(JNU))
            V11REL=MAXVAL(ABS(JV-JNU)/
     1         (ABS(JNU)+1.0E-4_LDP*MAX(T1,1.0E-30_LDP)))
            V11MAX=MAX(V11MAX,V11REL)
          END IF
          IF(MODE .EQ. 'dbl' .AND. IT .EQ. NIT)
     1         WRITE(305,'(I4,999ES16.8)')L,JNU(1:ND)
        END DO
        WRITE(*,'(A,I3,A,ES14.6,A,I8)')'IT=',IT,
     1   ' MOM_GS_MAXCHG=',MOM_GS_MAXCHG,' NCROW=',
     2   MOM_SW_NCROW+MOM_SW_NROW
      END DO
!
      IF(MODE .EQ. 'v11')
     1   WRITE(*,'(A,ES16.8E3)')'V11_MAXREL=',V11MAX
!
      IF(MODE .EQ. 'seed')THEN
!
! Inject the production failure state: at the upwind redder bin (L+1)
! J is ~0 inside the band while H is a negative flux growing in
! magnitude inward across the band -- so the lagged H-divergence
! HS(I)*H(I)-HS(I-1)*H(I-1) is a debit at every band row and the
! local source is zero. Then solve the transparent frequency L=1
! once and report the assembled RHS and solved field in the band.
!
        J_NEXT_FRQ(1:ND_NM,2)=0.0_LDP
        DO I=1,ND_NM-1
          T1=MIN(MAX((DBLE(I)/DBLE(ND_NM)-0.44_LDP)/0.22_LDP,
     1       0.0_LDP),1.0_LDP)
          H_NEXT_FRQ(I,2)=-2.0E-4_LDP*T1
        END DO
!
! Optically thick band at the seeded frequency too: the production
! negative-J rows were absorption-dominated (|TA|+|TC| << TB), so the
! diffusion coupling cannot lift the negative RHS back to J>0.
!
        CHI_L(1:ND)=CHI(1:ND)
        CHI_L(IB1:IB2)=CHI(IB1:IB2)*1.0E+4_LDP
        FREQ=FRQ(1); FIDX=1; NU_NXT=FRQ(2)
        CALL SET_MOM_NU_NEXT(NU_NXT)
        INIT=.FALSE.; NEW_FREQ=.TRUE.
        CALL MOM_J_CMF_NONMONO_V1(ETA,CHI_L,ESEC,
     1       V,SIG,SIGT,R,
     2       JNU,RSQH,
     3       VDOP,VDOP_FRAC,
     4       FREQ,DLN,
     5       INBM6,DBB,IC,IBSTAB,
     6       NTYP6,HCHK6,METH,COHR,OUTBC,
     7       INIT,NEW_FREQ,0,1,ND,FIDX,NCF)
        WRITE(*,'(A,ES16.8E3)')'SEED_XMSAVE_MIN=',MINVAL(XM_SAVE_NM)
        WRITE(*,'(A,ES16.8E3)')'SEED_RSQJ_MIN=',MINVAL(RSQJNU_NM)
        WRITE(*,'(A,I10)')'SEED_XMNEG=',MOM_XMNEG
      END IF
!
      IF(MODE .EQ. 'rfd' .OR. MODE .EQ. 'nurf' .OR.
     1   MODE .EQ. 'nuw')THEN
!
! Test D: refresh-order invariance. The sweep above leaves
! J_NEXT_FRQ populated and KN_M_READY set; re-drive the captured
! state through MOM_KNOP_REFRESH_V1 under forced descending and
! ascending pass orders from the same lagged store. Both orders
! converge the same lagged map, so the committed stores must agree
! to the refresh exit tolerance -- a direction-dependent fixed
! point would mean the traversal order itself is doing physics.
!
! nurf adds the MOM_NUSOLV arm: a third refresh with the implicit
! frequency-chain pass enabled. It must land on the same fixed
! point -- the captured equations are the swept row's own, so the
! composite map's solution is identical -- which the store
! comparison measures directly.
!
        ALLOCATE(JRF0(ND_NM,NCF),HRF0(ND_NM,NCF),
     1   JDSC(ND_NM,NCF),JASC(ND_NM,NCF),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        JRF0(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
        HRF0(1:ND_NM,1:NCF)=H_NEXT_FRQ(1:ND_NM,1:NCF)
        CALL SET_MOM_RF_DIR_OPT('DESC',.TRUE.,.TRUE.,.TRUE.)
        CALL MOM_KNOP_REFRESH_V1()
        JDSC(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
        IF(MODE .NE. 'rfd')
     1   WRITE(*,'(A,L2)')'RF_DESC_CONV=',MOM_RF_CONV
        J_NEXT_FRQ(1:ND_NM,1:NCF)=JRF0
        H_NEXT_FRQ(1:ND_NM,1:NCF)=HRF0
        IF(MODE .EQ. 'rfd')THEN
          CALL SET_MOM_RF_DIR_OPT('ASC ',.TRUE.,.TRUE.,.TRUE.)
          CALL MOM_KNOP_REFRESH_V1()
          JASC(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
          T1=MAXVAL(ABS(JRF0))
          WRITE(*,'(A,ES16.8E3)')'RF_ASCDESC_MAXDIFF=',
     1      MAXVAL(ABS(JASC-JDSC))/MAX(T1,1.0E-30_LDP)
          WRITE(*,'(A,ES16.8E3)')'RF_ASCDESC_MAXREL=',
     1      MAXVAL(ABS(JASC-JDSC)/
     2      (ABS(JDSC)+1.0E-4_LDP*MAX(T1,1.0E-30_LDP)))
        ELSE
!
! MOM_NUSOLV arm ('nurf' mode): the chain solve runs after each pass
! inside the same refresh driver. JDSC is the NUSOLV=F fixed point;
! the NUSOLV=T store must agree -- the chain solve satisfies the
! captured equations exactly, so any bulk disagreement is a
! capture/consistency defect, not physics.
!
          ALLOCATE(JNUV(ND_NM,NCF),STAT=IOS)
          IF(IOS .NE. 0)GOTO 999
          CALL SET_MOM_NUSOLV(.TRUE.)
          CALL SET_MOM_RF_DIR_OPT('ALT ',.TRUE.,.TRUE.,.TRUE.)
          CALL MOM_KNOP_REFRESH_V1()
          JNUV(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
          T1=MAXVAL(ABS(JRF0))
          WRITE(*,'(A,ES16.8E3)')'RF_NUSOLV_MAXDIFF=',
     1      MAXVAL(ABS(JNUV-JDSC))/MAX(T1,1.0E-30_LDP)
          WRITE(*,'(A,ES16.8E3)')'RF_NUSOLV_MAXREL=',
     1      MAXVAL(ABS(JNUV-JDSC)/
     2      (ABS(JDSC)+1.0E-4_LDP*MAX(T1,1.0E-30_LDP)))
          WRITE(*,'(A,L2)')'RF_NUSOLV_CONV=',MOM_RF_CONV
          WRITE(*,'(A,ES16.8E3)')'NU_PASS_CHG=',MOM_NU_CHG_LAST
          WRITE(*,'(A,I10)')'NU_NRING_T=',MOM_NU_NRING_T
          WRITE(*,'(A,I10)')'NU_NSING_T=',MOM_NU_NSING_T
          WRITE(*,'(A,I10)')'NU_NFRZ_T=',MOM_NU_NFRZ_T
          WRITE(*,'(A,I10)')'NU_NNEG_T=',MOM_NU_NNEG_T
          WRITE(*,'(A,F4.1)')'NU_FINISHED=',
     1      MERGE(1.0_LDP,0.0_LDP,.NOT. MOM_KN_REFRESHING)
          CALL SET_MOM_NUSOLV(.FALSE.)
          DEALLOCATE(JNUV)
!
! Block-GMRES arm: same restart, same refresh driver, but the captured
! rows are solved as one wavelength-depth system. JASC is reused as
! the block-solve store. The fixed point must again agree with the
! swept system; this arm exercises the production capture, inactive-
! column pinning, GMRES write-back and precedence wiring together.
!
          J_NEXT_FRQ(1:ND_NM,1:NCF)=JRF0
          H_NEXT_FRQ(1:ND_NM,1:NCF)=HRF0
          CALL SET_MOM_BLKGMRES(.TRUE.)
          CALL SET_MOM_NUSOLV(.TRUE.)
          CALL SET_MOM_RF_DIR_OPT('ALT ',.TRUE.,.TRUE.,.TRUE.)
          CALL MOM_KNOP_REFRESH_V1()
          JASC(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
          T1=MAXVAL(ABS(JRF0))
          WRITE(*,'(A,ES16.8E3)')'RF_BLK_MAXDIFF=',
     1      MAXVAL(ABS(JASC-JDSC))/MAX(T1,1.0E-30_LDP)
          WRITE(*,'(A,ES16.8E3)')'RF_BLK_MAXREL=',
     1      MAXVAL(ABS(JASC-JDSC)/
     2      (ABS(JDSC)+1.0E-4_LDP*MAX(T1,1.0E-30_LDP)))
          WRITE(*,'(A,L2)')'RF_BLK_CONV=',MOM_RF_CONV
          WRITE(*,'(A,ES16.8E3)')'BLK_PASS_CHG=',MOM_BLK_CHG_LAST
          WRITE(*,'(A,I10)')'BLK_NITER=',MOM_BLK_NITER
          WRITE(*,'(A,ES16.8E3)')'BLK_RESID=',MOM_BLK_RESID
          WRITE(*,'(A,I10)')'BLK_NFRZ=',MOM_BLK_NFRZ
          WRITE(*,'(A,I10)')'BLK_NNEG=',MOM_BLK_NNEG
          WRITE(*,'(A,I10)')'BLK_NCAP=',MOM_BLK_NCAP
          WRITE(*,'(A,F4.1)')'BLK_FINISHED=',
     1      MERGE(1.0_LDP,0.0_LDP,.NOT. MOM_KN_REFRESHING)
          CALL SET_MOM_BLKGMRES(.FALSE.)
          CALL SET_MOM_NUSOLV(.FALSE.)
        END IF
        WRITE(*,'(A,F4.1)')'RF_FINISHED=',
     1    MERGE(1.0_LDP,0.0_LDP,.NOT. MOM_KN_REFRESHING)
        DEALLOCATE(JRF0,HRF0,JDSC,JASC)
      END IF
!
      IF(MODE .EQ. 'rfst' .OR. MODE .EQ. 'rfj')THEN
!
! Compressing-branch stall fixture (checklists/use_lam_es_eq_f.md
! WP-M1). Production depth structure: a thin, scattering
! compressing band (jump cell sigt=-60; arg10 scales CHI inside it,
! arg13 sets its albedo) inside an absorption-anchored column (arg9
! albedo ~0). At tight log-frequency spacing (arg14 ~1E-4) and high
! outer velocity (arg15 ~1E5) the compressing rows reach the
! production marginal-map regime: the nu-census a-value (max
! |lagged J coupling|/|diagonal|) measures ~0.96 here vs 9.947E-01
! in the m15 diagnosis, with expanding rows absorption-anchored and
! contractive. One ALT refresh -- the production failing-arm
! schedule, no implicit arms -- must print the production signature
! "Compressing direction stalled - dropping it" while the expanding
! branch converges and carries the refresh to a converged exit.
!
        IF(RFBLK .GE. 0.0_LDP)
     1    CALL SET_MOM_RF_BULK(RFBLK .GT. 0.5_LDP)
        IF(RFAA .GE. 0.0_LDP)
     1    CALL SET_MOM_RF_AA(RFAA .GT. 0.5_LDP)
        IF(RFMAXP .GT. 0.0_LDP .OR. RFRLX .GT. 0.0_LDP)
     1    CALL SET_MOM_KN_OPTS(KN_M_TOL,
     2      MERGE(NINT(RFMAXP),KN_M_MAXPASS,RFMAXP .GT. 0.0_LDP),
     3      MERGE(RFRLX,JNEXT_RELAX,RFRLX .GT. 0.0_LDP),
     4      MOM_AA_M,MOM_AA_OUT,MOM_AA_OUT_M)
!
! MOM_AUG_RES is capture-only (the solve path is unchanged): arming
! it makes the per-pass nu-census accumulate over the ALT march, so
! RFST_AMAX below is the same lagged-coupling ratio the production
! diagnosis reported (nu-census amax 9.947E-01).
!
        CALL SET_MOM_AUG_RES(.TRUE.)
        IF(RFGS .GE. 0.0_LDP)
     1    CALL SET_MOM_RF_OPTS(RFGS .GT. 0.5_LDP,MOM_RF_KEEPNEG)
        IF(RFRLX2 .GE. 0.0_LDP)
     1    CALL SET_MOM_RF_RLX(RFRLX2)
        IF(RFDIR .EQ. 1.0_LDP)THEN
          CALL SET_MOM_RF_DIR_OPT('DESC',.TRUE.,.TRUE.,.TRUE.)
        ELSE IF(RFDIR .EQ. 2.0_LDP)THEN
          CALL SET_MOM_RF_DIR_OPT('ASC ',.TRUE.,.TRUE.,.TRUE.)
        ELSE
          CALL SET_MOM_RF_DIR_OPT('ALT ',.TRUE.,.TRUE.,.TRUE.)
        END IF
!
! Pre-refresh J snapshot: on exit, the retained field's per-cell
! distance back to this store localizes which rows the march could
! not settle -- the compressing stall should live on the marginal
! (census a -> 1) band cells.
!
        ALLOCATE(JRF0(ND_NM,NCF),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        JRF0(1:ND_NM,1:NCF)=J_NEXT_FRQ(1:ND_NM,1:NCF)
        CALL MOM_KNOP_REFRESH_V1()
        WRITE(*,'(A,F4.1)')'RFST_CONV=',
     1    MERGE(1.0_LDP,0.0_LDP,MOM_RF_CONV)
        WRITE(*,'(A,ES16.8E3)')'RFST_DCHG=',MOM_RF_DCHG
        WRITE(*,'(A,I10)')'RFST_N5=',MOM_RF_N5
        WRITE(*,'(A,I10)')'RFST_N50=',MOM_RF_N50
        WRITE(*,'(A,F4.1)')'RFST_FINISHED=',
     1    MERGE(1.0_LDP,0.0_LDP,.NOT. MOM_KN_REFRESHING)
        WRITE(*,'(A,ES16.8E3)')'RFST_AMAX=',NUC_AMAX
        WRITE(*,'(A,I10)')'RFST_NA9=',NUC_NA9
        WRITE(*,'(A,I10)')'RFST_NA99=',NUC_NA99
        WRITE(*,'(A,I10)')'RFST_XSNEG=',NUC_XSNEG
        WRITE(*,'(A,I10)')'RFST_NROW=',NUC_NROW
!
! Net drift cells (I,L): distance from the pre-refresh store to the
! retained field, ranked by the relative metric |D|/max(|J_prev|,
! 1E-4*column peak). Converged expanding rows legitimately moved
! far; a stalled compressing cell sits mid-range. The depth index
! tells whether the moving cells live on the compressing band.
!
        DO NCASE=1,5
          DMAX=0.0_LDP; IB1=0; IB2=0
          DO L=1,NCF
            T1=MAXVAL(ABS(JRF0(1:ND_NM,L)))
            DO I=1,ND_NM
              T2=ABS(J_NEXT_FRQ(I,L)-JRF0(I,L))/
     1          MAX(ABS(JRF0(I,L)),1.0E-4_LDP*T1)
              IF(T2 .GT. DMAX)THEN
                DMAX=T2; IB1=I; IB2=L
              END IF
            END DO
          END DO
          IF(IB1 .EQ. 0)EXIT
          WRITE(*,'(A,I2,A,2I8,A,ES12.4)')'RFST_TOP',NCASE,
     1      '_CELL=',IB1,IB2,' DRIFT=',DMAX
          JRF0(IB1,IB2)=J_NEXT_FRQ(IB1,IB2)
        END DO
        DEALLOCATE(JRF0)
        CALL SET_MOM_AUG_RES(.FALSE.)
      END IF
!
      IF(MODE .EQ. 'ckpt')THEN
        ALLOCATE(CKJ(ND_NM,NCF),CKH(ND_NM,NCF),CKQ(ND_NM,NCF),STAT=IOS)
        IF(IOS .NE. 0)GOTO 999
        CKJ=J_NEXT_FRQ(1:ND_NM,1:NCF)
        CKH=H_NEXT_FRQ(1:ND_NM,1:NCF)
! Seed a nonzero quarantine-age pattern so the VER5 age record is
! exercised: the round-trip must restore it and a rejected restore
! must reset it alongside the field.
        DO L=1,NCF
          DO I=1,ND_NM
            MOM_QUAR_AGE(I,L)=MOD(I+L,5)
          END DO
        END DO
        CKQ=MOM_QUAR_AGE(1:ND_NM,1:NCF)
        CALL SET_MOM_GS_TAG(7)
        CALL DUMP_MOM_JNEXT_GS(7,R,V,SIG,ND,291)
        J_NEXT_FRQ=-7.0_LDP
        H_NEXT_FRQ=-7.0_LDP
        MOM_QUAR_AGE=0_2
        CALL RESTORE_MOM_JNEXT_GS(R,V,ND)
        DMAX=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,1:NCF)-CKJ)),
     1       MAXVAL(ABS(H_NEXT_FRQ(1:ND_NM,1:NCF)-CKH)))
        WRITE(*,'(A,L2)')'CKPT_OK=',MOM_J_GS_OK
        WRITE(*,'(A,ES16.8)')'CKPT_MAXDIFF=',DMAX
        WRITE(*,'(A,L2)')'CKPT_AGE_OK=',
     1       ALL(MOM_QUAR_AGE(1:ND_NM,1:NCF) .EQ. CKQ)
!
! Fingerprint mismatch: perturbed V endpoint must reject and zero.
!
        J_NEXT_FRQ=1.0_LDP
        V(ND)=V(ND)*1.01_LDP
        CALL RESTORE_MOM_JNEXT_GS(R,V,ND)
        WRITE(*,'(A,L2)')'CKPT_FP_OK=',MOM_J_GS_OK
        WRITE(*,'(A,ES16.8)')'CKPT_FP_STORE=',
     1   MAXVAL(ABS(J_NEXT_FRQ))
        WRITE(*,'(A,L2)')'CKPT_FP_AGE0=',
     1   ALL(MOM_QUAR_AGE .EQ. 0_2)
!
! Tag mismatch: reject.
!
        V(ND)=V(ND)/1.01_LDP
        J_NEXT_FRQ=1.0_LDP
        CALL SET_MOM_GS_TAG(8)
        CALL RESTORE_MOM_JNEXT_GS(R,V,ND)
        WRITE(*,'(A,L2)')'CKPT_TAG_OK=',MOM_J_GS_OK
        WRITE(*,'(A)')'RESULT_OK'
        GOTO 999
      END IF
!
! Result counters.
!
      WRITE(*,'(A,I6)')'ND_NM=',ND_NM
      WRITE(*,'(A,ES14.6)')'GS_MAXCHG=',MOM_GS_MAXCHG
      WRITE(*,'(A,I10)')'GS_NBAD5=',MOM_GS_NBAD5
      WRITE(*,'(A,I10)')'GS_NBAD50=',MOM_GS_NBAD50
      WRITE(*,'(A,I10)')'GS_NTOT=',MOM_GS_NTOT
      WRITE(*,'(A,I10)')'GS_NQUAR=',MOM_GS_NQUAR
      WRITE(*,'(A,I10)')'NEGJ_LRG=',MOM_NEGJ_LRG
      WRITE(*,'(A,I10)')'NEGJ_T=',MOM_N_NEGJ_T
      WRITE(*,'(A,I10)')'N_DGTSV=',MOM_N_DGTSV
      WRITE(*,'(A,I10)')'N_NONF=',MOM_N_NONF
      WRITE(*,'(A,I10)')'N_RUNAWAY=',MOM_N_RUNAWAY
      WRITE(*,'(A,I10)')'TB_NONPOS=',MOM_TB_NONPOS
      WRITE(*,'(A,I10)')'INEXT_CAP=',MOM_INEXT_CAPPED
      WRITE(*,'(A,ES14.6)')'MIN_DR=',MOM_MIN_DR
      WRITE(*,'(A,ES14.6)')'EQ_RES=',MOM_EQ_RES
      WRITE(*,'(A,ES14.6)')'LAGB_MIN=',MOM_LAGB_MIN
      WRITE(*,'(A,I10)')'LAGB_CNT=',MOM_LAGB_CNT
      WRITE(*,'(A,I10)')'EPS_CAP=',MOM_EPS_CAP
      WRITE(*,'(A,ES14.6)')'JMIN=',MINVAL(JNU(1:ND))
      WRITE(*,'(A,ES14.6)')'JMAX=',MAXVAL(JNU(1:ND))
!
! Localize the drift argmax: the metric is relative to
! max(|J_old|, 1E-4*column peak), so an argmax sitting on a
! near-empty cell means the reported drift is a metric floor, not a
! physically significant field error.
!
      IF(ALLOCATED(J_NEXT_FRQ) .AND. MOM_GS_ARGI .GT. 0 .AND.
     1   MOM_GS_ARGL .GT. 0)THEN
        WRITE(*,'(A,2I8)')'ARGCELL=',MOM_GS_ARGI,MOM_GS_ARGL
        WRITE(*,'(A,ES14.6)')'ARGJ=',
     1   J_NEXT_FRQ(MOM_GS_ARGI,MOM_GS_ARGL)
        WRITE(*,'(A,ES14.6)')'ARGCOLMAX=',
     1   MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,MOM_GS_ARGL)))
      END IF
      WRITE(*,'(A,ES14.6)')'STORE_MINJ=',MINVAL(J_NEXT_FRQ)
      WRITE(*,'(A,ES14.6)')'STORE_MAXJ=',MAXVAL(J_NEXT_FRQ)
!
! Dumps: internal coefficient arrays on the expanded grid.
!
      OPEN(UNIT=301,FILE=TRIM(OUTTAG)//'_nm.txt',STATUS='REPLACE')
      WRITE(301,'(A)')'# I R_NM V_NM SIGT_NM K_ON_J NMJ NMH '//
     1 'CHI_NM ETA_NM ESEC_NM CHIHAT GAMS Q_NM W_NM EPS_NM '//
     2 'PSI PSIPREV HU HL HS TA TB TC DTAU DTAUONQ QF RSQJ RSQH '//
     3 'XMSAVE RSQJP RSQHP'
      DO I=1,ND_NM
        T1=0.0_LDP; T2=0.0_LDP
        IF(I .LE. ND_NM-1)THEN
          T1=NMID_ON_J(I); T2=NMID_ON_HMID(I)
        END IF
        WRITE(301,'(I6,31ES25.16)')I,R_NM(I),V_NM(I),SIGMA_TRUE_NM(I),
     1   K_ON_J(I),T1,T2,CHI_NM(I),ETA_NM(I),ESEC_NM(I),
     2   CHI_HAT_NM(I),GAM_SIGNED_NM(I),Q_NM(I),
     3   W_NM(MIN(I,ND_NM-1)),EPS_NM(MIN(I,ND_NM-1)),
     4   PSI_NM(I),PSIPREV_NM(I),
     5   HU_NM(MIN(I,ND_NM-1)),HL_NM(MIN(I,ND_NM-1)),
     6   HS_NM(MIN(I,ND_NM-1)),
     7   TA_NM(I),TB_NM(I),TC_NM(I),
     8   DTAU_NM(MIN(I,ND_NM-1)),DTAUONQ_NM(I),QF_NM(I),
     9   RSQJNU_NM(I),RSQHNU_NM(MIN(I,ND_NM-1)),
     1   XM_SAVE_NM(I),RSQJNU_PREV_NM(I),RSQHNU_PREV_NM(I)
      END DO
      CLOSE(301)
!
      OPEN(UNIT=302,FILE=TRIM(OUTTAG)//'_jnext.txt',
     1     STATUS='REPLACE')
      DO I=1,ND_NM
        WRITE(302,'(I6,16ES25.16)')I,
     1   (J_NEXT_FRQ(I,L),L=1,NCF)
      END DO
      CLOSE(302)
!
      OPEN(UNIT=303,FILE=TRIM(OUTTAG)//'_hnext.txt',
     1     STATUS='REPLACE')
      DO I=1,ND_NM-1
        WRITE(303,'(I6,16ES25.16)')I,
     1   (H_NEXT_FRQ(I,L),L=1,NCF)
      END DO
      CLOSE(303)
!
      OPEN(UNIT=304,FILE=TRIM(OUTTAG)//'_sm.txt',STATUS='REPLACE')
      WRITE(304,'(A)')'# I R V SIGT CHI JNU J_INDX'
      DO I=1,ND
        WRITE(304,'(I6,5ES25.16,I6)')I,R(I),V(I),SIGT(I),CHI(I),
     1   JNU(I),J_INDX_NM(I)
      END DO
      CLOSE(304)
!
      WRITE(*,'(A)')'RESULT_OK'
999   CONTINUE
      END
