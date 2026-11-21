! Replay driver for the KNRDUMP ray deck (phase 2.22 stall audit).
! Reads KNRDUMP_L<LS>_r<ref>.bin written by KN_RAY_DUMP inside
! KNOP_RAY_SOLVE_V1 (plane/fg_j_cmf_nonmono_v1.f): the deck carries
! every module input the production solve consumes -- ray geometry,
! GAM profile, captured CHI/ETA cubics, frequency tables, boundary
! constants, both lagged store columns and the seeded P field --
! plus the production outcome as a verification record.  The driver
! pokes the state into the same module globals and calls the
! production KNOP_RAY_SOLVE_V1 unchanged, so experiments run the
! exact ordered march and controller, not a reimplementation.
!
! Usage: drv_knop_replay <deck.bin> [mode] [value]
!   base            one production solve with the deck's own
!                   controller settings; RES/IXIT/NPO must match
!                   the appended production record.
!   loop <n>        n one-pass solves (KN_RAY_MAXPASS=1 per call);
!                   after each call the exported post-pass P field
!                   is written back into the store, replaying the
!                   seed->pass->seed chain under caller control.
!   om <v>          force KN_RAY_OMEGA=v before the run (combine
!                   with loop for a truly fixed-omega schedule --
!                   OMEGA_RAY is re-initialised from KN_RAY_OMEGA
!                   at every solve entry, so in-loop adaptation
!                   cannot drift it).
!   btol <x>        set KN_RAY_BULK_TOL=x.
!   stall <x>       set KN_RAY_STALL=x.
!   maxp <n>        override KN_RAY_MAXPASS.
!   seed cold       zero both store columns before the run.
!   seed scale <f>  multiply both store columns by f.
!   verify          extra unrelaxed pass on the final state
!                   (omega=1) -> the independent map defect.
!   mapd <file>     one unrelaxed, undamped pass on an arbitrary
!                   full-field candidate (KN_REPLAY_END format:
!                   committed P then M). Reports the raw map
!                   defect R(X)-X on the production DCELL scale
!                   and writes KN_MAPD_RX.bin for offline
!                   fixed-mask analysis.
!   p2=<v>          replay-only override of FG_ORD_P2DAMP,
!                   v in {0,1,on,off,true,false,yes,no,t,f,y,n}.
!                   May appear in any argument slot >=3 and
!                   composes with every mode (incl. omf). The
!                   production default and VADAT semantics are
!                   untouched; a malformed p2 token stops the run.
!
      PROGRAM DRV_KNOP_REPLAY
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      USE SET_KIND_MODULE
      USE FG_J_CMF_NONMONO_MOD_V1
      USE ITER_LEDGER_MOD, ONLY: LEDG_N_RF
      USE CONTROL_VARIABLE_MOD, ONLY: GS_FOREIGN_SEED
      IMPLICIT NONE
!
      INTEGER NI,NDX,NCF,NDD,NCPD,NXD,KNC,IMETH,FBITS,KMAXP,KCNT
      INTEGER NADD,NCPS,KNPD
      INTEGER MGX,LSX,NPO_I,IXIT_I,IOS,UF,JJ,KK,L
      INTEGER NSTORE,NPASS,NCAP,NNEG,NDIV,NMAX,NXIT
      INTEGER*8 NBAD8,NTOT8,NMAT8,NCMP8,NP28,NBB8,NTB8,NCLS8,NNLS8
      REAL(KIND=LDP) RESID,RHO_M,RHO_P,OMG_OUT,RES_OUT
      REAL(KIND=LDP) KR_TOL,KR_BTOL,KR_OMG,KR_STALL,NCAP_C
      REAL(KIND=LDP) CMT_NC,CMT_JMP,DIV_CAP,INR,VP2N,VP2R,VHP,VHT
      REAL(KIND=LDP) OARG,FSEED
      REAL(KIND=LDP), ALLOCATABLE :: FLDP_SEED(:,:),FLDM_POST(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: FLDP_END(:,:),TMP(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CAND_P(:,:),CAND_M(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: CAND_VP(:),CAND_VM(:)
      REAL(KIND=LDP), ALLOCATABLE :: OMGSEQ(:)
      REAL(KIND=LDP) T1,T2,T3,T4
      INTEGER NSEQ,IRF_D,ILS_D,IPO_D,IDIR_D,IPV_D,NOSC_D
      CHARACTER(LEN=12) IMETH_S
      CHARACTER(LEN=256) DECK
      CHARACTER(LEN=32) MODE
      character(len=4096) sval, omega_line
      INTEGER NLOOP,IV,JP2,K,UH,NARG,IARG,UC,NCMG,KNI_C,KCF_C
      INTEGER*8 NDF8
      LOGICAL LCOLD,LSCALE,LV1
      LOGICAL LP2_SET,LP2_VAL,LP2_DECK
      CHARACTER(LEN=4096) ATMP
      REAL(KIND=LDP) DMA,DMI,DTMP
      INTEGER IMA,JMA,IMI,JMI
      EXTERNAL KNOP_RAY_SOLVE_V1
!
      CALL GETARG(1,DECK)
      CALL GETARG(2,MODE)
      CALL GETARG(3,SVAL)
      IF(LEN_TRIM(DECK) .EQ. 0)THEN
        WRITE(*,'(A)')'RESULT usage: deck [mode] [value]'
        STOP
      END IF
      IF(MODE .EQ. '')MODE='base'
      OARG=-1.0_LDP
      IF(SVAL .NE. '')READ(SVAL,*,IOSTAT=IOS)OARG
      NLOOP=0
      LCOLD=.FALSE.; LSCALE=.FALSE.; FSEED=1.0_LDP
!
! p2=<v> (phase 2.30, replay only): scan every argument slot for a
! keyed P2 token before the deck is opened so a malformed override
! fails fast instead of silently falling back to the production
! default. Slot 3 doubles as the mode value; a schedule path or
! seed value must therefore never begin with 'p2'. The override
! itself is applied after the deck restores its recorded flag.
!
      LP2_SET=.FALSE.; LP2_VAL=.FALSE.
      NARG=COMMAND_ARGUMENT_COUNT()
      DO IARG=3,NARG
        CALL GETARG(IARG,ATMP)
        IF(ATMP(1:2) .EQ. 'p2' .OR. ATMP(1:2) .EQ. 'P2')THEN
          IF(ATMP(3:3) .NE. '=')ERROR STOP 'p2: malformed override'
          IF(LP2_SET)ERROR STOP 'p2: duplicate override'
          SELECT CASE(TRIM(ATMP(4:)))
          CASE('0','f','F','false','FALSE','off','OFF','n','N', &
      'no','NO')
            LP2_VAL=.FALSE.
          CASE('1','t','T','true','TRUE','on','ON','y','Y','yes', &
      'YES')
            LP2_VAL=.TRUE.
          CASE DEFAULT
            ERROR STOP 'p2: malformed override'
          END SELECT
          LP2_SET=.TRUE.
        END IF
      END DO
!
      OPEN(NEWUNIT=UF,FILE=TRIM(DECK),FORM='UNFORMATTED', &
      STATUS='OLD',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT deck_open_fail=1'
        STOP
      END IF
      READ(UF)MGX,LSX,NI,NDD,NXD,NCF,KNC,NCPD,NADD,NCPS,IMETH, &
      KMAXP,KCNT,FBITS,KNPD
! Deck format versions are date-styled literals chosen by the
! writer, not runtime stamps. v1 (20261122, live d50dcacd build)
! predates the KN_CHI_EXT/KN_ETA_EXT and JMP_SIB records; v2
! (20261123, current source) carries them under FBITS 9/14.
      LV1=MGX .EQ. 20261122
      IF(MGX .NE. 20261123 .AND. .NOT. LV1)THEN
        WRITE(*,'(A,I10,A)')'RESULT deck_version_unsupported=',MGX, &
      ' reader_supports=20261122,20261123'
        STOP
      END IF
! v1 never wrote FBITS entries above bit 8; clear them so a v1
! deck cannot alias into the v2-only record branches.
      IF(LV1)FBITS=IAND(FBITS,Z'1FF')
      READ(UF)KR_TOL,KR_BTOL,KR_OMG,KR_STALL,NCAP_C,CMT_NC,CMT_JMP, &
      DIV_CAP,INR,VP2N,VP2R,VHP,VHT
      READ(UF)IMETH_S
!
      KN_ND=NDD; KN_NP=1; KN_NCF=NCF
! The replay renumbers the dumped ray as LS=1, so the solver's core
! tests (LS .LE. KN_NC, LS .LE. NC_PNT_SRCE_L) must be re-expressed
! for that index. Passing the deck's KN_NC through made every
! replayed ray a core ray: a non-core ray then took the DIFFUSION
! I_CORE branch with SQRT(R_ND**2-P**2) of a negative number, and
! row NI of the P march was sanitized to 1E-30 instead of IM_NI.
      KN_NC=MERGE(1,0,LSX .LE. KNC)
! NC_PNT_SRCE_L is a solver PARAMETER and cannot be remapped here;
! it only matters for KN_INNER_METH='PNT_SRCE' decks.
      IF(IMETH .EQ. 3 .AND. LSX .GT. NCPS)WRITE(*,'(A)') &
      'RESULT warn_pnt_srce_core_index_not_remapped=1'
      KN_ND_EXT=NXD; ND_ADD=NADD; N_COMP=NCPD
      KN_INNER_METH=IMETH_S
      IF(IMETH .EQ. 1)KN_INNER_METH='DIFFUSION'
      IF(IMETH .EQ. 2)KN_INNER_METH='ZERO_FLUX'
      IF(IMETH .EQ. 3)KN_INNER_METH='PNT_SRCE'
!
      ALLOCATE(R_EXT(NXD),LOG_R_EXT(NXD),KN_R(NDD),KN_P(1), &
      Z(NI,1),R_RAY(NI,1),RAY_PNT(NI,1),NI_RAY(1), &
      J_PNT(NDD,1),COMP_DEPTH_IDX(NDD), &
      GAM(NI,1),GAM_MAG(NI,1),dGAMdR(NI,1), &
      KN_CHI_COEF(NXD,4,NCF),KN_ETA_COEF(NXD,4,NCF), &
      KN_FREQ(NCF),KN_DLN_P(NCF),KN_DLN_N(NCF), &
      KN_IC(NCF),KN_DBB(NCF), &
      I_P_DEPTH_NEXT(NCPD,1,NCF-1), &
      I_M_DEPTH_NEXT(NCPD,1,NCF-1), &
      FG_RAY_SEEDED(1),KN_CS_LS(1),KN_CS_REA(1), &
      FLDP_SEED(NI,NCF),FLDM_POST(NI,NCF),FLDP_END(NI,NCF), &
      KN_RES_PASS(2*(KMAXP+INT(VHP)+64)),KN_HOT_CNT(NCF), &
      KR_DBG_FLDM(NI,NCF),KR_DBG_FLDP(NI,NCF),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(*,'(A)')'RESULT alloc_fail=1'
        STOP
      END IF
!
      READ(UF)R_EXT
      READ(UF)LOG_R_EXT
      READ(UF)KN_R
      READ(UF)KN_P(1)
      READ(UF)Z(1:NI,1)
      READ(UF)R_RAY(1:NI,1)
      READ(UF)RAY_PNT(1:NI,1)
      READ(UF)J_PNT(1:NDD,1)
      READ(UF)COMP_DEPTH_IDX(1:NDD)
      READ(UF)GAM(1:NI,1)
      READ(UF)GAM_MAG(1:NI,1)
      READ(UF)dGAMdR(1:NI,1)
      READ(UF)KN_CHI_COEF
      READ(UF)KN_ETA_COEF
      IF(LV1 .OR. .NOT. BTEST(FBITS,9))THEN
! v1 (or v2 without the EXT flag): no KN_CHI_EXT/KN_ETA_EXT
! records. BUILD_FREQ reads the EXT node values unguarded for
! on-node points, so they are reconstructed here -- not
! fabricated. KN_CHI_COEF is the MON_INT_FUNS_V2 monotone cubic
! built on LOG(EXT): the interpolant passes through the node
! values by construction, so EXP(coef(K,4,L)) is exactly the
! tabulated node value the writer captured (the solver's JPT
! branch only skips the FP round-trip).
        ALLOCATE(KN_CHI_EXT(NXD,NCF),KN_ETA_EXT(NXD,NCF),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(*,'(A)')'RESULT alloc_fail=1'
          STOP
        END IF
        DO JJ=1,NCF
          DO KK=1,NXD
            KN_CHI_EXT(KK,JJ)=EXP(KN_CHI_COEF(KK,4,JJ))
            KN_ETA_EXT(KK,JJ)=EXP(KN_ETA_COEF(KK,4,JJ))
          END DO
        END DO
        WRITE(*,'(A)')'RESULT deck_ext_mode=coef_node'
      ELSE
        ALLOCATE(KN_CHI_EXT(NXD,NCF),KN_ETA_EXT(NXD,NCF),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(*,'(A)')'RESULT alloc_fail=1'
          STOP
        END IF
        READ(UF)KN_CHI_EXT
        READ(UF)KN_ETA_EXT
      END IF
      READ(UF)KN_FREQ
      READ(UF)KN_DLN_P
      READ(UF)KN_DLN_N
      READ(UF)KN_IC
      READ(UF)KN_DBB
      READ(UF)FLDM_POST
      READ(UF)((I_M_DEPTH_NEXT(JJ,1,KK),JJ=1,NCPD),KK=1,NCF-1)
      READ(UF)((I_P_DEPTH_NEXT(JJ,1,KK),JJ=1,NCPD),KK=1,NCF-1)
      IF(LV1 .OR. .NOT. BTEST(FBITS,14))THEN
! v1 decks carry no sibling maps; JMP_SIB only feeds diagnostic
! sibref prints, so zero-fill keeps replay numerics unchanged.
        ALLOCATE(JMP_SIB_M(NCF),JMP_SIB_P(NCF),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(*,'(A)')'RESULT alloc_fail=1'
          STOP
        END IF
        JMP_SIB_M=0.0_LDP; JMP_SIB_P=0.0_LDP
        WRITE(*,'(A)')'RESULT deck_sib_mode=zero'
      ELSE
        ALLOCATE(JMP_SIB_M(NCF),JMP_SIB_P(NCF),STAT=IOS)
        IF(IOS .NE. 0)THEN
          WRITE(*,'(A)')'RESULT alloc_fail=1'
          STOP
        END IF
        READ(UF)JMP_SIB_M
        READ(UF)JMP_SIB_P
      END IF
      READ(UF)FLDP_SEED
      READ(UF,IOSTAT=IOS)RES_OUT,NPO_I,IXIT_I,OMG_OUT
      IF(IOS .EQ. 0)READ(UF,IOSTAT=IOS)FLDP_END
      CLOSE(UF)
!
! Flags: FBITS 0=FG_FREQ_MATCH 1=P2DAMP 2=ORD_SWEEP 3=AITKEN 4=AA
! 5=probe 6=BULK_CONTR 7=ADAPT_OMEGA_V2 8=PLATEAU_DET
!
      FG_FREQ_MATCH=BTEST(FBITS,0)
      FG_ORD_P2DAMP=BTEST(FBITS,1)
      FG_ORD_SWEEP=BTEST(FBITS,2)
      FG_ORD_AITKEN=BTEST(FBITS,3)
      FG_ORD_AA=BTEST(FBITS,4)
      FG_BULK_CONTR=BTEST(FBITS,6)
      FG_ADAPT_OMEGA_V2=BTEST(FBITS,7)
      FG_PLATEAU_DET=BTEST(FBITS,8)
      KN_NOCP_PROBE_LS=-1
      IF(BTEST(FBITS,5))KN_NOCP_PROBE_LS=1
!
      KN_RAY_TOL=KR_TOL
      KN_RAY_BULK_TOL=KR_BTOL
      KN_RAY_OMEGA=KR_OMG
      KN_RAY_STALL=KR_STALL
      NONMONO_INEXT_CAP=NCAP_C
      FG_COMMIT_NC=CMT_NC
      FG_COMMIT_JUMP=CMT_JMP
      FG_DIV_HEAL_CAP=DIV_CAP
      INEXT_RELAX=INR
      FG_ORD_P2N=INT(VP2N); FG_ORD_P2REL=INT(VP2R)
      KN_RAY_HOTPASSES=INT(VHP)
      KN_HOT_TOTAL_PASSES=INT(VHT)
      KN_RAY_MAXPASS=KMAXP
      KN_CONT_NP=KCNT
      NI_RAY(1)=NI
! Match the refresh start: uncommitted certification status and the
! same seeded-flag state the dumped ray carried, so the first-commit
! jump veto behaves as in production.
      FG_SEED_COMMIT=BTEST(FBITS,10)
      FG_RAY_SEEDED(1)=MERGE(1,0,BTEST(FBITS,11))
      KN_Q_RESID=BTEST(FBITS,12)
      KN_CERT_CONV=BTEST(FBITS,13)
      KN_CS_LS=KN_CS_INV; KN_CS_REA=5
      KN_RES_PASS=0.0_LDP
      KN_HOT_CNT=0
      KN_HOT_SPENT=0_8
      KN_DBG_LS=1
      KN_SEED_ARM=-1
      KN_PSWEEP_DIAG=.FALSE.
      FREEZE_I_NEXT=.FALSE.
      FREEZE_I_PREV=.FALSE.
      KNOP_DIAG_NO_REGRID=.TRUE.
      GS_FOREIGN_SEED=.FALSE.
      KN_RAY_PACK_MB=2048.0_LDP
      FG_COMMIT_CONV=.TRUE.
      KN_ADAPT_RELAX=.TRUE.
      LEDG_N_RF=1
!
! p2= apply (phase 2.30): the deck-restored value is recorded for
! the RESULT echo, then the replay-only override is the single
! assignment that can change FG_ORD_P2DAMP in this driver. No
! other module state is touched, so an A/B pair differs only in
! the P2 midpoint commit.
!
      LP2_DECK=FG_ORD_P2DAMP
      IF(LP2_SET)FG_ORD_P2DAMP=LP2_VAL
      WRITE(*,'(A,L1,A,L1,A,L1)')'RESULT p2_deck=',LP2_DECK, &
      ' p2_eff=',FG_ORD_P2DAMP,' p2_ovr=',LP2_SET
!
! Variant overrides.
!
      IF(MODE .EQ. 'om' .AND. OARG .GT. 0.0_LDP)KN_RAY_OMEGA=OARG
      IF(MODE .EQ. 'btol' .AND. OARG .GE. 0.0_LDP)KN_RAY_BULK_TOL=OARG
      IF(MODE .EQ. 'stall' .AND. OARG .GT. 0.0_LDP)KN_RAY_STALL=OARG
      IF(MODE .EQ. 'maxp' .AND. OARG .GT. 0.0_LDP) &
      KN_RAY_MAXPASS=INT(OARG)
      IF(MODE .EQ. 'loop' .OR. MODE .EQ. 'hst') THEN
        NLOOP=MAX(1,INT(OARG))
        KN_RAY_MAXPASS=1
      END IF
      if (mode == 'omf' .or. mode == 'stage') then
        if (.not. fg_ord_sweep .or. ni < 3) error stop 'replay: ordered path required'
        kn_stage_enabled = .true.
        kn_stage_ls = 1
        kn_stage_ray = lsx
        kn_diag_refresh = 1_8
      end if
! omf keeps the production multi-pass lifecycle. The checked schedule
! overrides only the applied omega at ordered pass entry; the controller
! still updates its own residual, oscillation and stopping histories.
! Both channels are required because M supplies the P inner boundary.
! A canonical schedule names its original ray and refresh explicitly.
! No missing entry may use a base-omega fallback.
      IF(MODE .EQ. 'omf')THEN
        call getarg(4, omega_line)
        if (len_trim(omega_line) == 0) error stop 'omf: require refresh argument'
        read(omega_line, *, iostat=ios) irf_d
        if (ios /= 0) error stop 'omf: invalid refresh'
        call kn_replay_load(trim(sval), lsx, irf_d)
        kn_replay_ls = 1
        ledg_n_rf = irf_d
        kn_diag_refresh = int(irf_d, 8)
! KNOMG2 requires declared counts for both channels and IEEE bit checks.
! Legacy KNOMG records cannot establish those semantics on their own.
! Preserve input refresh identity instead of silently selecting a block.
! Applied direction 0 denotes ordered row-wise frequency marching.
! The original deck remains unchanged and is still read as v1/v2.
! Production loop state is retained by the single normal ray call below.
! The schedule is not a request to relax certification or stopping tests.
! A trajectory that terminates early fails the consumption check below.
! Extra requested passes likewise fail at the first unavailable value.
! Channel 1 is P, channel 2 is M; this is distinct from KN_PA indexing.
      END IF
! impose (phase 2.27G/E2): substitute the deck IM_NI boundary row for
! the replay-recomputed one at the M-to-P handoff inside the same
! uninterrupted solve, isolating the downstream P contribution from
! the M reconstruction gap. Stage rows record both values.
      IF(MODE .EQ. 'impose')THEN
        ALLOCATE(KN_IMPOSE_ROW(NCF),STAT=IOS)
        IF(IOS .NE. 0)ERROR STOP 'impose: alloc'
        KN_IMPOSE_ROW(1:NCF)=FLDM_POST(NI,1:NCF)
        KN_IMPOSE_IMNI=.TRUE.
        KN_STAGE_ENABLED=.TRUE.; KN_STAGE_LS=1; KN_STAGE_RAY=LSX
        KN_DIAG_REFRESH=1_8
      END IF
      IF(MODE .EQ. 'seed')THEN
        IF(SVAL .EQ. 'cold')THEN
          LCOLD=.TRUE.
        ELSE
          READ(SVAL,*,IOSTAT=IOS)FSEED
          LSCALE=(IOS .EQ. 0)
        END IF
      END IF
      IF(MODE .EQ. 'verify')THEN
! Seed the P store with the production final iterate and take one
! unrelaxed pass: RES is then the one-pass map defect measured at
! the stalled state itself, not along the seed transient.
        KN_RAY_OMEGA=1.0_LDP
        KN_RAY_MAXPASS=1
        DO K=1,NDD
          IF(COMP_DEPTH_IDX(K) .GE. 1)THEN
            JP2=J_PNT(K,1)
            IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
              DO L=2,NCF
                I_P_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &
     &            FLDP_END(JP2,L)
              END DO
            END IF
          END IF
        END DO
      END IF
      IF(MODE .EQ. 'mapd')THEN
! mapd (phase 2.30): independent raw map-defect evaluation. arg3
! is a KN_REPLAY_END-format candidate (magic, NI, NCF, committed
! P, committed M). Both lagged store columns are seeded from the
! candidate, P2 and relaxation are forced off, and the single
! unrelaxed pass then yields R(X) on both channels inside one
! production solve, so the internal M-to-P boundary hand-off is
! the production one. The candidate file is read-only. p2= is
! rejected here because this mode measures the UNDAMPED map.
        IF(LP2_SET)ERROR STOP 'mapd: p2 override not allowed'
        IF(LEN_TRIM(SVAL) .EQ. 0) &
      ERROR STOP 'mapd: require candidate file'
        ALLOCATE(CAND_P(NI,NCF),CAND_M(NI,NCF),STAT=IOS)
        IF(IOS .NE. 0)ERROR STOP 'mapd: alloc'
        OPEN(NEWUNIT=UC,FILE=TRIM(SVAL),FORM='UNFORMATTED', &
      STATUS='OLD',IOSTAT=IOS)
        IF(IOS .NE. 0)ERROR STOP 'mapd: candidate open'
        READ(UC,IOSTAT=IOS)NCMG,KNI_C,KCF_C
        IF(IOS .NE. 0 .OR. NCMG .NE. 20261124 .OR. &
      KNI_C .NE. NI .OR. KCF_C .NE. NCF) &
      ERROR STOP 'mapd: candidate shape'
        READ(UC,IOSTAT=IOS)CAND_P
        IF(IOS .NE. 0)ERROR STOP 'mapd: candidate P'
        READ(UC,IOSTAT=IOS)CAND_M
        IF(IOS .NE. 0)ERROR STOP 'mapd: candidate M'
        CLOSE(UC)
        KN_RAY_OMEGA=1.0_LDP
        KN_RAY_MAXPASS=1
        FG_ORD_P2DAMP=.FALSE.
        DO K=1,NDD
          IF(COMP_DEPTH_IDX(K) .GE. 1)THEN
            JP2=J_PNT(K,1)
            IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
              DO L=2,NCF
                I_P_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &
      CAND_P(JP2,L)
                I_M_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &
      CAND_M(JP2,L)
              END DO
            END IF
          END IF
        END DO
      END IF
      IF(LCOLD)THEN
        I_P_DEPTH_NEXT=0.0_LDP
        I_M_DEPTH_NEXT=0.0_LDP
      ELSE IF(LSCALE)THEN
        I_P_DEPTH_NEXT=I_P_DEPTH_NEXT*FSEED
        I_M_DEPTH_NEXT=I_M_DEPTH_NEXT*FSEED
      END IF
!
      NSTORE=0; NPASS=0; NCAP=0; NNEG=0; NDIV=0
      NMAX=0; NXIT=0; NBAD8=0_8; NTOT8=0_8; NMAT8=0_8
      NCMP8=0_8; NP28=0_8; NBB8=0_8; NTB8=0_8; NCLS8=0_8; NNLS8=0_8
      RESID=0.0_LDP; RHO_M=0.0_LDP; RHO_P=0.0_LDP
!
      IF(NLOOP .LE. 1)THEN
        CALL KNOP_RAY_SOLVE_V1(1,RESID,NSTORE,NPASS,NCAP,NNEG, &
      NBAD8,NTOT8,NDIV,NMAT8,NCMP8,NMAX,NXIT,NP28,NBB8,NTB8, &
      RHO_M,RHO_P,NCLS8,NNLS8)
        if (mode == 'omf') then
          if (any(kn_replay_used /= kn_replay_count)) &
            call kn_replay_fail(trim(sval), lsx, 0, 'unused scheduled passes')
        end if
        WRITE(*,'(A,ES13.6,A,I3,A,I4,A,F6.3)') &
      'RESULT res=',RESID,' ixit=',NXIT,' npo=',NPASS, &
      ' deck_omega=',OMG_OUT
        WRITE(*,'(A,I4,A,I3)')'RESULT prod_npo=',NPO_I, &
      ' prod_ixit=',IXIT_I
        WRITE(*,'(A,ES13.6)')'RESULT prod_res=',RES_OUT
        WRITE(*,'(A,ES13.6)')'RESULT resid_arg=',RESID
! Committed end state for offline trajectory work (phase 2.24):
! maxp=n runs reproduce production passes 1..n, so per-n dumps
! recover the production-trajectory committed fields.
        OPEN(NEWUNIT=UH,FILE='KN_REPLAY_END.bin', &
      FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=IOS)
        IF(IOS .EQ. 0)THEN
          WRITE(UH)20261124,NI,NCF
          WRITE(UH)KR_DBG_FLDP
          WRITE(UH)KR_DBG_FLDM
          CLOSE(UH)
        END IF
! Replay-vs-production field agreement; the test binary compiles
! -O0 against production -O3 objects, so FP-level drift is expected
! and only a large diff would mean the deck missed live state.
        IF(IOS .EQ. 0)WRITE(*,'(A,ES13.6)') &
      'RESULT replay_vs_prod_maxdiff=', &
      MAXVAL(ABS(KR_DBG_FLDP-FLDP_END))
        IF(MODE .EQ. 'mapd')THEN
! Raw defect D(i,l)=|R(X)-X|/MAX(|X|,1E-4*VMAXF(L)) on the
! candidate-derived column scale (the production pass scale
! recomputes VMAXF from its own pre-pass field, which is exactly
! the seeded candidate). KN_MAPD_RX.bin carries X, R(X) and the
! column scales so the fixed-baseline-mask and significance-
! migration checks live in the analysis layer, not the driver.
          ALLOCATE(CAND_VP(NCF),CAND_VM(NCF),STAT=IOS)
          IF(IOS .NE. 0)ERROR STOP 'mapd: scale alloc'
          DO L=1,NCF
            CAND_VP(L)=MAX(MAXVAL(CAND_P(1:NI,L)),1.0E-25_LDP)
            CAND_VM(L)=MAX(MAXVAL(CAND_M(1:NI,L)),1.0E-25_LDP)
          END DO
          DMA=0.0_LDP; DMI=0.0_LDP
          IMA=0; JMA=0; IMI=0; JMI=0; NDF8=0_8
          DO L=1,NCF
            DO K=1,NI
              DTMP=ABS(KR_DBG_FLDP(K,L)-CAND_P(K,L))/ &
      MAX(ABS(CAND_P(K,L)),1.0E-4_LDP*CAND_VP(L))
              IF(DTMP .GT. DMA)THEN
                DMA=DTMP; IMA=K; JMA=L
              END IF
              IF(DTMP .GT. KN_RAY_TOL)NDF8=NDF8+1_8
              DTMP=ABS(KR_DBG_FLDM(K,L)-CAND_M(K,L))/ &
      MAX(ABS(CAND_M(K,L)),1.0E-4_LDP*CAND_VM(L))
              IF(DTMP .GT. DMI)THEN
                DMI=DTMP; IMI=K; JMI=L
              END IF
              IF(DTMP .GT. KN_RAY_TOL)NDF8=NDF8+1_8
            END DO
          END DO
          WRITE(*,'(A,ES13.6,A,I6,A,I6)') &
      'RESULT mapd_dmax_p=',DMA,' at_i=',IMA,' l=',JMA
          WRITE(*,'(A,ES13.6,A,I6,A,I6)') &
      'RESULT mapd_dmax_m=',DMI,' at_i=',IMI,' l=',JMI
          WRITE(*,'(A,I10)')'RESULT mapd_nbad=',NDF8
          OPEN(NEWUNIT=UC,FILE='KN_MAPD_RX.bin', &
      FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=IOS)
          IF(IOS .EQ. 0)THEN
            WRITE(UC)20261201,NI,NCF
            WRITE(UC)CAND_P
            WRITE(UC)KR_DBG_FLDP
            WRITE(UC)CAND_M
            WRITE(UC)KR_DBG_FLDM
            WRITE(UC)CAND_VP
            WRITE(UC)CAND_VM
            CLOSE(UC)
          ELSE
            ERROR STOP 'mapd: KN_MAPD_RX.bin open failed'
          END IF
        END IF
      ELSE
        IF(MODE .EQ. 'hst')THEN
! Signed per-cell pass history (phase 2.24): KN_RAY_MAXPASS=1 makes
! each solve call exactly one march pass, so KR_DBG_FLDP after the
! call is the committed iterate x(n). The seed field is written as
! x(0); D1/D2 tests, denominator floors and the DS-corner lag can
! all be evaluated offline without solver instrumentation.
          OPEN(NEWUNIT=UH,FILE='KN_REPLAY_HST.bin', &
      FORM='UNFORMATTED',STATUS='REPLACE',IOSTAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(*,'(A)')'RESULT hst_open_fail=1'
            STOP
          END IF
          WRITE(UH)20261124,NI,NCF
          WRITE(UH)0
          WRITE(UH)FLDP_SEED
        END IF
        DO IV=1,NLOOP
          CALL KNOP_RAY_SOLVE_V1(1,RESID,NSTORE,NPASS,NCAP,NNEG, &
      NBAD8,NTOT8,NDIV,NMAT8,NCMP8,NMAX,NXIT,NP28,NBB8,NTB8, &
      RHO_M,RHO_P,NCLS8,NNLS8)
          WRITE(*,'(A,I3,A,ES13.6,A,I3)') &
      'RESULT pass=',IV,' res=',RHO_P,' ixit=',NXIT
          IF(MODE .EQ. 'hst')THEN
            WRITE(UH)IV
            WRITE(UH)KR_DBG_FLDP
          END IF
! Write the post-pass P field back into the store so the next
! one-pass call seeds from it -- the caller-side equivalent of
! the production pass loop, minus the in-sweep omega drift.
          DO K=1,NDD
            IF(COMP_DEPTH_IDX(K) .GE. 1)THEN
              JP2=J_PNT(K,1)
              IF(JP2 .GE. 1 .AND. JP2 .LE. NI)THEN
                DO L=2,NCF
                  I_P_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &
      KR_DBG_FLDP(JP2,L)
                END DO
              END IF
            END IF
          END DO
          KN_HOT_CNT=0
          KN_RES_PASS=0.0_LDP
        END DO
        IF(MODE .EQ. 'hst')CLOSE(UH)
      END IF
      STOP
      END
