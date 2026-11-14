!
! P8.2 (diagnostic only): electron-scattering redistribution
! convergence instrumentation for CMF_FLUX_SUB_V5.
!
! Tracks ES_J_CONV file generations across the incoherent-ES sweep
! loop and writes zero-safe per-depth and global residuals between
! consecutive convolved fields to ES_RESID_DIAG.  Every call site in
! cmf_flux_sub_v5.f is gated on [ES_RESID_DIAG] (default F); with the
! key off this module is never entered and the production output is
! unchanged.  The module holds the previous-generation buffer and the
! generation counter so the call sites only forward sweep/phase
! labels - no production array is touched.
!
! Generation bookkeeping: a pre-sweep-1 ES_J_CONV (the OLD_J warm
! start) is generation 0; each end-of-sweep COMP_J_CONV_V2 write of
! ES_J_CONV bumps the counter.  The render folds ESEC*conv(J_last)
! built by the in-place convolve of RJ_CMF_ST - by construction the
! same generation as the last file write.  ES_RESID_RENDER verifies
! the render field against the stored last generation rather than
! assuming it.  Under ACCURATE (extended depth grid) the file field
! is on NDEXT while the render field is on ND: the residual tracking
! still works (file vs file, same grid across generations) but the
! elementwise render check is only run when the grids match.
!
      MODULE ES_RESID_DIAG_MOD
      USE SET_KIND_MODULE
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE
!
      REAL(KIND=LDP), ALLOCATABLE :: RJ_PREV(:,:)
      INTEGER :: PREV_GEN=-1          !generation label of stored field
      INTEGER :: GEN_COUNT=0          !completed ES_J_CONV generations
      INTEGER :: N_SW_UPD=0           !UPDATE calls (sweep writes seen)
      INTEGER :: N_INVAL_SW=0         !sweeps marked invalid (non-finite)
      INTEGER :: N_NEG_TOT=0          !cumulative negative-node count
      INTEGER :: N_NONFIN_TOT=0       !cumulative non-finite-node count
      INTEGER :: LU_DIAG=-1
      LOGICAL :: DIAG_OPEN=.FALSE.
      LOGICAL :: HAS_PREV=.FALSE.
      LOGICAL :: DEPTH_ROWS=.TRUE.    ![ES_DIAG_DEPTHS] per-depth rows
      LOGICAL :: RENDER_DONE=.FALSE.
      LOGICAL :: RENDER_MATCH=.FALSE.
      LOGICAL :: RENDER_GRID_OK=.FALSE.
      INTEGER :: RENDER_GEN=-1        !generation the render consumed
      REAL(KIND=LDP) :: RENDER_DIFF=0.0_LDP
      CHARACTER(LEN=16) :: PATHWAY_LAB='UNKNOWN'
      REAL(KIND=LDP) :: AFLOOR_USE=1.0E-12_LDP
!
      CONTAINS
!
!
!
      SUBROUTINE ES_RESID_RESET()
!
! Clear all tracked state.  Called once per CMF_FLUX invocation so a
! stale buffer from a previous run or phase cannot be compared
! against the new field.
!
      IF(ALLOCATED(RJ_PREV))DEALLOCATE(RJ_PREV)
      PREV_GEN=-1
      GEN_COUNT=0
      N_SW_UPD=0
      N_INVAL_SW=0
      N_NEG_TOT=0
      N_NONFIN_TOT=0
      HAS_PREV=.FALSE.
      RENDER_DONE=.FALSE.
      RENDER_MATCH=.FALSE.
      RENDER_GRID_OK=.FALSE.
      RENDER_GEN=-1
      RENDER_DIFF=0.0_LDP
      RETURN
      END SUBROUTINE ES_RESID_RESET
!
!
!
      SUBROUTINE ES_RESID_OPEN_OUT(LUER)
      INTEGER LUER,IOS
!
      IF(DIAG_OPEN)RETURN
      CALL GET_LU(LU_DIAG,'ES_RESID_DIAG file')
      OPEN(UNIT=LU_DIAG,FILE='ES_RESID_DIAG',STATUS='REPLACE',
     1     FORM='FORMATTED',IOSTAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(LUER,*)'ES_RESID_DIAG: open failed, IOS=',IOS
        LU_DIAG=LUER
      END IF
      DIAG_OPEN=.TRUE.
      RETURN
      END SUBROUTINE ES_RESID_OPEN_OUT
!
!
!
      SUBROUTINE ES_RESID_READ_ESJ(FNAME,FIELD,NUV,ND_OUT,NCF_OUT,
     1                             IOS,LUER)
!
! Read one whole convolved field from an EDDFACTOR-layout direct
! access file.  The <FNAME>_INFO sidecar supplies ND and the system
! dependent record length; the record-3 header then carries
! INIT_REC/NCF/ND, so the reader follows whatever depth grid the
! writer used (ND or NDEXT) instead of assuming the core grid.
!
      CHARACTER(LEN=*) FNAME
      REAL(KIND=LDP), ALLOCATABLE :: FIELD(:,:)
      REAL(KIND=LDP), ALLOCATABLE :: NUV(:)
      INTEGER ND_OUT,NCF_OUT,IOS,LUER
!
      INTEGER INIT_REC,NCF_RD,ND_RD,IRECL,LU_RD,ML,K
      CHARACTER(LEN=20) FILE_DATE
!
      ND_OUT=0
      NCF_OUT=0
      CALL GET_LU(LU_RD,'ES_RESID read')
      CALL READ_DIRECT_INFO_V3(ND_RD,IRECL,FILE_DATE,FNAME,
     1                         LU_RD,IOS)
      IF(IOS .NE. 0)RETURN
!
      OPEN(UNIT=LU_RD,FILE=FNAME,FORM='UNFORMATTED',ACCESS='DIRECT',
     1     STATUS='OLD',ACTION='READ',RECL=IRECL,IOSTAT=IOS)
      IF(IOS .NE. 0)RETURN
!
! Header record: INIT_REC, NCF, ND (written by COMP_J_CONV_V2 from
! the input file's own header).
!
      READ(LU_RD,REC=3,IOSTAT=IOS)INIT_REC,NCF_RD,ND_OUT
      IF(IOS .NE. 0)THEN
        CLOSE(LU_RD)
        ND_OUT=0
        RETURN
      END IF
      IF(ND_OUT .NE. ND_RD)THEN
        CLOSE(LU_RD)
        WRITE(LUER,*)'ES_RESID_DIAG: ',TRIM(FNAME),
     1    ' INFO/header ND mismatch ',ND_RD,' vs ',ND_OUT
        ND_OUT=0
        RETURN
      END IF
!
      IF(ALLOCATED(FIELD))DEALLOCATE(FIELD)
      IF(ALLOCATED(NUV))DEALLOCATE(NUV)
      ALLOCATE(FIELD(ND_RD,NCF_RD),STAT=IOS)
      IF(IOS .EQ. 0)ALLOCATE(NUV(NCF_RD),STAT=IOS)
      IF(IOS .NE. 0)THEN
        CLOSE(LU_RD)
        WRITE(LUER,*)'ES_RESID_DIAG: alloc failed in read'
        ND_OUT=0
        NCF_OUT=0
        RETURN
      END IF
!
      DO ML=1,NCF_RD
        READ(LU_RD,REC=INIT_REC+ML-1,IOSTAT=IOS)
     1       (FIELD(K,ML),K=1,ND_RD),NUV(ML)
        IF(IOS .NE. 0)THEN
          CLOSE(LU_RD)
          WRITE(LUER,*)'ES_RESID_DIAG: short file ',TRIM(FNAME),
     1                 ' at record ',INIT_REC+ML-1
          ND_OUT=0
          NCF_OUT=0
          RETURN
        END IF
      END DO
      CLOSE(LU_RD)
      ND_OUT=ND_RD
      NCF_OUT=NCF_RD
      IOS=0
      RETURN
      END SUBROUTINE ES_RESID_READ_ESJ
!
!
!
      SUBROUTINE ES_RESID_COMPARE(RJ_NEW,RJ_OLD,ND,NCF,AFLOOR,
     1    R_K,D_K,S_K,RREL_K,NACT_K,NMASK_K,NNONF_K,NNEG_K,
     2    R_GLOB,K_GLOB,RREL_GLOB,D_GLOB,TOT_ACT,TOT_INVAL,TOT_NEG,
     3    BAD_K,BAD_M,BAD_WHICH,VALID)
!
! Zero-safe residual between two consecutive RJ_ES generations on
! the same grid (P8.2 spec sec 4).  This is the routine the unit
! tests drive directly: it takes both fields as arguments and holds
! no state, so manufactured cases exercise the production code path.
!
! Per depth K:
!   S_K   = max_M max(|NEW|,|OLD|)           (finite nodes only)
!   D_K   = max_M |NEW-OLD|
!   R_K   = D_K/S_K; S_K=0 & D_K=0 -> 0; S_K=0 & D_K!=0 -> 1
!   active mask: max(|NEW|,|OLD|) > AFLOOR*S_K
!   RREL_K = max |d|/max(|NEW|,|OLD|) over active nodes
!
! Non-finite nodes are counted, not skipped silently: VALID=F marks
! the whole comparison invalid so no false convergence is reported.
! Negative values are counted and kept - no clamping.
!
      INTEGER ND,NCF
      REAL(KIND=LDP) RJ_NEW(ND,NCF),RJ_OLD(ND,NCF)
      REAL(KIND=LDP) AFLOOR
      REAL(KIND=LDP) R_K(ND),D_K(ND),S_K(ND),RREL_K(ND)
      INTEGER NACT_K(ND),NMASK_K(ND),NNONF_K(ND),NNEG_K(ND)
      REAL(KIND=LDP) R_GLOB,RREL_GLOB,D_GLOB
      INTEGER K_GLOB,TOT_ACT,TOT_INVAL,TOT_NEG
      INTEGER BAD_K,BAD_M,BAD_WHICH
      LOGICAL VALID
!
      INTEGER K,M
      REAL(KIND=LDP) XN,XO,AMAX,AD,REL,SC
!
      DO K=1,ND
        S_K(K)=0.0_LDP
        D_K(K)=0.0_LDP
        NNONF_K(K)=0
        NNEG_K(K)=0
        DO M=1,NCF
          XN=RJ_NEW(K,M)
          XO=RJ_OLD(K,M)
          IF(.NOT. IEEE_IS_FINITE(XN) .OR.
     1       .NOT. IEEE_IS_FINITE(XO))THEN
            NNONF_K(K)=NNONF_K(K)+1
            IF(BAD_K .EQ. 0)THEN
              BAD_K=K
              BAD_M=M
              IF(.NOT. IEEE_IS_FINITE(XN))THEN
                BAD_WHICH=1           !new/current generation
              ELSE
                BAD_WHICH=-1          !old/previous generation
              END IF
            END IF
            CYCLE
          END IF
          AMAX=MAX(ABS(XN),ABS(XO))
          AD=ABS(XN-XO)
          IF(AMAX .GT. S_K(K))S_K(K)=AMAX
          IF(AD .GT. D_K(K))D_K(K)=AD
          IF(XN .LT. 0.0_LDP .OR. XO .LT. 0.0_LDP)
     1       NNEG_K(K)=NNEG_K(K)+1
        END DO
        IF(S_K(K) .GT. 0.0_LDP)THEN
          R_K(K)=D_K(K)/S_K(K)
        ELSE IF(D_K(K) .EQ. 0.0_LDP)THEN
          R_K(K)=0.0_LDP
        ELSE
          R_K(K)=1.0_LDP
        END IF
      END DO
!
! Second pass for the active-frequency mask - it needs the depth
! scale S_K computed over the whole row first.
!
      DO K=1,ND
        NACT_K(K)=0
        NMASK_K(K)=0
        RREL_K(K)=0.0_LDP
        SC=S_K(K)
        DO M=1,NCF
          XN=RJ_NEW(K,M)
          XO=RJ_OLD(K,M)
          IF(.NOT. IEEE_IS_FINITE(XN) .OR.
     1       .NOT. IEEE_IS_FINITE(XO))CYCLE
          AMAX=MAX(ABS(XN),ABS(XO))
          IF(AMAX .GT. AFLOOR*SC)THEN
            NACT_K(K)=NACT_K(K)+1
            IF(AMAX .GT. 0.0_LDP)THEN
              REL=ABS(XN-XO)/AMAX
              IF(REL .GT. RREL_K(K))RREL_K(K)=REL
            END IF
          ELSE
            NMASK_K(K)=NMASK_K(K)+1
          END IF
        END DO
      END DO
!
      R_GLOB=0.0_LDP
      K_GLOB=0
      RREL_GLOB=0.0_LDP
      D_GLOB=0.0_LDP
      TOT_ACT=0
      TOT_INVAL=0
      TOT_NEG=0
      DO K=1,ND
        TOT_ACT=TOT_ACT+NACT_K(K)
        TOT_INVAL=TOT_INVAL+NNONF_K(K)
        TOT_NEG=TOT_NEG+NNEG_K(K)
        IF(D_K(K) .GT. D_GLOB)D_GLOB=D_K(K)
        IF(RREL_K(K) .GT. RREL_GLOB)RREL_GLOB=RREL_K(K)
        IF(R_K(K) .GT. R_GLOB)THEN
          R_GLOB=R_K(K)
          K_GLOB=K
        END IF
      END DO
      VALID=(TOT_INVAL .EQ. 0)
      RETURN
      END SUBROUTINE ES_RESID_COMPARE
!
!
!
      SUBROUTINE ES_RESID_SEED(ND,NCF,KNOP,AFLOOR,DEPTH_RWS,
     1                         NUM_ES,LUER)
!
! Phase-1 entry: snapshot any pre-existing ES_J_CONV (OLD_J warm
! start) as generation 0 so the first UPDATE measures the update
! from the seed instead of skipping the comparison.
!
      INTEGER ND,NCF,NUM_ES,LUER
      LOGICAL KNOP,DEPTH_RWS
      REAL(KIND=LDP) AFLOOR
!
      REAL(KIND=LDP), ALLOCATABLE :: FLD(:,:),NUV(:)
      INTEGER ND_F,NCF_F,IOS
      LOGICAL EX
!
      CALL ES_RESID_RESET()
      CALL ES_RESID_OPEN_OUT(LUER)
      AFLOOR_USE=AFLOOR
      DEPTH_ROWS=DEPTH_RWS
      IF(KNOP)THEN
        PATHWAY_LAB='KNOP'
      ELSE
        PATHWAY_LAB='SOB'
      END IF
!
      WRITE(LU_DIAG,'(A)')'ES_RESID_DIAG v1'
      WRITE(LU_DIAG,'(A,A,A,I0,A,I0,A,I0,A,ES10.3)')
     1  'HEADER PATHWAY=',TRIM(PATHWAY_LAB),
     1  ' ND=',ND,' NCF=',NCF,' CONFIGURED_ES_SWEEPS=',NUM_ES,
     2  ' ACTIVE_FLOOR=',AFLOOR_USE
!
      INQUIRE(FILE='ES_J_CONV',EXIST=EX)
      IF(EX)THEN
        CALL ES_RESID_READ_ESJ('ES_J_CONV',FLD,NUV,ND_F,NCF_F,
     1                         IOS,LUER)
        IF(IOS .EQ. 0)THEN
          CALL MOVE_ALLOC(FLD,RJ_PREV)
          HAS_PREV=.TRUE.
          PREV_GEN=0
          WRITE(LU_DIAG,'(A,I0,A,I0,A)')
     1      'SEED GEN=0 FIELD_ND=',ND_F,' FIELD_NCF=',NCF_F,
     2      ' SOURCE=ES_J_CONV (warm-start field)'
        ELSE
          WRITE(LU_DIAG,'(A)')
     1      'SEED GEN=NONE - ES_J_CONV present but unreadable'
        END IF
      ELSE
        WRITE(LU_DIAG,'(A)')
     1    'SEED GEN=NONE - no ES_J_CONV (cold start, sweep 1 coherent)'
      END IF
      RETURN
      END SUBROUTINE ES_RESID_SEED
!
!
!
      SUBROUTINE ES_RESID_UPDATE(ND,NCF,SWEEP,NUM_ES,LUER)
!
! Called after each end-of-sweep COMP_J_CONV_V2 write of ES_J_CONV:
! reads the fresh generation, compares it against the stored
! previous generation, and writes per-depth + global rows.
!
      INTEGER ND,NCF,SWEEP,NUM_ES,LUER
!
      REAL(KIND=LDP), ALLOCATABLE :: RJ_NEW(:,:),NUV(:)
      REAL(KIND=LDP), ALLOCATABLE :: R_K(:),D_K(:),S_K(:),RREL_K(:)
      INTEGER, ALLOCATABLE :: NACT_K(:),NMASK_K(:),NNONF_K(:),NNEG_K(:)
      INTEGER ND_F,NCF_F,IOS,K
      REAL(KIND=LDP) R_GLOB,RREL_GLOB,D_GLOB
      INTEGER K_GLOB,TOT_ACT,TOT_INVAL,TOT_NEG
      INTEGER BAD_K,BAD_M,BAD_WHICH
      LOGICAL VALID
!
      GEN_COUNT=GEN_COUNT+1
      N_SW_UPD=SWEEP
!
      CALL ES_RESID_READ_ESJ('ES_J_CONV',RJ_NEW,NUV,ND_F,NCF_F,
     1                       IOS,LUER)
      IF(IOS .NE. 0)THEN
        WRITE(LU_DIAG,'(A,I0,A,I0,A)')
     1    'SWEEP ES_SWEEP=',SWEEP,' REFRESH_GENERATION=',GEN_COUNT,
     2    ' VALID=F READ_FAIL'
        N_INVAL_SW=N_INVAL_SW+1
        RETURN
      END IF
!
      IF(.NOT. HAS_PREV)THEN
        WRITE(LU_DIAG,'(A,I0,A,I0,A)')
     1    'GLOBAL ES_SWEEP=',SWEEP,' REFRESH_GENERATION=',GEN_COUNT,
     2    ' FIRST_GENERATION - no previous field, residual skipped'
      ELSE IF(ND_F .NE. SIZE(RJ_PREV,1) .OR.
     1        NCF_F .NE. SIZE(RJ_PREV,2))THEN
        WRITE(LU_DIAG,'(A,I0,A,I0,A)')
     1    'SWEEP ES_SWEEP=',SWEEP,' REFRESH_GENERATION=',GEN_COUNT,
     2    ' VALID=F GRID_CHANGE'
        N_INVAL_SW=N_INVAL_SW+1
      ELSE
        ALLOCATE(R_K(ND_F),D_K(ND_F),S_K(ND_F),RREL_K(ND_F),
     1           NACT_K(ND_F),NMASK_K(ND_F),NNONF_K(ND_F),NNEG_K(ND_F))
        BAD_K=0
        BAD_M=0
        BAD_WHICH=0
        CALL ES_RESID_COMPARE(RJ_NEW,RJ_PREV,ND_F,NCF_F,AFLOOR_USE,
     1    R_K,D_K,S_K,RREL_K,NACT_K,NMASK_K,NNONF_K,NNEG_K,
     2    R_GLOB,K_GLOB,RREL_GLOB,D_GLOB,TOT_ACT,TOT_INVAL,TOT_NEG,
     3    BAD_K,BAD_M,BAD_WHICH,VALID)
!
        N_NONFIN_TOT=N_NONFIN_TOT+TOT_INVAL
        N_NEG_TOT=N_NEG_TOT+TOT_NEG
        IF(.NOT. VALID)N_INVAL_SW=N_INVAL_SW+1
!
        IF(DEPTH_ROWS)THEN
          DO K=1,ND_F
            WRITE(LU_DIAG,'(A,I0,A,I0,A,ES12.5,A,ES12.5,A,ES12.5,'//
     1           'A,ES12.5,A,I0,A,I0,A,I0,A,I0)')
     2        'DEPTH SWEEP=',SWEEP,' K=',K,
     3        ' R=',R_K(K),' RREL=',RREL_K(K),' D=',D_K(K),
     4        ' S=',S_K(K),' NACT=',NACT_K(K),' NMASK=',NMASK_K(K),
     5        ' NNONFIN=',NNONF_K(K),' NNEG=',NNEG_K(K)
          END DO
        END IF
!
        IF(VALID)THEN
          WRITE(LU_DIAG,'(A,I0,A,I0,A,ES12.5,A,I0,A,ES12.5,A,ES12.5,'//
     1         'A,I0,A,I0,A,I0)')
     2      'GLOBAL ES_SWEEP=',SWEEP,' REFRESH_GENERATION=',GEN_COUNT,
     3      ' GLOBAL_MAX_R=',R_GLOB,' DEPTH=',K_GLOB,
     4      ' GLOBAL_MAX_REL=',RREL_GLOB,' GLOBAL_MAX_ABS=',D_GLOB,
     5      ' TOTAL_ACTIVE=',TOT_ACT,' TOTAL_INVALID=',TOT_INVAL,
     6      ' TOTAL_NEG=',TOT_NEG
        ELSE
          WRITE(LU_DIAG,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')
     1      'GLOBAL ES_SWEEP=',SWEEP,' REFRESH_GENERATION=',GEN_COUNT,
     2      ' VALID=F NONFINITE_AT K=',BAD_K,' M=',BAD_M,
     3      ' BAD_GENERATION=',BAD_WHICH,' TOTAL_INVALID=',TOT_INVAL
        END IF
        DEALLOCATE(R_K,D_K,S_K,RREL_K,
     1             NACT_K,NMASK_K,NNONF_K,NNEG_K)
      END IF
!
      IF(ALLOCATED(RJ_PREV))DEALLOCATE(RJ_PREV)
      CALL MOVE_ALLOC(RJ_NEW,RJ_PREV)
      HAS_PREV=.TRUE.
      PREV_GEN=GEN_COUNT
      RETURN
      END SUBROUTINE ES_RESID_UPDATE
!
!
!
      SUBROUTINE ES_RESID_RENDER(RJ,ND,NCF,LUER)
!
! Called after the in-place convolve of RJ_CMF_ST and before the
! ESEC fold: verifies the render field is the latest completed
! generation.  The render field is conv(J_last) on the core grid;
! the stored generation is conv(J_last) on the file grid - the
! elementwise check only runs when the grids coincide (always true
! without ACCURATE, which is the only mode the non-mono renders
! support).  A mismatch means the render consumed a stale or wrong
! field - reported, not asserted away.
!
      INTEGER ND,NCF,LUER
      REAL(KIND=LDP) RJ(ND,NCF)
!
      REAL(KIND=LDP) DMAX
!
      RENDER_DONE=.TRUE.
      RENDER_GEN=GEN_COUNT
!
      IF(.NOT. HAS_PREV)THEN
        RENDER_GRID_OK=.FALSE.
        RENDER_MATCH=.FALSE.
        WRITE(LU_DIAG,'(A,I0,A)')
     1    'RENDER GEN_LAST=',GEN_COUNT,
     2    ' FIELD_MATCH=F NO_STORED_GENERATION'
        RETURN
      END IF
!
      IF(SIZE(RJ_PREV,1) .NE. ND .OR. SIZE(RJ_PREV,2) .NE. NCF)THEN
        RENDER_GRID_OK=.FALSE.
        RENDER_MATCH=.FALSE.
        WRITE(LU_DIAG,'(A,I0,A,I0,A,I0,A)')
     1    'RENDER GEN_LAST=',GEN_COUNT,
     2    ' GRID_MISMATCH file_nd=',SIZE(RJ_PREV,1),
     3    ' render_nd=',ND,
     4    ' - elementwise check skipped'
        RETURN
      END IF
!
      RENDER_GRID_OK=.TRUE.
      DMAX=MAXVAL(ABS(RJ-RJ_PREV))
      RENDER_DIFF=DMAX
      RENDER_MATCH=(DMAX .EQ. 0.0_LDP)
      WRITE(LU_DIAG,'(A,I0,A,L1,A,ES12.5)')
     1  'RENDER GEN_LAST=',GEN_COUNT,' FIELD_MATCH=',RENDER_MATCH,
     2  ' MAX_ABS_DIFF=',DMAX
      RETURN
      END SUBROUTINE ES_RESID_RENDER
!
!
!
      SUBROUTINE ES_RESID_FINAL(NUM_ES,LUER)
!
! Summary block with the generation audit the spec asks for:
! configured vs completed sweeps, final refresh generation, and the
! generation the render actually consumed.
!
      INTEGER NUM_ES,LUER
!
      WRITE(LU_DIAG,'(A,I0,A,I0,A,I0,A,I0,A,L1,A,I0,A,I0,A,I0)')
     1  'FINAL CONFIGURED_ES_SWEEPS=',NUM_ES,
     2  ' COMPLETED_ES_SWEEPS=',N_SW_UPD,
     3  ' FINAL_REFRESH_GENERATION=',GEN_COUNT,
     4  ' FINAL_RENDER_INPUT_GENERATION=',RENDER_GEN,
     5  ' RENDER_MATCH=',RENDER_MATCH,
     6  ' INVALID_SWEEPS=',N_INVAL_SW,
     7  ' TOTAL_NONFINITE=',N_NONFIN_TOT,
     8  ' TOTAL_NEGATIVE=',N_NEG_TOT
      IF(.NOT. RENDER_DONE)THEN
        WRITE(LU_DIAG,'(A)')
     1    'FINAL NOTE - render check never ran (coherent or gated)'
      END IF
!
      IF(ALLOCATED(RJ_PREV))DEALLOCATE(RJ_PREV)
      HAS_PREV=.FALSE.
      IF(DIAG_OPEN .AND. LU_DIAG .NE. LUER)THEN
        CLOSE(LU_DIAG)
        DIAG_OPEN=.FALSE.
      END IF
      RETURN
      END SUBROUTINE ES_RESID_FINAL
!
!
!
      END MODULE ES_RESID_DIAG_MOD
