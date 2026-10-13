!
! Data module for MOM_J_CMF_NONMONO_V1.
! Based on MOD_MOM_J_V11 with non-monotonic extensions.
!
! Non-monotonic moment solver using the Knop et al. (2009, A&A 501, 813)
! generalized-opacity split:
!
!   chi_hat = chi + |GAM|/dLOG_NU   (always positive generalized opacity)
!   Q      = (GAM/dLOG_NU) / chi_hat (signed coupling coefficient)
!   S'     = SOURCE + |Q| * J_upwind (modified source; the implicit
!            diagonal PSI = |Q|*DTAUONQ already carries the sign split,
!            so the explicit upwind term enters with positive |coeff|
!            in both expanding and compressing regions)
!
! where GAM = (V/R) * c^{-1} * (1 + SIGMA_TRUE * K_ON_J) is the
! angle-integrated comoving-frame frequency gradient.
!
! Sign-dependent upwind (Knop Paper II, Gauss-Seidel):
!   GAM >= 0 (expanding):  upwind = J_PREV (previous/bluer frequency)
!   GAM <  0 (compressing): upwind = J_NEXT (next/redder frequency)
!
! J_NEXT is stored from the previous great iteration and converges
! over CMFGEN great iterations, mirroring the I_NEXT approach in
! FG_J_CMF_NONMONO_V1.
!
! The standard MOM_J_CMF_V11 remains unchanged for monotonic velocity
! fields.  This routine is only called when DO_NONMONO_KNOP=.TRUE.
!
! Ported to dev7 19-Sep-2026.  The Jacobian variant of the earlier
! non-monotonic development line (VAR_MOM_J_CMF_NONMONO_V1) and its
! J_NEXT_MD / STORE_J_NEXT_MD model-grid store are not carried over:
! that variant linearized the moment equation with the clamped SIGMA,
! which is inconsistent with the non-monotonic forward solve and
! produced an all-NaN BA matrix at the first BA step.  DO_VAR_CONT
! instead uses a fixed-J linearization (as under USE_LAM_ES) plus the
! diagonal dJ/dCHI, dJ/dETA self-coupling accumulated by
! FG_J_CMF_NONMONO_V1's own approximate-Lambda operator
! (GET_FG_ALO_DIAG), so the linearization matches the forward solve.
! GET_MOM_GS_MATURITY supports the GS warm-up gate in SOLVE_FOR_POPS.
!
! References:
!   Knop, Hauschildt & Baron 2009, A&A 501, 813 (Paper II)
!   Knop, Hauschildt & Baron 2009, A&A 496, 295 (Paper I)
!   Baron & Hauschildt 2004, A&A 427, 987 (wavelength boundary-value
!   problem and local upwind discretisation; Paper II's "Paper I")
!   Ceulemans et al. 2025, JQSRT (numerically stable CMF solver)
!
! edited by Sean Brennan (15-Sep-2026)
!
	MODULE MOD_MOM_J_NONMONO_V1
	USE SET_KIND_MODULE
	IMPLICIT NONE
!
! To be dimensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL(KIND=LDP), ALLOCATABLE :: LOG_R_SM_NM(:)
!
! V_SM_SV: velocity fingerprint at the last grid build.  The expanded
! sub-grid resolution (INT(|V_SM(I)-V_SM(I+1)|/T2)) and all of
! V_NM/SIGMA_NM/SIGMA_TRUE_NM are computed only in the build block,
! so a V change with R unchanged previously left the whole coupling
! (GAM, W, EPS, PSI, Q) computed from the stale velocity field.
! Compared each INIT to trigger NEW_R_GRID.
!
	REAL(KIND=LDP), ALLOCATABLE :: V_SM_SV(:)
!
! Dimensioned ND_SM,4
!
	REAL(KIND=LDP), ALLOCATABLE :: V_COEF_NM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_COEF_NM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_TRUE_COEF_NM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_COEF_NM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_COEF_NM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_COEF_NM(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
!
	REAL(KIND=LDP), ALLOCATABLE :: R_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: V_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SIGMA_TRUE_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: ESEC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: CHI_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_NM(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: RSQJNU_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: RSQHNU_NM(:)
!
	REAL(KIND=LDP), ALLOCATABLE :: RSQJNU_PREV_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: RSQHNU_PREV_NM(:)
!
! J_NEXT storage for Gauss-Seidel convergence (Knop Paper II).
! J_NEXT_FRQ(ND_NM,NCF) stores R^2*J at each frequency, written in-place
! at the end of each frequency solve.  For compressing zones (GAM < 0)
! the upwind field is the NEXT (redder) frequency, so the solve for
! frequency k reads slot k+1, which still holds the previous great
! iteration's value (identical read-before-write pattern to the
! I_M_DEPTH_NEXT store in FG_J_CMF_NONMONO_V1).
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: J_NEXT_FRQ(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: H_NEXT_FRQ(:,:)
	LOGICAL, SAVE :: J_NEXT_VALID=.FALSE.
!
! H_NEXT_FRQ(ND_NM,NCF): lagged R^2*H store, the moment-equation
! partner of J_NEXT_FRQ.  At a compressing midpoint the explicit
! W/(1+W)*H and EPS*J terms of the H equation must read the redder
! frequency (slot k+1), not RSQ{H,J}NU_PREV_NM (the bluer k-1
! field): the lagged map's fixed point must satisfy the correct
! upwind equations, otherwise a wrong-signed coupling is injected
! into rows whose upwind weight approaches 1.
!
! Rebuild-preserved GS state: when NEW_R_GRID rebuilds the expanded
! grid, the old J_NEXT_FRQ (R^2*J on the old R_NM grid) is stashed in
! JNEXT_OLD/R_NM_OLD and remapped onto the new grid by interpolation
! rather than zeroed -- zeroing would force a fresh LAMBDA warm-up for
! the compressing-zone coupling after every grid rebuild, and the Knop
! Gauss-Seidel state would never converge.
!
	LOGICAL, SAVE :: HAVE_JNEXT_OLD=.FALSE.
	INTEGER, SAVE :: ND_NM_OLD=0
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: R_NM_OLD(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: JNEXT_OLD(:,:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: HNEXT_OLD(:,:)
!
! Gauss-Seidel restart/checkpoint state (added 16-Sep-2026).
! J_NEXT_FRQ is process-local and was lost on every SLURM requeue;
! it is now checkpointed to MOMJNXT_GS by DUMP_MOM_JNEXT_GS (called
! from SCR_RITE_V2, always a complete-sweep boundary) and restored by
! RESTORE_MOM_JNEXT_GS after allocation.  MOM_J_GS_OK is .TRUE. only
! when this process restored a validated generation; MOM_GS_FILE_SEEN
! records that a checkpoint file existed at all (used by the warm-up
! gate to distinguish "mom solver not run" from "restore failed").
! MOM_GS_TAG_EXP is the NITSF restored by SCR_READ_V2 (SET_MOM_GS_TAG).
!
! MOM_GS_MAXCHG  : largest |new-old|/old relative change written into
!                  J_NEXT_FRQ during the current sweep (reset at
!                  FREQ_INDX=1); measures how converged the lagged
!                  field itself is.
! MOM_INEXT_CAPPED: count of non-finite XM_NM values zeroed this sweep.
! MOM_GS_CHG_FILE / MOM_GS_CAP_FILE: the same quantities recorded in
!                  the restored checkpoint (-1 = no file); the warm-up
!                  gate treats a checkpoint that is complete but still
!                  changing by >50% (or with cap events) as immature.
!
	LOGICAL, SAVE :: MOM_J_GS_OK=.FALSE.
	LOGICAL, SAVE :: MOM_J_RAN=.FALSE.
	INTEGER, SAVE :: MOM_GS_TAG_EXP=-1
	INTEGER, SAVE :: MOM_GS_NCF=-1
	REAL(KIND=LDP), SAVE :: MOM_GS_MAXCHG=0.0_LDP
	INTEGER*8, SAVE :: MOM_INEXT_CAPPED=0
	REAL(KIND=LDP), SAVE :: MOM_GS_CHG_FILE=-1.0_LDP
	INTEGER*8, SAVE :: MOM_GS_CAP_FILE=-1
!
! MOM_NEGJ_LRG  : count of SIGNIFICANT negative-J resets this sweep
!                 (|XM| > 1D-6*max|XM|; ignores the ~1D-68 boundary
!                 noise which is benign).  A real sign flip in the
!                 lagged field marks the GS state immature: the
!                 warm-up gate in SOLVE_FOR_POPS forces
!                 LAMBDA for that iteration (audit 17-Sep-2026).
!
	INTEGER*8, SAVE :: MOM_NEGJ_LRG=0
!
! Argmax of MOM_GS_MAXCHG: (depth, freq) of the largest relative
! J_NEXT_FRQ change this sweep -- localizes the residual-dominant
! cell for instability audits.
!
	INTEGER, SAVE :: MOM_GS_ARGI=0, MOM_GS_ARGL=0
!
! Bulk convergence metrics (Oct-2026): the max-norm MOM_GS_MAXCHG is
! pinned by a handful of roaming pathological cells even when the
! field is converged in bulk.  MOM_GS_NBAD5/NBAD50 count cells whose
! relative J_NEXT_FRQ change exceeds 5%/50% this sweep, and
! MOM_GS_NTOT is the cell count the metric covered -- the MOM
! analogue of the FG side's GS_NBAD/GS_NTOT bulk-maturity inputs.
!
	INTEGER*8, SAVE :: MOM_GS_NBAD5=0
	INTEGER*8, SAVE :: MOM_GS_NBAD50=0
	INTEGER*8, SAVE :: MOM_GS_NTOT=0
!
! MOM_GS_NHUGE: cells whose relative drift exceeded 1E6 this sweep.
! MAXCHG is a single-cell statistic, so one bounded-but-bad
! resonance cell (GAM~0 at the 1+sigma*K angle cancellation) pins
! the catastrophic UNSAFE veto forever; the count lets the gate
! distinguish a handful of locally-ill-conditioned cells from a
! genuine field collapse.
!
	INTEGER*8, SAVE :: MOM_GS_NHUGE=0
!
! MOM_GS_NQUAR: cells quarantined this sweep under MOM_RF_KEEPNEG --
! their stored value is frozen, so they cannot converge and must not
! count toward the residual metrics; the count stays visible so a
! growing quarantined fraction is still diagnosable.
!
	INTEGER*8, SAVE :: MOM_GS_NQUAR=0
!
! MOM_QUAR_AGE: consecutive sweeps each (depth,frequency) cell has
! spent quarantined under MOM_RF_KEEPNEG.  A cell that fails the
! post-solve guards for a structural reason is re-solved on every
! call and keeps the same stale iterate indefinitely; the per-sweep
! NQUAR cannot tell that apart from transient failures, so cells at
! MOM_QUAR_PERSIST sweeps or older are counted separately in
! MOM_GS_NPQUAR.  edited by Devin, 29-Sep-2026.
!
	INTEGER(KIND=2), ALLOCATABLE :: MOM_QUAR_AGE(:,:)
	INTEGER*8, SAVE :: MOM_GS_NPQUAR=0
	INTEGER, PARAMETER :: MOM_QUAR_PERSIST=3
!
! MOM_OWN_CMP: per (depth,frequency) ownership map under MOM_RF_OWNW
! -- TRUE marks rows owned by descending passes (compressing rows,
! plus the I=1 boundary row when A_BND<0 under OUT_BC_TYPE<=1).
! Ownership is not recoverable at pass end because Q_NM/A_BND are
! overwritten by each frequency solve, so the drift loop records it.
! MOM_KNOP_REFRESH_V1 uses it to merge each direction's best-pass
! rows into the retained composite field.  MOM_OWN_CMPH is the H-row
! equivalent keyed on GAMH_SIGNED_NM -- midpoint and cell-centre
! sigma signs differ near the transition, so J and H masks are not
! interchangeable there.  [DO_NONMONO_KNOP]
! edited by Devin, 30-Sep-2026.
!
	LOGICAL, ALLOCATABLE :: MOM_OWN_CMP(:,:)
	LOGICAL, ALLOCATABLE :: MOM_OWN_CMPH(:,:)
!
! Per-sweep failure/conditioning audit counters.  The convergence
! problem these expose is which subsystem fails first when the
! population change stalls near ~30%: a DGTSV failure, a loss of
! diagonal dominance, a negative/runaway J repair, or a large
! ESEC/CHI overshoot leave different signatures, and silently
! repaired fields cannot be distinguished from converged ones
! without these counts.  Reset at FREQ_INDX=1, printed by
! DUMP_MOM_JNEXT_GS.
!
	INTEGER*8, SAVE :: MOM_N_DGTSV=0
	INTEGER*8, SAVE :: MOM_N_NEGJ_T=0
	INTEGER*8, SAVE :: MOM_N_NONF=0
	INTEGER*8, SAVE :: MOM_N_RUNAWAY=0
	INTEGER*8, SAVE :: MOM_TB_NONPOS=0
!
! Per-sweep budget on negJ rows written to unit 48.  A cold-start
! sweep can reject a large fraction of the ND_NM*NCF cells, and an
! uncapped per-cell dump turns the CMF solve into a pure I/O job
! (tens of GB per sweep).  Counts are exact; only the writes stop.
!
	INTEGER*8, SAVE :: MOM_NEGJ_DUMP=0
	INTEGER*8, PARAMETER :: MOM_NEGJ_DUMP_MAX=20000
!
! MOM_MIN_DR is the SIGNED ratio TB/(|TA|+|TC|); it answers "is the
! diagonal positive and dominant" and goes negative when TB<=0 (which
! MOM_TB_NONPOS counts).  MOM_MIN_DRA is the ABS ratio
! |TB|/(|TA|+|TC|); it answers "how dominant is the diagonal in
! magnitude" and stays positive when TB flips sign.  Both are needed:
! a negative-signed but large-magnitude row is ill-conditioned in a
! different way than a small-magnitude one.
!
	REAL(KIND=LDP), SAVE :: MOM_MIN_DR=1.0E+30_LDP
	INTEGER, SAVE :: MOM_MIN_DR_I=0
	INTEGER, SAVE :: MOM_MIN_DR_L=0
	REAL(KIND=LDP), SAVE :: MOM_MIN_DRA=1.0E+30_LDP
	INTEGER, SAVE :: MOM_MIN_DRA_I=0
	INTEGER, SAVE :: MOM_MIN_DRA_L=0
	REAL(KIND=LDP), SAVE :: MOM_EQ_RES=0.0_LDP
	INTEGER, SAVE :: MOM_EQ_RES_I=0
	INTEGER, SAVE :: MOM_EQ_RES_L=0
	REAL(KIND=LDP), SAVE :: MOM_MAX_COH=0.0_LDP
	REAL(KIND=LDP), SAVE :: MOM_MAX_GAM=0.0_LDP
	REAL(KIND=LDP), SAVE :: MOM_MAX_SIGT=0.0_LDP
!
! Sweep-level maturity metrics: MOM_KNOP_REFRESH_V1 overwrites the
! MOM_GS_* counters with its per-pass map residuals, so the values the
! sweep produced are copied here at refresh entry.  The sweep drift
! (fresh XM vs stored J_NEXT) is the field change the populations
! actually see between great iterations; the refresh residual is the
! lagged map's self-consistency error, which asymptotes on marginally
! contractive compressing cells and is not a maturity signal.
!
	REAL(KIND=LDP), SAVE :: MOM_SW_MAXCHG=0.0_LDP
	INTEGER*8, SAVE :: MOM_SW_NBAD5=0
	INTEGER*8, SAVE :: MOM_SW_NBAD50=0
	INTEGER*8, SAVE :: MOM_SW_NTOT=0
	INTEGER*8, SAVE :: MOM_SW_NEGJ=0
	INTEGER*8, SAVE :: MOM_SW_NQUAR=0
	INTEGER*8, SAVE :: MOM_SW_CAPPED=0
	INTEGER*8, SAVE :: MOM_SW_NHUGE=0
!
! MOM_RF_*: post-refresh state of the retained J_NEXT_FRQ field,
! measured as |post-refresh - pre-refresh| on the same normalization
! as the sweep metric.  MOM_SW_* records the sweep drift (fresh solve
! vs stored field, what the populations saw this iteration); MOM_RF_*
! records how far that sweep field was from the converged refresh
! fixed point.  When the refresh converged, MOM_RF_* is the maturity
! signal for the retained field: a converged refresh whose retained
! store still moved means the populations linearized against an
! immature field this iteration but will see the fixed point next
! sweep.  MOM_RF_RAN latches that the refresh executed at all this
! great iteration; reset at FREQ_INDX=1 on non-refresh calls.
!
	LOGICAL, SAVE :: MOM_RF_RAN=.FALSE.
	LOGICAL, SAVE :: MOM_RF_CONV=.FALSE.
	REAL(KIND=LDP), SAVE :: MOM_RF_DCHG=0.0_LDP
	INTEGER*8, SAVE :: MOM_RF_N5=0
	INTEGER*8, SAVE :: MOM_RF_N50=0
!
! MOM_NONMM_CNT: interior rows this sweep with TC_NM>0, i.e. cells
! where |EPS|>HU and the tridiagonal loses its M-matrix structure
! (positive superdiagonal).  MOM_EPS_FWD: midpoints where
! sigma*NMID_ON_J<0 took the forward EPS difference.  MOM_EPS_RAT:
! max |EPS|/HU over the sweep.  Swept up per frequency call, reset
! at FREQ_INDX=1, printed by DUMP_MOM_JNEXT_GS.
!
	INTEGER*8, SAVE :: MOM_NONMM_CNT=0
	INTEGER*8, SAVE :: MOM_EPS_FWD=0
	REAL(KIND=LDP), SAVE :: MOM_EPS_RAT=0.0_LDP
!
! MOM_BND_HBC_FL: frequencies this sweep where the outer boundary
! Eddington factor measured by the formal solve was negative or NaN
! and the boundary row fell back to the zero-flux closure (HBC_ROW
! block in the OUT_BC_TYPE<=1 row).  BND_FL_OPEN latches the lazy
! open of MOM_BND_HBC_FLOOR (append-mode per-event log).
!
	INTEGER*8, SAVE :: MOM_BND_HBC_FL=0
	LOGICAL, SAVE :: BND_FL_OPEN=.FALSE.
!
! Under-relaxation factor for the lagged J_NEXT_FRQ Gauss-Seidel
! write, matching INEXT_RELAX in FG_J_CMF_NONMONO_V1.  With strong
! compression the undamped update (straight replacement of the
! upwind J fed into the next sweep's compressing source) amplifies
! the SE<->RT loop; writing OLD + JNEXT_RELAX*(NEW-OLD) damps the
! fixed-point map without moving its fixed point.  Cells with no
! prior value (OLD<=0) are filled directly.
!
	REAL(KIND=LDP), SAVE :: JNEXT_RELAX=0.5_LDP
!
! Per-cell relaxation schedule for the refresh write.  Cells whose
! update returns near its same-parity-pass value while still making
! large steps are cycling (a period-2 mode of the alternating-
! direction map); their factor is halved toward MOM_DAMP_MIN.  The
! damped multiplier 1-r*(1-lambda) contracts for r < 2/(1-lambda),
! so a floor of 0.02 covers lambda down to ~-100; truly divergent
! cells (lambda > 1) cannot be damped and stay caught by the
! runaway/abort guards.  Persisted across great iterations because
! the marginal cells track the fixed compression zone; the schedule
! is not checkpointed, so a requeue re-derives it over ~3 passes.
!
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: MOM_DAMP_NM(:,:)
	REAL(KIND=LDP), PARAMETER :: MOM_DAMP_MIN=0.02_LDP
!
! Vector Anderson acceleration on the composite refresh map.  At fixed
! ETA/CHI the pass map is affine in the lagged field, so AA on the
! residual sequence is equivalent to GMRES on (I-A)x=b and does not
! require the map to contract (Walker & Ni 2011 equivalence).  Because
! the pass direction alternates, the iterated object is the 2-pass
! composite map; AA history is therefore taken at even pass
! boundaries.  MOM_AA_M secant pairs are kept.  When enabled it
! supersedes the per-cell Aitken step (double extrapolation would
! overshoot); the cycle-detection damp schedule still runs since AA
! needs a bounded raw iterate sequence to mix from.
!
	LOGICAL, SAVE :: MOM_RF_AA=.TRUE.
	INTEGER, SAVE :: MOM_AA_M=3
!
! Outer (call-level) Anderson acceleration state -- applied at the
! exit of MOM_KNOP_REFRESH_V1.  The iterates are the committed lagged
! stores and the residuals their successive differences, so the
! accelerated map is the whole sweep+refresh outer composition,
! including the J<->Eddington coefficient drift that the inner
! (frozen-coefficient) pass AA cannot see.  RETJ_PREV/RETH_PREV hold
! the field the last outer map was applied to.  Gated on
! [MOM_AA_OUT] (default T inside the DO_NONMONO_KNOP path);
! [MOM_AA_OUT_M] sets the history depth (default 8).
!
	LOGICAL, SAVE :: MOM_AA_OUT=.TRUE.
	INTEGER, SAVE :: MOM_AA_OUT_M=8
	INTEGER, SAVE :: NHIST_OUT=-1
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: XHIST_OUT(:,:,:),
	1    FHIST_OUT(:,:,:),RETJ_PREV(:,:),RETH_PREV(:,:)
!
! MOM_FAC_RLX (VADAT real, default 1.0 = off): under-relax the
! Eddington-factor update between great iterations.  With pops held,
! J -> formal solve -> K_ON_J/NMID factors -> coefficients is the
! only channel that moves the lagged map's fixed point between
! iterations, and on this model the unrelaxed snap displaces the map
! by >50% at ~2E5 cells per sweep while the refresh already reaches
! a fixed point at the old coefficients.  Blending the update into
! the previous factor store damps that composite map; the fixed
! point itself is unchanged.  The store is per-frequency, seeded
! with the first call's factors and not checkpointed -- a requeue
! re-seeds from the snap and re-relaxes.  1.0 keeps the raw formal
! snap (all prior behavior).
!
	REAL(KIND=LDP), SAVE :: MOM_FAC_RLX=1.0_LDP
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: FACR_KONJ(:,:),
	1    FACR_NMJ(:,:),FACR_NMH(:,:)
	LOGICAL, ALLOCATABLE, SAVE :: FACR_SEED(:)
!
! MOM_FREQ_MATCH (VADAT key DO_NONMONO_FREQ_MATCH, default FALSE;
! set by SET_MOM_FREQ_MATCH): at compressing cells the lagged
! upwind frequency index is rematched to the bin containing the
! characteristic foot, n = |GAM|*dR/dln(nu), and the frequency
! derivative spans n*dln(nu).  An explicit one-bin upwind is a
! Courant violation in frequency space once the per-cell comoving
! drift exceeds the grid spacing -- the oscillatory-eigenvalue
! regime of Ceulemans et al. 2025.  NM_J/NM_H hold the matched
! index offsets (0 = off-grid, drop coupling); DLN_J/DLN_H the
! matched log spacings; MOM_NU_GRID caches the frequency grid so
! the matched span uses true log nu, with the uniform-spacing
! fallback n*DLNM.
!
	LOGICAL, SAVE :: MOM_FREQ_MATCH=.FALSE.
	INTEGER, ALLOCATABLE, SAVE :: MOM_NM_J(:)
	INTEGER, ALLOCATABLE, SAVE :: MOM_NM_H(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: MOM_DLN_J(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: MOM_DLN_H(:)
	REAL(KIND=LDP), ALLOCATABLE, SAVE :: MOM_NU_GRID(:)
!
! Fraction of cells newly re-damped in one pass that invalidates the
! AA secant history.  edited by Devin, 28-Sep-2026: damping a handful
! of the ~12M cells does not change the global map enough to matter,
! so a restart on any single damped cell would keep NHIST_AA at ~1
! forever and AA would never mix.  1E-4 (~1200 cells) resets only on
! a genuine map change.
!
	REAL(KIND=LDP), PARAMETER :: MOM_AA_DMP_FRAC=1.0E-4_LDP
!
! MOM_AA_GCAP: per-cell cap on the weighted residual/secant admitted
! to the AA Gram system.  Cells beyond ~100% relative change per
! composite step are outside the affine-map regime (quarantine,
! limiters, dampers are active there) and are fitted by the cycle
! machinery, not the global secant solve.
!
	REAL(KIND=LDP), PARAMETER :: MOM_AA_GCAP=1.0_LDP
!
	REAL(KIND=LDP), ALLOCATABLE :: CON_GAM_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: CON_GAMH_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: AV_SIGMA_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: AV_SIGMA_TRUE_NM(:)
!
! Non-monotonic generalized opacity and coupling arrays.
!
	REAL(KIND=LDP), ALLOCATABLE :: CHI_HAT_NM(:)    !chi + |GAM|/dLOG_NU (always > 0)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_MAG_NM(:)    !|GAM| (always >= 0)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_SIGNED_NM(:) !signed GAM (can be < 0)
	REAL(KIND=LDP), ALLOCATABLE :: Q_NM(:)          !(GAM/dLOG_NU)/chi_hat (signed)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH_MAG_NM(:)   !|GAMH| for H-equation
	REAL(KIND=LDP), ALLOCATABLE :: GAMH_SIGNED_NM(:)!signed GAMH
	REAL(KIND=LDP), ALLOCATABLE :: WH_MAG_NM(:)     !|W| = |GAMH|/dLOG_NU (always > 0)
!
	REAL(KIND=LDP), ALLOCATABLE :: DTAU_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: DTAUONQ_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: TA_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: TB_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: TC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: QF_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: XM_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: XM_SAVE_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: SOURCE_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: VB_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: VC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: HU_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: HL_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: HS_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: COH_VEC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH_NM(:)
!
! GAM_FAC_NM / GAMH_FAC_NM are the sign-preserved, floored factors
! 1 + sigma*K_on_J (depth points) and 1 + AV_SIGMA*NMID_ON_HMID
! (midpoints) used by the Knop sign split, chi_hat, Q, W, PSI and the
! explicit upwind source.  They are filled every call in the
! coefficient block below; module scope is required because the
! allocation block only runs on grid rebuild while the factor arrays
! are consumed on every call.
!
	REAL(KIND=LDP), ALLOCATABLE :: GAM_FAC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH_FAC_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAM_FAC_PREV_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: GAMH_FAC_PREV_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: W_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: WPREV_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: PSI_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: PSIPREV_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS_PREV_NM(:)
!
	REAL(KIND=LDP) VDOP_FRAC_SAVE_NM
!
	INTEGER ND_NM
!
	INTEGER, ALLOCATABLE :: R_PNT_NM(:)
!
	INTEGER, ALLOCATABLE :: J_INDX_NM(:)
	INTEGER, ALLOCATABLE :: H_INDX_NM(:)
!
	LOGICAL VERBOSE_NM
	LOGICAL FIRST_TIME_NM
	DATA FIRST_TIME_NM/.TRUE./
	DATA VDOP_FRAC_SAVE_NM/-10001.1_LDP/    !Absurd value
!
! ----------------------------------------------------------------------
! Per-frequency Knop refresh store (Oct-2026), moment-equation analog
! of the FG_J_CMF_NONMONO_V1 per-ray refresh.  The moment solver has
! the same one-bin-per-great-iteration lag: compressing depths read
! J_NEXT_FRQ(:,k+1) written during the PREVIOUS sweep.  The refresh
! converges the coupling by re-driving MOM_J_CMF_NONMONO_V1 itself in
! alternating frequency order on the captured per-frequency inputs
! (same discretisation, DGTSV solve and fallbacks by construction).
!
! MOM_NU_NEXT: nu(k+1) for the current call, set by SET_MOM_NU_NEXT so
! compressing depths use the forward spacing ln(nu_k/nu_{k+1}).
! MOM_KN_REFRESHING: set while MOM_KNOP_REFRESH_V1 re-drives the
! routine -- skips GET_MOMS_NON_REL (whose EDDFACTOR-backed inputs
! only exist for the last sweep frequency) and the capture write.
!
	REAL(KIND=LDP) MOM_NU_NEXT
	LOGICAL MOM_KN_REFRESHING
!
! KN_M_DESC: set by MOM_KNOP_REFRESH_V1 per pass -- TRUE while the
! pass runs in descending frequency order (LDI=-1).
!
! MOM_RF_GS_DESC (VADAT key, default FALSE): on a descending refresh
! pass the compressing-depth read of slot L+1 uses the live
! J_NEXT_FRQ (already updated this pass) instead of the frozen
! JNEXT_RF_READ snapshot -- true Gauss-Seidel propagation along the
! red-to-blue compressing characteristic.  The frozen read limits
! that propagation to one bin per pass, which cannot converge a
! coupling spanning more bins than KN_M_MAXPASS.  The earlier
! cascade divergence (job 11887664) predates the M-matrix/EPS,
! COH-floor and HBC-floor fixes; the abort guard plus the
! pre-refresh restore contain a recurrence if it is still unstable.
!
! MOM_RF_KEEPNEG (VADAT key, default TRUE): cells whose XM hit a
! post-solve repair (significant negJ, non-finite or runaway
! rejection) keep their previous J_NEXT_FRQ iterate instead of
! storing the repaired value -- the repair is a display-level fix,
! not a fixed-point output, and storing it injects discontinuous
! noise the lagged map then chases.  In the repair itself XM_FB_NM
! (the previous valid field) is already written into XM_NM, so
! this flag only controls whether the freeze shows up in J_NEXT_FRQ
! as an un-written cell vs the relaxed write of the fallback.
!
! MOM_RF_BULK_EXIT (VADAT key, default TRUE): the refresh max-norm
! asymptotes on a few marginally contractive compressing cells
! (|Q|/(chi+|Q|)->1), so the refresh otherwise burns passes until
! the stall exit retains the same bulk-converged field.  Once the
! per-pass drift satisfies the bulk bounds the gate applies to the
! sweep (fewer than 0.1% of cells moving >50%, no capped writes,
! negJ and quarantined fractions under 10%), exit as converged.
! CONVERGED=T routes the gate to the retained-vs-sweep delta
! (MOM_RF_*), which still rejects a fixed point that sits far from
! the field the populations saw.
!
! MOM_RF_ACTIVE (VADAT key, default FALSE): active-set
! front-tracking over frequencies in the refresh.  A frequency's
! solve only depends on the lagged store at slot L-1 (expanding
! cells) and slots L+MOM_NM_{J,H}(I) (compressing cells), so a
! converged frequency can be skipped while its inputs are quiet.
! Each pass records the per-frequency max drift; the next pass
! keeps L active when it was hot, its L-1 input was hot, or any
! slot in (L,L+max NM] was hot -- a conservative dilation of the
! hot set over every offset any compressing cell reads.  The mask
! resets to all-active each refresh because the map's frozen
! coefficients are rebuilt between calls.
!
	LOGICAL KN_M_DESC
	LOGICAL MOM_RF_GS_DESC
	LOGICAL MOM_RF_KEEPNEG
	LOGICAL MOM_RF_BULK_EXIT
	LOGICAL MOM_RF_ACTIVE
!
! MOM_RF_DIR (VADAT key, default 'AUTO'): refresh pass ordering.
! Knop Paper II's GS solve is physically guided -- the pass follows
! the direction information actually propagates in (r,nu) space.  On
! a grid where nearly every coupling row is compressing (m15_2p82d:
! sigma < -1 at all 70 depths) an ascending pass reads only the
! frozen pass-start snapshot -- a pure Jacobi step that re-injects a
! full pass of error -- so the alternating composite map stalls
! (job 12045050: negJ ~2.4M/12M cells, refresh residual pinned at
! ~1E4).  'AUTO' measures the compressing-row fraction on the sweep
! (MOM_SW_NCROW/MOM_SW_NROW) and runs only the dominant direction
! once one side holds >=90% of rows; 'ALT' forces the original
! alternating order, 'ASC'/'DESC' force a single direction.
!
! MOM_DOM_CMP / MOM_DOM_EXP: the resolved dominance state, set at the
! FREQ_INDX=1 reset from the previous sweep's counters and refreshed
! at MOM_KNOP_REFRESH entry.  The dominant-compressing state also
! suppresses the in-sweep J_NEXT_FRQ/H_NEXT_FRQ write at compressing
! rows: the ascending sweep's compressing upwind read is a
! one-iteration-old Jacobi value, so writing the solve output back
! at JNEXT_RELAX erodes the refresh-converged store every great
! iteration.  The descending refresh owns the compressing component
! of the lagged field; expanding rows are unaffected (their in-sweep
! read is fresh).  Suppression requires MOM_RF_ENABLED (the
! DO_KNOP_REFRESH setting): with no refresh running the suppressed
! cells would freeze at seed values forever.
!
	CHARACTER(LEN=4) MOM_RF_DIR
	LOGICAL MOM_DOM_CMP,MOM_DOM_EXP,MOM_RF_ENABLED
	INTEGER*8, SAVE :: MOM_SW_NCROW=0
	INTEGER*8, SAVE :: MOM_SW_NROW=0
	REAL(KIND=LDP), SAVE :: MOM_SW_FCMP=0.0_LDP
!
! MOM_EPS_LIM (VADAT key, default TRUE): donor-cell positivity
! limiter on the lagged EPS midpoint transfer.  The explicit
! lagged-frequency term scatters +-2*EPS_PREV*J_up between adjacent
! rows; at advection-dominated cells (|EPS|/HU ~ 20) the debit can
! exceed the donor row's accumulated positive RHS, returning
! negative J from a well-conditioned M-matrix -- the persistent negJ
! population of job 12045050.  The transfer is capped at the donor
! row's positive balance and the recipient is credited only what was
! debited, so the pair stays conservative.  MOM_EPS_CAP counts
! limited transfers.
!
	LOGICAL MOM_EPS_LIM
	INTEGER*8, SAVE :: MOM_EPS_CAP=0
	INTEGER*8, SAVE :: MOM_EPS_CAP_SW=0
!
! MOM_RF_OWNW (VADAT key, default TRUE): pass-direction ownership of
! lagged-store writes during MOM_KNOP_REFRESH.  On a descending pass
! only compressing rows (upwind slot L+1 already visited) see a fresh
! read; expanding rows read the not-yet-visited L-1 store, i.e. a
! Jacobi step, and vice versa on ascending passes.  Writing the
! stale-read result back every other pass dilutes the fresh-read
! subset and is the residual floor of the alternating map
! (residual ~1E4 plateau, job 12046369).  With OWNWrite each row only
! ever takes Gauss-Seidel steps in its own direction: the composite
! map is the per-direction GS operator instead of GS*Jacobi.
! MOM_RF_NSUPP counts suppressed writes per refresh pass;
! MOM_SW_NSUPP counts main-sweep dominant-mode suppressions.
!
	LOGICAL MOM_RF_OWNW
	INTEGER*8, SAVE :: MOM_RF_NSUPP=0
	INTEGER*8, SAVE :: MOM_SW_NSUPP=0
!
! MOM_GAM_PIN (VADAT key, default TRUE): evaluate the lagged
! comoving-gradient factors (GAM_FAC_PREV, GAMH_FAC_PREV, the
! Eddington boundary PSIPREV) at the CURRENT frequency's Eddington
! factors instead of the bluer neighbour's.  With independent
! evaluations a sign or magnitude drift across a frequency boundary
! near the 1+sigma*K=0 resonance lets the explicit upwind weight
! exceed the implicit diagonal (PSIPREV/PSI > 1), and the lagged
! fixed-point map then amplifies the field it reads.  Pinning makes
! the explicit weight identical to the implicit one -- the
! compressing branch's existing PSI convention -- so neither upwind
! channel is expansive by construction.  edited by Devin,
! 12-Oct-2026.
!
	LOGICAL MOM_GAM_PIN
!
! MOM_LAG_BND (VADAT key, default TRUE): per-row bound on the
! lagged frequency source.  The lagged fixed-point map amplifies
! wherever a row's lagged coupling weight -- upwind J coefficient,
! HS*H_up, and the two EPS_PREV midpoint transfers -- exceeds its
! diagonal excess TB-|TA|-|TC|.  Shock cells violate this: PSI
! collapses at the 1+sigma*K=0 resonance while the H and EPS
! channels still feed the row (m15 probes, 12-Oct-2026).  The
! row's whole lagged contribution is scaled by
! min(1, excess/weight); rows whose unbounded map was already
! contractive are untouched.  edited by Devin, 12-Oct-2026.
!
	LOGICAL MOM_LAG_BND
!
! MOM_ABS_FLR (VADAT key, default 0.0 = off): lower bound on the
! absorption fraction 1-COH in the tridiagonal diagonal.  At
! ESEC>=CHI cells the MAX(0,1-COH) clip removes the only
! J-independent diagonal weight; with PSI~0 near the 1+sigma*K=0
! resonance such rows have no diagonal excess, so MOM_LAG_BND nulls
! their entire lagged source and the row sloshes between its clipped
! and unclipped states instead of converging (KNOP probes,
! 13-Oct-2026: mom maxchg floored ~1E4 at sigma_t~-100 cells).  A
! small positive floor keeps a physical absorption anchor in the row
! so the bound can leave a nonzero lagged weight.  0.0 restores the
! exact clip; ~1D-3 retains 0.1% absorption.  edited by Devin,
! 13-Oct-2026.
!
	REAL(KIND=LDP) MOM_ABS_FLR
!
! MOM_H_IMPL (VADAT key, default FALSE): implicitize the lagged
! upwind-H divergence in the J row by adding |VB|+|VC| (the
! magnitude of the HS*H_up midpoint couplings) to the diagonal.
! At shock cells the MOM_LAG_BND excess TB-|TA|-|TC| is ~0, so
! the lagged source is nulled and the row oscillates between its
! clipped local solution and the lagged field -- the seed delta
! (~1E1 at the retained-vs-sweep metric) and the ~150 capped
! writes per refresh pass both live on these rows.  Adding the
! coupling magnitude to the diagonal damps the row's response to
! both neighbours and lagged terms, and inside MOM_LAG_BND it
! raises the excess by the H weight so less of the lagged source
! is clipped.  Consistency cost: the signed lagged term
! -HS(I-1)*H(I-1)+HS(I)*H(I) is a difference, while the added
! diagonal is a sum, so smooth rows carry a bounded extra
! diagonal mass (~2*HS <= 2).  Default F keeps the proven
! (a)+(c) behaviour; enable only for the KNOP stall test.
! edited by Devin, 13-Oct-2026.
!
	LOGICAL MOM_H_IMPL
!
! LAGB diagnostics: rows scaled this sweep, and the smallest scale
! applied (1.0 = the bound was inactive everywhere).
!
	INTEGER*8, SAVE :: MOM_LAGB_CNT=0
	REAL(KIND=LDP), SAVE :: MOM_LAGB_MIN=1.0_LDP
!
! XM_FB_NM: previous valid R^2*J at the current frequency, assembled
! before DGTSV overwrites XM_NM.  It is the fallback target for a
! failed solve or a rejected cell: the lagged GS iterate
! J_NEXT_FRQ(:,FREQ_INDX) where positive, else the local-equilibrium
! estimate SOURCE*R^2.  The bluer-frequency RSQJNU_PREV_NM is
! deliberately excluded -- see the assembly block.  The old code
! restored XM_SAVE_NM (the assembled RHS)
! after a DGTSV failure -- a RHS is not a radiation field, so that
! injected source-scale values into R^2*J and propagated them through
! J_NEXT_FRQ, RSQJNU_PREV and BA.
!
	REAL(KIND=LDP), ALLOCATABLE :: XM_FB_NM(:)
!
! XM_BAD_NM: marks cells repaired by the post-solve guards this call.
!
	LOGICAL, ALLOCATABLE :: XM_BAD_NM(:)
!
! MOM_CELL_QUAR: XM_BAD_NM mapped onto the SM (model) grid for the
! last non-refresh call; read by DO_VAR_CONT through
! GET_MOM_CELL_QUAR to suppress the ALO diagonal at frozen cells.
!
	LOGICAL, ALLOCATABLE :: MOM_CELL_QUAR(:)
!
! Diagonal response of the moment solution on the model grid, for the
! USE_LAM_ES=F linearization.  Neighbours and the lagged upwind
! (J_PREV or J_NEXT) are held fixed, so a compressing row, whose
! diagonal carries PSI = DTAUONQ*|Q|, reports only its local source
! response.  edited by Sean Brennan (25-Sep-2026)
!
	REAL(KIND=LDP), ALLOCATABLE :: MOM_ALO_DCHI(:)
	REAL(KIND=LDP), ALLOCATABLE :: MOM_ALO_DETA(:)
	LOGICAL :: MOM_ALO_VALID=.FALSE.
	INTEGER :: MOM_ALO_ND=0
!
! JNEXT_RF_READ: pass-start snapshot of J_NEXT_FRQ used while
! MOM_KN_REFRESHING.  The in-sweep compressing read is pure Jacobi
! (the previous great iteration's field); without the snapshot a
! descending refresh pass consumes this-pass writes -- a GS cascade
! that amplified residuals to ~1E10 (job 11887664, 20-Sep-2026).
!
	REAL(KIND=LDP), ALLOCATABLE :: JNEXT_RF_READ(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: HNEXT_RF_READ(:,:)
!
! Captured per-frequency inputs.  SM-grid arrays are the routine's
! arguments (ETA_SM/CHI_SM/ESEC_SM); NM-grid arrays are the
! post-GET_MOMS Eddington factors.  All double precision: ETA/CHI are
! log-interpolated in the re-solve, and valid double-precision values
! below ~1E-38 (cold faint cells) would underflow to exact 0.0 in
! single precision, turning LOG(ETA) into -Inf and NaN-ing the solve.
!
	REAL(KIND=LDP), ALLOCATABLE :: KN_M_CHISM(:,:),KN_M_ETASM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_M_ESECSM(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MKONJ(:,:),KN_MNMJ(:,:),
	1     KN_MNMH(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MHBC(:),KN_MNBC(:),KN_MINHBC(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MJPIB(:),KN_MJMIB(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MHIB(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MDBB(:),KN_MIC(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MDLN(:),KN_MNUNXT(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_MFREQ(:)
!
! OUT_BC_TYPE is per-frequency: COMP_J_BLANK forces 1 above
! OUT_BC_PARAM_ONE and uses RD_OUT_BC_TYPE below it.  Stored per
! frequency so the refresh replays the same boundary type the sweep
! used; a single scalar was only correct for RD_OUT_BC_TYPE=1.
!
	INTEGER, ALLOCATABLE :: KN_M_OUTBC(:)
!
! Sweep-fixed arguments captured once.
!
	REAL(KIND=LDP), ALLOCATABLE :: KN_M_RSM(:),KN_M_VSM(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_M_SIGSM(:),KN_M_SIGTSM(:)
	REAL(KIND=LDP), ALLOCATABLE :: KN_M_VDOPV(:)
	REAL(KIND=LDP) KN_M_IBSTAB,KN_M_VDOPF
	INTEGER KN_M_NC,KN_M_NP,KN_M_NDSM,KN_M_NDNM,KN_M_NCF
	CHARACTER(LEN=20) KN_M_NTYPE,KN_M_HCHK,KN_M_INBM
	CHARACTER(LEN=6) KN_M_METHOD
	LOGICAL KN_M_COH
!
	LOGICAL KN_M_READY
	LOGICAL KN_M_RF_WARNED
	LOGICAL KN_M_ANYCMP
	REAL(KIND=LDP) KN_M_TOL
	INTEGER KN_M_MAXPASS
!
! KN_NM_MAX_INS: per-original-cell insertion cap for the expanded
! grid.  |dV|/VDOP can request arbitrarily many sub-points in one
! interval on a strong shock; the cap keeps ND_NM bounded (and so
! the ND_NM x NCF lagged stores) and is reported, not silent.
!
	INTEGER KN_NM_MAX_INS
!
	DATA MOM_NU_NEXT /-1.0_LDP/
	DATA MOM_KN_REFRESHING /.FALSE./
	DATA KN_M_DESC /.FALSE./
	DATA MOM_RF_GS_DESC /.FALSE./
	DATA MOM_RF_KEEPNEG /.TRUE./
	DATA MOM_RF_BULK_EXIT /.TRUE./
	DATA MOM_RF_ACTIVE /.FALSE./
	DATA MOM_RF_DIR /'AUTO'/
	DATA MOM_RF_OWNW /.TRUE./
	DATA MOM_DOM_CMP /.FALSE./
	DATA MOM_DOM_EXP /.FALSE./
	DATA MOM_RF_ENABLED /.FALSE./
	DATA MOM_EPS_LIM /.TRUE./
	DATA MOM_GAM_PIN /.TRUE./
	DATA MOM_LAG_BND /.TRUE./
	DATA MOM_ABS_FLR /0.0_LDP/
	DATA MOM_H_IMPL /.FALSE./
	DATA KN_M_READY /.FALSE./
	DATA KN_M_RF_WARNED /.FALSE./
	DATA KN_M_ANYCMP /.FALSE./
	DATA KN_M_TOL /2.0E-03_LDP/
	DATA KN_M_MAXPASS /60/
	DATA KN_NM_MAX_INS /256/
!
	END MODULE MOD_MOM_J_NONMONO_V1
!
!
!
! MOM_J_CMF_NONMONO_V1: Non-monotonic moment solver for the mean
! intensity J at a single frequency in the comoving frame.
!
! Based on MOM_J_CMF_V11 with the Knop et al. (2009) generalized-opacity
! split for arbitrary (non-monotonic) velocity fields.
!
! Key differences from MOM_J_CMF_V11:
!   1. Uses SIGMA_TRUE (unclamped) instead of clamped SIGMA
!   2. Generalized opacity chi_hat = chi + |GAM|/dLOG_NU (always > 0)
!   3. Velocity coupling moved to source: Q*J_upwind (signed Q)
!   4. Sign-dependent upwind: J_PREV (GAM>=0), J_NEXT (GAM<0)
!   5. J_NEXT stored from previous great iteration (Gauss-Seidel)
!   6. EPS*J_mid (sigma*NMID_ON_J*J) coupling: depth-upwind midpoint
!      evaluation (23-Sep-2026).  The centered (J_m+J_{m+1}) form put a
!      positive off-diagonal on one row for EITHER sign once
!      |EPS|>HU/HL -- a non-M-matrix whose sign-flipped solution seeded
!      the lagged frequency couplings and blew up the outer luminosity
!      over 2-3 great iterations (verified with THOMAS and pivoted
!      DGTSV, jobs 11948451/11950204/11951526).  The midpoint J is now
!      taken at the upwind node (outer node for EPS>=0, inner for
!      EPS<0): every EPS contribution lands sign-definite, the matrix
!      is an M-matrix again, and the term stays fully implicit.  A
!      deferred-correction split was also tried and abandoned: it
!      cannot converge while the deferred part exceeds the implicit
!      diagonal (|EPS|/HU~16-26).
!
! The tridiagonal system is well-conditioned in the diffusion limit:
!   - chi_hat > 0  =>  DTAU > 0  =>  DTAUONQ > 0
!   - PSI = DTAUONQ * |GAM|/dLOG_NU  (always positive, no SIGMA term)
!   - W   = |GAMH|/dLOG_NU          (always positive, no SIGMA term)
! The signed velocity coupling enters only through the RHS source
! term Q*J_upwind, which cannot make the tridiagonal system singular.
! With the depth-upwind EPS form the full matrix is an M-matrix even
! at strongly compressing depths; DGTSV (partial pivoting) is kept as
! the solve since it is agnostic to the stencil.
!
! edited by Sean Brennan (15-Sep-2026)
!
	SUBROUTINE MOM_J_CMF_NONMONO_V1(ETA_SM,CHI_SM,ESEC_SM,
	1                  V_SM,SIGMA_SM,SIGMA_TRUE_SM,R_SM,
	1                  JNU,RSQHNU_SM,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  FREQ,dLOG_NU,
	1                  INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1                  N_TYPE,H_CHK_OPTION,METHOD,COHERENT,OUT_BC_TYPE,
	1                  INIT,NEW_FREQ,NC,NP,ND_SM,FREQ_INDX,NCF_IN)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
	INTEGER NC,NP,ND_SM
	INTEGER FREQ_INDX,NCF_IN
	REAL(KIND=LDP) ETA_SM(ND_SM),CHI_SM(ND_SM),ESEC_SM(ND_SM)
	REAL(KIND=LDP) V_SM(ND_SM),SIGMA_SM(ND_SM),SIGMA_TRUE_SM(ND_SM)
	REAL(KIND=LDP) R_SM(ND_SM)
!
	REAL(KIND=LDP) JNU(ND_SM),RSQHNU_SM(ND_SM)
	REAL(KIND=LDP) VDOP_VEC(ND_SM),VDOP_FRAC
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL(KIND=LDP) MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
	INTEGER OUT_BC_TYPE
	REAL(KIND=LDP) DBB,IC,FREQ,dLOG_NU
	REAL(KIND=LDP) IB_STAB_FACTOR
	CHARACTER(LEN=*) INNER_BND_METH
	CHARACTER(LEN=*) H_CHK_OPTION
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
!
	LOGICAL INIT,COHERENT,NEW_FREQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	REAL(KIND=LDP) T1,T2,T3,T4,T5
	REAL(KIND=LDP) A_BND
	REAL(KIND=LDP) HBC_ROW,HBC_PREV_ROW
	REAL(KIND=LDP) NBC_ROW,NBC_PREV_ROW
	LOGICAL BND_HBC_FLOORED
	REAL(KIND=LDP)  DELTA_R
	REAL(KIND=LDP)  DLNM
	REAL(KIND=LDP) GAM_MOM, GAM_MOM_MAG, Q_MOM
!
! edited by Devin, 28-Sep-2026: regularize the Knop moment closure
! at the zero-GAM resonance 1 + sigma*K_on_J = 0, where the implicit
! diagonal PSI collapses and the fixed-point map stops contracting.
! The limit is dimensionless; the sign is preserved so the upwind
! branch is unchanged away from the resonance.
!
	REAL(KIND=LDP), PARAMETER :: MOM_GAM_FAC_MIN=1.0E-2_LDP
!
! Minimum log-frequency step used in the advection coefficients.
! The continuum grid inserts bound-free edge doublets at dLOG_NU =
! 2E-11 (m15 SOB: 1774 of 9665 points); every frequency-advection
! coefficient carries 1/DLN and reaches ~1E8x the physical terms
! there.  V11 survives because its implicit and explicit terms
! share coefficients and the row degenerates to J_k = J_{k-1}; the
! sign-split discretization (lagged upwind stores, EPS limiter)
! does not.  1E-6 bounds the coefficients while leaving every real
! spacing (min ~1E-5) untouched.
!
	REAL(KIND=LDP), PARAMETER :: MOM_DLN_MIN=1.0E-6_LDP
	REAL(KIND=LDP), ALLOCATABLE :: JUP_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: HUP_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: EPS_SCL_NM(:)
!
! XLAG_NM: the lagged part of each row's assembled RHS (everything
! except the DTAUONQ*SOURCE*R^2 local source).  LAG_SCL_NM: the
! per-row bound factor applied under MOM_LAG_BND, reused by the
! negative-J rebuild so the rebuild matches the equation that was
! actually solved.
!
	REAL(KIND=LDP), ALLOCATABLE :: XLAG_NM(:)
	REAL(KIND=LDP), ALLOCATABLE :: LAG_SCL_NM(:)
!
! Work arrays for the pivoted tridiagonal solve (DGTSV destroys its
! matrix arguments, so the coefficients are copied before each call).
!
	REAL(KIND=LDP), ALLOCATABLE :: GT_DL(:),GT_D(:),GT_DU(:)
	INTEGER IT1,I,J,K,N_CAP
	INTEGER, PARAMETER :: IONE=1
	INTEGER IOS
	INTEGER GT_INFO
	EXTERNAL DGTSV
	INTEGER ICNT
	INTEGER N_QBND
	INTEGER, SAVE :: N_QPRT=0
	INTEGER, SAVE :: N_XPRT=0
	REAL(KIND=LDP), ALLOCATABLE :: R_QR_CL(:)
	REAL(KIND=LDP) R_HI,R_LO
	LOGICAL NEW_R_GRID
	LOGICAL HAS_COMPRESSING
	LOGICAL L_SUPP
!
	NEW_R_GRID=.FALSE.
	MOM_ALO_VALID=.FALSE.
!
! LUER is only assigned inside conditional blocks below; the DGTSV
! failure report could reach it unset.  Fix it once per call.
!
	LUER=ERROR_LU()
	IF(INIT .AND. ALLOCATED(R_NM))THEN
	  DO I=1,ND_SM
! R is re-derived each iteration and can differ at the few-ulp level
! between consecutive calls with no physical change.  Rebuild only on
! a real drift (relative change > 1e-8): bit-level noise was forcing a
! rebuild every INIT call, discarding and remapping the GS store
! unnecessarily.
	    IF(ABS(R_SM(I)-EXP(LOG_R_SM_NM(I))) .GT.
	1       1.0D-8*R_SM(I))THEN
	      NEW_R_GRID=.TRUE.
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Updating Rgrid in MOM_J_CMF_NONMONO_V1'
	      WRITE(LUER,'(A,I5,A,2ES16.8)')
	1      '  trigger: depth',I,'  R_SM vs stored logR:',
	1      R_SM(I),EXP(LOG_R_SM_NM(I))
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAVE_NM)THEN
	    NEW_R_GRID=.TRUE.
	    LUER=ERROR_LU()
	    WRITE(LUER,'(A,2ES16.8)')
	1    '  trigger: VDOP_FRAC change:',VDOP_FRAC_SAVE_NM,VDOP_FRAC
	  END IF
!
! Velocity drift also invalidates the expanded grid: the sub-point
! count is INT(|V_SM(I)-V_SM(I+1)|/T2), and V_COEF_NM/SIGMA_COEF_NM/
! SIGMA_TRUE_COEF_NM are computed only in the build block below.
! Without this check a V update with R unchanged left GAM, W, EPS,
! PSI and Q evaluated on the old velocity field.  Scale the
! tolerance by the peak |V| so near-zero cells do not fire on
! bit-level noise; 1e-6 of peak is far below the T2 sub-grid
! sensitivity (T2 ~ VDOP_FRAC*|VDOP|, order km/s).
!
	  IF(.NOT. NEW_R_GRID .AND. ALLOCATED(V_SM_SV))THEN
	    IF(SIZE(V_SM_SV) .EQ. ND_SM)THEN
	      T1=1.0D-6*MAX(MAXVAL(ABS(V_SM(1:ND_SM))),1.0E-10_LDP)
	      DO I=1,ND_SM
	        IF(ABS(V_SM(I)-V_SM_SV(I)) .GT. T1)THEN
	          NEW_R_GRID=.TRUE.
	          WRITE(LUER,'(A,I5,A,2ES16.8)')
	1            '  trigger: V_SM drift at depth',I,
	2            '  V_SM vs stored:',V_SM(I),V_SM_SV(I)
	          EXIT
	        END IF
	      END DO
	    END IF
	  END IF
	END IF
	IF(FIRST_TIME_NM)THEN
          OPEN(UNIT=48,FILE='MOM_J_NONMONO_ERRORS',STATUS='UNKNOWN')
          CALL GET_VERBOSE_INFO(VERBOSE_NM)
	ELSE IF(INIT)THEN
	  REWIND(48)
	END IF
!
! Deallocate all arrays if we have changed R grid or VDOP_FRAC.
!
	IF(ALLOCATED(R_NM) .AND. NEW_R_GRID)THEN
!
! Stash the GS frequency-coupling store before the old grid goes away;
! it is remapped onto the new expanded grid below.
!
	  IF(ALLOCATED(J_NEXT_FRQ) .AND. J_NEXT_VALID)THEN
	    ND_NM_OLD=SIZE(J_NEXT_FRQ,1)
	    IF(ALLOCATED(JNEXT_OLD))DEALLOCATE(JNEXT_OLD)
	    IF(ALLOCATED(R_NM_OLD))DEALLOCATE(R_NM_OLD)
	    ALLOCATE(JNEXT_OLD(ND_NM_OLD,SIZE(J_NEXT_FRQ,2)),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE(R_NM_OLD(ND_NM_OLD),STAT=IOS)
	    IF(IOS .EQ. 0)THEN
	      R_NM_OLD(1:ND_NM_OLD)=R_NM(1:ND_NM_OLD)
	      JNEXT_OLD(1:ND_NM_OLD,:)=J_NEXT_FRQ(1:ND_NM_OLD,:)
	      IF(ALLOCATED(H_NEXT_FRQ))THEN
	        IF(ALLOCATED(HNEXT_OLD))DEALLOCATE(HNEXT_OLD)
	        ALLOCATE(HNEXT_OLD(ND_NM_OLD,SIZE(H_NEXT_FRQ,2)),
	1          STAT=IOS)
	        IF(IOS .EQ. 0)HNEXT_OLD=H_NEXT_FRQ
	      END IF
	      HAVE_JNEXT_OLD=.TRUE.
	    END IF
	  END IF
	  DEALLOCATE ( R_NM )
	  DEALLOCATE ( R_PNT_NM )
	  DEALLOCATE ( LOG_R_SM_NM )
	  IF(ALLOCATED(V_SM_SV))DEALLOCATE ( V_SM_SV )
	  DEALLOCATE ( V_COEF_NM )
	  DEALLOCATE ( SIGMA_COEF_NM )
	  DEALLOCATE ( SIGMA_TRUE_COEF_NM )
	  DEALLOCATE ( ETA_COEF_NM )
	  DEALLOCATE ( ESEC_COEF_NM )
	  DEALLOCATE ( CHI_COEF_NM )
	  DEALLOCATE ( V_NM )
	  DEALLOCATE ( SIGMA_NM )
	  DEALLOCATE ( SIGMA_TRUE_NM )
	  DEALLOCATE ( ETA_NM )
	  DEALLOCATE ( ESEC_NM )
	  DEALLOCATE ( CHI_NM )
	  DEALLOCATE ( RSQJNU_NM )
	  DEALLOCATE ( RSQHNU_NM )
	  DEALLOCATE ( RSQJNU_PREV_NM )
	  DEALLOCATE ( RSQHNU_PREV_NM )
	  DEALLOCATE ( CON_GAM_NM )
	  DEALLOCATE ( CON_GAMH_NM )
	  DEALLOCATE ( AV_SIGMA_NM )
	  DEALLOCATE ( AV_SIGMA_TRUE_NM )
	  DEALLOCATE ( CHI_HAT_NM )
	  DEALLOCATE ( GAM_MAG_NM )
	  DEALLOCATE ( GAM_SIGNED_NM )
	  DEALLOCATE ( Q_NM )
	  DEALLOCATE ( GAMH_MAG_NM )
	  DEALLOCATE ( GAMH_SIGNED_NM )
	  DEALLOCATE ( WH_MAG_NM )
	  DEALLOCATE ( DTAU_NM )
	  DEALLOCATE ( DTAUONQ_NM )
	  DEALLOCATE ( TA_NM )
	  DEALLOCATE ( TB_NM )
	  DEALLOCATE ( TC_NM )
	  DEALLOCATE ( QF_NM )
	  DEALLOCATE ( XM_NM )
	  DEALLOCATE ( XM_SAVE_NM )
	  DEALLOCATE ( XM_FB_NM )
	  DEALLOCATE ( XM_BAD_NM )
	  DEALLOCATE ( SOURCE_NM )
	  DEALLOCATE ( VB_NM )
	  DEALLOCATE ( VC_NM )
	  DEALLOCATE ( HU_NM )
	  DEALLOCATE ( HL_NM )
	  DEALLOCATE ( HS_NM )
	  DEALLOCATE ( COH_VEC_NM )
	  DEALLOCATE ( GAM_NM )
	  DEALLOCATE ( GAMH_NM )
	  DEALLOCATE ( GAM_FAC_NM )
	  DEALLOCATE ( GAMH_FAC_NM )
	  DEALLOCATE ( GAM_FAC_PREV_NM )
	  DEALLOCATE ( GAMH_FAC_PREV_NM )
	  DEALLOCATE ( W_NM )
	  DEALLOCATE ( WPREV_NM )
	  DEALLOCATE ( PSI_NM )
	  DEALLOCATE ( PSIPREV_NM )
	  DEALLOCATE ( EPS_NM )
	  DEALLOCATE ( EPS_PREV_NM )
	  DEALLOCATE ( J_INDX_NM )
	  DEALLOCATE ( H_INDX_NM )
	  IF(ALLOCATED(J_NEXT_FRQ))DEALLOCATE(J_NEXT_FRQ)
	  IF(ALLOCATED(H_NEXT_FRQ))DEALLOCATE(H_NEXT_FRQ)
	  IF(ALLOCATED(MOM_QUAR_AGE))DEALLOCATE(MOM_QUAR_AGE)
	  IF(ALLOCATED(MOM_OWN_CMP))DEALLOCATE(MOM_OWN_CMP)
	  IF(ALLOCATED(MOM_OWN_CMPH))DEALLOCATE(MOM_OWN_CMPH)
	  IF(ALLOCATED(KN_M_CHISM))THEN
	    DEALLOCATE(KN_M_CHISM,KN_M_ETASM,KN_M_ESECSM,
	1      KN_MKONJ,KN_MNMJ,KN_MNMH,
	2      KN_MHBC,KN_MNBC,KN_MINHBC,
	3      KN_MJPIB,KN_MJMIB,KN_MHIB,
	4      KN_MDBB,KN_MIC,KN_MDLN,KN_MNUNXT,KN_MFREQ,
	5      KN_M_RSM,KN_M_VSM,KN_M_SIGSM,KN_M_SIGTSM,KN_M_VDOPV,
	6      KN_M_OUTBC)
	  END IF
	  KN_M_READY=.FALSE.
	  KN_M_ANYCMP=.FALSE.
	END IF
	VDOP_FRAC_SAVE_NM=VDOP_FRAC
!
! On the very first entry, define the expanded R grid and allocate arrays.
!
	IF( FIRST_TIME_NM .OR. .NOT. ALLOCATED(R_NM) )THEN
!
! Determine the number of points for the expanded R grid.
! We insert an EVEN number of points between grid points to resolve
! the Doppler shift.  For non-monotonic flows, |V(I)-V(I+1)| is used
! (not V(I)-V(I+1)) so that velocity reversals are handled correctly.
!
          K=1
          T2=VDOP_FRAC*MINVAL(ABS(VDOP_VEC(1:ND_SM)))
!
! T2=0 (all-zero VDOP) would divide by zero; a negative T2
! (VDOP_FRAC<0 or a negative VDOP entry) would silently skip every
! insertion.  Neither is a recoverable grid -- stop loudly.
!
	  IF(.NOT.(T2 .GT. 0.0_LDP))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'MOM_J_CMF_NONMONO_V1: invalid Doppler grid'
	1   ,' scale T2=',T2,' (VDOP_FRAC=',VDOP_FRAC,')'
	    STOP
	  END IF
!
	  N_CAP=0
          DO I=1,ND_SM-1
            IT1=INT( ABS(V_SM(I)-V_SM(I+1))/T2 )
	    IF(IT1 .GT. KN_NM_MAX_INS)THEN
	      IT1=KN_NM_MAX_INS;  N_CAP=N_CAP+1
	    END IF
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)K=K+IT1
            K=K+1
          END DO
          ND_NM=K
	  LUER=ERROR_LU()
	  WRITE(LUER,'(A,I5,A,I6,A,ES10.3,A,I5)')
	1   '  MOM grid:',ND_SM,' -> ',ND_NM,
	2   '  T2(km/s)=',T2,'  capped intervals=',N_CAP
	  IF(N_CAP .GT. 0)THEN
	    WRITE(LUER,*)'MOM_J_CMF_NONMONO_V1: insertion cap',
	1   KN_NM_MAX_INS,' hit in',N_CAP,' of',ND_SM-1,
	2   ' intervals (ND_NM=',ND_NM,')'
	  END IF
!
	  ALLOCATE ( R_NM(ND_NM),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( R_PNT_NM(ND_NM),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate R_NM in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
!
! Build the expanded R grid.  Use ABS(V_SM(I)-V_SM(I+1)) for the
! sub-point count so that velocity reversals (V(I) < V(I+1)) produce
! the same sub-grid resolution as the corresponding expansion.
!
          K=1
	  R_NM(1)=R_SM(1)
          R_PNT_NM(1)=1
	  DO I=1,ND_SM-1
            IT1=INT( ABS(V_SM(I)-V_SM(I+1))/T2 )
	    IF(IT1 .GT. KN_NM_MAX_INS)IT1=KN_NM_MAX_INS
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)THEN
              DELTA_R=(R_SM(I+1)-R_SM(I))/(IT1+1)
              DO J=1,IT1
                K=K+1
                R_NM(K)=R_NM(K-1)+DELTA_R
                R_PNT_NM(K)=I
              END DO
            END IF
            K=K+1
            R_NM(K)=R_SM(I+1)
            R_PNT_NM(K)=I
	  END DO
!
	  ALLOCATE ( LOG_R_SM_NM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_SM_SV(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_TRUE_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ESEC_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF_NM(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate COEF in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
!
	  ALLOCATE ( V_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( SIGMA_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( SIGMA_TRUE_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( ETA_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( ESEC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( CHI_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( RSQJNU_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( RSQHNU_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( RSQJNU_PREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( RSQHNU_PREV_NM(ND_NM),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate V_BLOCK in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  RSQJNU_NM(1:ND_NM)=0.0_LDP; RSQJNU_PREV_NM(1:ND_NM)=0.0_LDP
	  RSQHNU_NM(1:ND_NM)=0.0_LDP; RSQHNU_PREV_NM(1:ND_NM)=0.0_LDP
!
	  ALLOCATE ( CON_GAM_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( CON_GAMH_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( AV_SIGMA_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( AV_SIGMA_TRUE_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( CHI_HAT_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAM_MAG_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAM_SIGNED_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( Q_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAMH_MAG_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAMH_SIGNED_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( WH_MAG_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( DTAU_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( DTAUONQ_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( TA_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( TB_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( TC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( QF_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( XM_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( XM_SAVE_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( XM_FB_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( XM_BAD_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( SOURCE_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( VB_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( VC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( HU_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( HL_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( HS_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( COH_VEC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAM_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAMH_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAM_FAC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAMH_FAC_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAM_FAC_PREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( GAMH_FAC_PREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( W_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( WPREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( PSI_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( PSIPREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( EPS_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( EPS_PREV_NM(ND_NM),STAT=IOS )
	  ALLOCATE ( J_INDX_NM(ND_SM) )
	  ALLOCATE ( H_INDX_NM(ND_SM) )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate work arrays in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
!
! J_NEXT storage for Gauss-Seidel convergence (Knop Paper II).
! Allocated to (ND_NM,NCF_IN): the compressing-zone upwind field is the
! NEXT frequency's R^2*J from the previous great iteration, so the store
! must be frequency-indexed (the old 1-D J_NEXT_DEPTH held the previous
! -- bluer -- frequency of the current sweep, i.e. the wrong field).
! Converges over great iterations.
!
	  ALLOCATE ( J_NEXT_FRQ(ND_NM,NCF_IN),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate J_NEXT_FRQ in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  J_NEXT_FRQ(1:ND_NM,1:NCF_IN)=0.0_LDP
	  J_NEXT_VALID=.FALSE.
!
! H_NEXT_FRQ: lagged R^2*H partner of J_NEXT_FRQ for the
! compressing-midpoint H/EPS upwind reads.  Always allocated
! alongside J_NEXT_FRQ; its validity shares the J_NEXT_VALID flag.
!
	  ALLOCATE ( H_NEXT_FRQ(ND_NM,NCF_IN),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate H_NEXT_FRQ in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  H_NEXT_FRQ(1:ND_NM,1:NCF_IN)=0.0_LDP
!
! Quarantine persistence ages: 2-byte integers are enough (a cell
! quarantined for >32767 consecutive sweeps has long stopped
! mattering).  Allocated with J_NEXT_FRQ so a grid rebuild resets
! the ages together with the store they describe.
!
	  ALLOCATE ( MOM_QUAR_AGE(ND_NM,NCF_IN),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate MOM_QUAR_AGE in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  MOM_QUAR_AGE(1:ND_NM,1:NCF_IN)=0_2
!
	  ALLOCATE ( MOM_OWN_CMP(ND_NM,NCF_IN),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate MOM_OWN_CMP in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  MOM_OWN_CMP(1:ND_NM,1:NCF_IN)=.FALSE.
	  ALLOCATE ( MOM_OWN_CMPH(ND_NM,NCF_IN),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate MOM_OWN_CMPH in MOM_J_CMF_NONMONO_V1'
	     STOP
	  END IF
	  MOM_OWN_CMPH(1:ND_NM,1:NCF_IN)=.FALSE.
!
! Fresh lagged stores mean the call-level AA history describes a
! different field; a rebuilt grid invalidates it either way.
!
	  NHIST_OUT=-1
!
! Per-frequency Knop refresh store (Oct-2026): captured inputs for
! MOM_KNOP_REFRESH_V1's alternating-direction re-solve.  Allocated
! unconditionally; the refresh self-guards when no compressing
! depths exist, but the capture is cheap and keeps the logic simple.
!
	  ALLOCATE (KN_M_CHISM(ND_SM,NCF_IN),
	1     KN_M_ETASM(ND_SM,NCF_IN),KN_M_ESECSM(ND_SM,NCF_IN),
	2     KN_MKONJ(ND_NM,NCF_IN),KN_MNMJ(ND_NM,NCF_IN),
	3     KN_MNMH(ND_NM,NCF_IN),
	4     KN_MHBC(NCF_IN),KN_MNBC(NCF_IN),KN_MINHBC(NCF_IN),
	5     KN_MJPIB(NCF_IN),KN_MJMIB(NCF_IN),KN_MHIB(NCF_IN),
	6     KN_MDBB(NCF_IN),KN_MIC(NCF_IN),
	7     KN_MDLN(NCF_IN),KN_MNUNXT(NCF_IN),KN_MFREQ(NCF_IN),
	8     KN_M_RSM(ND_SM),KN_M_VSM(ND_SM),KN_M_SIGSM(ND_SM),
	9     KN_M_SIGTSM(ND_SM),KN_M_VDOPV(ND_SM),
	1     KN_M_OUTBC(NCF_IN),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate KN_M refresh store'
	     STOP
	  END IF
	  KN_M_CHISM=0.0; KN_M_ETASM=0.0; KN_M_ESECSM=0.0
	  KN_MKONJ=0.0; KN_MNMJ=0.0; KN_MNMH=0.0
	  KN_MHBC=0.0_LDP; KN_MNBC=0.0_LDP; KN_MINHBC=0.0_LDP
	  KN_MJPIB=0.0_LDP; KN_MJMIB=0.0_LDP; KN_MHIB=0.0_LDP
	  KN_MDBB=0.0_LDP; KN_MIC=0.0_LDP
	  KN_MDLN=1.0_LDP; KN_MNUNXT=-1.0_LDP; KN_MFREQ=0.0_LDP
	  KN_M_OUTBC=1
	  KN_M_RSM(1:ND_SM)=R_SM(1:ND_SM)
	  KN_M_VSM(1:ND_SM)=V_SM(1:ND_SM)
	  KN_M_SIGSM(1:ND_SM)=SIGMA_SM(1:ND_SM)
	  KN_M_SIGTSM(1:ND_SM)=SIGMA_TRUE_SM(1:ND_SM)
	  KN_M_VDOPV(1:ND_SM)=VDOP_VEC(1:ND_SM)
	  KN_M_IBSTAB=IB_STAB_FACTOR
	  KN_M_VDOPF=VDOP_FRAC
	  KN_M_NC=NC; KN_M_NP=NP; KN_M_NDSM=ND_SM
	  KN_M_NDNM=ND_NM; KN_M_NCF=NCF_IN
	  KN_M_NTYPE=N_TYPE; KN_M_HCHK=H_CHK_OPTION
	  KN_M_METHOD=METHOD; KN_M_INBM=INNER_BND_METH
	  KN_M_COH=COHERENT
	  KN_M_READY=.FALSE.
	  KN_M_ANYCMP=.FALSE.
!
! Grid-rebuild remap: if the old expanded grid carried a valid GS
! store, interpolate R^2*J from R_NM_OLD onto the new R_NM for every
! frequency.  This preserves the compressing-zone upwind field across
! grid rebuilds; without it each rebuild zeroed J_NEXT_FRQ and the GS
! iteration never converged.
!
	  IF(HAVE_JNEXT_OLD)THEN
! The RDINR drift moves the model R each iteration, so the rebuilt
! expanded grid can sit just outside the old grid's range.  MON_INTERP
! STOPs on out-of-range queries; clamp query points to the old
! endpoints (nearest-value extension -- a tiny boundary extrapolation,
! vastly better than zeroing the whole GS store).
	    ALLOCATE(R_QR_CL(ND_NM),STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Unable to allocate R_QR_CL in MOM_J_CMF_NONMONO_V1'
	      STOP
	    END IF
	    R_HI=MAX(R_NM_OLD(1),R_NM_OLD(ND_NM_OLD))
	    R_LO=MIN(R_NM_OLD(1),R_NM_OLD(ND_NM_OLD))
	    DO I=1,ND_NM
	      R_QR_CL(I)=MIN(MAX(R_NM(I),R_LO),R_HI)
	    END DO
	    DO K=1,MIN(NCF_IN,SIZE(JNEXT_OLD,2))
	      CALL MON_INTERP(J_NEXT_FRQ(1,K),ND_NM,IONE,R_QR_CL,ND_NM,
	1          JNEXT_OLD(1,K),ND_NM_OLD,R_NM_OLD,ND_NM_OLD)
	    END DO
	    IF(ALLOCATED(HNEXT_OLD))THEN
	      DO K=1,MIN(NCF_IN,SIZE(HNEXT_OLD,2))
	        CALL MON_INTERP(H_NEXT_FRQ(1,K),ND_NM,IONE,R_QR_CL,
	1            ND_NM,HNEXT_OLD(1,K),ND_NM_OLD,R_NM_OLD,ND_NM_OLD)
	      END DO
	      DEALLOCATE(HNEXT_OLD)
	    END IF
	    DEALLOCATE(R_QR_CL)
	    J_NEXT_VALID=.TRUE.
	    DEALLOCATE(JNEXT_OLD,R_NM_OLD)
	    HAVE_JNEXT_OLD=.FALSE.
	    LUER=ERROR_LU()
	    WRITE(LUER,'(A)')' MOM_J_CMF_NONMONO_V1: remapped J_NEXT_FRQ'
	1      //' onto rebuilt expanded grid'
	  END IF
!
! GS restart (16-Sep-2026): save NCF_IN for the checkpoint dump
! (SCR_RITE_V2 cannot see it), then attempt to restore a validated
! complete J_NEXT_FRQ generation from MOMJNXT_GS.  On any validation
! failure the array is left untouched, so a remapped store survives.
! The disk restore is skipped when the in-process remap already
! supplied a valid store (fresher than the checkpoint).
!
	  MOM_GS_NCF=NCF_IN
	  IF(.NOT. J_NEXT_VALID)THEN
! RESTORE_MOM_JNEXT_GS checks the checkpoint's physical-grid fingerprint
! (ND_F, R/V endpoints written by DUMP_MOM_JNEXT_GS against the SM grid);
! the expanded ND_NM belongs only to the J_NEXT_FRQ dimension check.
	    CALL RESTORE_MOM_JNEXT_GS(R_SM,V_SM,ND_SM)
	  END IF
!
! Interpolate V, SIGMA, and SIGMA_TRUE onto the expanded grid.
! SIGMA_TRUE is interpolated (not clamped) so the non-monotonic
! velocity gradient is preserved on the sub-grid.
!
	  LOG_R_SM_NM(1:ND_SM)=LOG(R_SM(1:ND_SM))
	  V_SM_SV(1:ND_SM)=V_SM(1:ND_SM)
          CALL MON_INT_FUNS_V2(V_COEF_NM,V_SM,LOG_R_SM_NM,ND_SM)
          CALL MON_INT_FUNS_V2(SIGMA_COEF_NM,SIGMA_SM,LOG_R_SM_NM,ND_SM)
          CALL MON_INT_FUNS_V2(SIGMA_TRUE_COEF_NM,SIGMA_TRUE_SM,
	1                       LOG_R_SM_NM,ND_SM)
          DO I=1,ND_NM
            K=R_PNT_NM(I)
            T1=LOG(R_NM(I)/R_SM(K))
            V_NM(I)=((V_COEF_NM(K,1)*T1+V_COEF_NM(K,2))*T1+
	1            V_COEF_NM(K,3))*T1+V_COEF_NM(K,4)
            SIGMA_NM(I)=((SIGMA_COEF_NM(K,1)*T1+SIGMA_COEF_NM(K,2))*T1+
	1            SIGMA_COEF_NM(K,3))*T1+SIGMA_COEF_NM(K,4)
            SIGMA_TRUE_NM(I)=((SIGMA_TRUE_COEF_NM(K,1)*T1+
	1            SIGMA_TRUE_COEF_NM(K,2))*T1+
	1            SIGMA_TRUE_COEF_NM(K,3))*T1+SIGMA_TRUE_COEF_NM(K,4)
          END DO
!
! Build J_INDX and H_INDX mappings (same as V11).
!
	  J_INDX_NM(1:ND_SM)=0
	  H_INDX_NM(1:ND_SM)=0
	  K=1
	  DO I=1,ND_SM
!
! edited by Devin, 29-Sep-2026: bound the march in-loop, as in
! FG_J_CMF_NONMONO_V1.  The expanded R_NM grid is built with inserted
! points for the non-monotonic field, so a depth that brackets no
! interval would otherwise walk K past ND_NM and read out of bounds
! in an unchecked build.
!
	    DO WHILE(J_INDX_NM(I) .EQ. 0)
	      IF(R_SM(I) .LE. R_NM(K) .AND. R_SM(I) .GE. R_NM(K+1))THEN
	        IF( (R_NM(K)-R_SM(I)) .LT. (R_SM(I)-R_NM(K+1)) )THEN
	          J_INDX_NM(I)=K
	        ELSE
	          J_INDX_NM(I)=K+1
	        END IF
	      ELSE
	        K=K+1
	        IF(K+1 .GT. ND_NM)THEN
	          LUER=ERROR_LU()
	          WRITE(LUER,*)'Error setting J_INDX_NM in'//
	1          ' MOM_J_CMF_NONMONO_V1 -- no bracket for depth'
	          WRITE(LUER,*)'Depth=',I,'R_SM(I)=',R_SM(I)
	          WRITE(LUER,*)'R_NM(1)=',R_NM(1),'R_NM(ND_NM)=',R_NM(ND_NM)
	          STOP
	        END IF
	      END IF
	    END DO
	  END DO
!
	  K=1
	  DO I=1,ND_SM-1
	    T1=0.5_LDP*(R_SM(I)+R_SM(I+1))
	    DO WHILE(H_INDX_NM(I) .EQ. 0)
	      IF(T1 .LT. R_NM(K) .AND. T1 .GT. R_NM(K+1))THEN
	        H_INDX_NM(I)=K
	      ELSE
	        K=K+1
	        IF(K+1 .GT. ND_NM)THEN
	          LUER=ERROR_LU()
	          WRITE(LUER,*)'Error setting H_INDX_NM in'//
	1          ' MOM_J_CMF_NONMONO_V1 -- no bracket for midpoint'
	          WRITE(LUER,*)'Depth=',I,'T1=',T1
	          WRITE(LUER,*)'R_NM(1)=',R_NM(1),'R_NM(ND_NM)=',R_NM(ND_NM)
	          STOP
	        END IF
	      END IF
	    END DO
	  END DO
!
	  IF(N_TYPE .NE. 'G_ONLY' .AND. N_TYPE .NE. 'N_ON_J')THEN
	    IF(ND_NM .NE. ND_SM)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Warning in MOM_J_CMF_NONMONO_V11'
	      WRITE(LUER,*)'N_TYPE=MIXED with inserted points may be inaccurate'
	    END IF
	  END IF
!
	  FIRST_TIME_NM=.FALSE.
	END IF
!
! Reset the per-sweep accumulators at the start of a sweep, BEFORE
! this frequency's solve contributes to them.  The reset used to sit
! in the J_NEXT update section below, so frequency 1's contributions
! were erased every sweep -- and on descending refresh passes (which
! end on FREQ_INDX=1) the whole pass's physics stats were wiped, which
! is why DUMP_MOM_JNEXT_GS printed sentinel/zero blocks on ~half of
! the great iterations.  Skipped during refresh: the driver resets
! the GS metrics per frequency and the physics stats per pass.
!
	IF(FREQ_INDX .EQ. 1 .AND. .NOT. MOM_KN_REFRESHING)THEN
!
! Compression-dominance resolution for MOM_RF_DIR.  The lagged
! coupling direction is a property of the velocity field, not of any
! single frequency solve, so the flags are set from the previous
! sweep's accumulated row census before it is cleared.  On the first
! sweep the counters are zero and both flags stay false -- the
! alternating order and unsuppressed writes run until the census
! exists.
!
	  IF(MOM_RF_DIR .EQ. 'DESC')THEN
	    MOM_DOM_CMP=.TRUE.;  MOM_DOM_EXP=.FALSE.
	  ELSE IF(MOM_RF_DIR .EQ. 'ASC ')THEN
	    MOM_DOM_CMP=.FALSE.;  MOM_DOM_EXP=.TRUE.
	  ELSE IF(MOM_RF_DIR .EQ. 'ALT ')THEN
	    MOM_DOM_CMP=.FALSE.;  MOM_DOM_EXP=.FALSE.
	  ELSE IF(MOM_SW_NROW .GT. 0_8)THEN
	    MOM_SW_FCMP=DBLE(MOM_SW_NCROW)/DBLE(MOM_SW_NROW)
	    MOM_DOM_CMP=(10_8*MOM_SW_NCROW .GE. 9_8*MOM_SW_NROW)
	    MOM_DOM_EXP=(10_8*MOM_SW_NCROW .LE. MOM_SW_NROW)
	  END IF
	  MOM_SW_NCROW=0;  MOM_SW_NROW=0
	  MOM_EPS_CAP=0;  MOM_EPS_CAP_SW=0
	  MOM_SW_NSUPP=0
	  MOM_GS_MAXCHG=0.0_LDP
	  MOM_INEXT_CAPPED=0
	  MOM_NEGJ_LRG=0
	  MOM_GS_NBAD5=0;  MOM_GS_NBAD50=0
	  MOM_GS_NHUGE=0
	  MOM_GS_NQUAR=0
	  MOM_GS_NPQUAR=0
	  MOM_GS_NTOT=INT(ND_NM,8)*INT(NCF_IN,8)
	  MOM_GS_ARGI=0;  MOM_GS_ARGL=0
	  MOM_SW_NTOT=0
	  MOM_RF_RAN=.FALSE.
	  MOM_RF_CONV=.FALSE.
	  MOM_RF_DCHG=0.0_LDP
	  MOM_RF_N5=0;  MOM_RF_N50=0
	  MOM_NONMM_CNT=0
	  MOM_EPS_FWD=0
	  MOM_EPS_RAT=0.0_LDP
	  MOM_BND_HBC_FL=0
	  MOM_N_DGTSV=0
	  MOM_N_NEGJ_T=0
	  MOM_N_NONF=0
	  MOM_N_RUNAWAY=0
	  MOM_TB_NONPOS=0
	  MOM_NEGJ_DUMP=0
	  MOM_MIN_DR=1.0E+30_LDP
	  MOM_MIN_DR_I=0;  MOM_MIN_DR_L=0
	  MOM_MIN_DRA=1.0E+30_LDP
	  MOM_MIN_DRA_I=0;  MOM_MIN_DRA_L=0
	  MOM_EQ_RES=0.0_LDP
	  MOM_EQ_RES_I=0;  MOM_EQ_RES_L=0
	  MOM_MAX_COH=0.0_LDP
	  MOM_MAX_GAM=0.0_LDP
	  MOM_MAX_SIGT=0.0_LDP
	  MOM_LAGB_CNT=0
	  MOM_LAGB_MIN=1.0_LDP
	END IF
!
!
!
! Get Eddington factors on the revised grid.
!
! MOM_KN_REFRESHING: the refresh driver pre-loads the per-frequency
! Eddington state from the captured store; GET_MOMS_NON_REL would
! rebuild it from EDDFACTOR records that only exist for the last
! frequency of the main sweep.
!
	IF(.NOT. MOM_KN_REFRESHING)THEN
!
!	  Request the real midpoint g = N/H under N_ON_J: the compressing
!	  branches of this routine key on GAMH_SIGNED = sign(1+sigma*g),
!	  which is dead at g=0 (the historical fill).  Opt-in only; stock
!	  callers keep the zero fill.
!
	  GMNR_NMJ_G=.TRUE.
	  CALL GET_MOMS_NON_REL(R_NM, V_NM, FREQ, N_TYPE, ND_NM)
!
! Sanitize the interpolated Eddington factors before they enter
! GAM/EPS/W.  GET_MOMS_NON_REL can return NaN or unphysical values
! at edges where the store is degenerate; a NaN K_ON_J or a
! non-finite NMID_ON_HMID then poisons the tridiagonal coefficients.
! Bounds: f=K/J in [1E-4,1] (isotropic 1/3 to free-streaming 1);
! NMID_ON_J keeps its sign (it is a signed flux-like ratio) but
! non-finite is zeroed.  NMID_ON_HMID is g=N/H -- genuinely signed
! and unbounded now that N_ON_J fills it for the direction split;
! only non-finite and |g|>1E12 are regularized, since this routine
! always takes |1+sigma*g| rather than 1+sigma*g itself.
!
	  DO I=1,ND_NM
	    IF(.NOT. (ABS(K_ON_J(I)) .LE. HUGE(K_ON_J(I))))
	1      K_ON_J(I)=1.0_LDP/3.0_LDP
	    K_ON_J(I)=MIN(1.0_LDP,MAX(1.0E-4_LDP,K_ON_J(I)))
	  END DO
	  DO I=1,ND_NM-1
	    IF(.NOT. (ABS(NMID_ON_J(I)) .LE. HUGE(NMID_ON_J(I))))
	1      NMID_ON_J(I)=0.0_LDP
	    IF(.NOT. (ABS(NMID_ON_HMID(I)) .LE. HUGE(NMID_ON_HMID(I))))
	1      NMID_ON_HMID(I)=0.0_LDP
	    NMID_ON_HMID(I)=SIGN(MIN(ABS(NMID_ON_HMID(I)),
	1       1.0E+12_LDP),NMID_ON_HMID(I))
	  END DO
!
! MOM_FAC_RLX blend: this call's sanitized factors are folded into
! the per-frequency relaxed store and the store is fed back, so the
! coefficients below see a damped factor trajectory.  First call per
! slot seeds directly.  Runs on the main sweep only -- refresh
! re-calls read the captured KN_M factor state (already relaxed) and
! skip the GET_MOMS_NON_REL block entirely.
!
	  IF(MOM_FAC_RLX .LT. 1.0_LDP .AND. FREQ_INDX .GE. 1 .AND.
	1     FREQ_INDX .LE. NCF_IN)THEN
	    IF(.NOT. ALLOCATED(FACR_KONJ) .OR.
	1       SIZE(FACR_KONJ,1) .NE. ND_NM .OR.
	2       SIZE(FACR_KONJ,2) .NE. NCF_IN)THEN
	      IF(ALLOCATED(FACR_KONJ))DEALLOCATE(FACR_KONJ,FACR_NMJ,
	1        FACR_NMH,FACR_SEED)
	      ALLOCATE(FACR_KONJ(ND_NM,NCF_IN),FACR_NMJ(ND_NM,NCF_IN),
	1        FACR_NMH(ND_NM,NCF_IN),FACR_SEED(NCF_IN),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'MOM_FAC_RLX: factor store alloc failed;',
	1        ' relax disabled'
	        MOM_FAC_RLX=1.0_LDP
	      ELSE
	        FACR_SEED=.FALSE.
	      END IF
	    END IF
	    IF(MOM_FAC_RLX .LT. 1.0_LDP)THEN
	      IF(FACR_SEED(FREQ_INDX))THEN
	        DO I=1,ND_NM
	          FACR_KONJ(I,FREQ_INDX)=FACR_KONJ(I,FREQ_INDX)+
	1            MOM_FAC_RLX*(K_ON_J(I)-FACR_KONJ(I,FREQ_INDX))
	          K_ON_J(I)=FACR_KONJ(I,FREQ_INDX)
	        END DO
	        DO I=1,ND_NM-1
	          FACR_NMJ(I,FREQ_INDX)=FACR_NMJ(I,FREQ_INDX)+
	1            MOM_FAC_RLX*(NMID_ON_J(I)-FACR_NMJ(I,FREQ_INDX))
	          NMID_ON_J(I)=FACR_NMJ(I,FREQ_INDX)
	          FACR_NMH(I,FREQ_INDX)=FACR_NMH(I,FREQ_INDX)+
	1            MOM_FAC_RLX*(NMID_ON_HMID(I)-FACR_NMH(I,FREQ_INDX))
	          NMID_ON_HMID(I)=FACR_NMH(I,FREQ_INDX)
	        END DO
	      ELSE
	        FACR_KONJ(1:ND_NM,FREQ_INDX)=K_ON_J(1:ND_NM)
	        FACR_NMJ(1:ND_NM,FREQ_INDX)=NMID_ON_J(1:ND_NM)
	        FACR_NMH(1:ND_NM,FREQ_INDX)=NMID_ON_HMID(1:ND_NM)
	        FACR_SEED(FREQ_INDX)=.TRUE.
	      END IF
	    END IF
	  END IF
	END IF
!
	IF(ND_NM .GT. ND_SM)THEN
!
! Interpolate opacities onto the expanded grid (log interpolation).
!
! Guard the LOG arguments: a zero, negative or NaN opacity at one SM
! node previously made LOG(0)=-Inf, which MON_INT_FUNS_V2 then spread
! over the whole interval's coefficients (NaN poisoning the solve).
! MERGE with a >0 test is used instead of MAX because Fortran does
! not define MAX on NaN operands; the predicate rejects both
! non-positive and non-finite inputs.  The floor maps back to ~0
! after EXP, which is the intended limit of a vanishing opacity.
!
	  DO I=1,ND_SM
	    TA_NM(I)=LOG(MERGE(CHI_SM(I),1.0E-290_LDP,
	1        CHI_SM(I) .GT. 0.0_LDP))
	  END DO
	  CALL MON_INT_FUNS_V2(CHI_COEF_NM,TA_NM,LOG_R_SM_NM,ND_SM)
	  DO I=1,ND_SM
	    TA_NM(I)=LOG(MERGE(ESEC_SM(I),1.0E-290_LDP,
	1        ESEC_SM(I) .GT. 0.0_LDP))
	  END DO
	  CALL MON_INT_FUNS_V2(ESEC_COEF_NM,TA_NM,LOG_R_SM_NM,ND_SM)
	  DO I=1,ND_SM
	    TA_NM(I)=LOG(MERGE(ETA_SM(I),1.0E-290_LDP,
	1        ETA_SM(I) .GT. 0.0_LDP))
	  END DO
	  CALL MON_INT_FUNS_V2(ETA_COEF_NM,TA_NM,LOG_R_SM_NM,ND_SM)
	  DO I=1,ND_NM
	    K=R_PNT_NM(I)
	    T1=LOG(R_NM(I)/R_SM(K))
	    T2=((CHI_COEF_NM(K,1)*T1+CHI_COEF_NM(K,2))*T1+
	1      CHI_COEF_NM(K,3))*T1+CHI_COEF_NM(K,4)
	    CHI_NM(I)=EXP(T2)
	    T2=((ESEC_COEF_NM(K,1)*T1+ESEC_COEF_NM(K,2))*T1+
	1      ESEC_COEF_NM(K,3))*T1+ESEC_COEF_NM(K,4)
	    ESEC_NM(I)=EXP(T2)
	    T2=((ETA_COEF_NM(K,1)*T1+ETA_COEF_NM(K,2))*T1+
	1      ETA_COEF_NM(K,3))*T1+ETA_COEF_NM(K,4)
	    ETA_NM(I)=EXP(T2)
	  END DO
!
	ELSE
	  ESEC_NM(1:ND_NM)=ESEC_SM(1:ND_SM)
	  CHI_NM(1:ND_NM)=CHI_SM(1:ND_SM)
	  ETA_NM(1:ND_NM)=ETA_SM(1:ND_SM)
	END IF
!
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0_LDP
	  END DO
	  MOM_ERR_CNT=0
	END IF
!
!*****************************************************************************
!
! Compute source function and coherent scattering vector.
!
! CHI_NM can underflow to exact 0.0 from the floored LOG
! interpolation above; a bare division would then give 0/0 NaN in
! SOURCE_NM or Inf in COH_VEC_NM.
!
	DO I=1,ND_NM
	  SOURCE_NM(I)=ETA_NM(I)/MAX(CHI_NM(I),1.0E-290_LDP)
	END DO
	IF(COHERENT)THEN
	  DO I=1,ND_NM
	    COH_VEC_NM(I)=ESEC_NM(I)/MAX(CHI_NM(I),1.0E-290_LDP)
	  END DO
	ELSE
	  DO I=1,ND_NM
	    COH_VEC_NM(I)=0.0_LDP
	  END DO
	END IF
!
! Compute Q factors from F (Eddington factor K/J) and optical depth.
! NB: QF_NM is the geometric Q factor (not to be confused with Q_NM
! the velocity coupling coefficient).
!
	CALL QFROMF(K_ON_J,QF_NM,R_NM,TA_NM,TB_NM,ND_NM)
!
! Bound QF_NM to a finite positive range.  The formal-solver Eddington
! factors can be degenerate at boundary depths (F ~ 0 or F <= 0 where J
! is tiny), in which case the (3F-1)/(R F) integral in QFROMF underflows
! to 0 or overflows to +Inf.  QF_NM=0 divides DTAUONQ by zero; QF_NM=Inf
! overflows the K_ON_J*QF_NM coefficient products (TC(1)=-Inf).  Both
! ends are the saturated asymptotic limit of the moment equations
! (decoupled / optically thick), so clamping preserves the limiting
! behaviour while keeping the tridiagonal coefficients finite.
!
	N_QBND=0
	DO I=1,ND_NM
	  IF(.NOT. (QF_NM(I) .GT. 1.0E-20_LDP))THEN
	    QF_NM(I)=1.0E-20_LDP
	    N_QBND=N_QBND+1
	  ELSE IF(.NOT. (QF_NM(I) .LT. 1.0E+30_LDP))THEN
	    QF_NM(I)=1.0E+30_LDP
	    N_QBND=N_QBND+1
	  END IF
	END DO
	IF(N_QBND .GT. 0 .AND. N_QPRT .LT. 15)THEN
	  N_QPRT=N_QPRT+1
	  WRITE(6,'(A,ES14.6,A,I5,A)')' MOM_J_CMF_NONMONO_V1: '//
	1     'QF_NM bounded at',FREQ,' :',N_QBND,' depths'
	END IF
!
	DO I=1,ND_NM
	  TA_NM(I)=CHI_NM(I)*QF_NM(I)
	END DO
	CALL DERIVCHI(TB_NM,TA_NM,R_NM,ND_NM,METHOD)
!
!*****************************************************************************
! Non-monotonic generalized opacity (Knop et al. 2009, Paper I Eq 7).
!
! For the moment equations, the comoving-frame frequency gradient is:
!   GAM_MOM = (V/R) * c^{-1} * (1 + SIGMA_TRUE * K_ON_J)
!
! where K_ON_J is the Eddington factor (angle-averaged mu^2).
! This is the angle-integrated analogue of the ray-based GAM in
! FG_J_CMF_NONMONO_V1.
!
! The generalized opacity is:
!   chi_hat = chi + |GAM_MOM| / dLOG_NU   (always positive)
!
! The signed coupling coefficient is:
!   Q = (GAM_MOM / dLOG_NU) / chi_hat     (signed)
!
! The modified source is:
!   S' = SOURCE + Q * J_upwind
!
! where J_upwind = J_PREV when GAM_MOM >= 0 (expanding) and
!   J_NEXT when GAM_MOM < 0 (compressing).
!
! This split guarantees:
!   1. chi_hat > 0  =>  DTAU > 0  (well-conditioned tridiagonal)
!   2. PSI = DTAUONQ * |GAM_MOM|/dLOG_NU  (always positive diagonal)
!   3. W   = |GAMH_MOM|/dLOG_NU          (always positive H-coupling)
!   4. The signed Q*J_upwind enters the RHS only (cannot destabilize)
!
!*****************************************************************************
!
	IF(NEW_FREQ)THEN
	  K_ON_J_PREV(1:ND_NM)=K_ON_J_SAVE(1:ND_NM)
	  NMID_ON_HMID_PREV(1:ND_NM)=NMID_ON_HMID_SAVE(1:ND_NM)
	  NMID_ON_J_PREV(1:ND_NM)=NMID_ON_J_SAVE(1:ND_NM)
	  RSQJNU_PREV_NM(1:ND_NM)=RSQJNU_NM(1:ND_NM)
	  RSQHNU_PREV_NM(1:ND_NM)=RSQHNU_NM(1:ND_NM)
	  HBC_PREV=HBC_SAVE;  IN_HBC_PREV=IN_HBC_SAVE
	  NBC_PREV=NBC_SAVE
	END IF
!
	IF(INIT)THEN
!
! At the bluest frequency no lagged frequency fields exist, so the
! explicit (PREV) coupling arrays are zeroed.  The CURRENT-frequency
! coefficients (CHI_HAT, GAM, W, EPS, Q) are physical properties of
! this frequency and are computed unconditionally below -- zeroing
! them here (as the V11-derived structure did) leaves CHI_HAT=0 on
! the first call, giving DTAU=0, Inf matrix rows and a non-finite
! tridiagonal solution at every depth, and on later calls silently
! reuses the reddest frequency's CHI_HAT -- the wrong opacity.
!
	  DO I=1,ND_NM
	    WPREV_NM(I)=0.0_LDP
	    PSIPREV_NM(I)=0.0_LDP
	    RSQJNU_PREV_NM(I)=0.0_LDP
	    RSQHNU_PREV_NM(I)=0.0_LDP
	    EPS_PREV_NM(I)=0.0_LDP
	  END DO
	  HBC_PREV=0.0_LDP;  IN_HBC_PREV=0.0_LDP; NBC_PREV=0.0_LDP
!
! The lagged Eddington arrays are never filled at FREQ_INDX=1
! (NEW_FREQ is false), but PSIPREV is assembled unconditionally below
! and compressing points still need the redder-neighbour J_next
! coupling, so set PREV = current values for the boundary frequency.
!
	  K_ON_J_PREV(1:ND_NM)=K_ON_J(1:ND_NM)
	  NMID_ON_J_PREV(1:ND_NM)=NMID_ON_J(1:ND_NM)
	  NMID_ON_HMID_PREV(1:ND_NM)=NMID_ON_HMID(1:ND_NM)
	END IF
!
! Current-frequency coefficients are required at every call,
! including the bluest frequency (INIT).
!
!
! Compute GAM_MOM, GAMH_MOM and their magnitudes.
! Use SIGMA_TRUE (unclamped) for the velocity gradient.
!
! GAM_MOM = (V/R) * c^{-1} * (1 + SIGMA_TRUE * K_ON_J)
! GAMH_MOM = (V_mid/R_mid) * c^{-1} * (1 + AV_SIGMA_TRUE * NMID_ON_HMID)
!
! The (1 + SIGMA * angle_factor) is the FULL comoving frequency gradient.
! The "1" is the V/R divergence (aberration); SIGMA * angle_factor is
! the velocity gradient correction.  Both must be present.
!
! c^{-1} = 3.33564E-06 (converts km/s/cm to dimensionless)
!
	  HAS_COMPRESSING=.FALSE.
	  DO I=1,ND_NM-1
	    CON_GAMH_NM(I)=2.0_LDP*3.33564E-06_LDP*
	1       (V_NM(I)+V_NM(I+1))/(R_NM(I)+R_NM(I+1))
	    AV_SIGMA_NM(I)=0.5_LDP*(SIGMA_NM(I)+SIGMA_NM(I+1))
	    AV_SIGMA_TRUE_NM(I)=0.5_LDP*(SIGMA_TRUE_NM(I)+SIGMA_TRUE_NM(I+1))
	    CON_GAM_NM(I)=3.33564E-06_LDP*V_NM(I)/R_NM(I)
!
! Signed GAMH_MOM (mid-point comoving frequency gradient).
! Uses SIGMA_TRUE (unclamped) and NMID_ON_HMID (angle factor).
! The (1 + ...) is the full velocity-gradient factor.
!
! edited by Devin, 28-Sep-2026: regularize the factor 1+sigma*K at the
! zero-GAM resonance so the sign split, chi_hat and the explicit/implicit
! weights stay away from zero.  The sign is preserved and the magnitude is
! limited from below; large-|GAM| branches are unchanged.
!
	    T1=1.0_LDP+AV_SIGMA_TRUE_NM(I)*NMID_ON_HMID(I)
	    GAMH_FAC_NM(I)=SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)
	    IF(T1 .EQ. 0.0_LDP)GAMH_FAC_NM(I)=0.0_LDP
	    GAM_SIGNED_NM(I)=CON_GAMH_NM(I)*GAMH_FAC_NM(I)
	    GAMH_SIGNED_NM(I)=GAM_SIGNED_NM(I)
	    GAM_MAG_NM(I)=ABS(GAM_SIGNED_NM(I))
	    IF(GAM_SIGNED_NM(I) .LT. 0.0_LDP)HAS_COMPRESSING=.TRUE.
	    IF(MOM_GAM_PIN)THEN
!
! Pin the lagged factor to the current-frequency evaluation: an
! independent 1+sigma*g_prev can carry the opposite sign or a larger
! magnitude near the zero-GAM resonance, which lets the lagged H
! weight exceed the implicit (1+W) normalization.  [MOM_GAM_PIN]
!
	      GAMH_FAC_PREV_NM(I)=GAMH_FAC_NM(I)
	    ELSE
	      T1=1.0_LDP+AV_SIGMA_TRUE_NM(I)*NMID_ON_HMID_PREV(I)
	      GAMH_FAC_PREV_NM(I)=SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)
	      IF(T1 .EQ. 0.0_LDP)GAMH_FAC_PREV_NM(I)=0.0_LDP
	    END IF
	  END DO
	  GAMH_FAC_NM(ND_NM)=1.0_LDP
	  GAMH_FAC_PREV_NM(ND_NM)=1.0_LDP
	  CON_GAM_NM(ND_NM)=3.33564E-06_LDP*V_NM(ND_NM)/R_NM(ND_NM)
!
! Frequency-match bookkeeping (VADAT DO_NONMONO_FREQ_MATCH, default
! F; set by SET_MOM_FREQ_MATCH).  At compressing cells the lagged
! frequency upwind is rematched from the adjacent bin k+1 to the
! bin k+n containing the characteristic foot, n = |GAM|*ds/DLN in
! grid bins, and the frequency derivative spans ln(nu_k/nu_{k+n}).
! Once the per-cell comoving drift exceeds one bin the +1 stencil
! extrapolates the derivative -- the oscillatory-eigenvalue regime
! of Ceulemans et al. (2025, JQSRT 343).  NM_J/NM_H hold the
! matched offsets at depth/mid points; DLN_J/DLN_H the matched log
! spans.  The arrays are filled every call so the coefficient loops
! read a single spacing table; with the flag off (or at expanding
! points) the fills reproduce the pre-match code exactly --
! compressing cells read k+1 over the forward spacing ln(nu_k/
! nu_{k+1}), expanding cells keep dLOG_NU.
!
	IOS=0
	IF(.NOT. ALLOCATED(MOM_NM_J))THEN
	  ALLOCATE(MOM_NM_J(ND_NM),MOM_NM_H(ND_NM),
	1      MOM_DLN_J(ND_NM),MOM_DLN_H(ND_NM),STAT=IOS)
	ELSE IF(SIZE(MOM_NM_J) .NE. ND_NM)THEN
	  DEALLOCATE(MOM_NM_J,MOM_NM_H,MOM_DLN_J,MOM_DLN_H)
	  ALLOCATE(MOM_NM_J(ND_NM),MOM_NM_H(ND_NM),
	1      MOM_DLN_J(ND_NM),MOM_DLN_H(ND_NM),STAT=IOS)
	END IF
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)
	1    'MOM_J_CMF_NONMONO_V1: freq-match alloc failed'
	  STOP
	END IF
	MOM_NM_J(1:ND_NM)=1
	MOM_DLN_J(1:ND_NM)=dLOG_NU
	MOM_NM_H(1:ND_NM)=1
	MOM_DLN_H(1:ND_NM)=dLOG_NU
!
! The frequency-grid cache fills one entry per call (plus the
! next-frequency entry supplied by SET_MOM_NU_NEXT), so a matched
! span uses true log(nu) once a full sweep has run; before that the
! uniform-spacing estimate n*DLN applies.
!
	IF(.NOT. ALLOCATED(MOM_NU_GRID))THEN
	  ALLOCATE(MOM_NU_GRID(NCF_IN),STAT=IOS)
	  IF(IOS .EQ. 0)MOM_NU_GRID=0.0_LDP
	ELSE IF(SIZE(MOM_NU_GRID) .NE. NCF_IN)THEN
	  DEALLOCATE(MOM_NU_GRID)
	  ALLOCATE(MOM_NU_GRID(NCF_IN),STAT=IOS)
	  IF(IOS .EQ. 0)MOM_NU_GRID=0.0_LDP
	END IF
	IF(ALLOCATED(MOM_NU_GRID))THEN
	  MOM_NU_GRID(FREQ_INDX)=FREQ
	  IF(FREQ_INDX .LT. NCF_IN .AND. MOM_NU_NEXT .GT. 0.0_LDP)
	1    MOM_NU_GRID(FREQ_INDX+1)=MOM_NU_NEXT
	END IF
!
	DO I=1,ND_NM-1
	  IF(GAMH_SIGNED_NM(I) .LT. 0.0_LDP .AND.
	1     MOM_NU_NEXT .GT. 0.0_LDP .AND.
	2     MOM_NU_NEXT .LT. FREQ)THEN
	    MOM_DLN_H(I)=LOG(FREQ/MOM_NU_NEXT)
	    IF(MOM_FREQ_MATCH .AND. NCF_IN-FREQ_INDX .GE. 1)THEN
	      DELTA_R=R_NM(I)-R_NM(I+1)
	      T4=ABS(GAMH_SIGNED_NM(I))*ABS(DELTA_R)/MOM_DLN_H(I)
	      IF(T4 .GE. DBLE(NCF_IN-FREQ_INDX))THEN
	        N_CAP=NCF_IN-FREQ_INDX
	      ELSE
!
! CEILING keeps the characteristic foot inside the matched
! stencil; NINT could round down and extrapolate up to half
! a bin past the foot.
!
	        N_CAP=MAX(1,CEILING(T4))
	      END IF
	      MOM_NM_H(I)=N_CAP
	      MOM_DLN_H(I)=N_CAP*MOM_DLN_H(I)
	      IF(ALLOCATED(MOM_NU_GRID))THEN
	        IF(MOM_NU_GRID(FREQ_INDX+N_CAP) .GT. 0.0_LDP .AND.
	1           MOM_NU_GRID(FREQ_INDX+N_CAP) .LT. FREQ)
	1          MOM_DLN_H(I)=LOG(FREQ/MOM_NU_GRID(FREQ_INDX+N_CAP))
	      END IF
	    END IF
	  END IF
	END DO
!
! Compute chi_hat at depth points (for DTAU and GAM); midpoint
! values are the average of the depth-point values below.  The
! earlier midpoint loop was removed: it read GAM_MAG_NM before the
! depth loop populated it and its result was never used.
! Knop et al. 2009 Paper I Eq 7: |GAM| in opacity keeps it positive.
!
! chi_hat(I) = chi(I) + |CON_GAM(I) * (1 + SIGMA_TRUE(I) * K_ON_J(I))|/dLOG_NU
!
	  DO I=1,ND_NM
!
! edited by Devin, 28-Sep-2026: same zero-GAM regularization as the
! midpoint factor above, so the depth-point sign split, chi_hat, Q and
! the explicit source all use a consistent non-zero factor.
!
	    T1=1.0_LDP+SIGMA_TRUE_NM(I)*K_ON_J(I)
	    GAM_FAC_NM(I)=SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)
	    IF(T1 .EQ. 0.0_LDP)GAM_FAC_NM(I)=0.0_LDP
	    GAM_SIGNED_NM(I)=CON_GAM_NM(I)*GAM_FAC_NM(I)
	    GAM_MAG_NM(I)=ABS(GAM_SIGNED_NM(I))
	    GAM_NM(I)=GAM_SIGNED_NM(I)
! The compressing branch keys off this depth-point sign (Q_NM<0), so the
! refresh gate must see it too -- the midpoint test above can stay
! positive when NMID_ON_HMID sits at its 0.01 floor.
	    IF(GAM_SIGNED_NM(I) .LT. 0.0_LDP)HAS_COMPRESSING=.TRUE.
	    IF(MOM_GAM_PIN)THEN
!
! Same pin as the midpoint factor: PSIPREV = DTAUONQ*|GAM*GAM_FAC|
! now equals PSI, so the expanding explicit upwind weight cannot
! exceed the implicit diagonal (the compressing branch already uses
! PSI as its coupling coefficient).  [MOM_GAM_PIN]
!
	      GAM_FAC_PREV_NM(I)=GAM_FAC_NM(I)
	    ELSE
	      T1=1.0_LDP+SIGMA_TRUE_NM(I)*K_ON_J_PREV(I)
	      GAM_FAC_PREV_NM(I)=SIGN(MAX(ABS(T1),MOM_GAM_FAC_MIN),T1)
	      IF(T1 .EQ. 0.0_LDP)GAM_FAC_PREV_NM(I)=0.0_LDP
	    END IF
	    IF(GAM_SIGNED_NM(I) .LT. 0.0_LDP .AND.
	1       MOM_NU_NEXT .GT. 0.0_LDP .AND.
	2       MOM_NU_NEXT .LT. FREQ)THEN
	      MOM_DLN_J(I)=LOG(FREQ/MOM_NU_NEXT)
	      IF(MOM_FREQ_MATCH .AND. NCF_IN-FREQ_INDX .GE. 1)THEN
	        IF(I .EQ. 1)THEN
	          DELTA_R=R_NM(1)-R_NM(2)
	        ELSE IF(I .EQ. ND_NM)THEN
	          DELTA_R=R_NM(ND_NM-1)-R_NM(ND_NM)
	        ELSE
	          DELTA_R=0.5_LDP*(R_NM(I-1)-R_NM(I+1))
	        END IF
	        T4=GAM_MAG_NM(I)*ABS(DELTA_R)/MOM_DLN_J(I)
	        IF(T4 .GE. DBLE(NCF_IN-FREQ_INDX))THEN
	          N_CAP=NCF_IN-FREQ_INDX
	        ELSE
	          N_CAP=MAX(1,CEILING(T4))
	        END IF
	        MOM_NM_J(I)=N_CAP
	        MOM_DLN_J(I)=N_CAP*MOM_DLN_J(I)
	        IF(ALLOCATED(MOM_NU_GRID))THEN
	          IF(MOM_NU_GRID(FREQ_INDX+N_CAP) .GT. 0.0_LDP .AND.
	1             MOM_NU_GRID(FREQ_INDX+N_CAP) .LT. FREQ)
	1            MOM_DLN_J(I)=
	1              LOG(FREQ/MOM_NU_GRID(FREQ_INDX+N_CAP))
	        END IF
	      END IF
	    END IF
	    CHI_HAT_NM(I)=CHI_NM(I)+GAM_MAG_NM(I)/
	1       MAX(MOM_DLN_J(I),MOM_DLN_MIN)
	    IF(CHI_HAT_NM(I) .LT. 1.0E-30_LDP)CHI_HAT_NM(I)=1.0E-30_LDP
	  END DO
!
! Compute GAMH (dimensionless, like V11 but with chi_hat instead of chi).
! GAMH = CON_GAMH / dLOG_NU / (CHI_HAT + CHI_HAT+1)
!
	  DO I=1,ND_NM-1
	    DLNM=MAX(MOM_DLN_H(I),MOM_DLN_MIN)
	    GAMH_NM(I)=CON_GAMH_NM(I)/DLNM/
	1       (CHI_HAT_NM(I)+CHI_HAT_NM(I+1))
	  END DO
!
! Compute W and WPREV using |1 + AV_SIGMA_TRUE * NMID_ON_HMID| (always >= 0).
! V11: W = GAMH * (1 + AV_SIGMA * NMID_ON_HMID)
! NONMONO: W = GAMH * |1 + AV_SIGMA_TRUE * NMID_ON_HMID| (always positive)
! The signed part is moved to the source via Q_H * J_upwind.
!
	  DO I=1,ND_NM-1
	    W_NM(I)=GAMH_NM(I)*ABS(GAMH_FAC_NM(I))
	    WPREV_NM(I)=GAMH_NM(I)*ABS(GAMH_FAC_PREV_NM(I))
	  END DO
!
! Compute GAM (dimensionless, like V11 but with chi_hat instead of chi).
! GAM = CON_GAM / CHI_HAT / dLOG_NU
!
	  DO I=1,ND_NM
	    DLNM=MAX(MOM_DLN_J(I),MOM_DLN_MIN)
! GAM_NM held the physical signed value above; it is deliberately
! overwritten here with the dimensionless matrix coefficient.  The
! signed copy lives in GAM_SIGNED_NM -- read that for the physics.
	    GAM_NM(I)=CON_GAM_NM(I)/CHI_HAT_NM(I)/DLNM
	  END DO
!
! Compute Q_J at depth points (signed, for J-equation source).
! Q_J = GAM * (1 + SIGMA_TRUE * K_ON_J) (signed)
! This is the Knop coupling coefficient for the J-equation.
!
	  DO I=1,ND_NM
	    Q_NM(I)=GAM_NM(I)*GAM_FAC_NM(I)
	  END DO
!
! Coupling-direction census for MOM_RF_DIR='AUTO' and the
! compression-dominant write suppression: count rows whose upwind
! reads the lagged redder bin.  Refresh re-calls are excluded --
! they re-solve the same sweep and would double-count it.
!
	IF(.NOT. MOM_KN_REFRESHING)THEN
	  MOM_SW_NCROW=MOM_SW_NCROW+
	1    INT(COUNT(Q_NM(1:ND_NM) .LT. 0.0_LDP),8)
	  MOM_SW_NROW=MOM_SW_NROW+INT(ND_NM,8)
	END IF
!
! EPS for N_ON_J formulation.
! V11: EPS = GAMH * AV_SIGMA * NMID_ON_J / (1 + W)
! NONMONO: same formula but with GAMH using chi_hat and AV_SIGMA_TRUE.
! EPS uses just the gradient part (AV_SIGMA_TRUE * NMID_ON_J), not
! the full (1 + AV_SIGMA_TRUE * NMID_ON_HMID), matching V11.
! The sign of EPS_NM itself now selects the depth-upwind node in the
! matrix assembly; EPS_SGN_NM (the abandoned forward-difference
! variant) was removed 23-Sep-2026.
!
! EPS_PREV denominator (verified against stock V11, which uses the
! identical /(1+W) form): eliminating H from the moment pair gives a
! single row whose coefficients are all normalized by the CURRENT
! frequency's (1+W) factor -- that factor comes from the current
! row's own diagonal elimination, not from the lagged field.  So
! both EPS_NM and EPS_PREV_NM carry 1/(1+W_NM); WPREV_NM enters the
! row only as the numerator of the lagged-H coefficient
! HS_NM = WPREV_NM/(1+W_NM) below.  Under MOM_GAM_PIN the lagged
! angle factor in EPS_PREV_NM is NMID_ON_J itself (pinned); with
! the flag off it is NMID_ON_J_PREV -- in both cases the lagged
! content is the field the term multiplies (RSQJNU_PREV_NM).  GAMH
! and AV_SIGMA_TRUE are not lagged because they carry no frequency
! index in this discretization.
!
!
	  IF(N_TYPE .EQ. 'G_ONLY')THEN
	    DO I=1,ND_NM-1
	      EPS_NM(I)=0.0_LDP
	      EPS_PREV_NM(I)=0.0_LDP
	    END DO
	  ELSE
	    DO I=1,ND_NM-1
	      EPS_NM(I)=GAMH_NM(I)*AV_SIGMA_TRUE_NM(I)*
	1         NMID_ON_J(I)/(1.0_LDP+W_NM(I))
	      IF(MOM_GAM_PIN)THEN
!
! Pin the lagged angle factor to the current-frequency NMID_ON_J:
! same lagged/current mismatch class as GAM_FAC_PREV -- an
! independent sigma*n_prev evaluation lets the explicit midpoint
! transfer outweigh the implicit EPS in TA/TB/TC.  With the pin the
! EPS_PREV<0 fixup below is already satisfied.  [MOM_GAM_PIN]
!
	        EPS_PREV_NM(I)=EPS_NM(I)
	      ELSE
	        EPS_PREV_NM(I)=GAMH_NM(I)*AV_SIGMA_TRUE_NM(I)*
	1           NMID_ON_J_PREV(I)/(1.0_LDP+W_NM(I))
	      END IF
	    END DO
	  END IF
!
! Compute DTAU using chi_hat (generalized opacity, always > 0).
! This replaces CHI*Q with chi_hat*Q in the optical depth.
!
	DO I=1,ND_NM
	  TA_NM(I)=CHI_HAT_NM(I)*QF_NM(I)
	END DO
	CALL DERIVCHI(TB_NM,TA_NM,R_NM,ND_NM,METHOD)
	CALL NORDTAU(DTAU_NM,TA_NM,R_NM,R_NM,TB_NM,ND_NM)
!
! Compute PSI and PSIPREV.
! V11: PSI = DTAUONQ * GAM * (1 + SIGMA * K_ON_J)
! NONMONO: PSI = DTAUONQ * |GAM * (1 + SIGMA_TRUE * K_ON_J)| (always >= 0)
! The chi in GAM (CON_GAM/CHI_HAT/dLOG_NU) cancels with chi_hat in DTAUONQ,
! giving PSI = dr * |CON_GAM * (1 + SIGMA_TRUE * K_ON_J)| / dLOG_NU,
! same as V11 where chi cancels.
!
	DO I=2,ND_NM-1
	  DTAUONQ_NM(I)=0.5_LDP*(DTAU_NM(I)+DTAU_NM(I-1))/QF_NM(I)
!
! PSI = DTAUONQ * |GAM * (1 + SIGMA_TRUE * K_ON_J)| (always >= 0, implicit)
! PSIPREV = DTAUONQ * |GAM * GAM_FAC_PREV| (explicit); under
!   MOM_GAM_PIN GAM_FAC_PREV = GAM_FAC so PSIPREV = PSI, with the
!   flag off it uses K_ON_J_PREV (independent lagged evaluation).
!
! V11: PSI = DTAUONQ * GAM * (1 + SIGMA * K_ON_J) (always positive in V11
!   because SIGMA is clamped to >= -0.999 and K_ON_J ~ 1/3)
! NONMONO: PSI uses |...| for positivity; PSIPREV is signed to handle
!   compressing zones via the upwind source term.
!
! PSIPREV is now kept positive (|GAM*GAM_FAC_PREV|): Paper I's
! construction puts a strictly POSITIVE weight on the upwind field
! (-xi*a*p^+/- over chi_hat > 0).  With unclamped SIGMA_TRUE the
! lagged factor GAM_FAC_PREV can carry the opposite sign to the
! current row's Q_NM near the zero-GAM resonance, which injected a
! negative upwind source onto expanding rows that the Knop split
! can never produce.  The magnitude convention (lagged Eddington
! factor) is unchanged; only the sign is regularized.
!
	  PSI_NM(I)=DTAUONQ_NM(I)*ABS(Q_NM(I))
	  PSIPREV_NM(I)=DTAUONQ_NM(I)*ABS(GAM_NM(I)*GAM_FAC_PREV_NM(I))
	END DO
!
! Compute H-equation vectors using |W| (always positive).
!
	DO I=1,ND_NM-1
	  HU_NM(I)=K_ON_J(I+1)*QF_NM(I+1)/(1.0_LDP+W_NM(I))/DTAU_NM(I)
	  HL_NM(I)=K_ON_J(I)*QF_NM(I)/(1.0_LDP+W_NM(I))/DTAU_NM(I)
	  HS_NM(I)=WPREV_NM(I)/(1.0_LDP+W_NM(I))
	END DO
!
! Compressing-midpoint fixup (the advected-frequency direction
! reverses when GAMH_SIGNED < 0).  The explicit lagged H term then
! reads the REDDER neighbour field, so its coefficient must carry
! the current-frequency weights: HS becomes W/(1+W) (the same
! convention PSI_NM uses for the J coupling).  The sigma*n*J piece
! is keyed on its own sign below.  Overwriting here makes the RHS
! assembly, the outer boundary row, the negative-J rebuild and the
! H reconstruction take the consistent compressing form without
! touching the expanding branch.
!
	DO I=1,ND_NM-1
	  IF(GAMH_SIGNED_NM(I) .LT. 0.0_LDP)
	1    HS_NM(I)=W_NM(I)/(1.0_LDP+W_NM(I))
!
!	  The sigma*n*J piece carries its own direction (flux split): its
!	  frequency upwind is sign(EPS) = sign(sigma*n), not the H-part
!	  characteristic sign(1+sigma*g).  Where it reads the lagged
!	  redder field the lagged coefficient takes the current weights,
!	  matching the PSI_NM convention for the J coupling.
!
	  IF(EPS_NM(I) .LT. 0.0_LDP)EPS_PREV_NM(I)=EPS_NM(I)
	END DO
!
!*****************************************************************************
! Build the tridiagonal system.
!
! The key difference from V11: the velocity-coupling is split into
!   implicit (|GAM*(1+SIGMA*K_ON_J)| in PSI, always >= 0) and
!   explicit (PSIPREV * J_upwind in RHS, signed).
!
! J_upwind = RSQJNU_PREV when Q >= 0 (expanding, standard CMFGEN)
! J_upwind = J_NEXT_FRQ(:,FREQ_INDX+1) when Q < 0 (compressing, Knop
! Gauss-Seidel: next/redder frequency lagged from previous iteration)
!
!*****************************************************************************
!
! Absorption fraction floor (24-Sep-2026).  ESEC and CHI are
! interpolated onto the expanded grid independently, so at near
! pure-scattering cells the ratio COH_VEC=ESEC/CHI can exceed 1.
! (1-COH)*DTAUONQ then subtracts from the diagonal; verified at a
! negJ cell (TB=0.004 vs |TA|+|TC|=0.055) that it can cancel ~99%
! of PSI+HL+HU, leaving a near-singular row whose solution flips
! to a large negative J that then drags neighbours down.
! Absorption cannot be negative: floor the fraction at pure
! scattering.  COH_VEC itself is left raw so the negJ record still
! shows the ESEC>CHI crossing.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND_NM-1
	    TA_NM(I)=-HL_NM(I-1)
	    TC_NM(I)=-HU_NM(I)
	    TB_NM(I)=DTAUONQ_NM(I)*
	1           MAX(MOM_ABS_FLR,1.0_LDP-COH_VEC_NM(I)) + PSI_NM(I) +
	2           HL_NM(I) + HU_NM(I-1)
	    VB_NM(I)=-HS_NM(I-1)
	    VC_NM(I)=HS_NM(I)
	    XM_NM(I)=DTAUONQ_NM(I)*SOURCE_NM(I)*R_NM(I)*R_NM(I)
	  END DO
	ELSE
	  DO I=2,ND_NM-1
! Depth-upwind sigma*NMID_ON_J*J midpoint term (23-Sep-2026).  The
! centered (J_m+J_{m+1}) average puts a positive off-diagonal on one
! of the two rows for EITHER sign of EPS once |EPS|>HU or HL, so the
! tridiagonal is not an M-matrix and the solution flips sign at the
! shock cells.  The midpoint value is instead taken at the upwind
! node: EPS>=0 reads the outer node (row m gets +E*J_m on the
! diagonal, row m+1 gets -E*J_m on the subdiagonal), EPS<0 reads the
! inner node (row m gets E*J_{m+1}=-|E| on the superdiagonal, row
! m+1 gets -E*J_{m+1}=+|E| on the diagonal).  Sign-definite for both
! branches -- a first-order flux evaluation where the term is
! advection-dominated; identical to the centered form in smooth
! zones where EPS~0 anyway.  The factor 2 keeps the V11 sum
! convention: Jsum=(J_m+J_{m+1})=2*J_mid, so J_mid=J_up gives
! 2*E*J_up -- matching the converged magnitude, not half of it.
	    TA_NM(I)=-HL_NM(I-1)
	    TC_NM(I)=-HU_NM(I)
	    TB_NM(I)=DTAUONQ_NM(I)*
	1           MAX(MOM_ABS_FLR,1.0_LDP-COH_VEC_NM(I)) + PSI_NM(I) +
	2           HL_NM(I) + HU_NM(I-1)
	    IF(EPS_NM(I-1) .GE. 0.0_LDP)THEN
	      TA_NM(I)=TA_NM(I)-2.0_LDP*EPS_NM(I-1)
	    ELSE
	      TB_NM(I)=TB_NM(I)-2.0_LDP*EPS_NM(I-1)
	    END IF
	    IF(EPS_NM(I) .GE. 0.0_LDP)THEN
	      TB_NM(I)=TB_NM(I)+2.0_LDP*EPS_NM(I)
	    ELSE
	      TC_NM(I)=TC_NM(I)+2.0_LDP*EPS_NM(I)
	    END IF
	    VB_NM(I)=-HS_NM(I-1)
	    VC_NM(I)=HS_NM(I)
	    XM_NM(I)=DTAUONQ_NM(I)*SOURCE_NM(I)*R_NM(I)*R_NM(I)
	  END DO
	END IF
!
! MOM_H_IMPL: pay for the lagged upwind-H divergence with its own
! coupling magnitude on the diagonal.  The row keeps the same RHS
! (the signed -HS(I-1)*H_up(I-1)+HS(I)*H_up(I) stays in XLAG_NM),
! but the solve sees the H channel's weight as diagonal mass: a
! clipped shock row stops oscillating between its clipped and
! unclipped states, and the MOM_LAG_BND excess below grows by the
! same |VB|+|VC| so less of the lagged source is scaled away.
! Boundary row 1 assembles its own bound and is unchanged.
!
	IF(MOM_H_IMPL)THEN
	  DO I=2,ND_NM-1
	    TB_NM(I)=TB_NM(I)+ABS(VB_NM(I))+ABS(VC_NM(I))
	  END DO
	END IF
!
! Sweep diagnostics: rows where the EPS coupling overwhelms the
! spatial superdiagonal (TC>0 -> non-M-matrix; zero under the
! depth-upwind form since the EPS<0 superdiagonal is -|E|), the
! count of EPS<0 midpoints (inner-node upwind), and the worst
! |EPS|/HU ratio.
!
	DO I=2,ND_NM-1
	  IF(TC_NM(I) .GT. 0.0_LDP)MOM_NONMM_CNT=MOM_NONMM_CNT+1
	END DO
	DO I=1,ND_NM-1
	  IF(EPS_NM(I) .LT. 0.0_LDP)MOM_EPS_FWD=MOM_EPS_FWD+1
	  IF(HU_NM(I) .GT. 0.0_LDP)MOM_EPS_RAT=MAX(MOM_EPS_RAT,
	1    ABS(EPS_NM(I))/HU_NM(I))
	END DO
!
! The velocity-coupling source term is handled by PSIPREV * J_upwind
! in the explicit frequency-coupling section below (matching V11's
! structure).  No separate Q*J_upwind term is needed.
!
! Inner boundary condition (depth ND).  All four method rows divide
! by DTAU_NM(ND_NM-1); a nonpositive increment would flip the row
! sign or produce Inf, so reject the frequency outright.
!
	IF(.NOT.(DTAU_NM(ND_NM-1) .GT. 0.0_LDP))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid inner-boundary DTAU:',DTAU_NM(ND_NM-1),
	1     ' FREQ_INDX=',FREQ_INDX
	  STOP
	END IF
	TA_NM(ND_NM)=-K_ON_J(ND_NM-1)*QF_NM(ND_NM-1)/DTAU_NM(ND_NM-1)
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  TB_NM(ND_NM)=K_ON_J(ND_NM)/DTAU_NM(ND_NM-1)
	  XM_NM(ND_NM)=DBB*R_NM(ND_NM)*R_NM(ND_NM)/3.0_LDP/CHI_HAT_NM(ND_NM)
	ELSE IF(INNER_BND_METH .EQ. 'PNT_SRCE')THEN
	  TB_NM(ND_NM)=K_ON_J(ND_NM)/DTAU_NM(ND_NM-1)
	  XM_NM(ND_NM)=HNU_AT_IB*R_NM(ND_NM)*R_NM(ND_NM)
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  TB_NM(ND_NM)=K_ON_J(ND_NM)/DTAU_NM(ND_NM-1)
	  XM_NM(ND_NM)=IB_STAB_FACTOR*TB_NM(ND_NM)*R_NM(ND_NM)*R_NM(ND_NM)*
	1      (JPLUS_IB+JMIN_IB)
	  TB_NM(ND_NM)=(1.0_LDP+IB_STAB_FACTOR)*TB_NM(ND_NM)
	ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	  TB_NM(ND_NM)=K_ON_J(ND_NM)/DTAU_NM(ND_NM-1)+IN_HBC
	  XM_NM(ND_NM)=R_NM(ND_NM)*R_NM(ND_NM)*IC*(0.25_LDP+0.5_LDP*IN_HBC)
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid INNER_BND_METH:', INNER_BND_METH
	  STOP
	END IF
	TC_NM(ND_NM)=0.0_LDP
	VB_NM(ND_NM)=0.0_LDP
	VC_NM(ND_NM)=0.0_LDP
	PSIPREV_NM(ND_NM)=0.0_LDP
!
! Add the explicit frequency-coupling terms to the RHS.
! V11: XM(I) = ... + PSIPREV(I) * RSQJNU_PREV(I)
! NONMONO: XM(I) = ... + coeff * J_upwind
!   where J_upwind = RSQJNU_PREV if Q >= 0 (expanding, coeff=PSIPREV)
!                   J_NEXT_FRQ(:,k+1) if Q < 0 (compressing, Knop GS,
!                   coeff=PSI_NM = the diagonal |Q_k| weight, matching
!                   FG's compressing source S + |Q_k|*I_NEXT).
! Equal explicit/implicit weights keep the lagged J_NEXT map
! non-expansive; a mismatched (K_ON_J_PREV-weighted) coefficient can
! exceed the diagonal weight and amplify J every great iteration.
!
! During refresh re-calls (MOM_KN_REFRESHING) the compressing upwind
! column is taken from the pass-start snapshot JNEXT_RF_READ, not the
! live J_NEXT_FRQ: on a descending refresh pass slot k+1 holds a
! this-pass write, turning the read into a GS cascade that amplified
! residuals to ~1E10 (job 11887664).  The snapshot restores the
! in-sweep pure-Jacobi semantics (previous generation's field).
!
	IF(.NOT. ALLOCATED(JUP_NM))THEN
	  ALLOCATE(JUP_NM(ND_NM))
	ELSE IF(SIZE(JUP_NM) .NE. ND_NM)THEN
	  DEALLOCATE(JUP_NM)
	  ALLOCATE(JUP_NM(ND_NM))
	END IF
	JUP_NM(1:ND_NM)=0.0_LDP
!
! LAGB workspaces (MOM_LAG_BND): XLAG_NM accumulates the lagged part
! of each row's RHS; LAG_SCL_NM holds the bound factor so the
! negative-J rebuild scales the same terms.
!
	IF(.NOT. ALLOCATED(XLAG_NM))THEN
	  ALLOCATE(XLAG_NM(ND_NM),LAG_SCL_NM(ND_NM))
	ELSE IF(SIZE(XLAG_NM) .NE. ND_NM)THEN
	  DEALLOCATE(XLAG_NM,LAG_SCL_NM)
	  ALLOCATE(XLAG_NM(ND_NM),LAG_SCL_NM(ND_NM))
	END IF
	XLAG_NM(1:ND_NM)=0.0_LDP
	LAG_SCL_NM(1:ND_NM)=1.0_LDP
	IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
!
! Per-point matched read: the compressing upwind is the bin
! FREQ_INDX+MOM_NM_J(I) holding the characteristic foot (offset 1
! unless DO_NONMONO_FREQ_MATCH rematched it).  The refresh snapshot
! vs live-GS choice below is unchanged; on descending passes slot
! k+n has already been visited, so the live read propagates the
! coupling just as the k+1 read did.
!
	  IF(MOM_KN_REFRESHING .AND. ALLOCATED(JNEXT_RF_READ) .AND.
	1     .NOT. (KN_M_DESC .AND.
	2     (MOM_RF_GS_DESC .OR. MOM_DOM_CMP)))THEN
!
! Frozen pass-start snapshot.  On ascending passes slot L+1 has not
! been visited yet, so the snapshot is the only consistent read; on
! descending passes it is the Jacobi fallback used when the live-GS
! option is off.
!
	    DO I=1,ND_NM
	      IF(MOM_NM_J(I) .GE. 1 .AND.
	1         FREQ_INDX+MOM_NM_J(I) .LE. NCF_IN)
	1        JUP_NM(I)=JNEXT_RF_READ(I,FREQ_INDX+MOM_NM_J(I))
	    END DO
	  ELSE
!
! Descending pass with MOM_RF_GS_DESC: slot L+1 was solved this pass
! -- a live Gauss-Seidel read that propagates the compressing
! (red->blue) coupling across the whole band in one pass.
!
	    DO I=1,ND_NM
	      IF(MOM_NM_J(I) .GE. 1 .AND.
	1         FREQ_INDX+MOM_NM_J(I) .LE. NCF_IN)
	1        JUP_NM(I)=J_NEXT_FRQ(I,FREQ_INDX+MOM_NM_J(I))
	    END DO
	  END IF
	END IF
!
! HUP_NM: the upwind R^2*H for the explicit HS*H_up midpoint terms.
! Expanding midpoints read the bluer chain (RSQHNU_PREV_NM);
! compressing midpoints read the lagged redder store H_NEXT_FRQ with
! the same refresh-snapshot semantics as JUP_NM.
!
	IF(.NOT. ALLOCATED(HUP_NM))THEN
	  ALLOCATE(HUP_NM(ND_NM))
	ELSE IF(SIZE(HUP_NM) .NE. ND_NM)THEN
	  DEALLOCATE(HUP_NM)
	  ALLOCATE(HUP_NM(ND_NM))
	END IF
	HUP_NM(1:ND_NM)=0.0_LDP
	DO I=1,ND_NM-1
	  IF(GAMH_SIGNED_NM(I) .GE. 0.0_LDP)THEN
	    HUP_NM(I)=RSQHNU_PREV_NM(I)
	  ELSE IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN .AND.
	1     FREQ_INDX+MOM_NM_H(I) .LE. NCF_IN)THEN
	    IF(MOM_KN_REFRESHING .AND. ALLOCATED(HNEXT_RF_READ) .AND.
	1       .NOT. (KN_M_DESC .AND.
	2       (MOM_RF_GS_DESC .OR. MOM_DOM_CMP)))THEN
	      T1=HNEXT_RF_READ(I,FREQ_INDX+MOM_NM_H(I))
	    ELSE
	      T1=H_NEXT_FRQ(I,FREQ_INDX+MOM_NM_H(I))
	    END IF
!
! A non-finite lagged R^2*H (stale restart entry or a failed earlier
! sweep) would enter the tridiagonal RHS unchecked; J_NEXT_FRQ gets
! the same finite test.  Zeroing drops this cell's advected H
! contribution rather than propagating NaN through the solve.
!
	    IF(.NOT. (ABS(T1) .LE. HUGE(T1)))THEN
	      HUP_NM(I)=0.0_LDP
	      MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+1
	    ELSE
	      HUP_NM(I)=T1
	    END IF
	  END IF
	END DO
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND_NM-1
	    T3=VB_NM(I)*HUP_NM(I-1) + VC_NM(I)*HUP_NM(I)
	    IF(Q_NM(I) .GE. 0.0_LDP)THEN
	      T3=T3 + PSIPREV_NM(I)*RSQJNU_PREV_NM(I)
	    ELSE
! GAM<0: the explicit term on the LHS is GAM*J_{nu+1}/dlnnu, i.e.
! -|coeff|*J_next; moved to the RHS it enters with a POSITIVE sign,
! +PSI_NM*J_next (matching FG's SOURCE_PRIME = S + |Q_k|*I_NEXT).
	      IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)T3=T3 +
	1         PSI_NM(I)*JUP_NM(I)
	    END IF
	    XM_NM(I)=XM_NM(I)+T3
	    XLAG_NM(I)=XLAG_NM(I)+T3
	  END DO
	ELSE
	  DO I=2,ND_NM-1
	    T3=VB_NM(I)*HUP_NM(I-1) + VC_NM(I)*HUP_NM(I)
	    IF(Q_NM(I) .GE. 0.0_LDP)THEN
	      T3=T3 + PSIPREV_NM(I)*RSQJNU_PREV_NM(I)
	    ELSE
	      IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)T3=T3 +
	1         PSI_NM(I)*JUP_NM(I)
	    END IF
	    XM_NM(I)=XM_NM(I)+T3
	    XLAG_NM(I)=XLAG_NM(I)+T3
	  END DO
!
! Lagged EPS*(sigma n J) midpoint terms.  The advected quantity uses
! the term's own frequency direction (flux split): EPS >= 0 reads the
! bluer RSQJNU_PREV_NM, EPS < 0 the lagged redder JUP_NM.  (Under
! N_ON_J before the NMID_ON_HMID fill, GAMH_SIGNED was always >= 0,
! so this term was read downwind wherever sigma*n < 0.)  The
! coefficient sign still selects the depth-upwind node.  Scattered
! per midpoint: midpoint I contributes +2E*J_up to row I and
! -2E*J_up to row I+1.  Row 1's share is added by the boundary block
! below and row ND is replaced by the inner BC, so neither is
! written here.
!
	IF(.NOT. ALLOCATED(EPS_SCL_NM))THEN
	  ALLOCATE(EPS_SCL_NM(ND_NM))
	ELSE IF(SIZE(EPS_SCL_NM) .NE. ND_NM)THEN
	  DEALLOCATE(EPS_SCL_NM)
	  ALLOCATE(EPS_SCL_NM(ND_NM))
	END IF
	EPS_SCL_NM(1:ND_NM)=1.0_LDP
	  DO I=1,ND_NM-1
	    IF(EPS_PREV_NM(I) .GE. 0.0_LDP)THEN
	      K=I
	    ELSE
	      K=I+1
	    END IF
	    IF(EPS_PREV_NM(I) .LT. 0.0_LDP)THEN
	      T4=JUP_NM(K)
	    ELSE
	      T4=RSQJNU_PREV_NM(K)
	    END IF
	    T3=2.0_LDP*EPS_PREV_NM(I)*T4
	    T5=T3
	    IF(MOM_EPS_LIM)THEN
!
!	      Donor-cell positivity limiter on the lagged transfer (a
!	      signed +T5 to row I, -T5 to row I+1).  The debit may not
!	      exceed the donor row's current positive RHS balance: this
!	      term is lagged-frequency auxiliary data, and at
!	      advection-dominated cells (|EPS|/HU ~ 20) the uncapped
!	      debit drove interior rows negative, returning negative J
!	      from a well-conditioned M-matrix -- the persistent negJ
!	      population of job 12045050.  The recipient is credited
!	      only the amount actually debited, so the pair stays
!	      conservative; boundary rows (not written here) are
!	      unaffected.  [MOM_EPS_LIM, default T]
!
	      IF(T5 .LT. 0.0_LDP .AND. I .GT. 1)THEN
	        T4=MAX(XM_NM(I),0.0_LDP)
	        IF(-T5 .GT. T4)THEN
	          T5=-T4
	          MOM_EPS_CAP=MOM_EPS_CAP+1
	          MOM_EPS_CAP_SW=MOM_EPS_CAP_SW+1
	        END IF
	      ELSE IF(T5 .GT. 0.0_LDP .AND. I .LT. ND_NM-1)THEN
	        T4=MAX(XM_NM(I+1),0.0_LDP)
	        IF(T5 .GT. T4)THEN
	          T5=T4
	          MOM_EPS_CAP=MOM_EPS_CAP+1
	          MOM_EPS_CAP_SW=MOM_EPS_CAP_SW+1
	        END IF
	      END IF
!
!	      The implicit half of the same transfer sits uncapped in
!	      TA/TB/TC: capping only this side turned a frequency
!	      derivative into a net source/sink wherever |EPS| >> HU
!	      (millions of events per sweep on m15).  Record the cap
!	      factor so the implicit half is scaled to match below --
!	      the pair stays a (smaller, consistent) derivative.
!
	      IF(T3 .NE. 0.0_LDP)EPS_SCL_NM(I)=T5/T3
	    END IF
	    IF(I .GT. 1)THEN
	      XM_NM(I)=XM_NM(I)+T5
	      XLAG_NM(I)=XLAG_NM(I)+T5
	    END IF
	    IF(I .LT. ND_NM-1)THEN
	      XM_NM(I+1)=XM_NM(I+1)-T5
	      XLAG_NM(I+1)=XLAG_NM(I+1)-T5
	    END IF
	  END DO
!
! Apply the cap factor to the implicit half: rows I and I+1 carry
! (EPS_SCL-1)*2*EPS on the same depth-upwind stencil the assembly
! used, then both EPS arrays are rescaled so the boundary row, the
! negJ rebuild and the H reconstruction see the capped coefficient.
!
	  DO I=1,ND_NM-1
	    IF(EPS_SCL_NM(I) .NE. 1.0_LDP)THEN
	      T3=2.0_LDP*EPS_NM(I)*(EPS_SCL_NM(I)-1.0_LDP)
	      IF(I .GT. 1)THEN
	        IF(EPS_NM(I) .GE. 0.0_LDP)THEN
	          TB_NM(I)=TB_NM(I)+T3
	        ELSE
	          TC_NM(I)=TC_NM(I)+T3
	        END IF
	      END IF
	      IF(I .LT. ND_NM-1)THEN
	        IF(EPS_NM(I) .GE. 0.0_LDP)THEN
	          TA_NM(I+1)=TA_NM(I+1)-T3
	        ELSE
	          TB_NM(I+1)=TB_NM(I+1)-T3
	        END IF
	      END IF
	      EPS_NM(I)=EPS_NM(I)*EPS_SCL_NM(I)
	      EPS_PREV_NM(I)=EPS_PREV_NM(I)*EPS_SCL_NM(I)
	    END IF
	  END DO
	END IF
!
! Per-row bound on the lagged frequency source [MOM_LAG_BND,
! default T].  A lagged fixed-point map amplifies wherever a row's
! lagged coupling weight exceeds the diagonal excess
! TB-|TA|-|TC| (for an M-matrix row, |A^-1|_inf <= 1/excess).
! Shock cells violate this: PSI collapses at the 1+sigma*K=0
! resonance while the lagged H and EPS_PREV transfers still feed
! the row, and the retained field then grows geometrically (m15
! probes, 12-Oct-2026).  Scale the row's whole lagged contribution
! (accumulated in XLAG_NM) by min(1, excess/weight); where the
! unbounded map was already contractive the factor is 1 and the
! fixed point is unchanged.  The implicit/explicit EPS pair match
! of MOM_EPS_LIM becomes approximate on scaled rows -- that is the
! documented consistency cost of this bound.  LAG_SCL_NM is kept
! so the negative-J rebuild below applies the same factor to the
! same terms; the H reconstruction is left unscaled, it solves the
! eliminated H equation rather than the bounded J-row map.
!
	IF(MOM_LAG_BND)THEN
	  DO I=2,ND_NM-1
	    IF(Q_NM(I) .GE. 0.0_LDP)THEN
	      T3=PSIPREV_NM(I)
	    ELSE IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	      T3=PSI_NM(I)
	    ELSE
	      T3=0.0_LDP
	    END IF
	    T4=T3+ABS(VB_NM(I))+ABS(VC_NM(I))+
	1      2.0_LDP*(ABS(EPS_PREV_NM(I-1))+ABS(EPS_PREV_NM(I)))
	    IF(T4 .GT. 0.0_LDP)THEN
	      T5=TB_NM(I)-ABS(TA_NM(I))-ABS(TC_NM(I))
	      T3=MIN(1.0_LDP,MAX(T5,0.0_LDP)/T4)
	    ELSE
	      T3=1.0_LDP
	    END IF
	    LAG_SCL_NM(I)=T3
	    IF(T3 .LT. 1.0_LDP)THEN
	      XM_NM(I)=XM_NM(I)-(1.0_LDP-T3)*XLAG_NM(I)
	      MOM_LAGB_CNT=MOM_LAGB_CNT+1
	      MOM_LAGB_MIN=MIN(MOM_LAGB_MIN,T3)
	    END IF
	  END DO
	END IF
!
!*****************************************************************************
! Outer boundary condition (depth 1).
!
! V11: PSI(1) = GAM(1) * (HBC + NBC * SIGMA(1))
! NONMONO: PSI(1) = |GAM(1) * (HBC + NBC * SIGMA_TRUE(1))| (always >= 0)
!   PSIPREV(1) = GAM(1) * (HBC_PREV + NBC_PREV * SIGMA_TRUE(1)) (signed)
!   where GAM(1) = CON_GAM(1) / CHI_HAT(1) / dLOG_NU
!
! A_BND is this row's own frequency-coupling coefficient,
! A_bnd = GAM(1)*(HBC + NBC*sigma_true(1)) -- the boundary Eddington
! analogue of the interior Q_NM = GAM*(1 + sigma*K_ON_J).  The
! expanding-vs-compressing upwind choice must key on sign(A_bnd):
! near sigma(1) = -1 (or with a negative boundary Eddington) the two
! coefficients differ in sign, and keying on Q_NM(1) classified the
! row as expanding while A_bnd < 0, injecting PSIPREV*J_prev < 0 into
! the RHS.  J(1) then tripped the negative-J rejection every sweep
! (originally the |XM|/10 fixer floor) and corrupted
! HFLUX_AT_OB = HBC_CMF(1)*RJ(1) (audit 24-Sep-2026).  The
! OUT_BC_TYPE > 1 branch still keys on Q_NM(1): that row is the
! negated interior equation, whose coupling coefficient IS Q_NM(1).
!
	BND_HBC_FLOORED=.FALSE.
	IF(OUT_BC_TYPE .LE. 1)THEN
!
! Boundary-Eddington floor (audit 24-Sep-2026, job 11960113): the
! formal solve returns HBC = HNU_AT_OB/JNU_STORE(1), which goes
! negative at near-zero-flux frequencies -- the thick-BC extension
! can emit more inward than the star emits outward (seen at
! ~1230-1270 A and ~22 um).  A negative HBC is an unphysical
! closure for an unilluminated atmosphere and makes this row
! degenerate under EITHER upwind branch: compressing gives the
! lagged J_next weight |A_bnd|/TB(1) > 1 (the diagonal drops below
! |A_bnd| when HBC < -K_ON_J(1)*QF(1)/DTAU(1)), expanding injects
! PSIPREV*J_prev < 0 and lands J(1) on the negative-J rejection.
!
! Floor the Eddington part of the closure at zero for this row
! only: H(1)=0 is the physical boundary condition where the
! measured closure is unusable.  HBC_ROW/HBC_PREV_ROW replace
! HBC/HBC_PREV in A_BND, the diagonal, and PSIPREV consistently so
! both branches stay contractive.  The .NOT.(.GE.0) test also
! catches NaN.  HBC_CMF is untouched, so HFLUX_AT_OB keeps
! reporting the measured flux.  Same defect class as the FEDD->1/3
! guard in FG_J_CMF_NONMONO_V1: a degenerate Eddington measurement
! is replaced by a physical closure.  Floored frequencies are
! logged to MOM_BND_HBC_FLOOR and counted in MOM_BND_HBC_FL.
!
	  HBC_ROW=HBC
	  IF(.NOT.(HBC_ROW .GE. 0.0_LDP))HBC_ROW=0.0_LDP
	  IF(.NOT.(HBC_ROW .LE. 1.0_LDP))HBC_ROW=0.0_LDP
	  HBC_PREV_ROW=HBC_PREV
	  IF(.NOT.(HBC_PREV_ROW .GE. 0.0_LDP))HBC_PREV_ROW=0.0_LDP
	  IF(.NOT.(HBC_PREV_ROW .LE. 1.0_LDP))HBC_PREV_ROW=0.0_LDP
!
! NBC enters A_BND and PSIPREV with the same weight as HBC.  An
! unphysical N/J (negative, non-finite, or above the |mu|<=1 bound
! NBC<=1) leaves the row ill-conditioned even with HBC floored, so
! it is validated with the same counted-floor policy.
!
	  NBC_ROW=NBC
	  IF(.NOT.(NBC_ROW .GE. 0.0_LDP))NBC_ROW=0.0_LDP
	  IF(.NOT.(NBC_ROW .LE. 1.0_LDP))NBC_ROW=0.0_LDP
	  NBC_PREV_ROW=NBC_PREV
	  IF(.NOT.(NBC_PREV_ROW .GE. 0.0_LDP))NBC_PREV_ROW=0.0_LDP
	  IF(.NOT.(NBC_PREV_ROW .LE. 1.0_LDP))NBC_PREV_ROW=0.0_LDP
	  BND_HBC_FLOORED=(HBC_ROW .NE. HBC) .OR.
	1                 (HBC_PREV_ROW .NE. HBC_PREV) .OR.
	2                 (NBC_ROW .NE. NBC) .OR.
	3                 (NBC_PREV_ROW .NE. NBC_PREV)
	  IF(BND_HBC_FLOORED)THEN
	    IF(.NOT. BND_FL_OPEN)THEN
	      OPEN(UNIT=49,FILE='MOM_BND_HBC_FLOOR',STATUS='UNKNOWN',
	1           POSITION='APPEND',ACTION='WRITE')
	      BND_FL_OPEN=.TRUE.
	    END IF
	    WRITE(49,'(A,I8,A,ES16.8)')
	1     '  bnd HBC floor: FREQ_INDX=',FREQ_INDX,' FREQ=',FREQ
	    WRITE(49,'(A,4ES13.4)')
	1     '    HBC,NBC,HBC_PREV,NBC_PREV=',
	2     HBC,NBC,HBC_PREV,NBC_PREV
	    WRITE(49,'(A,3ES13.4)')
	1     '    sigma_true(1),JNU(1),HNU_OB=',
	2     SIGMA_TRUE_NM(1),JNU_STORE(1),HNU_AT_OB
	  END IF
	  A_BND=GAM_NM(1)*(HBC_ROW+NBC_ROW*SIGMA_TRUE_NM(1))
	  PSI_NM(1)=ABS(A_BND)
	  IF(MOM_GAM_PIN)THEN
!
! Same pin as the interior PSIPREV: the boundary coupling weight is
! this row's own A_bnd, not the bluer frequency's Eddington factors
! (HBC_PREV/NBC_PREV), whose drift past the sign flip let the
! explicit upwind source outweigh the implicit diagonal.
! [MOM_GAM_PIN]
!
	    PSIPREV_NM(1)=PSI_NM(1)
	  ELSE
	    PSIPREV_NM(1)=ABS(GAM_NM(1)*
	1       (HBC_PREV_ROW+NBC_PREV_ROW*SIGMA_TRUE_NM(1)))
	  END IF
	  IF(.NOT.(DTAU_NM(1) .GT. 0.0_LDP))THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Invalid outer-boundary DTAU:',DTAU_NM(1),
	1       ' FREQ_INDX=',FREQ_INDX
	    STOP
	  END IF
	  TC_NM(1)=-K_ON_J(2)*QF_NM(2)/DTAU_NM(1)
	  TB_NM(1)=K_ON_J(1)*QF_NM(1)/DTAU_NM(1) + PSI_NM(1) + HBC_ROW
	  TA_NM(1)=0.0_LDP
	  VB_NM(1)=0.0_LDP
	  VC_NM(1)=0.0_LDP
!
! Same row bound as the interior [MOM_LAG_BND]: this row's lagged
! weight is the upwind-J coefficient itself, its diagonal excess is
! TB(1)-|TC(1)|.
!
	  LAG_SCL_NM(1)=1.0_LDP
	  IF(MOM_LAG_BND)THEN
	    IF(A_BND .GE. 0.0_LDP)THEN
	      T3=PSIPREV_NM(1)
	    ELSE IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	      T3=PSI_NM(1)
	    ELSE
	      T3=0.0_LDP
	    END IF
	    IF(T3 .GT. 0.0_LDP)THEN
	      T4=TB_NM(1)-ABS(TC_NM(1))
	      LAG_SCL_NM(1)=MIN(1.0_LDP,MAX(T4,0.0_LDP)/T3)
	      IF(LAG_SCL_NM(1) .LT. 1.0_LDP)THEN
	        MOM_LAGB_CNT=MOM_LAGB_CNT+1
	        MOM_LAGB_MIN=MIN(MOM_LAGB_MIN,LAG_SCL_NM(1))
	      END IF
	    END IF
	  END IF
	  IF(A_BND .GE. 0.0_LDP)THEN
	    XM_NM(1)=LAG_SCL_NM(1)*PSIPREV_NM(1)*RSQJNU_PREV_NM(1)
	  ELSE
	    IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	      XM_NM(1)=LAG_SCL_NM(1)*PSI_NM(1)*JUP_NM(1)
	    ELSE
	      XM_NM(1)=0.0_LDP
	    END IF
	  END IF
	ELSE
	  PSI_NM(1)=ABS(GAM_NM(1)*GAM_FAC_NM(1))
	  PSIPREV_NM(1)=ABS(GAM_NM(1)*GAM_FAC_PREV_NM(1))
	  T1=0.25_LDP*(CHI_HAT_NM(2)+CHI_HAT_NM(1))*(R_NM(2)-R_NM(1))
	  TC_NM(1)=HU_NM(1)/T1
! This row is the negated interior equation, so the same ESEC>CHI
! defect enters as +(COH-1); MIN gives -(floored absorption).  The
! MOM_ABS_FLR cap on the MIN argument keeps at least that fraction of
! absorption in the row, matching the interior MAX floor.
	  TB_NM(1)= (MIN(COH_VEC_NM(1),1.0_LDP-MOM_ABS_FLR)-1.0_LDP) -
	1   (HL_NM(1)+HBC)/T1 - PSI_NM(1)
! Midpoint 1 sees row 1 as the outer node; this row is the negated
! interior row, so the depth-upwind EPS term lands negated on the
! diagonal (EPS>=0) or negated on the superdiagonal (EPS<0).  Factor
! 2 as in the interior: J_mid=J_up in the Jsum=2*J_mid convention.
	  IF(EPS_NM(1) .GE. 0.0_LDP)THEN
	    TB_NM(1)=TB_NM(1)-2.0_LDP*EPS_NM(1)/T1
	  ELSE
	    TC_NM(1)=TC_NM(1)-2.0_LDP*EPS_NM(1)/T1
	  END IF
!
! Same row bound as the interior [MOM_LAG_BND].  This row is the
! negated interior equation, so the implicit excess is
! |TB(1)|-|TC(1)| and the lagged weight is the J-coupling
! coefficient plus (HS + 2|EPS_PREV|)/|T1|.
!
	  LAG_SCL_NM(1)=1.0_LDP
	  IF(MOM_LAG_BND)THEN
	    IF(Q_NM(1) .GE. 0.0_LDP)THEN
	      T3=PSIPREV_NM(1)
	    ELSE IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	      T3=PSI_NM(1)
	    ELSE
	      T3=0.0_LDP
	    END IF
	    T4=T3+(HS_NM(1)+2.0_LDP*ABS(EPS_PREV_NM(1)))/ABS(T1)
	    IF(T4 .GT. 0.0_LDP)THEN
	      T5=ABS(TB_NM(1))-ABS(TC_NM(1))
	      LAG_SCL_NM(1)=MIN(1.0_LDP,MAX(T5,0.0_LDP)/T4)
	      IF(LAG_SCL_NM(1) .LT. 1.0_LDP)THEN
	        MOM_LAGB_CNT=MOM_LAGB_CNT+1
	        MOM_LAGB_MIN=MIN(MOM_LAGB_MIN,LAG_SCL_NM(1))
	      END IF
	    END IF
	  END IF
	  IF(Q_NM(1) .GE. 0.0_LDP)THEN
	    XM_NM(1)=-SOURCE_NM(1)*R_NM(1)*R_NM(1) -
	1     LAG_SCL_NM(1)*(HS_NM(1)*HUP_NM(1)/T1 +
	2     PSIPREV_NM(1)*RSQJNU_PREV_NM(1))
	  ELSE
! This boundary row is the negated interior equation (note -SOURCE and
! -PSIPREV*J_prev for expansion), so the compressing upwind term enters
! as -PSI_NM*J_next, with the same diagonal |Q_k| weight as interior.
	    IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	      XM_NM(1)=-SOURCE_NM(1)*R_NM(1)*R_NM(1) -
	1       LAG_SCL_NM(1)*(HS_NM(1)*HUP_NM(1)/T1 +
	2       PSI_NM(1)*JUP_NM(1))
	    ELSE
	      XM_NM(1)=-SOURCE_NM(1)*R_NM(1)*R_NM(1) -
	1       LAG_SCL_NM(1)*HS_NM(1)*HUP_NM(1)/T1
	    END IF
	  END IF
	  IF(N_TYPE .NE. 'G_ONLY')THEN
! Lagged midpoint-1 EPS term on the negated boundary row, at the
! same depth-upwind node as the implicit part (outer node J_1 for
! EPS_PREV>=0, inner node J_2 for EPS_PREV<0) and the same frequency
! upwind as the interior scatter (bluer J_prev expanding, redder
! JUP_NM compressing).
	    IF(EPS_PREV_NM(1) .GE. 0.0_LDP)THEN
	      K=1
	    ELSE
	      K=2
	    END IF
	    IF(EPS_PREV_NM(1) .LT. 0.0_LDP)THEN
	      T4=JUP_NM(K)
	    ELSE
	      T4=RSQJNU_PREV_NM(K)
	    END IF
	    T5=2.0_LDP*EPS_PREV_NM(1)*T4/T1
!
!	    Same limiter as the interior scatter: the lagged debit may
!	    not exceed row 1's accumulated positive RHS.  Row 1 is the
!	    negated equation, so the transfer lands with opposite sign.
!	    The LAG_SCL factor applies the row bound after the cap.
!
	    IF(MOM_EPS_LIM .AND. T5 .GT. 0.0_LDP)THEN
	      T4=MAX(XM_NM(1),0.0_LDP)
	      IF(T5 .GT. T4)THEN
	        T5=T4
	        MOM_EPS_CAP=MOM_EPS_CAP+1
	        MOM_EPS_CAP_SW=MOM_EPS_CAP_SW+1
	      END IF
	    END IF
	    XM_NM(1)=XM_NM(1)-LAG_SCL_NM(1)*T5
	  END IF
	END IF
!
!*****************************************************************************
! Solve the tridiagonal system for R^2*J.
!
! DGTSV (LAPACK, partial pivoting) is retained as the solve: it was
! introduced while the centered EPS stencil made the matrix
! non-M-matrix, and stays because it is stencil-agnostic.  The
! coefficient arrays are copied because DGTSV overwrites them during
! factorization.
!
	IF(.NOT. ALLOCATED(GT_DL))THEN
	  ALLOCATE(GT_DL(ND_NM),GT_D(ND_NM),GT_DU(ND_NM))
	ELSE IF(SIZE(GT_DL) .NE. ND_NM)THEN
	  DEALLOCATE(GT_DL,GT_D,GT_DU)
	  ALLOCATE(GT_DL(ND_NM),GT_D(ND_NM),GT_DU(ND_NM))
	END IF
!
! Repair mask reset moved ahead of the solve so a DGTSV failure can
! mark the whole frequency as rejected and the post-solve guards do
! not erase the mark.
!
	XM_BAD_NM(1:ND_NM)=.FALSE.
!
! Previous valid field at this frequency, the fallback target for a
! failed solve or a rejected cell: the lagged GS iterate
! J_NEXT_FRQ(:,FREQ_INDX) where positive, else the local-equilibrium
! estimate SOURCE*R^2.  The bluer-frequency RSQJNU_PREV_NM is not a
! valid fallback here: on a cold sweep every failed cell would
! inherit its neighbour's value, which ratchets R^2*J upward along
! the frequency sweep and drives the iteration-1 luminosity to 1e20.
! XM_SAVE_NM below is the assembled RHS -- restoring it after a
! failure put source-scale values into R^2*J and propagated them
! through J_NEXT_FRQ, RSQJNU_PREV and BA.
!
	DO I=1,ND_NM
	  XM_FB_NM(I)=0.0_LDP
	  IF(SOURCE_NM(I) .GT. 0.0_LDP)
	1    XM_FB_NM(I)=SOURCE_NM(I)*R_NM(I)*R_NM(I)
	END DO
	IF(ALLOCATED(J_NEXT_FRQ))THEN
	  IF(FREQ_INDX .GE. 1 .AND.
	1     FREQ_INDX .LE. SIZE(J_NEXT_FRQ,2))THEN
	    DO I=1,ND_NM
	      IF(J_NEXT_FRQ(I,FREQ_INDX) .GT. 0.0_LDP)
	1        XM_FB_NM(I)=J_NEXT_FRQ(I,FREQ_INDX)
	    END DO
	  END IF
	END IF
!
! Matrix-quality audit: minimum |TB|/(|TA|+|TC|) over interior rows,
! non-positive diagonals, and the ESEC/CHI overshoot.  Swept up
! across frequencies; the sweep report shows whether the ~30% stall
! correlates with conditioning.
!
	DO I=2,ND_NM-1
	  T1=ABS(TA_NM(I))+ABS(TC_NM(I))
	  IF(T1 .GT. 0.0_LDP)THEN
	    T2=TB_NM(I)/T1
	  ELSE
	    T2=1.0E+30_LDP
	  END IF
	  IF(T2 .LT. MOM_MIN_DR)THEN
	    MOM_MIN_DR=T2; MOM_MIN_DR_I=I; MOM_MIN_DR_L=FREQ_INDX
	  END IF
	  IF(ABS(T2) .LT. MOM_MIN_DRA)THEN
	    MOM_MIN_DRA=ABS(T2)
	    MOM_MIN_DRA_I=I;  MOM_MIN_DRA_L=FREQ_INDX
	  END IF
	  IF(TB_NM(I) .LE. 0.0_LDP)MOM_TB_NONPOS=MOM_TB_NONPOS+1
	  IF(COH_VEC_NM(I) .GT. 1.0_LDP)
	1    MOM_MAX_COH=MAX(MOM_MAX_COH,COH_VEC_NM(I)-1.0_LDP)
	END DO
	MOM_MAX_GAM=MAX(MOM_MAX_GAM,MAXVAL(ABS(GAM_NM(1:ND_NM))))
	MOM_MAX_SIGT=MAX(MOM_MAX_SIGT,
	1    MAXVAL(ABS(SIGMA_TRUE_NM(1:ND_NM))))
!
	XM_SAVE_NM=XM_NM
	GT_DL(1:ND_NM-1)=TA_NM(2:ND_NM)
	GT_DU(1:ND_NM-1)=TC_NM(1:ND_NM-1)
	GT_D(1:ND_NM)=TB_NM(1:ND_NM)
	CALL DGTSV(ND_NM,IONE,GT_DL,GT_D,GT_DU,XM_NM,ND_NM,GT_INFO)
	IF(GT_INFO .NE. 0)THEN
	  MOM_N_DGTSV=MOM_N_DGTSV+1
	  WRITE(LUER,*)'DGTSV failed in MOM_J_CMF_NONMONO_V1: INFO=',
	1    GT_INFO,' FREQ_INDX=',FREQ_INDX
	  IF(GT_INFO .GE. 1 .AND. GT_INFO .LE. ND_NM)
	1    WRITE(LUER,*)'  TA,TB,TC at bad row:',TA_NM(GT_INFO),
	2    TB_NM(GT_INFO),TC_NM(GT_INFO),
	3    '  diag range:',MINVAL(TB_NM),MAXVAL(TB_NM)
!
! Whole-frequency quarantine: the tridiagonal couples every depth, so
! no part of a failed solve is trustworthy.  Restore the previous
! valid field (not the RHS) and mark every cell bad so the GS metric
! counts a rejection instead of a fake converged field, and the
! maturity gate forces Lambda iteration.
!
	  XM_NM(1:ND_NM)=XM_FB_NM(1:ND_NM)
	  XM_BAD_NM(1:ND_NM)=.TRUE.
	  MOM_GS_MAXCHG=MAX(MOM_GS_MAXCHG,1.0E+30_LDP)
	  MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+ND_NM
	END IF
!
! Equation residual |A.x - b| / |b| on the post-solve field before
! repairs.  Skipped on failure: XM_NM then holds the fallback field,
! not a solution of this system.
!
	IF(GT_INFO .EQ. 0)THEN
	  T2=MAX(MAXVAL(ABS(XM_SAVE_NM(1:ND_NM))),1.0E-30_LDP)
	  DO I=2,ND_NM-1
	    T1=ABS(TA_NM(I)*XM_NM(I-1)+TB_NM(I)*XM_NM(I)+
	1         TC_NM(I)*XM_NM(I+1)-XM_SAVE_NM(I))/
	2         MAX(ABS(XM_SAVE_NM(I)),1.0E-6_LDP*T2)
	    IF(T1 .GT. MOM_EQ_RES)THEN
	      MOM_EQ_RES=T1; MOM_EQ_RES_I=I; MOM_EQ_RES_L=FREQ_INDX
	    END IF
	  END DO
	END IF
!
! Guard against extremely negative solutions (same logic as V11).
!
	IF(MINVAL(XM_NM) .LE. -1.0E+20_LDP)THEN
	  DO I=2,ND_NM-1
	    ICNT=0
	    T1=1.0_LDP
	    IF(XM_NM(I) .LE. 0.0_LDP .AND. XM_SAVE_NM(I) .LE. 0.0_LDP)THEN
              T1=T1*0.5_LDP
! Rebuild the explicit upwind coupling with the same sign-dependent
! upwind field as the main solve: J_prev (expanding) or the
! previous-iteration J_next with the POSITIVE |coeff| (compressing).
	      IF(Q_NM(I) .GE. 0.0_LDP)THEN
	        T2=PSIPREV_NM(I)*RSQJNU_PREV_NM(I)
	      ELSE IF(J_NEXT_VALID .AND. FREQ_INDX .LT. NCF_IN)THEN
	        T2=PSI_NM(I)*JUP_NM(I)
	      ELSE
	        T2=0.0_LDP
	      END IF
	      IF(N_TYPE .EQ. 'G_ONLY')THEN
	        XM_NM(I)=DTAUONQ_NM(I)*SOURCE_NM(I)*R_NM(I)*R_NM(I) +
	1          LAG_SCL_NM(I)*(T2 +
	2          T1*(VB_NM(I)*HUP_NM(I-1)+VC_NM(I)*HUP_NM(I)))
	      ELSE
! Rebuild the lagged EPS coupling at the same depth-upwind nodes
! and the same frequency upwind as the main RHS (J_prev expanding,
! JUP_NM compressing).
	        T3=0.0_LDP
	        DO J=I-1,I
	          IF(EPS_PREV_NM(J) .GE. 0.0_LDP)THEN
	            K=J
	          ELSE
	            K=J+1
	          END IF
	          IF(EPS_PREV_NM(J) .LT. 0.0_LDP)THEN
	            T4=JUP_NM(K)
	          ELSE
	            T4=RSQJNU_PREV_NM(K)
	          END IF
	          T4=2.0_LDP*EPS_PREV_NM(J)*T4
	          IF(J .EQ. I-1)THEN
	            T3=T3-T4
	          ELSE
	            T3=T3+T4
	          END IF
	        END DO
!
!	        Rejected-row rebuild mirrors the sweep's lagged-EPS
!	        limiter and the MOM_LAG_BND row bound: the scaled lagged
!	        debit may not drive the rebuilt balance negative, or the
!	        fallback re-seeds the same negative RHS that just failed.
!	        [MOM_EPS_LIM, MOM_LAG_BND]
!
	        T4=T2 + T1*(VB_NM(I)*HUP_NM(I-1)+VC_NM(I)*HUP_NM(I))
	        T5=DTAUONQ_NM(I)*SOURCE_NM(I)*R_NM(I)*R_NM(I)
	        IF(MOM_EPS_LIM .AND. T1*T3 .LT. 0.0_LDP .AND.
	1           T5+LAG_SCL_NM(I)*(T4+T1*T3) .LT. 0.0_LDP)THEN
	          T3=-(T5+LAG_SCL_NM(I)*T4)/(LAG_SCL_NM(I)*T1)
	          MOM_EPS_CAP=MOM_EPS_CAP+1
	          MOM_EPS_CAP_SW=MOM_EPS_CAP_SW+1
	        END IF
	        XM_NM(I)=T5 + LAG_SCL_NM(I)*(T4+T1*T3)
	      END IF
	      ICNT=ICNT+1
	    ELSE
	      XM_NM(I)=XM_SAVE_NM(I)
	    END IF
	  END DO
	  XM_NM(1)=XM_SAVE_NM(1); XM_NM(ND_NM)=XM_SAVE_NM(ND_NM)
	  GT_DL(1:ND_NM-1)=TA_NM(2:ND_NM)
	  GT_DU(1:ND_NM-1)=TC_NM(1:ND_NM-1)
	  GT_D(1:ND_NM)=TB_NM(1:ND_NM)
	  CALL DGTSV(ND_NM,IONE,GT_DL,GT_D,GT_DU,XM_NM,ND_NM,GT_INFO)
	  IF(GT_INFO .NE. 0)THEN
	    WRITE(LUER,*)'DGTSV (fallback) failed: INFO=',GT_INFO
!
! Same rule as the primary failure: restore the previous valid field,
! never the assembled RHS.
!
	    MOM_N_DGTSV=MOM_N_DGTSV+1
	    XM_NM(1:ND_NM)=XM_FB_NM(1:ND_NM)
	    XM_BAD_NM(1:ND_NM)=.TRUE.
	    MOM_GS_MAXCHG=MAX(MOM_GS_MAXCHG,1.0E+30_LDP)
	    MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+ND_NM
	  END IF
	END IF
!
! Check for negative mean intensities.
!
! XM_BAD_NM marks cells touched by any post-solve repair below so the
! J_NEXT_FRQ update can optionally (MOM_RF_KEEPNEG) keep the previous
! iterate there instead of storing the repaired value.  The mask was
! reset before the solve so a DGTSV failure mark survives.
!
	IF(MINVAL(XM_NM(1:ND_NM)) .LE. 0.0_LDP .AND. VERBOSE_NM)THEN
	  WRITE(48,'(/,A,E16.8)')'MOM_J_CMF_NONMONO_V1: negative J, Freq=',FREQ
	  TA_NM(1:ND_NM)=XM_NM(1:ND_NM)/R_NM(1:ND_NM)/R_NM(1:ND_NM)
	  CALL WRITE_VEC(TA_NM,ND_NM,'XM Vec',48)
	  CALL WRITE_VEC(K_ON_J,ND_NM,'K_ON_J (FEDD) Vec',48)
	  IF(N_TYPE .NE. 'G_ONLY')CALL WRITE_VEC(NMID_ON_J,ND_NM,'NMID_ON_J Vec',48)
	  IF(N_TYPE .NE. 'N_ON_J')CALL WRITE_VEC(NMID_ON_HMID,ND_NM,'NMID_ON_HMID (G) Vec',48)
	  CALL WRITE_VEC(ETA_NM,ND_NM,'ETA Vec',48)
	  CALL WRITE_VEC(ESEC_NM,ND_NM,'ESEC Vec',48)
	  CALL WRITE_VEC(CHI_NM,ND_NM,'CHI Vec',48)
	  CALL WRITE_VEC(CHI_HAT_NM,ND_NM,'CHI_HAT Vec',48)
	  CALL WRITE_VEC(SIGMA_TRUE_NM,ND_NM,'SIGMA_TRUE Vec',48)
	END IF
!
	RECORDED_ERROR=.FALSE.
!
! Significant-negative threshold, scale-local: T4 is 1E-6 of the
! column peak (no absolute floor -- an earlier MAX(...,1.0) meant a
! -1E-20 solution in a cold EUV column with R^2*J ~ 1E-20 counted as
! "tiny" even though it was 100% of the local field), and T5 raises
! it to 1E-3 of the cell's previous valid value so a sign flip of
! order the local field is always rejected while roundoff on a
! strong cell is not.  MAXVAL is masked to exclude non-finite
! entries.
!
	T4=MAXVAL(ABS(XM_NM(1:ND_NM)),
	1    MASK=(ABS(XM_NM(1:ND_NM)) .LE. HUGE(1.0_LDP)))
	T4=MAX(T4,1.0E-280_LDP)*1.0E-6_LDP
	DO I=1,ND_NM
	  IF(XM_NM(I) .LT. 0.0_LDP)THEN
	    T5=MAX(T4,1.0E-3_LDP*ABS(XM_FB_NM(I)))
	    IF(ABS(XM_NM(I)) .GT. T5)THEN
	      MOM_NEGJ_LRG=MOM_NEGJ_LRG+1
	    ELSE
	      MOM_N_NEGJ_T=MOM_N_NEGJ_T+1
	    END IF
	    IF(.NOT. VERBOSE_NM .AND.
	1    MOM_NEGJ_DUMP .LT. MOM_NEGJ_DUMP_MAX)THEN
	      MOM_NEGJ_DUMP=MOM_NEGJ_DUMP+1
	      WRITE(48,'(I5,ES16.8,12ES13.4)')I,FREQ,XM_NM(I),ETA_NM(I),
	1         CHI_NM(I),ESEC_NM(I),K_ON_J(I),
	1         NMID_ON_J(I),NMID_ON_HMID(I),XM_NM(MAX(1,I-2):MIN(I+2,ND_NM))
	      WRITE(48,'(T112,5ES13.4)')XM_SAVE_NM(MAX(1,I-2):MIN(I+2,ND_NM))
	      WRITE(48,'(A,T38,4ES13.4)')'  sigma_true,gam,chi,chi_hat=',
	1         SIGMA_TRUE_NM(I),GAM_NM(I),CHI_NM(I),CHI_HAT_NM(I)
!
! Row diagnostics for the negJ investigation (24-Sep-2026): cells that
! flip negative at sigma~0 with a positive RHS imply TB<=0, which can
! only come from the absorption term (1-COH)*DTAUONQ when the coherent
! fraction ESEC/CHI exceeds 1 at a thick cell.  TB/TA/TC survive the
! solve (DGTSV works on copies).
!
	      WRITE(48,'(A,T38,6ES13.4)')'  TB,TA,TC,absorp,COH,DTAUONQ=',
	1         TB_NM(I),TA_NM(I),TC_NM(I),
	2         DTAUONQ_NM(I)*(1.0_LDP-COH_VEC_NM(I)),COH_VEC_NM(I),
	3         DTAUONQ_NM(I)
	    END IF
	    IF(ABS(XM_NM(I)) .GT. T5)THEN
!
!	      Significant negative: a failed solve, not noise.  Reject the
!	      update to the previous valid field rather than |XM|/10,
!	      which is an arbitrary positive value the lagged map then
!	      chases.
!
	      XM_NM(I)=XM_FB_NM(I)
	      XM_BAD_NM(I)=.TRUE.
	    ELSE
!
!	      Roundoff-level negative (below both 1D-6 of the column peak
!	      and 1D-3 of the local previous value): keep the magnitude;
!	      the sign flip is solver noise on J~0.
!
	      XM_NM(I)=ABS(XM_NM(I))
	    END IF
	    IF(.NOT. RECORDED_ERROR)THEN
	      IF(MOM_ERR_CNT .EQ. 0)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      ELSE IF(MOM_ERR_CNT .GT. N_ERR_MAX)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	      ELSE IF(MOM_ERR_ON_FREQ(MOM_ERR_CNT) .NE. FREQ)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        IF(MOM_ERR_CNT .LT. N_ERR_MAX)MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      END IF
	      RECORDED_ERROR=.TRUE.
	    END IF
	  END IF
	END DO
!
! A non-finite XM_NM must not be stored: RSQJNU_NM and J_NEXT_FRQ are
! lagged fields read back as explicit couplings by later frequencies
! (and by later great iterations at compressing depths), so a stored
! NaN/Inf injects NaN into every subsequent tridiagonal.  A failed
! entry is rejected to the previous valid field XM_FB_NM -- zeroing
! dropped the lagged coupling entirely, and the old 1E40 clip stored
! a wrong field at a wrong scale that the next solve's RHS coupling
! (PSI*JUP etc) then amplified.
!
	N_QBND=0
!
! Runaway reference: the largest previous-valid value at this
! frequency.  When no valid history exists (cold start with zero
! SOURCE and zero RSQJNU_PREV) the check is skipped -- there is no
! scale to reject against, and a blanket floor would reject a
! legitimate first solution.  The 1E3 bound is a hard numerical
! rejection, not a damping choice: a single GS step that moves a
! cell three dex above the previous valid field is a failed solve,
! and the retained iterate keeps the map near a valid state.
!
	T3=MAX(MAXVAL(XM_FB_NM(1:ND_NM)),0.0_LDP)
	DO I=1,ND_NM
	  IF(.NOT. (ABS(XM_NM(I)) .LE. HUGE(1.0_LDP)))THEN
	    MOM_N_NONF=MOM_N_NONF+1
	    XM_NM(I)=XM_FB_NM(I)
	    XM_BAD_NM(I)=.TRUE.
	    N_QBND=N_QBND+1
	  ELSE IF(T3 .GT. 1.0E-20_LDP .AND. ABS(XM_NM(I)) .GT.
	1      1.0E+3_LDP*MAX(XM_FB_NM(I),1.0E-3_LDP*T3))THEN
	    MOM_N_RUNAWAY=MOM_N_RUNAWAY+1
	    XM_NM(I)=XM_FB_NM(I)
	    XM_BAD_NM(I)=.TRUE.
	    N_QBND=N_QBND+1
	  END IF
	END DO
	IF(N_QBND .GT. 0 .AND. N_XPRT .LT. 15)THEN
	  N_XPRT=N_XPRT+1
	  WRITE(6,'(A,ES14.6,A,I5,A)')' MOM_J_CMF_NONMONO_V1: '//
	1     'non-finite/runaway XM_NM rejected at',FREQ,' :',
	2     N_QBND,' depths'
! Localize the first failure: a single non-finite tridiagonal row or
! RHS entry propagates through the tridiagonal solve, so the first bad index
! (not the count) identifies the source.  Also scan the matrix/RHS so
! an upstream NaN is reported at its origin, not its first XM echo.
	  DO I=1,ND_NM
	    IF(.NOT. (ABS(XM_SAVE_NM(I)) .LE. HUGE(1.0_LDP)) .OR.
	1         .NOT. (ABS(TA_NM(I)) .LE. HUGE(1.0_LDP)) .OR.
	1         .NOT. (ABS(TB_NM(I)) .LE. HUGE(1.0_LDP)) .OR.
	1         .NOT. (ABS(TC_NM(I)) .LE. HUGE(1.0_LDP)))THEN
	      WRITE(6,'(A,I6,A,4ES13.4,A,5ES13.4)')'   first bad row I=',I,
	1       '  TA,TB,TC,XMrhs=',TA_NM(I),TB_NM(I),TC_NM(I),
	1       XM_SAVE_NM(I),'  DTAUONQ,SRC,R,CHI_HAT,K_ON_J=',
	1       DTAUONQ_NM(I),SOURCE_NM(I),R_NM(I),CHI_HAT_NM(I),
	1       K_ON_J(I)
	      EXIT
	    END IF
	  END DO
	END IF
!
! Quarantined cells already hold XM_FB_NM (the previous valid field,
! which is J_NEXT_FRQ where that was positive) from the post-solve
! guards above, so RSQJNU_NM, JNU and the relaxed J_NEXT_FRQ write
! below all see the same frozen value -- the old MOM_BA_FIX block
! that re-froze only under a flag is removed.  The unconditional
! DO_VAR_CONT suppression of the ALO diagonal at MOM_CELL_QUAR cells
! is the matching statement that a frozen cell has zero response.
!
! Save R^2*J for next frequency.
!
	RSQJNU_NM(1:ND_NM)=XM_NM(1:ND_NM)
!
! Update J_NEXT_FRQ for Gauss-Seidel convergence.
! Store the current R^2*J at all depths for this frequency; during the
! solve for frequency k the slot k+1 still holds the previous great
! iteration's value (read-before-write), giving compressing zones the
! lagged redder-frequency upwind field.
!
! MOM_GS_MAXCHG/MOM_INEXT_CAPPED accumulate the last-sweep convergence
! metrics written into the GS checkpoint; a complete-but-still-moving
! restored field does not skip the restart warm-up.  The accumulator
! reset now lives at the top of the routine (before this frequency's
! solve contributes), gated off during refresh.
!
	IF(FREQ_INDX .GE. 1 .AND. FREQ_INDX .LE. NCF_IN)THEN
	  MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+N_QBND
	  IF(BND_HBC_FLOORED)MOM_BND_HBC_FL=MOM_BND_HBC_FL+1
!
!	  Denominator floor: the larger of the cell's previous value and
!	  1e-4 of this frequency's column peak, so cells carrying a
!	  negligible fraction of the field cannot pin the max-norm (and
!	  the maturity gate) at irrelevant floor oscillations.
!
	  T3=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,FREQ_INDX))),
	1    1.0E-28_LDP)
	  DO I=1,ND_NM
	    IF(MOM_KN_REFRESHING .AND. MOM_RF_OWNW)THEN
!
!	      Non-owned rows are not part of this pass's map -- their
!	      store write is suppressed below, so the solve-vs-store
!	      difference persists by construction.  Counting it would
!	      pin the refresh residual at the inter-direction
!	      discrepancy forever; the row is measured on its owning
!	      pass instead.  [MOM_RF_OWNW]
!
	      L_SUPP=(KN_M_DESC .NEQV. (Q_NM(I) .LT. 0.0_LDP))
	      IF(I .EQ. 1 .AND. OUT_BC_TYPE .LE. 1)L_SUPP=
	1          (KN_M_DESC .NEQV. (A_BND .LT. 0.0_LDP))
!
!	      Record which pass direction owns the row for the
!	      composite best-field merge in MOM_KNOP_REFRESH_V1;
!	      Q_NM/A_BND are overwritten by the next frequency's
!	      solve, so ownership is not recoverable later.
!
	      MOM_OWN_CMP(I,FREQ_INDX)=(Q_NM(I) .LT. 0.0_LDP)
	      IF(I .EQ. 1 .AND. OUT_BC_TYPE .LE. 1)
	1          MOM_OWN_CMP(I,FREQ_INDX)=(A_BND .LT. 0.0_LDP)
	      IF(L_SUPP)CYCLE
	    END IF
	    T1=J_NEXT_FRQ(I,FREQ_INDX)
	    T2=XM_NM(I)
	    IF(XM_BAD_NM(I) .AND. MOM_RF_KEEPNEG)THEN
!
!	      Quarantined cell: the stored value is frozen by
!	      construction, so solve-vs-stored disagreement measures the
!	      repair, not field drift -- count separately and exclude
!	      from the residual metrics.
!
	      MOM_GS_NQUAR=MOM_GS_NQUAR+1
!
!	      edited by Devin, 29-Sep-2026: age the quarantine on the
!	      sweep calls only.  Refresh re-solves the same captured
!	      inputs many times per sweep and would inflate the age.
!
	      IF(.NOT. MOM_KN_REFRESHING)THEN
	        IF(MOM_QUAR_AGE(I,FREQ_INDX) .LT. 32767_2)
	1          MOM_QUAR_AGE(I,FREQ_INDX)=
	2            MOM_QUAR_AGE(I,FREQ_INDX)+1_2
	        IF(MOM_QUAR_AGE(I,FREQ_INDX) .GE. MOM_QUAR_PERSIST)
	1          MOM_GS_NPQUAR=MOM_GS_NPQUAR+1
	      END IF
	    ELSE IF(XM_BAD_NM(I))THEN
!
!	      Rejected cell under the default path: the stored field was
!	      restored to the previous iterate, so |T2-T1| is ~0 and the
!	      drift measure would read converged on a failed solve.
!	      Count the rejection as a >50% event.
!
	      MOM_GS_NBAD5=MOM_GS_NBAD5+1
	      MOM_GS_NBAD50=MOM_GS_NBAD50+1
	      IF(MOM_GS_MAXCHG .LT. 1.0_LDP)THEN
	        MOM_GS_MAXCHG=1.0_LDP
	        MOM_GS_ARGI=I;  MOM_GS_ARGL=FREQ_INDX
	      END IF
	    ELSE IF(.NOT.(ABS(T1) .LE. HUGE(T1)) .OR.
	1       .NOT.(ABS(T2) .LE. HUGE(T2)))THEN
!
!	      Non-finite entry in either field: flag the sweep as fully
!	      unconverged rather than letting NaN propagate through MAX.
!
	      MOM_GS_MAXCHG=MAX(MOM_GS_MAXCHG,1.0E+30_LDP)
	      MOM_GS_ARGI=I;  MOM_GS_ARGL=FREQ_INDX
	      MOM_GS_NBAD5=MOM_GS_NBAD5+1
	      MOM_GS_NBAD50=MOM_GS_NBAD50+1
	      MOM_GS_NHUGE=MOM_GS_NHUGE+1
	    ELSE IF(ABS(T1) .GT. 1.0E-28_LDP)THEN
!
!	      Denominator floor prevents a vanishing old value from
!	      producing an arbitrarily large (or non-finite) ratio.
!
	      T4=ABS(T2-T1)/MAX(ABS(T1),1.0E-4_LDP*T3)
	      IF(T4 .GT. 0.05_LDP)MOM_GS_NBAD5=MOM_GS_NBAD5+1
	      IF(T4 .GT. 0.5_LDP)MOM_GS_NBAD50=MOM_GS_NBAD50+1
	      IF(T4 .GT. 1.0E+6_LDP)MOM_GS_NHUGE=MOM_GS_NHUGE+1
	      IF(T4 .GT. MOM_GS_MAXCHG)THEN
	        MOM_GS_MAXCHG=T4
	        MOM_GS_ARGI=I;  MOM_GS_ARGL=FREQ_INDX
	      END IF
	    ELSE IF(ABS(T2) .GT. 1.0E-4_LDP*T3)THEN
	      MOM_GS_NBAD5=MOM_GS_NBAD5+1
	      MOM_GS_NBAD50=MOM_GS_NBAD50+1
	      IF(MOM_GS_MAXCHG .LT. 1.0_LDP)THEN
	        MOM_GS_MAXCHG=1.0_LDP
	        MOM_GS_ARGI=I;  MOM_GS_ARGL=FREQ_INDX
	      END IF
	    END IF
!
!	    A cell that solved cleanly (or was repaired under the default
!	    non-KEEPNEG path, which writes the repaired value and so is
!	    not frozen) is not quarantined this sweep -- its age restarts.
!
	    IF(.NOT. MOM_KN_REFRESHING .AND. .NOT.
	1      (XM_BAD_NM(I) .AND. MOM_RF_KEEPNEG))
	2      MOM_QUAR_AGE(I,FREQ_INDX)=0_2
	  END DO
!
! Under-relaxed write (JNEXT_RELAX): damps the GS fixed-point map on
! compressing grids where straight replacement amplifies the SE<->RT
! loop; fixed point unchanged.  Cells with no prior value get the new
! field directly.
!
	  DO I=1,ND_NM
	    IF(J_NEXT_FRQ(I,FREQ_INDX) .GT. 0.0_LDP)THEN
!
!	      Compression-dominant write suppression: on a dominant-DESC
!	      grid the ascending main sweep's compressing upwind read is
!	      a one-iteration-old Jacobi value, so relaxing the store
!	      toward that output erodes the refresh-converged field every
!	      great iteration (m15_2p82d, job 12045050).  The descending
!	      refresh passes own the compressing component of the lagged
!	      field; suppress the write at compressing rows instead.
!	      Row 1's direction keys on A_BND under the Eddington boundary
!	      (matching the RHS assembly); interior rows key on Q_NM.
!	      Requires MOM_RF_ENABLED (a refresh actually runs to own the
!	      store) and KN_M_READY (the seeding sweep still writes).
!
!	      edited by Devin, 01-Oct-2026: extended to mixed grids under
!	      MOM_RF_OWNW.  On an alternating schedule the refresh's
!	      descending passes still own the compressing rows, and the
!	      sweep's per-iteration Jacobi write re-injected the error
!	      the refresh just removed (job 12047200: ~2M cells moved
!	      >50% between retained field and next sweep while the
!	      refresh residuals stayed ~1E4).  DOM_EXP is exempted
!	      because no descending pass ever owns compressing rows.
!
	      L_SUPP=MOM_RF_ENABLED .AND. KN_M_READY .AND.
	1        .NOT. MOM_KN_REFRESHING .AND. Q_NM(I) .LT. 0.0_LDP
	2        .AND. (MOM_DOM_CMP .OR.
	3        (MOM_RF_OWNW .AND. .NOT. MOM_DOM_EXP))
	      IF(I .EQ. 1 .AND. OUT_BC_TYPE .LE. 1)L_SUPP=
	1        MOM_RF_ENABLED .AND. KN_M_READY .AND.
	2        .NOT. MOM_KN_REFRESHING .AND. A_BND .LT. 0.0_LDP
	3        .AND. (MOM_DOM_CMP .OR.
	4        (MOM_RF_OWNW .AND. .NOT. MOM_DOM_EXP))
	      IF(L_SUPP)MOM_SW_NSUPP=MOM_SW_NSUPP+1
!
!	      Refresh pass-direction ownership [MOM_RF_OWNW]: a row
!	      whose upwind was stale this pass (compressing on an
!	      ascending pass, expanding on a descending one) took a
!	      Jacobi step; suppressing its write leaves the row at its
!	      last fresh-read iterate, so every stored update is a
!	      Gauss-Seidel step.  Row 1 keys on A_BND under the Eddington
!	      boundary, matching the RHS assembly's upwind selector.
!
	      IF(MOM_KN_REFRESHING .AND. MOM_RF_OWNW)THEN
	        L_SUPP=(KN_M_DESC .NEQV. (Q_NM(I) .LT. 0.0_LDP))
	        IF(I .EQ. 1 .AND. OUT_BC_TYPE .LE. 1)L_SUPP=
	1          (KN_M_DESC .NEQV. (A_BND .LT. 0.0_LDP))
	        IF(L_SUPP)MOM_RF_NSUPP=MOM_RF_NSUPP+1
	      END IF
	      IF(XM_BAD_NM(I) .AND. MOM_RF_KEEPNEG)THEN
!
! Repaired cell: keep the previous iterate so the stored field holds
! only genuine solve outputs; sign-flipped/clamped repairs are not
! fixed-point values and reseed the lagged read with noise.
!
	        CONTINUE
	      ELSE IF(L_SUPP)THEN
	        CONTINUE
	      ELSE
!
!	        During refresh, cells that cycle under the global relax
!	        factor get a reduced per-cell factor (MOM_DAMP_NM); the
!	        fixed point is unchanged, only the map multiplier.
!
	        T5=JNEXT_RELAX
	        IF(MOM_KN_REFRESHING .AND. ALLOCATED(MOM_DAMP_NM))
	1          T5=MOM_DAMP_NM(I,FREQ_INDX)
	        J_NEXT_FRQ(I,FREQ_INDX)=J_NEXT_FRQ(I,FREQ_INDX)+
	1          T5*(XM_NM(I)-J_NEXT_FRQ(I,FREQ_INDX))
	      END IF
	    ELSE
	      J_NEXT_FRQ(I,FREQ_INDX)=XM_NM(I)
	    END IF
	  END DO
	END IF
!
! Record the quarantine mask on the model grid so DO_VAR_CONT can
! zero the claimed dJ response at frozen cells.  Kept from the last
! non-refresh call: the SE Jacobian is built against the sweep solve,
! and refresh re-calls would otherwise overwrite it with pass masks.
!
	IF(.NOT. MOM_KN_REFRESHING)THEN
	  IF(.NOT. ALLOCATED(MOM_CELL_QUAR))ALLOCATE(MOM_CELL_QUAR(ND_SM))
	  IF(SIZE(MOM_CELL_QUAR) .NE. ND_SM)THEN
	    DEALLOCATE(MOM_CELL_QUAR)
	    ALLOCATE(MOM_CELL_QUAR(ND_SM))
	  END IF
	  DO I=1,ND_SM
	    MOM_CELL_QUAR(I)=XM_BAD_NM(J_INDX_NM(I))
	  END DO
!
! Model-grid diagonal of this frequency's moment row.  The unknown is
! R^2*J and the local source piece of the RHS is
! DTAUONQ*SOURCE*R^2 with SOURCE=ETA/CHI.  Holding the off-diagonal
! depths and the lagged upwind fixed,
!   dJ/dETA = DTAUONQ / (TB * CHI)
!   dJ/dCHI = -SOURCE * dJ/dETA
! TB includes PSI = DTAUONQ*|Q|, so a compressing depth (Q<0, upwind
! = J_NEXT) has a smaller local response.  J_NEXT itself is not
! differentiated: it is the outer iteration.  Boundary rows are left
! at zero; they are closures, not local-source rows.
! edited by Sean Brennan (25-Sep-2026)
!
	  IF(.NOT. ALLOCATED(MOM_ALO_DCHI))THEN
	    ALLOCATE(MOM_ALO_DCHI(ND_SM),MOM_ALO_DETA(ND_SM))
	  ELSE IF(SIZE(MOM_ALO_DCHI) .NE. ND_SM)THEN
	    DEALLOCATE(MOM_ALO_DCHI,MOM_ALO_DETA)
	    ALLOCATE(MOM_ALO_DCHI(ND_SM),MOM_ALO_DETA(ND_SM))
	  END IF
	  MOM_ALO_DCHI(1:ND_SM)=0.0_LDP
	  MOM_ALO_DETA(1:ND_SM)=0.0_LDP
	  DO I=1,ND_SM
	    K=J_INDX_NM(I)
	    IF(K .GE. 2 .AND. K .LE. ND_NM-1)THEN
	      IF(TB_NM(K) .GT. 0.0_LDP .AND.
	1         CHI_NM(K) .GT. 0.0_LDP)THEN
	        MOM_ALO_DETA(I)=DTAUONQ_NM(K)/(TB_NM(K)*CHI_NM(K))
	        MOM_ALO_DCHI(I)=-SOURCE_NM(K)*MOM_ALO_DETA(I)
	      END IF
	    END IF
	  END DO
	  MOM_ALO_VALID=.TRUE.
	  MOM_ALO_ND=ND_SM
	END IF
	J_NEXT_VALID=.TRUE.
!
! Record that the moment solver ran this process: SOLVE_FOR_POPS
! then includes the MOM GS state (J_NEXT_FRQ convergence)
! in the restart warm-up gate instead of checking only the FG I_NEXT
! state.
!
	MOM_J_RAN=.TRUE.
!
! Compute H (flux) from J.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=1,ND_NM-1
	    RSQHNU_NM(I)=HU_NM(I)*XM_NM(I+1)-HL_NM(I)*XM_NM(I)+
	1       HS_NM(I)*HUP_NM(I)
	  END DO
	ELSE
	  DO I=1,ND_NM-1
	    RSQHNU_NM(I)=HU_NM(I)*XM_NM(I+1)-HL_NM(I)*XM_NM(I)+
	1       HS_NM(I)*HUP_NM(I)
! -d/dlnu(EPS*J_mid) evaluated at the same depth-upwind node as the
! tridiagonal stencil (outer node for a positive coefficient, inner
! for a negative one) and at the same frequency upwind as the RHS
! assembly: a compressing midpoint reads the redder JUP_NM field.
! The lagged (EPS_PREV) and implicit (EPS) parts take their own
! signs since a midpoint can flip direction between frequencies.
	    IF(EPS_PREV_NM(I) .GE. 0.0_LDP)THEN
	      K=I
	    ELSE
	      K=I+1
	    END IF
	    IF(EPS_PREV_NM(I) .LT. 0.0_LDP)THEN
	      T4=JUP_NM(K)
	    ELSE
	      T4=RSQJNU_PREV_NM(K)
	    END IF
	    RSQHNU_NM(I)=RSQHNU_NM(I)+2.0_LDP*EPS_PREV_NM(I)*T4
	    IF(EPS_NM(I) .GE. 0.0_LDP)THEN
	      RSQHNU_NM(I)=RSQHNU_NM(I)-2.0_LDP*EPS_NM(I)*XM_NM(I)
	    ELSE
	      RSQHNU_NM(I)=RSQHNU_NM(I)-2.0_LDP*EPS_NM(I)*XM_NM(I+1)
	    END IF
	  END DO
	END IF
!
! Enforce |H| <= J (causality check, same as V11).
!
	IF(H_CHK_OPTION .EQ. 'AV_VAL')THEN
	  DO I=1,ND_NM-1
	    T1=(RSQJNU_NM(I)+RSQJNU_NM(I+1))/2.0_LDP
	    IF(RSQHNU_NM(I) .GT. T1)THEN
	      RSQHNU_NM(I)=0.99999_LDP*T1
	    ELSE IF(RSQHNU_NM(I) .LT. -T1)THEN
	      RSQHNU_NM(I)=-0.99999_LDP*T1
	    END IF
	  END DO
	ELSE IF(H_CHK_OPTION .EQ. 'MAX_VAL')THEN
	  DO I=1,ND_NM-1
	    T1=MAX(RSQJNU_NM(I),RSQJNU_NM(I+1))
	    IF(RSQHNU_NM(I) .GT. T1)THEN
	      RSQHNU_NM(I)=0.99999_LDP*T1
	    ELSE IF(RSQHNU_NM(I) .LT. -T1)THEN
	      RSQHNU_NM(I)=-0.99999_LDP*T1
	    END IF
	  END DO
	END IF
!
! The H-clip comparisons do not catch NaN; RSQHNU_NM is read back as
! RSQHNU_PREV_NM by the next frequency, so non-finite entries are
! zeroed before they can propagate.
!
	DO I=1,ND_NM-1
	  IF(.NOT. (ABS(RSQHNU_NM(I)) .LE. HUGE(1.0_LDP)))
	1      RSQHNU_NM(I)=0.0_LDP
	END DO
!
! Update the lagged R^2*H store with the same relaxed-write rule as
! J_NEXT_FRQ: compressing midpoints read slot k+1 as their upwind
! field on the next sweep.  Drift is folded into MOM_GS_MAXCHG only:
! H crosses zero legitimately, so per-cell fraction counts would be
! inflated by sign-change cells and could hold the maturity gate on
! noise; the max-norm still catches a runaway H field.
!
	IF(FREQ_INDX .GE. 1 .AND. FREQ_INDX .LE. NCF_IN)THEN
	  T3=MAX(MAXVAL(ABS(H_NEXT_FRQ(1:ND_NM-1,FREQ_INDX))),
	1    1.0E-28_LDP)
	  DO I=1,ND_NM-1
	    T1=H_NEXT_FRQ(I,FREQ_INDX)
	    T2=RSQHNU_NM(I)
	    IF(ABS(T1) .GT. 1.0E-28_LDP)THEN
	      T4=ABS(T2-T1)/MAX(ABS(T1),1.0E-4_LDP*T3)
	    ELSE IF(ABS(T2) .GT. 1.0E-4_LDP*T3)THEN
	      T4=1.0_LDP
	    ELSE
	      T4=0.0_LDP
	    END IF
!
!	    Non-owned midpoint rows do not write the store this pass --
!	    exclude them from the drift metric for the same reason as
!	    the J loop above.  [MOM_RF_OWNW]
!
	    IF(MOM_KN_REFRESHING .AND. MOM_RF_OWNW .AND.
	1      (KN_M_DESC .NEQV. (GAMH_SIGNED_NM(I) .LT. 0.0_LDP)))
	2      T4=0.0_LDP
	    IF(T4 .GT. MOM_GS_MAXCHG)THEN
	      MOM_GS_MAXCHG=T4
	      MOM_GS_ARGI=I;  MOM_GS_ARGL=FREQ_INDX
	    END IF
	    T5=JNEXT_RELAX
	    IF(MOM_KN_REFRESHING .AND. ALLOCATED(MOM_DAMP_NM))
	1      T5=MOM_DAMP_NM(I,FREQ_INDX)
!
!	    Same suppression as the J_NEXT write above: a compressing
!	    midpoint's H store entry is owned by the descending refresh.
!	    Mixed grids under MOM_RF_OWNW included (see J block): the
!	    alternating refresh still owns the compressing midpoint.
!
	    L_SUPP=MOM_RF_ENABLED .AND. KN_M_READY .AND.
	1        .NOT. MOM_KN_REFRESHING .AND.
	2        GAMH_SIGNED_NM(I) .LT. 0.0_LDP .AND.
	3        (MOM_DOM_CMP .OR.
	4        (MOM_RF_OWNW .AND. .NOT. MOM_DOM_EXP))
	    IF(L_SUPP)MOM_SW_NSUPP=MOM_SW_NSUPP+1
!
!	    Refresh pass ownership on the same midpoint direction
!	    (GAMH_SIGNED) the H upwind read uses.
!
	    IF(MOM_KN_REFRESHING .AND. MOM_RF_OWNW)THEN
	      L_SUPP=(KN_M_DESC .NEQV. (GAMH_SIGNED_NM(I) .LT. 0.0_LDP))
!
!	      H-row ownership for the composite best-field merge -- the
!	      midpoint sign is the H upwind selector, and it can differ
!	      from Q_NM(I) at sigma~0 transition cells.
!
	      MOM_OWN_CMPH(I,FREQ_INDX)=
	1          (GAMH_SIGNED_NM(I) .LT. 0.0_LDP)
	      IF(L_SUPP)MOM_RF_NSUPP=MOM_RF_NSUPP+1
	    END IF
	    IF(.NOT. L_SUPP)H_NEXT_FRQ(I,FREQ_INDX)=T1+T5*(T2-T1)
	  END DO
	  H_NEXT_FRQ(ND_NM,FREQ_INDX)=0.0_LDP
	END IF
!
! Regrid J and H onto the small (input) grid.
!
	DO I=1,ND_SM
	  K=J_INDX_NM(I)
	  JNU(I)=RSQJNU_NM(K)/R_SM(I)/R_SM(I)
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX_NM(I)
	  RSQHNU_SM(I)=RSQHNU_NM(K)
	END DO
!
! Save Eddington factors for next frequency.
!
	K_ON_J_SAVE(1:ND_NM)=K_ON_J(1:ND_NM)
	NMID_ON_HMID_SAVE(1:ND_NM)=NMID_ON_HMID(1:ND_NM)
	NMID_ON_J_SAVE(1:ND_NM)=NMID_ON_J(1:ND_NM)
	HBC_SAVE=HBC
	IN_HBC_SAVE=IN_HBC
	NBC_SAVE=NBC
!
! Capture the per-frequency inputs MOM_KNOP_REFRESH_V1 needs to
! re-drive this routine.  Slot s holds frequency s (same convention
! as J_NEXT_FRQ and the FG lagged store).  The capture is skipped
! during refresh re-calls: the store already holds the values being
! replayed, and rewriting it would corrupt later slots on a
! descending pass.
!
	IF(ALLOCATED(KN_M_CHISM) .AND. .NOT. MOM_KN_REFRESHING)THEN
	  K=FREQ_INDX
	  IF(K .GE. 1 .AND. K .LE. KN_M_NCF)THEN
	    KN_M_CHISM(1:ND_SM,K)=CHI_SM(1:ND_SM)
	    KN_M_ETASM(1:ND_SM,K)=ETA_SM(1:ND_SM)
	    KN_M_ESECSM(1:ND_SM,K)=ESEC_SM(1:ND_SM)
	    KN_MKONJ(1:ND_NM,K)=K_ON_J(1:ND_NM)
	    KN_MNMJ(1:ND_NM,K)=NMID_ON_J(1:ND_NM)
	    KN_MNMH(1:ND_NM,K)=NMID_ON_HMID(1:ND_NM)
	    KN_MHBC(K)=HBC;  KN_MNBC(K)=NBC;  KN_MINHBC(K)=IN_HBC
	    KN_MJPIB(K)=JPLUS_IB;  KN_MJMIB(K)=JMIN_IB
	    KN_MHIB(K)=HNU_AT_IB
	    KN_MDBB(K)=DBB;  KN_MIC(K)=IC
	    KN_MDLN(K)=dLOG_NU
	    KN_MNUNXT(K)=MOM_NU_NEXT
	    KN_MFREQ(K)=FREQ
	    KN_M_OUTBC(K)=OUT_BC_TYPE
	    IF(HAS_COMPRESSING)KN_M_ANYCMP=.TRUE.
	    IF(K .EQ. KN_M_NCF)KN_M_READY=.TRUE.
	  END IF
	END IF
!
	RETURN
	END
!
! ==============================================================================
!  Gauss-Seidel J_NEXT_FRQ restart/checkpoint routines (added 16-Sep-2026).
!
!  Same design as the FG_J_CMF_NONMONO_V1 I_NEXT checkpoint: J_NEXT_FRQ holds
!  a complete GS generation only between frequency sweeps; it is checkpointed
!  at every SCR_RITE_V2 call (post-sweep by construction), tagged with the
!  NITSF written to the POINT files, and accepted on restart only when the
!  tag equals the NITSF SCR_READ_V2 restored, the dims match, and the R/V
!  grid fingerprint is identical.
!
!  File MOMJNXT_GS (unformatted sequential):
!    rec 1: MAGIC, VERSION, TAG
!    rec 2: ND_NM, NCF_IN, ND_SM
!    rec 3: R(1), R(ND), V(1), V(ND)
!    rec 4: R_NM, R_PNT_NM expanded-grid fingerprints
!    rec 5: MOM_GS_CHG_FILE, MOM_GS_CAP_FILE sweep metrics
!    rec 6: J_NEXT_FRQ (ND_NM*NCF_IN values)
!    rec 7: H_NEXT_FRQ (VER 4 only; VER 3 files end after rec 6)
!    rec 8: END_MAGIC
!  Any failure leaves J_NEXT_FRQ zeroed, J_NEXT_VALID .FALSE. and
!  MOM_J_GS_OK .FALSE. -> LAMBDA warm-up fallback.
! ==============================================================================
!
	SUBROUTINE SET_MOM_GS_TAG(N)
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_GS_TAG_EXP
	IMPLICIT NONE
	INTEGER N
	MOM_GS_TAG_EXP=N
	RETURN
	END
!
	SUBROUTINE DUMP_MOM_JNEXT_GS(NITSF,R,V,SIGMA,ND,LU)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1
	USE IFPORT
	IMPLICIT NONE
	INTEGER NITSF,ND,LU
	REAL(KIND=LDP) R(ND),V(ND),SIGMA(ND)
	INTEGER IOS,I,ND_NM_F,LUER,ERROR_LU
	CHARACTER(LEN=7) SW_ST,SV_ST
	INTEGER, PARAMETER :: MAGIC=20260916, VER=4
	LOGICAL REN_OK
	EXTERNAL ERROR_LU
!
	IF(.NOT. ALLOCATED(J_NEXT_FRQ))RETURN
!
! Report the per-sweep Gauss-Seidel convergence of the J_NEXT frequency
! coupling (Knop et al. 2009 Paper II pathway) to OUTGEN.
!
	LUER=ERROR_LU()
!
! Tabulated sweep report.  The status tag mirrors the maturity gate
! tiers in SOLVE_FOR_POPS: [BAD] = catastrophic tier (>10% capped,
! >1% drifting >50%, max change >1E6, >10% negJ or quarantined),
! [WARN] = marginal tier (above the strict thresholds), else [OK].
!
	SW_ST='[OK]   '
	IF(MOM_INEXT_CAPPED .GT. MAX(50_8,MOM_GS_NTOT/10) .OR.
	1   MOM_GS_NBAD50 .GT. MAX(50_8,MOM_GS_NTOT/100) .OR.
	2   MOM_GS_MAXCHG .GT. 1.0E+6_LDP .OR.
	3   MOM_NEGJ_LRG .GT. MOM_GS_NTOT/10 .OR.
	4   MOM_GS_NQUAR .GT. MOM_GS_NTOT/10)THEN
	  SW_ST='[BAD]  '
	ELSE IF(MOM_GS_NBAD50 .GT. MAX(50_8,MOM_GS_NTOT/1000) .OR.
	1   MOM_GS_MAXCHG .GT. 1.0E+1_LDP .OR.
	2   MOM_INEXT_CAPPED .GT. MAX(50_8,MOM_GS_NTOT/1000) .OR.
	3   MOM_GS_NBAD5 .GT. MOM_GS_NTOT/10 .OR.
	4   MOM_BND_HBC_FL .GT. 0)THEN
	  SW_ST='[WARN] '
	END IF
	SV_ST='[OK]   '
	IF(MOM_N_DGTSV .GT. 0 .OR. MOM_N_NONF .GT. 0 .OR.
	1   MOM_TB_NONPOS .GT. 0)SV_ST='[BAD]  '
!
	WRITE(LUER,'(/,A)')' ----------------------- MOM_J_CMF_NONMONO_V1: GS sweep ------------------------'
	WRITE(LUER,'(A)')'   cols: d = max rel. change in lagged'//
	1    ' J_NEXT; counts are cells this sweep'
	WRITE(LUER,'(4X,A10,2X,ES10.3,1X,A,I10,A,I10,1X,A,I10,2X,A)')
	1    'drift',MOM_GS_MAXCHG,'n50=',MOM_GS_NBAD50,'/',
	2    MOM_GS_NTOT,'n5=',MOM_GS_NBAD5,SW_ST
	WRITE(LUER,'(4X,A10,1X,A,I10,1X,A,I10,1X,A,I10,1X,A,I10,2X,A)')
	1    'repairs','negJ=',MOM_NEGJ_LRG,'tiny=',MOM_N_NEGJ_T,
	2    'quar=',MOM_GS_NQUAR,'pquar=',MOM_GS_NPQUAR,
	3    'significant/roundoff negJ; quarantined cells'
	WRITE(LUER,'(4X,A10,1X,A,I6,1X,A,I10,1X,A,I10,1X,A,I10,2X,A)')
	1    'solver','DGTSV=',MOM_N_DGTSV,'nonfinite=',
	2    MOM_N_NONF,'runaway=',MOM_N_RUNAWAY,'TB<=0=',
	3    MOM_TB_NONPOS,SV_ST
	WRITE(LUER,'(4X,A10,1X,A,I10,1X,A,I10,1X,A,ES10.3,1X,A,I10)')
	1    'structure','nonMM=',MOM_NONMM_CNT,'negEPS=',
	2    MOM_EPS_FWD,'max|EPS|/HU=',MOM_EPS_RAT,'eps_cap=',
	3    MOM_EPS_CAP_SW
	WRITE(LUER,'(4X,A10,1X,A,I10,1X,A,ES10.3)')
	1    'lagbnd','scaled=',MOM_LAGB_CNT,'mins=',MOM_LAGB_MIN
	WRITE(LUER,'(4X,A10,1X,A,ES10.3,A,I5,A,I5,A,1X,A,ES10.3,A,I5,A,I5,A)')
	1    'matrix','minDR(s)=',MOM_MIN_DR,' at(',MOM_MIN_DR_I,',',
	2    MOM_MIN_DR_L,')','minDR|a|=',MOM_MIN_DRA,' at(',
	3    MOM_MIN_DRA_I,',',MOM_MIN_DRA_L,')'
	WRITE(LUER,'(4X,A10,1X,A,ES10.3,A,I5,A,I5,A,1X,A,ES10.3)')
	1    'eq resid','max=',MOM_EQ_RES,' at(',MOM_EQ_RES_I,',',
	2    MOM_EQ_RES_L,')','max ESEC/CHI-1=',MOM_MAX_COH
	WRITE(LUER,'(4X,A10,1X,A,ES10.3,1X,A,ES10.3,1X,A,I10)')
	1    'params','max|GAM|=',MOM_MAX_GAM,'max|sigT|=',
	2    MOM_MAX_SIGT,'HBC floors=',MOM_BND_HBC_FL
	WRITE(LUER,'(4X,A10,1X,A,I6,A,I6,A,1X,A,I10)')
	1    'misc','argmax(d,f)=(',MOM_GS_ARGI,',',MOM_GS_ARGL,')',
	2    'NCF=',MOM_GS_NCF
	WRITE(LUER,'(4X,A10,1X,A,A4,1X,A,ES10.3,1X,A,L1,1X,A,L1,'//
	1    '1X,A,I10)')
	2    'ordering','dir=',MOM_RF_DIR,'cmp_frac=',MOM_SW_FCMP,
	3    'dom_cmp=',MOM_DOM_CMP,'dom_exp=',MOM_DOM_EXP,
	4    'nsupp=',MOM_SW_NSUPP
	IF(MOM_RF_RAN)THEN
	  WRITE(LUER,'(4X,A10,1X,A,L2,1X,A,ES10.3,1X,A,I10,1X,A,I10)')
	1    'refresh','conv=',MOM_RF_CONV,'delta=',MOM_RF_DCHG,
	2    'n5=',MOM_RF_N5,'n50=',MOM_RF_N50
	END IF
	WRITE(LUER,'(A)')' -------------------------------------------------------------------------------'
!
! Format version 2 (16-Sep-2026): J_NEXT_FRQ lives on the EXPANDED
! depth grid R_NM/R_PNT_NM, whose construction depends on the internal
! velocity profile -- two models can share ND_NM and R/V endpoints yet
! have different expanded grids, so the fingerprint is the full R_NM
! and R_PNT_NM vectors plus the input-grid R/V endpoints.  Written to
! MOMJNXT_GS.tmp then atomically renamed so a kill mid-write cannot
! destroy the previous valid checkpoint.
!
	ND_NM_F=SIZE(J_NEXT_FRQ,1)
	OPEN(UNIT=LU,FILE='MOMJNXT_GS.tmp',FORM='UNFORMATTED',
	1     STATUS='REPLACE',IOSTAT=IOS)
	IF(IOS .NE. 0)RETURN
	WRITE(LU,IOSTAT=IOS)MAGIC,VER,NITSF
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)ND_NM_F,MOM_GS_NCF,ND
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)R(1),R(ND),V(1),V(ND)
	IF(IOS .EQ. 0 .AND. ALLOCATED(R_NM) .AND.
	1   ALLOCATED(R_PNT_NM))THEN
	  WRITE(LU,IOSTAT=IOS)(R_NM(I),I=1,ND_NM_F)
	  IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)(R_PNT_NM(I),I=1,ND_NM_F)
	ELSE
	  IOS=-1
	END IF
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)MOM_GS_MAXCHG,
	1                               MOM_INEXT_CAPPED
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)J_NEXT_FRQ
! VER 4: the lagged R^2*H partner store travels with J_NEXT_FRQ so a
! restored generation keeps both upwind fields consistent.
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)H_NEXT_FRQ
	IF(IOS .EQ. 0)WRITE(LU,IOSTAT=IOS)-MAGIC
	CLOSE(LU)
!
	IF(IOS .EQ. 0)THEN
	  REN_OK=RENAMEFILEQQ('MOMJNXT_GS.tmp','MOMJNXT_GS')
	END IF
	RETURN
	END
!
	SUBROUTINE RESTORE_MOM_JNEXT_GS(R,V,ND)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1
	IMPLICIT NONE
	INTEGER ND
	REAL(KIND=LDP) R(ND),V(ND)
	INTEGER IOS,LUER,ERROR_LU
	INTEGER MAG,VER,TAG
	INTEGER ND_NM_F,NCF_F,ND_F,ENDM,I
	INTEGER, PARAMETER :: MAGIC=20260916
	REAL(KIND=LDP) R1_F,RND_F,V1_F,VND_F
	REAL(KIND=LDP), ALLOCATABLE :: RNM_F(:)
	INTEGER, ALLOCATABLE :: RPNT_F(:)
	EXTERNAL ERROR_LU
	LOGICAL FILE_THERE
!
	LUER=ERROR_LU()
	MOM_J_GS_OK=.FALSE.
	MOM_GS_CHG_FILE=-1.0_LDP
	MOM_GS_CAP_FILE=-1
!
	WRITE(LUER,'(A)')'GS RESTART STATE (MOM_J_CMF_NONMONO_V1):'
!
! Hard invariant: NEWMOD (or any state where SCR_READ_V2 did not
! return a positive NITSF) never restores GS state -- exit before the
! file is even opened, so a stale file from another model cannot be
! accepted on a coincidental tag/dimension match.
!
	IF(MOM_GS_TAG_EXP .LE. 0)THEN
	  WRITE(LUER,'(A)')'  cold start (NITSF<=0): no GS restore'
	  WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (J_NEXT=0)'
	  RETURN
	END IF
!
	INQUIRE(FILE='MOMJNXT_GS',EXIST=FILE_THERE)
	IF(.NOT. FILE_THERE)THEN
	  WRITE(LUER,'(A)')'  J_NEXT checkpoint present: NO'
	  WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (J_NEXT=0)'
	  RETURN
	END IF
	WRITE(LUER,'(A)')'  J_NEXT checkpoint present: YES'
!
	OPEN(UNIT=292,FILE='MOMJNXT_GS',FORM='UNFORMATTED',
	1     STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(292,IOSTAT=IOS)MAG,VER,TAG
! VER 3 files carry J_NEXT_FRQ only; VER 4 adds the H_NEXT_FRQ record.
! A VER 3 restore leaves the H store cold -- it is rebuilt from the
! first sweep, cheaper than discarding the whole J generation.
	IF(IOS .NE. 0 .OR. MAG .NE. MAGIC .OR.
	1   (VER .NE. 3 .AND. VER .NE. 4))GOTO 900
	READ(292,IOSTAT=IOS)ND_NM_F,NCF_F,ND_F
	IF(IOS .NE. 0)GOTO 900
	READ(292,IOSTAT=IOS)R1_F,RND_F,V1_F,VND_F
	IF(IOS .NE. 0)GOTO 900
!
	WRITE(LUER,'(A,I8)')'  GS tag (checkpoint iter): ',TAG
	WRITE(LUER,'(A,I8)')'  expected tag (NITSF):     ',MOM_GS_TAG_EXP
	IF(TAG .NE. MOM_GS_TAG_EXP)THEN
	  WRITE(LUER,'(A)')'  J_NEXT valid: NO (stale/foreign checkpoint)'
	  GOTO 800
	END IF
	IF(ND_NM_F .NE. SIZE(J_NEXT_FRQ,1) .OR.
	1   NCF_F .NE. SIZE(J_NEXT_FRQ,2) .OR. ND_F .NE. ND)THEN
	  WRITE(LUER,'(A)')'  J_NEXT valid: NO (dimension mismatch)'
	  GOTO 800
	END IF
	IF(R1_F .NE. R(1) .OR. RND_F .NE. R(ND) .OR.
	1   V1_F .NE. V(1) .OR. VND_F .NE. V(ND))THEN
	  WRITE(LUER,'(A)')'  J_NEXT valid: NO (grid fingerprint mismatch)'
	  GOTO 800
	END IF
!
! Full expanded-grid fingerprint: every element of R_NM and the
! R_PNT_NM depth map must be identical -- the expanded grid depends on
! the internal velocity profile, so endpoint matching is insufficient
! to guarantee J_NEXT_FRQ(I,K) refers to the same physical depth.
! Each unformatted READ consumes a whole record, so the file vectors
! are read into matching-shape locals, not elementwise.
!
	ALLOCATE(RNM_F(ND_NM_F),RPNT_F(ND_NM_F),STAT=IOS)
	IF(IOS .NE. 0)GOTO 900
	READ(292,IOSTAT=IOS)RNM_F
	IF(IOS .NE. 0)GOTO 900
	READ(292,IOSTAT=IOS)RPNT_F
	IF(IOS .NE. 0)GOTO 900
	IF(ANY(RNM_F .NE. R_NM) .OR.
	1   ANY(RPNT_F .NE. R_PNT_NM))GOTO 810
	DEALLOCATE(RNM_F,RPNT_F)
!
	READ(292,IOSTAT=IOS)MOM_GS_CHG_FILE,MOM_GS_CAP_FILE
	IF(IOS .NE. 0)GOTO 900
!
	READ(292,IOSTAT=IOS)J_NEXT_FRQ
	IF(IOS .NE. 0)GOTO 900
	IF(VER .EQ. 4)THEN
	  READ(292,IOSTAT=IOS)H_NEXT_FRQ
	  IF(IOS .NE. 0)GOTO 900
	END IF
	READ(292,IOSTAT=IOS)ENDM
	IF(IOS .NE. 0 .OR. ENDM .NE. -MAGIC)GOTO 900
!
! Value validation (audit 17-Sep-2026): a structurally intact file can
! still carry a corrupted field -- the convergence metric must be
! finite and non-negative, the event count non-negative, and the
! restored R^2*J finite and non-negative everywhere.  Any violation
! rejects the whole generation: the next compressing-frequency solve
! would otherwise use the corrupted values as explicit source terms.
!
	IF(.NOT.(ABS(MOM_GS_CHG_FILE) .LE. HUGE(1.0_LDP)) .OR.
	1   MOM_GS_CHG_FILE .LT. 0.0_LDP .OR.
	1   MOM_GS_CAP_FILE .LT. 0)GOTO 820
	IF(ANY(.NOT.(ABS(J_NEXT_FRQ) .LE. HUGE(1.0_LDP))) .OR.
	1   ANY(J_NEXT_FRQ .LT. 0.0_LDP))GOTO 820
! H is signed; only finiteness is required.
	IF(ANY(.NOT.(ABS(H_NEXT_FRQ) .LE. HUGE(1.0_LDP))))GOTO 820
	CLOSE(292)
!
	MOM_J_GS_OK=.TRUE.
	J_NEXT_VALID=.TRUE.
	WRITE(LUER,'(A)')'  J_NEXT valid: YES'
	WRITE(LUER,'(A,I8,A)')
	1 '  RESTORED COMPLETE J_NEXT GENERATION (tag=',TAG,')'
	WRITE(LUER,'(A,ES10.3,A,I6)')
	1 '  last-sweep max change=',MOM_GS_CHG_FILE,
	1 '  non-finite events=',MOM_GS_CAP_FILE
	RETURN
!
900	CONTINUE
	WRITE(LUER,'(A)')'  J_NEXT valid: NO (read/format failure)'
	GOTO 800
810	CONTINUE
	WRITE(LUER,'(A)')'  J_NEXT valid: NO (expanded-grid mismatch)'
	GOTO 800
820	CONTINUE
	WRITE(LUER,'(A)')'  J_NEXT valid: NO (non-finite/negative data)'
800	CONTINUE
	INQUIRE(UNIT=292,OPENED=FILE_THERE)
	IF(FILE_THERE)CLOSE(292)
	MOM_J_GS_OK=.FALSE.
	J_NEXT_VALID=.FALSE.
	MOM_GS_CHG_FILE=-1.0_LDP
	MOM_GS_CAP_FILE=-1
	IF(ALLOCATED(J_NEXT_FRQ))J_NEXT_FRQ=0.0_LDP
	IF(ALLOCATED(H_NEXT_FRQ))H_NEXT_FRQ=0.0_LDP
	WRITE(LUER,'(A)')'  FALLING BACK TO LAMBDA REBUILD (J_NEXT=0)'
	RETURN
	END
!
! ==============================================================================
!  GET_MOM_GS_MATURITY (added 19-Sep-2026, dev7 port).
!
!  Return the moment-solver GS state for the N_LAM_RWARM warm-up decision
!  in SOLVE_FOR_POPS.  Kept as a subroutine (not a module USE) because
!  new_main/mod_subs compiles before plane/ in the library build order,
!  so module visibility is not available there.
!
!    RAN      : MOM_J_CMF_NONMONO_V1 ran in this process
!    GS_OK    : a validated MOMJNXT_GS generation was restored
!    CHG_FILE : last-sweep max relative change recorded in the file
!               (-1 when no checkpoint was read)
!    CAP_FILE : capped/non-finite write count recorded in the file
!    MAXCHG   : max relative change written into J_NEXT_FRQ this sweep
!    CAPPED   : non-finite or runaway XM_NM rejections this sweep
!    NEGJ     : significant negative-J resets this sweep
!    NBAD5    : cells with relative J_NEXT_FRQ change above 5 percent
!    NBAD50   : cells with relative change above 50 percent
!    NTOT     : total cells the metric covered (ND_NM*NCF)
!    NQUAR    : cells quarantined this sweep under MOM_RF_KEEPNEG
!    NPQUAR   : quarantined cells at MOM_QUAR_PERSIST sweeps or older
!               (the structurally failing subset of NQUAR)
!    RF_RAN   : MOM_KNOP_REFRESH_V1 executed this great iteration
!    RF_CONV  : the refresh converged (retained field is the map's
!               fixed point)
!    RF_DCHG  : max relative change of the retained J_NEXT_FRQ
!               between pre- and post-refresh
!    RF_N5/N50: cells whose retained field moved >5% / >50%
! ==============================================================================
!
	SUBROUTINE GET_MOM_GS_MATURITY(RAN,GS_OK,CHG_FILE,CAP_FILE,
	1                             MAXCHG,CAPPED,NEGJ,NBAD5,NBAD50,
	2                             NTOT,NQUAR,SW_MAXCHG,SW_NBAD5,
	3                             SW_NBAD50,SW_NTOT,SW_NEGJ,SW_NQUAR,
	4                             SW_CAPPED,SW_BNDFL,RF_RAN,RF_CONV,
	5                             RF_DCHG,RF_N5,RF_N50,NPQUAR,SW_NHUGE)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     MOM_J_RAN,MOM_J_GS_OK,MOM_GS_CHG_FILE,MOM_GS_CAP_FILE,
	2     MOM_GS_MAXCHG,MOM_INEXT_CAPPED,MOM_NEGJ_LRG,
	3     MOM_GS_NBAD5,MOM_GS_NBAD50,MOM_GS_NTOT,MOM_GS_NQUAR,
	4     MOM_SW_MAXCHG,MOM_SW_NBAD5,MOM_SW_NBAD50,MOM_SW_NTOT,
	5     MOM_SW_NEGJ,MOM_SW_NQUAR,MOM_SW_CAPPED,MOM_BND_HBC_FL,
	6     MOM_RF_RAN,MOM_RF_CONV,MOM_RF_DCHG,MOM_RF_N5,MOM_RF_N50,
	7     MOM_GS_NPQUAR,MOM_GS_NHUGE,MOM_SW_NHUGE
	IMPLICIT NONE
	LOGICAL RAN,GS_OK,RF_RAN,RF_CONV
	REAL(KIND=LDP) CHG_FILE,MAXCHG,SW_MAXCHG,RF_DCHG
	INTEGER*8 CAP_FILE,CAPPED,NEGJ,NBAD5,NBAD50,NTOT,NQUAR,NPQUAR
	INTEGER*8 SW_NBAD5,SW_NBAD50,SW_NTOT,SW_NEGJ,SW_NQUAR,SW_CAPPED
	INTEGER*8 SW_BNDFL,SW_NHUGE
	INTEGER*8 RF_N5,RF_N50
	RF_RAN=MOM_RF_RAN
	RF_CONV=MOM_RF_CONV
	RF_DCHG=MOM_RF_DCHG
	RF_N5=MOM_RF_N5
	RF_N50=MOM_RF_N50
	RAN=MOM_J_RAN
	GS_OK=MOM_J_GS_OK
	CHG_FILE=MOM_GS_CHG_FILE
	CAP_FILE=MOM_GS_CAP_FILE
	MAXCHG=MOM_GS_MAXCHG
	CAPPED=MOM_INEXT_CAPPED
	NEGJ=MOM_NEGJ_LRG
	NBAD5=MOM_GS_NBAD5
	NBAD50=MOM_GS_NBAD50
	NTOT=MOM_GS_NTOT
	NQUAR=MOM_GS_NQUAR
	NPQUAR=MOM_GS_NPQUAR
	SW_BNDFL=MOM_BND_HBC_FL
!
! When no refresh ran this iteration the sweep-level slots stay at
! their last values; report the live counters so the gate always sees
! a sweep-drift estimate.
!
	IF(MOM_SW_NTOT .GT. 0)THEN
	  SW_MAXCHG=MOM_SW_MAXCHG
	  SW_NBAD5=MOM_SW_NBAD5
	  SW_NBAD50=MOM_SW_NBAD50
	  SW_NTOT=MOM_SW_NTOT
	  SW_NEGJ=MOM_SW_NEGJ
	  SW_NQUAR=MOM_SW_NQUAR
	  SW_CAPPED=MOM_SW_CAPPED
	  SW_NHUGE=MOM_SW_NHUGE
	ELSE
	  SW_MAXCHG=MOM_GS_MAXCHG
	  SW_NBAD5=MOM_GS_NBAD5
	  SW_NBAD50=MOM_GS_NBAD50
	  SW_NTOT=MOM_GS_NTOT
	  SW_NEGJ=MOM_NEGJ_LRG
	  SW_NQUAR=MOM_GS_NQUAR
	  SW_CAPPED=MOM_INEXT_CAPPED
	  SW_NHUGE=MOM_GS_NHUGE
	END IF
	RETURN
	END
!
! ==============================================================================
! SET_MOM_NU_NEXT passes nu(k+1) to the solver each call so the
! compressing depths can use the forward spacing ln(nu_k/nu_{k+1})
! (asymmetric-spacing fix, Oct-2026; same role as SET_FG_NU_NEXT in
! FG_J_CMF_NONMONO_V1).  MOM_NU_NEXT<=0 or >= FREQ falls back to the
! backward spacing dLOG_NU.
! ==============================================================================
!
	SUBROUTINE SET_MOM_NU_NEXT(NU_NEXT)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_NU_NEXT
	IMPLICIT NONE
	REAL(KIND=LDP) NU_NEXT
	MOM_NU_NEXT=NU_NEXT
	RETURN
	END
!
! SET_MOM_RF_OPTS: runtime wiring for the MOM refresh toggles
! MOM_RF_GS_DESC (live GS read on descending passes) and
! MOM_RF_KEEPNEG (quarantine repaired cells).  Called once from
! RD_CONTROL_VARIABLES after the VADAT keys are read.
!
	SUBROUTINE SET_MOM_RF_OPTS(GS_DESC,KEEPNEG)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     MOM_RF_GS_DESC,MOM_RF_KEEPNEG
	IMPLICIT NONE
	LOGICAL GS_DESC,KEEPNEG
	MOM_RF_GS_DESC=GS_DESC
	MOM_RF_KEEPNEG=KEEPNEG
	RETURN
	END
!
! SET_MOM_FREQ_MATCH / SET_MOM_RF_AA: runtime wiring for
! DO_NONMONO_FREQ_MATCH and MOM_RF_AA, called once from
! RD_CONTROL_VARIABLES after the VADAT keys are read.
!
	SUBROUTINE SET_MOM_FREQ_MATCH(DO_MATCH)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_FREQ_MATCH
	IMPLICIT NONE
	LOGICAL DO_MATCH
	MOM_FREQ_MATCH=DO_MATCH
	RETURN
	END
!
	SUBROUTINE SET_MOM_RF_AA(DO_AA)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_RF_AA
	IMPLICIT NONE
	LOGICAL DO_AA
	MOM_RF_AA=DO_AA
	RETURN
	END
!
! SET_MOM_RF_BULK: runtime wiring for MOM_RF_BULK_EXIT, called once
! from RD_CONTROL_VARIABLES after the VADAT key is read.
!
	SUBROUTINE SET_MOM_RF_BULK(DO_BULK)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_RF_BULK_EXIT
	IMPLICIT NONE
	LOGICAL DO_BULK
	MOM_RF_BULK_EXIT=DO_BULK
	RETURN
	END
!
! SET_MOM_RF_ACT: runtime wiring for MOM_RF_ACTIVE, called once
! from RD_CONTROL_VARIABLES after the VADAT key is read.
!
	SUBROUTINE SET_MOM_RF_ACT(DO_ACT)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_RF_ACTIVE
	IMPLICIT NONE
	LOGICAL DO_ACT
	MOM_RF_ACTIVE=DO_ACT
	RETURN
	END
!
! SET_MOM_RF_DIR_OPT: runtime wiring for MOM_RF_DIR (refresh pass
! ordering: AUTO/ALT/ASC/DESC), MOM_EPS_LIM (lagged-EPS positivity
! limiter) and MOM_RF_ENABLED (DO_KNOP_REFRESH -- the dominant-mode
! write suppression is only safe while a refresh actually owns the
! compressing store).  Called once from RD_CONTROL_VARIABLES after
! the VADAT keys are read; an unrecognized DIR falls back to AUTO.
!
	SUBROUTINE SET_MOM_RF_DIR_OPT(DIR,EPSLIM,RF_EN,OWNW)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     MOM_RF_DIR,MOM_EPS_LIM,MOM_RF_ENABLED,MOM_RF_OWNW
	IMPLICIT NONE
	CHARACTER(LEN=*) DIR
	LOGICAL EPSLIM,RF_EN,OWNW
	IF(DIR .EQ. 'AUTO' .OR. DIR .EQ. 'ALT' .OR.
	1   DIR .EQ. 'ASC' .OR. DIR .EQ. 'DESC')MOM_RF_DIR=DIR
	MOM_EPS_LIM=EPSLIM
	MOM_RF_ENABLED=RF_EN
	MOM_RF_OWNW=OWNW
	RETURN
	END
!
! SET_MOM_LAG_OPTS: runtime wiring for MOM_GAM_PIN, MOM_LAG_BND and
! MOM_ABS_FLR, called once from RD_CONTROL_VARIABLES after the VADAT
! keys are read.
!
	SUBROUTINE SET_MOM_LAG_OPTS(GAMPIN,LAGBND,ABSFLR,HIMPL)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     MOM_GAM_PIN,MOM_LAG_BND,MOM_ABS_FLR,MOM_H_IMPL
	IMPLICIT NONE
	LOGICAL GAMPIN,LAGBND,HIMPL
	REAL(KIND=LDP) ABSFLR
	MOM_GAM_PIN=GAMPIN
	MOM_LAG_BND=LAGBND
	MOM_ABS_FLR=MAX(ABSFLR,0.0_LDP)
	MOM_H_IMPL=HIMPL
	RETURN
	END
!
! SET_MOM_FAC_RLX: Eddington-factor under-relaxation between great
! iterations (VADAT MOM_FAC_RLX, default 1.0 = raw formal snap).
! Bounds mirror JNEXT_RELAX.
!
	SUBROUTINE SET_MOM_FAC_RLX(RLX)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_FAC_RLX
	IMPLICIT NONE
	REAL(KIND=LDP) RLX
	MOM_FAC_RLX=MIN(MAX(RLX,1.0E-3_LDP),1.0_LDP)
	RETURN
	END
!
! SET_MOM_KN_OPTS: runtime wiring for the previously compile-time
! MOM refresh knobs KN_M_TOL, KN_M_MAXPASS, JNEXT_RELAX and the AA
! history depth MOM_AA_M.  Bounds are defensive: a zero tolerance
! or an AA depth that outgrows the local Gram workspaces would be
! silently wrong rather than obviously wrong.
!
	SUBROUTINE SET_MOM_KN_OPTS(TOL,MAXP,RLX,AA_M,AA_OUT,AA_OUT_M)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     KN_M_TOL,KN_M_MAXPASS,JNEXT_RELAX,MOM_AA_M,
	2     MOM_AA_OUT,MOM_AA_OUT_M
	IMPLICIT NONE
	REAL(KIND=LDP) TOL,RLX
	INTEGER MAXP,AA_M,AA_OUT_M
	LOGICAL AA_OUT
	KN_M_TOL=MAX(TOL,1.0E-12_LDP)
	KN_M_MAXPASS=MAX(MAXP,1)
	JNEXT_RELAX=MIN(MAX(RLX,1.0E-3_LDP),1.0_LDP)
	MOM_AA_M=MAX(1,MIN(AA_M,16))
	MOM_AA_OUT=AA_OUT
	MOM_AA_OUT_M=MAX(1,MIN(AA_OUT_M,16))
	RETURN
	END
!
! GET_MOM_CELL_QUAR: return the quarantine mask (XM_BAD mapped to the
! SM grid) recorded by the last non-refresh MOM solve, for
! DO_VAR_CONT's unconditional ALO diagonal suppression at frozen
! cells.  VALID=F when the solver has not yet produced a mask on an
! N-point grid.
!
	SUBROUTINE GET_MOM_CELL_QUAR(Q,N,VALID)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY: MOM_CELL_QUAR
	IMPLICIT NONE
	INTEGER N
	LOGICAL Q(N),VALID
	VALID=.FALSE.
	Q(1:N)=.FALSE.
	IF(ALLOCATED(MOM_CELL_QUAR))THEN
	  IF(SIZE(MOM_CELL_QUAR) .EQ. N)THEN
	    Q(1:N)=MOM_CELL_QUAR(1:N)
	    VALID=.TRUE.
	  END IF
	END IF
	RETURN
	END
!
! GET_MOM_ALO_DIAG: diagonal dJ/dCHI, dJ/dETA of the moment row just
! solved, on the model grid.  VALID=F until a non-refresh solve has
! stored them, and across a grid-size change.
! edited by Sean Brennan (25-Sep-2026)
!
	SUBROUTINE GET_MOM_ALO_DIAG(DJDCHI,DJDETA,N,VALID)
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1, ONLY:
	1     MOM_ALO_DCHI,MOM_ALO_DETA,MOM_ALO_VALID,MOM_ALO_ND
	IMPLICIT NONE
	INTEGER N
	REAL(KIND=LDP) DJDCHI(N),DJDETA(N)
	LOGICAL VALID
	VALID=.FALSE.
	IF(ALLOCATED(MOM_ALO_DCHI) .AND. ALLOCATED(MOM_ALO_DETA))THEN
	  IF(MOM_ALO_VALID .AND. MOM_ALO_ND .EQ. N)THEN
	    DJDCHI(1:N)=MOM_ALO_DCHI(1:N)
	    DJDETA(1:N)=MOM_ALO_DETA(1:N)
	    VALID=.TRUE.
	  END IF
	END IF
	RETURN
	END
!
! ==============================================================================
! MOM_KNOP_REFRESH_V1 (Oct-2026)
!
! Post-sweep refresh of the lagged compressing-branch wavelength
! coupling in the moment solver -- the MOM analogue of
! FG_KNOP_REFRESH_V1.  In the normal sweep, compressing depths read
! J_NEXT_FRQ(:,k+1) written during the PREVIOUS great iteration, so
! red-to-blue information advances one frequency bin per iteration and
! the GS field cannot mature on strongly non-monotonic models.
!
! This routine re-drives MOM_J_CMF_NONMONO_V1 itself over the captured
! per-frequency inputs (KN_M_* store) in alternating frequency order:
!
!   descending pass (L=NCF..1): compressing depths read the freshly
!     updated J_NEXT_FRQ(:,L+1); expanding depths read the lagged
!     bluer-neighbour field seeded from J_NEXT_FRQ(:,L-1).
!   ascending pass (L=1..NCF):  expanding depths read the freshly
!     updated bluer field; compressing depths lag one pass.
!
! Because the refresh calls the production routine, the tridiagonal
! assembly, boundary conditions, DGTSV solve and negative/non-finite
! safeguards are identical to the in-sweep discretisation by
! construction.  MOM_KN_REFRESHING suppresses GET_MOMS_NON_REL (its
! EDDFACTOR-backed inputs exist only for the last sweep frequency)
! and the KN_M capture writes; the per-frequency Eddington factors
! and boundary data are restored from the captured store before each
! call.
!
! Only the lagged store J_NEXT_FRQ (and the lagged-H scratch) is
! refreshed; the per-frequency JNU/RSQHNU outputs of the re-calls go
! to scratch and are discarded, matching the FG refresh semantics.
! ==============================================================================
!
	SUBROUTINE MOM_KNOP_REFRESH_V1()
	USE SET_KIND_MODULE
	USE MOD_MOM_J_NONMONO_V1
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	REAL(KIND=LDP) OMP_GET_WTIME
	EXTERNAL OMP_GET_WTIME
!
	INTEGER L,I,IOS,IPASS,LB,LE,LDI
	INTEGER IT1,J,PASS_ARGI,PASS_ARGL
	INTEGER*8 PASS_CAP,PASS_NEG
	INTEGER*8 PASS_N5,PASS_N50,PASS_NQ,N5_SV,N50_SV,NQ_SV
	INTEGER*8 PASS_NCMP,PASS_NMATCH,PASS_NOFF
	INTEGER*8 BEST_N5,BEST_N50
	INTEGER*8 BEST_N5_D(2),BEST_N50_D(2),LAST_N50(2)
	INTEGER NSTALL_D(2),IDIR,PASS_NMAX,N_FLOOR
	INTEGER*8 NDMP,NACC
	INTEGER NHIST_AA,JJ,KK,MP_AA,MTRY,INFO_AA
	INTEGER*8 NAA_SKIP
	INTEGER IPIV_AA(MOM_AA_M)
	REAL(KIND=LDP) PASS_CHG,T0,BEST_CHG,MIN_RES,RESID0
	REAL(KIND=LDP) RES_PREV(2),RATIO_M,FCMP,BEST_CHG_D(2)
	REAL(KIND=LDP) PASS_NU,PASS_R,PASS_V,PASS_SIGT,PASS_GAM,PASS_Q
	REAL(KIND=LDP) T3,T4,STEP,CYC,TCOL,D1,D2,LAMC,T1X
	REAL(KIND=LDP) GRM_AA(MOM_AA_M,MOM_AA_M),BRHS_AA(MOM_AA_M)
	REAL(KIND=LDP) WGRM(MOM_AA_M,MOM_AA_M),WBR(MOM_AA_M)
	REAL(KIND=LDP) GAM_AA(MOM_AA_M),FNRM,DNRM,SC,XV,WV,DTR
	INTEGER IPIV_OUT(MOM_AA_OUT_M)
	REAL(KIND=LDP) GRM_OUT(MOM_AA_OUT_M,MOM_AA_OUT_M),
	1    BRHS_OUT(MOM_AA_OUT_M)
	REAL(KIND=LDP) WGRM_OUT(MOM_AA_OUT_M,MOM_AA_OUT_M),
	1    WBR_OUT(MOM_AA_OUT_M),GAM_OUT(MOM_AA_OUT_M)
	EXTERNAL DGESV
	LOGICAL LINIT,CONVERGED,BEST_VALID,KEEP_BEST
	LOGICAL BESTV_D(2),SKIP_D(2)
	CHARACTER(LEN=7) RF_ST
!
! Latest RSQHNU per frequency.  Double precision is required: the
! runaway guard limits XM_NM to ~1E3 times the previous valid field,
! which can still exceed the float32 range (~3.4E38) when the field
! itself is large; a REAL store overflowed to Inf and was seeded
! back as RSQHNU_PREV_NM, producing XMrhs=+/-Infinity rows with
! finite matrix coefficients (TypeIIn_KNOP_TRIDIAG, 22-Sep-2026).
! RSQJNU needs no scratch: J_NEXT_FRQ(:,L) already holds the latest
! relaxed R^2*J.
!
	REAL(KIND=LDP), ALLOCATABLE :: HRF(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: OUTJ(:),OUTH(:)
	REAL(KIND=LDP), ALLOCATABLE :: ETA_L(:),CHI_L(:),ESEC_L(:)
	REAL(KIND=LDP), ALLOCATABLE :: JNEXT_SAVE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: HNEXT_SAVE(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: JPREV1(:,:),JPREV2(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: JPREV3(:,:),JPREV4(:,:)
	REAL(KIND=LDP), ALLOCATABLE :: XHIST_AA(:,:,:),FHIST_AA(:,:,:)
!
! JBEST/HBEST hold the best-pass iterate so an unconverged refresh
! can retain the propagation it produced instead of reverting the
! lagged store to the pre-refresh field.
!
	REAL(KIND=LDP), ALLOCATABLE :: JBEST(:,:),HBEST(:,:)
!
! Active-set front-tracking (MOM_RF_ACTIVE): per-frequency mask,
! residual and maximum upwind offset.  SHOT_L(J) = any hot slot at
! index >= J, giving an O(1) "hot within the matched range" test.
!
	LOGICAL, ALLOCATABLE :: ACT_CUR(:),ACT_NEW(:),SHOT_L(:)
	REAL(KIND=LDP), ALLOCATABLE :: RES_L(:)
	INTEGER, ALLOCATABLE :: NMUP_L(:)
	INTEGER NACT
	LOGICAL DO_ACT
!
	LUER=ERROR_LU()
! This routine is also invoked on the FG path, where the MOM solver
! never ran and every guard is always false; only report the guard
! state when the MOM path is in use (KN_M_READY gates the same thing).
	IF(KN_M_READY)THEN
	WRITE(LUER,'(A,4L2,2I8,L2,2I8)')
	1    ' MOM_KNOP_REFRESH guard: RDY CMP CHISM JNXT | NDNM NDNMnow'//
	2    ' | VALID KONJsz ND_NM=',
	3    KN_M_READY,KN_M_ANYCMP,ALLOCATED(KN_M_CHISM),
	4    ALLOCATED(J_NEXT_FRQ),KN_M_NDNM,ND_NM,
	5    J_NEXT_VALID,
	6    MERGE(SIZE(K_ON_J),-1,ALLOCATED(K_ON_J)),ND_NM
	END IF
	IF(.NOT. KN_M_READY)RETURN
	IF(.NOT. KN_M_ANYCMP)RETURN
	IF(.NOT. ALLOCATED(KN_M_CHISM))RETURN
	IF(.NOT. ALLOCATED(J_NEXT_FRQ))RETURN
	IF(.NOT. J_NEXT_VALID)RETURN
	IF(KN_M_NCF .LE. 1)RETURN
	IF(.NOT. ALLOCATED(RSQJNU_NM))RETURN
	IF(.NOT. ALLOCATED(K_ON_J))RETURN
	IF(SIZE(K_ON_J) .NE. ND_NM)RETURN
	IF(KN_M_NDNM .NE. ND_NM)THEN
	  IF(.NOT. KN_M_RF_WARNED)THEN
	    WRITE(LUER,'(A,2I8)')' MOM_KNOP_REFRESH: ND_NM changed since'//
	1      ' capture -- refresh skipped',KN_M_NDNM,ND_NM
	    KN_M_RF_WARNED=.TRUE.
	  END IF
	  RETURN
	END IF
!
	ALLOCATE(HRF(ND_NM,KN_M_NCF),OUTJ(KN_M_NDSM),OUTH(KN_M_NDSM),
	1   ETA_L(KN_M_NDSM),CHI_L(KN_M_NDSM),ESEC_L(KN_M_NDSM),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  IF(.NOT. KN_M_RF_WARNED)THEN
	    WRITE(LUER,'(A)')' MOM_KNOP_REFRESH: allocation failed --'//
	1      ' refresh disabled for this run'
	    KN_M_RF_WARNED=.TRUE.
	  END IF
	  RETURN
	END IF
!
! Seed HRF with the in-sweep RSQHNU chain state: the only per-frequency
! H available is the final (reddest) one, but the first pass only reads
! HRF(:,L-1) for expanding depths, and the in-sweep ordering already
! produced a consistent blue-to-red H chain which is approximated by
! the lagged values written below as the pass proceeds.
!
	HRF=0.0_LDP
!
! One-shot store-integrity scan: a non-positive or non-finite captured
! ETA/CHI would NaN the LOG interp in the re-solve.  Report the first
! bad (I,L) so a poisoned slot can be traced to its capturing call.
!
	IF(.NOT. KN_M_RF_WARNED)THEN
	  IT1=0;  J=0
	  DO L=1,KN_M_NCF
	    DO I=1,KN_M_NDSM
	      IF(.NOT. (KN_M_ETASM(I,L) .GT. 0.0))THEN
	        IF(IT1 .EQ. 0)WRITE(LUER,'(A,2I8,ES13.4)')
	1          ' MOM_KNOP_REFRESH: first non-pos ETA at (I,L)=',
	2          I,L,KN_M_ETASM(I,L)
	        IT1=IT1+1
	      END IF
	      IF(.NOT. (KN_M_CHISM(I,L) .GT. 0.0))THEN
	        IF(J .EQ. 0)WRITE(LUER,'(A,2I8,ES13.4)')
	1          ' MOM_KNOP_REFRESH: first non-pos CHI at (I,L)=',
	2          I,L,KN_M_CHISM(I,L)
	        J=J+1
	      END IF
	    END DO
	  END DO
	  WRITE(LUER,'(A,I10,A,I10,A,I8,A,I8)')
	1    ' MOM_KNOP_REFRESH store scan: bad ETA=',IT1,' bad CHI=',J,
	2    '  NDSM=',KN_M_NDSM,'  NCF=',KN_M_NCF
	  KN_M_RF_WARNED=.TRUE.
	END IF
!
	MOM_KN_REFRESHING=.TRUE.
	T0=OMP_GET_WTIME()
!
	ALLOCATE(JNEXT_RF_READ(ND_NM,KN_M_NCF),
	1   JNEXT_SAVE(ND_NM,KN_M_NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(HNEXT_RF_READ(ND_NM,KN_M_NCF),
	1   HNEXT_SAVE(ND_NM,KN_M_NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(JPREV1(ND_NM,KN_M_NCF),
	1   JPREV2(ND_NM,KN_M_NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(JPREV3(ND_NM,KN_M_NCF),
	1   JPREV4(ND_NM,KN_M_NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(JBEST(ND_NM,KN_M_NCF),
	1   HBEST(ND_NM,KN_M_NCF),STAT=IOS)
	IF(IOS .EQ. 0 .AND. MOM_RF_AA)ALLOCATE(
	1   XHIST_AA(ND_NM,KN_M_NCF,MOM_AA_M+1),
	2   FHIST_AA(ND_NM,KN_M_NCF,MOM_AA_M+1),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  MOM_KN_REFRESHING=.FALSE.
	  DEALLOCATE(HRF,OUTJ,OUTH,ETA_L,CHI_L,ESEC_L)
	  IF(ALLOCATED(JNEXT_RF_READ))DEALLOCATE(JNEXT_RF_READ)
	  IF(ALLOCATED(JNEXT_SAVE))DEALLOCATE(JNEXT_SAVE)
	  IF(ALLOCATED(HNEXT_RF_READ))DEALLOCATE(HNEXT_RF_READ)
	  IF(ALLOCATED(HNEXT_SAVE))DEALLOCATE(HNEXT_SAVE)
	  IF(ALLOCATED(JPREV1))DEALLOCATE(JPREV1)
	  IF(ALLOCATED(JPREV2))DEALLOCATE(JPREV2)
	  IF(ALLOCATED(JPREV3))DEALLOCATE(JPREV3)
	  IF(ALLOCATED(JPREV4))DEALLOCATE(JPREV4)
	  IF(ALLOCATED(JBEST))DEALLOCATE(JBEST)
	  IF(ALLOCATED(HBEST))DEALLOCATE(HBEST)
	  IF(ALLOCATED(XHIST_AA))DEALLOCATE(XHIST_AA)
	  IF(ALLOCATED(FHIST_AA))DEALLOCATE(FHIST_AA)
	  RETURN
	END IF
!
! AA history restarts every call: each refresh iterates a different
! map (the source function was rebuilt from updated populations), so
! secants from the previous refresh describe a different operator.
!
	NHIST_AA=0
!
! Per-cell relax schedule: sized to the current grid, reset to the
! global factor on any resize.  Not checkpointed -- a requeue
! re-derives the schedule within the first refresh's passes.
!
	IF(.NOT. ALLOCATED(MOM_DAMP_NM) .OR.
	1   SIZE(MOM_DAMP_NM,1) .NE. ND_NM .OR.
	2   SIZE(MOM_DAMP_NM,2) .NE. KN_M_NCF)THEN
	  IF(ALLOCATED(MOM_DAMP_NM))DEALLOCATE(MOM_DAMP_NM)
	  ALLOCATE(MOM_DAMP_NM(ND_NM,KN_M_NCF),STAT=IOS)
	  IF(IOS .EQ. 0)MOM_DAMP_NM=JNEXT_RELAX
	END IF
!
! Slow recovery across calls: the schedule is decrease-only within a
! refresh (otherwise flagged cells bounce back and re-oscillate
! before the map can contract), so cells whose physics has moved
! drift back toward the global factor once per great iteration.
!
	IF(ALLOCATED(MOM_DAMP_NM))
	1   MOM_DAMP_NM=MIN(JNEXT_RELAX,1.15_LDP*MOM_DAMP_NM)
!
! Active-set mask: all frequencies on pass 1 -- the map's frozen
! coefficients were rebuilt since the last refresh, so a frequency
! cold under the old map cannot be trusted cold under the new one.
!
	DO_ACT=.FALSE.
	IF(MOM_RF_ACTIVE)THEN
	  ALLOCATE(ACT_CUR(KN_M_NCF),ACT_NEW(KN_M_NCF),
	1     SHOT_L(KN_M_NCF+1),RES_L(KN_M_NCF),NMUP_L(KN_M_NCF),
	2     STAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    ACT_CUR=.TRUE.
	    RES_L=HUGE(1.0_LDP)
	    NMUP_L=1
	    DO_ACT=.TRUE.
	  ELSE
	    WRITE(LUER,'(A)')' MOM_KNOP_REFRESH: active-set alloc'//
	1      ' failed -- full-frequency passes'
	    IF(ALLOCATED(ACT_CUR))DEALLOCATE(ACT_CUR)
	    IF(ALLOCATED(ACT_NEW))DEALLOCATE(ACT_NEW)
	    IF(ALLOCATED(SHOT_L))DEALLOCATE(SHOT_L)
	    IF(ALLOCATED(RES_L))DEALLOCATE(RES_L)
	    IF(ALLOCATED(NMUP_L))DEALLOCATE(NMUP_L)
	  END IF
	END IF
	IF(.NOT. ALLOCATED(MOM_DAMP_NM))THEN
	  MOM_KN_REFRESHING=.FALSE.
	  DEALLOCATE(HRF,OUTJ,OUTH,ETA_L,CHI_L,ESEC_L)
	  DEALLOCATE(JNEXT_RF_READ,JNEXT_SAVE)
	  DEALLOCATE(HNEXT_RF_READ,HNEXT_SAVE)
	  DEALLOCATE(JPREV1,JPREV2)
	  DEALLOCATE(JPREV3,JPREV4)
	  DEALLOCATE(JBEST,HBEST)
	  IF(ALLOCATED(XHIST_AA))DEALLOCATE(XHIST_AA)
	  IF(ALLOCATED(FHIST_AA))DEALLOCATE(FHIST_AA)
	  IF(ALLOCATED(ACT_CUR))DEALLOCATE(ACT_CUR)
	  IF(ALLOCATED(ACT_NEW))DEALLOCATE(ACT_NEW)
	  IF(ALLOCATED(SHOT_L))DEALLOCATE(SHOT_L)
	  IF(ALLOCATED(RES_L))DEALLOCATE(RES_L)
	  IF(ALLOCATED(NMUP_L))DEALLOCATE(NMUP_L)
	  RETURN
	END IF
!
! Divergence containment, mirroring KNOP_RAY_SOLVE_V1: a ray whose
! fixed-point iteration diverges never writes back, so the in-sweep
! field is retained.  Here the store is global, so keep a snapshot of
! the pre-refresh J_NEXT_FRQ plus a snapshot of the best-pass iterate
! (JBEST/HBEST).  An unconverged refresh restores the best-pass field
! when its residual improved on the first-pass residual -- retaining
! the red-to-blue propagation the refresh was written to provide --
! and falls back to the pre-refresh field when no pass improved on
! the seed (a map that grew from the start must not be kept: an
! earlier revision unconditionally restored the best pass and a
! 1E5-residual iterate containing blown-up values poisoned the next
! great iteration).
!
	BEST_CHG=HUGE(1.0_LDP)
	MIN_RES=HUGE(1.0_LDP)
	BEST_N5=HUGE(0_8);  BEST_N50=HUGE(0_8)
	RESID0=-1.0_LDP
	BEST_VALID=.FALSE.
	CONVERGED=.FALSE.
	RES_PREV=0.0_LDP
	NSTALL_D=0
	SKIP_D=.FALSE.
	BESTV_D=.FALSE.
	BEST_CHG_D=HUGE(1.0_LDP)
	BEST_N5_D=HUGE(0_8);  BEST_N50_D=HUGE(0_8)
	LAST_N50=0_8
	KN_M_DESC=.FALSE.
	N5_SV=MOM_GS_NBAD5;  N50_SV=MOM_GS_NBAD50
	NQ_SV=MOM_GS_NQUAR
!
! Snapshot the sweep-level maturity metrics before the pass loop
! resets them: the sweep drift (fresh solve vs stored field) is the
! maturity signal the populations see; the pass residuals the loop
! produces measure map self-consistency instead.
!
	MOM_SW_MAXCHG=MOM_GS_MAXCHG
	MOM_SW_NBAD5=MOM_GS_NBAD5
	MOM_SW_NBAD50=MOM_GS_NBAD50
	MOM_SW_NTOT=MOM_GS_NTOT
	MOM_SW_NEGJ=MOM_NEGJ_LRG
	MOM_SW_NQUAR=MOM_GS_NQUAR
	MOM_SW_CAPPED=MOM_INEXT_CAPPED
	MOM_SW_NHUGE=MOM_GS_NHUGE
	JNEXT_SAVE(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	HNEXT_SAVE(1:ND_NM,1:KN_M_NCF)=H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
!
! The retained composite starts as the seed: under MOM_RF_OWNW each
! direction's best pass merges only its owned rows, so a direction
! that never improves contributes the seed values rather than
! uninitialized storage.  [MOM_RF_OWNW]
!
	JBEST(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	HBEST(1:ND_NM,1:KN_M_NCF)=H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
!
! Re-resolve dominance from the sweep that just finished: the
! FREQ_INDX=1 flags were set from the previous sweep's census, which
! is one iteration stale.  The counters still hold this sweep's
! totals because refresh re-calls do not accumulate.
!
	IF(MOM_RF_DIR .EQ. 'DESC')THEN
	  MOM_DOM_CMP=.TRUE.;  MOM_DOM_EXP=.FALSE.
	ELSE IF(MOM_RF_DIR .EQ. 'ASC ')THEN
	  MOM_DOM_CMP=.FALSE.;  MOM_DOM_EXP=.TRUE.
	ELSE IF(MOM_RF_DIR .EQ. 'ALT ')THEN
	  MOM_DOM_CMP=.FALSE.;  MOM_DOM_EXP=.FALSE.
	ELSE IF(MOM_SW_NROW .GT. 0_8)THEN
	  MOM_SW_FCMP=DBLE(MOM_SW_NCROW)/DBLE(MOM_SW_NROW)
	  MOM_DOM_CMP=(10_8*MOM_SW_NCROW .GE. 9_8*MOM_SW_NROW)
	  MOM_DOM_EXP=(10_8*MOM_SW_NCROW .LE. MOM_SW_NROW)
	END IF
!
! Temporary DO_NONMONO_KNOP diagnostic table.  The refresh can stall
! with a falling bulk residual while a small compression region grows.
! One row per pass ties the same-direction contraction ratio to that
! region and measures whether DO_NONMONO_FREQ_MATCH actually selects
! non-adjacent upwind bins.  The stock path never calls this routine.
!
	WRITE(LUER,'(/,A)')' ------------------------------ MOM_KNOP_REFRESH -------------------------------'
	WRITE(LUER,'(A)')'  lagged-field re-solve, direction-ordered'//
	1    ' frequency passes'
	WRITE(LUER,'(A)')'   per pass: drift = max rel. change'//
	1    ' (+ same-direction ratio); n5/n50 = cells'//
	2    ' >5%/>50%; argmax = nu,R,V,sigT,GAM,Q at'//
	3    ' the worst cell; match/off = upwind reads'//
	4    ' frequency-matched / dropped off-grid'
!
! Report the resolved pass order.  On a dominant grid the minority
! direction is a pure Jacobi step that re-injects a full pass of
! error, so only the dominant (fresh-read) direction runs -- the
! physically-guided ordering of Knop Paper II Sect. 3.2.
!
	FCMP=MOM_SW_FCMP
	IF(MOM_DOM_CMP)THEN
	  WRITE(LUER,'(A,F6.3,A)')'   compressing rows=',FCMP,
	1    ' -- pass order DESC (dominant-direction only)'
	ELSE IF(MOM_DOM_EXP)THEN
	  WRITE(LUER,'(A,F6.3,A)')'   compressing rows=',FCMP,
	1    ' -- pass order ASC (dominant-direction only)'
	ELSE
	  WRITE(LUER,'(A,F6.3,A)')'   compressing rows=',FCMP,
	1    ' -- pass order alternating'
	END IF
	DO IPASS=1,KN_M_MAXPASS
	  PASS_CHG=0.0_LDP
	  PASS_CAP=0
	  PASS_NEG=0
	  PASS_N5=0;  PASS_N50=0;  PASS_NQ=0
	  PASS_NCMP=0;  PASS_NMATCH=0;  PASS_NOFF=0
	  PASS_ARGI=0;  PASS_ARGL=0;  PASS_NMAX=1
	  PASS_NU=0.0_LDP;  PASS_R=0.0_LDP;  PASS_V=0.0_LDP
	  PASS_SIGT=0.0_LDP;  PASS_GAM=0.0_LDP;  PASS_Q=0.0_LDP
	  NDMP=0_8;  NACC=0_8
!
! Per-pass reset of the physics accumulators -- the in-routine
! FREQ_INDX=1 reset is suppressed during refresh, so without this a
! descending pass ending on L=1 would leave the sweep-report stats at
! sentinel values (the "blank GS sweep" blocks in OUTGEN).
!
	  MOM_NONMM_CNT=0
	  MOM_EPS_FWD=0
	  MOM_EPS_RAT=0.0_LDP
	  MOM_BND_HBC_FL=0
	  MOM_N_DGTSV=0
	  MOM_N_NEGJ_T=0
	  MOM_N_NONF=0
	  MOM_N_RUNAWAY=0
	  MOM_TB_NONPOS=0
	  MOM_NEGJ_DUMP=0
	  MOM_EPS_CAP=0
	  MOM_RF_NSUPP=0
	  MOM_MIN_DR=1.0E+30_LDP
	  MOM_MIN_DR_I=0;  MOM_MIN_DR_L=0
	  MOM_MIN_DRA=1.0E+30_LDP
	  MOM_MIN_DRA_I=0;  MOM_MIN_DRA_L=0
	  MOM_EQ_RES=0.0_LDP
	  MOM_EQ_RES_I=0;  MOM_EQ_RES_L=0
	  MOM_MAX_COH=0.0_LDP
	  MOM_MAX_GAM=0.0_LDP
	  MOM_MAX_SIGT=0.0_LDP
	  MOM_LAGB_CNT=0
	  MOM_LAGB_MIN=1.0_LDP
!
! Jacobi update for the compressing upwind read: freeze the
! J_{k+1} field at pass start so descending passes do not consume
! this-pass writes (the GS cascade diverged to ~1E10 residuals on
! job 11887664).  The routine reads JNEXT_RF_READ while
! MOM_KN_REFRESHING is set; the write side (J_NEXT_FRQ) still gets
! the JNEXT_RELAX-damped update inside the solver.  When
! MOM_RF_GS_DESC is set, descending passes instead read the live
! J_NEXT_FRQ(:,L+1) -- KN_M_DESC tells the solver which pass
! direction is running.
!
	  JNEXT_RF_READ(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	  HNEXT_RF_READ(1:ND_NM,1:KN_M_NCF)=H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
!
! Pass order: on a dominant grid every pass runs in the direction
! whose upwind read is fresh (descending freshens the compressing
! J_NEXT read; ascending freshens the expanding J_PREV read).  The
! alternating schedule is kept for genuinely mixed fields, where each
! direction owns a different row subset.  A direction whose
! same-parity drift ratio has stalled (SKIP_D) is dropped from the
! schedule so the still-contracting direction keeps its passes --
! the row subsets are disjoint under MOM_RF_OWNW, so dropping a
! stalled direction cannot disturb the other.  The two-pass
! composite-map machinery below (cycle check, AA, same-direction
! stall) stays valid: in the single-direction modes both passes of a
! pair run the same map, so the composite is simply M^2.
!
	  IF(MOM_DOM_CMP)THEN
	    LB=KN_M_NCF;  LE=1;  LDI=-1     !all descending: fresh J_NEXT
	  ELSE IF(MOM_DOM_EXP)THEN
	    LB=1;  LE=KN_M_NCF;  LDI=1      !all ascending: fresh J_PREV
	  ELSE IF(SKIP_D(1) .AND. .NOT. SKIP_D(2))THEN
	    LB=KN_M_NCF;  LE=1;  LDI=-1
	  ELSE IF(SKIP_D(2) .AND. .NOT. SKIP_D(1))THEN
	    LB=1;  LE=KN_M_NCF;  LDI=1
	  ELSE IF(MOD(IPASS,2) .EQ. 1)THEN
	    LB=1;  LE=KN_M_NCF;  LDI=1      !alternating: freshens J_PREV
	  ELSE
	    LB=KN_M_NCF;  LE=1;  LDI=-1     !alternating: freshens J_NEXT
	  END IF
	  KN_M_DESC=(LDI .LT. 0)
!
	  DO L=LB,LE,LDI
!
! Inactive frequencies keep their last-pass field: the map output
! there was already below tolerance with quiet inputs, and the
! frozen coefficients cannot have changed mid-refresh.
!
	    IF(DO_ACT)THEN
	      IF(.NOT. ACT_CUR(L))CYCLE
	    END IF
	    MOM_GS_MAXCHG=0.0_LDP
	    MOM_INEXT_CAPPED=0
	    MOM_NEGJ_LRG=0
	    MOM_GS_NBAD5=0;  MOM_GS_NBAD50=0
	    MOM_GS_NHUGE=0
	    MOM_GS_NQUAR=0
	    MOM_GS_ARGI=0;  MOM_GS_ARGL=0
!
! Current-frequency Eddington/boundary state from the captured store.
!
	    K_ON_J(1:ND_NM)=KN_MKONJ(1:ND_NM,L)
	    NMID_ON_J(1:ND_NM)=KN_MNMJ(1:ND_NM,L)
	    NMID_ON_HMID(1:ND_NM)=KN_MNMH(1:ND_NM,L)
!
! Same sanitize as the GET_MOMS_NON_REL path: a store captured before
! the guards existed (or from a failed sweep) can hold NaN factors.
!
	    DO I=1,ND_NM
	      IF(.NOT. (ABS(K_ON_J(I)) .LE. HUGE(K_ON_J(I))))
	1        K_ON_J(I)=1.0_LDP/3.0_LDP
	      K_ON_J(I)=MIN(1.0_LDP,MAX(1.0E-4_LDP,K_ON_J(I)))
	    END DO
	    DO I=1,ND_NM-1
	      IF(.NOT. (ABS(NMID_ON_J(I)) .LE. HUGE(NMID_ON_J(I))))
	1        NMID_ON_J(I)=0.0_LDP
	      IF(.NOT. (ABS(NMID_ON_HMID(I)) .LE. HUGE(NMID_ON_HMID(I))))
	1        NMID_ON_HMID(I)=0.0_LDP
	      NMID_ON_HMID(I)=MIN(1.0_LDP,MAX(0.0_LDP,NMID_ON_HMID(I)))
	    END DO
	    HBC=KN_MHBC(L);  NBC=KN_MNBC(L);  IN_HBC=KN_MINHBC(L)
	    JPLUS_IB=KN_MJPIB(L);  JMIN_IB=KN_MJMIB(L)
	    HNU_AT_IB=KN_MHIB(L)
	    MOM_NU_NEXT=KN_MNUNXT(L)
!
	    LINIT=(L .EQ. 1)
	    IF(.NOT. LINIT)THEN
!
! Lagged bluer-neighbour state: the routine's NEW_FREQ block copies
! *_SAVE -> *_PREV and RSQ{J,H}NU_NM -> RSQ{J,H}NU_PREV_NM.  The
! Eddington PREV values are captured inputs; the J/H PREV fields are
! the iterated ones (latest available: this pass if L-1 was already
! visited, else the previous pass).
!
	      K_ON_J_SAVE(1:ND_NM)=KN_MKONJ(1:ND_NM,L-1)
	      NMID_ON_J_SAVE(1:ND_NM)=KN_MNMJ(1:ND_NM,L-1)
	      NMID_ON_HMID_SAVE(1:ND_NM)=KN_MNMH(1:ND_NM,L-1)
	      DO I=1,ND_NM
	        IF(.NOT. (ABS(K_ON_J_SAVE(I)) .LE. HUGE(K_ON_J_SAVE(I))))
	1          K_ON_J_SAVE(I)=1.0_LDP/3.0_LDP
	        K_ON_J_SAVE(I)=MIN(1.0_LDP,MAX(1.0E-4_LDP,K_ON_J_SAVE(I)))
	      END DO
	      DO I=1,ND_NM-1
	        IF(.NOT. (ABS(NMID_ON_J_SAVE(I)) .LE.
	1            HUGE(NMID_ON_J_SAVE(I))))NMID_ON_J_SAVE(I)=0.0_LDP
	        IF(.NOT. (ABS(NMID_ON_HMID_SAVE(I)) .LE.
	1            HUGE(NMID_ON_HMID_SAVE(I))))
	2          NMID_ON_HMID_SAVE(I)=0.0_LDP
	        NMID_ON_HMID_SAVE(I)=MIN(1.0_LDP,
	1            MAX(0.0_LDP,NMID_ON_HMID_SAVE(I)))
	      END DO
	      HBC_SAVE=KN_MHBC(L-1);  NBC_SAVE=KN_MNBC(L-1)
	      IN_HBC_SAVE=KN_MINHBC(L-1)
	      RSQJNU_NM(1:ND_NM)=J_NEXT_FRQ(1:ND_NM,L-1)
	      RSQHNU_NM(1:ND_NM)=HRF(1:ND_NM,L-1)
	    END IF
!
	    DO I=1,KN_M_NDSM
	      ETA_L(I)=KN_M_ETASM(I,L)
	      CHI_L(I)=KN_M_CHISM(I,L)
	      ESEC_L(I)=KN_M_ESECSM(I,L)
	    END DO
!
	    CALL MOM_J_CMF_NONMONO_V1(ETA_L,CHI_L,ESEC_L,
	1        KN_M_VSM,KN_M_SIGSM,KN_M_SIGTSM,KN_M_RSM,
	1        OUTJ,OUTH,KN_M_VDOPV,KN_M_VDOPF,
	1        KN_MFREQ(L),KN_MDLN(L),
	1        KN_M_INBM,KN_MDBB(L),KN_MIC(L),KN_M_IBSTAB,
	1        KN_M_NTYPE,KN_M_HCHK,KN_M_METHOD,KN_M_COH,KN_M_OUTBC(L),
	1        LINIT,.NOT. LINIT,KN_M_NC,KN_M_NP,KN_M_NDSM,
	1        L,KN_M_NCF)
!
	    HRF(1:ND_NM,L)=RSQHNU_NM(1:ND_NM)
	    PASS_NCMP=PASS_NCMP+INT(COUNT(GAM_SIGNED_NM(1:ND_NM)
	1    .LT. 0.0_LDP),8)
	    PASS_NMATCH=PASS_NMATCH+INT(COUNT(MOM_NM_J(1:ND_NM)
	1    .GT. 1),8)
	    PASS_NOFF=PASS_NOFF+INT(COUNT(MOM_NM_J(1:ND_NM)
	1    .EQ. 0),8)
	    PASS_NMAX=MAX(PASS_NMAX,MAXVAL(MOM_NM_J(1:ND_NM)))
	    IF(DO_ACT)THEN
	      RES_L(L)=MOM_GS_MAXCHG
	      NMUP_L(L)=MAX(MAXVAL(MOM_NM_J(1:ND_NM)),
	1        MAXVAL(MOM_NM_H(1:ND_NM)))
	    END IF
	    IF(MOM_GS_MAXCHG .GT. PASS_CHG)THEN
	      PASS_ARGI=MOM_GS_ARGI;  PASS_ARGL=MOM_GS_ARGL
	      IF(PASS_ARGL .GE. 1 .AND. PASS_ARGL .LE. KN_M_NCF)
	1      PASS_NU=KN_MFREQ(PASS_ARGL)
	      IF(PASS_ARGI .GE. 1 .AND. PASS_ARGI .LE. ND_NM)THEN
	        PASS_R=R_NM(PASS_ARGI)
	        PASS_V=V_NM(PASS_ARGI)
	        PASS_SIGT=SIGMA_TRUE_NM(PASS_ARGI)
	        PASS_GAM=GAM_SIGNED_NM(PASS_ARGI)
	        PASS_Q=Q_NM(PASS_ARGI)
	      END IF
	    END IF
	    PASS_CHG=MAX(PASS_CHG,MOM_GS_MAXCHG)
	    PASS_CAP=PASS_CAP+MOM_INEXT_CAPPED
	    PASS_NEG=PASS_NEG+MOM_NEGJ_LRG
	    PASS_N5=PASS_N5+MOM_GS_NBAD5
	    PASS_N50=PASS_N50+MOM_GS_NBAD50
	    PASS_NQ=PASS_NQ+MOM_GS_NQUAR
	  END DO
!
! Active-set update: keep L active next pass when its own residual
! was hot, its expanding input L-1 was hot, or any slot in
! (L,L+NMUP] -- the range every compressing cell's matched read
! can reach -- was hot.  SHOT_L(J)=any hot index >=J turns the
! range test into a two-point lookup.  L=1 stays active so the
! LINIT branch of the solver always runs.  Cold cells retain their
! last-pass residual, so a region that re-ignites upstream pulls
! its readers back in one pass.
!
	  IF(DO_ACT)THEN
	    SHOT_L(KN_M_NCF+1)=.FALSE.
	    DO L=KN_M_NCF,1,-1
	      SHOT_L(L)=(RES_L(L) .GT. KN_M_TOL) .OR. SHOT_L(L+1)
	    END DO
	    DO L=1,KN_M_NCF
	      ACT_NEW(L)=RES_L(L) .GT. KN_M_TOL
	      IF(L .GT. 1)THEN
	        IF(RES_L(L-1) .GT. KN_M_TOL)ACT_NEW(L)=.TRUE.
	      END IF
	      IF(L .LT. KN_M_NCF .AND. SHOT_L(L+1))THEN
	        JJ=MIN(KN_M_NCF+1,L+NMUP_L(L)+1)
	        IF(.NOT. SHOT_L(JJ))ACT_NEW(L)=.TRUE.
	      END IF
	    END DO
	    ACT_NEW(1)=.TRUE.
	    ACT_CUR(1:KN_M_NCF)=ACT_NEW(1:KN_M_NCF)
	    NACT=COUNT(ACT_CUR(1:KN_M_NCF))
	    IF(NACT .LT. KN_M_NCF)WRITE(LUER,'(A,I8,A,I8)')
	1      '    active freqs=',NACT,' of',KN_M_NCF
	  END IF
!
	  MOM_GS_MAXCHG=PASS_CHG
	  MOM_INEXT_CAPPED=PASS_CAP
	  MOM_NEGJ_LRG=PASS_NEG
!
! edited by Devin, 28-Sep-2026: best-pass selection runs here, on
! the raw map output, before the cycle check and Anderson update
! modify J_NEXT_FRQ -- JBEST then holds the field whose PASS_N5/
! PASS_N50 were actually measured.  Selection ranks by the bulk
! metrics (cells drifting >50%, then >5%) rather than the max-norm
! residual, which is dominated by a few marginal compressing cells
! and is non-monotone across passes.  MIN_RES keeps the smallest
! max-norm seen purely as the divergence-abort reference.
!
	  IF(IPASS .EQ. 1)RESID0=PASS_CHG
	  IF(PASS_CHG .LT. MIN_RES)MIN_RES=PASS_CHG
!
! edited by Devin, 30-Sep-2026: rank each pass direction separately.
! Under MOM_RF_OWNW the pass metrics count only that direction's
! owned rows, so a whole-field snapshot ranked by raw n50 is biased
! toward whichever direction owns fewer rows -- the other direction's
! rows could regress while its side of the snapshot won.  JBEST/HBEST
! is now a composite: each direction's least-bad owned rows, merged
! through the MOM_OWN_CMP map recorded by the drift loop.  Without
! OWNW every row is owned by every pass and the merge degenerates to
! the old whole-array snapshot.
!
	  IDIR=1
	  IF(LDI .LT. 0)IDIR=2
	  LAST_N50(IDIR)=PASS_N50
	  IF(.NOT. BESTV_D(IDIR) .OR. PASS_N50 .LT. BEST_N50_D(IDIR) .OR.
	1    (PASS_N50 .EQ. BEST_N50_D(IDIR) .AND.
	2     PASS_N5 .LT. BEST_N5_D(IDIR)))THEN
	    IF(MOM_RF_OWNW)THEN
!
!	      J rows merge through MOM_OWN_CMP (cell-centre Q_NM sign),
!	      H rows through MOM_OWN_CMPH (midpoint GAMH sign) -- the two
!	      disagree where sigma crosses zero.
!
	      DO L=1,KN_M_NCF
	        DO I=1,ND_NM
	          IF(MOM_OWN_CMP(I,L) .EQV. (IDIR .EQ. 2))THEN
	            JBEST(I,L)=J_NEXT_FRQ(I,L)
	          END IF
	          IF(MOM_OWN_CMPH(I,L) .EQV. (IDIR .EQ. 2))THEN
	            HBEST(I,L)=H_NEXT_FRQ(I,L)
	          END IF
	        END DO
	      END DO
	    ELSE
	      JBEST(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	      HBEST(1:ND_NM,1:KN_M_NCF)=H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	    END IF
	    BEST_CHG_D(IDIR)=PASS_CHG
	    BEST_N5_D(IDIR)=PASS_N5;  BEST_N50_D(IDIR)=PASS_N50
	    BESTV_D(IDIR)=.TRUE.
	    WRITE(LUER,'(A,I4,3A,I12,A,I12,A,ES12.4)')
	1      '    best pass updated: pass',IPASS,'  dir=',
	2      MERGE('D','A',IDIR.EQ.2),'  n50=',BEST_N50_D(IDIR),
	3      '  n5=',BEST_N5_D(IDIR),'  residual=',BEST_CHG_D(IDIR)
	  END IF
	  BEST_VALID=.FALSE.
	  IF(BESTV_D(1) .OR. BESTV_D(2))BEST_VALID=.TRUE.
	  BEST_N5=0_8;  BEST_N50=0_8;  BEST_CHG=0.0_LDP
	  IF(BESTV_D(1))THEN
	    BEST_N5=BEST_N5+BEST_N5_D(1)
	    BEST_N50=BEST_N50+BEST_N50_D(1)
	    BEST_CHG=MAX(BEST_CHG,BEST_CHG_D(1))
	  END IF
	  IF(BESTV_D(2))THEN
	    BEST_N5=BEST_N5+BEST_N5_D(2)
	    BEST_N50=BEST_N50+BEST_N50_D(2)
	    BEST_CHG=MAX(BEST_CHG,BEST_CHG_D(2))
	  END IF
!
! Same-parity cycle check: a cell that returns within 30 percent of
! its two-passes-ago value while still stepping more than 50 percent
! (the NBAD50 significance floor) is cycling under the alternating-
! direction map.  Pass parity matches the upwind direction, so
! comparing k to k-2 isolates the same composite map.  Flagged cells
! get their relax factor halved.  The schedule is decrease-only
! within a refresh -- letting settled cells recover here bounces
! them straight back into oscillation; recovery happens once per
! call instead (above).  Snapshots keep the last four passes so the
! same-parity subsequence (k, k-2, k-4) is available for both the
! cycle check and the drift extrapolation below.
!
	  IF(IPASS .EQ. 1)THEN
	    JPREV1(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	    JPREV2(1:ND_NM,1:KN_M_NCF)=JPREV1(1:ND_NM,1:KN_M_NCF)
	    JPREV3(1:ND_NM,1:KN_M_NCF)=JPREV1(1:ND_NM,1:KN_M_NCF)
	    JPREV4(1:ND_NM,1:KN_M_NCF)=JPREV1(1:ND_NM,1:KN_M_NCF)
	  ELSE
	    IF(IPASS .GT. 2)THEN
	      NDMP=0_8;  NACC=0_8
	      DO L=1,KN_M_NCF
	        TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1          1.0E-28_LDP)
	        DO I=1,ND_NM
	          STEP=ABS(J_NEXT_FRQ(I,L)-JPREV1(I,L))
	          IF(STEP .GT. 0.5_LDP*MAX(ABS(JPREV1(I,L)),
	1            1.0E-4_LDP*TCOL))THEN
	            CYC=ABS(J_NEXT_FRQ(I,L)-JPREV2(I,L))
	            IF(CYC .LT. 0.3_LDP*STEP)THEN
	              MOM_DAMP_NM(I,L)=MAX(MOM_DAMP_MIN,
	1                0.5_LDP*MOM_DAMP_NM(I,L))
	              NDMP=NDMP+1
	            ELSE IF(IPASS .GE. 5)THEN
!
! Same-direction drift of the composite map: D1=x_k-x_{k-2} and
! D2=x_{k-2}-x_{k-4} are successive same-parity steps, so LAMC
! estimates the composite-map multiplier.  0<LAMC<1 means marginal
! contraction, which under-relaxation only slows; Aitken instead
! extrapolates the geometric tail x* ~ x_k + D1*LAMC/(1-LAMC).
! LAMC is capped at 0.9 (correction <= 9x step).  edited by Devin,
! 01-Oct-2026: a cell whose same-parity step is not contracting
! (LAMC<=0 sign flip without a period-2 return, or LAMC>=1 growth)
! can neither converge nor be extrapolated; it is one of the
! marginal-mode cells the sweep retains and cycles between calls.
! Halve its relax factor toward the floor instead of leaving it
! unchecked when AA owns the Aitken branch.  D2=0 (a previously
! stationary cell freshly kicked) is exempted -- the kick reflects
! a changed map, not a divergent mode.
!
	              D1=J_NEXT_FRQ(I,L)-JPREV2(I,L)
	              D2=JPREV2(I,L)-JPREV4(I,L)
	              IF(D2 .NE. 0.0_LDP)THEN
	                LAMC=D1/D2
	                IF(LAMC .LE. 0.0_LDP .OR. LAMC .GE. 1.0_LDP)THEN
	                  MOM_DAMP_NM(I,L)=MAX(MOM_DAMP_MIN,
	1                    0.5_LDP*MOM_DAMP_NM(I,L))
	                  NDMP=NDMP+1
	                ELSE IF(.NOT. MOM_RF_AA .AND.
	1                  MOM_DAMP_NM(I,L) .GE. 0.99_LDP*JNEXT_RELAX
	2                  .AND. ABS(D1) .GT. 0.5_LDP*
	3                  MAX(ABS(JPREV2(I,L)),1.0E-4_LDP*TCOL))THEN
	                  LAMC=MIN(LAMC,0.9_LDP)
	                  J_NEXT_FRQ(I,L)=J_NEXT_FRQ(I,L)+
	1                    D1*LAMC/(1.0_LDP-LAMC)
	                  NACC=NACC+1
	                END IF
	              END IF
	            END IF
	          END IF
	        END DO
	      END DO
	      IF(NDMP .GT. 0_8)WRITE(LUER,'(A,I12)')
	1        '    newly damped this pass=',NDMP
	    END IF
!
! Anderson acceleration on the composite map at even pass boundaries.
! X_i = JPREV2 (the field two passes back, i.e. the previous
! composite iterate), G_i = J_NEXT_FRQ (raw map output).  History
! slot 1 is the newest pair (X_i, f_i=G_i-X_i); older pairs sit at
! slots j+1.  Secants use differences against slot 1, so the Gram
! system is (f_i - f_{i-j}) . (f_i - f_{i-k}) and the correction is
! Sum_j gamma_j (G_i - G_{i-j}) = Sum_j gamma_j (J_NEXT - X_{i-j} -
! f_{i-j}).  The update is capped at 10||f_i|| and falls back to the
! raw map output on cells the extrapolation drives non-positive.
!
	    IF(MOM_RF_AA .AND. MOD(IPASS,2) .EQ. 0 .AND.
	1       IPASS .GE. 2 .AND. ALLOCATED(XHIST_AA))THEN
!
! The cycle check above halves MOM_DAMP_NM on flagged cells, which
! changes the map itself: secants stored before that change belong
! to a different operator and corrupt the AA mix.  Restart the
! history only when a non-trivial fraction of cells was re-damped
! this pass -- a handful of damped cells among ~12M does not move
! the global map enough to matter (MOM_AA_DMP_FRAC).
!
	      IF(DBLE(NDMP) .GT. MOM_AA_DMP_FRAC*
	1        DBLE(ND_NM)*DBLE(KN_M_NCF))THEN
	        NHIST_AA=0
	        WRITE(LUER,'(A,I12)')
	1          '    AA history restarted (damped=',NDMP
	      END IF
	      DO JJ=MOM_AA_M+1,2,-1
	        IF(JJ-1 .LE. NHIST_AA)THEN
	          XHIST_AA(1:ND_NM,1:KN_M_NCF,JJ)=
	1            XHIST_AA(1:ND_NM,1:KN_M_NCF,JJ-1)
	          FHIST_AA(1:ND_NM,1:KN_M_NCF,JJ)=
	1            FHIST_AA(1:ND_NM,1:KN_M_NCF,JJ-1)
	        END IF
	      END DO
	      XHIST_AA(1:ND_NM,1:KN_M_NCF,1)=
	1        JPREV2(1:ND_NM,1:KN_M_NCF)
	      FHIST_AA(1:ND_NM,1:KN_M_NCF,1)=
	1        J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)-JPREV2(1:ND_NM,1:KN_M_NCF)
	      NHIST_AA=MIN(NHIST_AA+1,MOM_AA_M+1)
	      MP_AA=NHIST_AA-1
	      IF(MP_AA .GE. 1)THEN
!
! The Gram system is built on residuals weighted by the same
! normalisation as the convergence metric (1/MAX(|X|,1E-4*TCOL)):
! the unweighted L2 norm across ~100 decades of J lets the largest
! cells alone choose gamma, leaving the small-J cells -- which the
! max-norm metric judges -- effectively unrepresented.
!
	        GRM_AA=0.0_LDP;  BRHS_AA=0.0_LDP
	        NAA_SKIP=0_8
	        DO L=1,KN_M_NCF
	          TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1            1.0E-28_LDP)
	          DO I=1,ND_NM
	            WV=1.0_LDP/MAX(ABS(XHIST_AA(I,L,1)),
	1              1.0E-4_LDP*TCOL)
	            SC=FHIST_AA(I,L,1)*WV
!
!	    Non-affine cell exclusion (edited by Devin, 30-Sep-2026):
!	    a cell whose weighted residual or secant exceeds ~100% per
!	    composite step is outside the affine-map regime -- the
!	    quarantine, limiter and cycle damper own those cells.  Their
!	    secants otherwise dominate the Gram system and the mixed
!	    gamma extrapolates the runaway rather than the marginal bulk:
!	    with all cells included, ||corr||>2||f|| vetoed every AA step
!	    on job 12046695 while the contracting bulk went unaccelerated.
!
	            T1X=ABS(SC)
	            DO JJ=1,MP_AA
	              T1X=MAX(T1X,ABS(FHIST_AA(I,L,1)-
	1              FHIST_AA(I,L,JJ+1))*WV)
	            END DO
	            IF(T1X .GT. MOM_AA_GCAP)THEN
	              NAA_SKIP=NAA_SKIP+1
	              CYCLE
	            END IF
	            DO JJ=1,MP_AA
	              XV=(FHIST_AA(I,L,1)-FHIST_AA(I,L,JJ+1))*WV
	              BRHS_AA(JJ)=BRHS_AA(JJ)+XV*SC
	              DO KK=1,JJ
	                GRM_AA(KK,JJ)=GRM_AA(KK,JJ)+
	1                 XV*(FHIST_AA(I,L,1)-FHIST_AA(I,L,KK+1))*WV
	              END DO
	            END DO
	          END DO
	        END DO
	        DO JJ=1,MP_AA
	          DO KK=JJ+1,MP_AA
	            GRM_AA(JJ,KK)=GRM_AA(KK,JJ)
	          END DO
	        END DO
!
! Tikhonov regularisation: nearly dependent secants make the Gram
! matrix near-singular and DGESV returns gamma ~ O(1/rcond), which
! the 10*||f|| cap then rescales globally while still allowing huge
! per-cell steps.  A small relative ridge keeps gamma bounded by
! the data scale.
!
	        DTR=0.0_LDP
	        DO JJ=1,MP_AA
	          DTR=MAX(DTR,GRM_AA(JJ,JJ))
	        END DO
	        DTR=1.0E-8_LDP*MAX(DTR,1.0E-30_LDP)
	        DO JJ=1,MP_AA
	          GRM_AA(JJ,JJ)=GRM_AA(JJ,JJ)+DTR
	        END DO
	        WGRM(1:MP_AA,1:MP_AA)=GRM_AA(1:MP_AA,1:MP_AA)
	        WBR(1:MP_AA)=BRHS_AA(1:MP_AA)
!
! If the cap excluded almost every cell the Gram is built on
! residual noise -- an unrepresentative fit is worse than none.
!
	        IF(ND_NM*INT(KN_M_NCF,8)-NAA_SKIP .LT.
	1           MAX(100_8,ND_NM*INT(KN_M_NCF,8)/1000))THEN
	          INFO_AA=-999
	          WRITE(LUER,'(A,I12)')
	1            '    AA skipped -- affine subset empty, excl=',
	2            NAA_SKIP
	        ELSE
	          CALL DGESV(MP_AA,1,WGRM,MOM_AA_M,IPIV_AA,WBR,MOM_AA_M,
	1            INFO_AA)
	        END IF
	        IF(INFO_AA .EQ. 0)THEN
	          GAM_AA(1:MP_AA)=WBR(1:MP_AA)
	          FNRM=0.0_LDP;  DNRM=0.0_LDP
	          DO L=1,KN_M_NCF
	            TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1              1.0E-28_LDP)
	            DO I=1,ND_NM
	              WV=1.0_LDP/MAX(ABS(XHIST_AA(I,L,1)),
	1                1.0E-4_LDP*TCOL)
	              IF(ABS(FHIST_AA(I,L,1))*WV .GT. MOM_AA_GCAP)
	1                CYCLE
	              FNRM=FNRM+(FHIST_AA(I,L,1)*WV)**2
	              SC=0.0_LDP
	              DO JJ=1,MP_AA
	                SC=SC+GAM_AA(JJ)*(J_NEXT_FRQ(I,L)-
	1                 XHIST_AA(I,L,JJ+1)-FHIST_AA(I,L,JJ+1))*WV
	              END DO
	              DNRM=DNRM+SC*SC
	            END DO
	          END DO
	          FNRM=SQRT(FNRM);  DNRM=SQRT(DNRM)
!
! edited by Devin, 29-Sep-2026: skip the AA update when the mixed
! correction outgrows the weighted residual.  For an affine map the
! optimal AA step cancels the residual, ||corr|| ~ ||f||; anything
! larger overshoots the data and is discarded -- the observed
! corr/f sequence (0.11, 0.51, 1.79, 4.3, ...) shows the deep
! windows fitting the non-affine kinks (limiter/quarantine/damping
! switching cells between regimes), where extrapolation amplifies
! rather than contracts.  The raw map output is kept on a skip.
!
	          IF(DNRM .GT. FNRM)THEN
	            WRITE(LUER,'(A,I2,A,ES12.4,A,ES12.4,A,I12)')
	1            '    AA m=',MP_AA,'  ||f||=',FNRM,
	2            '  ||corr||=',DNRM,'  skipped (>||f||) excl=',
	3            NAA_SKIP
	          ELSE
	            SC=1.0_LDP
	          DO L=1,KN_M_NCF
	            TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1              1.0E-28_LDP)
	            DO I=1,ND_NM
!
! Cells excluded from the Gram fit (weighted residual > MOM_AA_GCAP)
! keep the raw map output: an extrapolation fitted without them has
! no warrant on their scale.
!
	              WV=1.0_LDP/MAX(ABS(XHIST_AA(I,L,1)),
	1                1.0E-4_LDP*TCOL)
	              IF(ABS(FHIST_AA(I,L,1))*WV .GT. MOM_AA_GCAP)
	1                CYCLE
	              XV=0.0_LDP
	              DO JJ=1,MP_AA
	                XV=XV+GAM_AA(JJ)*(J_NEXT_FRQ(I,L)-
	1                 XHIST_AA(I,L,JJ+1)-FHIST_AA(I,L,JJ+1))
	              END DO
	              T4=J_NEXT_FRQ(I,L)
	              XV=T4-SC*XV
!
! Per-cell cap on the extrapolation: a cell may not move further
! than ten times its own (floored) residual, so one bad gamma
! cannot overshoot a small-J cell far above its scale.
!
	              T3=1.0E+1_LDP*MAX(ABS(FHIST_AA(I,L,1)),
	1                1.0E-4_LDP*TCOL)
	              STEP=T4-XV
	              IF(ABS(STEP) .GT. T3)
	1                XV=T4-SIGN(T3,STEP)
	              IF(XV .LE. 0.0_LDP .AND. T4 .GT. 0.0_LDP)
	1                XV=T4
	              J_NEXT_FRQ(I,L)=XV
!
! edited by Devin, 28-Sep-2026: drag H along by the same ratio as
! J so the retained moment pair keeps its local H/J (Eddington)
! content; the ratio is 1 on every fallback branch above, so H is
! untouched whenever J was not actually corrected.
!
	              IF(T4 .GT. 0.0_LDP .AND. XV .GT. 0.0_LDP)
	1                H_NEXT_FRQ(I,L)=H_NEXT_FRQ(I,L)*XV/T4
	            END DO
	          END DO
	          WRITE(LUER,'(A,I2,A,ES12.4,A,ES12.4,A,I12)')
	1            '    AA m=',MP_AA,'  ||f||=',FNRM,
	2            '  ||corr||=',DNRM*SC,'  excl=',NAA_SKIP
	          END IF
	        ELSE
	          WRITE(LUER,'(A,I3,A,I2)')
	1            '    AA solve failed (DGESV info=',INFO_AA,
	2            ') -- raw map output kept'
	        END IF
	      END IF
	    END IF
	    JPREV4(1:ND_NM,1:KN_M_NCF)=JPREV3(1:ND_NM,1:KN_M_NCF)
	    JPREV3(1:ND_NM,1:KN_M_NCF)=JPREV2(1:ND_NM,1:KN_M_NCF)
	    JPREV2(1:ND_NM,1:KN_M_NCF)=JPREV1(1:ND_NM,1:KN_M_NCF)
	    JPREV1(1:ND_NM,1:KN_M_NCF)=J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	  END IF
!
	  IDIR=1
	  IF(LDI .LT. 0)IDIR=2
	  RATIO_M=-1.0_LDP
	  IF(RES_PREV(IDIR) .GT. 0.0_LDP)
	1  RATIO_M=PASS_CHG/RES_PREV(IDIR)
	  WRITE(LUER,'(/,A,I3,3A)')' ------------- MOM refresh pass',IPASS,
	1    '  dir=',MERGE('D','A',LDI.LT.0),
	2    ' --------------------------------------'
	  WRITE(LUER,'(4X,A11,2(1X,ES10.3),1X,A,I10,1X,A,I10)')
	1    'drift',PASS_CHG,RATIO_M,'n5=',PASS_N5,'n50=',PASS_N50
	  WRITE(LUER,'(4X,A11,1X,A,I5,A,I7,A,6(1X,ES10.3))')
	1    'argmax(d,f)','(',PASS_ARGI,',',PASS_ARGL,')',
	2    PASS_NU,PASS_R,PASS_V,PASS_SIGT,PASS_GAM,PASS_Q
	  WRITE(LUER,'(4X,A11,1X,A,I10,1X,A,I10,1X,A,I10,1X,A,I11,A,I10,'//
	1    '1X,A,I6,1X,A,I6)')
	2    'repairs','negJ=',PASS_NEG,'quar=',PASS_NQ,
	3    'cap=',PASS_CAP,'match=',PASS_NMATCH,'/',PASS_NCMP,
	4    'off=',PASS_NOFF,'nmax=',PASS_NMAX
	  IF(MOM_EPS_CAP .GT. 0_8)WRITE(LUER,'(4X,A11,1X,A,I10)')
	1    'eps-limit','capped=',MOM_EPS_CAP
	  IF(MOM_RF_NSUPP .GT. 0_8)WRITE(LUER,'(4X,A11,1X,A,I10)')
	1    'owned-skip','suppressed=',MOM_RF_NSUPP
	  IF(ALLOCATED(MOM_DAMP_NM))THEN
	    NDMP=COUNT(MOM_DAMP_NM(1:ND_NM,1:KN_M_NCF) .LT.
	1      0.99_LDP*JNEXT_RELAX)
	    IF(NDMP .GT. 0 .OR. NACC .GT. 0)WRITE(LUER,
	1      '(4X,A11,1X,A,I10,1X,A,ES10.3,1X,A,I10)')
	2      'damp/accel','damped=',NDMP,'min relax=',
	3      MINVAL(MOM_DAMP_NM(1:ND_NM,1:KN_M_NCF)),
	4      'accel=',NACC
	  END IF
	  WRITE(LUER,'(A)')' -------------------------------------------------------------------------------'
	  IF(PASS_CHG .LT. KN_M_TOL)THEN
	    CONVERGED=.TRUE.
	    WRITE(LUER,'(A,I4,A,ES10.3,A)')
	1      ' MOM_KNOP_REFRESH: converged at pass',IPASS,
	2      '  residual=',PASS_CHG,'  [OK]'
	    EXIT
	  END IF
!
!	  Bulk-converged exit (MOM_RF_BULK_EXIT): the max-norm above
!	  asymptotes on marginally contractive compressing cells, so
!	  without this the refresh burns passes until the stall exit
!	  retains the same bulk-settled field.  The bounds match the
!	  gate's sweep criteria exactly; IPASS>=4 guarantees both sweep
!	  directions ran at least twice.  CONVERGED=T routes the gate to
!	  the retained-vs-sweep delta, which still rejects a fixed point
!	  sitting far from the field the populations saw.
!
	  IF(MOM_RF_BULK_EXIT .AND. IPASS .GE. 4 .AND.
	1     PASS_CAP .EQ. 0_8 .AND.
	2     LAST_N50(1)+LAST_N50(2) .LE.
	3     MAX(50_8,MOM_SW_NTOT/1000) .AND.
	4     PASS_NEG .LE. MOM_SW_NTOT/10 .AND.
	5     PASS_NQ .LE. MOM_SW_NTOT/10)THEN
	    CONVERGED=.TRUE.
	    WRITE(LUER,'(A,I4,A,I12,A)')
	1      ' MOM_KNOP_REFRESH: bulk-converged at pass',IPASS,
	2      '  cells>50%=',LAST_N50(1)+LAST_N50(2),'  [OK]'
	    EXIT
	  END IF
!
!	  Abort once the residual has exploded orders of magnitude above
!	  the best achieved -- further passes only amplify the unstable
!	  mode (observed: contraction for ~20 passes then a sudden
!	  4-dex blow-up on job 11888720).  The reference is MIN_RES, the
!	  smallest max-norm seen, since BEST_CHG now tracks the
!	  bulk-chosen pass.
!
	  IF(IPASS .GE. 3 .AND. PASS_CHG .GT. 1.0E+3_LDP*MIN_RES)THEN
	    WRITE(LUER,'(A,I4,A,ES12.4,A,ES12.4,A)')
	1      ' MOM_KNOP_REFRESH: diverging -- abort at pass',IPASS,
	2      ', residual=',PASS_CHG,' vs min=',MIN_RES,'  [BAD]'
	    EXIT
	  END IF
!
!	  Same-direction stall check (KNOP_RAY_SOLVE_V1 pattern): the two
!	  upwind directions alternate, so compare against the residual of
!	  the same-direction pass two iterations back.  Two consecutive
!	  same-direction ratios above 0.90 -- once past the transient --
!	  mean the map is cycling rather than contracting; exit instead
!	  of burning the remaining passes.  The field is still restored
!	  below because CONVERGED stays false.
!
!	  edited by Devin, 30-Sep-2026: the transient floor grows with
!	  the AA history depth (2*MOM_AA_M+4) when MOM_RF_AA is on.  AA
!	  only sees one composite-map iterate per even pass boundary, so
!	  an exit at pass ~8-10 starves the Krylov window it needs to
!	  collapse the marginal (spectral-radius ~1) eigenspace that the
!	  stall check detects -- the detector cannot see contraction AA
!	  has not yet had history to produce.
!
	  IDIR=1
	  IF(LDI .LT. 0)IDIR=2
	  N_FLOOR=8
	  IF(MOM_RF_AA)N_FLOOR=MAX(8,2*MOM_AA_M+4)
	  IF(RES_PREV(IDIR) .GT. 0.0_LDP)THEN
	    RATIO_M=PASS_CHG/RES_PREV(IDIR)
	    IF(IPASS .GE. N_FLOOR)THEN
	      IF(RATIO_M .GT. 0.90_LDP)THEN
	        NSTALL_D(IDIR)=NSTALL_D(IDIR)+1
	      ELSE
	        NSTALL_D(IDIR)=0
	      END IF
	    END IF
	  END IF
	  RES_PREV(IDIR)=PASS_CHG
!
!	  Per-direction stall bookkeeping: under MOM_RF_OWNW the two
!	  directions own disjoint row subsets, and the previous shared
!	  counter let the stalled direction abort the whole refresh
!	  while the other was still contracting.  A stalled direction is
!	  dropped from the schedule (its best rows are already in the
!	  retained composite); the refresh exits only when every active
!	  direction has stalled.  [MOM_RF_OWNW]
!
	  IF(NSTALL_D(IDIR) .GE. 2)THEN
	    IF(.NOT. SKIP_D(IDIR))WRITE(LUER,'(3A,I4,A)')
	1      ' MOM_KNOP_REFRESH: direction stalled -- dropping dir=',
	2      MERGE('D','A',IDIR.EQ.2),' after pass',IPASS,'  [WARN]'
	    SKIP_D(IDIR)=.TRUE.
	  END IF
	  IF(MOM_DOM_CMP)THEN
	    IF(SKIP_D(2))THEN
	      WRITE(LUER,'(A,I4,A,ES12.4,A)')
	1        ' MOM_KNOP_REFRESH: stalled -- exit at pass',IPASS,
	2        ', residual=',PASS_CHG,'  [WARN]'
	      EXIT
	    END IF
	  ELSE IF(MOM_DOM_EXP)THEN
	    IF(SKIP_D(1))THEN
	      WRITE(LUER,'(A,I4,A,ES12.4,A)')
	1        ' MOM_KNOP_REFRESH: stalled -- exit at pass',IPASS,
	2        ', residual=',PASS_CHG,'  [WARN]'
	      EXIT
	    END IF
	  ELSE IF(SKIP_D(1) .AND. SKIP_D(2))THEN
	    WRITE(LUER,'(A,I4,A,ES12.4,A)')
	1      ' MOM_KNOP_REFRESH: stalled -- exit at pass',IPASS,
	2      ', residual=',PASS_CHG,'  [WARN]'
	    EXIT
	  END IF
	END DO
!
	MOM_KN_REFRESHING=.FALSE.
	KN_M_DESC=.FALSE.
!
! If the refresh did not converge, keep the best-pass iterate when
! its residual improved on the first pass: the stored lagged field
! then retains the wavelength propagation the refresh produced
! instead of reverting to the pre-refresh field and advancing one
! frequency bin per great iteration.  A pass sequence that never
! improved on the seed restores the pre-refresh field instead -- the
! iterated field is then strictly worse than what the sweep saw.
! Every write inside the passes is already per-cell capped and
! NaN-guarded, so the best-pass field is bounded by the same guards
! as the in-sweep store.  Report the best residual as the sweep
! metric -- it measures how far the iterated field drifted from the
! field actually retained.
!
	IF(.NOT. CONVERGED)THEN
!
! edited by Devin, 28-Sep-2026: the residual-improvement test alone
! (BEST_CHG < RESID0) is nearly always true since first-pass
! residuals run 1E4-3E6, and it actively misfires once the best pass
! is chosen by bulk counts -- the bulk-best pass can carry a larger
! max-norm than pass 1 while being far better in aggregate.  The
! candidate is accepted on the same lexicographic order the selector
! uses: n50 (the worst-cell count) must strictly improve, or hold
! while n5 improves.  Requiring n5 non-regression jointly discarded
! passes that fixed the worst cells while stirring the quiet bulk
! slightly -- the store then reverted unchanged (delta=0) and the
! same stall repeated every iteration.  An n5 regression is still
! bounded at 1.5x so a pass cannot trade a few bad cells for broad
! churn, and the max-norm must stay below the same 1E3*MIN_RES
! runaway bound the abort uses.
!
	  KEEP_BEST=BEST_VALID .AND. (BEST_N50 .LT. N50_SV .OR.
	1    (BEST_N50 .LE. N50_SV .AND. BEST_N5 .LE. N5_SV)) .AND.
	2    BEST_N5 .LE. N5_SV+N5_SV/2 .AND.
	3    BEST_CHG .LE. 1.0E+3_LDP*MIN_RES
	  IF(BEST_VALID .AND. .NOT. KEEP_BEST)THEN
	    IF(BEST_N50 .GT. N50_SV)THEN
	      WRITE(LUER,'(A,2I12,A,2I12,2A)')
	1        ' MOM_KNOP_REFRESH: best-pass field rejected'//
	2        ' (n50 regressed: n5,n50=',BEST_N5,BEST_N50,
	3        ' vs sweep',N5_SV,N50_SV,')','  [WARN]'
	    ELSE IF(BEST_N5 .GT. N5_SV+N5_SV/2)THEN
	      WRITE(LUER,'(A,2I12,A,2I12,2A)')
	1        ' MOM_KNOP_REFRESH: best-pass field rejected'//
	2        ' (n5 churn >1.5x: n5,n50=',BEST_N5,BEST_N50,
	3        ' vs sweep',N5_SV,N50_SV,')','  [WARN]'
	    ELSE IF(BEST_N5 .GT. N5_SV)THEN
	      WRITE(LUER,'(A,2I12,A,2I12,2A)')
	1        ' MOM_KNOP_REFRESH: best-pass field rejected'//
	2        ' (no bulk improvement: n5,n50=',BEST_N5,BEST_N50,
	3        ' vs sweep',N5_SV,N50_SV,')','  [WARN]'
	    ELSE
	      WRITE(LUER,'(A,ES12.4,A,ES12.4,2A)')
	1        ' MOM_KNOP_REFRESH: best-pass field rejected'//
	2        ' (residual runaway: ',BEST_CHG,' vs 1E3*min=',MIN_RES
	3        ,')','  [WARN]'
	    END IF
	  END IF
	  IF(KEEP_BEST)THEN
	    DO L=1,KN_M_NCF
	      DO I=1,ND_NM
	        IF(.NOT. (ABS(JBEST(I,L)) .LT. 1.0E40_LDP) .OR.
	1          .NOT. (ABS(HBEST(I,L)) .LT. 1.0E40_LDP))
	2          KEEP_BEST=.FALSE.
	      END DO
	    END DO
	    IF(.NOT. KEEP_BEST)WRITE(LUER,'(A)')
	1      ' MOM_KNOP_REFRESH: best-pass field rejected'//
	2      ' (non-finite or >1E40 entry)  [WARN]'
	  END IF
	  IF(KEEP_BEST)THEN
	    J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)=
	1      JBEST(1:ND_NM,1:KN_M_NCF)
	    H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)=
	1      HBEST(1:ND_NM,1:KN_M_NCF)
	    WRITE(LUER,'(A,ES12.4,A,2I12,A)')
	1      ' MOM_KNOP_REFRESH: unconverged -- retained best-pass'//
	2      ' field, residual=',BEST_CHG,'  n5,n50=',BEST_N5,BEST_N50,
	3      '  [WARN]'
	  ELSE
	    J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)=
	1      JNEXT_SAVE(1:ND_NM,1:KN_M_NCF)
	    H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)=
	1      HNEXT_SAVE(1:ND_NM,1:KN_M_NCF)
	    WRITE(LUER,'(A,ES12.4,A)')
	1      ' MOM_KNOP_REFRESH: unconverged -- retained pre-refresh'//
	2      ' field, best residual=',BEST_CHG,'  [WARN]'
	  END IF
	  MOM_GS_MAXCHG=BEST_CHG
	  MOM_GS_NBAD5=N5_SV
	  MOM_GS_NBAD50=N50_SV
	  MOM_GS_NQUAR=NQ_SV
!
!	  PASS_CAP counts capped writes made by passes whose field was
!	  discarded; the metric must describe the retained field.  Recount
!	  non-finite or extreme (>=1E40, legacy-poisoned) entries in the
!	  retained field so a genuinely poisoned store still trips the
!	  fatal check.
!
	  MOM_INEXT_CAPPED=0
	  DO L=1,KN_M_NCF
	    DO I=1,ND_NM
	      IF(.NOT. (ABS(J_NEXT_FRQ(I,L)) .LT. 1.0E40_LDP))
	1        MOM_INEXT_CAPPED=MOM_INEXT_CAPPED+1
	    END DO
	  END DO
	END IF
!
! edited by Devin, 01-Oct-2026: call-level (outer) Anderson
! acceleration.  The committed lagged store is the iterate of the
! outer fixed-point map (sweep+refresh composition, including the
! J<->Eddington coefficient drift); the inner pass AA cannot see
! that drift because it accelerates a frozen-coefficient map.  The
! observed failure it answers is the outer limit cycle: on job
! 12047456 the retained-vs-entry delta grew to ~2E6 cells>50% while
! the inner residual floored at ~2E4 and FG maxchg plateaued at
! ~2E-3 -- the marginal-mode subspace wanders between calls because
! no accelerator acts on the committed iterate sequence.  Residual
! f_{n-1} = committed_out - committed_in is paired with the field
! the map was applied to (RETJ_PREV); the secant solve, weighting,
! GCAP exclusion, ridge, ||corr||>||f|| rejection, per-cell cap and
! H-ratio drag are the same as the pass-level AA.  When the raw map
! returned the input unchanged (KEEP_BEST=F) f=0 and the update is
! a no-op.  [MOM_AA_OUT] gates it (default T); history depth
! [MOM_AA_OUT_M] (default 8).
!
	IF(MOM_RF_AA .AND. MOM_AA_OUT)THEN
	  IF(.NOT. ALLOCATED(XHIST_OUT) .OR.
	1     SIZE(XHIST_OUT,1) .NE. ND_NM .OR.
	2     SIZE(XHIST_OUT,2) .NE. KN_M_NCF .OR.
	3     SIZE(XHIST_OUT,3) .NE. MOM_AA_OUT_M+1)THEN
	    IF(ALLOCATED(XHIST_OUT))DEALLOCATE(XHIST_OUT)
	    IF(ALLOCATED(FHIST_OUT))DEALLOCATE(FHIST_OUT)
	    IF(ALLOCATED(RETJ_PREV))DEALLOCATE(RETJ_PREV)
	    IF(ALLOCATED(RETH_PREV))DEALLOCATE(RETH_PREV)
	    ALLOCATE(XHIST_OUT(ND_NM,KN_M_NCF,MOM_AA_OUT_M+1),
	1      FHIST_OUT(ND_NM,KN_M_NCF,MOM_AA_OUT_M+1),
	2      RETJ_PREV(ND_NM,KN_M_NCF),RETH_PREV(ND_NM,KN_M_NCF),
	3      STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      IF(ALLOCATED(XHIST_OUT))DEALLOCATE(XHIST_OUT)
	      IF(ALLOCATED(FHIST_OUT))DEALLOCATE(FHIST_OUT)
	      IF(ALLOCATED(RETJ_PREV))DEALLOCATE(RETJ_PREV)
	      IF(ALLOCATED(RETH_PREV))DEALLOCATE(RETH_PREV)
	      WRITE(LUER,*)'Unable to allocate outer AA history'//
	1        ' in MOM_KNOP_REFRESH_V1'
	    END IF
	    NHIST_OUT=-1
	  END IF
	  IF(ALLOCATED(XHIST_OUT))THEN
!
! RETJ_PREV is valid only after the first call commits a field
! (NHIST_OUT=-1 before that): the first push then pairs (x_1, f_1)
! with no phantom seed pair.
!
	    IF(NHIST_OUT .GE. 0)THEN
	      DO JJ=MIN(NHIST_OUT+1,MOM_AA_OUT_M+1),2,-1
	        XHIST_OUT(1:ND_NM,1:KN_M_NCF,JJ)=
	1          XHIST_OUT(1:ND_NM,1:KN_M_NCF,JJ-1)
	        FHIST_OUT(1:ND_NM,1:KN_M_NCF,JJ)=
	1          FHIST_OUT(1:ND_NM,1:KN_M_NCF,JJ-1)
	      END DO
	      XHIST_OUT(1:ND_NM,1:KN_M_NCF,1)=
	1        RETJ_PREV(1:ND_NM,1:KN_M_NCF)
	      FHIST_OUT(1:ND_NM,1:KN_M_NCF,1)=
	1        J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)-
	2        RETJ_PREV(1:ND_NM,1:KN_M_NCF)
	      NHIST_OUT=MIN(NHIST_OUT+1,MOM_AA_OUT_M+1)
	      MP_AA=NHIST_OUT-1
	      IF(MP_AA .GE. 1)THEN
	        GRM_OUT(1:MP_AA,1:MP_AA)=0.0_LDP
	        BRHS_OUT(1:MP_AA)=0.0_LDP
	        NAA_SKIP=0_8
	        DO L=1,KN_M_NCF
	          TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1            1.0E-28_LDP)
	          DO I=1,ND_NM
	            WV=1.0_LDP/MAX(ABS(XHIST_OUT(I,L,1)),
	1              1.0E-4_LDP*TCOL)
	            SC=FHIST_OUT(I,L,1)*WV
	            T1X=ABS(SC)
	            DO JJ=1,MP_AA
	              T1X=MAX(T1X,ABS(FHIST_OUT(I,L,1)-
	1                FHIST_OUT(I,L,JJ+1))*WV)
	            END DO
	            IF(T1X .GT. MOM_AA_GCAP)THEN
	              NAA_SKIP=NAA_SKIP+1
	              CYCLE
	            END IF
	            DO JJ=1,MP_AA
	              XV=(FHIST_OUT(I,L,1)-FHIST_OUT(I,L,JJ+1))*WV
	              BRHS_OUT(JJ)=BRHS_OUT(JJ)+XV*SC
	              DO KK=1,JJ
	                GRM_OUT(KK,JJ)=GRM_OUT(KK,JJ)+
	1                 XV*(FHIST_OUT(I,L,1)-
	2                 FHIST_OUT(I,L,KK+1))*WV
	              END DO
	            END DO
	          END DO
	        END DO
	        DO JJ=1,MP_AA
	          DO KK=JJ+1,MP_AA
	            GRM_OUT(JJ,KK)=GRM_OUT(KK,JJ)
	          END DO
	        END DO
	        DTR=0.0_LDP
	        DO JJ=1,MP_AA
	          DTR=MAX(DTR,GRM_OUT(JJ,JJ))
	        END DO
	        DTR=1.0E-8_LDP*MAX(DTR,1.0E-30_LDP)
	        DO JJ=1,MP_AA
	          GRM_OUT(JJ,JJ)=GRM_OUT(JJ,JJ)+DTR
	        END DO
	        WGRM_OUT(1:MP_AA,1:MP_AA)=GRM_OUT(1:MP_AA,1:MP_AA)
	        WBR_OUT(1:MP_AA)=BRHS_OUT(1:MP_AA)
	        IF(ND_NM*INT(KN_M_NCF,8)-NAA_SKIP .LT.
	1           MAX(100_8,ND_NM*INT(KN_M_NCF,8)/1000))THEN
	          INFO_AA=-999
	          WRITE(LUER,'(A,I12)')
	1            '    outer-AA skipped -- affine subset empty, excl=',
	2            NAA_SKIP
	        ELSE
	          CALL DGESV(MP_AA,1,WGRM_OUT,MOM_AA_OUT_M,
	1            IPIV_OUT,WBR_OUT,MOM_AA_OUT_M,INFO_AA)
	        END IF
	!
	! Contractive-retry: if the full-depth correction is rejected
	! (||corr||>||f||, near-singular Gram) re-solve on the leading m-1
	! window -- the secant columns are newest-first, so this drops the
	! oldest, least relevant direction.  A bad deep solve degrades to
	! the largest contractive window instead of wasting the call.
	!
	        MTRY=MP_AA
	        DO WHILE(INFO_AA .EQ. 0 .AND. MTRY .GE. 1)
	          IF(MTRY .LT. MP_AA)THEN
	            WGRM_OUT(1:MTRY,1:MTRY)=GRM_OUT(1:MTRY,1:MTRY)
	            WBR_OUT(1:MTRY)=BRHS_OUT(1:MTRY)
	            CALL DGESV(MTRY,1,WGRM_OUT,MOM_AA_OUT_M,
	1              IPIV_OUT,WBR_OUT,MOM_AA_OUT_M,INFO_AA)
	            IF(INFO_AA .NE. 0)EXIT
	          END IF
	          GAM_OUT(1:MTRY)=WBR_OUT(1:MTRY)
	          FNRM=0.0_LDP;  DNRM=0.0_LDP
	          DO L=1,KN_M_NCF
	            TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1              1.0E-28_LDP)
	            DO I=1,ND_NM
	              WV=1.0_LDP/MAX(ABS(XHIST_OUT(I,L,1)),
	1                1.0E-4_LDP*TCOL)
	              IF(ABS(FHIST_OUT(I,L,1))*WV .GT. MOM_AA_GCAP)
	1                CYCLE
	              FNRM=FNRM+(FHIST_OUT(I,L,1)*WV)**2
	              SC=0.0_LDP
	              DO JJ=1,MTRY
	                SC=SC+GAM_OUT(JJ)*(J_NEXT_FRQ(I,L)-
	1                 XHIST_OUT(I,L,JJ+1)-FHIST_OUT(I,L,JJ+1))*WV
	              END DO
	              DNRM=DNRM+SC*SC
	            END DO
	          END DO
	          FNRM=SQRT(FNRM);  DNRM=SQRT(DNRM)
	          IF(DNRM .LE. FNRM)EXIT
	          IF(MTRY .EQ. MP_AA)THEN
	            WRITE(LUER,'(A,I2,A,ES12.4,A,ES12.4,A,I12)')
	1              '    outer-AA m=',MTRY,'  ||f||=',FNRM,
	2              '  ||corr||=',DNRM,'  retry at m-1, excl=',
	3              NAA_SKIP
	          END IF
	          MTRY=MTRY-1
	          INFO_AA=0
	        END DO
	        IF(INFO_AA .EQ. 0 .AND. MTRY .GE. 1)THEN
	          DO L=1,KN_M_NCF
	            TCOL=MAX(MAXVAL(ABS(J_NEXT_FRQ(1:ND_NM,L))),
	1              1.0E-28_LDP)
	            DO I=1,ND_NM
	              WV=1.0_LDP/MAX(ABS(XHIST_OUT(I,L,1)),
	1                1.0E-4_LDP*TCOL)
	              IF(ABS(FHIST_OUT(I,L,1))*WV .GT. MOM_AA_GCAP)
	1                CYCLE
	              XV=0.0_LDP
	              DO JJ=1,MTRY
	                XV=XV+GAM_OUT(JJ)*(J_NEXT_FRQ(I,L)-
	1                   XHIST_OUT(I,L,JJ+1)-FHIST_OUT(I,L,JJ+1))
	              END DO
	              T4=J_NEXT_FRQ(I,L)
	              XV=T4-XV
	              T3=1.0E+1_LDP*MAX(ABS(FHIST_OUT(I,L,1)),
	1                  1.0E-4_LDP*TCOL)
	              STEP=T4-XV
	              IF(ABS(STEP) .GT. T3)
	1                  XV=T4-SIGN(T3,STEP)
	              IF(XV .LE. 0.0_LDP .AND. T4 .GT. 0.0_LDP)
	1                  XV=T4
	              J_NEXT_FRQ(I,L)=XV
	              IF(T4 .GT. 0.0_LDP .AND. XV .GT. 0.0_LDP)
	1                  H_NEXT_FRQ(I,L)=H_NEXT_FRQ(I,L)*XV/T4
	            END DO
	          END DO
	          WRITE(LUER,'(A,I2,A,ES12.4,A,ES12.4,A,I12)')
	1              '    outer-AA m=',MTRY,'  ||f||=',FNRM,
	2              '  ||corr||=',DNRM,'  excl=',NAA_SKIP
	        ELSE IF(INFO_AA .EQ. 0)THEN
	          WRITE(LUER,'(A,ES12.4,A,I12)')
	1            '    outer-AA all depths rejected, ||f||=',FNRM,
	2            '  excl=',NAA_SKIP
	        ELSE IF(INFO_AA .NE. -999)THEN
	          WRITE(LUER,'(A,I3,A)')
	1            '    outer-AA solve failed (DGESV info=',INFO_AA,
	2            ') -- raw map output kept'
	        END IF
	      END IF
	    END IF
	    RETJ_PREV(1:ND_NM,1:KN_M_NCF)=
	1      J_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	    RETH_PREV(1:ND_NM,1:KN_M_NCF)=
	1      H_NEXT_FRQ(1:ND_NM,1:KN_M_NCF)
	    IF(NHIST_OUT .LT. 0)NHIST_OUT=0
	  END IF
	END IF
!
! Post-refresh maturity of the RETAINED field: the distance between
! the field the populations saw (the pre-refresh store, which is the
! sweep's relaxed output) and the store actually carried forward.
! The normalization matches the sweep metric: per-frequency column
! peak with a 1E-4 floor so negligible cells cannot pin the max.
! When the refresh did not converge and the pre-refresh field was
! retained, the delta is identically zero; a retained best-pass
! field instead reports the real distance it moved from the field
! the populations saw.  MOM_RF_CONV=F tells SOLVE_FOR_POPS to keep
! gating on the sweep drift either way.
!
	MOM_RF_RAN=.TRUE.
	MOM_RF_CONV=CONVERGED
	MOM_RF_DCHG=0.0_LDP
	MOM_RF_N5=0
	MOM_RF_N50=0
	DO L=1,KN_M_NCF
	  T3=MAX(MAXVAL(ABS(JNEXT_SAVE(1:ND_NM,L))),1.0E-28_LDP)
	  DO I=1,ND_NM
!
! edited by Devin, 28-Sep-2026: the floor test must see the field on
! both sides of the delta -- a cell tiny in JNEXT_SAVE but grown
! orders of magnitude in the retained field was skipped entirely,
! hiding exactly the blow-up this metric is meant to expose.
!
	    IF(MAX(ABS(JNEXT_SAVE(I,L)),ABS(J_NEXT_FRQ(I,L))) .GT.
	1       1.0E-4_LDP*T3)THEN
	      T4=ABS(J_NEXT_FRQ(I,L)-JNEXT_SAVE(I,L))/
	1         MAX(ABS(JNEXT_SAVE(I,L)),1.0E-4_LDP*T3)
	      IF(T4 .GT. 0.05_LDP)MOM_RF_N5=MOM_RF_N5+1
	      IF(T4 .GT. 0.5_LDP)MOM_RF_N50=MOM_RF_N50+1
	      MOM_RF_DCHG=MAX(MOM_RF_DCHG,T4)
	    END IF
	  END DO
	END DO
	RF_ST='[WARN] '
	IF(MOM_RF_DCHG .LT. 0.5_LDP .AND.
	1   MOM_RF_N50 .LE. MAX(50_8,MOM_SW_NTOT/1000))RF_ST='[OK]   '
	WRITE(LUER,'(A,ES12.4,A,I12,A,I12,2X,A)')
	1  ' MOM_KNOP_REFRESH retained-field delta: max=',MOM_RF_DCHG,
	2  '  cells>5%=',MOM_RF_N5,'  cells>50%=',MOM_RF_N50,RF_ST
!
! Restore the post-sweep module state (frequency NCF) so the next
! great iteration starts from the same external view as a normal
! sweep exit.
!
	L=KN_M_NCF
	RSQJNU_NM(1:ND_NM)=J_NEXT_FRQ(1:ND_NM,L)
	RSQHNU_NM(1:ND_NM)=HRF(1:ND_NM,L)
	K_ON_J(1:ND_NM)=KN_MKONJ(1:ND_NM,L)
	NMID_ON_J(1:ND_NM)=KN_MNMJ(1:ND_NM,L)
	NMID_ON_HMID(1:ND_NM)=KN_MNMH(1:ND_NM,L)
	DO I=1,ND_NM
	  IF(.NOT. (ABS(K_ON_J(I)) .LE. HUGE(K_ON_J(I))))
	1      K_ON_J(I)=1.0_LDP/3.0_LDP
	  K_ON_J(I)=MIN(1.0_LDP,MAX(1.0E-4_LDP,K_ON_J(I)))
	END DO
	DO I=1,ND_NM-1
	  IF(.NOT. (ABS(NMID_ON_J(I)) .LE. HUGE(NMID_ON_J(I))))
	1      NMID_ON_J(I)=0.0_LDP
	  IF(.NOT. (ABS(NMID_ON_HMID(I)) .LE. HUGE(NMID_ON_HMID(I))))
	1      NMID_ON_HMID(I)=0.0_LDP
	  NMID_ON_HMID(I)=MIN(1.0_LDP,MAX(0.0_LDP,NMID_ON_HMID(I)))
	END DO
	K_ON_J_SAVE(1:ND_NM)=K_ON_J(1:ND_NM)
	NMID_ON_J_SAVE(1:ND_NM)=NMID_ON_J(1:ND_NM)
	NMID_ON_HMID_SAVE(1:ND_NM)=NMID_ON_HMID(1:ND_NM)
	HBC=KN_MHBC(L);  NBC=KN_MNBC(L);  IN_HBC=KN_MINHBC(L)
	HBC_SAVE=HBC;  NBC_SAVE=NBC;  IN_HBC_SAVE=IN_HBC
	JPLUS_IB=KN_MJPIB(L);  JMIN_IB=KN_MJMIB(L)
	HNU_AT_IB=KN_MHIB(L)
	MOM_NU_NEXT=KN_MNUNXT(L)
!
	WRITE(LUER,'(/,A,I4,A,F8.2,A)')
	1  ' ------------- MOM_KNOP_REFRESH done:',MIN(IPASS,KN_M_MAXPASS),
	2  ' passes,',OMP_GET_WTIME()-T0,' s --------------------'
!
	DEALLOCATE(HRF,OUTJ,OUTH,ETA_L,CHI_L,ESEC_L)
	DEALLOCATE(JNEXT_RF_READ,JNEXT_SAVE)
	DEALLOCATE(HNEXT_RF_READ,HNEXT_SAVE)
	DEALLOCATE(JPREV1,JPREV2)
	DEALLOCATE(JPREV3,JPREV4)
	DEALLOCATE(JBEST,HBEST)
	IF(ALLOCATED(XHIST_AA))DEALLOCATE(XHIST_AA)
	IF(ALLOCATED(FHIST_AA))DEALLOCATE(FHIST_AA)
	IF(ALLOCATED(ACT_CUR))DEALLOCATE(ACT_CUR)
	IF(ALLOCATED(ACT_NEW))DEALLOCATE(ACT_NEW)
	IF(ALLOCATED(SHOT_L))DEALLOCATE(SHOT_L)
	IF(ALLOCATED(RES_L))DEALLOCATE(RES_L)
	IF(ALLOCATED(NMUP_L))DEALLOCATE(NMUP_L)
	RETURN
	END
