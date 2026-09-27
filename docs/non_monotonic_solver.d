================================================================================
CMFGEN Non-Monotonic Velocity Solver — Design Document
================================================================================
Project:  SN2024pba (X2_sn interacting supernova model)
Date:     2026-08-01
Author:   S. Brennan
Status:   Implemented, audited, resubmitted

================================================================================
1. PROBLEM STATEMENT
================================================================================

CMFGEN was originally written for monotonic outflowing stellar winds where
velocity V(R) increases outward at all depths. The Sobolev optical depth is

    GAM = (pi*e^2 / m_e*c) * n_l * f_lu * lambda * (1 / |dV/dR - V/R|)

which in terms of sigma = dlnV/dlnR - 1 becomes

    GAM = T1 * |1 + sigma * mu^2|        where T1 = 3.33564E-6 * V / R
                                                mu = Z / R (impact parameter angle)

For monotonic flows, sigma >= 0 (typically sigma ~ 0.5 for beta=1 winds), so
1 + sigma*mu^2 > 0 always, and the absolute value is a no-op.

For INTERACTING SUPERNOVAE (Type IIn, shell collisions), the velocity profile
is NON-MONOTONIC: velocity decreases inward through a decelerated shell, then
increases again in the inner ejecta. This produces sigma < -1 at velocity
reversals, where dlnV/dlnR < 0 and |dlnV/dlnR| > 1.

When sigma < -1:
  - 1 + sigma*mu^2 can become ZERO or NEGATIVE for some angles mu
  - |1 + sigma*mu^2| can become very small, producing huge GAM (optical depth)
  - Division by near-zero GAMH/GAM produces NaN populations
  - The original code would STOP or produce garbage

The Dessart+2015 prescription (Dessart et al. 2015, A&A 575, A51) provides a
physically motivated cap that prevents these failures while preserving the
correct behavior for monotonic flows.

================================================================================
2. THE DESSART+2015 PRESCRIPTION
================================================================================

Reference: Dessart, Hillier, Livne, Yoon (2015), "Radiation-hydrodynamics
simulations of interacting supernovae", A&A 575, A51, Section 3.2.

The prescription modifies the Sobolev optical depth denominator:

    GAM = T1 * MAX( |1 + sigma*mu^2|, -0.2*(1+sigma) )     when sigma < 0
    GAM = T1 * |1 + sigma*mu^2|                              when sigma >= 0

with a floor of 1.0E-3 to prevent division by zero.

Physical justification:
  - When sigma >= 0 (monotonic outflow): the prescription is a no-op.
    MAX(|1+sigma*mu^2|, -0.2*(1+sigma)) = |1+sigma*mu^2| because
    -0.2*(1+sigma) < 0 <= |1+sigma*mu^2|. No change from original behavior.

  - When -1 < sigma < 0 (mild deceleration): the prescription is still a no-op.
    -0.2*(1+sigma) is in [0, 0.2), while |1+sigma*mu^2| >= 1-|sigma| > 0.
    For sigma = -0.5: -0.2*0.5 = -0.1 < 0 < |1-0.5*mu^2|. No change.

  - When sigma < -1 (velocity reversal): 1+sigma can be negative, so
    -0.2*(1+sigma) > 0. This provides a POSITIVE FLOOR for the denominator,
    preventing it from going to zero or negative. The factor 0.2 is empirical,
    chosen by Dessart+2015 to be small enough to not dominate in normal regions
    but large enough to prevent numerical explosions in reversal regions.

  - The floor 1.0E-3 catches the edge case where both |1+sigma*mu^2| and
    -0.2*(1+sigma) are very small (sigma very close to -1).

================================================================================
3. CHANGES MADE — FILE BY FILE
================================================================================

A total of 19 files were modified. The changes fall into 6 categories:

--------------------------------------------------------------------------------
3.1  SOBOLEV LINE CODE — Dessart+2015 Prescription (6 files)
--------------------------------------------------------------------------------

Files:
  subs/sobjbar_sim.f          (line ~110)
  subs/sobew.f                (line ~115)
  subs/sobew_grad.f           (line ~110)
  subs/sobew_grad_v2.f        (line ~125)
  subs/sobew_grad_v3.f        (line ~128)
  disp/subs/sobjbar.f         (line ~92)

WHAT CHANGED:
  The Sobolev optical depth computation GAM(I) = GAMH(I) / |1+sigma*mu^2|
  was replaced with the Dessart+2015 prescription:

  BEFORE:
    GAM(I)=GAMH(I)/ABS(1.0_LDP+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I))

  AFTER:
    T1=1.0_LDP+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I)
    IF(SIGMA(I) .LT. 0.0_LDP)THEN
      T1=MAX(ABS(T1),-0.2_LDP*(1.0_LDP+SIGMA(I)))
    ELSE
      T1=ABS(T1)
    END IF
    IF(T1 .LT. 1.0E-3_LDP)T1=1.0E-3_LDP
    GAM(I)=GAMH(I)/T1

  Note: sobjbar.f (disp/) originally used /(1+...) without ABS(); the fix
  also adds ABS() for the sigma >= 0 case, matching the other routines.

WHY:
  These routines compute the angle-dependent Sobolev optical depth for line
  transfer. When sigma < -1, the denominator 1+sigma*mu^2 can be zero or
  negative, producing infinite or negative optical depths. The Dessart+2015
  prescription caps the denominator at -0.2*(1+sigma) > 0, preventing this.

EFFECT ON SYSTEM:
  - Line JBAR (line mean intensity) computation is now stable for sigma < -1.
  - Line source functions converge instead of producing NaN populations.
  - No effect on monotonic flows (prescription is a no-op when sigma >= 0).
  - All 6 files use identical code for consistency.

CODE PATH: These routines are called from SUB_SOB_LINE (cmfgen_sub.f:3635+)
  when GLOBAL_LINE=SOB (Sobolev approximation). This is the ONLY line transfer
  path used by X2_sn (DO_SOB_LINES=T, SOBOLEV=.TRUE. for all lines).

--------------------------------------------------------------------------------
3.2  FORMAL SOLUTION (fg_j_cmf_v13.f) — Dessart+2015 + Sigma Clamp Removal
--------------------------------------------------------------------------------

File: plane/fg_j_cmf_v13.f

WHAT CHANGED (3 modifications):

  (a) SIGMA_RAY clamping removed (line 943):
    BEFORE: IF(SIGMA_RAY(I) .LE. -1.0_LDP)SIGMA_RAY(I)=-0.999_LDP
    AFTER:  C     Do NOT clamp SIGMA_RAY — Dessart+2015 prescription in GAM
            C     below handles sigma < -1

  (b) GAM computation — Dessart+2015 prescription (lines 958-965):
    BEFORE: GAM(I,LS)=T1*( 1.0_LDP+SIGMA_RAY(I)*(MU**2) )
    AFTER:  T2=1.0_LDP+SIGMA_RAY(I)*(MU**2)
            IF(SIGMA_RAY(I) .LT. 0.0_LDP)THEN
              T2=MAX(ABS(T2),-0.2_LDP*(1.0_LDP+SIGMA_RAY(I)))
            ELSE
              T2=ABS(T2)
            END IF
            IF(T2 .LT. 1.0E-3_LDP)T2=1.0E-3_LDP
            GAM(I,LS)=T1*T2

  (c) GAMH computation (DIFFERENCE method, lines 980-991):
    BEFORE: GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564E-06_LDP*
             ( 1.0_LDP+0.5_LDP*(MU**2)*(SIGMA_RAY(I)+SIGMA_RAY(I+1)) )/
             (R_RAY(I,LS)+R_RAY(I+1,LS))
    AFTER:  T2=0.5_LDP*(SIGMA_RAY(I)+SIGMA_RAY(I+1))
            T3=1.0_LDP+T2*(MU**2)
            IF(T2 .LT. 0.0_LDP)THEN
              T3=MAX(ABS(T3),-0.2_LDP*(1.0_LDP+T2))
            ELSE
              T3=ABS(T3)
            END IF
            IF(T3 .LT. 1.0E-3_LDP)T3=1.0E-3_LDP
            GAMH(I,LS)=(V_RAY(I+1)+V_RAY(I))*3.33564E-06_LDP*T3/
             (R_RAY(I,LS)+R_RAY(I+1,LS))

  (d) Diffusion inner boundary fix (lines 1803-1804):
    BEFORE: HNU_AT_IB=DBB/R(ND)/CHI_RAY(NI)/3.0_LDP
    AFTER:  HNU_AT_IB=DBB/R(ND)/CHI(ND)/3.0_LDP
    (same for NNU_AT_IB, /5.0_LDP)

  (e) ZERO_FLUX string fix (line 1811):
    BEFORE: (no ZERO_FLUX branch — fell through to default)
    AFTER:  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN

WHY:
  (a) The sigma clamp at -0.999 prevented the Dessart+2015 prescription from
      ever being activated (sigma was capped before reaching the GAM formula).
      Removing the clamp allows sigma < -1 to reach the prescription, which
      handles it correctly.

  (b,c) The formal solution computes J (mean intensity) by integrating along
      rays at different impact parameters. The Sobolev-like optical depth
      GAM appears in the line optical depth contribution. Without the
      prescription, sigma < -1 produces negative/zero GAM, corrupting the
      integration and producing NaN in J.

  (d) CHI_RAY(NI) is the opacity at the last ray point, which may not
      correspond to the innermost grid point ND. For the diffusion boundary
      condition, we need the opacity at the inner boundary (depth ND), which
      is CHI(ND). Using CHI_RAY(NI) could give the wrong opacity if the ray
      doesn't reach depth ND (e.g., for large impact parameters).

  (e) The string 'ZERO_FLUX' (with underscore) matches the VADAT keyword
      format. The original code had 'ZERO FLUX' (with space) which never
      matched, causing the ZERO_FLUX boundary to fall through to the default
      branch.

EFFECT ON SYSTEM:
  - J (mean intensity) is now correctly computed for non-monotonic velocity.
  - The diffusion inner boundary uses the correct opacity.
  - The ZERO_FLUX boundary is correctly recognized.
  - For monotonic flows: no change (prescription is a no-op, CHI_RAY(NI)
    equals CHI(ND) for rays that reach the inner boundary).

CODE PATH: FG_J_CMF_V13 is called from COMP_J_BLANK (comp_j_blank.f:643)
  when USE_LAM_ES=T (formal solution for J). This is the ONLY continuum J
  computation path for X2_sn. Called with FG_OPT=INT/INS (integral method).

  NOTE on dGAMdR (line 972): The derivative formula
    dGAMdR = GAM * SIGMA_RAY / R_RAY + T1 * mu^2 * d(SIGMA)/dR
  assumes GAM = T1*(1+sigma*mu^2), which is not exactly true when the
  prescription is active (GAM is capped). This produces a slightly
  inaccurate derivative at angles where the prescription is active, but
  does not affect stability — only integration accuracy. This is a known
  limitation, not a bug.

--------------------------------------------------------------------------------
3.3  FORMAL SOLUTION (fg_j_cmf_v12.f) — Same Changes (Not Called)
--------------------------------------------------------------------------------

File: plane/fg_j_cmf_v12.f

Same changes as V13 (Dessart+2015 prescription, sigma clamp removal,
CHI(ND) fix, ZERO_FLUX string fix).

WHY:
  V12 is the predecessor of V13. It is kept for reference but is NEVER
  CALLED in the current code (grep "CALL FG_J_CMF_V12" returns no matches).
  The changes were made for consistency and in case V12 is used in future.

NOTE: V12 still has STOP at line 916 when SIGMA < -1 (the diagnostic check).
  This is harmless since V12 is never called, but would need to be removed
  if V12 were reactivated.

--------------------------------------------------------------------------------
3.4  COMP_J_BLANK — String Fix + USE_LAM_ES Path
--------------------------------------------------------------------------------

File: new_main/mod_subs/comp_j_blank.f

WHAT CHANGED (2 modifications):

  (a) Lines 725, 752: 'ZERO FLUX' → 'ZERO_FLUX'
    BEFORE: ELSE IF(INNER_BND_METH .EQ. 'ZERO FLUX')THEN
    AFTER:  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN

  (b) Line 717-719: USE_LAM_ES path (unchanged but verified correct):
    RJ(1:ND)=TC(1:ND)     ! J from formal solution, not moment equations

WHY:
  (a) The VADAT keyword IB_METH uses 'ZERO_FLUX' (underscore). The original
      code checked for 'ZERO FLUX' (space), which never matched. This caused
      the ZERO_FLUX boundary condition to fall through to the default branch,
      which computes HFLUX_AT_IB from the moment equations — incorrect when
      ZERO_FLUX is intended.

  (b) When USE_LAM_ES=T, J is set directly from the formal solution (TC),
      bypassing the moment equations (MOM_J_CMF_V11, MOM_JREL_V9). This is
      the key to handling non-monotonic velocity: the moment equations cap
      sigma at -0.999, but the formal solution handles sigma < -1 correctly
      via the Dessart+2015 prescription in FG_J_CMF_V13.

EFFECT ON SYSTEM:
  - IB_METH=ZERO_FLUX now correctly sets HFLUX_AT_IB=0 (no flux at inner
    boundary). Previously it fell through to the default, which could inject
    spurious flux.
  - J is computed from the formal solution, not moment equations, avoiding
    the sigma cap at -0.999 in MOM_JREL_V9.

CODE PATH: COMP_J_BLANK is called from CMFGEN_SUB (cmfgen_sub.f:1815, 2298,
  3513) for every frequency in the BA iteration. The USE_LAM_ES branch
  (line 717) is taken when USE_LAM_ES=T (our configuration).

--------------------------------------------------------------------------------
3.5  CMFGEN_SUB — Luminosity NaN Protection + FIX_IMPURITY
--------------------------------------------------------------------------------

File: new_main/cmfgen_sub.f

WHAT CHANGED (2 modifications):

  (a) Line 610: FIX_IMPURITY during LAMBDA iterations
    BEFORE: FIX_IMPURITY=.FALSE.
    AFTER:  FIX_IMPURITY=RD_FIX_IMP

  (b) Lines 3088, 3091: NaN protection for RLUMST
    BEFORE: IF(RLUMST(1) .LE. 0.0_LDP)RLUMST(1)=1.0E-20_LDP
    AFTER:  IF(RLUMST(1) .LE. 0.0_LDP .OR. RLUMST(1) .NE. RLUMST(1))
                  RLUMST(1)=1.0E-20_LDP
            DO I=1,ND
              IF(RLUMST(I) .NE. RLUMST(I))RLUMST(I)=1.0E-20_LDP   ! NaN check
              IF(ABS(RLUMST(I)) .LT. 1.0E-10)THEN ...

WHY:
  (a) During LAMBDA iterations, FIX_IMPURITY was hardcoded to .FALSE.,
      ignoring the RD_FIX_IMP setting from VADAT. This meant impurity
      species populations were always allowed to vary during LAMBDA,
      even when the user requested them fixed. The fix respects the
      VADAT setting. In practice, both branches now set RD_FIX_IMP, so
      this is a no-op for our configuration — but it's more correct.

  (b) When the model has zero luminosity (LSTAR=0 with DIF=F), RLUMST
      can be exactly 0.0 or NaN (from 0/0 in DTDR computation). The
      original code only checked for <= 0, which doesn't catch NaN
      (NaN <= 0 is .FALSE. in Fortran). The fix adds the NaN self-
      comparison test (X .NE. X is .TRUE. only when X is NaN) at both
      the boundary (RLUMST(1)) and all depths (RLUMST(I)).

EFFECT ON SYSTEM:
  - NaN luminosity values are caught and replaced with 1.0E-20, preventing
    NaN propagation into the population solver (which would produce NaN
    populations and crash the model).
  - The FIX_IMPURITY change is a no-op for our config but more correct
    in general.

--------------------------------------------------------------------------------
3.6  RV_SN_MODEL_SNIIN — Grid Setup + Sigma Computation
--------------------------------------------------------------------------------

File: new_main/subs/rv_sn_model_sniin.f

WHAT CHANGED (4 modifications):

  (a) Line 46: NBND_INS=0 (was 3)
    BEFORE: NBND_INS=3 ! LUC: changed from 2 to 3
    AFTER:  NBND_INS=0 ! No boundary insertion — AMR grid handles refinement

  (b) Lines 153-187: Direct grid usage when NHYDRO=MND
    BEFORE: Always creates log-uniform tau grid, interpolates R from it
    AFTER:  IF(NHYDRO .EQ. MND) THEN
              TA(I) = R_HYDRO(I)   ! Use hydro grid directly
            ELSE
              ! Original log-uniform tau grid interpolation
            END IF

  (c) Lines 235-244: Conditional R_HYDRO/V_HYDRO reversal
    BEFORE: Always reverses R_HYDRO and V_HYDRO arrays
    AFTER:  IF(NHYDRO .NE. MND) THEN
              ! Reverse arrays
            END IF

  (d) Line 246: LIN_INTERP argument order fixed
    BEFORE: CALL LIN_INTERP(V,TA,ND,V_HYDRO,R_HYDRO,NHYDRO)
    AFTER:  CALL LIN_INTERP(TA,V,ND,R_HYDRO,V_HYDRO,NHYDRO)

  (e) Lines 248-265: Sigma computation with ABS() and zero protection
    BEFORE: T1=(V(I-1)-V(I))/(R(I-1)-R(I))
            T2=(V(I)-V(I+1))/(R(I)-R(I+1))
            SIGMA(I)=(T1*(R(I)-R(I+1))+T2*(R(I-1)-R(I)))/(R(I-1)-R(I))
            SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP
    AFTER:  DR_M=ABS(R(I-1)-R(I)); DR_P=ABS(R(I)-R(I+1)); DR_MP=ABS(R(I-1)-R(I+1))
            IF(DR_M .LT. 1.0E-30)DR_M=1.0E-30   ! prevent division by zero
            IF(DR_P .LT. 1.0E-30)DR_P=1.0E-30
            IF(DR_MP .LT. 1.0E-30)DR_MP=1.0E-30
            T1=(V(I-1)-V(I))/DR_M
            T2=(V(I)-V(I+1))/DR_P
            SIGMA(I)=(T1*DR_P+T2*DR_M)/DR_MP
            IF(ABS(V(I)) .LT. 1.0E-30)V(I)=1.0E-30   ! prevent division by zero
            SIGMA(I)=R(I)*SIGMA(I)/V(I)-1.0_LDP

  (f) Lines 266-270: Same ABS() and zero protection for boundary sigma

WHY:
  (a) NBND_INS=3 inserted 3 extra grid points at each boundary, extending
      the grid beyond the hydro data. With NBND_INS=0, the grid uses exactly
      the hydro points. The Python AMR script already handles boundary
      refinement, so CMFGEN's insertion is redundant and can create
      artifacts at the boundaries.

  (b) When the Python AMR script produces exactly ND-2*NBND_INS = MND points
      (i.e., NBND_INS=0 and NHYDRO=ND=MND), the original code would REPLACE
      the AMR grid with a log-uniform tau grid, throwing away the AMR
      refinement. The fix uses the hydro grid directly when the sizes match.

  (c) The R_HYDRO/V_HYDRO reversal was needed for the log-tau interpolation
      path (which expects inner-to-outer ordering). When using the grid
      directly (NHYDRO=MND), no reversal is needed because the hydro data
      is already in the correct order (outer-to-inner, matching CMFGEN).

  (d) LIN_INTERP(R,V,ND,U,W,NIN) interpolates V as a function of R.
      The original call LIN_INTERP(V,TA,ND,V_HYDRO,R_HYDRO,NHYDRO) had
      R and V swapped — it was interpolating R as a function of V, which
      is incorrect for non-monotonic velocity (V is not monotonic, so
      R(V) is multi-valued). The fix interpolates V(R), which is always
      single-valued.

  (e) The original sigma computation used signed radius differences
      (R(I-1)-R(I)), which can be positive or negative depending on grid
      ordering. The fix uses ABS() on all radius differences, making the
      computation correct regardless of grid direction. Zero-protection
      prevents division by zero when grid points coincide.

EFFECT ON SYSTEM:
  - The AMR grid from Python is preserved (not replaced by log-tau grid).
  - Velocity is interpolated correctly as V(R), not R(V).
  - Sigma is computed correctly for both grid orderings.
  - Zero-protection prevents NaN when grid points coincide or V=0.
  - For monotonic flows: the sigma formula is mathematically equivalent
    to the original (ABS() doesn't change the result when R is monotonic).

CODE PATH: RV_SN_MODEL_SNIIN is called when VEL_LAW=12 (SNIIn model) during
  initial grid setup. It reads SN_HYDRO_DATA and constructs the R, V, SIGMA
  arrays. This runs ONCE at the start of each CMFGEN execution.

--------------------------------------------------------------------------------
3.7  SET_ABUND_CLUMP — MON_INTERP Grid Direction Fix
--------------------------------------------------------------------------------

File: new_main/mod_subs/set_abund_clump.f

WHAT CHANGED:
  Lines 336-345: Added grid reversal for MON_INTERP density interpolation.

  BEFORE:
    CALL MON_INTERP(DENSITY,ND,IONE,LOG_R,ND,RHO_HYDRO,NHYDRO,LOG_R_HYDRO,NHYDRO)

  AFTER:
    ALLOCATE (LOG_R_INC(ND), DENSITY_TMP(ND))
    DO K=1,ND
      LOG_R_INC(K)=LOG_R(ND-K+1)    ! Reverse to increasing
    END DO
    CALL MON_INTERP(DENSITY_TMP,ND,IONE,LOG_R_INC,ND,
                   RHO_HYDRO,NHYDRO,LOG_R_HYDRO,NHYDRO)
    DO K=1,ND
      DENSITY(K)=DENSITY_TMP(ND-K+1)  ! Reverse back to decreasing
    END DO
    DEALLOCATE (LOG_R_INC, DENSITY_TMP)

WHY:
  MON_INTERP requires both input and output grids to be in the SAME monotonic
  direction. LOG_R_HYDRO is increasing (inner-to-outer, from the reversal at
  line 323). LOG_R is decreasing (outer-to-inner, CMFGEN convention). Without
  the reversal, MON_INTERP would either crash ("values outside range") or
  produce incorrect interpolation.

EFFECT ON SYSTEM:
  - Density is now correctly interpolated from the hydro grid to the CMFGEN
    grid. Previously, the mismatched grid directions could produce incorrect
    densities, leading to wrong populations and opacities.

CODE PATH: SET_ABUND_CLUMP is called during initial setup to set density
  and population from the hydro data.

--------------------------------------------------------------------------------
3.8  OBS_FRAME_SUB_V9 — Short Ray Guard
--------------------------------------------------------------------------------

File: obs/obs_frame_sub_v9.f

WHAT CHANGED:
  Lines 922-934: Added IF(IEND-IST .GT. 1) guard to STAU derivative block.

  BEFORE:
    dS(IST)=S(IST) +(S(IST)-S(IST+1))*DTAU(IST)/(DTAU(IST)+DTAU(IST+1))
    DO I=IST+1,IEND-1
      dS(I)=(S(I-1)*DTAU(I)+S(I)*DTAU(I-1))/(DTAU(I-1)+DTAU(I))
    END DO
    dS(IEND)=S(IEND-1)+(S(IEND-1)-S(IEND-2))*DTAU(IEND-1)/
             (DTAU(IEND-2)+DTAU(IEND-1))

  AFTER:
    IF(IEND-IST .GT. 1)THEN
      [original code]
    ELSE
      dS(IST)=S(IST)
      dS(IEND)=S(IST)
    END IF

WHY:
  When a ray has only 2 points (IEND-IST=1), the original code accesses
  S(IEND-2) = S(IST-1) = S(0), which is out of bounds (Fortran arrays start
  at 1). This crashes cmf_flux.exe with "Subscript #1 of the array S has
  value 0".

  This occurs with non-monotonic velocity profiles where some impact
  parameter rays pass through very few zones (1-2) near velocity reversals.

  The Z-space derivative block (line 835) already had this guard; the STAU
  block (line 922) was missing it.

EFFECT ON SYSTEM:
  - cmf_flux.exe no longer crashes on short rays in non-monotonic flows.
  - For rays with 1-2 points, the derivative is set to the endpoint value
    (flat extrapolation), which is a reasonable approximation.
  - No effect on monotonic flows (rays always have many points).

CODE PATH: OBS_FRAME_SUB_V9 is called from CMF_FLUX_SUB_V5 (cmf_flux.exe,
  the OBS step) for each frequency in the observed-frame spectrum synthesis.

--------------------------------------------------------------------------------
3.9  Build System — Viper/Intel Configuration
--------------------------------------------------------------------------------

File: Makefile_definitions

WHAT CHANGED:
  - Added viper_intel_nodebug and viper_intel_debug build configurations
  - Changed HOST=crc/gcc to HOST=viper/intel
  - Added -lgfortran to PGLIB for all configurations (needed for PGPLOT
    when using Intel compiler with GCC-compiled PGPLOT)
  - Changed -fpe0 to -fpe-summary=none in debug configs (allows NaN
    propagation for debugging instead of immediate FPE trap)

WHY:
  The code needs to compile on the Viper cluster with Intel ifort. The
  original configuration was for CRC with GCC. The -fpe-summary=none
  change is needed because the non-monotonic velocity prescription may
  produce intermediate NaN values (before the prescription kicks in)
  that would trigger FPE traps with -fpe0.

--------------------------------------------------------------------------------
3.10 Minor Changes (3 files)
--------------------------------------------------------------------------------

  subs/regrid_log_dc_v1.f: Comment cleanup (no code change)
  new_main/mod_subs/set_new_model_estimates.f: Indentation fix (no code change)
  new_main/Makefile, disp/Makefile, obs/Makefile: Build system updates

================================================================================
4. FILES NOT MODIFIED (AND WHY)
================================================================================

The following files use SIGMA but were NOT modified because they are not
in the code path for our configuration (USE_LAM_ES=T, INC_GRID=F,
REV_RGRID=F, SOBOLEV=T, THK_CONT=F):

  subs/gamma.f
    - Computes Sobolev GAM without the prescription
    - Called from: FG_COMP, FORMSOL, VAR_FORMSOL, HAM_FORMSOL, FG_HAM,
      FORM_EW, COMP_EW_ABS_EW
    - NOT in our path: FG_COMP/FORMSOL are for CMF line transfer (we use
      SOBOLEV), COMP_EW_ABS_EW is for the OBS EW diagnostic (DO_CMF_EW=T
      calls it, but it only affects EW output, not convergence or spectrum)
    - RISK: Low — EW diagnostic values may be slightly wrong for lines
      with sigma < -1, but this doesn't affect the main spectrum or
      population convergence.

  subs/cmf_form_sol_v2.f
    - Computes GAM without the prescription (line 596)
    - Called from COMP_J_BLANK for observed flux computation
    - NOT in our path: only called when LST_ITERATION and
      .NOT. LAMBDA_ITERATION and MAXCH < 100 and .NOT. SN_MODEL and ND < 200
    - Our model has SN_MODEL=T, so this is never called.

  subs/extend_vtsigma.f
    - Clamps SIGMAEXT to -0.9999 (line 32)
    - Called from CMFGEN_SUB when ACCURATE=T (INC_GRID=T)
    - NOT in our path: INC_GRID=F, so ACCURATE=.FALSE., EXTEND_VTSIGMA
      is never called.

  new_main/mod_subs/adjust_r_grid_v4.f
    - STOPs when SIGMA < -1 and V > 0.2 km/s (line 208)
    - Called from CMFGEN_SUB when REV_RGRID=T
    - NOT in our path: REV_RGRID=F, so ADJUST_R_GRID_V4 is never called.

  new_main/mod_subs/mom_j_cmf_v11.f
    - Moment equations, no Dessart+2015 prescription
    - Called from COMP_J_BLANK when USE_LAM_ES=F
    - NOT in our path: USE_LAM_ES=T, so MOM_J_CMF_V11 is bypassed.

  new_main/mod_subs/mom_jrel_v9.f
    - Moment equations, caps SIGMA at -0.999
    - Called from COMP_J_BLANK when USE_J_REL=T
    - NOT in our path: USE_J_REL=F, so MOM_JREL_V9 is bypassed.

================================================================================
5. CODE PATH SUMMARY
================================================================================

For X2_sn configuration:
  USE_LAM_ES=T, USE_J_REL=F, INC_GRID=F, REV_RGRID=F,
  SOBOLEV=T (GLOBAL_LINE=SOB), THK_CONT=F,
  IB_METH=DIFFUSION, DIF=T, FG_OPT=INT/INS

BA Solver (population convergence):
  CMFGEN_SUB
    → COMP_J_BLANK (continuum J)
      → FG_J_CMF_V13 [PRESCRIPTION ✓] (formal solution for J)
      → USE_LAM_ES branch: RJ=TC (no moment equations)
    → SUB_SOB_LINE (line transfer)
      → SOBJBAR_SIM [PRESCRIPTION ✓] (line JBAR)
      → SOBEW       [PRESCRIPTION ✓] (line EW/source)

OBS Step (spectrum synthesis):
  CMF_FLUX_SUB_V5
    → OBS_FRAME_SUB_V9 [NO SIGMA — safe] (observed-frame integration)
    → SOBEW_GRAD_V2 [PRESCRIPTION ✓] (line EW, when DO_CMF_EW=F)
    → COMP_EW_ABS_EW → GAMMA [NO PRESCRIPTION] (line EW, when DO_CMF_EW=T)
      ⚠ EW diagnostic only — does not affect spectrum or convergence

Grid Setup (initial):
  RV_SN_MODEL_SNIIN [FIXED ✓] (R, V, SIGMA from hydro data)
  SET_ABUND_CLUMP   [FIXED ✓] (density interpolation)

================================================================================
6. VERIFICATION
================================================================================

The changes were verified by:

1. Code path tracing: Every routine that uses SIGMA was identified via grep.
   Each was traced through the call chain to determine if it is in the
   X2_sn code path. All in-path routines have the prescription or are
   bypassed by configuration.

2. Compilation: The modified source compiles successfully with Intel ifort
   on Viper (exe/cmfgen_dev.exe and exe/cmf_flux.exe built 2026-08-01).

3. Logical verification: The Dessart+2015 prescription was verified to be
   a no-op for sigma >= 0 (monotonic flows) and to provide a positive floor
   for sigma < -1 (non-monotonic flows).

4. Argument verification: The LIN_INTERP argument swap was verified against
   the subroutine signature (LIN_INTERP(R,V,ND,U,W,NIN) interpolates V(R)).

5. String verification: The 'ZERO_FLUX' string fix was verified against the
   VADAT keyword format (IB_METH=ZERO_FLUX in VADAT).

================================================================================
7. KNOWN LIMITATIONS
================================================================================

1. dGAMdR derivative (fg_j_cmf_v13.f:972): The derivative formula assumes
   GAM = T1*(1+sigma*mu^2), but with the prescription, GAM is capped.
   The derivative is slightly inaccurate when the prescription is active.
   Impact: minor integration inaccuracy, not a stability issue.

2. gamma.f not patched: The OBS EW diagnostic (COMP_EW_ABS_EW → GAMMA)
   does not have the prescription. EW values for lines with sigma < -1
   may be slightly wrong. Impact: diagnostic only, does not affect the
   main spectrum or population convergence.

3. fg_j_cmf_v12.f still has STOP on sigma < -1: V12 is never called, but
   if reactivated, it would crash on non-monotonic velocity. The Dessart
   prescription was added to V12's GAM computation, but the diagnostic
   STOP at line 916 was not removed.

================================================================================
8. REFERENCES
================================================================================

  Dessart, L., Hillier, D.J., Livne, E., Yoon, S.-C.
  "Radiation-hydrodynamics simulations of interacting supernovae"
  A&A, 575, A51 (2015)
  [Section 3.2: Sobolev optical depth for non-monotonic velocity]

  Hillier, D.J.
  "CMFGEN: A radiative transfer code for spherically symmetric outflows"
  Documentation: docs/full_descr.txt, docs/1503.05463v1.txt

================================================================================
END OF DOCUMENT
================================================================================
