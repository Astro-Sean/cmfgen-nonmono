# CMF_FLUX Observer-Frame Pipeline Audit

Source-level audit of the observed-frame flux pipeline in `obs/`,
performed against HEAD `ce82a9f0` (branch `main`, clean except untracked
`drv_knop_ray.i90`, `unit_tests/systems/_pv.txt`). All findings below are
from reading the current source; where a claim rests on runtime behavior
that has not been measured, it is marked UNVERIFIED. Prior audit
(`CMF_FLUX_SPECTRAL_AUDIT.md`) covered the same files at HEAD
`8224a00b`; this document supersedes its pipeline description but keeps
its open validation items.

Reference set: Hillier & Miller (1998, ApJ 496, 407); Dessart, Audit &
Hillier (2015, arXiv:1503.05463); Knop et al. (2009, Paper I); Mihalas
(1978); Rybicki & Lightman (1979).

---

## 1. Execution path

```
cmf_flux_v5.f                       driver: reads RVTJ (R,V,SIGMA,ED,T,CLUMP_FAC,DENSITY),
                                    reads controls via rd_cmf_flux_controls.f (OBS keys)
        |
cmf_flux_sub_v5.f                   per-frequency ML loop, SECTION='CONTINUUM' (CONT_VEL=T)
        |
        +-- OPACITIES_V5.INC        unclumped continuum CHI, ETA at FL:
        |      ESOPAC(ESEC), RAYLEIGH_SCAT(CHI_RAY), COMPUTE_DUST_OPACITY
        |      (CHI_ABS_DUST, CHI_SCAT_DUST), GENOPAETA_V10 free-free/bf,
        |      TWO_PHOT_OPAC_V3, XOPAC_V4.  CHI_SCAT = ESEC+CHI_RAY+CHI_SCAT_DUST.
        |      ESEC = 6.65E-15 * ED   (driver, cmf_flux_sub_v5.f)
        |
        +-- COMP_JCONT_V4.INC       CMF J solve on CLUMPED opacity
        |      CHI_CLUMP=CHI*CLUMP_FAC; ETA_CLUMP=ETA*CLUMP_FAC;
        |      ESEC_CLUMP, CHI_RAY_CLUMP, CHI_SCAT_CLUMP, CHI_SCAT_DUST_CLUMP
        |      (COMP_JCONT_V4.INC:123-128)
        |      DBB boundary gradient; HFLUX_AT_IB = DBB/CHI(ND)/3 (:114)
        |      fork: DBB floored to 1E-99 under nonmono when DBB<0 (:109-113)
        |      f-iteration loop (max 10 ACCURATE / 20 non-ext):
        |        FG formal solve (TC = ray-by-ray J):
        |          DO_NONMONO_KNOP -> FG_J_CMF_NONMONO_V1 (SIGMA_TRUE, :816)
        |          else           -> FG_J_CMF_V12 (clamped SIGMA, :825)
        |          ACCURATE grid: V12 only; DO_NONMONO_KNOP -> STOP (:401-404)
        |        MOM solve (RJ):
        |          USE_LAM_ES=T      -> RJ = TC (formal J, :899)
        |          DO_NONMONO_SOB+LAM_ES=F -> MOM_J_CMF_NONMONO_V1 (:920)
        |          else             -> MOM_J_CMF_V8 (:929)
        |          (USE_J_REL -> MOM_JREL_V9; USE_DJDT_RTE -> MOM_J_DDT_V6;
        |           USE_FORMAL_REL -> CMF_FORMAL_REL_V4/V5; PP modes)
        |      RJ -> LU_EDD record; non-coherent ES convolved J -> RJ_ES
        |      via COMP_J_CONV_V2 between NUM_ES_ITERATIONS sweeps.
        |
        +-- store, per ML (unclumped convention):
        |      ETA_CMF_ST=ETA, CHI_CMF_ST=CHI, RJ_CMF_ST=RJ, RJ_ES_ST,
        |      IB_J_CMF_ST=RJ(ND), IB_H_CMF_ST=HFLUX_AT_IB
        |      RAY_ST=CHI_RAY, DUST_ST=CHI_SCAT_DUST (stored unclumped)
        |
        +-- COMP_OBS_V3 (CONT_VEL only, per ML, subs/comp_obs_v3.f)
        |      stores I+(p) CMF boundary intensities; Doppler-transforms
        |      (FREQ_CONV_FAC), integrates over p with FLUX_WGHTS=mu*dmu,
        |      fills OBS_FLUX -> written to OBSFLUX after the ES loop.
        |      FLUX_CONST = 2*PI*(1E18*RMAX/PARSEC)^2 (U) or 4*PI*.. (I+)
        |
        +-- phase-1 render (always):
        |      ETA_RND_ST = ETA_CMF_ST + RJ_ES*ESEC (+ RJ*RAY_ST if
        |      INCL_RAY_SCAT) (+ RJ*DUST_ST if INCL_DUST)     [unclumped]
        |      CHI_RND_ST = CHI_CMF_ST                        [unclumped]
        |      DO_CLUMP_MODEL -> ETA_RND_ST*=CLUMP_FAC,
        |                          CHI_RND_ST*=CLUMP_FAC      (once)
        |      OBS_FRAME_SUB_V9(...)  -> OBS_FLUX -> OBSFRAME (phase 1)
        |
        +-- DO_SOBOLEV_LINES block: SUB_SOB_LINE_V3 per blend group
        |      -> CHIL_MAT/ETAL_MAT (Sobolev line opacity/emissivity on
        |      SIGMA_TRUE under DO_NONMONO_SOB), CHECK_LINE_OPAC /
        |      NEG_OPACITY maser census, EWDATA, MEANOPAC
        |
        +-- phase 2 (DO_SOB_FLUX_OBS=T only):
               ADD_SOB_LINE_PROF deposits each 'SOB' line's Doppler core
               into ETA_CMF_ST/CHI_CMF_ST UNCLUMPED
               (production routine, cmf_flux_sub_v5.f:3157-3379);
               ES rescatter pass; OBS_FRAME_SUB_V9 re-render ->
               OBSFRAME overwritten (continuum+BLANK+SOB lines).
               OBSFLUX keeps the phase-1 continuum+BLANK product.
```

`obs_frame_sub_v9.f` itself, per p-ray `LS`:
build `P`, `Z_RAY`, `R_RAY`, `VMU_RAY` (= z*v/(r*c)), `GAM_RAY`,
`NU_ON_NU0_RAY` = nu_obs/nu_cmf = 1/[(1-VMU)*GAM]; insert points on
projected-velocity and `OBS_ES_DTAU` criteria (`N_INS_OBS`); per observer
frequency ML inside `!$OMP PARALLEL ... !$OMP DO` (:785-818): hunt the
CMF bracket `K` on the band `SM_FREQ_CMF` (STRT_INDX_CMF..END_INDX_CMF),
interpolate ETA/CHI in the stored arrays, apply observer-frame factors,
integrate along the ray (`OBS_INT_METHOD` = STAU or ETAZ), then
`OBS_FLUX(ML) += PAR_FLUX * HQW_AT_RMAX(LS)`. After all rays:
`OBS_FLUX = OBS_FLUX * 6.59934 * R(1)**2` -> Jy at 1 kpc.

## 2. Answers to the audit questions

1. **CMF radiation field solver.** `COMP_JCONT_V4.INC` dispatches per
   frequency: formal ray J `TC` from `FG_J_CMF_V12` (default), or
   `FG_J_CMF_NONMONO_V1` under `DO_NONMONO_KNOP` (with `SET_FG_NU_NEXT`,
   `SIGMA_TRUE`); moment J `RJ` from `MOM_J_CMF_V8` (default),
   `MOM_J_CMF_NONMONO_V1` under `DO_NONMONO_SOB .AND. .NOT.USE_LAM_ES`,
   or `RJ=TC` under `USE_LAM_ES=T`. All solves run on clumped
   `CHI_CLUMP/ETA_CLUMP/ESEC_CLUMP` (`COMP_JCONT_V4.INC:123-128`).
   The `DO_NONMONO_SOB` default render thus uses the `FG_J_CMF_V12`
   formal J on clamped SIGMA (the 10-Nov-2026 decision recorded in the
   file header).

2. **Opacities/emissivities into the observer-frame solve.** Render
   inputs are `ETA_RND_ST`, `CHI_RND_ST`: stored unclumped ETA plus
   coherent-scattering emissivity `RJ*CHI_SCAT`-equivalent terms
   (`RJ_ES_ST*ESEC` non-coherent; `RJ_CMF_ST*RAY_ST`, `*DUST_ST`),
   then a single `*CLUMP_FAC` fold under `DO_CLUMP_MODEL`. Phase 2 adds
   unclumped Sobolev line deposits to `ETA_CMF_ST/CHI_CMF_ST` before the
   fold.

3. **Scattering.** Electron scattering opacity `ESEC=6.65E-15*ED`
   (unclumped ED from RVTJ; `CLUMP_FAC` multiplies it at the same fold
   as every other opacity, the standard CMFGEN `chi_eff=f*chi`
   convention, population-side included). Non-coherent ES: `RJ_ES`
   from `COMP_J_CONV_V2` between `NUM_ES_ITERATIONS` sweeps.
   Rayleigh `CHI_RAY` coherent; dust scattering `CHI_SCAT_DUST`
   isotropic unless `USE_HEN_GREEN` (then only `CMF_FORMAL_REL_V5`).
   Blanketed lines enter through `ETA/CHI` of `OPACITIES_V5`/line
   blanketing; Sobolev lines only through the phase-2 deposit.

4. **Impact-parameter integration.** `PAR_FLUX` per ray weighted by
   `HQW_AT_RMAX` = `mu*dmu` quadrature at `R(1)` (Gauss-Legendre
   `GAULEG` for plane-parallel, `NPEXT` ray set otherwise);
   `OBS_FLUX = sum_LS PAR_FLUX*HQW * 6.59934 * R(1)^2`. With R in
   1e10 cm and D=1 kpc, `6.59934 ~ 2*pi*(1e10/3.086e21)^2/1e-23`
   -> Jy. Constant verified dimensionally; the exact digits depend on
   the PARSEC() constant used at build time (UNVERIFIED to the 4th
   decimal; check `FUN_PI`/`PARSEC` source before quoting).

5. **Flux normalization / Jy.** See 4. `COMP_OBS_V3` separately writes
   `OBSFLUX` via `FLUX_CONST*(I+ or U)` (`subs/comp_obs_v3.f:129-131`).

6. **Observable spectrum files.** `OBSFRAME` = `OBS_FRAME_SUB_V9`
   formal ray integral (phase 1 continuum+BLANK; phase 2 adds SOB
   lines when `DO_SOB_FLUX_OBS=T`). `OBSFLUX` = `COMP_OBS_V3`
   Doppler-transform product (continuum+BLANK only) - a *different
   numerical method*, not a component of OBSFRAME. Diagnostics:
   `IP_DATA` (I(p)), `RTAU_DATA`, `ZTAU_DATA`, `dFR_DATA`,
   `FLUX_FILE`, `EWDATA`, `TRANS_INFO`, `MEANOPAC`, `J_COMP`,
   `EDDFACTOR`, `ES_J_CONV`, `PLANCK_KAPPA_MEAN`, `PHOTOSPHERIC_RADIUS`,
   `CMF_FORCE_DATA`/`SOB_FORCE_DATA`, `OUTGEN`.

7. **Unused / uninitialized / overwritten quantities.**
   - `HFLUX_AT_IB` (= `IB_H_CMF_ST`) is the analytic estimate
     `DBB/(3*CHI(ND))` (`COMP_JCONT_V4.INC:114`), not the solver's H:
     `FG_J_CMF_V12`, `FG_J_CMF_NONMONO_V1`, `MOM_J_CMF_V8`,
     `MOM_J_CMF_NONMONO_V1` never return an H at the inner boundary;
     only `MOM_JREL_V9`/`MOM_J_DDT_V6` update it. `IB_J_CMF_ST=RJ(ND)`
     is the solved J. The injected pair (J solved, H analytic) is only
     consistent in the diffusion limit.
   - `CLUMP_FAC` argument of production `ADD_SOB_LINE_PROF`
     (`cmf_flux_sub_v5.f:3175`) is declared `INTENT(IN)` and never
     referenced; `DROP_SUM`/`INJ_SUM` are documented by caller comments
     as "clumped emissivity" but accumulate unclumped `ETAL_MAT`.
   - `sob_line_es_update.f` is compiled (`obs/Makefile:34`) but called
     only by `unit_tests/f90/drv_sob_line_es.f` - dead in production.
     It assumes `ETA_ST` already contains clumped scattering terms
     (`:86-94`), a different storage convention than the production
     arrays it mimics.
   - `COMP_OBS_V2` exists in `subs/` and is called by `cmfgen_sub.f`
     and `lte_sub.f`, not by `cmf_flux_sub_v5` (which calls V3).
   - `OBS_FLUX` is written by `COMP_OBS_V3` during the ML loop, zeroed
     and refilled by `OBS_FRAME_SUB_V9` afterwards; the two uses are
     temporally disjoint (OBSFLUX is written before the render), so no
     stale-data leak, but the shared array is a readability hazard.

## 3. Confirmed defects (source level)

**D1. ETAZ inner-boundary asymmetry (Pillar 3). FIXED**
`OBS_FRAME_SUB_V9` with `DO_IB_DIFF_FLUX` injects the diffusion
boundary intensity `I_nu(mu) = J_nu + 3*mu*H_nu` (relativistically
transformed `*NU_ON_NU0^3`) only in the `STAU` integration branch
(`obs_frame_sub_v9.f:~1191-1210`): when the ray foot reaches the inner
boundary (`SM_NRAY == NRAY`) `PAR_FLUX` is initialized from
`IB_J_CMF/IB_H_CMF` interpolated at the local comoving frequency. The
`ETAZ` branch (`:~1243-1263`) initialized `PAR_FLUX=0` for the same
geometry - the boundary-originating intensity was silently dropped, so
`INT_METH=ETAZ` under `DIF` lost the inner-boundary flux on every ray
that intersects the core. Fix: the ETAZ branch now seeds
`PAR_FLUX = I_IB * exp(-TAU(IEND))` for rays with
`DO_IB_DIFF_FLUX .AND. SM_NRAY==NRAY .AND. LS<=NC`, i.e. the same
interpolated J+3muH seed attenuated by the whole-ray optical depth
built from the same corrected DTAU. Manufactured tests show thin-model
STAU/ETAZ agreement at 5e-9 and convergence of the thick-model
difference with ray resolution (9.6% -> 3.9% at 4x insertions).

**D2. ADD_SOB_LINE_PROF production vs .INC divergence (Pillar 4/5). FIXED**
Two implementations existed with different clumping conventions and
different signatures:
- Production (`cmf_flux_sub_v5.f:3157-3379`): deposits
  `CHIL_MAT(I,S)*WT`, `ETAL_MAT(I,S)*WT` - UNCLUMPED, matching the
  stored-array convention (render fold applies CLUMP_FAC once).
  Signature includes `ARAD_VEC`, `C4_VEC`, `ED_IN` for
  `DO_SOB_DEP_WINGS`/`DO_SOB_DEP_FLOOR`; has the per-line >10% maser
  drop report.
- `.INC` copy (`obs/ADD_SOB_LINE_PROF.INC:127-138`): deposited
  `CHIL_MAT*CLUMP_FAC`, `ETAL_MAT*CLUMP_FAC` - CLUMPED at deposit,
  which would have double-clumped under the render fold; lacked the
  wings/floor arguments and the per-line report.
The in-file comment said "keep the two in sync"; they were not.
`test_sob_flux_obs` even asserted the clumped-deposit expectation.
FIX: `obs/ADD_SOB_LINE_PROF.INC` now holds the production routine
verbatim and `cmf_flux_sub_v5.f` INCLUDEs it - a single copy compiled
by production and tests. `test_sob_flux_obs` now asserts the
unclumped convention (deposit independent of CLUMP_FAC, integral =
CHIL*1e-15) plus static pins that the .INC has no CLUMP_FAC factor
and the render fold exists exactly once.

**D3. Blue-edge guard asymmetry at loc 1 (Pillar 6, latent). FIXED**
The loc-1 frequency hunt (`obs_frame_sub_v9.f:824-866`) only hunted K
upwards (red direction). If `T1 >= SM_FREQ_CMF(1)` the loop was
skipped, `K=1`, and `T2=(T1-SM_FREQ_CMF(1))/(SM_FREQ_CMF(0)-...)`
read `SM_FREQ_CMF(0)` - out of bounds. The I>1 path has the
symmetric blue-edge clamp (:899-924); loc 1 did not. In practice
`VMU_RAY(1) = z*v(1)/(r(1)*c) > 0` on every ray (front surface,
outward flow) so `T1 < OBS_FREQ <= SM_FREQ_CMF(1)` and the hole did
not fire on supernova models; it was unguarded, not unreachable. Fix:
loc 1 now mirrors the I>1 clamp (K=2, T1 clamped to the band edge,
census increment, abort beyond the 2% margin). Side effect: every
stored `INDX` is now in [2,SM_NCF], so the boundary lookup
`K=INDX(SM_NRAY)` cannot read element 0 (closes the P2.4 hole).

## 4. Potential defects / numerical risks

**P1. OpenMP `INDX` semantics (Pillar 2).** `INDX` is `PRIVATE` in the
`!$OMP PARALLEL` region (`:785-787`); the serial fill `INDX(1:NRAY)=1`
(:778) does not propagate to threads - ifort gives each thread a
private copy with inherited bounds but *undefined contents* (the
:789-799 comment documents the verified allocation behavior).
Per-iteration analysis: `INDX(1)` is written (:826, :865) before any
read; for I>1, `K=INDX(I-1)` (:895) reads the value written at
`INDX(I)=K` (:980) earlier in the same ML iteration. Write-before-read
appears to hold for all I within an iteration *provided the I loop
always starts at 2 and executes serially* - which it does. Residual
risks: (a) if a code path ever reads `INDX(I)` for I not yet written
(e.g. `K=INDX(SM_NRAY)` at :1200 - written at :980 during the I loop,
safe unless the loop early-exited before I=SM_NRAY; SM_NRAY can be
reduced by TAU_MAX truncation - CHECK whether :1200 is reached only
when SM_NRAY==NRAY full count, i.e. boundary-intersecting rays, in
which case INDX(SM_NRAY)=INDX(NRAY) was written); (b) no explicit
per-thread reinitialization of the `K=1` start (line :825 hardcodes
`K=1`, fine). Verdict: probably correct by write-before-read, but an
explicit `INDX(1:NRAY)=1` per-thread after :804 (or FIRSTPRIVATE) is
cheap insurance and removes reliance on control-flow analysis.
VERIFIED AND FIXED: per-thread `INDX=1` init added (commit
`d1e09a31`); a 16-frequency driver case gives bitwise-identical output
at OMP_NUM_THREADS=1/2/4/8 (x2 reps), and a bounds/uninit/fpe0 build
runs clean.

**P2. `K=INDX(SM_NRAY)` validity at the inner boundary (:1200).**
Read of `SM_FREQ_CMF(K-1)` with `K=INDX(SM_NRAY)`; K=1 would read
element 0 (same class as D3). RESOLVED: with the loc-1 blue-edge
clamp every stored INDX lies in [2,SM_NCF], so K>=2 and
`J-1 = STRT_INDX_CMF+K-2 >= STRT_INDX_CMF >= 1`. No element-0 read
reachable.

**P3. `IB_H_CMF_ST` provenance.** As above, H is analytic
`DBB/(3*chi(ND))` while J is solved. Under DIF this is the diffusion
estimate; when the J solve returns a boundary J inconsistent with the
local Planck gradient (line scattering dominating at the photosphere),
`J+3*mu*H` mixes a converged J with a first-order H. Physics choice,
upstream-consistent in the DIF limit; flag for manufactured test D.

**P4. Edge clamps substitute boundary opacity.** Three clamp sites
(:841-855 red loc1, :909-924 blue, :955-970 red) use boundary opacity
within a 2% window and abort beyond. Census exists
(`N_EDGE_CLAMP`, `CLAMP_FRAC_MAX`, worst indices) but there is no
strict production gate requiring zero unexplained clamps, and the
clamp means a genuinely under-covered Doppler reach degrades the
spectrum silently up to 2%. Sobolev profile truncation at band edges
is not separately tracked.

**P5. `sob_line_es_update.f` convention.** Dead in production (unit
test only); assumes clumped scattering inside `ETA_ST` and computes
`IB_H_ST=DBB*CLUMP_FAC(ND)/(3*CHI_ST(ND,ML))` (:203) - consistent only
if `CHI_ST` is the unclumped stored opacity (then it reduces to the
COMP_JCONT value times CLUMP_FAC(ND); i.e. it *differs* from the
production `IB_H_CMF_ST` by exactly `CLUMP_FAC(ND)` under the same
inputs). If ever wired into production it changes the boundary H.
Document or align before use.

**P6. Two different "observed spectrum" methods coexist.** `OBSFLUX`
(COMP_OBS_V3 Doppler transform of CMF boundary intensities) vs
`OBSFRAME` (direct observer-frame formal integral). Their agreement is
itself a consistency check that is not currently exercised in tests.

**P7. OpenMP `K` bounds in hunts rely on 10-cell window margin.** The
`STRT_INDX_CMF/END_INDX_CMF` window (built from
`MAX_VMU=MAXVAL(ABS(VMU_RAY))` plus margin) is what makes `K-1>=1` and
`K+1<=SM_NCF` true in practice; there is no assert. A margin shrink
turns P2/D3 into live reads.

**P8. Relativistic powers.** Verified by inspection: emissivity
`*NU_ON_NU0^2`, opacity `/NU_ON_NU0`, boundary intensity
`*NU_ON_NU0^3` with `NU_ON_NU0 = nu_obs/nu_cmf` and
`nu_cmf = nu_obs*(1-VMU)*GAM` - the standard `eta/nu^2`, `chi*nu`,
`I/nu^3` invariants (Mihalas 1978; R&L 1979). Angular `mu` used in
`J+3*mu*H` is the ray-frame z/r at the foot - the CMF-to-observer mu
aberration is *not* applied to the boundary moment direction (H is a
CMF moment injected with observer mu; at v/c<<1 the error is O(v/c);
DO_RELO runs should use the transformed mu - UNVERIFIED whether this
matters at SN velocities, ~1% level).
`GAM_RAY` is built from `VMU_RAY` (the projected speed?) - check
whether gamma uses total or projected v (line ~570-700 region;
UNVERIFIED - if `GAM=1/sqrt(1-VMU^2)` uses the projected component the
Lorentz factor is wrong for non-radial ray segments; check
`VMU_RAY`/`GAM_RAY` construction lines 560-700 before finalizing).
No superluminal input rejection observed (radicand floor? check).
Zero-velocity limit reduces to `NU_ON_NU0=1`, `GAM=1` - trivially
consistent.

**P9. Ray insertion criterion.** Insertion keys on the *preceding*
cell's projected velocity change (`MAX_DEL_V_RES_ZONE`,
`OBS_ES_DTAU`); on velocity reversals the Doppler excursion across a
cell is not bounded by the upwind cell width alone - flagged by prior
audit as a coverage risk for non-monotonic fields; resolution
convergence (Pillar 11) untested.

**P10. Two-phase contamination.** Phase-2 deposits write into
`ETA_CMF_ST/CHI_CMF_ST` in place; the arrays are not restored, so a
post-phase-2 re-render would double-deposit if ever invoked. `OBSFLUX`
is written before phase 2, so phase-2 cannot contaminate it - verified
by ordering (:2379 write precedes phase-2 render). UNVERIFIED: no
guard prevents a second `DO_SOB_FLUX_OBS` pass.

**P11. TRANS_TYPE values matching no spectral channel.** Every line
lands in exactly one of: 'BLANK' -> inserted into the CMF grid by
`INS_LINE_V5`/`INS_LINE_OBS_V4` and blanketed; 'SOB' -> excluded
from that grid, reaching the observed spectrum only via the
`DO_SOB_FLUX_OBS` deposit (flag off = EW estimate only, upstream
semantics). Any other `TRANS_TYPE` - e.g. 'CMF' from `TRANS_*='CMF'`
or `GLOBAL_LINE='CMF'`, an accepted VADAT value - matches neither
channel and was previously invisible in the observed spectrum with
no diagnostic (a longstanding upstream quirk, verified in
`cmfgen_reference/cur_cmf` - every inserter filters strictly on
'BLANK'). ADDRESSED by census: `cmf_flux_sub_v5.f` now counts
BLANK/SOB/other after `VEC_TRANS_TYPE` is finalized, names every
unmatched transition, and warns when SOB lines exist with
`DO_SOB_FLUX_OBS=F`. Diagnostic only (LUER lines); the upstream
quirk itself is preserved - aborting on 'CMF' was rejected since
CMF is a legitimate population-solve transfer choice whose
observed-frame meaning is simply "not blanketed here".

## 5. Test coverage status

- `unit_tests/` drivers exist for Sobolev machinery; `drv_sob_line_es.f`
  exercises `SOB_LINE_ES_UPDATE` (orphan routine, P5). The Sobolev
  deposit tests now exercise the production routine via the shared
  `.INC` (D2 fixed).
- Manufactured tests now exist for: STAU/ETAZ boundary agreement and
  boundary moment scaling (drv_obs_frame_ib, 10 checks), OpenMP
  reproducibility at 1/2/4/8 threads, clump-independence of the
  deposit (CLUMP_FAC in {0.5,1,1.5,2}), nonuniform-NU deposit
  integrals, maser/blend-repair drops, off-grid and non-SOB skips,
  negative-CHI STAU->ETAZ fallback, reversing-velocity ETAZ.
- Still no manufactured tests for: angular integration vs
  `pi*I*(R/D)^2`, frequency/ray-resolution convergence sweeps,
  superluminal rejection, `DO_SOB_DEP_WINGS` Voigt path.
- Prior audit open items still apply: no live numerical comparison of
  OBSFRAME vs OBSFLUX, no non-monotonic runtime validation of the
  render path, KNOP maturity gate absent (`GET_FG_GS_MATURITY` exists
  solver-side but no OBS-side acceptance gate).

## 6. Limitations of this audit

Source-level only. No numerical reproduction yet. `VMU_RAY`/`GAM_RAY`
construction (lines ~560-700) and the exact STAU/ETAZ quadrature
weights were inspected; equations above are quoted from code, not
re-derived on paper. `subs/comp_obs_v2.f` inspected by header only.
Population-side clumping consistency (whether `comp_j_blank.f` applies
`ESEC_CLUMP` identically) assumed from the CMFGEN-wide convention -
spot-verify before closing Pillar 4.
