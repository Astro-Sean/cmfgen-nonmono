# P8.1 - Electron-scattering redistribution audit

Source-level audit of the electron-scattering (e.s.) contribution to
the CMF flux pipeline and the observer-frame formal integral. No
production source was modified. All references are to this tree.

Scope: `obs/cmf_flux_sub_v5.f`, `obs/COMP_JCONT_V4.INC`,
`obs/obs_frame_sub_v9.f`, `new_main/mod_subs/comp_j_blank.f`,
`new_main/mod_subs/do_var_cont.f`, `new_main/cmfgen_sub.f`,
`subs/comp_j_conv_v2.f`, `new_main/mod_subs/comp_opac.f`,
`obs/OPACITIES_V4/V5.INC`, `obs/rd_cmf_flux_controls.f`,
`new_main/mod_subs/radiation_mod.f`.

## 1. Definitions and units

| Symbol | File:line | Physical meaning | Units |
|---|---|---|---|
| `ESEC` | `subs/esopac.f:14` | electron-scattering opacity `sigma_es = sigma_T * ED` | opacity (per length, code units: `6.65E-15*ED`, ED in 1e10 cm^-3 -> sigma_T=6.65e-25) |
| `RJ` | `radiation_mod` | comoving mean intensity J_nu at the current frequency | J units |
| `RJ_ES` | `radiation_mod.f:11` | "Convolution of RJ with e.s. R(v'v')" - the redistribution-convolved mean intensity `Jbar_es(nu) = int R(nu',nu) J(nu') dnu'` | J units (NOT emissivity) |
| `CHI_SCAT` | `comp_opac.f:181`, `OPACITIES_V4/5.INC` | total scattering opacity: `ESEC + CHI_RAY (+ CHI_SCAT_DUST obs side)` | opacity |
| `CHI_RAY` | `comp_opac.f:175-179` | Rayleigh scattering opacity; `=0` unless `INC_RAY=T` (`rd_control_variables.f:2921`) | opacity |
| `CHI_SCAT_DUST` | `OPACITIES_V5.INC:47-49` | dust scattering opacity; `=0` unless `INCL_DUST` | opacity |
| `NUM_ES_ITERATIONS` | `cmf_flux_cntrl_var_mod.f:151`, read `rd_cmf_flux_controls.f:105` (`NUM_ES`) | OBS-side e.s. redistribution iteration count | count |
| `RD_COHERENT_ES`/`COHERENT_ES` | `cmf_flux_cntrl_var_mod.f:297`, read `rd_cmf_flux_controls.f:107` (`COH_ES`) | requested / active coherent-e.s. flag | logical |
| `USE_OLDJ_FOR_ES` | `cmf_flux_cntrl_var_mod.f:306`, read `rd_cmf_flux_controls.f:109` (`OLD_J`) | warm-start Jbar seeding from an old J file | logical |

`RJ_ES` is a mean-intensity-like quantity, not an emissivity. The
scattered emissivity is `ESEC * RJ_ES`; `ESEC` is folded in at
emissivity construction, never into `RJ_ES` itself. `RJ_ES` is a
read (not iterated) input to the source function
(`comp_j_blank.f:790`).

## 2. Physics of the construction

The emissivity entering every formal solve is

    eta_tot = eta_thermal + sigma_es * Jbar_es
              + (chi_Rayleigh + chi_dust) * J

where `Jbar_es` is the e.s.-convolved mean intensity (coherent
limit `R -> delta`: `Jbar_es = J`). The source function used by
the moment/variable-Eddington solvers is

    S = eta/chi + theta_scat * Jbar_es,   theta_scat = chi_scat/chi

which is the standard two-stream scattering source
(Mihalas 1978; Hubeny & Mihalas 2014). Incoherent e.s. obeys
`eta_es(nu) = sigma_es int R(nu',nu) J(nu') dnu'` - the Rybicki &
Hummer (1992, ALI-I) form of angle-averaged redistribution.

`subs/comp_j_conv_v2.f` computes `Jbar_es`:

- Two-term exponential-in-ln-nu fit to the dipole e.s.
  redistribution function (fit coefficients `ACOEF`, `BCOEF`,
  `comp_j_conv_v2.f:114-120`, citing Rybicki & Hummer A&A 290,553 -
  i.e. RH 1994, not the 1992 ALI paper; the doc line at :8 says
  "Rybicki and Hummer" without a full cite).
- Each fit term solves a tridiagonal system over the frequency
  grid (`A x_{i-1} - [H+A+C] x_i - C x_{i+1} = D`, :66-73,
  :282-294) - the integral approximated as a second-derivative
  operator inverse, applied to `J(nu)` and to `B_nu(T)`.
- Planck-function preservation: wavelength shift on the Wien
  side + scaling `J_ES *= exp(logB - logB_es)` (:296-354) so that
  `Jbar_es[B] = B` at depth - the correction compensates for the
  missing Compton and stimulated-e.s. terms in the RH formalism
  (:9-13).
- Output is stored in `J_STORE` (overwriting J, :362) and either
  returned in the call or written to the `ES_J_CONV` direct-access
  file in EDDFACTOR format (:371-390).

## 3. Call chain

### 3.1 Model executable (cmfgen_sub + comp_j_blank)

```
cmfgen_sub iteration loop (NUM_ITS_TO_DO great-iterations)
  -> continuum frequency loop, per frequency ML:
       COMP_OPAC -> CHI_SCAT = ESEC + CHI_RAY          (comp_opac.f:171-185)
       comp_j_blank.f:
         if .NOT. COHERENT_ES:
           READ ES_J_CONV rec -> RJ_ES                (:781-800)
         TA = ETA_CLUMP + CHI_SCAT_CLUMP * RJ_ES      (:843 non-coherent;
                                                       :841 coherent uses RJ)
         formal solve (FG_J_CMF_V13 / NONMONO_V1 /
                       PP_FORM_CMF_V2 / CMF_FORMAL_REL_V*)  -> RJ
         MOM_* call: TA rebuilt, same convention        (:932-936)
         SOURCE = ZETA + THETA * RJ_ES                (:1244)
         WRITE EDDFACTOR rec <- RJ                    (:1251-1252)
       STEQ_T += FQW * ESEC * (RJ - RJ_ES)            (:2697 heating term)
  -> end of iteration:
       COMP_J_CONV_V2(EDDFACTOR RJ -> ES_J_CONV)      (:4783-4788)
       (FLUX_CAL_ONLY variant: convolve then GOTO 20000
        for another iteration, :4540-4558)
```

The model-side e.s. iteration cadence is the main iteration loop:
`RJ_ES` is refreshed once per great-iteration. `NUM_ITS_TO_DO`, not
`NUM_ES_ITERATIONS`, controls it. Within a great-iteration the
`RJ_ES` read is fixed - the classic Lambda-iteration lag.

### 3.2 Observer-frame executable (cmf_flux_sub_v5)

```
rd_cmf_flux_controls: NUM_ES_ITERATIONS, COH_ES, OLD_J
DO_NONMONO_KNOP guard: NUM_ES>=2 and COH_ES=F else STOP (:654-662)

DO SOB_PHASE = 1, NSOB_PHASE                       (:1556; 2 only under
                                                     DO_SOB_FLUX_OBS)
  warm start: USE_OLDJ_FOR_ES reads OLD_J_FILE ->
              ES_J_CONV via COMP_J_CONV_V2            (:1528-1534)
  DO ES_COUNTER = 1, NUM_ES_ITERATIONS               (:1564)
    LST_ITERATION = (ES_COUNTER == NUM_ES)           (:1568)
    CONTINUUM frequency loop, per ML:
      opacities: iteration 1 -> OPACITIES include;
        later iterations -> reload ETA/CHI/RAY/DUST
        from *_CMF_ST stores, rebuild
        CHI_SCAT = ESEC + RAY + DUST                 (:1890-1902)
      obs/COMP_JCONT_V4.INC per frequency:
        coherent:   TA = ETA_CLUMP + CHI_SCAT_CLUMP*RJ       (:749)
        non-coh:    TA = ETA_CLUMP + ESEC_CLUMP*RJ_ES
                         + (CHI_RAY_CLUMP+CHI_SCAT_DUST_CLUMP)*RJ (:751-752)
        USE_FORMAL_REL branch rebuilds TA identically (:774-777)
        solver: FG_J_CMF_V12 / PP_FORM_CMF_V2 /
                CMF_FORMAL_REL_V* / FG_J_CMF_NONMONO_V1 (KNOP) -> RJ
        RJ_ES read from ES_J_CONV when non-coherent   (:731-740)
      store RJ_CMF_ST(:,ML) = RJ every iteration      (:2253)
      ETA/CHI_CMF_ST stored on ES_COUNTER==1 only     (:2247-2251)
      IB_J_CMF_ST=RJ(ND), IB_H_CMF_ST=HFLUX_AT_IB     (:2254-2255)
    KNOP: if DO_NONMONO_KNOP .and. ES_COUNTER<NUM_ES
          FG_KNOP_REFRESH_V1()                        (:2280-2283)
    SOB moment path: MOM_KNOP_REFRESH_V1 mirror       (:2293-2296)
    end of sweep:
      COHERENT_ES=RD_COHERENT_ES                      (:2361)
      if .NOT. COHERENT_ES .and. .NOT. USE_FIXED_J
        .and. SOB_PHASE==1:
        COMP_J_CONV_V2(EDDFACTOR -> ES_J_CONV)        (:2367-2375)
  END DO ES_COUNTER                                  (:2379)
  render fold-in (once per phase):
    ETA_RND += RAY_CMF_ST * RJ_CMF_ST   (raw J)      (:2553)
    ETA_RND += DUST_CMF_ST * RJ_CMF_ST  (raw J)      (:2561)
    if non-coherent: RJ_CMF_ST <- conv(RJ_CMF_ST) in place (:2568-2573)
    ESEC = 6.65e-15*ED                                (:2575)
    ETA_RND += ESEC * RJ_CMF_ST         (convolved)  (:2577-2579)
    ETA_RND,CHI_RND *= CLUMP_FAC if DO_CLUMP_MODEL   (:2597-2602)
  OBS_FRAME_SUB_V9(ETA_RND_ST, CHI_RND_ST, ...)      (:2688 phase 1)
  ... SOB line deposit loop (phase 1) ...
  OBS_FRAME_SUB_V9(ETA_RND_ST, CHI_RND_ST, ...)      (:3108 phase 2)
END DO SOB_PHASE                                     (:3132)
```

### 3.3 obs_frame_sub_v9 itself

`obs_frame_sub_v9.f` recomputes `ESEC = 6.65e-15*ED` (:419) and
uses it only for the `TAU_ES`/`DTAU_ES` diagnostics (:426, :537)
reported through `OBS_ES_DTAU`. It adds no scattering emissivity
or redistribution term of its own: the e.s. contribution arrives
already folded into `ETA_RND_ST`. One fold total.

## 4. ESEC / RJ_ES folding census

Every multiplication or assignment involving ESEC and RJ_ES in the
production chain:

| # | Site | Operation | Quantity folded | Times applied |
|---|---|---|---|---|
| 1 | `comp_opac.f:181` | `CHI_SCAT = ESEC + CHI_RAY` | opacity assembly, model side | once per frequency |
| 2 | `OPACITIES_V4.INC:41-42` | `CHI_SCAT = ESEC + CHI_RAY`; `CHI = CHI_SCAT` | opacity assembly, obs V4 | once |
| 3 | `OPACITIES_V5.INC:55-56` | `CHI_SCAT = ESEC+CHI_RAY+DUST`; `CHI = CHI_SCAT + CHI_ABS_DUST` | opacity assembly, obs V5 | once |
| 4 | `comp_j_blank.f:843` | `TA += CHI_SCAT_CLUMP*RJ_ES` | emissivity, model non-coherent | once per freq per inner it. |
| 5 | `comp_j_blank.f:841` | `TA += CHI_SCAT_CLUMP*RJ` | emissivity, model coherent | once |
| 6 | `comp_j_blank.f:935` | `TA += CHI_SCAT_CLUMP*RJ_ES` | emissivity for MOM_* calls, non-coherent | once |
| 7 | `comp_j_blank.f:1292` | `TA += CHI_SCAT_CLUMP*RJ_ES` | emissivity for CMF_FORM_SOL obs-render candidate (model side, restricted predicate) | once |
| 8 | `comp_j_blank.f:573,1244` | `SOURCE = ZETA + THETA*RJ_ES` | source function (= eta/chi + chi_scat/chi * Jbar) | once |
| 9 | `do_var_cont.f:314` | `TA = ETA_CLUMP + ESEC_CLUMP*RJ_ES` | Jacobian-side emissivity | once |
| 10 | `do_var_cont.f:603` | `TB = RJ_ES` | convolved J handed to the linearization | once |
| 11 | `do_var_cont.f:756` | `ESECEXT*RJEXT_ES` | extended-grid emissivity for VAR_MOM_JEXT_CMF_V3 | once |
| 12 | `cmfgen_sub.f:2697,6427` | `STEQ_T += FQW*ESEC*(RJ-RJ_ES)` | net e.s. energy exchange in temperature equation | once per freq |
| 13 | `obs/COMP_JCONT_V4.INC:751-752` | `TA += ESEC_CLUMP*RJ_ES + (CHI_RAY+DUST)*RJ` | emissivity, obs non-coherent | once per freq per ES it. |
| 14 | `obs/COMP_JCONT_V4.INC:749` | `TA += CHI_SCAT_CLUMP*RJ` | emissivity, obs coherent | once |
| 15 | `obs/COMP_JCONT_V4.INC:775-777` | USE_FORMAL_REL branch: same convention rebuilt inside | emissivity | once |
| 16 | `obs/COMP_JCONT_V4.INC:553-556` | `SOURCE = ZETA + (ESEC*RJ_ES + CHI_RAY*RJ)/CHI` | source function | once |
| 17 | `cmf_flux_sub_v5.f:2577-2579` | `ETA_RND += ESEC * conv(RJ_CMF_ST)` | render emissivity, e.s. part | once per phase |
| 18 | `cmf_flux_sub_v5.f:2553,2561` | `ETA_RND += (RAY or DUST)*RJ_CMF_ST` | render emissivity, coherent parts | once per phase |
| 19 | `cmf_flux_sub_v5.f:1976-1997` | `CHI` floors at `0.1*ESEC`-relative guards | positivity guard on opacity, not emissivity | guard |
| 20 | `cmf_flux_sub_v5.f:2155,2437` | `LINE_FLUXMEAN/ESEC` | line-force diagnostic normalization | diagnostic |
| 21 | `obs_frame_sub_v9.f:419-426,537` | `TAU_ES`,`DTAU_ES` from local `ESEC=6.65e-15*ED` | e.s. optical-depth diagnostic | diagnostic only |

No site folds ESEC into RJ_ES more than once; no site multiplies
RJ_ES by ESEC twice in the same emissivity. The census shows
exactly-one-fold structure throughout.

## 5. Coherent versus non-coherent convention

| Path | Coherent (R=delta) | Non-coherent |
|---|---|---|
| model `comp_j_blank` emissivity | `eta + CHI_SCAT*RJ` | `eta + CHI_SCAT*RJ_ES` |
| model `do_var_cont` emissivity | same via CHI_SCAT | `eta + ESEC*RJ_ES` |
| obs `COMP_JCONT` emissivity | `eta + CHI_SCAT*RJ` | `eta + ESEC*RJ_ES + (RAY+DUST)*RJ` |
| obs render `ETA_RND` | `eta + CHI_SCAT*RJ` (via fold-ins) | `eta + RAY*RJ + DUST*RJ + ESEC*Jbar` |

The e.s. part is consistent everywhere: `sigma_es` multiplies the
convolved J exactly once. The difference is Rayleigh (and model-side
dust, which does not exist there):

- `comp_j_blank` non-coherent uses `CHI_SCAT*RJ_ES` where
  `CHI_SCAT = ESEC + CHI_RAY`: Rayleigh photons are convolved with
  the electron-scattering redistribution function.
- The obs side and `do_var_cont` use `ESEC*RJ_ES`, keeping Rayleigh
  (and obs-side dust) on the raw coherent `RJ`.

With `INC_RAY=F` (the default; `CHI_RAY=0`) the two conventions are
identical. With `INC_RAY=T` they differ: the model treats Rayleigh
as incoherent-by-e.s., the render treats it as coherent. Physically
Rayleigh scattering IS nearly coherent (much smaller frequency
redistribution than e.s.), so the obs/do_var_cont convention is the
better physics; the model-side form is the questionable one. Impact
is limited to models run with `INC_RAY=T` - check before calling
this a production defect.

## 6. NUM_ES_ITERATIONS refresh timing

- `RD_COHERENT_ES=T` forces `NUM_ES_ITERATIONS=1` at run time
  (`cmf_flux_sub_v5.f:1539`) - no redistribution to iterate.
- Non-coherent: each `ES_COUNTER` sweep solves J over all
  frequencies with `RJ_ES` read from `ES_J_CONV` (written after the
  previous sweep). Sweep k uses `Jbar_{k-1}` - the intended
  Lambda lag. After every sweep `COMP_J_CONV_V2` rewrites
  `ES_J_CONV` (`:2367-2375`), including after the last sweep (the
  call is gated on `SOB_PHASE==1` and `COHERENT_ES`, not on
  `ES_COUNTER`; the extra write is harmless - it leaves a fresh
  warm-start file for subsequent calls/`USE_OLDJ_FOR_ES`).
- The render uses `Jbar_N` (the freshest convolution), while the
  last sweep's emissivity used `Jbar_{N-1}`. If the iteration has
  converged these coincide; if NUM_ES is too small, the rendered
  spectrum is evaluated at the better J but the last solve was not
  self-consistent. There is **no convergence gate**: NUM_ES is a
  fixed count, and nothing checks `|Jbar_N - Jbar_{N-1}|`. That is
  the P8.2 maturity question, not a double-fold defect.
- `USE_OLDJ_FOR_ES` (`:1528-1534`): when no `ES_J_CONV` exists, an
  old J file is convolved to seed `ES_J_CONV`, letting iteration 1
  start from a finite Jbar rather than the `FIRST_FREQ` zero-fill.
- Model side: refresh is once per great-iteration
  (`cmfgen_sub.f:4783-4788`), so `RJ_ES` is always one
  great-iteration stale by design - consistent with the
  Lambda-iteration structure; the model loop runs to convergence
  anyway.

Verified per-sweep call sequence (OBS path, non-coherent):

1. `ES_J_CONV` open `STATUS='OLD'` (`:1603-1608`); if the open
   fails, `COHERENT_ES=T` for that sweep - a cold start, not a
   stale read.
2. Per frequency `ML` (continuum loop to `10000 CONTINUE`,
   `:2257`): `COMP_JCONT` reads record `ML` of `ES_J_CONV` into
   `RJ_ES` - sweep k consumes `Jbar_{k-1}` - builds
   `TA = eta + ESEC*RJ_ES + (CHI_RAY+DUST)*RJ`, solves, stores
   `RJ_CMF_ST(:,ML)=RJ` (`:2253`) so the array always holds the
   CURRENT sweep's J.
3. After the frequency loop, non-final sweeps only:
   `FG_KNOP_REFRESH_V1` (`:2280-2283`, DO_NONMONO_KNOP) /
   `MOM_KNOP_REFRESH_V1` (`:2293-2296`, DO_NONMONO_SOB without
   USE_LAM_ES) - the lagged I_NEXT/J_NEXT matures between
   sweeps.
4. Then, still inside sweep k: `COMP_J_CONV_V2(WM,...,
   'EDDFACTOR',...)` (`:2367-2375`) convolves the J just written
   to EDDFACTOR records during this sweep and rewrites
   `ES_J_CONV` - producing `Jbar_k` for sweep k+1. Refresh
   therefore precedes the convolve-write in the same
   inter-sweep block; both run inside the ES loop before
   `END DO` at `:2379`.
5. After the loop: the render fold-in convolves `RJ_CMF_ST`
   (=J_N) in place (`:2568-2573`) and folds `ESEC*Jbar_N` into
   `ETA_RND_ST` (`:2577-2579`) - the freshest convolution
   reaches `OBS_FRAME_SUB_V9`.

No stale-RJ_ES path to the render exists in the covered chain:
the only file consumer is the next sweep's `COMP_JCONT`, and
`USE_FIXED_J` skips both the convolve write (`:2367` gate) and
the iteration it would feed. `DO_SOB_FLUX_OBS` phase 2 forces
coherent and convolves in memory (item 4 of section 7).

## 7. KNOP / non-monotonic interactions

1. The KNOP formal solver `FG_J_CMF_NONMONO_V1` receives `TA`
   already containing `ESEC_CLUMP*RJ_ES` (`comp_j_blank.f:897`,
   obs include `:816`) - same convention, one fold.
2. `FG_KNOP_REFRESH_V1` runs only between ES sweeps and is skipped
   on the last iteration (`cmf_flux_sub_v5.f:2280-2283`). The
   lagged I_NEXT field therefore needs `NUM_ES>=2` to mature at
   all on the OBS path. This is already enforced fail-closed:
   `rd_cmf_flux_controls.f:654-662` STOPs on
   `DO_NONMONO_KNOP .and. (NUM_ES<2 .or. RD_COHERENT_ES)`.
3. `MOM_KNOP_REFRESH_V1` mirrors the refresh for the
   `DO_NONMONO_SOB .and. .NOT. USE_LAM_ES` moment path
   (`:2293-2296`).
4. `DO_SOB_FLUX_OBS` phase 2 forces `COHERENT_ES=.TRUE.` per sweep
   (`:1591`) and skips the file convolve (`:2367` gate on
   `SOB_PHASE==1`); the line-bearing arrays reuse the in-memory
   `RJ_CMF_ST` convolve at render (`:2568-2573`), so phase-2 wing
   light is redistributed through memory, not the stale file.
   Documented in the :1541-1552 comment block.
5. `DO_ES_ALI` (comp_j_blank `:829-834`) is an
   Olson-Auer-Buchler accelerator gated on `COHERENT_ES` - it is
   inert under non-coherent e.s. and does not alter the RJ_ES
   convention.
6. `SOB_PHASE==2` rebuilds `CHI_SCAT` per frequency from stored
   arrays with the deposited line opacity (`:1890-1902`), so the
   phase-2 e.s. opacity is the same `ESEC`; deposited line opacity
   lands in `CHI_CMF_ST`, not in `CHI_SCAT` - no double counting.

## 8. Confirmed defects versus unresolved questions

**Confirmed defects: none in the folding chain.** Every ESEC x
RJ_ES site applies the product exactly once as the scattered
emissivity; the observer-frame transfer receives the folded
emissivity and adds nothing further; the coherent limit reduces
to `sigma*J` consistently on both sides. The folding algebra
is now backed by the executable battery in section 8.5
(87 checks, all passing at the measured tolerances) plus the
static call-chain verification in sections 6-7.

**Unresolved questions / audit items:**

1. (Q-RAY) RESOLVED quantitatively 15-Nov-2026 (manufactured,
   `drv_es_conv qray` + `drv_es_render` through the real
   `OBS_FRAME_SUB_V9`; `test_es_conv.py`). `comp_j_blank`
   non-coherent `CHI_SCAT*RJ_ES` convolves Rayleigh with the e.s.
   redistribution, while the obs side and `do_var_cont` keep
   Rayleigh coherent on `RJ`. Measured two-regime bound:
     - smooth continuum, thermal-dominated ETA: `max|DTA/TA|`
       <=1.8e-8, `max|dF/F|` <=7.2e-11 - negligible;
     - scattering-dominated ETA (~1% absorption) with a narrow
       emission feature inside the Rayleigh band:
       `max|DTA/TA|` = 1.5e-4 / 5.6e-4 / 8.6e-4 and
       `max|dF/F|` = 1.6e-3 / 6.0e-3 / 9.2e-3 at
       chi_ray/sigma_es = 0.1 / 0.5 / 1.0 - measurable, scaling
       roughly linearly with the Rayleigh fraction.
   The convention difference is consequential exactly where the
   physics says it must be: Rayleigh a large fraction of the
   scattering opacity AND a J field with redistribution contrast
   at feature scales (narrow line light, sharp edges - the
   SOB-wing regime). Inert while `INC_RAY=F`; the obs-side form
   (Rayleigh coherent) is the physically motivated convention,
   but no source change is proposed here - the model-side
   expression feeds the Jacobian through `comp_j_blank`, so a
   convention fix must first map the solver implications
   (which derivatives assume `CHI_SCAT*RJ_ES`).
2. (Q-MATURITY) No acceptance gate ties `NUM_ES_ITERATIONS` to a
   measured `|Jbar_N - Jbar_{N-1}|` residual. Convergence is
   assumed from a fixed count. This is the P8.2 item.
3. (Q-WASTE) `COMP_J_CONV_V2` also runs after the last ES sweep
   (`:2367` not gated on `ES_COUNTER`): an extra full-grid
   tridiagonal convolve + file write whose only consumer is the
   next invocation's warm start. Harmless, slight waste.
4. (Q-STEQ-SIGN) `STEQ_T += FQW*ESEC*(RJ - RJ_ES)`
   (`cmfgen_sub.f:2697,6427`): the sign of the e.s. term in the
   temperature equation depends on the `STEQ_T` convention
   (heating-positive vs net-cooling-positive). Not verified here;
   flagged for the P8 energy-budget check rather than asserted
   wrong.
5. (Q-DEAD-CODE) `do_var_cont.f:275-309` contains commented-out
   ES_COH_VEC variants for a `RJ vs RJ_ES` correction - dead code
   noting an abandoned alternative; document only.
6. (P7-CARRYOVER) The linear-Eddington inner-boundary law
   `I(mu)=J+3*mu*H` requires `3*|H| < J` for positivity over the
   full comoving hemisphere. The manufactured tests satisfy this
   by construction; production moments are unverified. Record as a
   production-model validity diagnostic to evaluate during P8 or
   final validation. No clamp without evidence of a real
   violation and a physical derivation.
7. (Q-LEGACY) `disp/` carries an older e.s. path (`disp/maingen.f`,
   `disp/OPACITIES.INC:102`) outside the production chain; not
   audited here.
8. (Q-CONSUMERS) `spec_plt/plt_rj.f`, `spec_plt/new_plt_rj.f`,
   `spec_plt/subs/do_clump_shift.f` read `RJ_ES`/`ES_J_CONV` as
   downstream consumers; they do not feed the flux chain.
9. (Q-CITE) `comp_j_conv_v2.f:114` cites "Rybicki and Hummer A&A
   290,553" - verified: Rybicki, G. B. & Hummer, D. G. 1994,
   'An Accelerated Lambda Iteration Method for Multilevel
   Radiative Transfer. III. Noncoherent Electron Scattering',
   A&A 290, 553-562 (arXiv astro-ph/9404019). That paper IS the
   coefficient source (the two-term exponential fit in ln nu
   and the ACOEF/BCOEF constants) and the physical basis for
   the line-wing/continuum-edge redistribution treatment. The
   code comment is correct; the additional Planck-preserving
   shift/rescale in the same routine is a CMFGEN extension
   beyond the published fit.

## 8.5 Measured operator characterization (P8.1 acceptance tests)

`unit_tests/f90/drv_es_conv.f` drives the real `COMP_J_CONV_V2`
in 'J PASSED VIA CALL' mode; `drv_es_render.f` renders
manufactured emissivities through the real `OBS_FRAME_SUB_V9`
(V=0 static atmosphere, ND=4); `drv_es_seq.f` exercises the ES
iteration ordering; `drv_knop_guard.f` drives the production
store reader + `RD_CMF_FLUX_CONTROLS`. Harness:
`unit_tests/python/test_es_conv.py`, 87/87 PASS. Production
tree sources link as the exact `.o` files the production
executables use (`compile_and_link_prod`, ifort `-qopenmp`,
MKL sequential); only drivers and the `obs_v9_adapter` shim
compile fresh.

The production operator is NOT a discrete kernel matrix. Per
depth and per Rybicki & Hummer fit term L it solves the
Helmholtz-type tridiagonal system
`(1 - (beta/BCOEF_L)^2 * d^2/dln nu^2) D = J` and sums
`J_ES = sum_L ACOEF_L * D_L` (`comp_j_conv_v2.f:249-294`),
followed by a Wien-side wavelength shift and a `B/B_es` rescale
(`:302-354`). `beta = 1.84E-3*sqrt(T)` (`:269`) sets the
smearing width; (beta/B_L)^2 ~ 1e-6 in ln nu at 1e4 K, i.e.
the operator is a weak second-derivative smoother, narrower
than any resolvable feature on production grids.

Three implementation layers must be distinguished when
interpreting these measurements: the physical redistribution
integral (frequency-integrated, roughly conservative); the
Rybicki & Hummer (1994, A&A 290, 553; arXiv astro-ph/9404019)
two-term fitting approximation of that integral; and the
production tridiagonal solve plus the ADDITIONAL
Planck-preserving wavelength shift/rescale, which is a CMFGEN
extension, not part of the published fit. Properties of the
underlying integral (e.g. photon conservation) must not be
asserted for the implementation; the measured integrated-excess
drift below is the honest characterization.

Boundary treatment (from source, `comp_j_conv_v2.f:249-258,
302-380` post-fix): rows 1 and NCF of the tridiagonal operator
use the one-sided stencil `-2/(ln(nu_1/nu_2))^2` (and the
symmetric form at NCF); the wavelength shift is derived only
on the Wien side (`NU(I) > 1.5*TEMP(K)`), applied to all lower
frequencies as a uniform wavelength offset; the interpolation
back onto the grid runs `ML=2,NCF-1`; nodes whose shifted
bracket still contains the forced `PLANCK_NU(1)=NU(1)` value
(final `I=1`) and the upper endpoint itself receive NO shift
and NO rescale - `PLANCK_ES=PLANCK_FN` (factor unity) and the
raw convolved `J_ES` stands, because no genuine shifted
`PLANCK_ES` exists there and the `B/B_es` ratio is then not a
calibration (sec 8.5.1). The lower endpoint keeps its raw
`PLANCK_ES` (a benign ~unity factor on the smooth
Rayleigh-Jeans side). Interior edge points extrapolate if
needed; the whole log-space shift/interp/rescale block is
gated on `MINVAL(PLANCK_FN)>0 .AND. MINVAL(PLANCK_ES)>0`
(`:302`), so a depth whose scaled Planck underflows to
exactly 0 escapes the rescale entirely. This is the treatment
the `bnd` test measures.

Feature coordinate note: the manufactured line feature is a
Gaussian in nu, while the operator acts in ln nu. The local
approximation sigma_ln_nu ~ sigma_nu/nu0 = 0.200 understates
the actual ln(nu) width; the measured second-moment width is
0.242 (NCF=256) / 0.245 (NCF=512), i.e. 6.4 / 13.0 grid points
per actual sigma. The kernel width itself ((beta/B_L)^2 ~ 1e-6
in ln nu at 1e4 K) is far below the grid cell; this battery
resolves the FEATURE, not the kernel, and no kernel-width
claim is made.

Measured battery results (87/87 PASS):

| case | result |
| --- | --- |
| planck `J=B(T)` | max rel 1.46e-5, rms 1.28e-6 (significant-J support; Wien-underflow masked) |
| const | Jbar >= 0; B/B_es rescale documented; 0 zeros at Wien edge |
| linepair 256 | peak 3.9986e-4 -> 3.9978e-4; core -7.5e-8 / wing +2.1e-8; max|dexc|/pk 2.25e-4; full-domain integrals below; no ringing |
| linepair 512 | peak 3.9983e-4 -> 3.9975e-4; core -7.5e-8 / wing +2.3e-8; max|dexc|/pk 2.26e-4; resolution-consistent |
| bnd 256/512 | interior vs extended max 1.30e-4/1.34e-4; lower edge max 1.0e-4/2.0e-4, rms 4.2e-5/8.1e-5; upper edge endpoint 8.3e-3/1.6e-2, feature preserved (sec 8.5.1) |
| wcp 256/512 | Wien-dead non-Planck continuum: top nodes vs ext max 2.2e-4/1.8e-2; out/in bounded [1.00,1.52] |
| swp nu0=299.5..200 | feature peak ratio 0.983-0.999 at all positions incl. ~0.5 cells from endpoint |
| edge (step) | max rel 9.1e-4, no ringing |
| tempu 5e3/1e4/5e4 K (uniform T) | max 1.55e-4 / 1.82e-13 / 4.46e-14; beta/B_L reported per case; same mask for max and rms |
| pos (zeroed interior band) | Jbar >= 0 |
| grid NCF 128/256/512 | max err 1.32e-5 / 1.46e-5 / 1.66e-5 - saturates (operator width is physical) |

Full-domain line-excess integrals (linepair, all nodes - not
restricted to the input feature mask; ascending-nu trapezoid):

| measure | NCF | in | out | drift |
| --- | --- | --- | --- | --- |
| int(exc) dlnnu | 256 | 2.0978e-4 | 2.0978e-4 | -3.96e-6 |
| int(exc) dlnnu | 512 | 2.0978e-4 | 2.0978e-4 | -3.45e-6 |
| int(exc) dnu   | 256 | 4.0112e-4 | 4.0112e-4 | +4.47e-7 |
| int(exc) dnu   | 512 | 4.0105e-4 | 4.0105e-4 | +1.09e-6 |

On the extended grid (linext) the dlnnu drift is at machine
precision (~1e-15), and the production-vs-extended full-domain
output integral agrees to +2.3e-5 relative. The measured
drift is reported, not elevated to a conservation law: the
interior column sums of the discrete operator are -H = 1 but
the one-sided endpoint rows do not share that property, and
the pointwise B/B_es rescale is not an integral-preserving
step. The output excess outside the input >5%-of-peak
footprint is 2.7% (256) / 2.2% (512) of int(exc_out) dlnnu -
the redistributed wing content. FWHM_ln(nu): 0.415 -> 0.452
(256), 0.452 -> 0.452 (512); the change is below one grid
cell at both resolutions and is reported as unresolved (the
assertion is that the feature does not NARROW by more than
two cells - a smearing operator must not sharpen).

Extended-domain (linext) boundary sensitivity of the smeared
feature: exc_out restricted to the shared production nodes
differs from the production-grid exc_out by max 8.1e-4 (256) /
1.5e-3 (512) on the significant mask - the boundary influence
on an interior feature is ~1e-3, an order below the feature
itself.

### 8.5.1 Upper-edge Wien-rescale suppression (diagnosed, corrected)

The live-feature boundary test exposed a real production
defect, confirmed against the extended-domain reference and
corrected in `comp_j_conv_v2.f` (this commit):

- With a significant field alive at the upper production bound
  (manufactured Gaussian at nu~260, amplitude ~0.5), the
  pre-fix production-grid Jbar at the high-nu end was
  suppressed to <=1e-2 of the extended-domain result at
  T >= 1e4 K; at T = 2e4 K and 1e4 K the output was crushed
  to ~1e-82 and ~1e-273 respectively (max rel discrepancy 1.0).
- The suppression was one-sided (no amplification, no
  ringing), and the depth whose scaled Planck underflows to
  exactly 0 (T=0.5e4) escaped the rescale per the
  `MINVAL(PLANCK_FN)>0` gate and agreed to <=2e-4.

Mechanism (source-traced, pre-fix `comp_j_conv_v2.f:302-354`):
in the Wien-dead zone the Helmholtz-solved `B_es` at the
endpoint is enormously larger than `B` (measured B_es ~ 1e-38
vs B ~ 1e-57 at T=5e4; ~1e62 vs ~1e-316 at T=1e4). Two
coupled defects followed:

1. `PLANCK_NU(1)=NU(1)` is FORCED, not derived - the endpoint
   anchor violates the `ln B_es(NU) = ln B(PLANCK_NU)`
   relation the interior map satisfies by construction.
2. The log-interpolation runs `ML=2,NCF-1`; the endpoint and
   the nodes whose shifted bracket `[PLANCK_NU(I),
   PLANCK_NU(I+1)]` still contains the forced value (final
   `I=1`, i.e. `NU(ML) > PLANCK_NU(2)`) therefore keep or
   interpolate against the unshifted `ln B_es`, while the
   final `EXP(J_ES + PLANCK_FN - PLANCK_ES)` applies the
   `B/B_es` factor to ALL nodes. The corrupted zone is
   contiguous at the top (measured: nu >= 248, 6 nodes at
   NCF=256; nu >= 283, 4 nodes at NCF=512 - identical to the
   `|logF|>0.01` set). For a Planck input the endpoint factor
   is algebraically exact (J_es = B_es divides back out), so
   the defect is invisible to the LTE test but destroys any
   non-Planck field by ~1e-20 to ~1e-280.

Options evaluated against the extended-domain reference:

- Extend the shift interpolation to the endpoints (rejected):
  the endpoint's correct shifted position is off the grid.
  Deriving PLANCK_NU(1) gives nu~252 (NCF=256, T=5e4) - the
  smeared B_es is smallest AT the boundary, so no on-grid
  node supplies B_es(NU1) = B(NU1). Naive extension
  extrapolates the log-interp with weight T1 = -9.7 and still
  suppresses the live field by a factor ~360 - uncontrolled
  extrapolation, exactly the failure the boundary condition
  warns about.
- Restrict the rescale where the shifted evaluation is not
  genuine (selected): at the endpoint and the I=1-bracket
  nodes no valid shifted `PLANCK_ES` exists, so the `B/B_es`
  factor is not a calibration - it is the convolution's own
  boundary error applied to an arbitrary J. Setting
  `PLANCK_ES = PLANCK_FN` (factor unity) and `J_ES = A(ML)`
  (the raw convolved value) there reproduces the
  extended-domain treatment: on the extended grid the same
  nodes sit inside the domain and, when B underflows at the
  extended top, the whole block is skipped and the output is
  the raw J_ES.
- Consistent endpoint treatment on both fields (equivalent
  statement of the selected fix): where the map cannot
  produce a genuine shifted evaluation, no calibration data
  exists and the operator output stands.

Corrected measurements (post-fix):

| case | before | after |
| --- | --- | --- |
| bnd 256 upper edge | max rel ~1.0 (1e-20..1e-273 crush) | max 8.3e-3 (T=5e4), 3.4e-3/1.7e-3/8.4e-4 per T; endpoint 8.3e-3 vs ext |
| bnd 512 upper edge | same suppression | endpoint 1.6e-2; corrupted-zone nodes raw (rel <=5e-8..1.6e-2); deeper Wien-side nodes carry the designed shift (zone bounded factor 0.69-1.02, not orders) |
| wcp continuum (top nodes) | n/a | max 2.2e-4 (256) / 1.8e-2 (512) vs extended; out/in in [1.00,1.52] - bounded designed shift, no loss |
| swp positions nu0=299.5..200 | n/a | peak ratio 0.983-0.999 at all positions incl. 0.5 cells from the endpoint; endpoint node ratio ~1.0 where significant |
| planck sig support | 1.46e-5 | 1.46e-5 (unchanged; masked dead nodes now emit raw B_es) |
| const Wien zeros | 3 | 0 |

Residual, documented: at the endpoint a Planck input now
returns the convolved B_es (~1e-38 at T=5e4) instead of the
algebraically-exact B (~1e-57) - LTE invariance is traded for
non-thermal preservation only where both values are denormal
(emitted field ~1e-38 vs ~1e-57, both << any real field). The
interior Wien-side wavelength-shift distortion of non-Planck
fields (up to ~50% on a steep ramp at nu~110, T=1e4) is the
stock CMFGEN correction and is unchanged.

Latent risk noted (pre-existing, unchanged): the PLANCK_NU
match loops `DO WHILE(PLANCK_ES(I) .GT. PLANCK_FN(ML))` /
`.LT. PLANCK_FN(ML-1)` at `:310-315` carry no `ML <= NCF`
bound; if `PLANCK_ES(I)` lies outside the entire `PLANCK_FN`
range the index overruns the array. Not triggered on the
manufactured or validation grids; recorded for P8.2.

Test-harness lesson (driver fix, this round): the first
linext run showed a spurious ~6.4e-4 excess band near nu~0.2
on the extended grid only. Cause: `drv_es_conv` allocated
`J_STORE`/`J_IN` at the pre-extension NCF while the extended
grid used NCF+NCF/2, so the convolver walked past `J_STORE`
into `J_IN` (second convolve read corrupted input). Fixed by
resolving the extended NCF before ALLOCATE; the anomaly
vanished and the boundary results changed at the ~1e-4 level
(lower edge improved ~6x where the polluted zone overlapped).
The clean numbers above are post-fix.

ES-iteration sequencing (T5, `drv_es_seq.f`): scripted
two-sweep harness mirroring the production ordering on the
real convolver - sweep k consumes `Jbar_{k-1}`, stores J_k,
refresh slot, convolve stored J_k to `Jbar_k`, render folds
`ESEC*Jbar_2`. Perturbation probes show both scripted
generations propagate to the rendered TA; broken-ordering
modes are detected exactly where they must be silent:
'stale' (sweep 2 fed `Jbar_0`) is bitwise silent to the J_1
generation; 'preconv' (convolve before the store) identically;
'oldrend' (render folds `Jbar_1`) is bitwise silent to the
J_2 generation. REFRESH_COUNT=1 between sweeps. This harness
mirrors the audited ordering; it is NOT a link of
`CMF_FLUX_SUB_V5` and does not by itself cover the file-slot
mechanics (`ES_J_CONV` direct-access records), which remain
source-verified.

KNOP ES guard (T6, `drv_knop_guard.f`): the real
`RD_OPTIONS_INTO_STORE` + `RD_CMF_FLUX_CONTROLS` on a
manufactured control file (production VADAT key lines + OBS
keys). `KNOP=T, NUM_ES=1, COH_ES=F` and `KNOP=T, NUM_ES=2,
COH_ES=T` both STOP with the named
'DO_NONMONO_KNOP=T requires NUM_ES>=2 and COH_ERENT_ES=F'
diagnostic; `KNOP=T, NUM_ES=2, COH_ES=F` returns cleanly
(GUARD_PASS). Verified by message content, not exit status.

Emissivity/render checks (manufactured fields, real
`OBS_FRAME_SUB_V9`):

- Both emissivity conventions rebuilt elementwise to <=1e-12
  of the emissivity scale; `DTA = CHI_RAY*CLUMP*(RJ_ES-RJ)`
  exact to machine noise.
- Coherent-vs-non-coherent emergent flux differs by
  8.4e-3..1.4e-2 - the convolved fold measurably propagates to
  the observed flux (T4 discriminator).
- All render flux differences trace to the supplied TA variant
  alone: the transfer adds no second scattering term (T7).

## 9. Equations verified against references

- Emissivity `eta = eta_th + sigma_es*Jbar_es` (+ coherent
  Rayleigh/dust) matches the scattering source-function structure
  of Mihalas (1978) and Hubeny & Mihalas (2014): scattering
  contributes `sigma*J` coherent, `sigma*int R J` incoherent.
- `S = zeta + theta*Jbar`, `theta = chi_scat/chi` is the standard
  two-term source used by the Eddington-factor and moment solvers.
- `RJ_ES = int R J` as a *field* input iterated between sweeps is
  the Lambda-iteration/ALI structure of Rybicki & Hummer (1992):
  the updated mean intensity feeds back into the emissivity until
  self-consistency.
- The Planck-preserving shift/scale in `comp_j_conv_v2` is an
  implementation device (documented in-source) to recover LTE at
  depth; it is not itself a published formula and is not
  attributed to one in the code.

## 10. Proposed implementation plan

See `docs/audits/CMF_FLUX_P8_IMPLEMENTATION_PLAN.md` for the
manufactured tests, numerical acceptance criteria, and the P8.2
maturity-gate proposal.
