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
to `sigma*J` consistently on both sides.

**Unresolved questions / audit items:**

1. (Q-RAY) `comp_j_blank` non-coherent `CHI_SCAT*RJ_ES` convolves
   Rayleigh with the e.s. redistribution, while the obs side and
   `do_var_cont` keep Rayleigh coherent on `RJ`. Inert while
   `INC_RAY=F`. Needs: an `INC_RAY=T` census of `CHI_RAY/ESEC`
   magnitudes in the production models, and a decision on the
   correct physical convention (Rayleigh is near-coherent, so the
   obs-side form is likely the intended one).
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
   290,553" (the 1994 redistribution-function fit); the task's
   listed RH-1992 ALI reference is the correct basis for the
   iteration scheme but is not the coefficient source. The code
   comment is correct as far as it goes; no formula attribution
   issue found.

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
