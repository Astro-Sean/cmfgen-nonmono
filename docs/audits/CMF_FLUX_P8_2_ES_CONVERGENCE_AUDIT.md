# P8.2 - ES redistribution convergence audit

Diagnostic-only phase: instrument the electron-scattering redistribution
iteration in CMF_FLUX to measure whether `RJ_ES` and the emergent
spectrum converge with `NUM_ES_ITERATIONS` sweeps, and whether the
final render consumes the freshest convolved field. No production
physics is changed; the fixed-sweep calculation remains the reference.

Scope and conventions follow the P8.1 audit
(`CMF_FLUX_P8_ES_REDISTRIBUTION_AUDIT.md`); line numbers verified
against the tree on 15-Nov-2026.

## 1. Production iteration sequence (verified from source)

### 1.1 File and field roles

| Symbol / file | Role | Verified |
|---|---|---|
| `EDDFACTOR` (LU_EDD) | per-frequency record `REC=ACCESS_F` holds `RJ(1:ND), FL` - the sweep's solved J | `COMP_JCONT_V4.INC:978,1141` write |
| `ES_J_CONV` (LU_ES) | same EDDFACTOR layout; record `ACCESS_F` holds `RJ_ES(1:ND), NU(ML)` - the convolved field consumed by the next sweep | `COMP_JCONT_V4.INC:731-739` read; `comp_j_conv_v2.f:531-537` write |
| `RJ_CMF_ST(ND,NCF)` | in-memory copy of the current sweep's J, stored per frequency every sweep | `cmf_flux_sub_v5.f:2271` |
| `ETA/CHI/RAY/DUST_CMF_ST` | stored on sweep 1 only | `:2265-2270` |
| `EDD_CONT_REC=3`, `INITIAL_ACCESS_REC=15`, `FINISH_REC=5` | record indices | `eddfac_rec_defs_mod.f` |

`ACCESS_F` is reset to 5 at the head of each `SOB_PHASE` iteration
(`:1569`) and marched to `INITIAL_ACCESS_REC` inside the record
machinery; frequency `ML` maps to record `INITIAL_ACCESS_REC+ML-1`.

### 1.2 Sweep-level sequence - all OBS pathways (phase 1)

```
rd_cmf_flux_controls: NUM_ES, COH_ES, OLD_J
DO_NONMONO_KNOP guard: NUM_ES>=2 and COH_ES=F else STOP
                       (rd_cmf_flux_controls.f:654-662)
IF(RD_COHERENT_ES) NUM_ES_ITERATIONS=1              (:1545)

warm start (only when USE_OLDJ_FOR_ES and no live ES_J_CONV):
    COMP_J_CONV_V2(OLD_J_FILE -> ES_J_CONV)          (:1534-1540)
    -> field generation G0 = conv(J_old)

DO SOB_PHASE = 1, NSOB_PHASE                       (:1562; 2 only
                                                     under DO_SOB_FLUX_OBS)
DO ES_COUNTER = 1, NUM_ES_ITERATIONS               (:1582)
  LST_ITERATION = (ES_COUNTER == NUM_ES)            (:1586)
  ES_J_CONV open STATUS='OLD'; open failure ->
      COHERENT_ES=T for this sweep (cold start)     (:1618-1626)
  phase 2 forces COHERENT_ES=T, COMPUTE_EDDFAC=T    (:1609,1617)

  CONTINUUM frequency loop, per ML:
    RJ_ES(1:ND)  <- READ ES_J_CONV rec ACCESS_F     (COMP_JCONT:732)
                    sweep k consumes generation k-1
    TA = ETA_CLUMP + ESEC_CLUMP*RJ_ES
         + (CHI_RAY_CLUMP + CHI_SCAT_DUST_CLUMP)*RJ (COMP_JCONT:751-752)
    formal solve -> RJ(1:ND):
        DO_NONMONO_KNOP: FG_J_CMF_NONMONO_V1        (COMP_JCONT:816)
        default (SOB render / stock): FG_J_CMF_V12  (:825)
        DO_NONMONO_SOB + USE_LAM_ES=F:
            MOM_J_CMF_NONMONO_V1                    (:920)
    WRITE EDDFACTOR rec ACCESS_F <- RJ              (:978/:1141)
    RJ_CMF_ST(:,ML) = RJ                            (:2271)

  inter-sweep refresh, skipped on last sweep:
    DO_NONMONO_KNOP:        FG_KNOP_REFRESH_V1      (:2300)
    DO_NONMONO_SOB + !LAM:  MOM_KNOP_REFRESH_V1     (:2313)

  end-of-sweep convolution, phase 1 only,
  gated .NOT.COHERENT_ES .AND. .NOT.USE_FIXED_J:
    COMP_J_CONV_V2(EDDFACTOR -> ES_J_CONV)          (:2388-2393)
    -> generation G_k = conv(J_k), written after EVERY
       sweep including the last

  [render fold-in, once per phase]
    ETA_RND += RAY*RJ_CMF_ST  (raw J)               (:2580)
    ETA_RND += DUST*RJ_CMF_ST (raw J)               (:2588)
    in-place COMP_J_CONV_V2 on RJ_CMF_ST            (:2597)
        -> RJ_CMF_ST becomes Jbar_N = conv(J_N)
    ETA_RND += ESEC * Jbar_N                        (:2616)
    OBS_FRAME_SUB_V9                                (:2726 phase 1,
                                                      :3146 phase 2)
```

### 1.3 Generation genealogy

| Generation | Produced by | Consumed by |
|---|---|---|
| G0 = conv(J_old) | OLDJ block :1534 | sweep 1 (if ES_J_CONV exists) |
| G_k = conv(J_k) | end-of-sweep convolve :2390, k=1..N | sweep k+1 reads G_k |
| Jbar_N = conv(J_N) | file write :2390 (gen N) AND in-place convolve :2597 | render ESEC fold :2616 |

Three verified facts define the convergence question:

1. Sweep `k` consumes generation `k-1`: the lambda-lag is structural.
2. The final file generation `G_N` is written at `:2390` after the last
   sweep but read by **no** further sweep - it exists only as a
   warm-start seed for any later invocation (P8.1 Q-WASTE).
3. The render does **not** read `ES_J_CONV`. It convolves
   `RJ_CMF_ST = J_N` in place (`:2597`) - the identical routine on the
   identical field. The render therefore uses the same generation-N
   field the file holds. `FINAL_RENDER_INPUT_GENERATION = N` by
   construction; the diagnostic verifies it by comparing the
   post-convolve `RJ_CMF_ST` against the saved generation-N field.

### 1.3b Model-side counterparts (inspection note)

The spec-listed model files use the same mechanism and were
verified as context:

- `new_main/mod_subs/comp_j_blank.f:782` reads `RJ_ES` from
  `ES_J_CONV` record `ACCESS_F` inside the population-solve formal
  J, with a frequency-vs-record guard (`T1 .NE. FL` -> named error,
  `:783-787`). Same read pattern as `COMP_JCONT_V4.INC:731`.
- `new_main/mod_subs/do_var_cont.f` does not reference `ES_J_CONV`
  or `COMP_J_CONV_V2` - it is the J-variation (`USE_FIXED_J`)
  machinery, not part of the redistribution chain. The model-side
  `ES_J_CONV` write lives in the great-iteration driver
  (`new_main/cmfgen_sub.f`), outside the OBS pathway audited here.

### 1.4 SOB_FLUX_OBS phase 2

`DO_SOB_FLUX_OBS` adds `SOB_PHASE=2`: each phase-2 sweep is forced
`COHERENT_ES=T` (`:1609`), the file convolve is skipped (gate at
`:2388`), and the second render folds `ESEC*conv(J_phase2)` from the
in-place convolve. Phase-2 sweeps carry **no** `RJ_ES` generations -
they iterate a coherent formal solve only. The residual diagnostic is
meaningful in phase 1 only.

### 1.5 KNOP versus SOB differences (for equivalent-generation comparison)

| | SOB render / stock | KNOP render |
|---|---|---|
| formal solver | `FG_J_CMF_V12` (clamped SIGMA) | `FG_J_CMF_NONMONO_V1` (unclamped SIGMA_TRUE) |
| inter-sweep refresh | none (V12 carries no lagged field) or `MOM_KNOP_REFRESH_V1` when `DO_NONMONO_SOB .and. .NOT.USE_LAM_ES` | `FG_KNOP_REFRESH_V1` on sweeps 1..N-1 |
| RJ_ES generations | G0..GN via same file mechanism | G0..GN via same file mechanism |
| extra convergence layer | none | lagged `I_NEXT` (compressing upwind intensity) matures across sweeps; census already printed by `GET_FG_GS_MATURITY` (:2419) |
| controls | none beyond NUM_ES/COH_ES | requires NUM_ES>=2 and COH_ES=F (STOP enforced) |

Equal `NUM_ES_ITERATIONS` does **not** imply equivalent convergence:
under KNOP the solved J additionally depends on `I_NEXT` maturity, so
`Delta_JES` alone understates the sweep-to-sweep change. The RJ_ES
residual is still the correct field comparison for both pathways;
the I_NEXT maturity census (:2419) is the companion KNOP diagnostic.
Sweep k of one pathway is generation-equivalent to sweep k of the
other - compare like-for-like generations only.

## 2. Diagnostic design

### 2.1 New production source: `obs/es_resid_diag.f`

Module `ES_RESID_DIAG_MOD` (fixed-form `.f`, matching
`obs/cmf_flux_cntrl_var_mod.f` convention) holding module state:

- `RJ_PREV(:,:)` allocatable - previous-generation field buffer
- `HAS_PREV` - whether a stored generation exists
- `GEN_COUNT` - convolve-generation counter
- `N_SWEEP_REC` - completed-sweep counter for the summary
- `N_NEG_CUM`, `N_NONFIN_CUM` - cumulative anomaly counters
- `DIAG_LU` - logical unit for `ES_RESID_DIAG` output

Four public entry points:

| Routine | Call site | Action |
|---|---|---|
| `ES_RESID_SEED` | after `ACCESS_F=5`, before `DO ES_COUNTER`; gated `SOB_PHASE==1`, diag on | read existing `ES_J_CONV` (OLDJ seed or stale file); store as G0, or mark cold start |
| `ES_RESID_UPDATE` | inside convolve block after `CLOSE(LU_ES)` (:2393, call at :2400); diag on | read fresh `ES_J_CONV` -> gen=GEN_COUNT+1; compare vs stored; write per-depth + global rows; save field |
| `ES_RESID_RENDER` | after in-place convolve `RJ_CMF_ST` (:2597, call at :2609), before ESEC fold (:2616); phase 1, diag on | compare `RJ_CMF_ST` (=Jbar_N) vs stored gen N; write render-generation verdict |
| `ES_RESID_FINAL` | after render, end of phase-1 block; diag on | summary block: CONFIGURED/COMPLETED sweeps, FINAL_REFRESH_GENERATION, FINAL_RENDER_INPUT_GENERATION, match flag, cumulative anomalies |

Two-layer split for testability: the file I/O and state management sit
in `ES_RESID_*` wrappers; the pure residual computation lives in
`ES_RESID_COMPARE(RJ_NEW,RJ_OLD,...)` so manufactured drivers call the
production comparison directly without file plumbing.

### 2.2 Residual definitions (spec section 4)

Per depth K over frequency nodes M (same grid, consecutive generations):

- `S_K = max_M max(|NEW|,|OLD|)`; `D_K = max_M |NEW-OLD|`
- `R_K = D_K/S_K`; explicit handling: `S_K>0` divide; `S_K==0,D_K==0`
  -> `R_K=0`; `S_K==0,D_K!=0` -> `R_K=1`
- `R_GLOBAL = max_K R_K` plus the depth index
- active mask: `max(|NEW|,|OLD|) > ACTIVE_FLOOR * S_K`,
  `ACTIVE_FLOOR=1e-12` default (diagnostic mask, not a tolerance);
  `R_REL(K,M)=|dJ|/max(|NEW|,|OLD|)` on active nodes; masked/active
  counts reported; "no active nodes" reported explicitly
- non-finite check on BOTH fields first: index (K,M), which
  generation, result marked invalid - never reported as converged
- negatives: counted and located per depth; no clamping

### 2.3 Output: dedicated file `ES_RESID_DIAG`

Formatted text, separate from OBSFLUX/OBSFRAME. Per sweep:

```
ES_RESID_DIAG v1
HEADER PATHWAY=<SOB|KNOP> ND=<nd> NCF=<ncf> CONFIGURED_ES_SWEEPS=<n>
       ACTIVE_FLOOR=<f>
SEED GEN=<0|NONE> FIELD_ND=<..> FIELD_NCF=<..> SOURCE=ES_J_CONV
DEPTH SWEEP=<k> K=<k> R=<rk> RREL=<r'> D=<dk> S=<sk>
      NACT=<n> NMASK=<n> NNONFIN=<n> NNEG=<n>
GLOBAL ES_SWEEP=<k> REFRESH_GENERATION=<g> GLOBAL_MAX_R=<r>
       DEPTH=<k*> GLOBAL_MAX_REL=<r'> GLOBAL_MAX_ABS=<d>
       TOTAL_ACTIVE=<n> TOTAL_INVALID=<n> TOTAL_NEG=<n>
RENDER GEN_LAST=<g> FIELD_MATCH=<T|F> MAX_ABS_DIFF=<d>
FINAL CONFIGURED_ES_SWEEPS=<n> COMPLETED_ES_SWEEPS=<m>
      FINAL_REFRESH_GENERATION=<g> FINAL_RENDER_INPUT_GENERATION=<g'>
      RENDER_MATCH=<T|F> INVALID_SWEEPS=<i> TOTAL_NONFINITE=<n>
      TOTAL_NEGATIVE=<n>
```

plus the `FIRST_GENERATION` marker on a cold-start sweep 1. Rows are
machine-readable (CSV-friendly fixed keywords).

### 2.4 Control keys (rd_cmf_flux_controls.f, cmf_flux_cntrl_var_mod.f)

| Key | Type | Default | Purpose |
|---|---|---|---|
| `ES_RESID_DIAG` | logical | F | master switch; F = zero work |
| `ES_DIAG_AFLOOR` | real | 1.0D-12 | active-frequency mask floor |
| `ES_DIAG_DEPTHS` | logical | T | per-depth rows vs global-only |

All defaults preserve production behaviour; a run without the keys is
bit-identical to the pre-diagnostic run.

### 2.5 Safety properties

- Diagnostic loops run serially in the post-sweep sequential region;
  no OMP directives added -> no races by construction. The unit test
  runs the same manufactured case under 1 and N threads and requires
  byte-identical output.
- `RJ_PREV` allocated on first use at `(ND_RES,NCF_RES)`; reallocated
  only when dimensions change; deallocated in `ES_RESID_FINAL`.
- First sweep with no seed: `HAS_PREV=F` -> reports "first
  generation, no residual" rather than comparing uninitialised data.
- No writes into any production array; reads are read-only.
- The routine is called from inside `IF(.NOT.COHERENT_ES ...)` blocks,
  so coherent/phase-2 sweeps cannot produce spurious generations.

## 3. Insertion-point rationale

The spec requires "immediately after a new ES convolution has
completed". The only point where a complete new `RJ_ES` generation is
available as a whole field is after `COMP_J_CONV_V2` has rewritten
`ES_J_CONV` and the unit is closed (`:2393`). Reading the file there is
the least invasive possible probe: zero coupling to production arrays,
no signature changes, and the record layout is already defined by the
production writer (`comp_j_conv_v2.f:531-537`).

The alternative - passing `J_STORE` out of `COMP_J_CONV_V2` - would
touch a production signature for diagnostic-only benefit and was
rejected. Per-frequency interception inside `COMP_JCONT` would give
the *consumed* field, not the *produced* one - useful only as a
secondary probe, not needed since the produced field is what the spec
defines `Delta_JES` on.

## 4. Spectral convergence experiment design

Controlled NUM_ES sweep on one physical model per pathway:

- model: `m15_2p82d` (existing `cmfgen_validation/m15_2p82d_*` cells;
  converged populations reused so only the OBS runs vary)
- `NUM_ES` = 1,2,3,4,6,8,12 per pathway; KNOP floor NUM_ES>=2 already
  enforced by `rd_cmf_flux_controls.f:654-662` -> KNOP runs 2..12
- identical model input, warm-start disabled (`OLD_J=F` and a clean
  obs dir, so no `ES_J_CONV` survives) - each run cold-starts on
  sweep 1 by the same rule
- identical grids enforced by shared `MIN_CF=2.4D-3` +
  `OBS_LAM_MAX=908117`: the P4 coverage guard requires the CMF grid
  below `NU_MIN_OBS*X_MIN_RAY` while the default `NU_MIN_OBS` tracks
  `MIN_CF` (see the implementation plan) - the pin reproduces the
  reference OBS_FREQ grid and decouples the requirement
- analysis tools `unit_tests/python/es_sweep_compare.py` (spectral
  pairs) and `unit_tests/python/es_resid_collect.py` (diag tree):
  files, verifies identical wavelength grids, masks non-finite,
  computes D_F/F_SCALE residuals on consecutive sweep counts plus
  continuum-window, line-profile (core/wings/EW/centroid) and
  integrated differences; emits CSV + PNG per spec section 6.

## 5. Known limitations (pre-registered)

1. `RJ_ES` file generations are compared; `RJ_CMF_ST` is the render's
   copy - identical by construction and verified, not assumed.
2. Phase-2 (DO_SOB_FLUX_OBS) sweeps are coherent and carry no RJ_ES
   generation; the diagnostic is phase-1 only by design.
3. Under KNOP, `Delta_JES` does not capture `I_NEXT` lag convergence;
   the existing `GET_FG_GS_MATURITY` census is the companion, and the
   audit reports both.
4. The file read at `:2400` costs one `ES_J_CONV` scan per sweep
   (ND*NCF reads) - negligible vs the sweep's formal solve; disabled
   cost is one logical test.
5. A stale `ES_J_CONV` from an earlier batobs run is not possible:
   `obs_template/batobs.sh` deletes generated files at setup; the
   seed read therefore only ever sees the current run's OLDJ output.

## 6. Results

Experiment status: the NUM_ES sweep matrix on
`m15_2p82d_SOB_l14_lte_adp2_tauto/cmfgen/2.82days` is configured
(SOB 1,2,3,4,6,8,12; KNOP 2,3,4,6,8,12; diag-off controls at NUM_ES=2)
and queued. The coverage configuration required for the current
executable is in `CMF_FLUX_P8_2_IMPLEMENTATION_PLAN.md` (the
`MIN_CF=2.4D-3` + `OBS_LAM_MAX=908117` pair that satisfies the P4
red-edge guard while reproducing the reference OBS_FREQ grid).

### 6.1 Residual histories

First dataset: `SN_A_11dh_p2p85r200_mc170_H10_SOB_d32/cmfgen/32.00days`
(obs_*_esN dirs; sob N=1,2,3,4,6,8,12 complete; knop N=2,3,4 complete,
6,8,12 in progress). `OLD_J=F`, cold start, `ES_RESID_DIAG=T`.

SOB pathway, `GLOBAL_MAX_R` per completed generation (ND=100,
NCF=4650; identical across runs of different NUM_ES - deterministic):

| gen | R_GLOBAL | argmax depth |
|---|---|---|
| 2 | 2.29749E-03 | K=100 |
| 3 | 9.72360E-04 | K=100 |
| 4 | 4.56571E-04 | K=88 |
| 5 | 2.59774E-04 | - |
| 6 | 1.53856E-04 | - |
| 7 | 9.40516E-05 | - |
| 8 | 6.00725E-05 | - |
| 9 | 3.90902E-05 | - |
| 10 | 2.57399E-05 | - |
| 11 | 1.71605E-05 | - |
| 12 | 1.13668E-05 | - |

Monotone geometric contraction, ratio ~0.55-0.65 per sweep, no
stall or oscillation through gen 12. Slowest depth sits at the
inner boundary (K=88..100): the thermalization floor dominates
the residual, as expected for a scattering-dominated field.
Zero non-finite and zero negative nodes in every sweep.

KNOP pathway (ND=100, NCF=102254 - 22x denser frequency
coverage than the SOB grid because the non-mono formal solve
needs line-resolution frequencies):

| gen | R_GLOBAL | argmax depth | R_REL max |
|---|---|---|---|
| 2 | 2.77769E-02 | K=88 | 1.96E-01 |
| 3 | 6.18909E-03 | K=71 | 4.23E-02 |
| 4 | 1.97009E-03 | K=71 | 1.67E-02 |

Faster per-sweep contraction (~0.22-0.32) but from a 12x larger
first residual; by gen 4 the absolute R is still ~4x the SOB
value at equal count. TOTAL_ACTIVE ~9.84M nodes per sweep.
Note (section 1.5): KNOP J additionally depends on I_NEXT lag
maturity, so Delta_JES understates the true sweep-to-sweep
change - the maturity census is the companion diagnostic.

[pending m15: `p82_residual_history.csv` per-run histories]

### 6.2 Final-generation render verification

Every completed SN_A run writes:

    RENDER GEN_LAST=N FIELD_MATCH=T MAX_ABS_DIFF=0.00000E+00
    FINAL CONFIGURED_ES_SWEEPS=N COMPLETED_ES_SWEEPS=N
        FINAL_REFRESH_GENERATION=N FINAL_RENDER_INPUT_GENERATION=N
        RENDER_MATCH=T INVALID_SWEEPS=0
        TOTAL_NONFINITE=0 TOTAL_NEGATIVE=0

Verified for sob N=1..12 and knop N=2,3,4: the render consumes
the last completed RJ_ES generation bit-for-bit. The
manufactured stale-render and cold-start cases are covered by
unit tests 15-17 in `test_es_resid.py`.

[pending m15: `p82_final_audit.csv`]

### 6.3 NUM_ES vs observer-frame spectrum

SN_A first dataset (`es_sweep_compare.py`, consecutive-pair
residuals on the 33865-point grid, optical continuum windows
4100-4300, 5300-5600, 6100-6300, 7000-7400 AA):

SOB pathway - `obs_fin_10` is NUM_ES-invariant:

    obs_fin_10 (sob):  BIT-IDENTICAL for N=1,2,3,4,6,8,12

    Mechanism (source-verified): the phase-2 render
    (OBS_FRAME_SUB_V9 at :3146 -> OBSFRAME -> obs_fin_10)
    runs with COHERENT_ES forced T (:1605) and a recomputed
    EDDFACTOR seed, so the phase-1 incoherent RJ_ES history
    never reaches the line-bearing render. Under coherent
    ES the outer sweep re-solves an already-converged
    system - extra sweeps are no-ops.

    The ES-sensitive SOB spectral product is OBSFLUX
    (COMP_OBS_V3 phase-1 continuum+BLANK render, same grid):

    pair    GLOBAL_MAX_R    CONT_REL 5300-5600   LINE_DEW Hb
    1->2    2.935e-04       1.529e-04            -2.13e-02 AA
    2->3    2.201e-04       7.645e-05            -8.96e-03 AA
    3->4    1.468e-04       7.819e-05            -
    4->6    1.468e-04       -                    -
    6->8    7.338e-05       -                    -
    8->12   7.338e-05       -                    -

    The OBSFLUX residual tracks the field history through
    ~gen 4 then flattens near 7e-5 - the printed-format floor
    (~4 significant digits in WRITV_V2 output), not a physics
    plateau; the field itself is still contracting at gen 12.

KNOP pathway - `obs_fin_10` responds to NUM_ES (single-phase,
incoherent fold reaches the emergent spectrum):

    pair    GLOBAL_MAX_R    CONT_REL 5300-5600   max-rel (active)
    2->3    1.900e-03       1.900e-03            8.93e-03
    3->4    8.249e-04       -                    -

    Global argmax sits at the extreme-UV grid end (weak nodes);
    the optical-window and active-mask numbers are the
    meaningful metrics.

[pending m15 + knop N=6,8,12: `p82_spectrum_metrics.csv`]

### 6.4 Line-profile convergence

SN_A knop 2->3 (representative windows, es_sweep_compare.py):
Hb core 1.58e-3, dEW -0.041 AA, dcentroid 0.031 AA;
HeI 5876 core 1.15e-3, dEW -0.075 AA, dcentroid 0.026 AA;
Ha core 2.02e-3, dEW +0.231 AA. Line metrics shrink roughly
in step with the flux residual (3->4 values pending; the
spectrum is still measurably changing at gen 4).

[pending: `p82_line_metrics.csv` per-window table]

### 6.5 Runtime and diagnostic overhead

[pending: `p82_runtimes.csv` - wall time vs NUM_ES; NUM_ES=2
diag-on vs diag-off pairs on m15]

### 6.6 SOB versus KNOP

SN_A first dataset, like-for-like generation comparison:

| | SOB | KNOP |
|---|---|---|
| frequency coverage NCF | 4650 | 102254 |
| gen-2 R_GLOBAL | 2.30e-3 | 2.78e-2 |
| gen-4 R_GLOBAL | 4.57e-4 | 1.97e-3 |
| contraction/sweep | ~0.55 | ~0.25 |
| emergent spectrum ES-sensitive | only via phase-1 OBSFLUX; obs_fin_10 invariant (coherent phase 2) | yes, obs_fin_10 direct fold |
| render-gen verified | MATCH all N | MATCH all N |

Equal NUM_ES does not imply equivalent convergence (section
1.5): the fields differ in grid density, the knop J also
carries I_NEXT lag, and the emergent-spectrum sensitivity is
structurally different between the two renders.

## 7. Conclusions

[pending - m15 sweep completion and knop high-N data]
