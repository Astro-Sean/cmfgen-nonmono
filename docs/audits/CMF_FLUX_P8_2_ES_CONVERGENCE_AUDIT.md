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
`m15_2p82d_SOB_l14_lte_adp2_tauto/cmfgen/2.82days` is complete
for SOB (N=1,2,3,4,6,8,12 plus a NUM_ES=2 diag-off control) and
in progress for KNOP (N=2..12 at 4 CPU - KNOP needs more than
the interactive node's 30G at 8 CPU; see 6.6). The coverage
configuration required for the current executable is in
`CMF_FLUX_P8_2_IMPLEMENTATION_PLAN.md` (the `MIN_CF=2.4D-3` +
`OBS_LAM_MAX=908117` pair that satisfies the P4 red-edge guard
while reproducing the reference OBS_FREQ grid).

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

| gen | R_GLOBAL | argmax depth |
|---|---|---|
| 2 | 2.77769E-02 | K=88 |
| 3 | 6.18909E-03 | K=71 |
| 4 | 1.97009E-03 | K=71 |
| 5 | 7.70972E-04 | K=71 |
| 6 | 3.51556E-04 | K=12 |
| 7 | 1.59749E-04 | K=12 |
| 8 | 7.39900E-05 | K=10 |
| 9 | 3.46684E-05 | K=7 |
| 10 | 1.62949E-05 | K=4 |
| 11 | 7.70055E-06 | K=1 |
| 12 | 3.69700E-06 | K=1 |

Monotone geometric contraction ~0.45/sweep through gen 12.
Unlike SOB (worst depth pinned at the inner boundary), the
KNOP argmax depth marches outward as deep layers converge
first - K=88 -> 71 -> 12 -> 7 -> 4 -> 1 - leaving the outer
atmosphere as the residual's last refuge. Zero non-finite and
zero negative nodes. TOTAL_ACTIVE ~9.84M nodes per sweep.
Note (section 1.5): KNOP J additionally depends on I_NEXT lag
maturity, so Delta_JES understates the true sweep-to-sweep
change - the maturity census is the companion diagnostic.

Second dataset: `m15_2p82d_SOB_l14_lte_adp2_tauto/cmfgen/2.82days`
(obs_p82_* dirs; sob N=1,2,3,4,6,8,12 complete). `OLD_J=F`, cold
start, `ES_RESID_DIAG=T`, `MIN_CF=2.4D-3` + `OBS_LAM_MAX=908117`.

SOB pathway, `GLOBAL_MAX_R` per completed generation (ND=70,
NCF=171207; identical across runs of different NUM_ES - the
redistribution is deterministic):

| gen | R_GLOBAL | argmax depth |
|---|---|---|
| 2 | 1.97678E-01 | K=11 |
| 3 | 9.09502E-02 | K=15 |
| 4 | 5.55618E-02 | K=15 |
| 5 | 3.79413E-02 | K=12 |
| 6 | 2.75945E-02 | K=12 |
| 7 | 2.09139E-02 | K=12 |
| 8 | 1.64012E-02 | K=12 |
| 9 | 1.32295E-02 | K=12 |
| 10 | 1.12186E-02 | K=52 |
| 11 | 1.03017E-02 | K=46 |
| 12 | 9.55676E-03 | K=46 |

Monotone decrease but the contraction ratio itself decays -
0.46, 0.61, 0.68, 0.73, 0.76, 0.78, 0.81, 0.85, 0.92, 0.93 -
i.e. the residual is shrinking slower than a clean geometric
sequence and is still ~1% of the local field scale at gen 12.
The argmax depth migrates (K=11 -> 15 -> 12 -> 52 -> 46): the
residual's last refuge moves from the inner thermalization
layers to mid-atmosphere. Zero non-finite, zero negative nodes
in every sweep; TOTAL_ACTIVE ~11.8M nodes per sweep.
Unlike the SN_A dataset (ratio ~0.55-0.65 held to gen 12), this
model's per-sweep gain is still decaying at N=12 - production
NUM_ES=2 leaves a ~20% first-refresh delta un-iterated.

[pending m15 KNOP: see 6.6]

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

m15 second dataset: verified for sob N=1,2,3,4,6,8,12 - the
render consumes the last completed RJ_ES generation bit-for-bit
(RENDER_MATCH=T, MAX_ABS_DIFF=0 in every run;
`p82_final_audit.csv`).

Diagnostic-off control: `obs_p82_sob_n2_diagoff` (ES_RESID_DIAG=F)
produced no ES_RESID_DIAG file and its obs_fin_10 is
bit-identical to the diag-on `obs_p82_sob_n2` run (max|dF|=0) -
the diagnostic is observational only.

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

    pair    GLOBAL_MAX_R    argmax lam (AA)
    2->3    1.900e-03       5524  (He I 5876 region)
    3->4    8.249e-04       6516  (Ha wing)
    4->6    5.079e-04
    6->8    1.270e-04
    8->12   6.350e-05       (printed-format floor)

    Argmax lands on real line features, not UV weak nodes -
    line cores and wings move measurably between low sweeps.

m15 second dataset (GLOBAL_LINE=BLANK, so the incoherent ES
fold reaches obs_fin_10 on the SOB pathway as well - unlike the
SN_A coherent-phase-2 render, obs_fin_10 IS ES-sensitive here).
Consecutive-pair residuals on the shared 158443-point grid:

    pair    GLOBAL_MAX_R    GLOBAL_MAX_REL   argmax lam (AA)
    1->2    6.540e-02       1.879e-01        1241
    2->3    2.871e-02       1.053e-01        -
    3->4    1.860e-02       6.852e-02        -
    4->6    2.310e-02       8.563e-02        -
    6->8    1.399e-02       5.592e-02        -
    8->12   1.621e-02       7.019e-02        -

The spectrum-level delta shrinks overall but not monotonically
at the tail (4->6 and 8->12 uptick over their predecessors):
the emergent spectrum is still drifting ~1.5% between N=8 and
N=12, consistent with the field residual still being ~1% at
gen 12. INTEGRATED_DF stays in the 6e7..1e8 Jy*AA band.
H-alpha metrics (LINE_DEW -7..-10 AA per pair, LINE_CORE
~1.2-3.7e-2) confirm the line region keeps moving.

[pending m15 KNOP: `p82_spectrum_metrics.csv` knop rows]

### 6.4 Line-profile convergence

SN_A knop, es_sweep_compare.py representative windows:

    pair    metric            Hb        HeI5876    Ha        CaII
    2->3    LINE_CORE         1.58e-3   1.15e-3    2.02e-3   -
            LINE_DEW (AA)     -0.041    -0.075     +0.231    +0.373
            LINE_DCENT (AA)   0.031     0.026      -         -
    3->4    LINE_DEW (AA)     +3.6e-4   -0.027     +0.168    +0.047

EW shifts of tenths of an Angstrom between 2 and 3 sweeps;
line centres move ~0.03 AA. The line metrics shrink in step
with the flux residual.

m15 second dataset, SOB pathway (BLANK render - obs_fin_10 is
ES-sensitive; representative windows):

    pair    metric            Hb        HeI5876    Ha        CaII
    1->2    LINE_CORE         4.83e-2   1.18e-1    3.70e-2   2.24e-3
            LINE_DEW (AA)     +2.86     -7.51      -10.21    -0.31
            LINE_DCENT (AA)   +2.09     +3.95      +113.5    +0.56
    2->3    LINE_CORE         3.05e-2   -          2.06e-2   -
            LINE_DEW (AA)     +1.45     -          -8.75     -

EW shifts are Angstrom-scale on m15 (this model's lines are far
stronger than SN_A's); Ha EW moves ~10 AA between N=1 and N=2
and still ~8 AA between N=8 and N=12 - the line-bearing regions
track the still-non-converged field. LINE_DCENT=nan on some
pairs marks windows where the computed centroid denominator
vanishes - reported, not hidden.

### 6.5 Runtime and diagnostic overhead

m15 batobs wall times (shared interactive node, 8 CPU SOB /
4 CPU KNOP; `p82_runtimes.csv`):

    sob  N=1  273s   N=2  237s   N=3  245s   N=4  375s
         N=6  468s   N=8  850s   N=12 770s

Diagnostic overhead is within run-to-run noise on this shared
node: NUM_ES=2 diag-on 237s vs diag-off 259s (the diag-off
control ran slower, so the instrumentation cost is not
measurable at this granularity - per-update file IO is a few
hundred formatted writes total).

Runtime scales roughly linearly with NUM_ES (each sweep runs
the full convolve + EDDFACTOR rewrite); KNOP sweeps are
slower than SOB sweeps at equal N because the non-mono formal
solve precedes the refresh (see 6.6 note on 30G memory).

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

m15 second dataset (KNOP gens verified through n12 run in
progress; deterministic across runs - n2/n3/n12 share gen-1..3
values bit-for-bit):

| | SOB | KNOP |
|---|---|---|
| frequency coverage NCF | 171207 | 182712 |
| gen-2 R_GLOBAL | 1.98e-1 | 5.78e-1 |
| gen-3 R_GLOBAL | 9.10e-2 | 2.65e-1 |
| gen-4 R_GLOBAL | 5.56e-2 | 1.78e-1 |
| contraction 2->3->4 | 0.46 / 0.61 | 0.46 / 0.67 |
| obs_fin_10 ES-sensitive | yes (BLANK fold) | yes |

KNOP's first refresh moves J by ~58% - nearly 3x SOB's first
delta on the same model - and the per-depth rows show the
change is broad (R ~0.2-0.58 at most depths), not localized to
a boundary band. KNOP sweeps cost ~10min each at 8 CPU
(~15min at 4 CPU): the interactive node's 30G memory is
sufficient only at 4 CPU (8-CPU runs peaked at 31.2G and were
OOM-killed; the convolved-field workspace scales with threads
here).

## 7. BLANK render on non-monotonic populations

`obs_line_treatment: blank` in the pipeline yaml produces an
OBS-stage run with `GLOBAL_LINE=BLANK`, all `DO_NONMONO_*`
flags F: the populations are loaded from the converged model,
but the OBS solve runs the stock `FG_J_CMF_V12` +
`MOM_J_CMF_V8` + `N_ON_J` f-convergence machinery
(`obs/COMP_JCONT_V4.INC`). None of the fork's
non-monotonic alterations in that file apply to this branch.

### 7.1 Failure on SN_A mc170_H10

All five attempted configurations stop at the same guard,
`ABS(RJ(1)) .GT. 1.0E+30`, inside the per-frequency f
convergence loop, at nu ~ 0.053e15 Hz (lambda ~ 5.6 um,
far-IR continuum, first inner-boundary thermalisation zone):

| run | overrides | result |
|---|---|---|
| obs_blank | stock (N_ON_J, DELV_MOM=0.5, OLD_J=F) | FAIL ~2 min |
| obs_blank_dm | DELV_MOM=1.0 | FAIL ~2 min |
| obs_blank_gonly | N_TYPE=G_ONLY + DELV_MOM=1.0 | FAIL ~2 min |
| obs_blank_full | N_TYPE=MIXED + DELV_MOM=1.0 + INC_GRID=T | FAIL ~2 min |
| obs_blank_oldj | OLD_J=T (warm RJ_ES seed from model EDDFACTOR) | FAIL ~2 min |

The knobs suggested by the stock error message (N_TYPE,
DELV_MOM/DELV_FG, INC_GRID+NINS) do not rescue the run.

Per project directive, `N_TYPE=N_ON_J` is the fixed choice
for BLANK renders - it is the acceptable variant because the
alternatives (G_ONLY, MIXED) engage velocity-switch handling
this model does not need. The G_ONLY/MIXED rows above were
probes of the stock suggestion list only; no further N_TYPE
experiments are planned.

### 7.2 Counter-example: m15 SDEC render

`cmfgen_validation/m15_2p82d_SOB_l14_lte_adp2_tauto/cmfgen/
2.82days/SDEC/Cr` is a genuine plain-BLANK OBS run -
`DO_NONMONO_KNOP=F`, `DO_NONMONO_SOB=F`, `DO_SOB_FLUX_OBS=F`,
`GLOBAL_LINE=BLANK`, full frequency range (MIN_CF=3e-3 to
MAX_CF=200, 331778 continuum frequencies) - that completed
"CMF_FLUX has finished" on the same Sobolev-relaxed m15
population files. BLANK-on-SOB populations is therefore not
categorically blocked: the mc170_H10 divergence is a
model-specific numerical instability of the J/f iteration at
the far-IR boundary, not a structural rejection of the
population provenance.

### 7.3 USE_FIXED_J mechanism

`COMP_JCONT_V4.INC:132` carries a stock `USE_FIXED_J` branch:
under `CONT_VEL` (auto-true for SECTION='CONTINUUM') it calls
`RD_CONT_J_OBS`, which walks EDDFACTOR records and
interpolates the stored continuum J at the current frequency,
then renders without a J/f solve. This is the closest stock
mechanism to "BLANK transfer evaluated on populations/J from
a cheaper line treatment" - the same division of labour as
the deposited-SOB render, with BLANK chi/eta structure on a
stored field.

Constraints found while testing it:

- It needs an EDDFACTOR on a compatible grid. The model's own
  EDDFACTOR fails immediately: `COMPUTE_EDDFAC=F` then reads
  records index-aligned and checks `ABS(T1-FL)/FL .LE. 1e-14`
  (`COMP_JCONT_V4.INC:719-725`); the model grid (~3700
  records) differs from the OBS continuum grid (4650), so the
  run stops at ML=2 with "incorrect reading of the mean
  intensity" (job obs_blank_fixedj, 6 s).
- An OBS-grid EDDFACTOR from a completed run on the same
  model (obs_sob_es1) is NOT usable either: the sob and blank
  runs build different continuum grids (obs_blank_fixedj2
  stopped at FL=200.0 against a zero tag - the grids differ
  from the first record).  USE_FIXED_J seeding therefore has
  no viable source file on this model.

### 7.4 Source adaptation: DAMP_FJ and J_RESCUE

Because no parameter or stored-field path reaches a BLANK
render on mc170_H10, two guarded mechanisms were added to
`COMP_JCONT_V4.INC` (ND-path f loop; the ACCURATE/INC_GRID
branch is unchanged):

- `DAMP_FJ` (real, default 1.0): under-relaxes the f/J
  iteration - `FEDD <- FOLD + w*(FEDD-FOLD)` and
  `RJ <- RJ_old + w*(RJ-RJ_old)` - while the convergence test
  still measures the undamped residual.  At the default the
  update is bit-identical to stock.
- `J_RESCUE` (logical, default F): if the iteration exits
  with non-finite or exploding boundary J (>1e25 inside the
  loop triggers early exit; >1e30 or NaN at loop exit
  triggers the restore), the previous frequency's converged
  J/f/BC fields are restored, a warning is logged, and the
  run continues.  The EDDFACTOR record stores the fallback J,
  so later frequencies warm-start from it.  No effect at the
  first frequency - no previous field exists there, and the
  stock guard still stops the run.

Project directive: a non-monotonic model's OBS stage must
always return a blanketed spectrum.  These controls are the
mechanism; N_TYPE stays N_ON_J throughout (no velocity-switch
variants).

Validation runs on mc170_H10 (v1 implementation, RJ(1)-based
detection):

| run | overrides | result |
|---|---|---|
| obs_blank_damp | DAMP_FJ=0.3 + J_RESCUE=T | completed; ~32k rescues |
| obs_blank_resc | J_RESCUE=T (no damping) | completed; ~35k rescues |

DAMP_FJ=0.3 did not converge the infrared wall either - the
failure is not iteration overshoot.

### 7.5 Root cause: JNU_PREV poison chain, not a single stiff frequency

Reading the resc-run EDDFACTOR and ES_J_CONV showed what the
log had hidden:

- The first divergent frequency (FL=0.0530) exited the f loop
  via the L>20 path with J already at 9.7e24..6.4e53 across
  depth - but `ABS(RJ(1))` stayed under the 1e25 in-loop
  trigger and the 1e30 post-loop guard, so no rescue fired.
  The garbage J was written to EDDFACTOR.
- `JNU_PREV` is saved at the top of every frequency iteration
  from the previous frequency's stored RJ.  Once record
  0.0530 stored garbage, every later frequency "rescued" back
  to the same garbage - a self-perpetuating poison chain.
- The convolve pass then carried the garbage into RJ_ES
  (ES_J_CONV holds the identical 9.7e24..6.4e53 at IR
  frequencies), and the incoherent emissivity
  `ESEC*RJ_ES` injected it into sweep 2.

So the enormous RJ_ES is downstream contamination; the first
divergence is coherent-stage (sweep 1, `COH_ES=F`, no ESEC
term at all).

### 7.6 v2 J_RESCUE: field-wide detection + RJ_ES sanitize

Corrections to the rescue machinery (same file, same flags):

- Detection now uses `MAXVAL(ABS(RJ(1:ND)))` plus non-finite
  checks at all depths, for both the in-loop early exit and
  the post-loop restore.
- The bound is per-frequency:
  `RESC_CAP = MIN(1e10, MAX(1e4*max|JNU_PREV|, 1e-3))` -
  a 4-order jump between adjacent continuum frequencies is
  unambiguous runaway, and the relative test catches fields
  far below any absolute ceiling.  FIRST_FREQ uses the 1e10
  backstop alone (no previous field).
- New sanitize: when `.NOT. COHERENT_ES`, the read-back
  RJ_ES is tested against the same bound; a non-physical
  convolved field falls back to the local coherent J for
  that frequency only, so `ESEC*RJ_ES` stays sane in sweep 2.
- Message throttling (first 20 + every 500th) and a
  first-divergence depth dump: scalars (DBB, DIF, HBC_CMF,
  NBC_CMF, dLOG_NU) plus per-depth R, V, ED, CHI, ETA, ESEC,
  CHI_SCAT, ZETA, RJ, FEDD, JNU_PREV.

Expected effect: even where the underlying IR solve still
diverges, the stored J stays bounded and sane, the warm
start is unpoisoned, and sweep-2 emissivities are physical.
The first-divergence dump records which input field drives
the runaway.  obs_blank_resc2 (J_RESCUE=T, DAMP_FJ=1.0,
job 13157399) is the validation run.

### 7.7 First-divergence anatomy and the ratchet flaw

The NRESC=1 dump at FL=0.0532 localized the runaway:
J is sane over K=1..96, then spikes at the inner boundary
(K=97..100: 6.6e-6 -> 1.1e-4 -> 2.8e-3 -> 0.24).  FEDD
oscillates 0.35..0.84 non-monotonically in those cells;
DBB=-2.0e-12 (DTDR<0 from the cool inner zone);
ESEC~2e-5 dominates CHI~2e-6 - dense, pure-scattering
inner layers.  Extending the existing nonmono luminosity-
derived DTDR and the DBB>=0 floor to J_RESCUE made DBB
positive (+1.5e-11) but the runaway was identical
(obs_blank_resc3, job 13157606) - the sign of the
boundary coefficient is not the driver.

Reading resc2's stored EDDFACTOR exposed a subtler defect
in the v2 cap: Jmax vs frequency is sane outside
0.025..0.048 but inside the band ramps smoothly from ~1e-4
to ~3e8 with stored FEDD of 1e4..1e6.  A converged moment
solution cannot have f>1, so the band is a slow ratchet:
each frequency stores slightly-larger under-cap garbage,
JNU_PREV absorbs it, and the next iteration grows from
there.  RESC_CAP = 1e4*prev was too loose to catch it.

Two candidate gates were tested:

- v3 gated restore on a new F_DID_NOT_CONV flag (any L>20
  or early-exit).  obs_blank_resc4 (job 13157729) showed
  this is wrong on this model: stiff frequencies routinely
  exit L>20 with residuals ~1e-3 and a sane field - stock
  semantics keep them, sob_es1 has 153 such exits and a
  bounded spectrum.  Restoring on convergence status
  froze J at the first frequency's field; the run was
  cancelled.
- v4 gates on field validity instead: non-finite J or f,
  |RJ|max > RESC_CAP, or |FEDD|max > 2.  f <= 1 is
  required for a real radiation field, so f > 2 at any
  depth marks an invalid iterate; the ratchet fields had
  f ~ 1e4-1e6 while every converged field has f <= 1.
  Mild L>20 exits keep their fields (stock behaviour), so
  the warm-start chain is unbroken outside genuinely
  diverged bands.  obs_blank_resc5 (job 13157778)
  carries v4.

v4 result: obs_blank_resc5 produced an IR tail identical
to resc2 (max I ~ 6e18, 4109 points >1e10).  The FEDD>2
test fired nowhere - a smooth multiplicative ratchet
keeps f in bounds because f is the ratio K/J: if J and K
scale together the shape looks converged while the
amplitude drifts.  Deeper reading showed why no
iteration-level test can catch it: in the pure-scattering
IR regime (ETA ~ 1e-21 << scattering), the formal-solve
source TA = ETA + CHI_SCAT*RJ is nearly homogeneous in J,
so J = A.J + b.ETA has a continuum of near-fixed-points -
any amplitude satisfies the iteration.  The loop's only
convergence test is on FEDD, which is amplitude-blind.
KNOP's coupled moment solver does not share this defect
and converges the same frequencies cleanly: its stored
EDDFACTOR field is flat J ~ 1e-7-1e-5 across the whole
band, which is the reference answer.

The remaining defect was in RESC_CAP itself:
MAX(1e4*prev, 1e-3) is a floor - once JNU_PREV inflated
past ~1e-7 the relative term recalibrated the cap upward
and the ratchet walked underneath it.  v5 replaces the
floor with the smaller of two bounds, split by role:

- RESC_CAP (post-loop restore + RJ_ES sanitize):
  MIN(1e4*JNU_PREVmax, RESC_JMAX).  The relative bound
  catches single-step explosions and dark-region garbage;
  RESC_JMAX (new control, default 1e-3) is an absolute
  ceiling - legitimate J peaks at 4.8e-5 across the whole
  converged grid (KNOP and SOB agree), and the largest
  legitimate adjacent-frequency Jmax ratio measured is
  ~44x, so the 1e4x bound has ~200x margin and cannot
  false-fire on real transitions.
- RESC_CAP_IT (in-loop early exit only): the previous
  loose bound MIN(1e10, MAX(1e4*prev,1e-3)).  Iterates
  may transiently overshoot a tight cap, so the
  early-exit stays on runaway scale - it trims wasted
  iterations, it does not decide validity.

obs_blank_resc6 (job 13158209) carries v5.

v5 result on mc170_H10 - the objective is met.  The
BLANK render completes and returns a bounded blanketed
spectrum:

- Rescues fire only inside the eigenproblem band
  FL~0.031-0.053 (every diverged frequency restores the
  last converged field); frequencies below 0.031 converge
  on their own, so the deep tail is real, not fallback.
- obs_fin_10: max I = 2.3e6 Jy, zero points >1e10 - vs
  resc5's 5.9e18 and the converged KNOP reference at
  1.6e6.  Deep-IR edge 3.4e4 vs KNOP 2.6e4.
- Ratio to the fully-converged KNOP spectrum: median 1.00
  over FL 0.5-5 (optical/UV, exact), 1.36 below the band,
  ~4.3 median inside the restored band (max ~70x).  The
  band overshoot is the flat last-converged field sitting
  slightly above KNOP's local J - bounded, flagged, and
  ~10 orders below the previous garbage.
- 57 mild L>20 warnings kept (stock semantics, same as
  sob_es1's 153-exit pattern); no RJ_ES sanitize events.

The remaining ~4x band overshoot is acceptable for the
science goal (the band sits at lambda >~5.6 um, outside
every photometric band used in the light-curve fit), and
it is honest: every rescued frequency is logged.

Verification on the v5 objects: P8.1 battery 116/116,
P8.2 battery 20/20 - J_RESCUE=F default leaves the stock
path untouched.

Known limitation: all J_RESCUE machinery lives in the ND
(non-INC_GRID) branch of COMP_JCONT_V4.INC, which is the
path every production render on these models takes.  The
INC_GRID branch keeps stock behaviour under J_RESCUE=T.

## 8. Conclusions

[pending - m15 sweep completion and knop high-N data;
BLANK bounded-rescue outcome in section 7.7]
