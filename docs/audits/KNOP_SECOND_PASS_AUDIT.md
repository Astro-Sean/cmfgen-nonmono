# KNOP second-pass mathematical and physical audit - 11-Nov-2026

Scope: `DO_NONMONO_KNOP` formal solver `plane/fg_j_cmf_nonmono_v1.f`
(FG), population coupling `new_main/mod_subs/solveba_v13.f` /
`solve_for_pops.f`, observer path, and the lagged `I_NEXT` field
machinery.  This is a verification audit of the convergence-failure
mechanisms reported in `KNOP_PRODUCTION_AUDIT.md`; it changes no solver
equation and no production pathway.

Primary runtime evidence: `m15_2p82d_FRESH_knop` (STOP it124,
`GS_REL_NOPROG_LIM`), `m15_2p82d_FRESH_knop_d1a` (running, same
freeze/release signature), `m15_2p82d_VAL_knop_formal` (STOP it1489,
luminosity divergence).  All line numbers verified against the working
tree at the audit date.

Classification used below: **[DEFECT]** confirmed coding defect,
**[MATH]** confirmed mathematical inconsistency, **[APPROX]** physical
approximation / deliberate truncation, **[HYP]** unverified hypothesis,
**[MISSING]** missing runtime evidence.

---

## Pillar 1 - near-zero-GAM convergence interpretation

### The reviewer's limiting argument is correct

For a compressing cell the lagged-field source is

```
S' = S + |Q| * I_NEXT,     |Q| = |GAM| / (DLN * chi_hat)
chi_hat = CHI + |GAM| / DLN
```

so, for finite positive `CHI` and `DLN`,

```
|Q| = |GAM| / (CHI*DLN + |GAM|) = x/(1+x),   x = |GAM|/(CHI*DLN)
```

As `GAM -> 0` (i.e. `|1 + SIGMA*MU^2| -> 0`), `x -> 0`, `|Q| -> 0`,
the `|GAM|/DLN` term in `chi_hat` vanishes, and both sides of the
lagged coupling disappear.  A perturbation `dI_NEXT` produces
`dI_solved ~ |Q| * dI_NEXT -> 0`; the store update
`NEW = (1-omega)*OLD + omega*I_solved` then contracts at the
relaxation rate.  These cells are weakly coupled and contractive.
**[MATH - confirmed]**

### "Marginally contractive near-zero-GAM cells" is contradicted

The production audit's phrase is **wrong** and is retracted here.  The
marginal (`|Q| -> 1`) regime is the opposite limit:

```
|GAM| / (DLN * CHI) -> infinity   (strong compression / low opacity)
=>  |Q| -> 1^-
```

In the thick limit `I -> S' ~ |Q| * I_NEXT`, so `dI_new/dI_NEXT ->
|Q| -> 1`: gain approaches unity per frequency hop.  The solver's own
comments already name this mechanism - "marginally contractive
(|Q|->1) cells at strong compressions" (`fg_j_cmf_nonmono_v1.f`
~1021), "marginal |Q|~1 modes" (~10205), "near-unity |Q|/(chi+|Q|)
links" (~7343).  **[MATH - confirmed]**

### Runtime cell evidence (FRESH_knop, final iteration)

From the `KNOP FIELD` block and iteration ledger:

- Refresh residual `max = 3.0986E-02` vs `tol = 2.0E-03` (~15x).
- `unconverged cells = 138464 of 620909597` (`frac = 2.23E-04`).
- `matched/compress = 17174408 / 258714053` - **94% of ray cells are
  compressing** (`GAM < 0`), so the lagged `I_NEXT` coupling is the
  norm, not the exception.
- `argmax history (d,r,f)=v`: cells `(3,82,1249)`, `(5,80,1249)`,
  `(4,81,1249)`, `(9,76,1372)` all pinned at `v = 3.06E-02` across
  iterations - the *same* cells hold the residual, the signature of a
  non-contracting fixed point rather than a transient.
- `worst cell: i=109 l=1469`, `worst ray=12 ni=139`.
- `bad cells last/best pass = 138464 / 139173` - the best pass is
  *worse* than the last, i.e. the refresh is not healing the store.
- `release_no_progress: demand=3.843E+01 ... count=10` then
  `GS_REL_NOPROG_LIM`.

### Velocity-gradient structure

`NEW_SN_R_GRID` gives `SIGMA = d ln V / d ln R` on the model grid:
negative at **every** depth, from `-1.04` at the outer edge to `-157`
in the inner region (minimum `-155.9`).  A compressing/reversing
velocity field is global, so `|1 + SIGMA*MU^2| -> 0` is reachable for a
range of `MU` at essentially every depth, and the `|Q| -> 1` strong
coupling is widespread.  The argmax cells sit at depths 3-5; under the
`NEW_SN_R_GRID` row ordering these correspond to `SIGMA ~ -1.07` to
`-1.18` (mild compression), where `|Q| -> 1` requires
`|GAM| >> DLN*CHI` rather than `GAM -> 0`.

### Verdict on Pillar 1

The residual plateau is **not** caused by near-zero-`GAM` cells -
those decouple and contract.  It is caused by the `|Q| -> 1`
strong-coupling compressing cells that dominate the grid (94%
compressing) together with the outer population-field coupling
(`se_res ~ 1.65E+13`, concentrated in the worst columns; `se_res2 ~
3.8E+10` after excluding them).  Whether the discretised operator is
intrinsically non-contractive on these cells, or merely slow, is
**[HYP]** - it requires the finite-difference `I_NEXT` sensitivity
measurement described in `KNOP_SECOND_PASS_TODO.md` P1-T1.  What is
confirmed is that the *stated* mechanism (near-zero-`GAM`) is
mathematically incorrect.

---

## Pillar 2 - `KN_NO_GAIN` physical consistency

### What the flag modifies **[DEFECT-adjacent / MATH]**

`KN_NO_GAIN` is enforced in `new_main/cmfgen_sub.f` (~2605-2621 main
loop, ~6346-6362 final-spectrum pass).  For every depth where
`CHI(I) < CHI_SCAT(I)`:

```fortran
CHI(I) = CHI_SCAT(I)          ! net absorption clamped to zero
NEG_OPACITY(I) = .TRUE.
NEG_OPAC_FAC(I) = 0.0
```

`CHI < CHI_SCAT` is precisely a **population inversion** (negative net
absorption = maser gain).  The floor therefore:

- raises physical total extinction `CHI` to `CHI_SCAT` (albedo
  `THETA = CHI_SCAT/CHI` capped at 1);
- reduces `chi_hat = CHI + |GAM|/DLN`;
- reduces `ZETA = ETA/CHI`;
- caps the source-coefficient sum at 1;
- propagates into the captured `KN_CHI_EXT`/`PK_*` pack, so `Q` is
  normalised by floored opacities;
- sets `NEG_OPACITY`/`NEG_OPAC_FAC=0`, which suppresses line-opacity
  contributions in the BA assembly and crushes `CHIL_MAT` at flagged
  depths.

It does **not** touch `GAM`, `GAM_MAG`, `dLOG_NU`, `CHI_SCAT`, `ETA`,
or the pack arrays directly - but the floored `CHI` is what the solver
consumes.  **[MATH - confirmed]**

### Physical vs numerical gain **[APPROX]**

In floored cells the solved problem is **not** the original transfer
equation: `CHI` is replaced by `CHI_SCAT`, changing `chi_hat`, `ZETA`,
the source-coefficient sum and `Q`.  The flag removes **all** gain,
including physically real maser amplification.  The in-code rationale
("inversions sit on high-n Rydberg and intra-n far-IR levels ...
overpopulated by the recombination cascade on transient iterates",
~2597-2603) is a code-comment assertion, **not independently
verified**.  Whether those inversions are genuine physics or iterate
artifacts is **[HYP]**.

The floor is load-bearing: with `KN_NO_GAIN=F` the lagged store
recycled amplified intensity at a measured `8.6x/iter`, diverging
`L_out` to `~1E29` (documented at `rd_control_variables.f` ~3444-3452).
So the flag is a necessary stabiliser of a *deliberately modified*
operator, not a harmless clip.  **[APPROX]**

### Configuration and runtime census

- Under `DO_NONMONO_KNOP`, `KN_NO_GAIN` is **forced `.TRUE.`** with no
  escape hatch (`rd_control_variables.f` ~3443-3467); the only bypass
  is disabling KNOP.
- `FRESH_knop`: `floored 13538 cells (depths 1-56); max
  chi_scat/chi = 1.242E+01` at `nu=8.221E-03` (`lambda=3.647E+05 A`,
  far-IR) - a 12x opacity modification.
- `VAL_knop_formal`: `floored 8447 cells`, `max chi_scat/chi=5.1`.
- `NEG_OPAC_OPT` (`0.1*CHI_SCAT` floors) is a weaker pre-existing
  control applied earlier in the same loops; `KN_NO_GAIN` is the
  stronger KNOP-only floor at exactly `CHI_SCAT`.

**[MISSING]** No controlled `KN_NO_GAIN=T` vs `F` comparison of
emergent luminosity, continuum flux, energy balance or population
residuals exists, because `F` diverges.  Quantifying the physical cost
of the floor requires a run that survives without it (e.g. smaller
`INEXT_RELAX` or a gain-limited rather than gain-zeroing variant) -
see TODO P2-T1.

---

## Pillar 3 - ordered-march branch usage

### Dispatch **[confirmed]**

`FG_ORD_SWEEP` (default `.TRUE.`, set from VADAT `KN_RAY_ORDERED`)
selects `SWEEP_FLD_ORD` over `SWEEP_FLD` per ray, gated by
`ORD_SWEEP = FG_ORD_SWEEP .AND. HAVE_PACK .AND. NI .GE. 3`
(`fg_j_cmf_nonmono_v1.f` ~7591).  `SWEEP_FLD` alternates frequency
order by pass parity; `SWEEP_FLD_ORD` marches along the characteristic
and orders frequency by `sign(GAM)` (`GAM>=0` ascending, else
descending), using the `SPW` ring for the parabolic slope term.

For all three dead/live runs (`NI=70`, pack `64*NI*NCF` well under the
2 GB cap), the defaults imply `ORD_SWEEP=T` and `SWEEP_FLD_ORD` should
have executed for both `I_M` and `I_P` of every ray.

### Observability gap **[DEFECT - diagnostics]**

**No counter, ledger column or print anywhere records whether
`SWEEP_FLD` or `SWEEP_FLD_ORD` executed.**  The only `ORD_SWEEP` uses
are the dispatch lines; the `KNOP FIELD` block reports residual, pass
counts, exit census and AA/p2damp counters but nothing branch-specific.
Branch execution therefore **cannot** be confirmed from surviving
output - the earlier inference that ordered marching ran is plausible
but unproven.  **[MISSING]**

### A/B plan (TODO P3-T1)

Controlled `KN_RAY_ORDERED=T` vs `F` on an identical warm-started
checkpoint, holding model, atomic data, `I_NEXT` seed, refresh
tolerances, pass limits and OMP fixed.  Measure post-refresh residual,
per-ray convergence rate, pass count, failed rays, store drift and
wall-clock.  First add the `ORD_SWEEP`/`HAVE_PACK`/`NI` census line and
an ITER_LEDGER branch column so the comparison is auditable.

---

## Pillar 4 - exact transfer-equation derivation

### The `+3 a_nu I` term is missing by convention **[APPROX]**

The implemented equation is `dI/ds = eta - chi*I - a_nu*(dI/d(ln nu)
- 3*I)` in nu-form.  The proportional (`+3 a_nu I`) piece **does not
vanish**; the surviving implicit term is the `|GAM|/DLN` inside
`chi_hat`.  The detailed justification block (~3501-3537) states the
term is "deliberately not added" because it "would change results vs
stock V13" - i.e. an **inherited CMFGEN/V13 truncation**, not a Knop
term-by-term implementation.

**[MATH - confirmed inconsistency in comments]** The file header
(~89-101) claims the term "cancels identically ... so nothing is
dropped", which contradicts the detailed block.  One of the two
justifications is wrong; the detailed block is the defensible one and
the header should be corrected (documentation fix, no solver change).

A secondary sign inconsistency exists: the documented equation writes
`+ GAM*dI/d(ln nu)` (~3474) while another comment writes
`- GAM*dI/d(ln nu)` (~4023).  **[MATH]**

### Mapping to Knop and the monotonic limit **[APPROX / MATH]**

`chi_hat` corresponds to Knop Paper I Eq. 7 and the signed `Q` to Eq.
9; the frequency-coupling stencil and the `|GAM|/DLN` opacity increment
are the implemented pieces.  The `+3 a_nu I` term is absent.

**Monotonic limit verified**: for `GAM > 0` everywhere the KNOP
expanding branch reduces *exactly* to stock V13 - V13 computes
`Q=GAM/dLOG_NU; CHI_RAY+=Q; Q=Q/CHI_RAY; SOURCE_RAY=ETA/CHI_RAY`,
identical to the KNOP expanding branch, and V13 has no `S_tilde`
either.  The removed `S_tilde = 4*Q*dLOG_NU*I` was the wrong
(lambda-form) coefficient and double-counted the drift.

**Consequence**: the solver must not be described as term-by-term
equivalent to Knop Paper II - it omits a physical term by convention.
The impact of that omission on the non-monotonic solution is
**[HYP]** (bounded in the monotonic limit by construction, unquantified
in compression).  Quantification is TODO P4-T1.

---

## Pillar 5 - population runaway and null-space

### Mechanism **[MATH - confirmed on fast model, corroborated on production]**

Per depth, `GENERATE_FULL_MATRIX_V3` assembles the `NT x NT` SE block
then executes `C_MAT(I,J) *= POPS(J,L)` (`generate_full_matrix_v3.f`
~385-390), converting the Jacobian to derivatives w.r.t. the
*fractional* correction `dz_j/z_j`.  A level at `POPS ~ 1E-30` then
contributes a column of magnitude `~1E-30 x O(rates)` - an approximate
null column.

`CMF_BLKBAND_PIN_NULL_V1` pins only *exact* `0=0` rows (not approximate
nulls); `DGEEQU` sets `COL_SF(J) ~ 1/POPS(J) ~ 1E30-1E43`; `DGETRF`
succeeds (no exact-zero pivot) but the LU residual error is amplified
on unscale `x_j = COL_SF(J)*z_j`, yielding fractional-correction
demands of order `1E24-1E47%`.  Conservation rows do not restore the
direction because their entry in the floor column is also
`POPS`-scaled to `~1E-30`.

This is **approximate rank loss / ill-conditioning amplified by column
unscaling** - not an exact singularity and not a floor-clipping
artifact.  `BA_DEM_PFLOOR` exempts floor rows from the *demand metric*
only; the demand scans read raw `STEQ_RAW`, so masking cannot hide it.
The `1E7` in `maxraw` is the DECREASE sentinel, not the raw demand;
uncapped floor demand reaches `1E21-1E53` in `exm`.  **[MATH]**

### Evidence

- Fast m15 KNOP model (Item-1 census): `FeI:3d6_4s2_a5De[2]` @d63,
  `POPS~1E-31`, `rcond~1E-8`, `colmin_u~1E-43`, `BA_DEM_EXR~97%`.
- `VAL_knop_formal`: `exm` oscillates `1E21-1E53`, `maxraw` pinned at
  `1.0000E+07` from it1463, `maxapp=9.9E3` (MAX_LAM clip), `lum_ratio
  -> -3.5`, died it1489 on luminosity divergence - same signature.
- `FRESH_knop`: cold transient `maxraw=1E7, exm=1.8E17` (it1-8), then
  demand falls to `~38-500`, `exm->~20-100`, and it stalls in the GS
  freeze/release cycle.

### Important distinction **[MATH]**

`FRESH_knop`'s **terminal** stall is a *different* mode: the GS
freeze/release field-maturity limit cycle with moderate (~38%) demand,
not the `1E7` floor runaway.  H1 explains `VAL_knop_formal`'s death and
`FRESH_knop`'s cold-start transient, **not** `FRESH_knop`'s end state.
Two distinct kill mechanisms operate on this model.

### Missing runtime evidence **[MISSING]**

Neither dead run set `POP_AUDIT`, `POP_MASK_FLOOR`, `POP_TRACE_MAXFAC`
or `POP_RUNAWAY_*`; `DO_BA_QUAR=F` so quarantine is fully inert.
Per-depth `rcond`/`colmin_u`/`nflr` census therefore does not exist for
the production runs - the conditioning proof is inherited from the fast
model plus `ITER_HISTORY`.  `OUTGEN` is consumed by restart prep, so
the `RUNAWAY` channel lines do not survive in the run dirs.  Settling
H1 on the production model needs a `POP_AUDIT=T` warm-restart rerun
with OUTGEN preserved (TODO P5-T1).

Do **not** port `BA_DEM_PFLOOR`/`BA_EX_RELFRAC` semantics from SOB
until this is settled - they are demand-metric exemptions, not
conditioning fixes.

---

## Pillar 6 - observer validation

### D1-D4 status **[confirmed defects, fixes unapplied]**

The four observer defects from `CMF_FLUX_KNOP_CONVERGENCE_AUDIT.md`
stand: the observer render uses a **cold** `I_NEXT` store on the OBS
grid, which in compressing cells acts as an artificial sink; the
reported `~9.5x` continuum suppression is **not** a validated physical
SOB-vs-KNOP difference - it is at least partly the cold-store artefact.
Until a converged `I_NEXT` exists there is no meaningful KNOP
spectrum to validate.

### Controlled `NUM_ES` plan (TODO P6-T1)

Once a converged `I_NEXT` is available (currently blocked on Pillars
1/5), run `NUM_ES = 2, 3, 5` on identical populations and record the
post-refresh residual and integrated luminosity per case.  The final
observer calculation must use a sufficiently converged `I_NEXT`.
Adaptive refresh is **not** to be implemented until the existing
refresh diagnostics demonstrate it is required - the current evidence
(`passes<=12`, `hot-pass budget spent`, residual plateau) does not yet
isolate a refresh-count deficiency from the field non-contraction.

---

## Pillar 7 - documentation consistency

### Formal-solver audit vs production audit **[MATH - terminology]**

`KNOP_FORMAL_SOLVER_AUDIT.md` item 20 states the "end-to-end
integration test" is "satisfied by the model run itself"
(`DO_NONMONO_KNOP=T` through `run_model_grid.py`, run
`m15_2p82d_KNOP_w19`).  `KNOP_PRODUCTION_AUDIT.md` states no full
KNOP model has converged and "no full-iteration integration tests
exist" (F-23).  These are **reconcilable but were using the same
phrase for two different claims**:

- item 20 = *execution* integration - the VADAT->solver->output path
  runs to completion and produces files (true for w19);
- F-23 = *convergence-behaviour* integration - the coupled
  population/field iteration converges correctly (absent; no run has
  done so).

A model execution without convergence must not be classified as a
successful end-to-end *convergence* test.  The formal-solver audit's
item 19/20 wording should be read as execution-only; this audit
records the distinction explicitly.  **[MATH]**

### Future-dated source comments **[DEFECT - documentation, fixed]**

`plane/fg_j_cmf_nonmono_v1.f` carried `audit 20-Nov-2026` references at
the `KN_COUPLING_OK` gate (lines ~3747 and ~5483) and the refresh-side
mirror (~9814).  The field-drift audit fixes are dated **10-Nov-2026**
and the production audit **11-Nov-2026**; 20-Nov is in the future.  All
three comments now read `10-Nov-2026`.  Comment-only change; no solver
behaviour altered.  **[DEFECT - corrected]**

The header/detailed-block contradiction on the `+3 a_nu I` term
(Pillar 4) is recorded as an open documentation defect - the comment
fix is deferred to avoid touching solver-adjacent text during this
audit; it is listed in the TODO.

---

## Synthesis - what actually causes the KNOP convergence failures

The failures are **not** a gross coding defect in the transport solve.
The formal solver is term-by-term correct against Knop Paper I/II
(modulo the inherited `+3 a_nu I` truncation that matches V13), the
sweep/refresh operators agree to `2.6e-8`, and checkpoints are
bit-exact.  The convergence failures arise from an **inadequately
conditioned coupled system** with three compounding ingredients:

1. **Field side (suspected)**: `|Q| -> 1` strong-coupling compressing
   cells (94% of the grid) give near-unity gain per frequency hop -
   the leading candidate for a marginally contractive lagged map.
   This is the suspected mechanism, not a proven non-contractive
   iteration mode; P1-T1 (finite-difference `I_NEXT` sensitivity)
   remains outstanding.  (Pillar 1 - corrected mechanism.)
2. **Population side**: `C_MAT *= POPS` column scaling produces
   approximate null columns at floor populations, yielding
   `1E21-1E53` correction demands and, on `VAL_knop_formal`, a
   luminosity-divergent runaway.  (Pillar 5.)
3. **Deliberate modifications**: `KN_NO_GAIN` floors real maser gain
   (modified equation) and the `+3 a_nu I` term is omitted by V13
   convention (physical approximation).  Neither is a bug, but both
   mean the converged object - if reached - is not the literal Knop
   equation.  (Pillars 2, 4.)

Two distinct kill mechanisms are confirmed on the production model:
field-maturity stall (`FRESH_knop`, `GS_REL_NOPROG_LIM`) and
population/luminosity divergence (`VAL_knop_formal`).  Ordered
marching is the leading untested field-side remedy but is
runtime-unverified (Pillar 3).  This is consistent with the standing
verdict: **KNOP is not production-ready**; `DO_NONMONO_SOB` +
`USE_LAM_ES=T` remains the production route.

## Evidence classification summary

| Item | Class |
|---|---|
| `\|Q\|->0` as `GAM->0`; near-zero-GAM cells contractive | MATH confirmed; prior phrase retracted |
| `\|Q\|->1` strong-compression marginal contraction | MATH confirmed (mechanism) |
| Field residual plateau = intrinsic non-contraction | HYP (needs P1-T1 FD measurement) |
| `KN_NO_GAIN` modifies physical extinction / suppresses real gain | MATH + APPROX confirmed |
| `+3 a_nu I` omitted by V13 convention | APPROX confirmed |
| Header "cancels identically" comment | MATH defect (doc, open) |
| Population floor-pinned approximate null-space | MATH confirmed (fast model), corroborated |
| `FRESH_knop` terminal stall = field-maturity cycle | MATH confirmed (distinct from runaway) |
| Ordered march actually executed in dead runs | MISSING (no branch counter) |
| OBS D1-D4 / cold `I_NEXT` sink | DEFECT confirmed, fixes unapplied |
| `+3 a_nu I` impact on non-monotonic solution | HYP unquantified |

## Source references

- `plane/fg_j_cmf_nonmono_v1.f`: GAM/Q/chi_hat ~2815-2822, 3738-3857;
  S' ~4065-4152; store write ~4814-4819; KN_COUPLING_OK ~5509-5527;
  SWEEP_FLD/ORD dispatch ~7591-7669, 7913+, 8480+; KN_RF_COUPLING
  ~9604-9678, 9824-9829; WRITE_STORE ~10202-10216.
- `new_main/cmfgen_sub.f`: KN_NO_GAIN floor ~2555-2643, 6314-6362;
  census ~3599-3611.
- `new_main/mod_subs/solveba_v13.f`: STEQ_RAW ~1178-1187; demand scans
  ~1533-1608; RUNAWAY census ~2152-2237.
- `new_main/mod_subs/generate_full_matrix_v3.f`: `C_MAT *= POPS`
  ~385-390; conservation rows ~254-277.
- `new_main/mod_subs/cmf_blkband_v3.f`: PIN_NULL ~904-1056;
  POP_AUDIT census ~936-949, 1018-1033.
- Runtime: `m15_2p82d_FRESH_knop`, `FRESH_knop_d1a`,
  `VAL_knop_formal` OUTGEN / ITER_LEDGER / GS_STATE / NEW_SN_R_GRID.
