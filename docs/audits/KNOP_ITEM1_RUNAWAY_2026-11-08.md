# KNOP Item 1 - runaway-level diagnostic (2026-11-08)

Status: diagnostics only. No solver behaviour changed; default path
bit-identical. Hypothesis under test: **H1** - floor-pinned trace
populations give the per-depth SE block an approximate null space, so
`SOLVEBA_V13` demands order-1E7 corrections every iteration and the
exit gate can never pass. This note adds the instrumentation that
names the exact channels; it does not yet fix anything.

## What was added

Two opt-in VADAT keys (defaults preserve current behaviour):

- `POP_RUNAWAY_DEM` (REAL, default `1.0E3` = 1E3 percent). When the
  largest demanded fractional correction on any level/conservation
  row reaches this, the census prints. `<=0` disables.
- `POP_RUNAWAY_TOPN` (INTEGER, default `10`, clamped to `[0,200]`).
  Number of top-demand channels printed per call.

### Channel census - `RUNAWAY` lines in `SOLVEBA_V13`

`new_main/subs/solveba_v13.f`, inside the `POPS_OLD` block just before
it is deallocated (after the committed `MAXCH_APP` scan and any
`BA_MAXCH_FAIL` restore). For each call whose max demanded correction
on rows `J=1..NT-2` reaches `POP_RUNAWAY_DEM` percent:

- a header `RUNAWAY it=<counter> maxdem%=<> thresh%=<> topn=<>`;
- one line per top channel:
  `RUNAWAY d=<depth> eq=<J> <species>:<level> n=<n_old> ionf=<frac>
   dem%=<> app%=<> fl=<> cons=<>`

Field meanings:

- `d`, `eq`: depth index and global equation/variable index.
- `<species>:<level>`: `EQ_DESC(J)` label (`ION:levelname`,
  `ION:ion-eq`, `T`, `ED`).
- `n`: pre-update population `POPS_OLD(J,I)` (the state the demand
  was computed against).
- `ionf`: `n` as a fraction of the owning ion's level-row total at
  that depth (`-1` for non-level rows).
- `dem%`: demanded fractional correction `100*STEQ(J,I)` (signed;
  `+` = decrease). This is the *raw* demand, including rows the
  `BA_DEM_PFLOOR` metric exempts - that is the point of the census.
- `app%`: applied fractional change actually committed,
  `(POPS-POPS_OLD)/POPS_OLD`, so scale limits, the Stage-B trace
  hold, and any `BA_MAXCH_FAIL` restore are all reflected.
- `fl`: `1` when `n <= 1E-30` or `n < BA_DEM_PFLOOR` (the floor the
  JFNK audit found hosting the near-null SE columns).
- `cons`: `1` for non-level conservation/ion-eq rows.

A new contained helper `EQ_LOC(J,ID,LEV)` decodes the equation index
to (ion ID, 1-based level) using the same `ATM(ID)%EQXzV`/`NXzV`
ownership scan as `EQ_DESC`, so the ion-stage fraction and the
conservation flag do not re-parse the label string.

### Per-depth floor/column census - `RUNAWAY` lines in `CMF_DIAG_SOLVE`

`new_main/subs/cmf_blkband_v3.f`, gated by the existing `POP_AUDIT`
logical (now imported from `CONTROL_VARIABLE_MOD`). For every depth
`K` of the DIAG band solve, after `GENERATE_FULL_MATRIX_V3` +
`CMF_BLKBAND_PIN_NULL_V1` and again after `DGEEQU` equilibration:

- `nflr`: count of unknowns `J=1..N-2` with
  `POPS(J,K) <= max(1E-30, BA_DEM_PFLOOR)`;
- `colmin_u/colmax_u`, `colmin_e/colmax_e`: min/max over columns of
  `max_i |C_MAT(i,j)|` on the unequilibrated and equilibrated block.

One line per depth:
`RUNAWAY d=<K> nflr=<>  colmin_u=<> colmax_u=<> colmin_e=<> colmax_e=<>`.

This is the cheap column-norm proxy the JFNK audit (section F)
called for: a column whose norm collapses toward `1E-16..1E-22` of
the block scale while `nflr` is large is the floor-pinned null column
behind the demand cap. No SVD, no factorisation change; the solve is
untouched. TRI path is not instrumented (DIAG is the production
route); noted as a limitation.

## Files changed

- `new_main/mod_subs/control_variable_mod.f` - declare the two keys.
- `new_main/mod_subs/rd_control_variables.f` - `RD_STORE_DBLE` /
  `RD_STORE_INT` parsing + clamps.
- `new_main/subs/solveba_v13.f` - USE import, locals, `EQ_LOC`,
  channel census.
- `new_main/subs/cmf_blkband_v3.f` - USE import, locals, per-depth
  floor/column census.
- `unit_tests/f90/drv_vadat.f90` - echo the two keys.
- `unit_tests/python/test_vadat_keys.py` - default + explicit parse
  cases.

## Verification

- `ifort -extend-source` syntax compile of all four edited `.f`
  files: clean.
- `python3 unit_tests/python/test_vadat_keys.py`: **85 passed,
  0 failed** (includes the two new POP_RUNAWAY cases).
- Full `./build.sh` + `unit_tests/run_all_tests.py`: not yet run end
  to end in this session; targeted compiles above cover the edited
  translation units.

## What to collect from the next run

- `RUNAWAY it=` header + channel lines: which (species,ion,level,
  depth) channels carry the 1E7 demand, their `n`, `ionf`, `fl`,
  `cons`.
- `RUNAWAY d=` per-depth lines (need `POP_AUDIT=T`): `nflr` and the
  `colmin/colmax` column-norm extrema - confirms whether the weakest
  columns sit on floor populations (H1) before Item 2 masks them.
- Cross-check against the existing `pops demand/applied`,
  `worst margin (exm)`, and `BA_DEM_EXM/EXR` lines.

## Not done / next

- Item 2 (`POP_MASK_FLOOR`, `POP_TRACE_MAXFAC`) waits on this output
  to confirm the masked set matches the near-null columns.
- The census runs on every `SOLVEBA_V13` call whose max demand
  reaches the threshold; on a stalled run that is every iteration.
  Cost is two O(NT*ND) scans plus TOPN argmax passes - negligible
  next to the band solve. `POP_AUDIT` gates only the heavier
  per-depth column census.
