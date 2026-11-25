# kv-scan Phase 0 inventory - WRITE/PRINT multi-field output sites

2026-11-25. Scanner: `tools/scan_kv_writes.py` (stdlib-only).
Tokens: `tools/parser_tokens.txt` (Phase 1 parser strings).
CSV: `docs/audits/kv_scan/inventory.csv`. Re-run:

```bash
python3 tools/scan_kv_writes.py . \
    --tokens tools/parser_tokens.txt \
    --csv docs/audits/kv_scan/inventory.csv
```

Scans every `.f/.f90/.for/.inc` in the tree (1618 files), assembles
logical statements (fixed-form col-6 + tab rules, free-form `&`,
form-feed-safe line numbers), resolves labeled/inline/list-directed
formats, traces units through `OPEN`/`GEN_ASCI_OPEN`, tracks DO-loop
and single-line-IF context, and reports output statements emitting
>=2 label=value / label:value / 'label N' fields on one output line.

## Totals: 2426 sites

| kind | count | meaning |
|---|---|---|
| kv-line | 958 | several labelled scalars on one line - convert |
| table-row | 481 | inside a DO loop - table or leave, never vertical |
| data-file | 391 | data/restart/model file or unformatted write - do not touch |
| parsed | 97 | emits a string external consumers grep - protect |
| unclear | 499 | unresolved format or unit - needs human decision |
| namelist | 0 | none in tree |

### By directory x kind (all dirs)

| dir | kv-line | unclear | table-row | data-file | parsed |
|---|---|---|---|---|---|
| blas | 1 | 0 | 0 | 0 | 0 |
| disp | 43 | 70 | 23 | 12 | 0 |
| dust | 3 | 0 | 0 | 0 | 0 |
| lte_hydro | 8 | 6 | 4 | 4 | 0 |
| misc | 40 | 57 | 40 | 14 | 0 |
| new_main | 306 | 105 | 139 | 82 | 71 |
| newsubs | 6 | 6 | 17 | 2 | 0 |
| obs | 50 | 6 | 43 | 65 | 0 |
| opac | 12 | 0 | 3 | 8 | 0 |
| pgplt | 41 | 8 | 24 | 9 | 0 |
| plane | 148 | 69 | 104 | 94 | 19 |
| spec_plt | 49 | 72 | 23 | 17 | 0 |
| stark | 5 | 1 | 1 | 0 | 0 |
| subs | 83 | 36 | 26 | 63 | 6 |
| tests | 63 | 11 | 2 | 0 | 0 |
| tools | 11 | 7 | 1 | 0 | 0 |
| unit_tests | 90 | 45 | 31 | 21 | 0 |

### By program area

| area | kv-line | table-row | parsed | data-file | unclear |
|---|---|---|---|---|---|
| core (cmfgen: new_main, plane, subs, newsubs, lte_hydro, misc, tools, stark, opac, dust) | 621 | 335 | 97 | 267 | 287 |
| support programs (disp, obs, spec_plt, pgplt) | 183 | 113 | 0 | 103 | 156 |
| test drivers (unit_tests, tests) | 153 | 33 | 0 | 21 | 56 |

## kv-line field-count distribution (core area)

| fields | sites |
|---|---|
| 2 | 364 |
| 3 | 91 |
| 4 | 60 |
| 5 | 39 |
| 6 | 15 |
| 7 | 8 |
| 8 | 13 |
| 9 | 14 |
| 10 | 9 |
| 11 | 3 |
| 14 | 2 |
| 15 | 2 |
| 17 | 1 |

Per the prompt's edge-case rule, most of the 364 two-field lines stay
single-line (if under 100 columns). The definite vertical-block pool
is ~257 core sites with >=3 fields plus wide two-field stragglers.

## Top sites by field count

| fields | site | kind |
|---|---|---|
| 20 | new_main/subs/cmf_sub_sonic_hydro.f:930 DO_CMF_HYDRO_V3 | unclear |
| 18 | new_main/mod_subs/solve_for_pops.f:1857 SOLVE_FOR_POPS | parsed |
| 17 | plane/fg_j_cmf_nonmono_v1.f:11856 KNOP_RAY_SOLVE_V1 | kv-line |
| 17 | plane/jgrey_hub_ddt_v{3,4,5}.f (5 sites) | unclear |
| 15 | new_main/subs/ddt_work.f:328,342 DDT_WORK | kv-line |
| 14 | new_main/subs/cmf_sub_sonic_hydro.f:1544 CHECK_HYD | unclear |
| 14 | new_main/subs/smooth_pops_as_we_iterate.f:32 | kv-line |
| 14 | plane/fg_j_cmf_v{12,13}.f FG_J_CMF | table-row |
| 14 | plane/mom_j_cmf_nonmono_v1.f:5579 MOM_KNOP_REFRESH_V1 | kv-line |
| 14 | subs/hydro_terms_v5.f:185 HYDRO_TERMS_V5 | data-file |

(Full list in inventory.csv; repeated v2/v3/v4/v5 file families count
separately.)

## Phase 1: parsers found and how each is handled

Consumers of OUTGEN/stdout (workspace `scripts/` + run scripts):

| parser | strings depended on | handling |
|---|---|---|
| `check_cmfgen.py`, `cmfgen_convergence.py`, `plot_cmfgen_diagnostics.py`, `plot_maxch_comparison.py` | `Luminosity of star \(d=1,ND\)\(iteration N\) is:`, `Maximum % increase/decrease at depth D is V (...) --- iteration N` | `parsed`: keep wording verbatim |
| `iter_phase0_metrics.py` | `ITER_SUM: it= mode= -> ba= frozen= held=`, `ITER_TYPE it=`, `MODE-OUT it`, `Time taken for iteration N:`, `pops/se residual/ba/lum diag` metric rows, `next lam / comp_ba / fix_t`, `masked depth columns`, `soft rejects`, `BA_RES_ACC (PARTIAL) REJ it`, `Error performing NG acceleration`, `Averaging of oscilating` | `parsed` |
| `check_convergence.py` | `ITER_MODE ... reason=`, `Maximum %`, `--- iteration N`, `Max(imm|imum) changes as returned by SOLVEBA_V13`, `EPS_TERM` | `parsed` |
| `parse_knop_iso.py`, `parse_knop_seed.py` | `ITER_SUM: it=`, `Model Convergence: N %`, `Luminosity of star`, `Time taken for iteration`, `KNOP I_NEXT refresh` blocks, `matched/compress=`, `KN-JT/KN-NOP/KN-SEED` lines, `ray exits:`, `commit verdicts:`, `release_no_progress`, `adaptive INEXT_RELAX`, `q argmax/q ray chmax`, `FG-KNOP worst ray`, `unconverged cells`, `retained-field delta`, `residual max=`, `eod even/odd`, `solves M ord`, `fg_rf`, `POPULATION SOLVE` boundary | `parsed` |
| `run_model_grid.py` | `Number of frequencies is:`, `[KEY]` VADAT echo lines, `TcurA`/`TcurB`, `EKO`/`kin.energy` | `parsed` (echo lines are already 1-field) |
| `submit_cmfgen.sh` (generated), `jfnk_watch.sh` | `Maximum % increase at depth`, `Current great iteration count is`, `MAXCH demand / applied`, `NG Acceleration performed`, `depth-block columns`, `post-accel demand`, `JFNK:` count | `parsed` |
| `unit_tests/python/*.py` | static source checks pin literals: `' Current iteration count is'`, `' Time taken for iteration'`, `'Luminosity of star'`, `' POPULATION SOLVE'`, `'Model Convergence:'`, `'ITER_TYPE'`, `' DEC_VEC: '/' INC_VEC: '`, `' Maximum % increase/decrease at depth '`, `' RUNAWAY it='`, `'Maximum changes as returned by SOLVEBA_V13 is'`, `'Average number of calls to FG_J_CMF is'` | `parsed`; tests updated in same commit when a site converts |

97 sites are marked `parsed`. Where a parsed site is already
sentence-style (`Maximum % increase at depth 6 is 1.6E+10 (LAMBDA)`)
it stays verbatim; the `Maximum %` lines additionally appear inside a
sed `s/.*is *\([0-9.E+]*\).*/\1/` extraction in generated
submit_cmfgen.sh requeue logic - value must remain the last `is NUM`
on the line.

Single-field parsed lines (`Model Convergence:`, `Luminosity of
star`, `DEC_VEC:`/`INC_VEC:` vectors, `Time taken for iteration`)
are below the >=2-field threshold and are inherently out of scope;
they are listed here because the wording must survive any nearby
edit.

## unclear breakdown (499)

| reason | count | resolution path |
|---|---|---|
| multi-open | 197 | literal unit opened on different files in different programs; decide per program |
| unresolved-var | 156 | unit var never assigned/opened in file (dummy arg, GET_LU()) |
| unopened-literal | 53 | fort.N or unit never opened anywhere scanned |
| fmt-var | 50 | format held in a variable/expression; read the assignment |
| multi-file | 42 | unit var opens several files |

## Scanner limits (known approximations)

- `header_text` is a heuristic (first literal up to first `:`/`=`).
- `est_width` approximates field widths (A->8 when standalone, 13
  per list-directed item); ordering only.
- `n_fields` counts label->value positions, max over `/`-separated
  lines; `[x,y]` pairs under one label count once.
- String variables are probed for parser tokens via per-routine
  fragment collection (400-char cap); a var reused for unrelated
  text can over-flag `parsed`.
- Internal writes to CHARACTER variables are excluded (string
  builds, e.g. `WRITE(ITER_MSG,...)`, `WRITE(FMSG,...)`).
- `data-file` errs to the safe side: unformatted opens, `REC=`,
  no-FMT writes to file units, all-data multi-open names.
- `unclear` errs to the cautious side: runtime-computed file names,
  ambiguous units, variable formats.

## Notes for Phase 2 planning (not yet started)

- `new_main/mod_subs/iter_log_fmt_mod.f90` already exists and covers
  much of the requested helper module: `ilf_section`, `ilf_kv*` (30-
  wide `key : value` style), `ilf_table_*`, `ilf_metric_*`,
  `ilf_es` (ES11.4E2 with ES11.3E3 fallback, NaN handled), tagged
  diagnostics, issue census, `OUT_LEVEL`/`CONV_DIAG_LVL` verbosity
  0..3. Round-1/round-2 refactor work (ITER_SUM table, moment-
  refresh table, RUNAWAY census) is already landed.
- Gaps vs the prompt's spec: block-scope label-width padding,
  `label = right-aligned value` style (existing `ilf_kv*` uses
  `key : value`), logical typed variant, `LOG_COMPACT_KV` switch,
  dedicated edge-case unit tests (1E+99/1E-99/0/10-digit int/long
  labels/100-col cap).
- Decision needed: adopt the prompt's `=` + right-aligned numeric
  style as a second kv-block helper, or restyle `ilf_kv*` - the
  prompt's target format differs from the style rounds 1-2 shipped.

## Verdict

25-Nov-2026. PARTIAL - the inventory itself is VERIFIED for the sites
the scanner can resolve: line numbers confirmed against grep on the
ES_ALI site (`new_main/mod_subs/comp_j_blank.f:1401-1411`), the parsed
set covers every Phase-1 token including variable-built strings
(`Luminosity of star` via `STRING` at `cmfgen_sub.f:3623`), and direct
evidence was spot-checked across all five kinds. REJECTED on
completeness by automation alone: 499 `unclear` sites and the
table-vs-leave decision on 481 `table-row` sites need the human scope
review this document exists for. Evidence: `inventory.csv` (2426
rows), spot-check output, `tools/scan_kv_writes.py` header notes.

## Verdict (update)

25-Nov-2026 (phase 3 close). VERIFIED - core conversion complete.
Of the 122 real `nf >= 3` kv-line sites in core dirs (new_main 84,
plane 30, subs 4, stark 4), all are now either `ILF_KVB` blocks or
deliberately left: per-ray/per-Krylov-direction loop telemetry,
embedded-sentence messages (alloc-fail, DEFER-PROBE prose), parsed
tokens, and `rd_bs_strk_tab` expected-vs-read pairs (unclear default).
`rd_control_variables.f` (7 sites) and `mom_j_cmf_nonmono_v1.f` (3)
were skipped while carrying uncommitted edits. Helper:
`ilf_kvb`/`ilf_v{i,r,l,s}` in `iter_log_fmt_mod` with the temporary
`CMFGEN_COMPACT_KV` env hook (removal owed at refactor close).
Tests: `unit_tests/python/test_iter_log_kvb.py` (26 checks).
