# KNOP Refresh FLD / Residual Audit (P1-T1f, PART 6)

Status: **telemetry largely already present; the independent
fixed-point residual `q` is specified but NOT yet implemented.**
This report maps every PART 6 requirement to the existing mechanism
in `KNOP_RAY_SOLVE_V1` (`plane/fg_j_cmf_nonmono_v1.f`, lines
7319-7920) and isolates the single remaining gap.

Scope: `KNOP_RAY_SOLVE_V1`, its internal procedures `SEED_FLD`,
`SWEEP_FLD`, `SWEEP_FLD_ORD`, `MARCH`, `BUILD_FREQ`, `MATCH_FREQ`,
`SRC_PRIME`, and the commit writers `WRITE_STORE`/`WRITE_STORE_CAP`.
No production numerical equation was modified; no refresh convergence
criterion was changed.

Evidence classes: **[MATH]** code structure, **[MANUF]** manufactured
`drv_knop_ray` model, **[PROD]** real non-monotonic model.

---

## 1. PART 6 requirement -> existing mechanism

`KNOP_RAY_SOLVE_V1` is called once per ray inside the OMP-parallel
refresh loop in `FG_KNOP_REFRESH_V1`. Its internal procedures
host-associate the routine's locals, so every per-ray variable below
is private to that ray's solve (thread-safe by construction).

| Required record | Existing mechanism | Status |
|---|---|---|
| initial `FLD` seed | `SEED_FLD(FLD, I_M_DEPTH_NEXT)` / `SEED_FLD(FLD_P, I_P_DEPTH_NEXT)` fill from the lagged store (lines 7647, 7702) | PRESENT (reconstructable) |
| per-pass `RES` | `KN_RES_PASS` module array accumulates each pass's max residual (allocated line 6907) | PRESENT |
| `NBAD`, `NTOT` | returned via `NBAD`, `NTOT`, `NB_M/NT_M`, `NB_P/NT_P`, plus best-pass `NBBLK/NTBLK` | PRESENT |
| iteration count | `NPASS` return; per-direction `NP_M`, `NP_P` | PRESENT |
| `IXIT` exit code | `IXIT_M`, `IXIT_P` -> `NXIT` return; codes 0-6 documented at 7437-7444 | PRESENT |
| iteration-cap status | `IXIT=4` (pass cap), `IXIT=3` (stall), `IXIT=6` (plateau) distinguished | PRESENT |
| final internal `FLD` | `KR_DBG_FLDM`, `KR_DBG_FLDP` copied at 7830-7835 for ray `LS == KN_DBG_LS`, before `DEALLOCATE` - i.e. the converged internal field regardless of commit path | PRESENT (per-ray dump) |
| committed store | `WRITE_STORE(FLD/FLD_P, I_*_DEPTH_NEXT)` at 7811-7812 | PRESENT |
| `WRITE_STORE_CAP` activity | capped provisional commit at 7792-7793; verdict `FG_RAY_VC=3` | PRESENT |
| rejected commit | jump-veto `FG_RAY_VC=1` (7748-7753), nc-reject `FG_RAY_VC=2` (7807) | PRESENT |
| `HEAL_STORE` activity | diverged-ray heal at 7851-7857; verdict `FG_RAY_VC=4` | PRESENT |
| floor and cap activity | `NNEG`/`NNG_P` (negative floor), `NCAP`/`NCP_P` (`NONMONO_INEXT_CAP`), returned as `NCAP`, `NNEG`, `NCAP_LST`, `NNG_LST` | PRESENT |
| per-ray commit verdict | `FG_RAY_VC(LS)` in {0 committed,1 jump-rej,2 nc-rej,3 capped,4 div}; printed as the "commit verdicts" census at 7224-7228 | PRESENT |

The refresh already distinguishes commit from convergence: a ray can
report a small `RESID` yet be classed unconverged by `IXIT`, and the
commit verdict records whether the field was actually written. The
cap-vs-raw and floor-vs-cap channels are separate (`NCAP` vs `NNEG`).

## 2. The missing piece: independent residual `q`

PART 6 asks for

```
q = norm(FLD - REFRESH_OPERATOR(FLD)) / max(norm(FLD), scale_floor)
```

evaluated **without committing** to the production store.

What exists today is `RES_M`/`RES_P`: the **max-overwrite** residual
of the last completed sweep pass, `max_cells |VNEW-VOLD|/scale`, where
`DCELL` is computed per cell in `MARCH` (line 10050) and the ordered
path's inline march. That is an `L-inf`-style measure of
`||T(FLD)-FLD||`, **not** the requested `L2` operator residual, and it
is measured *during* the pass rather than as an independent
re-evaluation.

### Why it matters

`RES` and `IXIT` are self-reported by the same loop that produced
`FLD`. An independent `q` re-applies the operator to the converged
field and measures the fixed-point defect directly. It separates
"the loop thinks it converged" from "the field is actually a fixed
point of `T`" - the distinction PART 6 wants for telling real
convergence apart from cap/stall/tolerance artefacts.

### Implementation plan (NOT yet done)

The per-cell delta is already formed inside `MARCH` as
`VNEW-VOLD` (and `FLD(I+1,L)=VOLD+OMEGA_RAY*(VNEW-VOLD)` applies the
relaxed update). The ordered path retains `FLD_PRV`/`DLT_PRV`
(previous iterate and pass delta) for Anderson acceleration. A
surgical `q` needs, per direction field:

1. accumulate `sum((VNEW-VOLD)/scale)^2` and `sum(VOLD^2)` (or the
   applied `OMEGA*(VNEW-VOLD)` delta, matching whichever residual
   convention is chosen) over all marched cells during the **last**
   pass only;
2. expose `q = sqrt(sum_delta2)/max(sqrt(sum_fld2), floor)` through a
   new return or module diagnostic, gated on a flag default-off;
3. cover all three march paths: `MARCH` (unordered), the `SWEEP_FLD`
   inline residual block, and `SWEEP_FLD_ORD`'s ordered march.

Because `MARCH`/`SWEEP_FLD`/`SWEEP_FLD_ORD` are internal procedures of
`KNOP_RAY_SOLVE_V1`, a per-ray accumulator local to
`KNOP_RAY_SOLVE_V1` is thread-safe (each OMP ray gets fresh locals).
No OMP reduction is needed inside `MARCH` (it is sequential within a
ray); the accumulator is per-ray, summed over the ray's cells.

This is a hot-kernel change: it must not alter `RES`, `NBAD`, the
update `FLD=VOLD+OMEGA*(VNEW-VOLD)`, or any convergence gate. It adds
one multiply-add per cell on the diagnostic path only.

### Manufactured-test plan (PART 6, once `q` lands)

In `drv_knop_ray`, run the refresh cases (`rfsh`, `rf3`, `rfail`,
`scold`) across several seeds and progressively tighter
`KN_RAY_TOL`/pass caps; record per-ray `q` beside `RESID`, `IXIT`,
commit verdict. Expect `q ~ RES`-consistent at converged exits and a
divergence between `q` and the reported `RESID` exactly on the
cap/stall/plateau exits - that divergence is the diagnostic signal.

For a production checkpoint, prioritise rays with alternating
expanding/compressing bands (`HAS_COMPRESSING` with sign changes in
`GAM`) and the largest `KN_RF_RESID` / `FG_GS_NBAD` residuals.

## 3. What is and is not established

- **Established [MATH]:** every PART 6 telemetry item except the
  independent `L2` residual `q` already exists and is wired to the
  commit/convergence distinction.
- **NOT established:** `q` itself; refresh seed-independence on a
  production model; whether cap/stall exits mask unconverged fixed
  points. None is claimed.

## 4. Remaining work

- Implement `q` (section 2) behind a default-off diagnostic flag.
- Manufactured seed/tolerance sweep to validate `q` against `RES`/`IXIT`.
- Production run once a PART 1 checkpoint exists.
