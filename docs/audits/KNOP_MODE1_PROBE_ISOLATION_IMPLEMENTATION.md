# KNOP MODE 1 diagnostic-probe isolation - implementation record

Date: 20-Nov-2026. Base revision: `90602193`.
Companion audits: `KNOP_PROBE_CONTROLLER_AUDIT.md`,
`KNOP_ACTIVE_FIELD_VALIDITY_AUDIT.md`, `KNOP_THREE_SEED_EXPERIMENT.md`.

## 1. Scope

Implements MODE 1 isolation for `KN_NOCP_PROBE` (the non-committing
diagnostic probe ray, `KN_NOCP_PROBE_LS`, ray 83 in the running
experiment). The probe keeps its seed injection, ray solve, per-ray
residuals and all KN-NOP/KN-SEED telemetry. It no longer feeds any
quantity the global controller reads:

- the refresh residual maximum (`KN_RF_RESID`),
- the adaptive-omega sample pool (`KN_OMEGA_POOL_V1` / `RHO_P90`),
- the share censuses (`FG_GS_NBAD/NTOT/NBAD_BLK/NTOT_BLK/NDIV`,
  `FG_RF_CAPPED`, `FG_INEXT_CAPPED`, the exit-class bad-cell sums
  behind `FG_GS_NSTAL/NPLAT`),
- the pending-delta maturity input (`FG_RF_PEND`),
- the commit/rejection counters (`FG_COMMIT_REJ`, `FG_NC_REJ`,
  `FG_NC_CAP`, `FG_DIV_HEAL`).

Not changed: any solver equation, tolerance, threshold, omega
formula, `KN_RHO_MIN_SAM`, `INEXT_RELAX` bounds, the MARCH
recurrence, or the probe's own solve path. MODE 2 (committing
probe) is untouched; `KN_NOCP_PROBE=F` reproduces the old path
exactly.

## 2. Source changes (all in `plane/fg_j_cmf_nonmono_v1.f`)

### 2.1 `KN_PROD_RESID_V1(RES,NDV,XRAY)` - new pure helper (:2037)

Controller residual pool. `XRAY` is the armed probe index
(`0` = no probe; `KN_NOCP_PROBE_RAY` is the single source of the
verdict so the mask can never disagree with the commit guards).

- `XRAY=0` reproduces the old `MAXVAL(RES,MASK=NDV==0)` exactly.
- All production rays diverged -> maximum over production
  residuals (diverged values), not a silent zero.
- Degenerate no-production-ray run -> unmasked `MAXVAL(RES)`,
  preserving the pre-change value for that pathological case.

`KN_RF_RESID` now calls this helper (:8904) after
`KN_PROD_MASK`/`KN_PROD_XRAY` are built (:8889-8895).

### 2.2 `KN_OMEGA_POOL_V1` - optional `EXCL` arg (:2078)

Last positional arg, `INTEGER, OPTIONAL`. When present the ray at
`EXCL` is skipped at the sample-pool fill (:2093); `RHOM/RHOP`,
`NDV` filtering, positive-rho filtering, sort, percentile
convention, `KN_RHO_MIN_SAM` retain path and the omega formula are
all unchanged. Call site passes `KN_PROD_XRAY` (:8912). An absent
or `0` EXCL samples every ray - the probe-off path is the old code
byte for byte.

### 2.3 `KN_PROD_MASK` / `KN_PROD_XRAY` - refresh-locals

`LOGICAL KN_PROD_MASK(KN_NP)` allocated with the other per-ray
work arrays (:8249), filled at :8889-8895 from
`.NOT. KN_NOCP_PROBE_RAY(LS)`, deallocated at :9611.
`KN_PROD_XRAY` records the single excluded index for the
index-based reductions.

### 2.4 Census reductions masked (:8955-8970)

`FG_RF_CAPPED`, `FG_INEXT_CAPPED`, `FG_GS_NDIV`, `FG_GS_NBAD`,
`FG_GS_NTOT`, `FG_GS_NBAD_BLK`, `FG_GS_NTOT_BLK` are `SUM(...,
MASK=KN_PROD_MASK)`. The printed `diverged rays=`/`clamp events=`
lines use the same mask; a new `probe ray=` line prints the
probe's own resid/div/nbad next to them so the excluded
contribution is still visible in OUTGEN.

### 2.5 Exit-class bad-cell sums (:9008-9020)

`NX_BAD(NXT_LS(LS))` accumulates `NBD_LS` only for production
rays: those sums feed `FG_GS_NSTAL/NPLAT`, which are gate inputs.
`NX_CNT` (exit census) and `NX_RMX` (per-class max residual)
remain probe-inclusive - they are printed diagnostics, not gate
inputs.

### 2.6 `FG_RF_PEND` masked (:9518-9532)

`MAXVAL(FG_RAY_DCH, MASK=(VC==1 .OR. VC==2))` gains
`.AND. KN_PROD_MASK`. The probe's veto/reject verdicts carry real
deltas in `FG_RAY_DCH` (9.8E4 in the prbv test) and would
otherwise hold a phantom production debt on `RF_DCHG` ->
`GS_MAX_ALL`. Belt-and-braces with the later verdict-7 rewrite:
correct under either evaluation order.

### 2.7 Commit/rejection counters (:10220, :10303, :10315)

`FG_COMMIT_REJ`, `FG_NC_CAP`, `FG_NC_REJ` increments are wrapped
in `IF(.NOT. L_NOPROBE)` inside their OMP CRITICALs. The probe's
verdicts still record on `FG_RAY_VC`/`FG_RAY_DCH` before the
code-7 rewrite, so the `would-be verdict=` KN-NOP line reports
what the controller would have seen. `FG_DIV_HEAL` needed no
gate - it already lives inside the `.NOT. L_NOPROBE` heal block.

### 2.8 Ordinary sweep writes deliberately retained (:5750)

The in-sweep `I_P_DEPTH_NEXT`/`I_M_DEPTH_NEXT` update
(WRITE ... T2+OMEGA*(T1-T2)) is NOT gated on `L_NOPROBE`, and a
comment block states why: it reads the production `I_P` field and
the retained slot snapshot, so it is arm-independent production
state; gating it would change the model, not isolate the
diagnostic. `GS_SW_MAXCHG` keeps the probe's sweep drift for the
same reason.

## 3. What stays probe-inclusive (telemetry only)

`RES_LS`, `RHO_M_LS`, `RHO_P_LS`, all per-ray census arrays
(`NBD/NTT/NBB/NTB/NCL/NNL/NDV/NCP/NNG_LS`), `NX_CNT` exit counts,
`NX_RMX` per-class maxima, `FG_RAY_VC`, `FG_RAY_DCH`,
`FG-KNOP worst ray=`, all KN-NOP/KN-SEED records and the probe's
seed files.

## 4. Data flow after the change

    per-ray solves -> RES_LS/RHO/NBD/... (all rays incl. probe)
                           |
             KN_PROD_MASK = .NOT. KN_NOCP_PROBE_RAY
                           |
      +--------------------+--------------------+
      | production-only (mask applied)           | probe telemetry (unmasked)
      | KN_RF_RESID, omega pool, FG_GS_* census, | 'probe ray=' line, VC/DCH,
      | NX_BAD, FG_RF_PEND, REJ/CAP/HEAL counters| KN-NOP, KN-SEED, worst ray
      +------------------------------------------+

## 5. Tests

New driver cases in `unit_tests/f90/drv_knop_ray.f90` (built by
`unit_tests/python/test_knop_ray.py`'s recipe) and assertions in
`unit_tests/python/test_knop_mode1.py`. Baseline is `rf3` (three
identical rays); cases `prb`/`prbv`/`prbs` arm ray 3 as probe with
its store scaled x2.5 / x1E-4 / x2.5+KN_SEED_ARM=1 - the sibling
medians are unaffected (median of {s,s,k*s} = s), so rays 1-2 see
identical inputs and every production export must equal `rf3`.

`punit` exercises `KN_PROD_RESID_V1` (excluded argmax, XRAY=0
parity, below-max probe, all-diverged fallback, diverged-probe
fallback, degenerate all-probe fallback) and `KN_OMEGA_POOL_V1`
EXCL (probe above/below P90, pool size 40 vs 42, min-sample
retain, diagnostic rho preserved).

Results (this tree, 20-Nov-2026):

- `test_knop_mode1.py`: 6/6 PASS. `prb`: RF_RESID=9.1124E-15 =
  rf3 bitwise, NTOT/NTB 180 vs 270 (probe's 90 cells out), probe
  resid 5.68E-15 still reported, VC(3)=7, KN-NOP delta=0,
  production RSTORE bitwise identical to rf3. `prbv`: would-be
  verdict=1 (real veto) with JREJ/NCREJ/NCCAP/HEAL=0 and
  PEND=0 despite FG_RAY_DCH=9.84E4. `prbs`: arm=1 seed ran
  (KN-SEED record present), production exports = rf3.
- `test_knop_ray.py`: 43/43 PASS - the probe-off production path
  is unchanged (bitwise store comparisons included).

## 6. Production compatibility evidence

- `XRAY=0`/`EXCL=0` paths are the old expressions verbatim
  (verified in `punit`: `M1U_RES_P0=5.0`, `M1U_OM_INCL` pool of
  42 samples with omega identical because the old pool is used).
- `rf3` (no probe) output is unchanged and still serves as the
  exact baseline for the armed cases.
- No new branches exist outside `L_NOPROBE`/`KN_NOCP_PROBE_RAY`/
  `KN_PROD_*` scope; the probe-off code path has no probe
  branches at all (the mask is all-true and reduces identically).

## 7. Remaining indirect coupling / limitations

- The probe's ordinary sweep writes are real production writes
  (Section 2.8) - by design, not residual contamination.
- `KN-NOP worst ray=` and exit-census prints are probe-inclusive
  (diagnostics only); `RF_STAL/RF_PLAT` production sums exclude
  the probe through `NX_BAD`.
- `FG_INEXT_CAPPED` now excludes the probe's same-sweep cap/floor
  events; its sweep writes still land in the store and in
  `GS_SW_MAXCHG` as intended.
- A probe whose solve diverges leaves `NDV_LS(probe)=1` in
  telemetry but out of `FG_GS_NDIV`; the freeze/release gate sees
  production divergence only.
- The running A/B/A-A experiment still uses the pre-change
  executable (`exe_frozen/`); nothing here modifies those jobs.

## 8. Verification artifacts

- Build: `make DEBUG=nodebug all` clean (only pre-existing #6075
  warnings), `exe/cmfgen_dev.exe` sha256
  `03ccd243848d1f3fac1ec0121f8d69e70c13135876e1da204b12e9a8355e5a37`.
- Test driver: `unit_tests/bin/knray` (production objects + driver).
