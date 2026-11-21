# Census failure/heal instrumentation design

Status: DESIGN + PARTIAL INSTRUMENTATION. The 18-column KNRAY_CENSUS row
(writer KN_CENSUS_ROW, called once per ray per armed refresh under the
KN_RDUMP_LS* keys) is source-reviewed and compiled into phase227_io.
Runtime tests are prepared, not run. No production classification logic
was duplicated: the writer exports KN_CS_LS/KN_CS_REA as computed by
KN_CERTIFY itself; validation only checks the table's own consistency.

## Captured per refresh x ray

serial (KN_DIAG_REFRESH), ray, KN_CS class, KN_CS_REA reason, exit
class NXT, bad-cell count NBD, max residual, ray geometry KN_P/NI,
per-channel RES/IXIT/NPAS/OMEGA (KN_PA_*), production-mask flag.

## Derived offline (no solver change needed)

per-ray classification timeline, first failure refresh, first recovery
refresh, consecutive-FULL and consecutive-non-FULL runs, transition
matrix, repeated-switching vs progressive convergence. A final FULL
never erases earlier rows: history is the join over the serial key.

## Missing for the required analysis

1. Max-defect location and frequency: KN_PA_IX is the exit code, not
   a cell index. Per-cell defect detail exists only for armed rays in
   KNSTAGE/KNCELL dumps. Needs a census field or per-ray join to a
   defect-argmax export; NOT yet implemented.
2. First failure PASS within a refresh: row exports NPAS only; a ray
   that fails early, heals, fails again inside one solve is invisible.
   Per-pass residual history is captured for armed rays only.
3. Q residual / bulk residual / strict-cell count: KN_CERTIFY inputs
   are not exported per ray; only the verdict and NBD reach the row.
4. Probe handling: KN_PROD_MASK is exported per row; consumers must
   not assume ray 83 is always the probe.

## Refresh serial review

KN_DIAG_REFRESH increments once, after the eligibility guards and
before coefficient rebuild, inside FG_KNOP_REFRESH_V1 - only when an
armed KN_RDUMP_LS* key is set. Skipped refreshes (early returns before
the increment) cannot burn a serial, so no gap or duplicate arises
there. Unique per process; a restarted capture must use a fresh
directory. Runtime checks (multi-iteration monotonicity, early-return,
re-entry) are in the prepared phase227_io plan; NOT RUN.
