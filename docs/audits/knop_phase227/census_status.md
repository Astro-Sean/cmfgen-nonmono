# Per-ray census - BLOCKED

Inspection and provenance: parent Phase 2.27 report, 2026-11-21T05:55:59Z.
Phase 2.27B scheduler recheck: both original jobs remain PENDING
(JobHeldUser); neither was released or modified. No actual KNOMG or
KNRAY_CENSUS records were found. Existing OUTGEN aggregates do not identify
rays or transitions. No complete classification matrix can be produced.

The isolated candidate writer now calls KN_CENSUS_ROW, the same row
routine compiled into the synthetic Fortran format harness. Integer and
binary64-width boundary cases are prepared, not runtime-tested. The ledger
refresh argument uses its actual 64-bit kind. The 18-column census still
lacks Q/bulk and per-cell absolute/relative defects and locations required
for the final residual-history study. No synthetic rows count as model
census evidence, and no second classification formula was added.

Authoritative verdicts are KN_CS_LS and KN_CS_REA from KN_CERTIFY, with
KN_PROD_MASK applied to the aggregate. Export all ray rows including probes,
then count only mask-true rows. Confirm the actual mask and active count;
do not assume that ray 83 is excluded in every configuration or conflate
84 solver rays with 84 production-mask entries. Keep ITLIM and DIV codes
even if absent in the historical 59 FULL /16 BULK /8 STAG /1 INV snapshot.
Reconcile every new refresh to its own OUTGEN, not to that historical total.

The existing format at fg_j_cmf_nonmono_v1.f:9410-9416 is defective: after
three real fields it expects I4, but the next argument KN_PA_RES is real.
An armed run may terminate at its first census. Existing binary was not
replaced. Safe isolated repair/build is required before relying on capture.
The table also lacks Q/bulk residuals, maximum absolute/relative defects,
argmax depth/frequency, cap flags and per-cell histories of changing rays.

Required analysis remains: complete per-refresh ID uniqueness and counts;
classification matrix; per-ray timeline; ever-BULK/STAG/INV lists; six-class
transition matrix; residual summaries; at least three consecutive selected
cell states for each changing ray. Missing diagnostics must be marked
unavailable, never inferred from aggregate counts. Counts varying by three
do not demonstrate a fixed marginal set or identify its cause. BULK is not
FULL and does not satisfy FG_RF_ALL_CONV.

Current scripts collide in working directory, scratch, prepared inputs and
append files. User approved preparing isolation, not scheduler changes.
User subsequently approved `scontrol hold 13464217 13464311`. Both were
rechecked PENDING, held, and verified JobHeldUser. No production job changed.
No cancel/release/resubmission is approved. The extra census I4 descriptor
is now removed in source; ifort syntax check passed. Armed binary unchanged.
