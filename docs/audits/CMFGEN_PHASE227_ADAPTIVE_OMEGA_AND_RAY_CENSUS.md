# Phase 2.27 - adaptive omega, ray census and update-state attribution

Status: PARTIAL / BLOCKED. No new capture data are available. Prior claims
that the omega-history question, production contraction class or complete
census were resolved are withdrawn. No convergence fix is proposed.

## Phase 2.27B preparation update (21-Nov-2026)

Phase remains PARTIAL. Both original jobs remain JobHeldUser; no scheduler
mutation or new submission occurred in Phase 2.27B. Installed production
and armed diagnostic executable hashes remain the values below.

Confirmed preparation results:

- Original dirty source, tracked diff, untracked replay source, job
  references and hashes were preserved outside the repository. Pre-existing
  audit documents were preserved again before documentation updates.
- Isolated build_v4 normal and debug candidates compile and link. The
  executable, replay and actual-Fortran IO harness in each variant share
  one revised solver object. Debug checks do not cover the unchanged
  archived objects. No installed module, library or executable was replaced.
- omf now calls the normal ray solver ONCE and overrides applied omega at
  the ordered pass entry, preserving solve-local state and retained M
  boundary data across P passes. This is a source-reviewed architecture
  repair, NOT a demonstrated production trajectory match.
- New canonical KNOMG2 schedules require BOTH channels, declared counts,
  exact identities/order and an applied binary64 bit check. Missing or
  ambiguous records fail; Q verification does not consume production rows.
- KNAPPLY records the proposed and actual applied values before the field
  updates, including terminal passes. Existing end-of-pass KNOMG alone is
  not a sufficient source for canonical schedules.
- KN_DIAG_REFRESH supplies a process-local unique diagnostic serial.
  LEDG_N_RF resets each great iteration and cannot identify whole-run
  appended blocks. Stage/omega/cell/census output and v2 deck filenames
  now share the diagnostic serial; deck record framing remains unchanged.
- Stage diagnostics distinguish global refresh-entry retained stores,
  per-channel seeds, pass inputs, selected-cell raw formal and actual
  relaxed results, full march results, FLDP_END, Q scratch states,
  authoritative classification and retained sparse-store commit stages.
  Sparse stores and full ray fields are different coordinate spaces.

Strong evidence remains limited to the source-level lifecycle comparison.
The v2 deck's INEXT_RELAX is retained. Its omg_out is now printed as
`deck_omega`; it is not the replay terminal controller state. FLDP_END
precedes Q/certification/commit, so same-refresh commit capping still cannot
explain an FLDP_END-to-FLDP_END discrepancy.

Unresolved: input/restart read-set closure of the prepared model snapshot;
all actual-writer/reader runtime tests; the revised ray-84 control and
bytewise O0/O3 fields; production versus replay M state before P;
per-pass 81/82 metrics; shared refresh-budget/context equivalence; stage
coverage on failure/heal paths; and the full 84-slot census with probe mask.
The census row writer is now independently callable for format tests, but
its 18-column sidecar still lacks several requested Q/bulk/per-cell defect
and location fields. A complete residual-history analysis cannot be
claimed from those rows alone.

Hypotheses about DS-corner causal contributions, omega-distribution effects
and subsequent-refresh commit-state effects remain hypotheses. AA and P2
stay disabled. Free-matter submission and convergence changes are deferred.
No physical equation, convergence threshold, commit cap, damping default,
FG_RF_ALL_CONV rule, population or observer solver was changed.

See `knop_phase227/capture_isolation_design.md`, `replay_state_inventory.md`,
`replay_lifecycle_comparison.md`, `runtime_test_results.md`, the manifest
schema and completion_status.json. Exact hashes and generation commands
are in the dedicated build/capture manifests and phase227b_artifacts.json.
Historical line numbers below refer to the pre-2.27B inspected source, not
the newly instrumented file; use the hashed source and routine anchors in
the new inventory for current locations.

## Phase 2.27C update (21-Nov-2026)

Status stays PARTIAL. No submissions; both jobs still JobHeldUser.
captures_v4 preserved unchanged; captures_v5 fixes two manifest gaps found
in the independent review: the CMF_KN_RF_MB environment pin and a missing
work/data runtime-output directory (all ./data OPENs found are writes).

Read-set audit (read_set_audit.md/json): the 251-file manifest covers all
233 launcher link names plus 19 core model inputs. Every literal OPEN name
resolves to a manifest name, a conditional restart file, or an output.
Remaining unproven: observed read set (no runtime trace), fresh-vs-warm
restart assumption, and species-dependent dynamic-name coverage.

Equivalence checklist added to replay_lifecycle_comparison.md: M-once,
boundary retention, controller and stopping ordering are confirmed in
source only; no runtime equivalence is claimed. census_failure_heal_design.md
records that per-refresh rows preserve fail/heal history (a final FULL
erases nothing) but the row still lacks defect location, first-failure
pass, and per-ray Q/bulk residual fields. capture_approval_checklist.json
gives the 13 machine-readable gates: overall BLOCKED; six gates await an
approved allocation; the census gate is partial.

## Phase 2.27D update (21-Nov-2026)

Status stays PARTIAL; no submissions, jobs still JobHeldUser, production
exe unchanged. Source audit of the uninterrupted omf lifecycle is
complete: one KNOP_RAY_SOLVE_V1 call per ray, M-before-P order identical,
the KNOMG2 schedule overwrites OMEGA_RAY at each pass top (the exact
application point), and the L_QP verify pass is correctly outside the
iterated sequence. Solve-local controller state (RATIO_OM, NOSC,
RES_PREV, BEST_*, IXIT machinery) has no SAVE attribute and is preserved
across passes precisely because of the single call.

New finding on controller semantics: OMEGA_RAY doubles as the applied
relaxation AND the controller's multiplicative chain state, so each
injection re-anchors the chain - the replay controller can only diverge
through its own residuals (NOSC/RATIO/exit decisions), never by free
drift, and the only trajectory-level effect is a divergent stopping
decision. KNAPPLY's proposed-vs-applied columns are the built-in Mode C
measure; Mode D (controller-state restore) stays unimplemented.

Ranked discrepancy hypotheses (first_divergence_report.md): R1 replay
recomputes M from retained stores while the deck holds post-M state
(leading), R2 omega schedule, R3 commit/seed state, R4 FP noise, R5
controller side-channel. No runtime evidence exists for any; every gate
in capture_approval_checklist.json is unchanged and the candidate remains
BLOCKED pending allocation.

## Phase 2.27E update (21-Nov-2026)

PARTIAL; jobs held; production exe unchanged. M-to-P audit complete:
IM_NI = FLD(NI,:) is the ONLY M-produced quantity P consumes, and the
v2 deck already records production's live value (KN_RAY_DUMP phase-1
fires exactly post-M/pre-P). Deck lacks: full M field, RES_M/IXIT_M/
NP_M (census sidecar carries them), M terminal omega (P reseeds anyway),
M pass schedule (armed M KNOMG rows exist, IDIR=2), and KN_HOT_SPENT.

New divergence channel found: KN_HOT_SPENT is a refresh-global shared
hot-pass budget (480 total, 120/ray) consumed by whichever rays stall
first under OpenMP - production's spend at a given ray is scheduling-
dependent and unreconstructable; replay starts at zero. Under Mode A a
stalling replay ray gets the full budget and can outlive production's
pass count; under Mode B the strict schedule length bounds it. This
targets exactly the stalled-ray population (81/82).

Added M_SOLVE_END stage row (diagnostic-only, syntax-checked, build_v5,
captures_v6): records fld(ni:ni,:) = IM_NI with RES_M and terminal M
omega. Cheapest next runtime check: replay-IM_NI vs deck-IM_NI bitwise +
distributional, then per-pass M KNOMG IDIR=2 join. P attribution is
deliberately gated on the M verdict. All comparisons remain NOT_RUN.

## Phase 2.27F update (21-Nov-2026)

PARTIAL; jobs held; production exe unchanged. KN_HOT_SPENT audit
complete (hot_budget_source_audit.md, hot_budget_openmp_audit.md):

- Refresh-global INTEGER*8, reset per refresh before the ray loop;
  preserved across JFNK evals. Per-sweep cap NHOT<=120 (reset before
  M and before P: <=240/ray), global cap 480, overshoot charged.
- Sync class B: consume+verdict atomic inside CRITICAL(KN_HOT_SPD);
  TRY_HOT entry test is an unsynchronised read - stale admission
  burns one unit, zero extra marches. Bounded race, no fix applied.
- SCHEDULE(DYNAMIC) over rays => per-sweep hot allowance is
  scheduling-dependent; production is an ensemble over thread
  schedules for contended hot rays. Replay can only target the
  captured realisation.
- Deck stores the limits only. Mode B is already schedule-bound
  (recorded pass count includes hot passes). Mode C needs per-sweep
  executed hot counts - sufficient, since the consume is atomic and
  monotone; global order need not be replayed. Telemetry design:
  per-thread buffers appended inside the existing critical, flushed
  post-region; not implemented.

## Phase 2.27G update (21-Nov-2026)

PARTIAL; jobs held; production exe unchanged; no source changes.
Deep audit (hot_budget_deep_audit.md): global budget contention
ELIMINATED by existing output - archived OUTGEN reports
hot-pass budget spent=11-29/480 on every refresh; the cap never
bound, so no cross-ray scheduling channel exists in this run.
New source facts: IXIT=4 conflates itlim/per-sweep-cap/global-
reject into KN_CS_ITLIM; hot passes drive the full controller on
window-restricted residuals (only RHO_WIN/KN_RES_PASS excluded);
a rejected consume exits before march/omega-inject/KNOMG row, so
armed schedules carry executed passes only. Stale-admission race
is real in source but required 480 - never reached. Divergence
hypotheses re-ranked: M-state reconstruction leads; budget dead.
2.27H: read-only census/replay-IM_NI comparison recommended.

## Provenance and preservation

Inspection: 2026-11-21T05:55:59Z (UTC). Commands: `squeue -j
13464217,13464311`, `sacct -j 13464217,13464311 -X`, `scontrol show job
13464217`, `scontrol show job 13464311`, `git rev-parse HEAD`,
`sha256sum exe/cmfgen_dev.exe exe/cmfgen_dev_v2dump.exe`.

- HEAD: e28833fe9d1286c9b92f909f559d92d28f4f7cd1; dirty diagnostic source.
- Production SHA256: d50dcacd756088e3b61667edbaf10a0dfd80aab0ee5ac1a499e3c1605d18aba6.
- Armed diagnostic SHA256: 508e39c9e92729cd4383dd853505ae72bdd539bfbe5adfa7a2397dc24f9830ed.
- Both jobs: PENDING, AssocMaxJobsLimit, no start time, eight CPUs, 29G,
  partition interactive. Scheduler QOS is interactive for BOTH jobs;
  local scripts request normal/debug respectively. Neither is a debug-build
  comparison: both resolve `exe/cmfgen_dev_v2dump.exe` through batch_v2cap.sh.
- Both share the m15_2p82d_KNOP_diagdump_i/cmfgen/2.82days working directory,
  SCRATCH, prepared inputs, OUTGEN and append-only diagnostic files.
  This is an output-collision and input-mutation risk, not an isolated A/B.
- No job, raw output, executable or production setting was changed this pass.
  No KNOMG or KNRAY_CENSUS files were found in that directory.

After the initial inspection, the user explicitly approved holding both
pending jobs. Rechecked PENDING and ran `scontrol hold 13464217 13464311`;
subsequent squeue reports JobHeldUser for both. No production job changed.
Do not cancel or release either job without further approval.
Build and stage isolated diagnostic captures with
copied, hash-pinned executables and independent mutable inputs, scratch,
logs and output directories. Reconcile Slurm's spooled scripts with the
staged scripts before any release/replacement. Merely editing a local
submitted script does not update Slurm's stored copy. Input and checkpoint
hashes must be recorded before preparation and again at executable launch.

## Confirmed source findings

References below are to plane/fg_j_cmf_nonmono_v1.f in the inspected tree.

1. P seeding: SEED_FLD at 10565-10566. Deck phase 1 at 10579-10580
   records the P seed, old stores and the ALREADY-SOLVED inward field.
2. P solve: SWEEP_FLD_ORD at 10583. Applied in-cell relaxation is
   `VOLD + OMEGA_RAY*(VNEW-VOLD)` (e.g. 12367-12370).
3. Deck phase 2: 10607-10608 writes FLD_P and production exit data.
   This is BEFORE Q verification (10657 onward), KN_CERTIFY (10746),
   jump veto (10776) and WRITE_STORE_CAP (10854-10860).
4. Classification uses the same KN_CERTIFY verdict recorded in KN_CS_LS
   and KN_CS_REA. FULL commits through WRITE_STORE; non-FULL may take
   provisional WRITE_STORE_CAP. Classification precedes this commit.
5. WRITE_STORE sanitizes and bounds the candidate then writes it without
   INEXT_RELAX (14542-14584). WRITE_STORE_CAP additionally bounds change
   by FG_COMMIT_NC*max(old,1e-4*column_peak) (14656-14691).
6. KR_DBG_FLDP/FLDM copy the local solve fields at 10921-10925, not the
   modified retained stores. Commit capping changes STORE, not FLD_P.
7. The v2 scalar record DOES include INEXT_RELAX (11275-11279), read
   into INR by the driver. Earlier statements that it is absent are false.

Therefore a same-refresh cap cannot directly explain a metric comparing
KR_DBG_FLDP to deck FLDP_END: both are solve-stage fields. Caps can affect
future refresh seeds. The first differing operation is still unidentified.
The deck M field is post-M-solve, whereas replay runs M again: compare
that stage explicitly before attributing the P mismatch to omega or caps.

## Strong evidence, with limits

Historical printed baseline P maxdiffs: 81=2.656368e-4,
82=2.442209e-4, 84=3.487919e-8. Existing O0/O3 reports agree to printed
precision. A bytewise full-field comparison has not been demonstrated
here; compiler effects are not a credible explanation within those tested
outputs, but their global contribution is not measured as exactly zero.
Ray 84 does not bound floating-point sensitivity on rays 81/82.

The ray-82 seed-omega sweep leaves maxdiff nearly unchanged. This tests
four initial omega values, NOT a recorded per-pass history and NOT the
M-channel median/p90/p99. It does not establish schedule independence.
Moreover, the driver prints OMG_OUT read from the deck, not replay terminal
OMEGA_RAY. Its printed 0.450 values do not establish replay saturation;
that requires actual applied-omega telemetry.

The verify-mode scalar defects are historical diagnostics of its actual
reseed-and-solve operation. It recomputes M, maps P end values into store
slots, invokes normal solve/commit logic and can take more than one pass.
It is not yet a stage-validated unrelaxed application to the identical
full production state. Neither fixed-point error nor a cap cause follows
from those scalar values alone.

## Capture and replay defects blocking acceptance

- The armed binary's census format has an extra I4 edit descriptor before
  KN_PA_RES (a real). This can terminate a run at its first census.
  Removed that descriptor in source only; solver-source syntax check passed.
- Census omits per-channel Q, bulk residual, defect positions/frequencies,
  absolute defects, cap flags and selected-cell histories. Aggregate
  totals and every changing ray's three-refresh cell history cannot be
  supplied from that table alone.
- The armed binary's KNOMG uses ES14.7, insufficient for exact binary64
  round-trip omega. Source now uses ES26.17E3; syntax check passed, no
  executable replaced. Pass-top snapshots exist, but early exits precede the pass-end write;
  completeness must be checked, not presumed.
- omf performs repeated complete one-pass ray solves, resetting controller,
  M solution, seed/boundary and stopping state. That is not the uninterrupted
  production P pass loop. It cannot yet serve as a production-history replay.
- This pass hardens only replay input validation: longer paths, pass count
  from deck NPO, finite/range checks, contiguous P records, explicit failure
  on missing/ambiguous history, and no SIZE call on an unallocated array
  in non-omf modes. Syntax check passed with ifort 2021.11.1; two pre-existing
  unused-variable remarks (NDX, TMP). Existing binaries were not overwritten.
  Runtime malformed-history tests and full replay regressions remain pending
  isolated allocated execution, as required by project rules.

## Hypotheses and unresolved questions

Observed caller-loop D1 decay describes that replay operator only. It does
not prove production contraction, a spectral radius, non-normality or a
near-neutral bright tail. Finite gains above one can reflect nonlinear
switches and controller resets; they do not establish non-normality.

All per-refresh ray IDs, transition matrices, ever-BULK/STAG/INV lists and
cell-level causes remain unavailable. Do not force a new capture to match
the old final 59/16/8/1 totals; reconcile every refresh against its own
OUTGEN using KN_PROD_MASK, retaining all six class codes. Probe exclusion
must be taken from actual records, not inferred from ray number alone.

DS-corner source lag remains a causal hypothesis, not a measured residual
fraction or a proven lower bound. AA and period-2 damping remain disabled.
Free-matter experiment is NOT READY and requires explicit launch approval;
identical executables do not guarantee cancellation of nonlinear bias.
Mixed-sign active-cell validation remains unclaimed.

## Deferred work and acceptance

After safe isolation and a diagnostic-only repair/build: validate exact
application-point histories and uninterrupted-loop schedule injection;
record both channel fields/residuals at matched stages; measure maximum,
median, p90 and p99 with explicit floors and exclusions; independently
reconcile every census refresh; enumerate changing rays and obtain at least
three consecutive selected-cell records. Capture compiler/FMA/MKL/OpenMP,
input/checkpoint and executable provenance at launch. Run binary-field
O0/O3 comparisons, control-ray regression and independent v2 parser checks.
No acceptance item A-G is claimed complete from this pass. Restriction H
(no production physics/default/acceptance changes) is preserved.

References for interpretation: Knop, Hauschildt & Baron (2009), A&A 496,
295-298, doi:10.1051/0004-6361:200810290; Paper II, A&A 501, 813-820,
doi:10.1051/0004-6361/200911793; Saad (2003), Iterative Methods for Sparse
Linear Systems; Trefethen & Embree (2005), Spectra and Pseudospectra;
Hubeny & Mihalas (2014), Theory of Stellar Atmospheres, Ch. 13.
No unverified equation number is used as evidence.
