# Full-refresh capture design (Phase 2.27F)

Status: designed; NOT submitted; readiness gates unchanged.

## Why the full 84 rays

Budget allocation is a population property: which sweep gets the
tail of the 480 depends on which rays are concurrently hot under
DYNAMIC scheduling. A 3-ray capture removes the competitors and
cannot measure production's allocation to 81/82.

## Capture matrix

- ALL rays: hot-budget event buffer (Part 5), census rows
  (class/IXIT/NP/NB/RES per channel), termination status.
- Rays 81/82/84 only: KN_RAY_DUMP decks, M_SOLVE_END + stage rows,
  KNOMG armed-ray omega schedule, FLDP_END.

## Provenance to record

exe hash (build_v5/v6 candidate), source commit + diff hash, input/
restart hashes per manifest, OMP_NUM_THREADS, OMP_SCHEDULE/DYNAMIC,
OMP thread count N_TH, node, refresh serial(s), MKL/OMP env, ray
completion order (budget buffer TID+seq gives it).

## Gates (unchanged, all currently NOT_RUN / BLOCKED)

read-set closure; input hashes; exe provenance; isolated paths;
runner resolves exact exe; Fortran census IO test; omega schedule IO
test; lifecycle review; refresh-serial test; ray-84 control; O0/O3
comparison; failure/heal census fields; no production change.

## Post-capture experiments

A fresh-budget, B recorded-schedule, C per-sweep hot-allowance
enforced, D (only if C insufficient) full ordered event replay.
First-divergence attribution order: M reconstruction -> hot budget
-> formal solution -> omega/controller -> stopping -> next-refresh
seed.
