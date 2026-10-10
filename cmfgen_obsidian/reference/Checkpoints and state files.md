# Checkpoints and state files

What carries state across requeues, and what each file actually
contains. Critical because SLURM requeues kill the process - anything
not on disk is re-zeroed (and the lagged radiation store is the
classic requeue wreck).

## Warm-restart set

| File | Writer | Contents |
|---|---|---|
| `POINT1` / `POINT2` | `SCR_RITE_V2` / `SCR_READ_V2` | R, V, SIGMA, POPS (the classic checkpoint pair) |
| `FGINEXT_GS` | `DUMP_FG_INEXT_GS` (inside SCR_RITE_V2) | FG lagged `I_NEXT` store + fingerprint (VER=3) |
| `MOMJNXT_GS` | `DUMP_MOM_JNEXT_GS` | MOM lagged `J_NEXT_FRQ`/`H_NEXT_FRQ` (VER 3; older cleanly rejected -> Lambda warmup rebuild) |
| `BAMATPNT` (`SCRATCH/`) | `STORE_BA_DATA_V3` | stored BA matrix pointer (frozen-Jacobian reuse) |
| `EDDFACTOR` / `EDDFACTOR_INFO` | | persisted Eddington factors; the `COMP_F` flip on requeue reuses them |
| `POP<SPECIES>` | at CMFGEN completion | per-species populations - **required by cmf_flux**; only written at completion, so an unfinished run resubmitted OBS-only can lack them |
| `RVTJ` | | final R,V,T,J dump |
| `GS_STATE` | SOLVE_FOR_POPS each call | freeze state for the slurm wrapper (`n_frozen` first key, `n_trial_acc` persisted) |

## Fingerprinting

- `FGINEXT_GS` restore validates bitwise R/V/SIGMA + `COMP_DEPTH_IDX`
  match; deliberately strict. On grid change the store is remapped
  (checkpoint VER=3). A failed restore leaves the store zeroed -
  covered by `N_LAM_RWARM` warmup.
- `GS_RWARM_SKIP` requires BOTH the FG I_NEXT and (when the MOM solver
  ran) MOM J_NEXT checkpoints mature: store restored, last sweep <5%,
  no capped writes.

## History lessons (why the preserve list is what it is)

- `run_model_grid.py`'s resubmit wipe once dropped POINT1/2 + POP* ->
  forced cold starts and an INDEXX segfault on OBS (`N_LINE_FREQ=0`
  with no POP files). Now preserved.
- `FGINEXT_GS`/`MOMJNXT_GS` were initially missing from the preserve
  list - every resubmission cold-started the coupling and convergence
  could never accumulate through the workflow.
- POINT1/2 do **not** contain the GS radiation state (despite a code
  comment claiming so).
- During JFNK evals, `FREEZE_I_NEXT` + the ACCESS_F cursor reset +
  `JFNK_IN_EVAL` keep residual/JVP evals from mutating the lagged
  store or checkpoints.

## Related

- [[Model pipeline]] - requeue mechanics, COMP_F flip
- [[Field maturity and GS freeze]] - warmup rules around the store
- [[Build and git workflow]] - run products never enter the git repo
