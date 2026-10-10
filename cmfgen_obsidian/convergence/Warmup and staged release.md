# Warmup and staged release

Start-up sequencing: what runs before the model is allowed to take
aggressive steps, and how physics is released in stages rather than
all at once.

## Cold start

- The model starts from `LTE_EST` populations at the hydro `T_gas`.
  On the Type IIn model that start is itself the first instability:
  the LTE state is inconsistent with the radiation field of a shocked,
  compressing structure, and with unlimited Lambda corrections
  (`LAM_SCALE_OPT=NONE`, `MAX_LAM=1E10`) it1 already explodes
  (max change ~6e16%, [[convergence_issues]] Sect. 6). The it2
  luminosity jump is the same transient on the stock control run -
  a setup effect, not a solver defect.
- The current production recipe deliberately keeps
  `LAM_SCALE_OPT=NONE` and `MAX_LAM=1e10`: a small Lambda cap makes
  the LTE-to-NLTE population cascade impractically slow. Safety comes
  from residual acceptance, `BA_ENTER_RAW`/`BA_SAT_RAW`, the
  `BA_DAMP_CAP` leash and `BA_PROBE_LAM` probes. The older suggestion
  of `LIMIT`/`MAX_LAM~1e2` is superseded for this model.
- Starting from a converged SOB population checkpoint is still the
  preferred staged route; it removes most of the cold-start transient
  before Knop/JFNK runs.
- `LSTAR`/`IT_ON_T` gotcha: `LSTAR=0` with `IT_ON_T=T` collapses T to
  the grey floor at model setup - keep `IT_ON_T=F` on ejecta runs
  ([[VADAT keys]]).

## Warmup windows

| Window | Key | What it protects |
|---|---|---|
| Per-process Lambda warmup | `N_LAM_RWARM` (4) | `I_NEXT` store is process-local and re-zeroed on requeue; a BA step on a zeroed store is the classic requeue wreck |
| JFNK admission warmup | `JFNK_WARMUP` (hyb: 6) | Newton not armed during early transients |
| Checkpoint skip | `GS_RWARM_SKIP` | a validated mature `FGINEXT_GS` skips re-warmup |

The warmup count is driven by `ITS_SINCE_START`, not the checkpointed
`MAIN_COUNTER`, matching the process-local lifetime of the field
store.

## Staged physics release

- `FIX_T=T` holds the hydro temperature structure; `DO_T_AUTO`
  releases it when `MAXCH < 50` once a non-Lambda step has run.
  `N_FULL_TREL` / `T_HELD_PHASE` phase the release over accepted full
  steps rather than flipping it in one go.
- `USE_FIXED_J` models run the early iterations on a fixed J and
  switch to variable J at `MAXCH < 50`; `DO_GT_AUTO` optionally
  inserts a `GREY_T_ITERATE` at the switch.
- `D2_EQ_D1` slaves the depth-2 populations to depth 1 on non-Lambda
  iterations - an outer-boundary helper while the inner boundary is
  still settling.
- Grid events force a Lambda recovery iteration and re-checkpoint:
  `REVISE_R`/`R_GRID_REVISED`, `DONE_HYDRO_REVISION`, clump revision,
  and the NG kick all set `LAMBDA_ITERATION=T, FIXED_T=T` for the
  next call.

## Runtime scheduling knobs

- `IN_ITS` (re-read each iteration): `NUM_ITS`, `DO_LAM_IT`,
  `DO_LAM_AUTO`, `DO_GT_AUTO`, `DO_T_AUTO`, `D2_EQ_D1`. Edit it live
  to extend a run or force Lambda.
- `IT_SPECIFIER` (optional file; excluded from pipeline run
  directories by default): rotating directives `LAM n`, `FULL`,
  `FULL(FIXT)`, `FULL(FIXBA)`, `DO_NG`, `DO_AV`, `REVISE_R`, `NORM`.
  See `subs/specify_it_cycle_v3.f`.

## Restart

- `FGINEXT_GS`/`MOMJNXT_GS` checkpoints preserve the lagged field
  store across requeue, and `POINT1/2`/`BAMATPNT` preserve
  populations and the BA matrix.
- What is **not** checkpointed: the maturity flags. A warm restart
  re-earns field maturity over ~25 iterations and the run cannot
  restart in BA mode - several hours of wall time per resubmit
  ([[convergence_issues]] issue 1.6).
- `GS_STATE` carries the stall-watchdog state and `N_TRIAL_ACC`.

## Related

- [[Convergence strategy]] - the staged plan these serve
- [[Lambda iteration]] - warmup pins this mode
- [[Field maturity and GS freeze]] - `N_LAM_RWARM`, `GS_RWARM_SKIP`
- [[Checkpoints and state files]] - what survives a requeue
- [[Acceleration methods]] - `IT_SPECIFIER` NG/AV directives
