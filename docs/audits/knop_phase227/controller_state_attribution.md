# Controller-state attribution - applied omega vs controller state

Status: SOURCE-REVIEWED. No runtime comparison exists; all controller
divergence claims below are mechanism statements, not measurements.

## Where each quantity lives (SWEEP_FLD_ORD, fg_j_cmf_nonmono_v1.f)

- Pass top (NPO incremented): `omega_proposed = OMEGA_RAY` records the
  value carried in - i.e. the controller's current chain value, which is
  exactly what production would apply this pass (no separate proposal
  variable exists; the controller updates OMEGA_RAY in place at pass end).
- Injection: `call kn_replay_apply(ls, merge(2,1,inward), npo, omega_ray)`
  overwrites OMEGA_RAY with `kn_replay_omega(npo, channel)` - only when a
  KNOMG2 schedule is loaded, only for `ls == kn_replay_ls`, only
  `.not. l_qp` (the KN_Q_RESID verify pass is outside the iterated
  sequence and correctly receives no scheduled row).
- `omg_used = omega_ray` then snapshots the actually-applied value and
  `kn_omega_applied(...)` writes (proposed, applied) to KNAPPLY.
- Pass end: controller multiplies OMEGA_RAY in place:
  `MAX(0.45, OMEGA_RAY*0.7)` on growth/oscillation, `MIN(1.0, *1.15)` on
  contraction. KNOMG row writes (NPO, OMG_USED, OMEGA_RAY, RES, RATIO_OM,
  NOSC) - i.e. applied value, post-update chain value, and the residual
  state that drove the update.

## What injection does and does not control

OMEGA_RAY is BOTH the applied relaxation and the controller's
multiplicative chain state. Injecting at pass top therefore:

- sets this pass's applied value to production's recorded value - the
  intended semantics;
- re-anchors the chain each pass: the replay controller's next proposal
  becomes `injected * f(replay_residual)`, so the omega chain can only
  diverge through residual differences, never by free-running drift.

Controller state NOT governed by injection (driven by replay's own
residuals): RATIO_OM, NOSC, KN_SGN_SUM/PREV, RES_PREV, BEST_RES/NB/NT,
WRC_* counters, hot-pass accounting, and every stopping/exit decision
(IXIT, bulk/stall/cap/plateau exits). These CAN diverge from production
if the trajectory diverges for any other reason (e.g. the recomputed M
state). A divergent NOSC/RATIO cannot corrupt the applied omega (it is
overridden anyway) but CAN change an exit decision - the dominant
controller-side failure mode for Mode B is an IXIT production never took
(or skipped one it did), which then truncates or extends the pass count
and trips the strict `kn_replay_used /= kn_replay_count` check.

## Measurement plan (blocked on runtime)

Production KNAPPLY rows show proposed == applied. Replay Mode B rows show
the replay controller's free-running proposal vs the injected applied
value. Mode C therefore needs no new instrumentation:

- first pass where replay `proposed` diverges from production `applied`
  (beyond bit-identical equality) marks controller-state divergence;
- KNOMG (RES, RATIO_OM, NOSC) columns per pass give the first residual
  history divergence, which precedes any omega-chain divergence;
- IXIT/exit-class mismatch is the only channel through which divergent
  controller state can alter the trajectory.

Mode D (restoring production controller internals) is not implemented
and stays unimplemented unless Mode C shows a controller-state-driven
trajectory difference that omega injection alone cannot explain.
