# P-state attribution - gated on the M verdict

Status: BLOCKED-BY-ORDERING. Per the phase requirement, P-channel
interpretation waits until the M reconstruction verdict exists.

When M matches: first-divergence scan across P stages in order -
seeded field, PASS_START omega (proposed/applied), per-pass formal
result, relaxed field, march result, FLDP_END, verify, classification,
commit, next seed. Metrics per two_channel_replay_metrics.json; P and
M reported separately; KNAPPLY proposed-vs-applied supplies the
controller side-channel read (a proposed!=applied gap that leaves
applied==production is harmless unless it flips an exit decision).

If M does NOT match: P differences are inherited input error, not a P
defect; the fix is diagnosing why M recomputes differently (budget
channel vs missing module state vs FP), not P instrumentation. An
imposed-M experiment may then bound the downstream P contribution but
is labelled non-faithful.

FLDP_END stays pre-cap: no FLDP_END discrepancy may be attributed to
same-refresh commit capping.
