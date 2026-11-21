# Anderson acceleration - deferred

Provenance: source inspection and replay limitations in the parent Phase
2.27 report. No experiment or production setting change was made.

AA remains disabled until uninterrupted production-history replay and the
full per-refresh census are validated. Existing residual sequences do not
establish its safety or prove that the accelerated operator is the limiting
production operator. Earlier local-eigenvalue and guaranteed-acceleration
interpretations are withdrawn.

A later proposal must explicitly audit vector/residual definitions,
omega-dependent history rescaling, rank/conditioning safeguards, positivity,
rollback and independent Q verification. For the piecewise-affine update,
record changes in slope limiter branches, commit caps, active-cell masks
and certification classes. Specify when history is invalidated after such
switches and compare against a pinned non-AA baseline. Do not silently
change certification or cap policy for the experiment.

Acceptance must improve independent KN_Q_RESID, transfer residual, FULL
production-ray count, refresh cost and field stability, not merely reduce
committed-store motion. Population and observer checks remain separate.
No statistically supported period-2 history exists; KN_ORD_P2DAMP also
stays disabled. Any enablement requires a separately approved experiment.
