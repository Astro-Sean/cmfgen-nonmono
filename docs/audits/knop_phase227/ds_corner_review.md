# DS-corner review - hypothesis, not causal attribution

Provenance and inspected source: parent Phase 2.27 report. No diagnostic
measurement or equation change was made. Knop, Hauschildt & Baron (2009),
Paper II, A&A 501, 813-820, doi:10.1051/0004-6361/200911793, is the formal
solution reference; earlier equation-number assertions have not been
independently verified against the paper in this continuation.

The source preview KN_ORD_SPW_FILL and KN_ORD_DS ring secants can couple a
not-yet-marched neighbor to the current row. This supports investigating
mixed old/fresh source dependencies. It does not establish that a pure
previous-pass source applies at every frequency or that this term explains
the observed replay spread. Frequency ordering and ring overwrite times
must be retained explicitly.

Minimum observational diagnostic: save all ring values and limiter branch
inputs actually used at the row update, coefficients multiplying DSD/DSU,
pre-update intensity and raw VNEW; later evaluate neighbor source on the
post-march state using PRIVATE scratch ring/arrays. Keep all other inputs
fixed for the shadow slope evaluation. Record row, frequency, direction,
pass and refresh, then estimate local deltaV = A3*deltaDSD+A4*deltaDSU.
Record limiter branch changes because the slope is piecewise defined.
Neither scratch evaluation nor ring reconstruction may modify the live
march or shared module state. Confirm thread ownership from allocation
and OpenMP scope, not from variable names alone.

This is a local conditional sensitivity, NOT a causal fraction of the
complete iterative residual and NOT a guaranteed lower bound: propagation,
feedback, cancellation and limiter switches can increase or decrease it.
Quantitative contribution remains unknown. No corner/donor alteration or
new DS instrumentation is implemented in this continuation.
