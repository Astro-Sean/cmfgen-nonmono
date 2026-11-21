# Iteration-map characterization - limited replay evidence

Provenance: historical KN_REPLAY_HST_L081/L082/L084.bin caller-loop
histories; source stage audit and inspection commands in the parent report.
No new history was generated in this continuation. Existing numerical
summary in map_characterization.json is retained, not revalidated here.

Historical D1 relative-L2 first/last: 81 1.33e-2/1.71e-10;
82 1.23e-2/3.43e-10. Median gain about 0.64; max gains about 6.5 and 9.3.
These are finite-history measurements of repeated complete one-pass ray
solves with caller store writes. They are NOT measurements of uninterrupted
production sweeps, the complete refresh operator or population feedback.
The caller writes KR_DBG_FLDP into P store slots, bypassing the usual
capped-commit trajectory; each full solve recomputes M and resets local
controller state. KR_DBG fields are solve outputs, not committed stores.

D1 decay in that sequence is evidence of decreasing increments there,
not a proof of contraction of the production map. A finite gain above one
does not identify non-normality. Neither a spectral radius nor local
lambda=1 follows from these norms or sign counts. The earlier bright-tail
near-neutral and guaranteed linear-convergence claims are withdrawn.
A convincing period-2 history is absent; KN_ORD_P2DAMP stays disabled.

Pending production-history analysis must declare state/stage, full versus
bright subset, denominator floors, signed D1/D2 alignment and absolute L2
as well as relative norms. G ratios based on changing normalizations must
not be presented as an operator-norm estimate. Relate threshold crossings
to actual per-cell histories and classifier inputs before attributing stalls.

Interpretation: Saad (2003), Iterative Methods for Sparse Linear Systems,
for iterate increments versus residual convergence; Trefethen & Embree
(2005), Spectra and Pseudospectra, for transient amplification without
inferring spectral information from a short nonlinear history.
