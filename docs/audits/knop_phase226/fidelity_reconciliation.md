# Replay-vs-production discrepancy reconciliation (Phase 2.26 sec 2)

Comparisons are on DIFFERENT array pairs - do not merge them:

A. replay end field vs deck FLDP_END (production end of the SAME
   refresh): max ABS = 2.4E-4 (82P), 2.7E-4 (81P).
B. median RELATIVE diff of A over cells I>1E-9: 3.7% (82), 3.3% (81),
   roughly uniform across rows incl. the M-fed tangent row.
   A and B are consistent: an absolute-diff floor of ~1E-4 gives
   ~4% median relative on cells of typical I~1E-3.
C. seed: replay seeds from the deck's own FLDP_SEED record ->
   identical to the production starting state BY CONSTRUCTION
   (0 difference; verified: rec25 = production seed).
D. end vs end = A/B above.
E. v1 reconstruction: KN_CHI_EXT/KN_ETA_EXT absent (FBITS bit9=0);
   replay rebuilds node values EXP(COEF(K,4)) - exact IF coefs and
   EXT are same-generation. The M-channel solve itself already
   differs from production FLD(M) by up to ~1.7% rel (row3:
   2.2E-4/1.25E-2; row4: 8.8E-5/1.7E-2), and the M field depends on
   EXT/geometry only - so an input discrepancy, not the P iteration,
   produces the systematic offset. If the coefs are one refresh
   older than the EXT used on nodes, the replay solves a slightly
   different map: consistent with a uniform ~3-4% field offset and
   with residuals/ixit still matching (the stall is a map property,
   not a state property).

Remaining candidate causes, ranked:
 1. coef-vs-EXT generation mismatch (fits uniform offset + M-only
    evidence; UNVERIFIED - needs a v2 deck with real EXT).
 2. FP/compiler drift (-O0 driver vs -O3 production): explains
    scattered noise, NOT a uniform 3.5% median.
 3. Boundary handling: FIXED this phase (core-index bug); residual
    diff after fix is 2.4E-4 abs, so BC is no longer the suspect.
 4. Different initial conditions: excluded - seed is read from the
    deck (C=0).
 5. Iteration stage: replay runs a FRESH solve from the seed (not a
    continuation), so npo differs (25 vs 15); end-state comparison
    is still valid for the same seeded map.

Not the cause: geometry (identical arrays verified during parse),
deck corruption, driver record misalignment (all reads
size-verified: I_M/I_P_DEPTH_NEXT are REAL*8 coupling stores,
31.5MB each, read exactly).
