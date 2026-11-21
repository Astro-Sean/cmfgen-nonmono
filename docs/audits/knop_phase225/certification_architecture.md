# Certification chain (current source)

1 ray solver:   SWEEP_FLD_ORD exit IXIT (1 tol,2 bulk,3 stall,...) ->
                KN_CS_LS class FULL/BULK/ITLIM/STAG/DIV/INV
                (fg_j_cmf_nonmono_v1.f ~9365 census).
2 refresh:      FG_RF_ALL_CONV = .NOT.FG_RF_FAILED .AND.
                KN_NFULL.EQ.KN_NACT .AND. KN_NACT>0 (~9381).
                BULK rays do NOT count as FULL.
                FG_COMMIT_NC caps provisional store writes for
                non-FULL rays (WRITE_STORE_CAP ~14420).
3 field ready:  GS_RF_ALLCONV -> LAM_FIELD_NOT_READY=T when F
                (solve_for_pops.f ~3507-3527, reason 'rf_not_allconv').
4 populations:  CONV_EXIT_MET requires .NOT.LAM_FIELD_NOT_READY
                (cmfgen_sub.f ~5081, ~5513).
5 spectrum:     knop_spec_gate.f classifies GS_STATE/manifest
                before CMF_FLUX_SUB_V5 writes OBSFLUX.

Insertion point for a separate physical-significance verdict:
between 2 and 3 - a new per-refresh flag (e.g. FG_RF_PHYS_OK)
computed from measured rate sensitivity, consumed ONLY by the
LAM_FIELD_NOT_READY fold-in as an alternative to GS_RF_ALLCONV,
with KN_CS_* and FG_RF_ALL_CONV untouched and both flags written
to GS_STATE/OUTGEN. Not implemented (no source change authorized).

Latest diagdump_i census (OUTGEN line 22835): full=59 bulk=16
itlim=0 stag=8 div=0 inv=1, rf_all_conv=F, active=84; stag
fluctuates 8-10 between refreshes. Last PSD solve per channel:
ixit=3 only on 81P/82P/83P; ixit=2 on 24M-26M,68P-80P.
