# USE_LAM_ES=F WP-M5 - minimal fix and regression gate

Date: 27-Nov-2026. Evidence for `checklists/use_lam_es_eq_f.md`
WP-M5. Fixture: the committed WP-M1 `rfj` arm (NCF=32 build).
Every number below is from that fixture unless marked
production.

## The defect (WP-M2) and the fix pair

WP-M2 named the non-contractive ingredient: under
MOM_RF_OWNW+ALT, compressing rows write on descending passes but
read the frozen pass-start snapshot JNEXT_RF_READ
(:2838-2856, :6185-6186) - a pure Jacobi step every other pass.
The map's spectral radius is the census a (0.964 fixture /
0.995 production) and ~150-300 passes are needed against the
KN_M_MAXPASS=60 budget.

The fix is two existing/paired knobs, no algorithmic surgery:

1. **MOM_RF_GS_DESC=T** (existing VADAT key, default F): swaps
   the descending-pass compressing read from the frozen snapshot
   to the live J_NEXT_FRQ (:2838, :2906) - true Gauss-Seidel
   propagation along the red-to-blue compressing characteristic.
   The flag's own comment notes the cascade divergence that
   motivated the freeze predates the M-matrix/EPS/COH/HBC fixes.
2. **MOM_RF_RLX=1.0** (new VADAT key, default -1 = inherit
   JNEXT_RELAX): a refresh-only relax factor for the per-cell
   MOM_DAMP_NM schedule. JNEXT_RELAX=0.5 also damps the
   sweep-side J_NEXT write (:4280-4290, the SE<->RT loop), so it
   cannot be raised globally; under live-GS the contractive step
   is halved again by the Jacobi-calibrated under-relaxation.
   Sites: mom_j_cmf_nonmono_v1.f :398 (var), :5768 (resolved
   RLX_RF), :5956-5976 (init+cap), :6699 (Aitken exemption),
   :6972,7021 (damp diagnostics); setter SET_MOM_RF_RLX :5374;
   VADAT read rd_control_variables.f :1523-1532. Clamp (0,1]
   matches JNEXT_RELAX.

## Ablation on the committed fixture (ALT schedule, bulk-exit off)

| Arm | compressing | expanding | refresh | note |
|-----|-------------|-----------|---------|------|
| production defaults (ctrl) | STALL pass ~16 | survives | converged 54p | WP-M1 signature, amax 0.964 |
| MOM_RF_GS_DESC only | no stall | survives | NOT converged 60p, best-pass rejected | contracts ~0.8/cmp-pass, needs ~61+ |
| GS_DESC + MOM_RF_RLX=1.0 | no stall | survives | **converged 12p** | the fix pair |
| MOM_RF_RLX=1.0 only | STALL pass ~24 | stall | not converged | relax alone worsens the Jacobi map |
| GS_DESC + MOM_RF_RLX=0.7 | no stall | survives | converged 35p | margin exists |

Both converged arms store the identical field (RFST_DCHG =
6.02823747E-01 = the control's expanding-carried value) - the
fix reaches the same fixed point, not a different one.

## Secondary defects found along the way (recorded, not fixed)

- **Best-pass rank saturation**: once PASS_N50 and PASS_N5 both
  saturate at 0, the per-direction best-pass rank keeps the
  earliest pass (:6543-6560 ranks N50 then N5) - the GS_DESC-only
  run's composite "best" field carried an AA-spiked PASS_CHG of
  53.3 and hit the 1E3*MIN_RES runaway bound (:7204), discarding
  the whole refresh (kept prior, RFST_DCHG=0). The metric needs
  a PASS_CHG tiebreak when the bulk counters saturate.
- **Expanding-branch noise-floor stall**: the stall detector has
  no absolute floor (RATIO_M > 0.90 twice after N_FLOOR), so a
  settled field drifting ~0.01-1.8 eventually stalls too
  (AA-off runs; production never showed it). Candidate for a
  follow-up gating flag, not needed for the fix.
- GS_DESC alone leaves the compressing direction contracting but
  1 pass over budget on the fixture; production's rho~0.995
  needs the relax pairing, not a bigger budget.

## Regression gate (committed)

test_mom_solve.py, `rfj` arm, four assertions on
GS_DESC+RF_RLX=1 under the production-failing ALT schedule:
driver completes, neither direction stalls, RFST_CONV=1 with
RFST_DCHG > 1E-2 (a real converged field, not a best-pass
retention), and the refresh exits `converged [OK]` inside 60
passes. All 153 checks pass.

## What the fix does not claim

- It does not change the sweep write (JNEXT_RELAX stays 0.5
  there), the stall detector, the guard census, or the fixed
  point - same field, reached faster.
- It does not help the Jacobi path: MOM_RF_RLX alone makes the
  stall fire sooner. The pairing is the point: live reads make
  the map contractive, the refresh relax commits the step.
- Production efficacy is WP-M6's question: the fixture's a~0.96
  converged in 12 passes; production's a~0.995 has ~2x the
  contraction distance but the same mechanism, plus the
  sweep-side reseeding the refresh must outrun.

## Verdict

WP-M5 VERIFIED: the compressing march converges on the WP-M1
fixture under MOM_RF_GS_DESC=T + MOM_RF_RLX=1.0, gated behind
two flags with safe defaults (production path bit-identical:
GS_DESC stays F, MOM_RF_RLX=-1 inherits JNEXT_RELAX). The
implicitness improvement (live-GS read) is the fix; the
refresh-only relax is what lets it finish inside the pass
budget. No tolerance or budget loosening.
