# Phase 2.23 residual/metric definitions (source: plane/fg_j_cmf_nonmono_v1.f)

Four distinct quantities. Do not use one as a proxy for another.

## A. Per-pass transfer residual (RES / RES_FREQ / dmax in PSD records)
- Location: MARCH inner loops, ~12214-12229 (P chan), ~12336 (M chan);
  SWEEP_FLD_ORD accumulate, ~12789+; SWEEP_FLD equivalent.
- Def: DCELL = |VNEW - VOLD| / MAX(VOLD, 1E-4*VMAXF(L)); RES=MAX over
  cells in the pass. VMAXF(L)=MAX(MAXVAL(|FLD(1:NI,L)|),1E-25): live-field
  per-frequency column peak floored at 1E-25.
- VNEW sanitized before the defect: <=0/NaN -> 1E-30 (NNEG), >CAP -> CAP.
- Measured BEFORE the relaxed commit FLD=VOLD+OMEGA*(VNEW-VOLD): it is the
  defect of the pre-relaxation map output vs the previous iterate.
- Significance gate (NTOT/NBAD) applies AFTER RES: floor cells inflate RES
  but do not enter the bad-cell census.

## B. Independent fixed-point residual (KN_Q_RESID, q_M/q_P)
- Location: verify pass ~10596-10695; KN_Q_ACCUM ~12222.
- One extra pass on a scratch copy (L_QP=T), before commit: ||T(x)-x||
  measured with the same DCELL normalization. Never writes production
  stores. Runs BEFORE KN_RES_PASS / relaxation-history updates.
- q census: per-channel sigmax (max significant-cell defect), absmax
  (max |DV|), flrmax (max floor-domain defect). absmax<=1.5E-4 on A2:
  the map is converged in intensity units; the tail is denominator noise.

## C. Store-update metrics (FG_GS_MAXCHG / GS_RF_DCHG / FG_RF_PEND)
- Location: in-sweep store comparison ~6040-6189; retained-delta export
  ~9853-9891; argmax ring ~6189+.
- FG_GS_MAXCHG: |fresh - lagged_store| / MAX(store, 1E-4*ray_peak) over the
  sweep; argmax ring carries (d,r,f,value,class bndry|cycle|margnl).
- GS_RF_DCHG/FG_RF_DCHG: retained-store change returned by WRITE_STORE /
  WRITE_STORE_CAP per refresh (committed generation only).
- FG_RF_PEND: pending delta from rejected/non-converged rays.

## D. Provisional-store cap (FG_COMMIT_NC / WRITE_STORE_CAP)
- Location: ~14420. DVB = FG_COMMIT_NC * MAX(V2, 1E-4*VMAX_F(LS2));
  |V1-V2|>DVB -> V1 = V2 + SIGN(DVB, V1-V2).
- Saturating commit returns exactly NC=0.05. It bounds store motion;
  it is not a convergence tolerance. Verdict codes: 0 committed,
  1 jump-reject, 2 nonconv-reject, 3 capped commit, 4 diverged,
  -1 unsolved, 7 probe-suppressed.
- Capped writes enter the retained store -> consumed by next SEED_FLD
  and by radiation-dependent rates.

## Controller values (deck-confirmed, rays 81/82)
KN_RAY_TOL=2E-3, KN_RAY_BULK_TOL=0.06, KN_RAY_OMEGA=1.0,
KN_RAY_STALL=0.98, NONMONO_INEXT_CAP=1E20, FG_COMMIT_NC=0.05,
FG_COMMIT_JUMP=50, FG_DIV_HEAL_CAP=1E6, INEXT_RELAX=0.5,
FG_ORD_P2N=2, FG_ORD_P2REL=3, KN_INNER_METH=DIFFUSION.
