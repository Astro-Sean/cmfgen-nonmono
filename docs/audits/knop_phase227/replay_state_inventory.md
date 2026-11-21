# Phase 2.27B replay state inventory

Status: SOURCE-REVIEWED; lifecycle changes are not runtime-validated.
All source locations below are routine/declaration anchors in
plane/fg_j_cmf_nonmono_v1.f. Exact preserved source and pre-edit tracked
diff are under the external preservation_v2 directory named in
capture_isolation_design.md. Final artifact hashes belong in the build and
completion manifests, not inferred from the source commit alone.

## Scope and entry sequence

FG_KNOP_REFRESH_V1 resets refresh counters and hot-frequency spending,
builds coefficients and sibling scales, establishes KN_PROD_MASK, then
calls KNOP_RAY_SOLVE_V1 once per ray in its parallel ray loop. Each ray
builds geometry/donor lookups and optional static coefficient packs once,
seeds and solves M once, retains IM_NI=FLD(NI,:), seeds P once, solves P,
verifies copies of both fields, classifies, and commits/rejects stores.
The refresh pools residual samples, updates generation-level INEXT_RELAX,
and tallies the authoritative certification census and all-ray verdict.

The new omf path uses that SAME single ray call. It does not implement a
new state container, resume after termination, or copy arbitrary work
arrays between separate solves. Persistence follows lexical production
lifetimes. These facts are source-reviewed, not a demonstrated replay match.

## Semantically persistent state

| State / declaration anchor | Initialization | Updates and lifetime | Required across passes / old omf |
|---|---|---|---|
| FLD, FLD_P (ray allocatables) | SEED_FLD once per channel | in-place march, acceleration/rollback; through verification and commit | yes / old omf reconstructed P through sparse store and recomputed M |
| IM_NI (ray allocatable) | M sweep result at inner boundary | held fixed through all P passes | yes / old omf recomputed each call |
| OMEGA_RAY, RES_PREV(2), RATIO_OM (ray host) | each channel sweep head | controller at pass end; scalar omega applied at next pass's cell updates | yes / old omf reset |
| NPO, NSTALL, NOSC, NPLAT, L_PLAT, BEST_RES/BEST_NB/BEST_NT, RATIO_PV (SWEEP_FLD_ORD locals) | sweep head | pass accounting, best-pass and stopping tests | yes / old omf reset |
| KN_SGN_SUM, KN_SGN_PREV (ordered locals); SGN_ACC, SGN_PRV, NOSC_FL (ray host, alternating path) | sweep head | sign accumulation and oscillation detection | previous sum/counter yes; current sum resets / old omf reset all |
| RHO_WIN, RHO_SRT, NRHO/JW/KW (sweep locals) | ring initialized at sweep head | window persists; sorted scratch assembled for export | ring yes, sorting scratch no / old omf reset |
| W_LO, W_HI, NHOT (ray host); HOTMODE (sweep local) | channel setup/sweep head | hot-window activation, retained restriction and spending | yes / old omf reset |
| KN_HOT_CNT, KN_HOT_SPENT (module) | refresh head | shared frequency activity and budget; OMP protected | refresh-wide yes / old omf zeroed counters and lost cross-ray context |
| FLD_PRV, DLT_PRV (ordered allocatables) | allocated if Aitken enabled | previous iterate/delta history | conditional yes / old omf reset |
| AA_XH, AA_RH, NH_AA, AA_REF, AA_PEND, AA_OFF, DO_AA (ordered locals) | sweep head/lazy history allocation | secant history, acceptance/rollback and disable state | conditional yes / old omf reset |
| GRM_AA, WGR_AA, BRH_AA, WBR_AA, GAM_AA, IPIV_AA and norm/step scalars | AA update | recomputed solve scratch | no independent persistence requirement |
| P2_SGN, P2_FLIP, P2_DMP (ordered allocatables) | optional allocation at sweep head | per-cell mode and sign histories | conditional yes / old omf reset; P2 remains disabled |
| CLEAN_EXIT, DIVERGED, INVR, IXIT_M/P, RES_M/P, NP_M/P (ray host) | ray/channel entry | sticky invalid reason, solve verdicts and residuals | yes / old omf reinitialized |
| Q_M/P, FLD_Q, WRC_SV, NCLS_SV, RHOS_SV, flags/counter saves (ray host) | ray defaults; copies at verify entry | scratch verification saves/restores production state | verify-lifetime only / old omf invoked verification every full call |
| ICLS, IREA and KN_CS_LS/REA | INV defaults; KN_CERTIFY | authoritative verdict before commit | retained output; not a substitute for solver state |
| WRC_I/L/PI/PL/SAME, WRC_IM/LM/SM, WRC_ADEL | sweep/ray head | worst-cell location and persistence, preserved around Q | yes where previous-pass comparison uses them / old omf reset |
| NCAP/NNEG, NCP_P/NNG_P, NCL_M/P, NNL_M/P, NCP_C, NBAD_CUR/NTOT_CUR, NB_M/P, NT_M/P, NB_BM/BP, NT_BM/BP, NP2_CNT | ray or pass head as named | pass counters reset; cumulative/terminal outputs retained | cumulative and terminal counts yes; per-pass accumulators recomputed |
| I_M/P_DEPTH_NEXT, FG_RAY_SEEDED, FG_RAY_VC/DCH, JMP_SIB_M/P (module) | restart/refresh setup | commit policy, retained generations and immutable sibling reference | refresh context yes / old omf bypassed P cap with caller overwrite |

## Correctly recomputed state

- IKC1/2, JPR1/2, WINX, IGRD, NGRD, DRF_L: initialized in the ray's
  geometry loop; invariant for that ray solve.
- PK_DTAU/A0/A1/A2/A3/A4/SR/Q, PK_ETANI/CHINI: built once from static
  captured opacities, geometry and boundary inputs; shared by M and P.
  HAVE_PACK/ORD_SWEEP determine the actual branch. A requested ordered
  flag does not itself prove the ordered branch was taken.
- CHI_L/ETA_L/DCHR_L/Q_L/SR_L/SP_L/S_L/DS_L/DTAU_L/A*L/E*L, NFM_L,
  DLNM_L: coefficient/interpolation scratch; rebuilt or unpacked as the
  original implementation requires. No snapshot between passes needed.
- SPW four-row source ring: refilled by KN_ORD_SPW_FILL according to row
  order; retained within the live sweep, not externally reassembled.
- VMAXF, RES_FREQ, RES, NBAD_CUR/NTOT_CUR, KN_SGN_SUM: recomputed/reset
  per pass. VMAXF must use the live pass input, not an end-state peak.
- PSD_X/PSD_REC/PSD_N and NOP diagnostic state affect reporting, not the
  physical map. Diagnostic collection lifetimes still stay uninterrupted.

## Module/SAVE/COMMON audit boundary

Module variables are implicitly saved. The ray implementation has local
host work and contained sweeps; no new COMMON or SAVE state was introduced.
The new diagnostic schedule arrays/counts and selectors are module state,
read only for one explicitly selected replay ray. They are not a parallel
multi-ray schedule service. Production leaves schedules unallocated; stage
writes are restricted to existing dump-ray selections or explicit replay
selection. Existing probe/seed diagnostic arms must be off or separately
captured; arbitrary module-state equivalence is not established by v2.

## Solve-local controller state (2.27D audit, SWEEP_FLD_ORD)

All initialized at sweep entry (~12288-12310), live only for the call -
no SAVE attribute. Single-call replay preserves them exactly as
production does; repeated entry would silently reset every row below.

| Variable | Init | Update | Lifetime | Production | omf replay |
|---|---|---|---|---|---|
| OMEGA_RAY | =KN_RAY_OMEGA seed | *0.7/*1.15 clamped [0.45,1.0] at pass end | whole channel sweep | applied+chain | applied overridden per pass; chain re-anchored |
| RATIO_OM | 0 | residual ratio per pass | sweep | controller input | replay's own residuals |
| RES_PREV(2) | 0 | residual history | sweep | controller input | same structure, own values |
| NOSC | reset on monotone | +/-1 on sign-flip streak | sweep | stopping input | own values |
| KN_SGN_SUM/PREV | 0 | sign accumulator | sweep | oscillation detect | own values |
| BEST_RES/NB/NT | HUGE/0 | best-pass tracking | sweep | best-iterate commit | own values |
| NPO/NB_O/NT_O | 0 | pass/depth counters | sweep | termination | identical machinery |
| IXIT | - | exit class at stops | sweep | classification feed | identical machinery |
| L_CONT/NP_EFF | from KN_CONT_NP | continuation budget | sweep | pass budget | identical |
| NHOT/HOTMODE | per refresh | hot-pass spend | refresh-global | shared budget | replay: fresh (no refresh context) |
| AA_XH/AA_RH | lazy at NPO>=4 | Aitken history | sweep | AA disabled | disabled, unused |
| kn_replay_omega | deck schedule | read-only | driver scope | n/a | injected at pass top |

Module-level persistent inputs the driver pokes from the deck:
I_P_DEPTH_NEXT/I_M_DEPTH_NEXT stores, IM_NI, seeded P field, INEXT_RELAX,
KN_CONT_NP. Module init is verified; live-vs-deck equality of IM_NI and
the M solution is NOT - replay recomputes M.

## Diagnostic identity lifetime

KN_DIAG_REFRESH is an integer(kind=8) module variable, initialized to zero
per process and incremented only at eligible armed FG_KNOP_REFRESH_V1 entry.
It does not feed solver equations, policy or original ledger counters.
Stage/omega/cell/census records and deck filenames share this serial.
LEDG_N_RF is reset in LEDG_ITER_RESET (iter_ledger_mod.f90) and was not a
whole-run key. Old build_v3 and earlier appended telemetry could repeat
refresh identities across great iterations. Those candidates are retained,
not submitted. A resumed capture must use a fresh capture directory and
manifest; serials are unique within a process/capture, not across restarts.

## Known remaining fidelity gaps

The v2 deck starts immediately before P, after production M. Baseline
replay starts at ray entry and recomputes M from old stores. Its M output
must be compared against the deck's post-M field before interpreting P.
The deck does not provide production per-pass fields, shared hot-budget
history of other rays, or a complete general module-state checkpoint.
Single-call persistence repairs the old lifecycle reset; it does not by
itself close those input/context gaps. KN_CONT_NP and all stopping rules
are unchanged. Replays fail rather than fabricate additional scheduled
passes when the live solver terminates at a different pass count.
