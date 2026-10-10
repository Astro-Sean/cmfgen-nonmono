# Adaptive block Gauss-Seidel design for the Knop wavelength coupling

2026-09-28. Design doc for the next iteration stage of
FG_J_CMF_NONMONO_V1 / MOM_J_CMF_NONMONO_V1, following the review of
the refresh architecture and the observed residual floors.

## 1. What exists today (verified against the source)

FG side (plane/fg_j_cmf_nonmono_v1.f):

- FG_KNOP_REFRESH_V1 (line 4645): after the full frequency sweep,
  per-ray parallel solve.  Each ray LS gets KNOP_RAY_SOLVE_V1 (4816).
- KNOP_RAY_SOLVE_V1 owns a per-ray field FLD(NI,NCF) and iterates
  SWEEP_FLD for the I_M (inward) then I_P (outward) field.
- The frequency coupling is NEAREST-NEIGHBOUR (SRC_PRIME, 5376):
    GAM(I)>=0: SP_L(I) = SR_L + Q_L * FLD(I,L-1)   (expanding link)
    GAM(I)< 0: SP_L(I) = SR_L + |Q_L| * FLD(JPR,L+1) (compressing link)
  On an ascending pass the L-1 read is fresh (GS), the L+1 read is
  stale (Jacobi lag); the descending pass inverts this.  The
  alternating pair is symmetric GS on a chain graph.
- Q_L(I) = GAM/(DLN*CHI) after chi_hat = chi + |Q| absorption
  (BUILD_FREQ 5287): the per-link gain is ~|Q|/(chi+|Q|), so cells
  with |Q| ~ chi (strong compressions) carry near-unity propagation
  factors -- the observed marginal mode.
- SWEEP_FLD controls: KN_RAY_TOL=2e-3 max-overwrite residual;
  NSTALL>=2 same-direction stall exit (pass NPO>=8); bulk exit
  NBAD<=max(50,NTOT/1000) at NPO>=6; OMEGA_RAY adaptive damping
  (same-direction ratio keyed); divergence abort at 1e3*BEST_RES.
- WRITE_STORE merges FLD back into I_*_DEPTH_NEXT under INEXT_RELAX
  and counts >50% changes into FG_GS_NCONV (the gate metric).

MOM side (plane/mom_j_cmf_nonmono_v1.f): J_NEXT_FRQ(ND_NM,NCF) is
iterated globally per refresh pass; compressing reads use the
pass-start snapshot JNEXT_RF_READ (Jacobi) or live J_NEXT_FRQ on
descending passes (MOM_RF_GS_DESC).  Vector AA on the composite
same-parity map (commit 76d3534).

## 2. The actual bottleneck

The chain-coupling gain |Q|/(chi+|Q|) -> 1 at the compression ramp
cells.  On a chain, N sequential marginal links give effective
contraction ~(1-eps)^N per pass over the coupled run: with runs of
hundreds of strongly-coupled frequencies, the alternating global
pass pair propagates the fixed-point information only a few
frequencies per pass -- the pathological scaling Knop Paper II's
iterative GS coupling was designed to avoid, now reappearing
inside the refresh because the *interior* iteration is still a
global Jacobi-lagged pass pair.

Three consequences observed in OUTGEN:

  a. refresh residual floor ~2e3 (MOM) / store churn floor
     (LAM pre-refinement): a small set of compressing cells pins
     the global max residual.
  b. stall exit fires globally: a handful of marginal cells declare
     the whole 57k-frequency field "stalled", truncating passes that
     were still making progress elsewhere (the STALL vs NOT-YET-
     PROPAGATED distinction).
  c. the maturity gate reads the global counters, so one marginal
     pocket blocks BA indefinitely even when 99% of the field is
     converged.

## 3. Design: coupling-adaptive block GS (FG refresh first)

The frequency dependency graph is a chain (nearest neighbour), so
a block partition is mathematically exact: blocks interact only
through their boundary frequencies.  Per ray:

### 3.1 Coupling metric and block detection

Per-frequency coupling weight for the lagged direction:

    W_L = max over ray points I with GAM(I,LS)<0 of
          |PK_Q(I,L)| / max(PK_CHI threshold, tiny)

i.e. the largest upwind-coupling opacity ratio at frequency L on
this ray.  W_L ~ O(1) means the link is marginal; W_L << 1 means
the link is weak (error decays per hop).  PK_Q is already in the
pass-invariant pack, so W_L costs one max-reduction per frequency,
once per ray per refresh.

Blocks = maximal contiguous runs of {L: W_L > W_CRIT} dilated by
OVL frequencies on both sides (OVL ~ 4-8).  W_CRIT ~ 0.1-0.3.
Runs separated by < OVL merge.

Data-driven alternative/fallback: mark frequencies whose per-L max
cell residual exceeds KN_RAY_TOL after pass 2 (RES_FREQ tracking).
The W_L criterion is sharper because it predicts marginality before
the residuals prove it; RES_FREQ is the safety net for rays whose
W profile underestimates coupling.

### 3.2 Iteration schedule

Per ray, replacing the flat alternating pass pair:

    phase A: 2 global passes (ascend, descend) -- propagate boundary
             conditions, seed the field, build the hot-block set
    phase B: block refinement -- for each block (swept in the same
             alternating order), iterate ONLY the block's frequency
             range with local ascending/descending mini-passes until
             the block-local residual < KN_RAY_TOL or the block
             stalls; the block boundary read uses the neighbour's
             current FLD column (chain coupling => boundary column
             L=a-1 or L=b+1 is the only external input)
    phase C: 1-2 global passes -- re-verify weak links and update
             the global residual/conv counters the gate reads

Blocks iterate independently: nearest-neighbour coupling means
block interiors depend only on the two boundary columns, so inner
iterations are exact given fixed boundaries, and block-level GS
(alternating block order) propagates the between-block coupling.

### 3.3 Convergence criteria (two-tier, per user direction)

  KN_RAY_TOL   per-cell overwrite residual (unchanged, block-local
               now)
  KN_RF_PHYS   physical tolerance: exit when the cells that feed
               the store (J_PNT-mapped compressing points only)
               have converged below a looser threshold -- e.g.
               rel. change < 5e-3 -- since SE only ever reads those
               cells.  Cells with no store slot cannot influence
               populations; polishing them buys nothing.

Residual tracking gets a per-frequency array RES_FREQ(NCF) per ray
(~0.5 MB): max over cells at each L during MARCH.  Cost is one
extra MAX into a local per L inside the existing march loop.

### 3.4 What changes for I_NEXT / FLD / store / restart

  - I_*_DEPTH_NEXT store format: UNCHANGED.  WRITE_STORE is still
    the only writer; blocks change which passes run, not what is
    stored.
  - FLD layout: UNCHANGED (NI x NCF per ray per direction).  Blocks
    are index ranges, not copies; the boundary column read is a
    normal FLD(:,L) access.
  - FGINEXT_GS checkpoint format: UNCHANGED.  Block data is
    transient per-ray state; nothing crosses the requeue boundary.
  - MOM_J_CMF_NONMONO_V1: the same construction applies to
    J_NEXT_FRQ where a block is a frequency interval and "cells"
    are depths.  The refresh driver there is frequency-major with a
    pass-start snapshot; a block pass re-solves L in [a-OVL,b+OVL]
    against the snapshot at the boundary.  Deferred until the FG
    version proves the structure.

### 3.5 Why this helps (and what it does not do)

Block GS fixes allocation and measurement: work concentrates on
the strongly-coupled runs, and the stall/converged verdict becomes
local, so one marginal pocket can no longer pin the global field
state.  It also creates the natural container for stronger local
methods: a ~500-frequency x NI block is small enough for block
vector-AA or even a local Newton solve on the block residual if GS
stays marginal inside it.

It does NOT by itself change the contraction factor of the
compressing links -- it isolates them where they can be treated
(and measured) honestly.

## 4. Implementation stages

Stage 1 (instrumentation, low risk, unblocks everything):
  - Add RES_FREQ(NCF) per-ray residual tracking in MARCH/SWEEP_FLD.
  - Report top-K hot frequencies per ray + per-block NBAD in the
    refresh summary.  This alone tells us whether the marginal set
    is contiguous (chain blocks) or scattered.

Stage 2 (schedule, moderate):
  - W_L computation from PK_Q in the pack build.
  - Block detection + phase A/B/C schedule in SWEEP_FLD: the L-loop
    bounds become a block range for inner iterations; boundary
    columns read current FLD.
  - Local NSTALL per block; global stall counter only over weak
    links.

Stage 3 (criteria):
  - Store-feeding subset convergence test (KN_RF_PHYS) wired into
    the maturity getters, so the gate measures what the SE solve
    can actually see.

Stage 4 (optional, only if blocks still stall):
  - Per-block vector Anderson on the block's compressing-cell slice
    (reuse the AA machinery from mom_j_cmf_nonmono_v1.f, applied to
    a NI x B-sized vector -- small Gram systems).
  - Last resort: implicit block solve (dense Newton on the block's
    nearest-neighbour chain; NI*B unknowns).

## 5. MOM-side analogue

J_NEXT_FRQ's lagged compressing read is the same chain coupling in
frequency, averaged over rays.  The block idea transfers as
frequency-interval iteration of the moment solve over [a,b] with
JNEXT_RF_READ boundary values; per-block AA is already the right
shape (the AA Gram system is small).  No changes until the FG
prototype validates the block semantics.

## 6. Test plan

  - LAM 12007689 (running, refined ND=201 grid): the gate now prints
    churn diagnostics; compare store cells>50% decay vs the 0.33%
    floor of the 70-pt run.
  - After Stage 1 lands: one LAM run gives the RES_FREQ histogram ->
    confirms whether the marginal set is a contiguous block (expected:
    frequencies whose line cores cross the sigma=-120 ramp).
  - Stage 2 acceptance: refresh residual and FG_GS_NBAD per block;
    global stall must no longer fire while a hot block is still
    contracting.
