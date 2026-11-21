# M-to-P dependency audit (source-reviewed, KNOP_RAY_SOLVE_V1)

Traced at fg_j_cmf_nonmono_v1.f 10679-10820. M runs on FLD seeded from
I_M_DEPTH_NEXT (10693); P runs on FLD_P seeded from I_P_DEPTH_NEXT
(10772). KN_RAY_DUMP(1) fires at 10790 - exactly post-M/pre-P - so the
deck boundary already sits at the requested capture point.

| Quantity | Production producer | P consumer | In v2 deck? | Replay reconstruction | Verified |
|---|---|---|---|---|---|
| IM_NI(NCF) | FLD(NI,:) after M (10739) | P inner boundary in march | YES (phase-1, production value) | recomputed by replay M | runtime pending |
| P seed | I_P_DEPTH_NEXT store | SEED_FLD (10772) | YES | poked from deck | yes |
| M seed | I_M_DEPTH_NEXT store | SEED_FLD (10693) | YES | poked from deck | yes |
| March inputs (PK_*, GAM, WINX, JPR/IKC, CHI/ETA cubics, boundary consts, RES_FREQ) | refresh capture | both sweeps | YES | poked from deck | yes |
| FLD full M field | M sweep | only via IM_NI row NI | NO | recomputed | unverifiable from deck (row NI only) |
| RES_M/IXIT_M/NP_M | M sweep | census/KN_PA col 1; divergence gates P (10728) | NO (census sidecar has them) | recomputed | pending |
| M terminal omega | M sweep end | NOT consumed by P (P reseeds OMEGA_RAY=KN_RAY_OMEGA, 12308) | n/a | independent chain | yes (source) |
| RES_FREQ (M residuals) | M sweep | KN_HOT_CNT histogram only (10744-50); overwritten by P | n/a | write-only diagnostic | yes |
| KN_HOT_SPENT | refresh-global SAVE counter, shared across OpenMP rays | bounds hot-pass extension per sweep (12389-98) | NO (caps in deck, spend is not) | starts fresh in replay | DIVERGENCE CHANNEL |
| KN_HOT_CNT | per-ray RES_FREQ fold | summary print only | n/a | diagnostic | yes |
| WRC_* worst-cell | both sweeps | overwritten by P (10717) | n/a | diagnostic | yes |
| KN_PA_* | both sweeps | census only | n/a | diagnostic | yes |
| OMEGA seed KN_RAY_OMEGA | module | each sweep entry | YES (RCF(4), 8460) | poked | yes |

## Findings

1. The v2 deck ALREADY contains production's live IM_NI - the single
   most important M-derived P input. No new capture is needed for the
   boundary value itself; comparing replay-computed IM_NI against the
   deck value is possible on existing decks at runtime.
2. The deck does NOT contain: full M field, RES_M/IXIT_M/NP_M, M
   terminal omega, M per-pass schedule (M KNOMG rows exist in captures
   when armed: IDIR=2 via merge(2,1,INWARD)), or KN_HOT_SPENT.
3. NEW divergence channel: KN_HOT_SPENT (cap 480/refresh, 120/ray) is
   consumed by whichever rays stall first under OpenMP. Production's
   spend at this ray's start is scheduling-dependent; replay starts at
   zero, granting the full budget. Mode A replay can therefore extend
   M or P further than production did. Mode B is bounded by schedule
   length (extra pass request fails strictly). This channel affects
   STALLING rays - precisely 81/82's class - and cannot affect
   IM_NI when M exits before needing hot passes.
4. M and P have independent omega chains (reseeding at each entry).
5. Diverged M skips P and freezes both stores (10728-37): IXIT_M is a
   control-flow input to P's existence, absent from the deck (census
   sidecar carries it).

## Answer to the standing question

"Does the deck contain every M-derived quantity consumed by P?" - for
the FIELD handoff, yes: IM_NI is the only M quantity P's march reads,
and the deck records production's. The unverified part is whether
replay M reproduces it, and whether KN_HOT_SPENT state changes the M
(or P) pass count. An imposed-M experiment (poke deck IM_NI, still run
M for diagnostics) would isolate the downstream P contribution but is
explicitly not faithful replay.
