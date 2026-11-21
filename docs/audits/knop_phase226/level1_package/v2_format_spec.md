# KNRDUMP v2 format spec (magic 20261123)

Source of truth: `plane/fg_j_cmf_nonmono_v1.f`, `KN_RDUMP` writer,
lines 11205-11284 (build fd78d84a, source e28833f+dirty).

One logical block per ray per refresh, written in two calls on the
same append unit:

- PHASE=1 (pre-solve): records 1-25 (+3 optional)
- PHASE=2 (post-solve): records 26-27

## Record order (PHASE=1)

| # | contents | type | size |
|---|----------|------|------|
| 1 | 20261123, LS, NI, KN_ND, KN_ND_EXT, NCF, KN_NC, N_COMP, ND_ADD, NC_PNT_SRCE_L, IMETH_D, KN_RAY_MAXPASS, KN_CONT_NP, FBITS_D, KN_NP | i4 | 15 |
| 2 | KN_RAY_TOL, KN_RAY_BULK_TOL, KN_RAY_OMEGA, KN_RAY_STALL, NONMONO_INEXT_CAP, FG_COMMIT_NC, FG_COMMIT_JUMP, FG_DIV_HEAL_CAP, INEXT_RELAX, P2N, P2REL, HOTPASSES, HOT_TOTAL | r8 | 13 |
| 3 | KN_INNER_METH | char | str |
| 4 | R_EXT(1:KN_ND_EXT) | r8 | KN_ND_EXT |
| 5 | LOG_R_EXT(1:KN_ND_EXT) | r8 | KN_ND_EXT |
| 6 | KN_R(1:KN_ND) | r8 | KN_ND |
| 7 | KN_P(LS) | r8 | 1 |
| 8 | Z(1:NI,LS) | r8 | NI |
| 9 | R_RAY(1:NI,LS) | r8 | NI |
| 10 | RAY_PNT(1:NI,LS) | i4 | NI |
| 11 | J_PNT(1:KN_ND,LS) | i4 | KN_ND |
| 12 | COMP_DEPTH_IDX(1:KN_ND) | i4 | KN_ND |
| 13 | GAM(1:NI,LS) | r8 | NI |
| 14 | GAM_MAG(1:NI,LS) | r8 | NI |
| 15 | dGAMdR(1:NI,LS) | r8 | NI |
| 16 | KN_CHI_COEF(1:KN_ND_EXT,1:4,1:NCF) | r8 | KN_ND_EXT*4*NCF |
| 17 | KN_ETA_COEF(1:KN_ND_EXT,1:4,1:NCF) | r8 | KN_ND_EXT*4*NCF |
| 18* | KN_CHI_EXT(1:KN_ND_EXT,1:NCF) | r4 | fbit9 |
| 19* | KN_ETA_EXT(1:KN_ND_EXT,1:NCF) | r4 | fbit9 |
| 20 | KN_FREQ(1:NCF) | r8 | NCF |
| 21 | KN_DLN_P(1:NCF) | r8 | NCF |
| 22 | KN_DLN_N(1:NCF) | r8 | NCF |
| 23 | KN_IC(1:NCF) | r8 | NCF |
| 24 | KN_DBB(1:NCF) | r8 | NCF |
| 25 | FLD(1:NI,1:NCF) M-channel seed field | r8 | NI*NCF |
| 26 | I_M_DEPTH_NEXT(JJ,LS,KK), JJ=1..N_COMP, KK=1..NCF-1 | r8 | N_COMP*(NCF-1) |
| 27 | I_P_DEPTH_NEXT(JJ,LS,KK), JJ=1..N_COMP, KK=1..NCF-1 | r8 | N_COMP*(NCF-1) |
| 28* | JMP_SIB_M(1:NCF) | r8 | fbit14 |
| 29* | JMP_SIB_P(1:NCF) | r8 | fbit14 |
| 30 | FLD_P(1:NI,1:NCF) P-channel seed field | r8 | NI*NCF |

## Record order (PHASE=2)

| # | contents | type | size |
|---|----------|------|------|
| 26 | RES_P, NP_P, IXIT_P, OMEGA_RAY | r8,i4,i4,r8 | 4 items |
| 27 | FLD_P(1:NI,1:NCF) P-channel final field | r8 | NI*NCF |

## FBITS bit map (writer order, lines 11222-11239)

| bit | meaning |
|-----|---------|
| 0 | FG_FREQ_MATCH |
| 1 | FG_ORD_P2DAMP |
| 2 | FG_ORD_SWEEP |
| 3 | FG_ORD_AITKEN |
| 4 | FG_ORD_AA |
| 5 | L_NOPROBE |
| 6 | FG_BULK_CONTR |
| 7 | FG_ADAPT_OMEGA_V2 |
| 8 | FG_PLATEAU_DET |
| 9 | KN_CHI_EXT/KN_ETA_EXT present |
| 10 | FG_SEED_COMMIT |
| 11 | ray seeded |
| 12 | KN_Q_RESID |
| 13 | KN_CERT_CONV |
| 14 | JMP_SIB_M/JMP_SIB_P present |

## v1 deltas (magic 20261122)

- No records 18-19 (CHI_EXT/ETA_EXT); v1 reader reconstructs node
  values via EXP(KN_CHI_COEF(KK,4,JJ)) / eta analog - exact at
  tabulated nodes under the monotone cubic, interpolant values
  between nodes are NOT recoverable.
- No records 28-29 (JMP_SIB_*); reader zero-fills; sibling maps feed
  diagnostics only.
- No KN_DBB record in v1.
- fbits must be masked to bits 0-8 on v1 read.
- Deck filename: KNRDUMP_L%3.3d_r%4.4d.bin, refresh index =
  MOD(LEDG_N_RF,10000); appended across refreshes.
