# KNRDUMP deck format spec (v1=20261122, v2=20261123)

Writer: KN_RAY_DUMP in fg_j_cmf_nonmono_v1.f ~11200 (phase-1 header
per refresh, phase-2 trailer per solve: RES_P,NP_P,IXIT_P,OMEGA_RAY,
FLD_P). Reader: unit_tests/f90/drv_knop_replay.f90.

The magic numbers are date-styled literals the author picks at
writer-edit time, not runtime timestamps. 20261122/20261123 label
two consecutive writer revisions; the v1 tag predates the local
checkout date only because the literal was chosen ahead of it.

Record order (Fortran unformatted sequential):
 1  header ints: magic,LS,NI,KN_ND,KN_ND_EXT,NCF,KN_NC,N_COMP,
    ND_ADD,NC_PNT_SRCE_L,IMETH,KN_RAY_MAXPASS,KN_CONT_NP,FBITS,KN_NP
 2  controller dbls: TOL,BTOL,OMEGA,STALL,INEXT_CAP,FG_COMMIT_NC,
    FG_COMMIT_JUMP,FG_DIV_HEAL_CAP,INEXT_RELAX,P2N,P2REL,HOTP,HOTT
 3  KN_INNER_METH (char)
 4-15 geometry/profile: R_EXT,LOG_R_EXT,KN_R,KN_P,Z,R_RAY,RAY_PNT,
    J_PNT,COMP_DEPTH_IDX,GAM,GAM_MAG,dGAMdR
 16 KN_CHI_COEF(NXD,4,NCF)   17 KN_ETA_COEF(NXD,4,NCF)
 18 v2-only IF FBITS bit9: KN_CHI_EXT(NXD,NCF)
 19 v2-only IF FBITS bit9: KN_ETA_EXT(NXD,NCF)
 20 KN_FREQ  21 KN_DLN_P  22 KN_DLN_N  23 KN_IC  24 KN_DBB
 25 FLDM seed field (NI,NCF)
 26 I_M_DEPTH_NEXT(NCPD,NCF-1)
 27 I_P_DEPTH_NEXT(NCPD,NCF-1)
 28 v2-only IF FBITS bit14: JMP_SIB_M(NCF)
 29 v2-only IF FBITS bit14: JMP_SIB_P(NCF)
 30 FLDP_SEED (NI,NCF)
 31+ per-solve trailers: RES_P,NP_P,IXIT_P,OMEGA_RAY + FLD_P

v1 decks lack records 18/19 and 28/29 and never set FBITS bits 9+.
BUILD_FREQ reads KN_CHI_EXT unguarded for on-node points (NGRD>0),
which is why a bare v1 read segfaults. Reader fix (this phase):
explicit version table, v1 reconstruction EXT=EXP(COEF(K,4,:))
(exact node values of the monotone cubic built on LOG(EXT)),
JMP_SIB zero-filled (diagnostic prints only). Limitations: v1 bits
10-13 (seeded/cert flags) were not written by the v1 writer - the
reader clears them, so replay treats the ray as unseeded; the
first-commit jump veto can differ from production for a ray whose
seed would have been rejected there.
