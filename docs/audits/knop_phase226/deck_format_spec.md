# KNRDUMP deck format - verified layout (20-Nov-2026)

Sequential unformatted Fortran records, one 28-record block per
refresh, appended (POSITION='APPEND'). Verified empirically on
KNRDUMP_L082_r0000.bin (43 blocks) and L081 (43 blocks).

Magic: rec0 int32. v1=20261122 (live d50dcacd exe), v2=20261123
(current source). v1 = v2 minus the two optional record groups.

rec | bytes (L082: NI=4,NDD=NXD=70,NCF=56377,NCPD=70) | content
 0  60   header ints: magic,LS,NI,NDD,NXD,NCF,KNC,NCPD,NADD,NCPS,
          IMETH,KMAXP,KCNT,FBITS,KNPD
 1  104  controls: KN_RAY_TOL,BULK_TOL,OMEGA,STALL,... (mixed i32/f64)
 2  12   KN_INNER_METH char
 3  560  R_EXT(NXD f64)
 4  560  LOG_R_EXT
 5  560  KN_R(NDD)
 6  8    KN_P(LS)
 7  32   Z(NI)
 8  32   R_RAY(NI)
 9  16   RAY_PNT(NI i32)
10  280  J_PNT(NDD i32)
11  280  COMP_DEPTH_IDX(NDD i32)
12  32   GAM(NI)
13  32   GAM_MAG(NI)
14  32   dGAMdR(NI)
15  126284480  KN_CHI_COEF(NXD,4,NCF) f64
16  126284480  KN_ETA_COEF
[17,18 v2 only, FBITS bit9: KN_CHI_EXT(NXD,NCF), KN_ETA_EXT]
17  451016   KN_FREQ(NCF) - nu in 1E15 Hz
18  451016   KN_DLN_P
19  451016   KN_DLN_N
20  451016   KN_IC
21  451016   KN_DBB
22  1804064  FLD (I_M channel, NI x NCF) - post-solve M field
23  31570560 I_M_DEPTH_NEXT(NCPD,NCF-1) f64 - donor coupling store
24  31570560 I_P_DEPTH_NEXT(NCPD,NCF-1) f64
[25,26 v2 only, FBITS bit14: JMP_SIB_M(NCF), JMP_SIB_P(NCF)]
25  1804064  FLD_P seed (NI x NCF)
26  24       trailer: RES_P(f64),NP_P(i32),IXIT_P(i32),OMEGA_RAY(f64)
27  1804064  FLD_P end (NI x NCF)

FBITS: bit0 FG_FREQ_MATCH 1 P2DAMP 2 ORD_SWEEP 3 AITKEN 4 AA
 5 probe 6 BULK_CONTR 7 ADAPT_OMEGA_V2 8 PLATEAU_DET
 9 HAS_EXT (v2) 10 SEED_COMMIT 11 RAY_SEEDED 12 Q_RESID 14 HAS_SIB (v2)
v1 deck FBITS=85 (bits 0,2,4,6): FREQ_MATCH+ORD_SWEEP+AA+BULK_CONTR.
NOTE: bit4 (AA) is set in the deck header - flag was ON in the
captured run's config namespace; check before interpreting
(AA was never enabled per VADAT history - this bit may reflect a
different flag ordering in the v1 writer).

v1/v2 delta is EXACTLY the optional groups; all shared records
identical in position and size. I_M/I_P_DEPTH_NEXT are REAL(LDP)
coupling arrays, NOT index maps.
