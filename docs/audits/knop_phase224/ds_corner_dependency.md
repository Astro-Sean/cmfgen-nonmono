# DS-corner dependency (ordered P march, fg_j_cmf_nonmono_v1.f)

P-march loop ~12410: DO J=NI-1,1,-1; preview fill
KN_ORD_SPW_FILL(FLD,J-1) writes SPW slot for the not-yet-marched
row J-1 from its PRE-march FLD state (comment ~12171: "the only
not-fresh source this step reads", "lagged C-corner source of
Knop et al. 2009, Paper II Eqs. 31/32" ~11783).

Per cell (J,L):
  SPJ  = KN_ORD_SPV(FLD,J,L)              fresh source, current pass
  SPW(J)=SPJ                              ring fill
  DSU  = KN_ORD_DS(J+1,L)                 upwind slope, fresh rows
  DSD  = KN_ORD_DS(J,L)  -> SEC(J-1),SEC(J)
         SEC(J-1) = (SPW(J)-SPW(J-1))/DTAU(J-1)
         SPW(J-1) = PREVIEW of row J-1 = S' evaluated on
                   FLD(J-1) BEFORE this pass's march reached it
  VNEW = FLD(J+1)*PK_A0 + SPW(J+1)*PK_A1 + SPJ*PK_A2
         - DSU*PK_A3 - DSD*PK_A4

For J=2: DSD=DS(2) couples fresh S'(2) to PREVIOUS-ITERATE S'(1).
Row 1 is the vacuum-boundary row; its update happens after J=2 is
marched (J descends to 1), so DS(2) carries a one-pass Jacobi lag.
The lag vanishes at the fixed point (identical successive iterates)
- it changes the iteration sequence, not the map's fixed points.

Measured (replay, production trajectory, passes 13-16):
bad cells by row: j1=2030 j2=2764 j3=798 j4=0 - the boundary pair
plus a row-3 tail. j1/j2 committed deltas move in lockstep, same
sign, j2 attenuated (e.g. l=1434: j1 +1.46e-27,+9.93e-28,-7.36e-26;
j2 +6.74e-28,+4.43e-28,-1.65e-25). Rows 3-4 frozen at the argmax
columns. Consistent with the documented coupling; isolating the
lag's causal share needs per-cell source-term instrumentation.
