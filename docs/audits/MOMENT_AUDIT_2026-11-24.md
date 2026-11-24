# Moment-machinery audit: mom_j_cmf_nonmono_v1.f

User-directed 24-Nov-2026. Question: is the moment-equation
procedure correct and robust. Method: line-level static audit of
each component, cross-checked against the python-driven unit and
regression batteries. Verdicts: VERIFIED / VERIFIED-WITH-NOTES /
DEFECT / OPEN.

## 1. Row assembly (TA/TB/TC/XM), interior rows ~:2604-2622

```
TA(I) = -HL(I-1)            TC(I) = -HU(I)
TB(I) = DTAUONQ*MAX(MOM_ABS_FLR,1-COH) + PSI + HL(I) + HU(I-1)
```
EPS sign-split (per midpoint, sign of EPS_PREV):
- EPS(I-1) >= 0: TA -= 2*EPS(I-1)   (upwind donor to subdiagonal)
- EPS(I-1) <  0: TB -= 2*EPS(I-1)   (adds +2|EPS| to diagonal)
- EPS(I)   >= 0: TB += 2*EPS(I)     (diagonal)
- EPS(I)   <  0: TC += 2*EPS(I)     (superdiagonal -= 2|EPS|... 
  TC += 2*EPS with EPS<0 lowers TC)

VERIFIED-WITH-NOTES. The split is the standard donor-cell
upwind form: positive EPS (outward scatter) couples to the
subdiagonal/diagonal, negative EPS (inward, compressing) loads
the diagonal/superdiagonal. Conservation check: the sum of the
EPS contributions to (TA+TB+TC) is -2E(I-1)+2E(I) regardless of
signs, matching the divergence form of the scatter term - the
sign split only redistributes which entry carries each piece, so
the row sum (and hence the continuum limit) is sign-exact.
Note 1: TC can go positive when 2|EPS| > HU (the census counts
these as MOM_NONMM_CNT); the depth-upwind comment at :2645
describes the intended form. Positive TC breaks the M-matrix
sign pattern locally - this is physical (compressing scatter is
genuinely anti-diffusive) and is exactly the population LAG_BND
then descales.
Note 2: MOM_ABS_FLR floors the absorption coefficient - a small
positive floor kept for conditioning; documented in NONMONO_KEYS.

## 2. Frequency coupling assembly ~:2965-2990

Expanding rows (Q>=0): RHS gets PSIPREV*J(k-1) lagged.
Compressing rows: RHS gets PSI*J_next(k+NMJ) with PSI = diagonal
|Q_k| weight, matching FG's compressing source form. Equal
explicit/implicit weights keep the lagged map at unit spectral
radius rather than amplifying.

VERIFIED. The compressing coefficient uses the full diagonal
weight PSI (not PSIPREV), which is the correct Baron-style
splitting for the compressing direction: the implicit piece of
the coupling sits in TB, the explicit piece in the RHS, with
matching magnitude.

## 3. MOM_LAG_BND scaling ~:2985-3005

```
T4 = coupling + |VB| + |VC| + 2(|EPS(I-1)|+|EPS(I)|)
T5 = TB - |TA| - |TC|            (diagonal excess)
LAG_SCL = min(1, max(T5,0)/T4);  XM -= (1-LAG_SCL)*XLAG
```

VERIFIED. The scaling enforces the M-matrix property on the
lagged operator splitting: post-scaling lagged weight <= excess,
so every row's effective coupling ratio a < 1 strictly. Rows
already non-M-matrix (T5 <= 0) get LAG_SCL = 0 - full descaling
of lagged content, the correct conservative limit. This is the
bound that makes matrix-side amplification through the NU chain
solve impossible, confirmed empirically (NMARG=0 arm 13653300).
Note: the census (arm 13653300) shows 864K-1031K cells/pass at
the bound - a large fraction of the compressing grid lives at or
past the M-matrix edge, i.e. the diffusion-style closure is
marginal there by construction, not by defect.

## 4. MOM_XM_POS clipping ~:3020-3032

Negative XM rows are zeroed (debit limited to the balance), with
census MOM_XMWORST/MOM_XMNEG.

VERIFIED-WITH-NOTES. Mathematically local and documented. The
structural caveat is physical, not implementational: negative XM
arises only where lagged H-divergence/EPS/upwind content
outweighs a near-zero local source - compressing absorption-
dominated rows. Arm 13653300 censused thousands of such rows per
refresh with worst clipped magnitudes 1.9e5-4.5e5. The clip is
the correct guard for the solve, but its firing rate is the
census signature that the operator split itself is marginal in
exactly those regions (advection-dominated transport does not
decompose into positive sources). Not a defect; a diagnostic.

## 5. H reconstruction and lagged store ~:4160-4274

H is algebraically reconstructed per frequency from the solved
same-frequency J row plus upwind-H/EPS terms; causality-clipped
|H| <= J pointwise; non-finite entries zeroed; the lagged store
H_NEXT_FRQ is updated with the same relaxed/damped ownership
write as J_NEXT. H drift is folded into the max-norm
MOM_GS_MAXCHG only (sign-crossing cells would inflate fractional
metrics - correct justification in the comments).

VERIFIED-WITH-NOTES. The algebra is consistent by construction:
H is not an independent solution but a derived quantity of the J
row, so no separate tridiagonal inconsistency can arise. The
structural consequence for the iteration (phase-3 finding):
because H scales with the committed J field, the lagged H
content in X0 scales with whatever field was last committed -
including inflated commits. The feedback loop of the nusolv arms
runs through this store, not through the matrix inverse. The
store update logic itself is correct; the loop-gain defect is in
the lagged splitting it serves.

## 6. Boundary rows

Inner boundary ~:2660-2700: four INNER_BND_METH forms, all divide
by DTAU(ND-1) with a positivity guard rejecting the frequency on
nonpositive increments. TC(ND)=VB=VC=PSIPREV=0 correctly close
the system.
Outer boundary ~:3040+: sign(A_bnd) keys the upwind choice (not
Q_NM(1)) after the 24-Sep-2026 audit showed misclassification at
sigma(1) ~ -1; boundary-Eddington floor for negative-flux
frequencies.

VERIFIED. Both boundaries carry audit-hardened guards with the
physical reasoning documented inline.

## 7. Capture consistency for MOM_NUSOLV

The capture records TA/TB/TC and the post-scaling RHS after every
bound (LAG_SCL, EPS_SCL, XM_POS) is folded in, so the chain solve
solves exactly the equations DGTSV solved - fixed-point
consistent. Verified numerically: test_mom_nu 53/53 including
exact dense-solve agreement; the iterated probe reproduces
theoretical growth rates to 4+ digits (unguarded 1000.00/iter =
c/(1-a) exactly).

VERIFIED.

## 8. Write-back and damping

Damped non-negative update XOLD + RLX*(XNEW-XOLD) with rejection
of negative/nonfinite results, per-cell damping MOM_DAMP_NM during
refresh, ownership-aware store writes (MOM_RF_OWNW) so compressing
midpoints are written only by the sweep direction that owns them.

VERIFIED. Ownership logic is subtle but documented and unit-
covered through test_mom_solve refresh arms.

## 9. Python test verification (this audit)

Re-run of the full battery against the current tree:

| battery | result | coverage |
|---|---|---|
| test_mom_nu | 53/53 | chain solve exactness, freeze band, write cap, iterated dynamics |
| test_mom_solve | 143/143 | fixed-point regression, refresh path, nurf/nuw arms |
| test_mom_gs | 34/34 | GS ordering, drift metrics |
| test_mom_aa_state | 62/62 | AA state management |
| test_mom_sweep_table | 8/8 | sweep table logic |
| test_vadat_keys | 89/89 | key plumbing incl MOM_NUSOLV, MOM_H_IMPL |

## 10. Coverage gaps (open items)

1. No coefficient-level numeric probe isolates the EPS sign-split
   through the full routine (the conservation argument in section
   1 is analytic). A manufactured grid with controlled EPS signs
   could pin this - moderate value, not yet built.
2. H reconstruction is verified structurally but not against an
   independent manufactured solution with nonzero H.
3. Boundary-row sign logic is audit-hardened but has no dedicated
   unit driver coverage at sigma(1) near -1.

## Verdict

The moment machinery is correct as implemented - no assembly
defect found. The instability the programme is fighting is not
an implementation bug: it is the lagged operator splitting's
structural marginality in compressing, advection-dominated
regions, evidenced by the XM_POS/lagbnd censuses. The procedure
is robust in the guard sense (nothing propagates unbounded
without census), but the iteration it serves cannot contract
there. Dated: 24-Nov-2026. Verdict: VERIFIED per component;
structural finding stands as the phase-3 input.
