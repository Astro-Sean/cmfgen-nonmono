# EXPONX recommendations (no implementation authorized)

## 1. Diagnostic before any code change (P3, zero-risk)

- Add a one-line min(GAM) accumulator in SOBJBAR_SIM logged per
  (group,depth,iteration) — establishes whether the (-0.5,0) window is
  ever populated. This is the only datum currently missing.
- Extend NEG_OPAC to also record group sums in (-0.5,0) at LST_ITERATION.

## 2. If a safeguard is justified (P2 — requires evidence from step 1)

Options in order of physical defensibility:

a. Keep EXPONX untouched; add a finite-GAM guard in the callers
   (SOBJBAR_SIM/SOBEW): if GAM < GAM_FLOOR (e.g. -50), treat the cell
   like a flagged cell or cap at the unsaturated-gain ceiling.
   Physically motivated: F(-a)=(e^a-1)/a is correct but a real line's
   gain saturates; the mathematically honest fix is a saturation cap,
   not clipping to 1/X or |GAM| (both wrong — would turn amplification
   into absorption).

b. Extend the crush from summed-tau<-0.5 to summed-tau<0 (makes the
   window unreachable). Cost: excises all inversions entirely; safe
   direction but removes the (real, weak) gain the unsaturated
   expression currently models for mild cells.

c. d_EXPONX_dX must receive whatever treatment EXPONX gets — they are
   used in the same VB terms; fixing one without the other is
   unacceptable.

Do NOT: clip EXPONX itself (silently changes correct math for X>-1),
replace F(-a) with 1/X (mathematically wrong sign/behavior), or treat
convergence as proof the window is empty.

## 3. Required tests before any implementation

1. X=0 and ±epsilon limits (F→1, F'→-0.5).
2. X=+1,+10,+40,+100 against analytic reference.
3. X=-1,-10,-40 finite and matching (e^a-1)/a.
4. X=-700 finite; X=-710/-1000 -> documented Inf behavior.
5. d_EXPONX_dX vs finite differences for |X|<50.
6. Isolated inverted-line EW calculation sign check (negative EW).
7. Stock monotonic regression: bit-identical BETA/EW/rates.
8. Non-monotonic regression on the oxygen model.
9. Guard activation test: synthetic group with tau in (-0.5,0) and
   T1~1e-4 — verifies the window closes.
10. JFNK replay path test confirming the duplicate crush stays in sync
    (the two copies at ~4137 and ~7300 must not diverge).

## 4. Physics follow-up (separate, larger task)

The crushed cells are genuine recombination-cascade inversions
(HI/He2/O IV-VI high-n). Whether a saturation-capable treatment
(Gray 2012-style) matters for the emergent spectrum is open but bounded:
affected lines are far-IR/Rydberg at trace depths. Not blocking.
