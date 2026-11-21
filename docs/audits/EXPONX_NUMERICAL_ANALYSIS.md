# EXPONX numerical analysis (subs/exponx.f)

## The function and its branches

```
|X|<1e-3 : 1 - X(0.5 - X/6(1 - X/4))      (Taylor, 4 terms)
X<40     : (1 - EXP(-X))/X
else     : 1/X
```

The branch condition `X < 40` accepts ALL negative X — there is no
negative-X guard. For X<0, EXP(-X)=EXP(|X|) grows without bound.

## Verified behavior (f64, independent evaluation, audit 6 harness)

| X | EXPONX | d_EXPONX_dX | note |
|---|---|---|---|
| 1000 | 1.0e-3 | ~-1e-6 | 1/X branch, correct |
| 10 | 9.9995e-2 | | exact |
| 1 | 6.321e-1 | | exact |
| 0.1 | 9.516e-1 | | exact |
| 1e-6 | 1.0 | -0.5 | Taylor branch, correct |
| 0 | 1.0 | -0.5 | Taylor limit |
| -0.1 | 1.0517 | | correct (F>1, amplification) |
| -1 | 1.7183 | | correct |
| -10 | 2.203e3 | | correct |
| -100 | 2.7e41 | | correct |
| -700 | 1.4e301 | | correct but near f64 max |
| -710 | +Inf | Inf | EXP(+710) overflows |
| -1000 | +Inf | Inf | overflow |

Overflow threshold: X < ~-709.78 (f64 EXP limit). Between -700 and -709
the result is finite but ~1e300 — Inf-equivalent in any physical sum.

## Mathematical notes

- F(-a) = (exp(a)-1)/a grows exponentially; it is the correct
  unsaturated-amplification escape factor, NOT replaceable by 1/X or
  |X| tricks. Any "fix" changing its value for X>-1 would silently
  alter real weak-inversion physics.
- Small-X Taylor truncation: O(X^4/120) — error <1e-14 at |X|=1e-3. Sound.
- d_EXPONX_dX = [(X+1)exp(-X)-1]/X^2 implemented as (Y-(1-Y)/X)/X,
  Y=exp(-X). For X>=40 uses Y=0 giving -(1)/X^2 + ... — consistent
  with F's own 1/X branch (F' = -1/X^2 there; check: (0-(1-0)/X)/X=-1/X^2,
  correct). Taylor branch: -0.5+X/3-X²/8+X³/30 — correct expansion.
  Same overflow point at X<-709 (Y=Inf -> Inf/X² → -Inf... actually
  (Inf - Inf/X)/X → Inf/-700 = -Inf, sign correct).

## The latent window

EXPONX's only structural vulnerability is negative X with |X|>~50
(huge-but-finite) and |X|>~710 (Inf). Whether this is reached depends
entirely on the upstream guards — see CALLSITE_TABLE. The function
itself is mathematically faithful everywhere.
