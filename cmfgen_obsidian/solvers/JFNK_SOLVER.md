# JFNK solver: problem and the changes in this tree

This note covers the Jacobian-free Newton-Krylov population solve
used when VADAT `DO_JFNK=T` (`USE_JFNK` is a deprecated alias).  The
default CMFGEN path (`DO_JFNK` absent or `F`) is unchanged.  JFNK is
admitted under `DO_NONMONO_KNOP` or `DO_NONMONO_SOB`; under SOB only
the frozen-J residual tier runs, because the coupled eval sweep does
not recompute the Sobolev line rates.

The method follows Knoll & Keyes (2004, J. Comput. Phys. 193, 357)
and Arramy, de la Cruz Rodriguez & Leenaarts (2024, A&A 690, A12;
arXiv:2406.17234).  Arramy's tests are static plane-parallel solar
atmospheres with six-level atoms.  The speed-up they report does not
transfer to this supernova model.  The convergence requirements do.

## What is being solved

On a full iteration the residual F(n) is the statistical-equilibrium
equation at the trial populations.  With `JFNK_COUPLED=T` the mean
intensity is recomputed by the Knop formal solution at those
populations.  The lagged frequency-coupling store (I_NEXT) is held
fixed.  That is one inner step of the nested iteration: the radiation
field responds to the populations, and the redward intensity from the
previous great iteration stays a parameter.  The finite-difference
Jacobian-vector product uses the same residual, so it includes dJ/dn.
The depth-block BA matrix is only the preconditioner.

Admission also requires `USE_LAM_ES=T`.  With `USE_LAM_ES=F`, the
MOM formal-solver path updates `J_NEXT_FRQ` and `H_NEXT_FRQ` in place
and has no evaluation-time freeze.  JFNK stays deferred on that path
to prevent trial residuals from mutating the lagged radiation state.

## Where the solve fails

The Dessart m15 2.82 d run is the example.  Lambda iterations bring
the largest population correction down to about 40 percent.  The
first full iteration then asks for a raw correction of order 1e7
percent.  Three separate defects produce that step.

1. Newton is started outside its basin.  The linear model of F is
   valid only near the root.  Arramy et al. (2024) report the same
   stall when H I is started from LTE, and they name the missing
   pieces: a line search, a continuation method, and a few
   operator-split iterations before the first Newton step.  Switching
   on the fully coupled operator at a 50 percent correction is still
   far from that neighbourhood.  A hard factor-of-ten clip limits the
   step but does not follow the residual.

2. The unknown can cross a kink.  Populations used
   n = n0 * (1 + x) with a floor at 1e-100.  The floor makes F
   non-differentiable, and a forward-difference Jacobian that crosses
   it is not the derivative GMRES assumes.  Arramy et al. (2024,
   their Fig. 2) show that a plain forward difference is a poor
   Jacobian on this problem.  A complex-step difference is accurate
   and is not practical here: the formal sweep branches and clips.
   Logarithmic unknowns are the change that removes the floor.

3. The inner linear solve is asked to do too much with too local a
   preconditioner.  Arramy et al. found their block-Jacobi
   preconditioner, which is the same object as this BA matrix,
   inadequate once the problem has many frequencies or strong
   scattering.  Their useful GMRES relative tolerance sat between
   about 3e-4 and 3e-2.  Tighter than that spends formal solutions
   without improving the Newton step.  Janett et al. (2024, A&A 682,
   A48) show a simplified physical operator used as a preconditioner
   cutting Krylov iterations by a large factor, and that LTE itself
   is useless as a preconditioner.

The lagged I_NEXT store should stay lagged.  Differentiating through
it would make every matrix-vector product a second frequency
coupling.  That is not the nested iteration.  Once the inner solve
is reliable at fixed store, the outer update of the store is a
fixed-point map.  Anderson acceleration is the method for that map.
It is not a substitute for the inner Newton solve, and it needs the
map to be contractive first.

## What this tree now does

All three are on by default for a JFNK solve.  Each has a VADAT key.

`JFNK_LOGX` (default T).  The unknown for every level and for the
electron density is x = ln(n/n0), so n = n0 * exp(x).  Temperature
already used this map.  At x = 0 the derivative matches the old
fractional unknown, so the BA preconditioner is still the right
local Jacobian at the start of the step.  Away from x = 0 the
finite-difference product carries the extra factor exp(x), and the
preconditioner apply multiplies by exp(-x) so the two stay
consistent.  The committed step is still limited to
`MAX_JFNK` in population ratio, which is |x| <= ln(MAX_JFNK).
Inside a residual evaluation x is only clipped at +/-40, to keep
exp from overflowing; that clip is counted and the difference step
is retried, as before.

`JFNK_PTC` (default T).  Pseudo-transient continuation
(Knoll & Keyes 2004, JCP 193, 357 Sec 4.2; switched evolution
relaxation, Kelley & Keyes 1998, SIAM J. Numer. Anal. 35, 508).
The residual Newton sees is

    G(x) = F(x) + (x - x_k) * s / dtau

where x_k is the current Newton iterate (the anchor; re-anchored at
each accepted step) and s_i = max(|F_i| at the current baseline, a
small floor) so the added term has the same units as F.  The term
vanishes at x_k, so the stored residual is F; every residual the
Newton loop reads has the term stripped back out, and all norms,
ared/pred, the EW forcing term and convergence tests act on |F|
(06-Oct-2026 fix: previously the term was anchored at the baseline
x = 0 and the norms acted on |G|, so "converged" could mean
F(x) = -x*s/dtau /= 0).  dtau starts at 1.  After a
Newton iteration that lowers |F|, dtau is multiplied by
|F_old|/|F_new|, and it is capped at 1e8.  If |F| does not fall,
dtau is halved, down to 0.01.  Large dtau drops the extra term and
the iteration is Newton's method.  The GMRES operator keeps the
diagonal: it solves (J + D/dtau) s = -F(x_k).  The preconditioner rescales the
BA solve by `d'/(d'+sigma)`, where `d'` is the signed diagonal of
the applied preconditioner, including its ridge or identity-row repair,
and `sigma=s/dtau`.  This is scalar-exact for
either sign of the BA diagonal, tends to the BA solve as dtau grows,
and tends to the PTC diagonal solve as dtau shrinks. If `d'` is
zero, the input residual is divided by `sigma` directly; if
`d'+sigma` is near zero, the BA result is retained rather than
dividing by a singular value.

GMRES forcing floor.  The Eisenstat-Walker eta is not allowed below
3e-4.  That is the bottom of the range Arramy et al. found useful.
The solver still loosens eta when the residual is falling fast, so
a distant iterate does not pay for a tight inner solve.

`JFNK_EPS` is a lower bound on the perturbation norm
`h*||v||inf` in coupled mode. Clip retries preserve that floor, using
the actual perturbation norm rather than comparing the scalar h,
which also depends on `||v||inf`. Frozen residuals use deterministic
stored-array arithmetic and may continue shrinking on clip retries.

## What was left out on purpose

A full Chen & Shen (2006) update of the preconditioner from the
Krylov products.  The diagonal continuation blend above is the part
required to match G(x).  Rebuilding a nonlocal preconditioner from
the matvecs is the next preconditioner change, not this one.

Raising the GMRES restart dimension (`JFNK_M`, default 10).  Arramy's
many-frequency cases wanted a larger space.  Each extra vector is
another formal solution.  Leave the dimension alone until a run
shows GMRES stagnating inside the forcing floor.

Differentiating through I_NEXT, and Anderson acceleration of the
outer store update.  Both belong to the outer iteration.  They are
not part of making one Newton step a descent direction.

## OUTGEN

Once per JFNK solve the formulation block names the unknown, whether
pseudo-time is on, the starting dtau, and the forcing floor.  Once
per accepted Newton iteration a single line gives the old and new
dtau and |F|.  Per-frequency sweep text stays in `JFNK_EVAL.log`.


## Testing

`tests/jfnk_core/` is a standalone unit test of `newton_hook_v1` on a
stiff synthetic system; see its README.  Run:
`cd tests/jfnk_core && module load intel/2024.0 mkl && make test`.
