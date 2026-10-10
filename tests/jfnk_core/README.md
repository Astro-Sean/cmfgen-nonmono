# tests/jfnk_core

Standalone unit test of the Newton core -- `newton_hook_v1` and
`gmres_m_v1` from `new_main/mod_subs/jfnk_solver_v1.f90`, compiled
unchanged.  `jfnk_int_stub.f90` replaces `jfnk_cmfgen_int_mod` (the
production module carries the whole CMFGEN model state).

Problem: N=60 stiff synthetic residual

    F_i(x) = a_i*(exp(x_i) - exp(xs_i)) + b*(x_i - x_{i-1})
             + b*(x_i - x_{i+1})

with a_i log-spaced 1e-3..1e3, b = 0.3, known root xs_i =
1.5*sin(1.7*i).  The residual callback adds the pseudo-transient
term (x - x0)*s/dtau exactly as JFNK_RESIDUAL_V1 does, so the tests
exercise the 06-Oct-2026 SER anchor fix: with the old
baseline-anchored x*s/dtau form, the Newton norms measured |G| and a
returned "converged" iterate could still have F /= 0.

Tests:

- T1: PTC off, start xs+0.5 -> info=0, max|x-xs| < 1e-6.
- T2: PTC on, baseline-shifted unknown (production form), start at
  u=0; asserts |F_true| < tol at the returned iterate.
- T3: PTC on, far start xs+3.0 (exp overshoot exercises trust
  region/line search/Cauchy); prints jfnk_n_ls_used,
  jfnk_n_tr_rejected, jfnk_n_cauchy_used.

Run:

    module load intel/2024.0 mkl
    cd tests/jfnk_core
    make test

Exit code nonzero if any test FAILs.
