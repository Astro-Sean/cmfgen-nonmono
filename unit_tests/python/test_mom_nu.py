#!/usr/bin/env python3
"""Verification of MOM_NU_CHAIN_SOLVE (MOM_NUSOLV path).

drv_mom_nu.f seeds the NU_* capture arrays and J_NEXT_FRQ with
manufactured row equations, calls the production per-depth
frequency-chain solve, and dumps the output field.  This test
solves the same directed-graph systems independently (dense NCF
x NCF solve of x = r + a*x_upwind, per depth, with the same
Gauss-Seidel-in-depth ordering) and compares, plus checks the
ring/freeze/rejection census counters.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "new_main/mod_subs/iter_ledger_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "subs/thomas.f",
    "plane/mom_j_cmf_v11.f",
]

ND = 4
NCF = 10
AV0 = 0.5
RV0 = 0.5


def build():
    return compile_and_link_prod(
        "mom_nu",
        [TREE / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_nu.f"],
    )


def load(tag):
    d = np.loadtxt(WORKDIR / f"{tag}_xv.txt")
    i, k = d[:, 0].astype(int), d[:, 1].astype(int)
    jold = np.full((ND, NCF), np.nan)
    xnew = np.full((ND, NCF), np.nan)
    jold[i - 1, k - 1] = d[:, 2]
    xnew[i - 1, k - 1] = d[:, 3]
    return jold, xnew


def kv(out, key):
    for line in out.splitlines():
        if line.startswith(key + "="):
            return float(line.split("=", 1)[1])
    raise KeyError(key)


def graph_solve(r, a, t):
    """Independent dense solve of x(k) = r(k) + a(k)*x(t(k))."""
    n = len(r)
    mat = np.eye(n)
    for k in range(n):
        if t[k] > 0:
            mat[k, t[k] - 1] = -a[k]
    return np.linalg.solve(mat, r)


def expect_simple(ncup, ndir, r=RV0, a=AV0):
    """exp (ndir=-1) or cmp (ndir=+1) all-rows chain, boundary at the
    unlinked end."""
    r_v = np.full(NCF, r)
    a_v = np.zeros(NCF)
    t_v = np.zeros(NCF, int)
    if ndir == -1:
        a_v[1:] = a
        t_v[1:] = np.arange(1, NCF)
    else:
        a_v[:-1] = a
        t_v[:-1] = np.arange(2, NCF + 1)
    return np.tile(graph_solve(r_v, a_v, t_v), (ND - 1, 1))


def run_case(rep, exe, mode, tag):
    r = run(exe, mode, tag, cwd=WORKDIR)
    rep.check(f"{mode}: driver completed", "RESULT_OK" in r.stdout)
    _, xv = load(tag)
    return r.stdout, xv


def main():
    rep = Report("test_mom_nu")
    exe = build()

    # --- exp: pure expanding chain ------------------------------------
    # x(k) = r(1-a^k)/(1-a) = 1 - a^k for r=a=0.5.
    out, xv = run_case(rep, exe, "exp", "nu_exp")
    exp = expect_simple(NCF, -1)
    rep.check(
        "exp: solved field is the exact geometric chain",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )
    rep.check(
        "exp: frozen depth row untouched",
        np.allclose(xv[ND - 1], 1.0),
        "depth ND must not be written",
    )
    rep.check("exp: no rings", kv(out, "NU_NRING") == 0)

    # --- cmp: pure compressing chain ----------------------------------
    out, xv = run_case(rep, exe, "cmp", "nu_cmp")
    cmp_ = expect_simple(NCF, +1)
    rep.check(
        "cmp: solved field is the exact back-propagated chain",
        np.allclose(xv[: ND - 1], cmp_, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-cmp_).max():.2e}",
    )
    rep.check("cmp: no rings", kv(out, "NU_NRING") == 0)

    # --- ring2: expanding/compressing 2-cycle -------------------------
    # {2,3}: x2 = r2 + a2*x3, x3 = r3 + a3*x2; all other nodes
    # uncoupled (x=r).
    out, xv = run_case(rep, exe, "ring2", "nu_ring2")
    a2, a3 = 0.5, 0.3
    x2 = (RV0 + a2 * RV0) / (1.0 - a2 * a3)
    x3 = (RV0 + a3 * RV0) / (1.0 - a2 * a3)
    exp = np.full((ND - 1, NCF), RV0)
    exp[:, 1] = x2
    exp[:, 2] = x3
    rep.check(
        "ring2: 2-cycle solved exactly",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )
    rep.check(
        "ring2: one ring per depth",
        kv(out, "NU_NRING") == ND - 1,
        f"NRING={kv(out, 'NU_NRING')}",
    )
    rep.check("ring2: none singular", kv(out, "NU_NSING") == 0)

    # --- ring4: stride-3 edge closes a 4-member ring -------------------
    # 1->4 (cmp a=.5), 4->3, 3->2, 2->1 (exp a=.4); rest uncoupled.
    out, xv = run_case(rep, exe, "ring4", "nu_ring4")
    t_v = np.array([4, 1, 2, 3])  # 1->4? members: 1 reads 4; 4 reads 3...
    # walk order from 1: 1 -> t(1)=4 -> t(4)=3 -> t(3)=2 -> t(2)=1 ring
    r_v = np.full(4, RV0)
    a_v = np.array([0.5, 0.4, 0.4, 0.4])
    t_loc = np.array([4, 1, 2, 3])
    sol4 = graph_solve(r_v, a_v, t_loc)
    exp = np.full((ND - 1, NCF), RV0)
    exp[:, 0] = sol4[0]
    exp[:, 1] = sol4[1]
    exp[:, 2] = sol4[2]
    exp[:, 3] = sol4[3]
    rep.check(
        "ring4: 4-member ring solved exactly",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )
    rep.check(
        "ring4: one ring per depth",
        kv(out, "NU_NRING") == ND - 1,
        f"NRING={kv(out, 'NU_NRING')}",
    )

    # --- mix: chains feeding two 2-cycles ------------------------------
    # 1 bnd; 2 exp->1; 3 cmp->4; 4 exp->3 (ring); 5 cmp->6;
    # 6 cmp->7; 7 exp->6 (ring); 8 bnd.
    out, xv = run_case(rep, exe, "mix", "nu_mix")
    r_v = np.full(NCF, RV0)
    a_v = np.array([0.0, 0.4, 0.5, 0.3, 0.4, 0.5, 0.3, 0.0, 0.0, 0.0])
    t_v = np.array([0, 1, 4, 3, 6, 7, 6, 0, 0, 0])
    sol = graph_solve(r_v, a_v, t_v)
    exp = np.tile(sol, (ND - 1, 1))
    rep.check(
        "mix: chains + two rings solved exactly",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )
    rep.check(
        "mix: two rings per depth",
        kv(out, "NU_NRING") == 2 * (ND - 1),
        f"NRING={kv(out, 'NU_NRING')}",
    )

    # --- depth: lagged TA/TC terms ------------------------------------
    # I=2's rhs reads J(1) (GS-fresh after I=1 solve) and J(3) (old).
    out, xv = run_case(rep, exe, "depth", "nu_depth")
    exp = np.empty((ND - 1, NCF))
    r_v = np.full(NCF, RV0)
    a_v = np.zeros(NCF)
    t_v = np.zeros(NCF, int)
    a_v[1:] = AV0
    t_v[1:] = np.arange(1, NCF)
    exp[0] = graph_solve(r_v, a_v, t_v)
    r2 = 0.5 - 0.1 * exp[0] - 0.1 * 0.9
    exp[1] = graph_solve(r2, a_v, t_v)
    exp[2] = graph_solve(r_v, a_v, t_v)
    rep.check(
        "depth: lagged TA/TC reads use the GS-fresh store",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )

    # --- fresh: uncaptured columns freeze -------------------------------
    out, xv = run_case(rep, exe, "fresh", "nu_fresh")
    exp = expect_simple(NCF, -1)
    rep.check(
        "fresh: captured columns solved",
        np.allclose(xv[: ND - 1, :5], exp[:, :5], rtol=1e-12,
                    atol=1e-14),
    )
    rep.check(
        "fresh: skipped columns frozen at the seeded field",
        np.allclose(xv[:, 5:], 3.0),
        "NU_FRESH=F slots must be untouched",
    )

    # --- sing: A=1 ring freezes ----------------------------------------
    out, xv = run_case(rep, exe, "sing", "nu_sing")
    rep.check(
        "sing: incompatible ring counted",
        kv(out, "NU_NSING") == ND - 1,
        f"NSING={kv(out, 'NU_NSING')}",
    )
    rep.check(
        "sing: ring members keep the old iterate",
        np.allclose(xv[:, 1], 4.0) and np.allclose(xv[:, 2], 5.0),
        "frozen members must not move",
    )

    # --- marg: near-marginal ring freezes in the band -------------------
    # a2*a3 = 0.9999 -> |1-A|=1e-4 < MOM_NU_RING_BAND: both members
    # must freeze at the seeded field, counted NMARG not NSING.
    out, xv = run_case(rep, exe, "marg", "nu_marg")
    rep.check(
        "marg: near-marginal ring counted per depth",
        kv(out, "NU_NMARG") == ND - 1,
        f"NMARG={kv(out, 'NU_NMARG')}",
    )
    rep.check(
        "marg: ring is not counted singular",
        kv(out, "NU_NSING") == 0,
    )
    rep.check(
        "marg: ring members keep the old iterate",
        np.allclose(xv[:, 1], 4.0) and np.allclose(xv[:, 2], 5.0),
        "frozen members must not move",
    )
    rep.check(
        "marg: uncoupled columns still solve",
        np.allclose(xv[: ND - 1, 0], RV0) and np.allclose(
            xv[: ND - 1, 3:], RV0
        ),
        "columns outside the ring are ordinary boundaries",
    )

    # --- margok: ring just above the band still solves ------------------
    # a2*a3 = 0.9409 -> |1-A|=0.0591 > band: exact solve, no freeze.
    out, xv = run_case(rep, exe, "margok", "nu_margok")
    xsol = (RV0 + 0.97 * RV0) / (1.0 - 0.97 * 0.97)
    exp = np.full((ND - 1, NCF), RV0)
    exp[:, 1] = xsol
    exp[:, 2] = xsol
    rep.check(
        "margok: ring solved exactly above the band",
        np.allclose(xv[: ND - 1], exp, rtol=1e-12, atol=1e-14),
        f"max diff {np.abs(xv[:ND-1]-exp).max():.2e}",
    )
    rep.check(
        "margok: nothing marginal, nothing capped",
        kv(out, "NU_NMARG") == 0 and kv(out, "NU_NCAP") == 0,
        f"NMARG={kv(out, 'NU_NMARG')} NCAP={kv(out, 'NU_NCAP')}"
        f" (solved {xsol:.3f} is below the 1e2 write cap)",
    )

    # --- cap: trust region clamps an explosive write --------------------
    # Column 5 solves to 1e6 against an incoming scale of 1: the
    # write must clamp to MOM_NU_WRITECAP*max(|XOLD|, 1e-4*FSCL)
    # = 1e2*max(1, 1e-4) = 100. Other columns solve freely.
    out, xv = run_case(rep, exe, "cap", "nu_cap")
    rep.check(
        "cap: one capped write per depth",
        kv(out, "NU_NCAP") == ND - 1,
        f"NCAP={kv(out, 'NU_NCAP')}",
    )
    rep.check(
        "cap: explosive cell clamped at the cap value",
        np.allclose(xv[: ND - 1, 4], 100.0, rtol=1e-12),
        f"col5 = {xv[:ND-1, 4]}",
    )
    rep.check(
        "cap: ordinary columns unaffected",
        np.allclose(xv[: ND - 1, 0], RV0) and np.allclose(
            xv[: ND - 1, 5:], RV0
        ),
    )

    # --- neg: negative chain results rejected ---------------------------
    out, xv = run_case(rep, exe, "neg", "nu_neg")
    rep.check(
        "neg: every negative result rejected",
        kv(out, "NU_NNEG") == (ND - 1) * NCF,
        f"NNEG={kv(out, 'NU_NNEG')}",
    )
    rep.check(
        "neg: field keeps the previous iterate",
        np.allclose(xv, 0.7),
        "rejected cells must not write",
    )

    # --- nonfin: bad rows freeze ---------------------------------------
    out, xv = run_case(rep, exe, "nonfin", "nu_nonfin")
    rep.check(
        "nonfin: dead-diagonal + nonfinite RHS frozen",
        kv(out, "NU_NFRZ") == 2,
        f"NFRZ={kv(out, 'NU_NFRZ')}",
    )
    rep.check(
        "nonfin: frozen cells keep the old value",
        np.isclose(xv[1, 3], 2.5) and np.isclose(xv[0, 5], 2.5),
    )
    rep.check(
        "nonfin: downwind readers consume the frozen value",
        np.isclose(xv[1, 4], RV0 + AV0 * 2.5)
        and np.isclose(xv[0, 6], RV0 + AV0 * 2.5),
        "readers of a frozen node solve with it as a constant",
    )

    # --- iter: composite-map dynamics on a marginal ring ----------------
    # Explosive-regime emulation: the 2-cycle's exact solve amplifies
    # the lagged feedback by c/(1-a) = 1000 per iteration. Theory
    # rates: conf1 (freeze band) rides the damped Jacobi step at
    # ~1.025/iter; conf2 (cap only) is bound by the trust region at
    # 100/iter; conf3 (unguarded) amplifies at 1000/iter. Measuring
    # the realised rates pins both the amplification mechanism and
    # the guard mechanics under iteration - the static cases above
    # cannot see this (the m15 failure was dynamical).
    r = run(exe, "iter", "nu_iter", cwd=WORKDIR)
    rep.check("iter: driver completed", "RESULT_OK" in r.stdout)
    trajs = {1: [], 2: [], 3: []}
    cens = {1: [], 2: [], 3: []}
    for line in r.stdout.splitlines():
        if line.startswith("ITER"):
            p = line.split()
            conf, it = int(p[1]), int(p[2])
            trajs[conf].append(float(p[3]))
            cens[conf].append(tuple(int(v) for v in p[5:9]))
    rep.check(
        "iter: all three configurations ran 40 iterations",
        all(len(trajs[c]) == 40 for c in (1, 2, 3)),
    )

    def slope(xs, i0, i1):
        """Log-linear growth rate per iteration over xs[i0:i1]."""
        ly = np.log(np.array(xs[i0:i1]))
        return np.exp(np.polyfit(np.arange(len(ly)), ly, 1)[0])

    r1 = slope(trajs[1], 20, 40)
    rep.check(
        "iter conf1 (freeze): ring drifts at the damped Jacobi rate",
        abs(r1 / 1.025 - 1.0) < 0.15,
        f"measured {r1:.4f}/iter vs theory 1.025",
    )
    r2 = slope(trajs[2], 4, 14)
    rep.check(
        "iter conf2 (cap only): growth bound by the trust region",
        abs(r2 / 100.0 - 1.0) < 0.10,
        f"measured {r2:.2f}/iter vs theory 100",
    )
    r3 = slope(trajs[3], 2, 10)
    rep.check(
        "iter conf3 (unguarded): growth at the amplification rate",
        abs(r3 / 1000.0 - 1.0) < 0.10,
        f"measured {r3:.2f}/iter vs theory 1000",
    )
    rep.check(
        "iter: freeze band fires every iteration in conf1 only",
        all(c[0] == ND - 1 for c in cens[1][5:])
        and all(c[0] == 0 for c in cens[2])
        and all(c[0] == 0 for c in cens[3]),
        f"conf1 NMARG/it={cens[1][10][0]}, "
        f"conf2={cens[2][10][0]}, conf3={cens[3][10][0]}",
    )
    rep.check(
        "iter: write cap binds on both ring members in conf2",
        all(c[1] == 2 * (ND - 1) for c in cens[2][2:15]),
        f"NCAP/it={cens[2][10][1]}",
    )
    rep.check(
        "iter: unguarded run never caps",
        all(c[1] == 0 for c in cens[3]),
    )
    rep.check(
        "iter: ordering of outcomes matches the m15 phenomenology",
        trajs[1][-1] < trajs[2][-1] < trajs[3][-1],
        f"x9(40): freeze {trajs[1][-1]:.2e} < cap {trajs[2][-1]:.2e}"
        f" < unguarded {trajs[3][-1]:.2e} - freeze drifts (stall),"
        " cap explodes slower, unguarded explodes fastest",
    )

    rep.summary()


if __name__ == "__main__":
    main()
