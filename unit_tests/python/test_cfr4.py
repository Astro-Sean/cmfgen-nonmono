#!/usr/bin/env python3
"""CMF_FORMAL_REL_V4 coastal-ray + N_STORE guards
(plane/cmf_formal_rel_v4.f vs the frozen reference).

Alterations under test, compared directly against the pristine
reference source built with production flags:

1. N_STORE floor (INNER_BND_METH='HOLLOW').  The hollow-core
   frequency-wrap store is sized NINT(2*V(ND)); a coasting inner
   velocity (V(ND)<0.5 km/s - exactly the regime a stalled
   non-monotonic inner region reaches) truncates it to 0, and the
   first MOD(CUR_LOC+1,N_STORE) then faulted.  Prod floors N_STORE
   at 1.  Case 'coast' sets V(ND)=0.1: upstream dies in the wrap,
   prod completes.

2. Coastal SQRT guards.  MU at the boundaries is evaluated as
   SQRT(MAX(0,(R-P)*(R+P)))/R; a ray whose P sits numerically above
   R(1) or R(ND) gave upstream a NaN where prod takes the 0 branch.
   Case 'nudge' pushes the outermost tangent ray 1e-7 above R(1):
   prod completes with finite JNU/FEDD (the characteristics clamp
   and the SQRT guard work together); upstream propagates NaN.

Case 'ctl' (HOLLOW boundary, normal velocity law) must match the
reference bit-for-bit - the guards are inert on the normal path.

Note the coastal trigger is shared with characteristics_v2's own
clamp: 'nudge' exercises the coastal-guard family end to end rather
than isolating the SQRT line alone.
"""
import re
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F90DIR, REF, TREE, Report, compile_and_link_prod, run

DRV = F90DIR / "drv_cfr4.f90"
HELPERS = [
    TREE / "plane" / "define_grid_v2.f",
    TREE / "plane" / "mod_space_grid_v2.f",
    TREE / "plane" / "solve_cmf_formal_v2.f",
    TREE / "plane" / "optdepth_v4.f",
    TREE / "plane" / "rel_variables.f",
    TREE / "plane" / "j_weight.f",
    TREE / "plane" / "h_weight.f",
    TREE / "plane" / "k_weight.f",
    TREE / "plane" / "n_weight.f",
    TREE / "subs" / "phys_con.f",
    TREE / "subs" / "tuvghd_rh.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "tune.f",
    TREE / "tools" / "get_lu.f",
    TREE / "tools" / "set_kind_module.f",
    TREE / "plane" / "velocity_law.f",
    TREE / "plane" / "derivs.f",
    TREE / "plane" / "derivr.f",
    TREE / "plane" / "runge_kutta.f",
]
PROD_SOURCES = [
    TREE / "plane" / "mod_ray_mom_store.f",
    TREE / "plane" / "characteristics_v2.f",
    TREE / "plane" / "cmf_formal_rel_v4.f",
    TREE / "subs" / "mon_int_funs_v2.f",
]
REF_SOURCES = [
    # REF modules/interpolators first: their interfaces must shadow
    # the production .mod files when the REF body compiles.
    REF / "plane" / "mod_ray_mom_store.f",
    REF / "plane" / "characteristics_v2.f",
    REF / "plane" / "cmf_formal_rel_v4.f",
    REF / "subs" / "mon_int_funs_v2.f",
]


def floats(tag, out):
    m = re.search(tag + r"=\s*((?:[-+0-9.Ee]+\s*)+)", out)
    return [float(x) for x in m.group(1).split()]


def finite(vals):
    return all(v == v and abs(v) != float("inf") for v in vals)


DRV_G = F90DIR / "drv_cfr4g.f90"
HELPERS_G = [
    TREE / "plane" / "define_grid_v2_gam.f",
    TREE / "plane" / "mod_space_grid_v2_gam.f",
    TREE / "plane" / "solve_cmf_formal_v3_gam.f",
    TREE / "plane" / "optdepth_v4_gam.f",
    TREE / "plane" / "rel_variables.f",
    TREE / "plane" / "j_weight.f",
    TREE / "plane" / "h_weight.f",
    TREE / "plane" / "k_weight.f",
    TREE / "plane" / "n_weight.f",
    TREE / "subs" / "phys_con.f",
    TREE / "subs" / "tuvghd_rh.f",
    TREE / "subs" / "mon_interp.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "tune.f",
    TREE / "unix" / "set_line_buffering.f",
    TREE / "tools" / "get_lu.f",
    TREE / "tools" / "set_kind_module.f",
    TREE / "plane" / "velocity_law.f",
    TREE / "plane" / "derivs.f",
    TREE / "plane" / "derivr.f",
    TREE / "plane" / "runge_kutta.f",
    TREE / "new_main" / "mod_subs" / "gam_mu_mod.f90",
    TREE / "new_main" / "mod_subs" / "mod_rd_gamray_cntrl_variables.f90",
]
PROD_G = [
    TREE / "plane" / "mod_ray_mom_store.f",
    TREE / "plane" / "characteristics_v2_gam.f",
    TREE / "plane" / "cmf_formal_rel_v4_gam.f",
    TREE / "subs" / "mon_int_funs_v2.f",
]
REF_G = [
    REF / "plane" / "mod_ray_mom_store.f",
    REF / "plane" / "characteristics_v2_gam.f",
    REF / "plane" / "cmf_formal_rel_v4_gam.f",
    REF / "subs" / "mon_int_funs_v2.f",
]


def main():
    rep = Report("CMF_FORMAL_REL_V4 coastal + N_STORE guards")
    prod = compile_and_link_prod(
        "drv_cfr4", PROD_SOURCES + HELPERS + [DRV])
    ref = compile_and_link_prod(
        "drv_cfr4_ref", REF_SOURCES + HELPERS + [DRV])

    with tempfile.TemporaryDirectory() as d:
        p = run(prod, "ctl", cwd=d).stdout
        r = run(ref, "ctl", cwd=d).stdout
        rep.check("ctl: both binaries complete",
                  "DONE_CFR4" in p and "DONE_CFR4" in r)
        rep.check("ctl: N_STORE allocated (prod log)",
                  "N_STORE=" in p)
        pj = floats("JNU", p)
        rep.check("ctl: prod JNU finite", finite(pj))
        rep.check("ctl: prod JNU identical to ref (guards inert)",
                  pj == floats("JNU", r))
        rep.check("ctl: prod FEDD identical to ref",
                  floats("FEDD", p) == floats("FEDD", r))

        # V(ND)=0.1 truncates upstream's N_STORE to 0.  Prod prints
        # N_STORE=1 and moves on (the run later stops on an unrelated
        # frequency-range check for the unrealistic coast velocity -
        # what matters is that it passed the store); upstream dies
        # right at MOD(CUR_LOC+1,0) with an integer divide by zero.
        q = run(prod, "coast", cwd=d, check=False)
        rep.check("coast: prod floors N_STORE at 1",
                  "N_STORE=           1" in q.stdout)
        qr = run(ref, "coast", cwd=d, check=False)
        rep.check("coast: ref crashes in the wrap (divide by zero)",
                  "DONE_CFR4" not in qr.stdout
                  and ("divide by zero" in (qr.stderr + qr.stdout)
                       or qr.returncode != 0))

        # Coastal outer ray: the ray itself is degenerate (never
        # reaches storage) so JNU is unaffected, but the boundary-mu
        # diagnostic writes the raw SQRT result - 0.0 in prod, NaN
        # upstream.
        pd, rd = Path(d) / "nudge_p", Path(d) / "nudge_r"
        pd.mkdir()
        rd.mkdir()
        q = run(prod, "nudge", cwd=pd, check=False)
        qr = run(ref, "nudge", cwd=rd, check=False)
        rep.check("nudge: both binaries complete",
                  "DONE_CFR4" in q.stdout and "DONE_CFR4" in qr.stdout)
        rep.check("nudge: prod JNU identical to ref on the "
                  "degenerate ray",
                  floats("JNU", q.stdout) == floats("JNU", qr.stdout))
        p_mu = (pd / "MU_VALUE_CHK").read_text()
        r_mu = (rd / "MU_VALUE_CHK").read_text()
        rep.check("nudge: prod mu diagnostic is 0.0 (SQRT floored)",
                  "0.000000E+00" in p_mu)
        rep.check("nudge: ref mu diagnostic is NaN",
                  "NaN" in r_mu)

    # GAM variant: identical guards in the gamma-relativistic branch,
    # exercised through the same three cases.
    gam_prod = compile_and_link_prod(
        "drv_cfr4g", PROD_G + HELPERS_G + [DRV_G])
    gam_ref = compile_and_link_prod(
        "drv_cfr4g_ref", REF_G + HELPERS_G + [DRV_G])

    with tempfile.TemporaryDirectory() as d:
        # The gam grid writer touches ./data/GAM_MU_GRID.
        for sub in ("gctl_p", "gctl_r", "gc_p", "gc_r", "np", "nr"):
            (Path(d) / sub / "data").mkdir(parents=True)

        p = run(gam_prod, "ctl", cwd=Path(d) / "gctl_p").stdout
        r = run(gam_ref, "ctl", cwd=Path(d) / "gctl_r").stdout
        rep.check("gam ctl: both complete",
                  "DONE_CFR4G" in p and "DONE_CFR4G" in r)
        rep.check("gam ctl: prod JNU identical to ref",
                  floats("JNU", p) == floats("JNU", r))

        q = run(gam_prod, "coast", cwd=Path(d) / "gc_p", check=False)
        rep.check("gam coast: prod floors N_STORE at 1",
                  "N_STORE=           1" in q.stdout)
        qr = run(gam_ref, "coast", cwd=Path(d) / "gc_r", check=False)
        rep.check("gam coast: ref crashes in the wrap",
                  "DONE_CFR4G" not in qr.stdout
                  and ("divide by zero" in (qr.stderr + qr.stdout)
                       or qr.returncode != 0))

        pd, rd = Path(d) / "np", Path(d) / "nr"
        q = run(gam_prod, "nudge", cwd=pd, check=False)
        qr = run(gam_ref, "nudge", cwd=rd, check=False)
        rep.check("gam nudge: both complete",
                  "DONE_CFR4G" in q.stdout
                  and "DONE_CFR4G" in qr.stdout)
        rep.check("gam nudge: prod mu diagnostic 0.0, ref NaN",
                  "0.000000E+00" in (pd / "MU_VALUE_CHK").read_text()
                  and "NaN" in (rd / "MU_VALUE_CHK").read_text())

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
