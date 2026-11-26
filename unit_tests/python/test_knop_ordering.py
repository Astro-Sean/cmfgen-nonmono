#!/usr/bin/env python3
"""Ordering-stall reproduction for the Knop non-monotonic solver.

Regression gate for the WP-1..WP-4 remediation: the manufactured ray
sits in the stale-donor regime the m15 stall diagnosis assigns to
SWEEP_FLD_ORD. Every interior row is "inserted" at WINX=0.95 toward
the deeper node, so the ordered march's compressing donor
    (1-WINX)*FLD(JPR1,L+NFM) + WINX*FLD(JPR2,L+NFM)
reads JPR1=J (own row, fresh on the row's descending-L pass) and
JPR2=J+1 (deeper row, a full pass stale on the ascending J march).
With |Q|/(chi+|Q|)=0.999 parked on the unmarched row the solve
contracts at ~WINX*|Q| per pass and cannot reach KN_RAY_TOL inside
KN_RAY_MAXPASS -- the production signature (stalled production rays
exit IXIT=3/6 at residual ~3e-2 over np 10-37). SWEEP_FLD computes
the whole donor column inside one pass, so the same operator is
block lower triangular and converges exactly.

Cases (drv_knop_ordering.f90, NI=140, NCF=200):
  ord    ordered sweep, production exits armed, matched spans (NFM=3)
  uno    same physics, FG_ORD_SWEEP=F
  ordx/unox   early exits off; clean 60-pass contraction measure
  ordb/unob   unmatched spans: drift 30 DLN bins, chi_hat*dz~30,
              |Q|=0.999 -- the drift-dominated algebra literally
  ordm/unom   matched mixed-sign ray (expanding fraction arg)
  ordmb/unomb unmatched mixed-sign ray

Measured expectations on the current default (ifort -O0):
  ord   NXIT=6 plateau exit at NPASS=25, RESID=4.12e-2, RHO_M=0.949
  uno   NXIT=1 clean at NPASS=3, RESID=4.1e-16
  ordx  NXIT=4 pass-cap at 60, RESID=1.588e-2 (the stall floor)
  unox  NXIT=1 clean at NPASS=3
  ordb  NXIT=6 plateau at NPASS=17, RESID=0.86
  unob  NXIT=1 clean at NPASS=3
  ordm  stalls like ord at every expanding fraction (rho~0.95)
  unom  ALSO fails to converge on the matched grid: at |Q|~1 on
        both signs plus A0~0.4 spatial coupling the alternating
        pass cannot isolate the stale share. The spec's "<=6
        passes" estimate does not survive at near-unity coupling
        on both signs -- recorded, not asserted away.
  ordmb/unomb  both converge (~13 passes): expanding rows anchor
        the chain because their own-row L-1 donors stay fresh --
        the stall needs an unbroken compressing chain.

Assertions:
  * uno/unox/unob: NXIT=1, RESID<KN_RAY_TOL, NPASS<=8.
  * ord/ordx/ordb/ordm: no clean tolerance exit (NXIT in {3,4,6}
    or RESID>tol at cap).
  * unom: NXIT in {3,4,6} (measured non-contraction, see above).
  * ordmb/unomb: NXIT in {1,2}, RESID<3e-3, NPASS<=30.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, F90DIR, TREE, WORKDIR, Report, compile_and_link_prod, run

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
    "new_main/mod_subs/iter_log_fmt_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
]

KN_RAY_TOL = 2.0e-3


def _f(x):
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def parse_out(text):
    out = {}
    for key in ("RESID", "RHO_M", "RHO_P"):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = _f(m.group(1)) if m else None
    for key in ("NPASS", "NXIT", "NMAX", "NBAD", "NTOT"):
        m = re.search(rf"{key}=\s*(\S+)", text)
        out[key.lower()] = int(m.group(1)) if m else None
    m = re.search(r"CLS_R=\s*(\S+)\s+(\S+)", text)
    out["cls"] = (int(m.group(1)), int(m.group(2))) if m else None
    out["ok"] = "RESULT_OK" in text
    out["fresh_m"] = [
        _f(m[1]) for m in re.findall(r"FRESH_M\s+(\d+)\s+(\S+)", text)
    ]
    out["fresh_p"] = [
        _f(m[1]) for m in re.findall(r"FRESH_P\s+(\d+)\s+(\S+)", text)
    ]
    out["text"] = text
    return out


def stalled(d):
    """The stall signature: no clean tolerance exit."""
    return d["nxit"] in (3, 4, 5, 6) or (
        d["resid"] is not None and d["resid"] > KN_RAY_TOL
    )


def main():
    rep = Report("KNOP ordering stall reproduction")
    bin_ = compile_and_link_prod(
        "knord", [SRC / h for h in HELPERS] + [F90DIR / "stub_vadat.f", F90DIR / "drv_knop_ordering.f90"]
    )

    data = {}
    for case, arg in (
        ("ord", ""), ("uno", ""), ("ordx", ""), ("unox", ""),
        ("ordb", ""), ("unob", ""),
        ("ordm", "2"), ("unom", "2"), ("ordm", "5"), ("unom", "5"),
        ("ordm", "8"), ("unom", "8"),
        ("ordmb", "2"), ("unomb", "2"), ("ordmb", "5"), ("unomb", "5"),
        ("ordmb", "8"), ("unomb", "8"),
        ("m2u", ""), ("m2m", "2"), ("m2m", "5"), ("m2m", "8"),
    ):
        r = run(bin_, case, arg, cwd=WORKDIR, check=False)
        d = parse_out(r.stdout)
        data[(case, arg)] = d
        rep.check(
            f"{case}{arg}: driver completed",
            r.returncode == 0 and d["ok"],
            f"rc={r.returncode}",
        )
        sys.stdout.write(
            f"    {case+arg:<7s} resid={d['resid']:.3e} npass={d['npass']}"
            f" nxit={d['nxit']} rho_m={d['rho_m']:.3f} cls={d['cls']}\n"
        )

    # -- T1: frequency-major sweep is the exact solve on the
    #    uniform-sign ray (block lower triangular).
    for case in ("uno", "unox", "unob"):
        d = data[(case, "")]
        rep.check(
            f"{case}: clean exit <=8 passes, resid < KN_RAY_TOL",
            d["nxit"] == 1 and d["resid"] < KN_RAY_TOL and d["npass"] <= 8,
            f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
        )

    # -- T2: the ordered march reproduces the stall on the same
    #    physics: no clean tolerance exit on any arm.
    for case in ("ord", "ordx", "ordb"):
        d = data[(case, "")]
        rep.check(
            f"{case}: stall signature (no clean tol exit)",
            stalled(d),
            f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
        )

    # -- T3: the default arm's residual sits at the measured stall
    #    floor (production stalled rays: 2.9-8.8e-2).
    d = data[("ord", "")]
    rep.check(
        "ord: residual floor ~3e-2+ and marginal contraction",
        d["resid"] > 5.0 * KN_RAY_TOL and d["rho_m"] > 0.9,
        f"resid={d['resid']:.3e} rho_m={d['rho_m']:.3f}",
    )

    # -- T4: exits-off arm reaches the pass cap with the residual
    #    floor bounded well above tolerance.
    d = data[("ordx", "")]
    rep.check(
        "ordx: 60-pass cap, residual floor > 5*tol",
        d["nxit"] == 4 and d["resid"] > 5.0 * KN_RAY_TOL,
        f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
    )

    # -- T5: ordered march stalls on the matched mixed ray at all
    #    expanding fractions (the stale compressing chain persists).
    for f in ("2", "5", "8"):
        d = data[("ordm", f)]
        rep.check(
            f"ordm{f}: ordered march stalls on mixed ray",
            stalled(d),
            f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
        )

    # -- T6: measured contradiction -- at |Q|~1 on both signs the
    #    alternating sweep does not isolate the stale share either;
    #    this is the lagged-map non-contraction the nbulk audit
    #    measured, not an ordering artifact. Asserted as observed
    #    non-clean so the fixture records it.
    for f in ("2", "5", "8"):
        d = data[("unom", f)]
        rep.check(
            f"unom{f}: alternating sweep also non-contracting at |Q|~1",
            stalled(d),
            f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
        )

    # -- T7: on the unmatched (drift-dominated, spatially decoupled)
    #    grid, expanding rows anchor the chain and BOTH marches
    #    converge -- the stall needs an unbroken compressing chain.
    for f in ("2", "5", "8"):
        for case in ("unomb", "ordmb"):
            d = data[(case, f)]
            rep.check(
                f"{case}{f}: converges (broken stale chain)",
                d["nxit"] in (1, 2)
                and d["resid"] < 3.0e-3
                and d["npass"] <= 30,
                f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
            )

    # -- T8: donor-offset census -- NFM=3 on the matched arm,
    #    offset-1 on the unmatched (spec's NFM=3 sits in the stall
    #    regime; offsets >=5 converge).
    rep.check(
        "NFM=3 on matched, offset-1 on unmatched",
        data[("ord", "")]["nmax"] == 3 and data[("ordb", "")]["nmax"] == 1,
        f"nmax ord={data[('ord','')]['nmax']} ordb={data[('ordb','')]['nmax']}",
    )

    # -- WP-2 freshness contract. KN_FRESH_FR records, per pass,
    #    the share of |Q| donor weight the march read from a
    #    not-yet-updated cell. On this fixture WINX=0.95 parks the
    #    compressing share on the deeper row J+1.
    #
    #    ord (mode 0): the I_M march fronts rows ascending, so
    #    JPR2=J+1 is unmarched -- stale share ~WINX. On I_P rows
    #    descend and J+1 is done -- stale ~0. The asymmetry is the
    #    measured defect: the ordered march's lag lives entirely
    #    in the I_M WINX share, not in "never waits a pass".
    d = data[("ord", "")]
    rep.check(
        "ord fresh: pass1 I_M stale ~WINX, I_P ~0",
        len(d["fresh_m"]) >= 1
        and d["fresh_m"][0] > 0.85
        and len(d["fresh_p"]) >= 1
        and d["fresh_p"][0] < 0.10,
        f"fresh_m={d['fresh_m'][:3]} fresh_p={d['fresh_p'][:3]}",
    )

    #    uno (mode 1): pass 1 ascends, so every compressing donor
    #    column is pending -- stale ~1. Pass 2 descends and the
    #    whole column is fresh -- stale ~0. The alternating pair
    #    is why the frequency-major sweep converges at all here.
    d = data[("uno", "")]
    rep.check(
        "uno fresh: pass1 stale ~1, pass2 fresh ~0",
        len(d["fresh_m"]) >= 2
        and d["fresh_m"][0] > 0.90
        and d["fresh_m"][1] < 0.10,
        f"fresh_m={d['fresh_m'][:4]}",
    )

    #    m2u (mode 2): uniform-compressing first pass descends, so
    #    pass 1 reads fresh donors where mode 1 spent a full stale
    #    pass. Pass 2 is the stale ascending half of the pair.
    d = data[("m2u", "")]
    rep.check(
        "m2u: clean exit <=8 passes, resid < KN_RAY_TOL",
        d["nxit"] == 1 and d["resid"] < KN_RAY_TOL and d["npass"] <= 8,
        f"nxit={d['nxit']} npass={d['npass']} resid={d['resid']:.3e}",
    )
    rep.check(
        "m2u fresh: pass1 (descending) ~0, pass2 stale ~1",
        len(d["fresh_m"]) >= 2
        and d["fresh_m"][0] < 0.10
        and d["fresh_m"][1] > 0.90,
        f"fresh_m={d['fresh_m'][:4]}",
    )

    #    m2m (mode 2 on the mixed ray): the dominant-sign rule
    #    keeps the mode-1 ascending start, so the case must match
    #    unom bitwise at every expanding fraction.
    for f in ("2", "5", "8"):
        dm, du = data[("m2m", f)], data[("unom", f)]
        rep.check(
            f"m2m{f}: bitwise identical to unom{f} (mode 1 on mixed)",
            dm["npass"] == du["npass"]
            and dm["nxit"] == du["nxit"]
            and dm["resid"] == du["resid"],
            f"m2m npass={dm['npass']} nxit={dm['nxit']} resid={dm['resid']:.3e}"
            f" vs unom npass={du['npass']} nxit={du['nxit']} resid={du['resid']:.3e}",
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
