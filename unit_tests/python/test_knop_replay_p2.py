"""Regression tests for the replay-only p2 override and mapd mode
(phase 2.30).

The p2=<v> keyed token lets an A/B replay pair differ only in
FG_ORD_P2DAMP without touching the VADAT reader or the production
default. mapd evaluates the undamped production map defect R(X)-X
at an arbitrary full-field candidate. Both are driver-side
features; no KNRDUMP deck is required to exercise the argument
layer, and the production source is pinned untouched.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, Report, compile_and_link_prod, run

TREE = F90DIR.parents[1]
DRV = F90DIR / "drv_knop_replay.f90"
PROD = TREE / "plane/fg_j_cmf_nonmono_v1.f"
RDVADAT = TREE / "new_main/mod_subs/rd_control_variables.f"

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
    "unit_tests/f90/stub_vadat.f",
]


def main():
    rep = Report("knop_replay_p2")
    src = DRV.read_text()
    prod = PROD.read_text()
    rdv = RDVADAT.read_text()

    # ---- source pins: p2= override -------------------------------
    rep.check(
        "p2 token scan covers every arg >=3",
        "COMMAND_ARGUMENT_COUNT()" in src
        and "CALL GETARG(IARG,ATMP)" in src
        and "ATMP(1:2) .EQ. 'p2'" in src,
    )
    rep.check(
        "malformed p2 rejected (no '=' after key)",
        "ATMP(3:3) .NE. '='" in src and "ERROR STOP 'p2: malformed override'" in src,
    )
    rep.check(
        "duplicate p2 rejected",
        "ERROR STOP 'p2: duplicate override'" in src,
    )
    rep.check(
        "bad p2 value rejected",
        "CASE DEFAULT" in src,
    )

    # the override must be the ONLY write to FG_ORD_P2DAMP besides
    # the deck restore and the mapd forced-off; any further write
    # would let an arm silently change more than the damper.
    assigns = re.findall(r"FG_ORD_P2DAMP\s*=", src)
    rep.check(
        "exactly three FG_ORD_P2DAMP assignments in driver",
        len(assigns) == 3,
        f"found {len(assigns)}",
    )
    rep.check(
        "deck value recorded before override",
        "LP2_DECK=FG_ORD_P2DAMP" in src and "IF(LP2_SET)FG_ORD_P2DAMP=LP2_VAL" in src,
    )
    rep.check(
        "override applied after deck restore, before solve",
        src.find("LP2_DECK=FG_ORD_P2DAMP") > src.find("FG_ORD_P2DAMP=BTEST(FBITS,1)")
        and src.find("LP2_DECK=FG_ORD_P2DAMP") < src.find("CALL KNOP_RAY_SOLVE_V1"),
    )
    rep.check(
        "effective P2 echoed in RESULT",
        "RESULT p2_deck=" in src and "p2_eff=" in src and "p2_ovr=" in src,
    )

    # ---- source pins: mapd mode ----------------------------------
    rep.check(
        "mapd rejects p2 override",
        "ERROR STOP 'mapd: p2 override not allowed'" in src,
    )
    rep.check(
        "mapd requires candidate argument",
        "ERROR STOP 'mapd: require candidate file'" in src,
    )
    rep.check(
        "mapd forces undamped unrelaxed single pass",
        re.search(
            r"MODE \.EQ\. 'mapd'[\s\S]*?KN_RAY_OMEGA=1\.0_LDP"
            r"[\s\S]*?KN_RAY_MAXPASS=1"
            r"[\s\S]*?FG_ORD_P2DAMP=\.FALSE\.",
            src,
        )
        is not None,
    )
    rep.check(
        "mapd seeds both lagged store columns from candidate",
        "I_P_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &\n      CAND_P(JP2,L)" in src
        and "I_M_DEPTH_NEXT(COMP_DEPTH_IDX(K),1,L-1)= &\n      CAND_M(JP2,L)" in src,
    )
    rep.check(
        "mapd validates candidate magic/shape",
        "NCMG .NE. 20261124" in src
        and "KNI_C .NE. NI" in src
        and "KCF_C .NE. NCF" in src,
    )
    rep.check(
        "mapd reports defect and writes R(X) dump",
        "RESULT mapd_dmax_p=" in src
        and "RESULT mapd_dmax_m=" in src
        and "RESULT mapd_nbad=" in src
        and "KN_MAPD_RX.bin" in src,
    )
    rep.check(
        "mapd defect uses production DCELL normalization",
        re.search(
            r"ABS\(KR_DBG_FLDP\(K,L\)-CAND_P\(K,L\)\)[\s\S]*?"
            r"MAX\(ABS\(CAND_P\(K,L\)\),1\.0E-4_LDP\*CAND_VP\(L\)\)",
            src,
        )
        is not None,
    )

    # ---- production untouched ------------------------------------
    rep.check(
        "production default still F",
        "DATA FG_ORD_P2DAMP /.FALSE./" in prod,
    )
    rep.check(
        "no replay-only tokens in production solver",
        "mapd" not in prod and re.search(r"(?<![a-zA-Z])p2=", prod) is None,
    )
    rep.check(
        "VADAT reader untouched by replay modes",
        "mapd" not in rdv and "p2=" not in rdv,
    )

    # ---- argument-layer behaviour (no deck needed) ---------------
    exe = compile_and_link_prod(
        "knray_replay",
        [TREE / h for h in HELPERS] + [DRV],
    )
    r = run(exe, "no_such_deck.bin", "base", "x", "p2junk", check=False)
    rep.check(
        "bare p2* token rejected before deck open",
        r.returncode != 0 and "p2" in r.stderr,
        f"rc={r.returncode} err={r.stderr.strip()[-80:]}",
    )
    r = run(exe, "no_such_deck.bin", "base", "x", "p2=", check=False)
    rep.check(
        "empty p2 value rejected",
        r.returncode != 0 and "p2" in r.stderr,
        f"rc={r.returncode}",
    )
    r = run(exe, "no_such_deck.bin", "base", "x", "p2=maybe", check=False)
    rep.check(
        "non-boolean p2 value rejected",
        r.returncode != 0 and "p2" in r.stderr,
        f"rc={r.returncode}",
    )
    r = run(exe, "no_such_deck.bin", "base", "x", "p2", check=False)
    rep.check(
        "keyless p2 token rejected",
        r.returncode != 0 and "p2" in r.stderr,
        f"rc={r.returncode}",
    )
    r = run(
        exe,
        "no_such_deck.bin",
        "base",
        "x",
        "p2=0",
        "p2=1",
        check=False,
    )
    rep.check(
        "duplicate p2 rejected",
        r.returncode != 0 and "duplicate" in r.stderr,
        f"rc={r.returncode} err={r.stderr.strip()[-80:]}",
    )
    r = run(exe, "no_such_deck.bin", "base", "x", "p2=1", check=False)
    rep.check(
        "valid p2=1 passes parse (fails later at deck open)",
        "deck_open_fail" in r.stdout,
        f"rc={r.returncode} out={r.stdout.strip()[-80:]}",
    )
    r = run(exe, "no_such_deck.bin", "omf", "sched.txt", "3", "p2=1", check=False)
    rep.check(
        "p2= composes with omf arg layout",
        "deck_open_fail" in r.stdout,
        f"rc={r.returncode}",
    )

    for sweep in ("ordered", "frequency"):
        r = run(exe, "no_such_deck.bin", "sweep", sweep, check=False)
        rep.check(f"valid sweep {sweep} passes parse", "deck_open_fail" in r.stdout)
    for sweep in ("", "invalid"):
        r = run(exe, "no_such_deck.bin", "sweep", sweep, check=False)
        rep.check(
            f"invalid sweep {sweep!r} rejected before deck open",
            r.returncode != 0 and "sweep:" in r.stderr,
        )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
