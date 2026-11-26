"""WP-6 forcing-term refresh tolerance (KN_FORCING_TOL).

drv_knop_frc.f90 drives the production path the refresh sees:
SET_FG_FORCING_TOL installs the flag/knobs, SET_FG_FRC_RES hands
in the last population demand, KN_FRC_UPDATE runs one Eisenstat &
Walker 1996 (SIAM J. Sci. Comput. 17, 16-32) Choice-2 step:

  eta_k = gamma * (F_k/F_{k-1})^alpha, floored at
  gamma*eta_{k-1}^alpha when that exceeds 0.1, capped at 0.9;
  KN_RAY_TOL_EFF = clamp(eta*TMAX/0.9, TMIN, min(TMAX, CAP)).

The fixture uses the production defaults gamma=0.9, alpha=2,
TMIN=1e-4, TMAX=1e-2 and the I2 cap 0.1*GS_MATURE_TOL=5e-3 with
KN_RAY_TOL=2e-3 as the flag-off base.

Audit: docs/audits/KNOP_REMEDIATION_WP6_2026-11-26.md.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import SRC, F90DIR, WORKDIR, Report, compile_and_link_prod, run

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

REP = Report("test_knop_frc (WP-6 KN_FORCING_TOL)")

BIN = compile_and_link_prod(
    "knofrc",
    [SRC / h for h in HELPERS]
    + [F90DIR / "stub_vadat.f", F90DIR / "drv_knop_frc.f90"],
)


def _f(x):
    if "E" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


def steps(text):
    out = []
    for m in re.finditer(
        r"FRC_STEP\s+(\d+)\s+ETA=\s*(\S+)\s+RAT=\s*(\S+)\s+SG=\s*(\d+)\s+TOL=\s*(\S+)",
        text,
    ):
        out.append(
            (int(m.group(1)), _f(m.group(2)), _f(m.group(3)),
             int(m.group(4)), _f(m.group(5)))
        )
    return out


def scalar(text, key):
    m = re.search(rf"{key}=\s*(\S+)", text)
    return _f(m.group(1)) if m else None


def close(got, want, tag, tol=1e-9):
    REP.check(tag, got is not None and abs(got - want) <= tol * max(1.0, abs(want)),
              detail=f"got {got} want {want}")


# ---- seed: first armed call has no ratio; seeds FPRV, holds base ----
r = run(BIN, "seed")
REP.check("seed ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
close(scalar(r.stdout, "FRC_ETA"), -1.0, "seed eta unset")
close(scalar(r.stdout, "FRC_FPRV"), 1.0, "seed fprv stored")
close(scalar(r.stdout, "FRC_TOL"), 2.0e-3, "seed tol at base")

# ---- seq: contracting demand tightens eta and tol monotonically ----
# ratios 0.5, 0.5, 0.4 -> eta = 0.9*rat^2 = 0.225, 0.225, 0.144
# tol = eta*(1e-2/0.9) -> 2.5e-3, 2.5e-3, 1.6e-3 (cap 5e-3 not hit)
r = run(BIN, "seq")
REP.check("seq ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
st = steps(r.stdout)
REP.check("seq has 4 steps", len(st) == 4, detail=f"{len(st)}")
if len(st) == 4:
    # step 1 seeds: eta=-1, tol=base
    close(st[0][1], -1.0, "seq step1 eta seeds")
    close(st[0][4], 2.0e-3, "seq step1 tol base")
    close(st[1][2], 0.5, "seq step2 ratio")
    close(st[1][1], 0.225, "seq step2 eta")
    close(st[1][4], 2.5e-3, "seq step2 tol")
    close(st[3][1], 0.144, "seq step4 eta")
    close(st[3][4], 1.6e-3, "seq step4 tol")
    REP.check("seq ladder never inverts",
              all(st[i][4] <= st[i - 1][4] + 1e-12 for i in (2, 3)))

# ---- rise: flat demand -> eta hits 0.9; tol lands on the I2 cap ----
r = run(BIN, "rise")
REP.check("rise ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
st = steps(r.stdout)
if len(st) == 4:
    close(st[3][1], 0.9, "rise eta capped at 0.9")
    close(st[3][4], 5.0e-3, "rise tol lands on I2 cap not TMAX")

# ---- sflo: safeguard floor raises eta instead of dropping ----
# eta_prev=0.6 -> floor = 0.9*0.36 = 0.324 > 0.1; raw eta =
# 0.9*0.09 = 0.081 < floor -> eta = 0.324, SG=1.
r = run(BIN, "sflo")
REP.check("sflo ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
close(scalar(r.stdout, "FRC_ETA"), 0.324, "sflo eta on safeguard floor")
close(scalar(r.stdout, "FRC_RAT"), 0.3, "sflo ratio")
REP.check("sflo safeguard fired",
          re.search(r"FRC_SG=\s*1\b", r.stdout) is not None)
close(scalar(r.stdout, "FRC_TOL"), 0.324 * (1.0e-2 / 0.9), "sflo tol")

# ---- tmin: fast contraction hits the TMIN floor ----
r = run(BIN, "tmin")
REP.check("tmin ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
st = steps(r.stdout)
if len(st) == 4:
    # step2: rat=0.05 -> eta=0.9*2.5e-3=2.25e-3 -> tol=2.5e-5 < 1e-4
    close(st[1][1], 0.9 * 0.05**2, "tmin step2 eta")
    close(st[1][4], 1.0e-4, "tmin step2 tol on floor")
    # step4: rat=0.5 -> raw eta 0.225, but safeguard floor
    # 0.9*eta_3^2 binds when eta_3 small enough? eta_3 = 0.9*0.2^2
    # = 0.036 -> floor 0.9*0.036^2 = 1.1664e-3 < 0.1 -> no guard.
    close(st[3][1], 0.225, "tmin step4 eta")
    close(st[3][4], 2.5e-3, "tmin step4 tol recovers")

# ---- noup: no pending demand leaves the tolerance at base ----
r = run(BIN, "noup")
REP.check("noup ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
close(scalar(r.stdout, "FRC_TOL"), 2.0e-3, "noup tol stays base")
close(scalar(r.stdout, "FRC_ETA"), -1.0, "noup eta stays unset")

# ---- getr: GET_FG_FRC exports the last update verbatim ----
# sequence 1.0 -> 0.4 -> 0.16 -> 0.064: rat 0.4 each step, eta =
# 0.9*0.16 = 0.144, tol = 0.144*(1e-2/0.9) = 1.6e-3, no safeguard
# (floor 0.9*0.144^2 = 1.87e-2 < 0.1).
r = run(BIN, "getr")
REP.check("getr ran", r.returncode == 0 and "RESULT_OK" in r.stdout)
st = steps(r.stdout)
close(scalar(r.stdout, "GET_ETA"), 0.144, "getr eta")
close(scalar(r.stdout, "GET_RAT"), 0.4, "getr ratio")
REP.check("getr sg", re.search(r"GET_SG=\s*0\b", r.stdout) is not None)
close(scalar(r.stdout, "GET_TOL"), 1.6e-3, "getr tol")

sys.exit(REP.summary())
