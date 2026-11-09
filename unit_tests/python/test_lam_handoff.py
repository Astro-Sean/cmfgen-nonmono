#!/usr/bin/env python3
"""test_lam_handoff -- RD_LAMBDA -> full-iteration handoff hardening.

Compiles the production LAM_HANDOFF_MOD (new_main/mod_subs/
lam_handoff_mod.f90) with drv_lam_handoff and asserts on the printed
verdicts, then pins the production integration points by source
inspection:

  * LAM_FIELD_GATE: held/frozen block with or without a lagged
    field; raw immaturity and the N_MATURE/GS_REL_STABLE window only
    count on lagged-field routes (KNOP, SOB+LAM_ES=F, SOB_FG_KNOP).
  * LAM_STABLE_STEP: consecutive quiet calls earn the count; force,
    luminosity block, over-demand and field-not-ready reset it; a
    resacc-accepted probe pins the count only while the field is
    ready (the old latch re-armed after a rejected probe).
  * LAM_HANDOFF_OK: stock path keeps the single-iteration
    MAXCH_CTRL < 50 test untouched; non-mono adds stable count AND
    field readiness.
  * Source pins: the HAS_LAGGED_FIELD gate around the maturity
    getters, the SOB-mom route in HAS_LAGGED_FIELD, the production
    call sites, the BA_PROBE_LAM > 0 arm, the rejected-probe latch
    clear, the DO_LAM_AUTO T-hold arm, the blocking diagnostic and
    the check_param_consistency warning.

Run: python3 unit_tests/python/test_lam_handoff.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common

REP = common.Report("test_lam_handoff (RD_LAMBDA handoff hardening)")

WORK = common.WORKDIR / "lam_handoff"
WORK.mkdir(parents=True, exist_ok=True)

binpath = common.compile_and_link_prod(
    "drv_lam_handoff",
    [
        common.TREE / "tools" / "set_kind_module.f",
        common.TREE / "new_main" / "mod_subs" / "lam_handoff_mod.f90",
        common.F90DIR / "drv_lam_handoff.f90",
    ],
)
r = common.run(binpath, cwd=WORK)
REP.check("driver ran", r.returncode == 0)

vals = {}
for ln in r.stdout.splitlines():
    m = re.match(r"^\s*([A-Za-z0-9_]+)=(.+)$", ln.strip())
    if m:
        vals[m.group(1)] = m.group(2).strip()


def tf(name):
    return vals.get(name, "").upper().startswith("T")


def cnt(name):
    return int(vals.get(name, "-1").split()[0])


def prov(name):
    return "proven=T" in vals.get(name, "").replace(" ", "")


# ---- LAM_FIELD_GATE -----------------------------------------------
# No lagged field: only held/frozen block; lagged terms unjudged.
REP.check("gate: no-lagged clean passes", not tf("G_nolag_clean"))
REP.check("gate: no-lagged held blocks", tf("G_nolag_held"))
REP.check("gate: no-lagged frozen blocks", tf("G_nolag_freeze"))
REP.check(
    "gate: no-lagged ignores raw-immature input", not tf("G_nolag_rawimm")
)
# Lagged-field routes: every term blocks independently.
REP.check("gate: lagged held blocks", tf("G_lag_held"))
REP.check("gate: lagged frozen blocks", tf("G_lag_freeze"))
REP.check("gate: lagged raw-immature blocks", tf("G_lag_rawimm"))
REP.check("gate: window 0/2 blocks", tf("G_lag_window0of2"))
REP.check("gate: window 1/2 blocks", tf("G_lag_window1of2"))
REP.check("gate: window 2/2 clears", not tf("G_lag_window2of2"))
REP.check("gate: window 0/1 blocks (first call)", tf("G_lag_window0of1"))
REP.check("gate: window 1/1 clears", not tf("G_lag_window1of1"))

# ---- LAM_STABLE_STEP ----------------------------------------------
REP.check(
    "step: 5 consecutive quiet calls earn+prove",
    cnt("S_earn5") == 5 and prov("S_earn5"),
)
REP.check(
    "step: 4 quiet calls do not prove",
    cnt("S_earn4") == 4 and not prov("S_earn4"),
)
REP.check(
    "step: demand >= LAM_EXIT resets",
    cnt("S_reset_demand") == 0 and not prov("S_reset_demand"),
)
REP.check(
    "step: forced lambda resets",
    cnt("S_reset_force") == 0 and not prov("S_reset_force"),
)
REP.check(
    "step: luminosity block resets",
    cnt("S_reset_lum") == 0 and not prov("S_reset_lum"),
)
REP.check(
    "step: field-not-ready resets (SOB-mom fix)",
    cnt("S_reset_field") == 0 and not prov("S_reset_field"),
)
REP.check(
    "step: field-not-ready revokes a latched proof",
    cnt("S_unlatch_field") == 0 and not prov("S_unlatch_field"),
)
REP.check(
    "step: frozen call is neutral (no earn, no reset)",
    cnt("S_frozen_neutral") == 3 and not prov("S_frozen_neutral"),
)
REP.check(
    "step: accepted probe pins count at threshold",
    cnt("S_pin_clean") == 5 and prov("S_pin_clean"),
)
REP.check(
    "step: pin blocked by not-ready field",
    cnt("S_pin_field_blocked") == 0 and not prov("S_pin_field_blocked"),
)
REP.check(
    "step: pin blocked by forced lambda",
    cnt("S_pin_force_blocked") == 0 and not prov("S_pin_force_blocked"),
)
REP.check(
    "step: pin blocked by luminosity block",
    cnt("S_pin_lum_blocked") == 0 and not prov("S_pin_lum_blocked"),
)

# ---- LAM_HANDOFF_OK ------------------------------------------------
# Stock: MAXCH_CTRL < 50 alone decides (bit-identical to upstream).
REP.check("handoff: stock releases below 50", tf("H_stock_clean"))
REP.check("handoff: stock holds at exactly 50", not tf("H_stock_at50"))
REP.check("handoff: stock holds on jfnk-unconverged", not tf("H_stock_jfnk"))
REP.check(
    "handoff: stock ignores field/stable inputs", tf("H_stock_ignores")
)
REP.check("handoff: stock holds without RD_LAMBDA", not tf("H_stock_nordlam"))
REP.check("handoff: stock holds without DO_LAM_AUTO", not tf("H_stock_noauto"))
# Non-mono: stable AND field-ready required on top of stock terms.
REP.check("handoff: non-mono clean releases", tf("H_nm_clean"))
REP.check("handoff: non-mono unstable holds", not tf("H_nm_unstable"))
REP.check("handoff: non-mono field-not-ready holds", not tf("H_nm_field_nrdy"))
REP.check("handoff: non-mono holds at exactly 50", not tf("H_nm_at50"))
REP.check("handoff: non-mono holds on jfnk-unconverged", not tf("H_nm_jfnk"))
REP.check(
    "handoff: non-mono releases on quiet demand with ready field",
    tf("H_nm_quiet_dem"),
)

# ---- Production integration source pins ---------------------------
def src(rel):
    return (common.TREE / rel).read_text()


sfp = src("new_main/mod_subs/solve_for_pops.f")
cms = src("new_main/cmfgen_sub.f")
rcv = src("new_main/mod_subs/rd_control_variables.f")
cpc = src("new_main/mod_subs/check_param_consistency.f")
cvm = src("new_main/mod_subs/control_variable_mod.f")

REP.check(
    "pin: maturity getters gated on HAS_LAGGED_FIELD",
    "IF(HAS_LAGGED_FIELD)THEN\n\t  CALL GET_FG_GS_MATURITY" in sfp,
)
REP.check(
    "pin: HAS_LAGGED_FIELD covers the SOB moment route",
    re.search(
        r"HAS_LAGGED_FIELD=DO_NONMONO_KNOP\s*\.OR\.\s*DO_SOB_FG_KNOP\s*\.OR\."
        r"\s*\n\s*1\s*\(DO_NONMONO_SOB\s*\.AND\.\s*\.NOT\.\s*USE_LAM_ES\)",
        rcv,
    )
    is not None,
)
REP.check(
    "pin: field gate computed from production call state",
    "LAM_FIELD_NOT_READY=LAM_FIELD_GATE(HAS_LAGGED_FIELD" in sfp,
)
REP.check(
    "pin: RD_LAMBDA counter runs LAM_STABLE_STEP",
    "CALL LAM_STABLE_STEP(FORCE_LAMBDA,LUM_BLOCKED,GS_FREEZE," in sfp,
)
REP.check(
    "pin: cmfgen_sub handoff runs LAM_HANDOFF_OK",
    "ELSE IF(LAM_HANDOFF_OK(DO_LAMBDA_AUTO,RD_LAMBDA,MAXCH_CTRL," in cms,
)
REP.check(
    "pin: probe arm still requires BA_PROBE_LAM > 0",
    "BA_PROBE_LAM .GT. 0" in sfp and "BA_PROBE_LAM=0" not in sfp.split(
        "BA_PROBE_LAM .GT. 0"
    )[0][-200:],
)
REP.check(
    "pin: rejected probe clears the accept latch",
    "IF(RD_LAMBDA)BA_RDLAM_ACCEPT=.NOT. BA_RAW_REJECTED" in sfp,
)
REP.check(
    "pin: handoff clears the accept latch",
    re.search(r"RD_LAMBDA=\.FALSE\.[\s\S]{0,600}?BA_RDLAM_ACCEPT=\.FALSE\.", cms)
    is not None,
)
REP.check(
    "pin: handoff arms T-hold when N_FULL_TREL > 0",
    re.search(
        r"RD_LAMBDA=\.FALSE\.[\s\S]{0,900}?IF\(N_FULL_TREL \.GT\. 0\)THEN"
        r"[\s\S]{0,200}?T_HELD_PHASE=\.TRUE\.",
        cms,
    )
    is not None,
)
REP.check(
    "pin: blocked-handoff diagnostic names the blocking term",
    "LAM_HANDOFF: stable count met" in sfp
    and "lam_field_hold" in sfp,
)
REP.check(
    "pin: field-not-ready added to probe arm",
    ".NOT. LAM_FIELD_NOT_READY" in sfp,
)
REP.check(
    "pin: consistency warning for dead BA_PROBE_LAM",
    "BA_PROBE_LAM .GT. 0 .AND." in cpc and "cannot fire" in cpc,
)
REP.check(
    "pin: LAM_FIELD_NOT_READY / HAS_LAGGED_FIELD declared",
    "LOGICAL LAM_FIELD_NOT_READY" in cvm
    and "LOGICAL HAS_LAGGED_FIELD" in cvm,
)
REP.check(
    "pin: no_refresh token gated to FG routes",
    "IF(DO_NONMONO_KNOP .OR. DO_SOB_FG_KNOP)\n\t1      CALL APPEND_ITER_REASON(ITER_REASON,'no_refresh')"
    in sfp,
)

sys.exit(REP.summary())
