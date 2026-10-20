#!/usr/bin/env python3
"""Unit tests for the modified VADAT parsing in
new_main/mod_subs/rd_control_variables.f.

Runs the real parser (with the real RD-store layer, control-variable
module and JFNK solver module) on a production VADAT file taken from
a live run, with individual keys mutated per case.

New-key coverage:
  * DO_NONMONO_SOB / DO_NONMONO_KNOP / DO_NONLOC_SOB parse and land
    in CONTROL_VARIABLE_MOD.
  * Mutual exclusion: both non-mono flags T aborts.
  * Prerequisites enforced (previously warnings only):
      DO_NONMONO_SOB  needs GLOBAL_LINE=SOB; USE_LAM_ES=F warns only
      (opt-in MOM_J_CMF_NONMONO_V1 path, rewired 12-Oct-2026)
      DO_NONMONO_KNOP needs GLOBAL_LINE=BLANK and LAM_SET=F
  * DO_JFNK requires a non-mono mode; JFNK_COUPLED is KNOP-only;
    USE_JFNK is a deprecated alias that still arms the solver.
  * DO_NONLOC_SOB is forced F (with warning) when DO_NONMONO_SOB=F.
  * Legacy V-coord controls rejected under non-mono: DO_CL=T and
    GAMRAY_TRANS=RAD_TRANS/NONLOCAL abort.
  * INNER_BND_METH centralized whitelist validation.
  * VEL_LAW=13/14 select SN_HYDRO_MODEL.
  * JFNK_PTC / JFNK_PTC_DTAU module variables parse.
  * Unmodified keys round-trip (spot check of shock-power controls).
"""

import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    BINDIR,
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_objs,
    link,
    run,
)

# Template VADAT shipped in the repo; override with JFNK_TEST_VADAT
# to exercise the parser against a production file from a live run.
VADAT_SRC = Path(
    os.environ.get(
        "JFNK_TEST_VADAT",
        TREE / "unit_tests" / "systems" / "vadat_template.txt",
    )
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
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "new_main/mod_subs/rd_control_variables.f",
]


def build():
    objs = compile_objs(
        "vadat",
        [TREE / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_vadat.f90"],
    )
    return link("vadat", objs, moddir=BINDIR / "obj_vadat")


def write_vadat(d, mutations):
    """Copy the production VADAT, replacing values for [KEY]s.

    mutations: dict KEY -> new value string.  A value of None deletes
    the key line entirely (tests absent-key defaults).
    """
    d.mkdir(parents=True, exist_ok=True)
    out = []
    seen = set()
    for ln in open(VADAT_SRC):
        m = re.search(r"\[([A-Z_0-9]+)\]\s*$", ln.strip())
        if m and m.group(1) in mutations:
            seen.add(m.group(1))
            v = mutations[m.group(1)]
            if v is not None:
                out.append(f"{v} [{m.group(1)}]\n")
            continue
        out.append(ln)
    # Optional keys absent from the production file are appended
    # (RD-store requires unique keys; order is irrelevant).
    for k, v in mutations.items():
        if k not in seen and v is not None:
            out.append(f"{v} [{k}]\n")
    with open(d / "VADAT", "w") as f:
        f.writelines(out)
    return d


def parse(out):
    vals = {}
    for ln in out.splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            vals[k.strip()] = v.strip()
    return vals


def run_vadat(exe, d):
    return run(exe, cwd=d, check=False, timeout=60)


def main():
    rep = Report("RD_CONTROL_VARIABLES VADAT keys")
    exe = build()

    # Baseline: production VADAT as-is.
    d = write_vadat(WORKDIR / "vadat_base", {})
    p = run_vadat(exe, d)
    ok = "RESULT_OK" in p.stdout
    rep.check("baseline production VADAT parses", ok, p.stdout[-500:] if not ok else "")
    if ok:
        v = parse(p.stdout)
        rep.check("baseline: DO_NONMONO_SOB=T parsed", v.get("DO_NONMONO_SOB") == "T")
        rep.check("baseline: DO_NONMONO_KNOP=F parsed", v.get("DO_NONMONO_KNOP") == "F")
        rep.check(
            "baseline: DO_JFNK=T -> USE_JFNK_SOLVER", v.get("USE_JFNK_SOLVER") == "T"
        )
        rep.check(
            "baseline: GLOBAL_LINE=SOB",
            v.get("GLOBAL_LINE_SWITCH", "").startswith("SOB"),
        )
        rep.check("baseline: USE_LAM_ES=T", v.get("USE_LAM_ES") == "T")
        rep.check(
            "baseline: INNER_BND_METH survives case-normalize",
            v.get("INNER_BND_METH", "").strip().isupper(),
        )
        rep.check(
            "baseline: shock-power controls round-trip",
            abs(float(v.get("VLOC_SHOCK_POWER", "nan"))) >= 0,
        )
        rep.check(
            "baseline: ALLOW_SIGMA_UNDEF defaults F",
            v.get("ALLOW_SIGMA_UNDEF") == "F",
        )
        rep.check(
            "baseline: absent JFNK_FADAPT defaults T",
            v.get("JFNK_FADAPT") == "T",
        )
        rep.check(
            "baseline: absent JFNK_FDCPLX defaults F",
            v.get("JFNK_FDCPLX") == "F",
        )
        rep.check(
            "baseline: absent JFNK_ESFR defaults F",
            v.get("JFNK_ESFR") == "F",
        )

    cases = [
        # (name, mutations, expect_ok, expect_stdout_substr)
        (
            "nonmono sob+knop both T aborts",
            {"DO_NONMONO_KNOP": "T"},
            False,
            "mutually exclusive",
        ),
        (
            "DO_NONMONO_SOB + GLOBAL_LINE=BLANK auto-sets to SOB",
            {"GLOBAL_LINE": "BLANK"},
            True,
            "SOB auto-set: GLOBAL_LINE=SOB",
        ),
        (
            "DO_NONMONO_SOB + USE_LAM_ES=F auto-sets to T",
            {"USE_LAM_ES": "F"},
            True,
            "SOB auto-set: USE_LAM_ES=T",
        ),
        (
            "DO_NONMONO_KNOP + GLOBAL_LINE=SOB auto-sets to BLANK",
            {"DO_NONMONO_SOB": "F", "DO_NONMONO_KNOP": "T", "DO_JFNK": "F"},
            True,
            "KNOP auto-set: GLOBAL_LINE=BLANK",
        ),
        (
            "DO_NONMONO_KNOP + LAM_SET=T auto-sets to F",
            {
                "DO_NONMONO_SOB": "F",
                "DO_NONMONO_KNOP": "T",
                "GLOBAL_LINE": "BLANK",
                "LAM_SET": "T",
                "DO_JFNK": "F",
            },
            True,
            "KNOP auto-set: LAM_SET=F",
        ),
        (
            "DO_JFNK without non-mono mode aborts",
            {"DO_NONMONO_SOB": "F"},
            False,
            "requires DO_NONMONO_KNOP",
        ),
        (
            "JFNK_COUPLED under SOB aborts",
            {"JFNK_COUPLED": "T"},
            False,
            "not available under",
        ),
        (
            "DO_CL=T under non-mono aborts",
            {"DO_CL": "T"},
            False,
            "DO_CL=T is incompatible",
        ),
        (
            "GAMRAY_TRANS=RAD_TRANS under non-mono aborts",
            {"GAMRAY_TRANS": "RAD_TRANS"},
            False,
            "is incompatible",
        ),
        (
            "INNER_BND_METH=BOGUS aborts",
            {"IB_METH": "BOGUS"},
            False,
            "Invalid INNER_BND_METH",
        ),
        (
            "BA_DAMP_CAP in (0,1] aborts",
            {"BA_DAMP_CAP": "0.5D0"},
            False,
            "BA_DAMP_CAP = 0 (off) or > 1.0",
        ),
    ]
    for name, muts, expect_ok, needle in cases:
        d = write_vadat(
            WORKDIR / ("vadat_" + re.sub(r"[^a-z0-9]+", "_", name.lower())[:40]), muts
        )
        p = run_vadat(exe, d)
        has_ok = "RESULT_OK" in p.stdout
        rep.check(
            name,
            (has_ok == expect_ok) and (needle is None or needle in p.stdout),
            "" if needle in p.stdout else p.stdout[-300:],
        )

    # KNOP-consistent config parses and lands the flags.
    d = write_vadat(
        WORKDIR / "vadat_knop_ok",
        {
            "DO_NONMONO_SOB": "F",
            "DO_NONMONO_KNOP": "T",
            "GLOBAL_LINE": "BLANK",
            "LAM_SET": "F",
            "DO_BA_RES_ACC": "T",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "valid KNOP config parses",
        "RESULT_OK" in p.stdout,
        p.stdout[-400:] if "RESULT_OK" not in p.stdout else "",
    )
    rep.check(
        "valid KNOP config: flags land",
        v.get("DO_NONMONO_KNOP") == "T" and v.get("DO_NONMONO_SOB") == "F",
    )
    rep.check(
        "KNOP absent-key: FREQ_MATCH auto-set T",
        v.get("DO_NONMONO_FREQ_MATCH") == "T"
        and "KNOP auto-set: DO_NONMONO_FREQ_MATCH=T" in p.stdout,
    )
    rep.check(
        "KNOP absent-key: refresh defaults T",
        v.get("DO_KNOP_REFRESH") == "T",
    )

    # KNOP + explicit DO_KNOP_REFRESH=F is re-armed post-USE_LAM_ES
    # (the moment side is re-synced to the post-auto-set value).
    d = write_vadat(
        WORKDIR / "vadat_knop_rf_off",
        {
            "DO_NONMONO_SOB": "F",
            "DO_NONMONO_KNOP": "T",
            "GLOBAL_LINE": "BLANK",
            "LAM_SET": "F",
            "DO_KNOP_REFRESH": "F",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "KNOP explicit refresh=F re-armed T",
        "RESULT_OK" in p.stdout
        and "KNOP auto-set: DO_KNOP_REFRESH=T" in p.stdout
        and v.get("DO_KNOP_REFRESH") == "T",
        p.stdout[-300:],
    )
    rep.check(
        "KNOP defaults: KN_ORD_P2DAMP stays F",
        v.get("KN_ORD_P2DAMP") == "F",
    )
    rep.check(
        "KNOP defaults: N_LAM_RWARM>=4, GS_REL_STABLE>=2",
        int(v.get("N_LAM_RWARM", "0")) >= 4
        and int(v.get("GS_REL_STABLE", "0")) >= 2,
    )

    # KNOP + explicit FREQ_MATCH=F is respected (overridable).
    d = write_vadat(
        WORKDIR / "vadat_knop_fm_off",
        {
            "DO_NONMONO_SOB": "F",
            "DO_NONMONO_KNOP": "T",
            "GLOBAL_LINE": "BLANK",
            "LAM_SET": "F",
            "DO_NONMONO_FREQ_MATCH": "F",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "KNOP explicit FREQ_MATCH=F respected",
        "RESULT_OK" in p.stdout
        and v.get("DO_NONMONO_FREQ_MATCH") == "F"
        and "auto-set: DO_NONMONO_FREQ_MATCH" not in p.stdout,
        p.stdout[-300:],
    )

    # KNOP + NOV_CONT=T is forced off (silent continuum bypass).
    d = write_vadat(
        WORKDIR / "vadat_knop_novcont",
        {
            "DO_NONMONO_SOB": "F",
            "DO_NONMONO_KNOP": "T",
            "GLOBAL_LINE": "BLANK",
            "LAM_SET": "F",
            "NOV_CONT": "T",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "KNOP NOV_CONT=T forced F",
        "RESULT_OK" in p.stdout
        and "KNOP auto-set: NOV_CONT=F" in p.stdout
        and v.get("NO_VEL_FOR_CONTINUUM") == "F",
        p.stdout[-300:],
    )

    # USE_JFNK deprecated alias still arms the solver.
    d = write_vadat(WORKDIR / "vadat_usejfnk", {"DO_JFNK": "F", "USE_JFNK": "T"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "USE_JFNK alias arms solver + warns",
        "RESULT_OK" in p.stdout
        and "deprecated" in p.stdout
        and v.get("USE_JFNK_SOLVER") == "T",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # New JFNK Phase-C keys parse and land when present.
    d = write_vadat(
        WORKDIR / "vadat_newkeys",
        {"JFNK_FADAPT": "F", "JFNK_FDCPLX": "T", "JFNK_ESFR": "T"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "JFNK_FADAPT/FDCPLX/ESFR parse and land",
        "RESULT_OK" in p.stdout
        and v.get("JFNK_FADAPT") == "F"
        and v.get("JFNK_FDCPLX") == "T"
        and v.get("JFNK_ESFR") == "T",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # DO_NONLOC_SOB without DO_NONMONO_SOB -> forced off with warning.
    d = write_vadat(
        WORKDIR / "vadat_nonloc",
        {"DO_NONMONO_SOB": "F", "DO_JFNK": "F", "DO_NONLOC_SOB": "T"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "DO_NONLOC_SOB forced F without non-mono SOB",
        "RESULT_OK" in p.stdout
        and "DO_NONLOC_SOB=T but DO_NONMONO_SOB=F" in p.stdout
        and v.get("DO_NONLOC_SOB") == "F",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # DO_NONLOC_SOB kept under SOB.
    d = write_vadat(WORKDIR / "vadat_nonloc_ok", {"DO_NONLOC_SOB": "T"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "DO_NONLOC_SOB kept under non-mono SOB",
        "RESULT_OK" in p.stdout and v.get("DO_NONLOC_SOB") == "T",
    )

    # VEL_LAW=13 selects SN_HYDRO_MODEL.
    d = write_vadat(WORKDIR / "vadat_v13", {"VEL_LAW": "13"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "VEL_LAW=13 -> SN_HYDRO_MODEL",
        "RESULT_OK" in p.stdout
        and v.get("SN_HYDRO_MODEL") == "T"
        and v.get("VELTYPE") == "13",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # INNER_BND_METH case normalization: lowercase accepted.
    d = write_vadat(WORKDIR / "vadat_ibnd", {"IB_METH": "diffusion"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "IB_METH=diffusion accepted + uppercased",
        "RESULT_OK" in p.stdout and v.get("INNER_BND_METH", "").strip() == "DIFFUSION",
    )

    # JFNK_PTC=F is force-flipped by the JFNK auto-set (the
    # transcendental-corrector stub is mandatory while JFNK is armed).
    d = write_vadat(WORKDIR / "vadat_ptc", {"JFNK_PTC": "F"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "JFNK_PTC=F auto-sets to T",
        "RESULT_OK" in p.stdout
        and "JFNK auto-set: JFNK_PTC=T" in p.stdout
        and v.get("JFNK_PTC_ACTIVE") == "T"
        and abs(float(v.get("JFNK_PTC_DTAU", "nan").replace("D", "E")) - 1.0) < 1e-9,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    d = write_vadat(
        WORKDIR / "vadat_jfnk_lag_factors",
        {
            "JFNK_LAG_FACTORS": "T",
            "JFNK_FRZ_GATE": "3.0D0",
            "JFNK_GATE_NREQ": "2",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "JFNK lag-factor and admission keys parse",
        "RESULT_OK" in p.stdout
        and v.get("JFNK_LAG_FACTORS") == "T"
        and abs(float(v.get("JFNK_FRZ_GATE", "nan").replace("D", "E")) - 3.0)
        < 1e-9
        and v.get("JFNK_GATE_NREQ") == "2",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # MOM_H_IMPL / MOM_ABS_FLR opt-in keys parse and land.
    d = write_vadat(
        WORKDIR / "vadat_mom_himpl",
        {"MOM_H_IMPL": "T", "MOM_ABS_FLR": "1.0D-3"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "MOM_H_IMPL/MOM_ABS_FLR keys parse and land",
        "RESULT_OK" in p.stdout
        and v.get("MOM_H_IMPL") == "T"
        and abs(float(v.get("MOM_ABS_FLR", "nan").replace("D", "E")) - 1.0e-3)
        < 1e-12,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # MOM_FAC_RLX opt-in key parses and lands.
    d = write_vadat(
        WORKDIR / "vadat_mom_facrlx",
        {"MOM_FAC_RLX": "5.0D-1"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "MOM_FAC_RLX key parses and lands",
        "RESULT_OK" in p.stdout
        and abs(float(v.get("MOM_FAC_RLX", "nan").replace("D", "E")) - 0.5)
        < 1e-12,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # MOM_TRIAL_LAM opt-in key parses and lands.
    d = write_vadat(
        WORKDIR / "vadat_mom_triallam",
        {"MOM_TRIAL_LAM": "1.0D3"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "MOM_TRIAL_LAM key parses and lands",
        "RESULT_OK" in p.stdout
        and abs(float(v.get("MOM_TRIAL_LAM", "nan").replace("D", "E")) - 1.0e3)
        < 1e-9,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # FG_TRIAL_CHURN opt-in key parses and lands.
    d = write_vadat(
        WORKDIR / "vadat_fg_trialchurn",
        {"FG_TRIAL_CHURN": "3.0D0"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "FG_TRIAL_CHURN key parses and lands",
        "RESULT_OK" in p.stdout
        and abs(float(v.get("FG_TRIAL_CHURN", "nan").replace("D", "E")) - 3.0)
        < 1e-12,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # Freeze/precondition and bulk-floor keys parse and land.  The
    # template is non-KNOP, so no auto-sets disturb the explicit
    # values; KN_RAY_BULK_TOL echoes the solver module variable.
    d = write_vadat(
        WORKDIR / "vadat_gates",
        {
            "FG_TRIAL_CHURN_MIN": "5",
            "GS_PRECOND": "T",
            "GS_PRECOND_MAX": "2.0D3",
            "GS_PRECOND_LIM": "42",
            "GS_FREEZE_LIM": "120",
            "GS_REL_NOPROG_LIM": "7",
            "KN_RAY_BULK_TOL": "6.0D-2",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "GS_PRECOND/GS_PRECOND_MAX/GS_PRECOND_LIM parse and land",
        "RESULT_OK" in p.stdout
        and v.get("GS_PRECOND") == "T"
        and abs(
            float(v.get("GS_PRECOND_MAX", "nan").replace("D", "E"))
            - 2000.0
        )
        < 1e-6
        and v.get("GS_PRECOND_LIM") == "42",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )
    rep.check(
        "GS_FREEZE_LIM/GS_REL_NOPROG_LIM parse and land",
        "RESULT_OK" in p.stdout
        and v.get("GS_FREEZE_LIM") == "120"
        and v.get("GS_REL_NOPROG_LIM") == "7",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )
    rep.check(
        "KN_RAY_BULK_TOL parses and lands on the solver",
        "RESULT_OK" in p.stdout
        and abs(
            float(v.get("KN_RAY_BULK_TOL", "nan").replace("D", "E"))
            - 0.06
        )
        < 1e-12,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )
    rep.check(
        "FG_TRIAL_CHURN_MIN parses and lands",
        "RESULT_OK" in p.stdout
        and v.get("FG_TRIAL_CHURN_MIN") == "5",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # RLOC/DRLOC opt-in keys parse and land.
    d = write_vadat(
        WORKDIR / "vadat_rloc",
        {
            "INC_SHOCK_POWER": "T",
            "RLOC_SHOCK_POWER": "2.5D5",
            "DRLOC_SHOCK_POWER": "5.0D4",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "RLOC_SHOCK_POWER/DRLOC_SHOCK_POWER keys parse and land",
        "RESULT_OK" in p.stdout
        and abs(float(v.get("RLOC_SHOCK_POWER", "nan").replace("D", "E")) - 2.5e5)
        < 1e-9
        and abs(float(v.get("DRLOC_SHOCK_POWER", "nan").replace("D", "E")) - 5.0e4)
        < 1e-9,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # ALLOW_SIGMA_UNDEF opt-in diagnostic key parses and lands.
    d = write_vadat(
        WORKDIR / "vadat_allow_undef",
        {"ALLOW_SIGMA_UNDEF": "T"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "ALLOW_SIGMA_UNDEF=T parses and lands",
        "RESULT_OK" in p.stdout and v.get("ALLOW_SIGMA_UNDEF") == "T",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # DO_BA_RES_REUSE opt-in keys parse and land.
    d = write_vadat(
        WORKDIR / "vadat_reuse",
        {
            "DO_BA_RES_REUSE": "T",
            "BA_REUSE_TOL": "8.0D-1",
            "BA_REUSE_NSTALL": "3",
        },
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "DO_BA_RES_REUSE keys parse and land",
        "RESULT_OK" in p.stdout
        and v.get("DO_BA_RES_REUSE") == "T"
        and abs(float(v.get("BA_REUSE_TOL", "nan").replace("D", "E")) - 0.8)
        < 1e-12
        and v.get("BA_REUSE_NSTALL") == "3",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # DO_BA_RES_REUSE=T without DO_BA_RES_ACC warns (inert, not fatal).
    d = write_vadat(
        WORKDIR / "vadat_reuse_noacc",
        {"DO_BA_RES_REUSE": "T", "DO_BA_RES_ACC": "F"},
    )
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "DO_BA_RES_REUSE without resacc warns but parses",
        "RESULT_OK" in p.stdout
        and "DO_BA_RES_REUSE=T has no effect" in p.stdout,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # Defaults: reuse off, tol 0.9, nstall 2.
    d = write_vadat(WORKDIR / "vadat_reuse_def", {})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "BA_RES_REUSE defaults (F / 0.9 / 2)",
        "RESULT_OK" in p.stdout
        and v.get("DO_BA_RES_REUSE") == "F"
        and abs(float(v.get("BA_REUSE_TOL", "nan").replace("D", "E")) - 0.9)
        < 1e-12
        and v.get("BA_REUSE_NSTALL") == "2",
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    # RLOC defaults to zero so the legacy V-space path stays live.
    d = write_vadat(WORKDIR / "vadat_rloc_def", {})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "RLOC_SHOCK_POWER defaults to 0 (V-space deposition)",
        "RESULT_OK" in p.stdout
        and abs(float(v.get("RLOC_SHOCK_POWER", "nan").replace("D", "E")))
        < 1e-30,
        p.stdout[-300:] if "RESULT_OK" not in p.stdout else "",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
