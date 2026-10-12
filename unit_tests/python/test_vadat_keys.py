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

VADAT_SRC = Path(
    "/viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/"
    "m15_2p82d_SOB_jfnk_zlaw14_fresh/cmfgen/2.82days/VADAT"
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

    cases = [
        # (name, mutations, expect_ok, expect_stdout_substr)
        (
            "nonmono sob+knop both T aborts",
            {"DO_NONMONO_KNOP": "T"},
            False,
            "mutually exclusive",
        ),
        (
            "DO_NONMONO_SOB + GLOBAL_LINE=BLANK aborts",
            {"GLOBAL_LINE": "BLANK"},
            False,
            "requires GLOBAL_LINE=SOB",
        ),
        (
            "DO_NONMONO_SOB + USE_LAM_ES=F warns but parses",
            {"USE_LAM_ES": "F"},
            True,
            "USE_LAM_ES=F",
        ),
        (
            "DO_NONMONO_KNOP + GLOBAL_LINE=SOB aborts",
            {"DO_NONMONO_SOB": "F", "DO_NONMONO_KNOP": "T", "DO_JFNK": "F"},
            False,
            "requires GLOBAL_LINE=BLANK",
        ),
        (
            "DO_NONMONO_KNOP + LAM_SET=T aborts",
            {
                "DO_NONMONO_SOB": "F",
                "DO_NONMONO_KNOP": "T",
                "GLOBAL_LINE": "BLANK",
                "LAM_SET": "T",
                "DO_JFNK": "F",
            },
            False,
            "requires LAM_SET=F",
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

    # JFNK_PTC key lands in the solver module (dtau has no VADAT key:
    # it is fixed at 1.0 by the reader).
    d = write_vadat(WORKDIR / "vadat_ptc", {"JFNK_PTC": "F"})
    p = run_vadat(exe, d)
    v = parse(p.stdout)
    rep.check(
        "JFNK_PTC key lands in solver module",
        "RESULT_OK" in p.stdout
        and v.get("JFNK_PTC_ACTIVE") == "F"
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

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
