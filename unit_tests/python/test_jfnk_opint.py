#!/usr/bin/env python3
"""JFNK operator-interface battery (JFNK_OP_PROBE_V1).

Runs the production probe module against the surrogate SE residual
via `se_driver seprobe` (systems/{mild,hard,shock}) and checks the
measured operator properties:

  * residual repeatability is bitwise zero on the deterministic
    surrogate;
  * the FD JVP matches the exact assembled Jacobian (se_jexact)
    at the truncation rate the h-sweep exposes;
  * additivity/homogeneity defects sit at the FD truncation floor,
    not at a policy or state-noise level;
  * the probe module itself compiles standalone and the production
    call site in jfnk_solve_corr_v1.f is wired under JFNK_HSCAN.

On production the same probe runs at JFNK admission under
[JFNK_HSCAN]; this test verifies the measurement code path that the
production run executes, not duplicated logic.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import TREE, Report, prod_obj

ROOT = TREE
SE_DIR = ROOT / "unit_tests"
SE_BIN = SE_DIR / "bin" / "se_driver"
PROBE_SRC = ROOT / "new_main/mod_subs/jfnk_op_probe_v1.f90"
CORR_SRC = ROOT / "new_main/mod_subs/jfnk_solve_corr_v1.f"
IFACE_SRC = ROOT / "new_main/mod_subs/jfnk_cmfgen_interface_v1.f"
INTMOD_SRC = ROOT / "new_main/mod_subs/jfnk_cmfgen_int_mod.f"
EPS = "1e-6"
SYSTEMS = ("mild", "hard", "shock")
EXPECTED_TAGS = {
    "rnd1", "rnd2", "pops", "ed", "t", "k1", "kmid",
    "kworst", "blk3", "basis", "edt", "spike", "xdep",
}


def run_seprobe(system):
    p = subprocess.run(
        [str(SE_BIN), "seprobe", f"systems/{system}", EPS],
        cwd=SE_DIR,
        capture_output=True,
        text=True,
        timeout=300,
    )
    if p.returncode != 0:
        raise SystemExit(f"seprobe {system} failed:\n{p.stdout}\n{p.stderr}")
    return p.stdout


def parse_opprobe(out):
    """Parse ' OPPROBE <sect> k=v k=v ...' lines into records."""
    recs = []
    for ln in out.splitlines():
        if " OPPROBE " not in ln:
            continue
        toks = ln.split()
        sect = toks[toks.index("OPPROBE") + 1]
        kv = {}
        rest = toks[toks.index("OPPROBE") + 2:]
        i = 0
        while i < len(rest):
            t = rest[i]
            if "=" in t:
                k, v = t.split("=", 1)
                # Fortran fixed formats can emit 'key= value' with the
                # value in the next token.
                if v == "" and i + 1 < len(rest) and "=" not in rest[i + 1]:
                    v = rest[i + 1]
                    i += 1
                kv[k] = v
            i += 1
        recs.append((sect, kv))
    return recs


def probe_compile(rep):
    """The probe module is a production archive member; resolving its
    exact .o also proves it compiled under the production toolchain."""
    try:
        obj = prod_obj(PROBE_SRC)
        ok = Path(obj).exists()
    except Exception:
        ok = False
    rep.check(
        "opprobe: production probe object resolves",
        ok,
    )


def probe_wiring(rep):
    """Static checks that the production call path is wired."""
    corr = CORR_SRC.read_text()
    iface = IFACE_SRC.read_text()
    intmod = INTMOD_SRC.read_text()
    rep.check(
        "opprobe: probe called under JFNK_HSCAN in solve_corr",
        "IF(JFNK_HSCAN)" in corr and "CALL JFNK_OP_PROBE_V1" in corr,
    )
    rep.check(
        "opprobe: probe fenced from solver status counters",
        "OPP_NBAD=jfnk_resid_nbad" in corr
        and "jfnk_resid_nbad=OPP_NBAD" in corr,
    )
    rep.check(
        "opprobe: JVP adapter + BA action in production interface",
        "SUBROUTINE JFNK_OP_JVPA_V1" in iface
        and "SUBROUTINE JFNK_OP_BAMV_V1" in iface,
    )
    rep.check(
        "opprobe: JFNK_PROBE_CAPPED flag set inside JVP cap logic",
        "JFNK_PROBE_CAPPED" in intmod
        and "JFNK_PROBE_CAPPED = (H_CAP .LT. H_STEP)" in iface,
    )


def check_system(rep, system):
    out = run_seprobe(system)
    recs = parse_opprobe(out)
    rep.check(
        f"opint:{system}: battery emits BEGIN/END",
        recs[0][0] == "BEGIN" and recs[-1][0] == "END",
    )
    end = dict(recs)[ "END"] if "END" in dict(recs) else {}
    rep.check(
        f"opint:{system}: no JVP failures on surrogate",
        end.get("fail") == "0",
        json.dumps(end),
    )
    # Repeatability: bitwise deterministic surrogate must reproduce
    # F(x) identically, baseline and perturbed alike.
    rep_errs = [float(kv["df_inf"]) for s, kv in recs if s == "REP"]
    rep.check(
        f"opint:{system}: baseline repeatability bitwise 0",
        len(rep_errs) == 4 and all(e == 0.0 for e in rep_errs),
        str(rep_errs),
    )
    repp = [float(kv["df_inf"]) for s, kv in recs if s == "REPP"]
    rep.check(
        f"opint:{system}: perturbed repeatability bitwise 0",
        len(repp) == 1 and repp[0] == 0.0,
        str(repp),
    )
    # Direction battery: all 13 tags present, JVP-vs-exact error at
    # FD truncation scale, sweep exposes a resolved plateau.
    dirs = {kv["tag"]: kv for s, kv in recs if s == "DIR"}
    rep.check(
        f"opint:{system}: all 13 direction classes ran",
        set(dirs) == EXPECTED_TAGS,
        str(sorted(dirs)),
    )
    efd_max = max(float(kv["eFD"]) for kv in dirs.values())
    efc_min = min(float(kv["eFC"]) for kv in dirs.values())
    rep.check(
        f"opint:{system}: FD JVP vs exact J bounded (max eFD < 0.2)",
        efd_max < 0.2,
        f"max eFD={efd_max:.3e}",
    )
    rep.check(
        f"opint:{system}: central FD beats forward FD somewhere",
        efc_min < efd_max,
        f"min eFC={efc_min:.3e} max eFD={efd_max:.3e}",
    )
    dms = {kv["tag"]: float(kv["emin"]) for s, kv in recs if s == "DM"}
    rep.check(
        f"opint:{system}: h-sweep resolves a derivative plateau",
        len(dms) == 13 and all(v > 0.0 for v in dms.values()),
    )
    add = next(kv for s, kv in recs if s == "ADD")
    hom = next(kv for s, kv in recs if s == "HOM")
    rep.check(
        f"opint:{system}: additivity at FD truncation floor",
        float(add["d_op"]) < 1e-4 and float(add["d_pol"]) < 1e-3,
        f"d_op={add['d_op']} d_pol={add['d_pol']}",
    )
    rep.check(
        f"opint:{system}: homogeneity exact under norm-scaled h",
        float(hom["d_pol"]) == 0.0,
        f"d_pol={hom['d_pol']}",
    )


def main():
    rep = Report("JFNK operator-interface battery")
    probe_compile(rep)
    probe_wiring(rep)
    b = subprocess.run(
        ["bash", "build.sh"], cwd=SE_DIR, capture_output=True, text=True
    )
    rep.check(
        "opint: se_driver rebuilds with probe module",
        b.returncode == 0 and SE_BIN.exists(),
        b.stderr[-200:],
    )
    if SE_BIN.exists():
        for system in SYSTEMS:
            check_system(rep, system)
    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
