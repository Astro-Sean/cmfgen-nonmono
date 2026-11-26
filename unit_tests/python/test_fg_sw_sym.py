#!/usr/bin/env python3
"""test_fg_sw_sym -- symmetric FG sweep drift + bad-cell identity.

Round-2 phase-2 patch C: the in-sweep store drift metric was
one-sided (|fresh-stored| / max(stored, 1e-4*peak)), which reads
~1eN on a cell whose stored value was floor-sized -- a magnitude
blowup indistinguishable from genuine field churn. FG_GS_MAXSYM /
FG_SW_MAXSYM add the symmetric companion
|fresh-stored| / max(|fresh|,|stored|,1e-4*peak) with its own
argmax (depth, ray, frequency, P/M channel), snapshotted and
exported beside the one-sided metric. KN_BADID adds the cross-
refresh bad-cell identity census (sparse packed-id lists, serial
merge, Jaccard diff).

Driver checks run the production module procedures directly
(drv_badid) and the real ray solve end-to-end (drv_knop_ray cases
badid/bido, which arm the flag and verify the merged census equals
the solve's own bad-cell count). Source pins lock the census
sites, lifecycle (reset/snapshot/push/pop), and gating.

Run: python3 unit_tests/python/test_fg_sw_sym.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common

REP = common.Report("test_fg_sw_sym (symmetric drift + bad-cell id)")

WORK = common.WORKDIR / "fg_sw_sym"
WORK.mkdir(parents=True, exist_ok=True)

FG = (
    common.SRC / "plane" / "fg_j_cmf_nonmono_v1.f"
).read_text()
SFP = (
    common.SRC / "new_main" / "mod_subs" / "solve_for_pops.f"
).read_text()
CVM = (
    common.SRC / "new_main" / "mod_subs" / "control_variable_mod.f"
).read_text()
RDV = (
    common.SRC / "new_main" / "mod_subs" / "rd_control_variables.f"
).read_text()


def count(pat, text=FG):
    return len(re.findall(pat, text))


# ----------------------------------------------------------------------
# Source pins: symmetric metric.
# ----------------------------------------------------------------------

REP.check(
    "symmetric formula uses max(|fresh|,|stored|,floor)",
    "dval=abs(fresh-stored)/" in FG
    and "max(abs(fresh),abs(stored),floor)" in FG,
)
REP.check(
    "symmetric accumulator keeps argmax",
    "fg_gs_maxsym=dval" in FG
    and "fg_gs_symk=k" in FG and "fg_gs_symc=idir" in FG,
)
REP.check(
    "P channel accumulates with VMAX_IP, channel id 1",
    count(r"VMAX_IP\(LS\)\),K,LS,FREQ_IDX_ITER,1\)") == 2,
)
REP.check(
    "M channel accumulates with VMAX_IM, channel id 2",
    count(r"VMAX_IM\(LS\)\),K,LS,FREQ_IDX_ITER,2\)") == 2,
)
REP.check(
    "one-sided T4 metric still computed and gated",
    count(r"IF\(T4 \.GT\. FG_GS_MAXCHG\)THEN") == 2,
)
REP.check(
    "sweep snapshot moves both metrics",
    "FG_SW_MAXCHG=FG_GS_MAXCHG" in FG
    and "FG_SW_MAXSYM=FG_GS_MAXSYM" in FG,
)
REP.check(
    "snapshot copies the argmax cell",
    "FG_SW_SYMK=FG_GS_SYMK" in FG
    and "FG_SW_SYMR=FG_GS_SYMR" in FG
    and "FG_SW_SYML=FG_GS_SYML" in FG
    and "FG_SW_SYMC=FG_GS_SYMC" in FG,
)
REP.check(
    "INIT resets live and snapshot state",
    "FG_GS_MAXSYM=0.0_LDP" in FG
    and "FG_GS_SYMK=0" in FG
    and "FG_SW_MAXSYM=0.0_LDP" in FG
    and "FG_SW_SYMK=0" in FG,
)
REP.check(
    "eval push/pop saves and restores the sym state",
    "FG_EVAL_GSYM=FG_GS_MAXSYM" in FG
    and "FG_GS_MAXSYM=FG_EVAL_GSYM" in FG
    and "FG_EVAL_SYK=FG_GS_SYMK" in FG
    and "FG_GS_SYMK=FG_EVAL_SYK" in FG,
)
REP.check(
    "epilogue prints both maxima on one line",
    "fg_sw_one=" in FG and "fg_sw_sym=" in FG,
)
REP.check(
    "argmax reported under a positive-ray guard",
    "fg_sw_sym argmax" in FG and "FG_SW_SYMR .GT. 0" in FG,
)
REP.check(
    "GET_FG_SW_SYM exports snapshot + argmax",
    "SUBROUTINE GET_FG_SW_SYM(SYM,ARGK,ARGR,ARGL,ARGC)" in FG,
)
# The symmetric max must never feed a *_TOL comparison.
REP.check(
    "sym metric never compared to a tolerance",
    not re.search(r"FG_(GS|SW)_MAXSYM\s*\.[GL][ET]\.", FG)
    and not re.search(r"\.[GL][ET]\.\s*FG_(GS|SW)_MAXSYM", FG)
    and not re.search(
        r"(MATURE|UNSAFE|SETTLE)_TOL.{0,80}MAXSYM", FG, re.S
    ),
)

# ----------------------------------------------------------------------
# Source pins: solve_for_pops integration.
# ----------------------------------------------------------------------

REP.check(
    "solve_for_pops reads the sym metric once per maturity call",
    "CALL GET_FG_SW_SYM(GS_SW_MAXSYM,GS_SYK,GS_SYR,GS_SYL,GS_SYC)"
    in SFP,
)
REP.check(
    "fg summary row prints sw= and sym=",
    "' sw=',GS_SW_MAXCHG,' sym=',GS_SW_MAXSYM" in SFP,
)
REP.check(
    "GS_DRIFT ray line carries sym + sym argmax",
    "ray sw/res/rf/sym" in SFP
    and "ray sym arg(k,ray,frq,ch)" in SFP,
)

# ----------------------------------------------------------------------
# Source pins: bad-cell identity tracking.
# ----------------------------------------------------------------------

REP.check(
    "BID_PUSH packs fixed-stride ids",
    "(LS-1_8)*17179869184_8" in FG and "(LF-1_8)*8192_8" in FG
    and "(IPT-1_8)*2_8" in FG,
)
REP.check(
    "append count: MARCH x3 + ORD x4 census sites",
    count(r"CALL BID_PUSH") == 7,
)
REP.check(
    "identity lists reset at the pass tops",
    count(r"KN_BADRAY\(LS,MERGE\(1,2,INWARD\)\)%N=0_8") == 3,
)
REP.check(
    "verify passes never record (L_QP guard at every site)",
    count(r"IF\(\.NOT\. L_QP\)CALL BID_PUSH") == 7,
)
REP.check(
    "serial merge runs only over the production mask",
    "CALL KN_BADID_MERGE_V1(NBD_LS,KN_PROD_MASK,KN_NP,LUER)" in FG,
)
REP.check(
    "merge reports cur/prev/inter/new/gone/jac on one line",
    "bad-cell id: cur=" in FG and "jac=" in FG
    and "bad-cell ch: Mcur=" in FG,
)
REP.check(
    "cap + sentinel bound any one channel list",
    "KN_BADID_CAP" in FG and "%N=-1_8" in FG,
)
REP.check(
    "untracked rays counted, not silently dropped",
    "KN_BID_NUNTR" in FG and "untracked rays=" in FG,
)
REP.check(
    "flag off -> lists deallocated, state reset",
    "KN_BADID_ON=DO_B" in FG and "KN_BID_NP=-1_8" in FG,
)
REP.check(
    "VADAT key wired: declared, parsed, pushed",
    "LOGICAL KN_BADID" in CVM
    and "CALL RD_STORE_LOG(KN_BADID,'KN_BADID'" in RDV
    and "CALL SET_FG_BADID(KN_BADID)" in RDV,
)
REP.check(
    "criterion untouched: NBAD_CUR test is still DCELL>KN_RAY_TOL",
    "IF(DCELL .GT. KN_RAY_TOL)THEN" in FG
    or "IF(T1 .GT. KN_RAY_TOL)THEN" in FG,
)

# ----------------------------------------------------------------------
# Driver: analytic symmetric values + accumulator + sort + merge.
# The FG module's external references need the same helper set as
# the knray driver.
# ----------------------------------------------------------------------

HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/indexx.f",
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

binpath = common.compile_and_link_prod(
    "drv_badid",
    [common.SRC / h for h in HELPERS] + [common.F90DIR / "stub_vadat.f", common.F90DIR / "drv_badid.f90"],
)
r = common.run(binpath, cwd=WORK)
REP.check("drv_badid ran", r.returncode == 0 and "RESULT_OK" in r.stdout)

vals = {}
for ln in r.stdout.splitlines():
    m = re.match(r"^\s*([A-Za-z0-9_]+)=(.+)$", ln.strip())
    if m:
        vals[m.group(1)] = m.group(2).strip()


def f(name, default=None):
    try:
        return float(vals.get(name, "nan"))
    except ValueError:
        return default


def i(name, default=None):
    try:
        return int(vals.get(name, "-999"))
    except ValueError:
        return default


# Analytic: dsym(2,1,1e-4)=0.5; dsym(1e10,1,f)~1; floor dominance;
# sign flip can exceed 1; equal values -> 0.
REP.check("dsym(2,1,floor)=0.5", abs(f("DSYM_A") - 0.5) < 1e-12)
REP.check(
    "dsym(1e10,1,floor)~1 (bounded, not 1e10)",
    0.999 < f("DSYM_B") < 1.001,
)
REP.check(
    "dsym small when both near floor",
    f("DSYM_C") < 0.02,
)
REP.check(
    "dsym sign-flip >1 (one-sided pin would hide it)",
    f("DSYM_D") > 1.0,
)
REP.check("dsym(0,0)=0", f("DSYM_E") == 0.0)
REP.check(
    "dsym floor-bound pair < 1",
    0.0 < f("DSYM_F") < 1.0,
)

# Accumulator: max wins, argmax carries (k,ray,freq,channel).
REP.check("accumulator picked the 0.7 cell", abs(f("ACC_SYM") - 0.7) < 1e-12)
m_arg = re.search(
    r"ACC_ARG=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", r.stdout
)
REP.check(
    "accumulator argmax = (9,3,200,2)",
    m_arg is not None
    and (int(m_arg.group(1)), int(m_arg.group(2)),
         int(m_arg.group(3)), int(m_arg.group(4))) == (9, 3, 200, 2),
)
m_get = re.search(
    r"GET_ARG=\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", r.stdout
)
REP.check(
    "GET_FG_SW_SYM echoes snapshot state",
    abs(f("GET_SYM") - 0.42) < 1e-12
    and m_get is not None
    and tuple(int(g) for g in m_get.groups()) == (4, 7, 88, 1),
)
REP.check("int8 heapsort ordered", vals.get("SORTED") == "T")

# Merge: first census has no previous -> jac=-1, cur=5; the merge
# stores the new census for the next diff so NP reads 5 after.
REP.check("first census cur=5 stored as next prev",
          i("M1_NC") == 5 and i("M1_NP") == 5)
REP.check("first census all new, none gone",
          i("M1_NEW") == 5 and i("M1_GONE") == 0)
REP.check("first census jac=-1 sentinel", f("M1_JAC") < 0.0)
REP.check("channel split M=4 P=1", i("M1_MCUR") == 4 and i("M1_PCUR") == 1)
REP.check("no untracked rays", i("M1_UNTR") == 0)

# Identical second census: Jaccard = 1.
REP.check("identical census: inter=5 jac=1",
          i("M2_INT") == 5 and abs(f("M2_JAC") - 1.0) < 1e-12)
REP.check("identical census: no churn",
          i("M2_NEW") == 0 and i("M2_GONE") == 0)
REP.check("channel inter M=4 P=1",
          i("M2_MINT") == 4 and i("M2_PINT") == 1)

# Partial overlap: ray2 cell moved -> inter=4, new=1, gone=1,
# union=6, jac=4/6.
REP.check("moved cell: inter=4 new=1 gone=1",
          i("M3_INT") == 4 and i("M3_NEW") == 1 and i("M3_GONE") == 1)
REP.check("moved cell: jac=2/3", abs(f("M3_JAC") - 4.0 / 6.0) < 1e-6)

# Untracked ray: sentinel channel excludes the ray, keeps counts.
REP.check("untracked ray excluded: cur=4", i("M4_NC") == 4)
REP.check("untracked counted", i("M4_UNTR") == 1)

# ----------------------------------------------------------------------
# Driver: end-to-end through KNOP_RAY_SOLVE_V1 (cap config, tight
# tolerance -> nonzero last-pass bad cells). The merged census must
# equal the solve's own NBAD exactly -- every above-tolerance cell
# was recorded once. badid runs the unordered MARCH path, bido the
# ordered SWEEP_FLD_ORD path.
# ----------------------------------------------------------------------

binray = common.compile_and_link_prod(
    "knray_badid",
    [common.SRC / h for h in HELPERS] + [common.F90DIR / "stub_vadat.f", common.F90DIR / "drv_knop_ray.f90"],
)

for case in ("badid", "bido"):
    r2 = common.run(binray, case, cwd=WORK, check=False)
    d = {}
    for ln in r2.stdout.splitlines():
        m = re.match(r"^(\w+)=\s*(\S+)", ln)
        if m:
            try:
                d[m.group(1)] = float(m.group(2))
            except ValueError:
                pass
    REP.check(
        f"{case}: refresh ran with a nonempty last-pass census",
        "RESULT_OK" in r2.stdout and d.get("RF_NBAD", -1) > 0,
    )
    REP.check(
        f"{case}: identity list total equals FG_GS_NBAD exactly",
        int(d.get("BID_SUM", -1)) == int(d.get("RF_NBAD", -2)),
    )
    REP.check(
        f"{case}: merged census cur equals FG_GS_NBAD",
        int(d.get("BID_NC", -1)) == int(d.get("RF_NBAD", -2)),
    )
    REP.check(
        f"{case}: stored prev equals cur after first merge",
        int(d.get("BID_NP", -1)) == int(d.get("RF_NBAD", -2)),
    )
    REP.check(
        f"{case}: first census jac=-1, untracked=0",
        d.get("BID_JAC", 0) < 0 and int(d.get("BID_UNTR", 1)) == 0,
    )


def main():
    return REP.summary()


if __name__ == "__main__":
    sys.exit(0 if main() == 0 else 1)
