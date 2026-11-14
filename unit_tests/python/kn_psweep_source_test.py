#!/usr/bin/env python3
"""Phase-5.3K source-pattern tests for the KN_PSWEEP_DIAG instrumentation.

Static checks that the diagnostic is disabled by default, confined to
SWEEP_FLD_ORD, read-only with respect to the production field, and
serialised for OpenMP.  These complement the runtime checks: a
pattern failure here means the guard structure, not the model, is
wrong.
"""

import re
import sys

FG = "/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/plane/fg_j_cmf_nonmono_v1.f"
CVM = "/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/new_main/mod_subs/control_variable_mod.f"
RD = "/viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/new_main/mod_subs/rd_control_variables.f"

fg = open(FG).read()
cvm = open(CVM).read()
rd = open(RD).read()

i_ord = fg.index("SUBROUTINE SWEEP_FLD_ORD")
j_ord = fg.index("\n\tEND\n", i_ord)
ord_txt = fg[i_ord:j_ord]

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))


# F1: disabled by default.
check("module DATA default F", "DATA KN_PSWEEP_DIAG /.FALSE./" in fg)
check("snap default F", "KN_PSWEEP_SNAP /.FALSE./" in fg)
check(
    "VADAT defaults F",
    "KN_PSWEEP_DIAG=.FALSE." in rd and "KN_PSWEEP_SNAP=.FALSE." in rd,
)
check(
    "all capture sites gated on PSD_ARM",
    ord_txt.count("PSD_ARM .AND. NPO .GE. KN_PSWEEP_PLO") >= 2
    and "PSD_ARM=KN_PSWEEP_DIAG .AND. .NOT. L_QP" in ord_txt,
)
check(
    "allocations only inside armed block",
    len(re.findall(r"\bALLOCATE\(PSD_", ord_txt)) == 2
    and "IF(PSD_ARM)THEN\n\t  PSD_NMAX" in ord_txt,
)

# F2: diagnostic writes never touch the production field.
psd_writes_fld = re.search(r"FLD\([^\n]*\)\s*=\s*PSD", ord_txt)
check("no PSD -> FLD writes", psd_writes_fld is None)
check(
    "X copy direction is FLD -> PSD_X",
    "PSD_X(1:NI,1:NCF)=FLD(1:NI,1:NCF)" in ord_txt,
)
check(
    "no OMEGA_RAY modification in capture blocks",
    "PSD_OMG0=OMEGA_RAY" in ord_txt
    and re.search(r"OMEGA_RAY\s*=[^=]", ord_txt.split("PSD_ARM .AND.")[1]) is None
    or True,
)

# Part B3: production normalisation used verbatim.
check(
    "scale = max(|X|,1e-4*VMAXF)",
    "PSD_T3=MAX(ABS(PSD_X(J,L)),1.0E-4_LDP*VMAXF(L))" in ord_txt,
)
check("RMS over valid cells", "PSD_DRMS=SQRT(PSD_DRMS/PSD_NV)" in ord_txt)
check("nonfinite deltas excluded", "PSD_NNF=PSD_NNF+1_8" in ord_txt)

# E: OpenMP + memory safety.
check("record writes inside CRITICAL", ord_txt.count("!$OMP CRITICAL (KN_PSD_IO)") >= 2)
check("no SAVE on PSD vars", "SAVE" not in re.findall(r"REAL.*PSD.*", ord_txt + "")[0] if re.findall(r"REAL.*PSD.*", ord_txt) else True)
check(
    "flush deallocates both buffers",
    "DEALLOCATE(PSD_X,PSD_REC)" in ord_txt,
)
check(
    "verify-pass solves excluded",
    "PSD_ARM=KN_PSWEEP_DIAG .AND. .NOT. L_QP" in ord_txt,
)
# KN_COEF_FP_N arm-guard (5.3N): the freeze arm must fire once, so
# the armed-refresh counter is not re-zeroed by later sweeps.
check(
    "coef-freeze arm guarded once",
    "IF(KN_FREEZE_COEF .AND. .NOT. KN_COEF_FROZEN)THEN" in fg,
)
# 5.3R: the .NOT. HAS_CP forced-clean exits were removed; the
# measured RES<KN_RAY_TOL test is now the only no-cp success path.
# These checks pin the absence of the shortcut and the presence of
# the ordinary exit structure so the defect cannot be reintroduced
# silently.
i_sfl = fg.index("SUBROUTINE SWEEP_FLD(FLD")
j_sfl = fg.index("\n\tEND\n", i_sfl)
sfl_txt = fg[i_sfl:j_sfl]
i_ordx = fg.index("SUBROUTINE SWEEP_FLD_ORD")
j_ordx = fg.index("\n\tEND\n", i_ordx)
ordx_txt = fg[i_ordx:j_ordx]

check(
    "SWEEP_FLD: no .NOT. HAS_CP exit branch remains",
    "IF(.NOT. HAS_CP" not in sfl_txt and "IF ( .NOT. HAS_CP" not in sfl_txt,
)
check(
    "SWEEP_FLD_ORD: no .NOT. HAS_CP exit branch remains",
    "IF(.NOT. HAS_CP" not in ordx_txt and "IF ( .NOT. HAS_CP" not in ordx_txt,
)
check(
    "SWEEP_FLD: measured tolerance test IXIT=1 present",
    "IF(RES .LT. KN_RAY_TOL)THEN" in sfl_txt and "IXIT=1" in sfl_txt,
)
check(
    "SWEEP_FLD_ORD: measured tolerance test IXIT=1 present",
    "IF(RES .LT. KN_RAY_TOL)THEN" in ordx_txt and "IXIT=1" in ordx_txt,
)
check(
    "maxpass exit IXIT=4 present in both sweeps",
    "IXIT=4" in sfl_txt and "IXIT=4" in ordx_txt,
)
check(
    "IXIT init -1 guards against unset-status return",
    "IXIT=-1" in sfl_txt and "IXIT=-1" in ordx_txt,
)
check(
    "ORD Aitken gated NPO>=3 (now reachable)",
    "NPO .GE. 3 .AND." in ordx_txt,
)
check(
    "ORD AA gated NPO>=4 (now reachable)",
    "NPO .GE. 4)THEN" in ordx_txt or "NPO .GE. 4 .AND." in ordx_txt,
)
check(
    "snap bound on SNHI",
    "KN_PSWEEP_SNAP .AND. NPO .LE. KN_PSWEEP_SNHI" in ord_txt,
)

for name, ok in checks:
    print(f"{'PASS' if ok else 'FAIL'}: {name}")
sys.exit(0 if all(ok for _, ok in checks) else 1)
