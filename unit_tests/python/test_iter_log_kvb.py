#!/usr/bin/env python3
"""Unit tests for the ilf_kvb / ilf_v* vertical key/value block
helpers in new_main/mod_subs/iter_log_fmt_mod.f90.

The block format converts multi-field one-liners: header line kept
verbatim, one 'label = value' row per field, labels padded to the
longest in the block, values right-aligned (I10 integers, ES reals
at the caller's decimals with a two-digit exponent, T/F logicals).
LOG_COMPACT_KV collapses the block to the legacy one-liner.

Driven end-to-end: drv_ilf_kvb.f90 links the production module and
emits the edge cases (1E+99, 1E-99, 0, negatives, |x|>=1E100,
|x|<1E-99, NaN, +/-Inf, a 10-digit integer, an empty label list,
labels past 30 characters, a block that hits the 100-column cap).
Static pins cover the export list and the layout rules.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (
    BINDIR,
    F90DIR,
    PROD_MOD_INC,
    TREE,
    Report,
    compile_objs,
    prod_objs,
    link_prod,
    run,
)

ILF_SRC = TREE / "new_main" / "mod_subs" / "iter_log_fmt_mod.f90"
PROD = [
    TREE / "tools" / "set_kind_module.f",
    ILF_SRC,
]


def build():
    md = BINDIR / "obj_ilf_kvb"
    objs = prod_objs(PROD)
    objs += compile_objs(
        "ilf_kvb",
        [F90DIR / "drv_ilf_kvb.f90"],
        moddir=md,
        extra_flags=PROD_MOD_INC,
    )
    return link_prod("ilf_kvb", objs, moddir=md)


def cases(out):
    """Split driver stdout into {name: [lines]}."""
    blocks, name = {}, None
    for ln in out.splitlines():
        if ln.startswith("CASE "):
            name = ln.split(None, 1)[1]
            blocks[name] = []
        elif name is not None:
            blocks[name].append(ln)
    return blocks


def eq_col(rows):
    """Column of the '=' character in each row; None if absent."""
    cols = [r.find("=") for r in rows if "=" in r]
    return cols[0] if cols and all(c == cols[0] for c in cols) else None


def main():
    rep = Report("iter_log_kvb")
    ilf = ILF_SRC.read_text()

    rep.check(
        "module exports ilf_kvb + ilf_v* + log_compact_kv",
        all(
            re.search(rf"public ::[^!]*\b{name}\b", ilf)
            for name in (
                "ilf_kvb", "ilf_vi", "ilf_vr", "ilf_vl", "ilf_vs",
                "log_compact_kv",
            )
        ),
    )
    rep.check(
        "compact switch defaults .false.",
        re.search(r"log_compact_kv\s*=\s*\.false\.", ilf) is not None,
    )
    rep.check(
        "block printer present and uses the issue census for wraps",
        "subroutine ilf_kvb" in ilf and "ilf_issue_add" in ilf,
    )

    binary = build()
    r = run(binary)
    blocks = cases(r.stdout)
    # The 100-column budget binds block layout only; compact mode
    # reproduces legacy one-liners at their original width.
    allout = [
        ln
        for name, rows in blocks.items()
        if name != "compact"
        for ln in rows
    ]
    rep.check(
        "every emitted line <= 100 columns",
        all(len(ln) <= 100 for ln in allout),
        f"worst {max(len(ln) for ln in allout)}",
    )
    rep.check(
        "no '*'-overflow and no trailing whitespace anywhere",
        all("*" not in ln and ln == ln.rstrip() for ln in allout),
    )

    rows = blocks["spec"]
    rep.check("spec: header kept verbatim", rows[0] == " ES_ALI:")
    rep.check("spec: five field rows", len(rows[1:]) == 5)
    c = eq_col(rows[1:])
    rep.check("spec: '=' aligned (indent 4 + pad 17 + '  ')", c == 23)
    rep.check(
        "spec: I10 integer field right-aligned",
        rows[1].endswith("        98") and rows[3].endswith("         0"),
    )
    rep.check(
        "spec: ES11.2E2 real field keeps decimals",
        rows[2].endswith("  1.96E-03"),
    )
    rep.check(
        "spec: labels padded to the block maximum",
        rows[2][4:21] == "max final dJ/J   " and rows[1][4:21] == "freqs unconverged",
    )

    rows = blocks["edges"][1:]
    rep.check(
        "edges: 1E+99 / 1E-99 keep two-digit exponents",
        "1.00E+99" in rows[0] and "1.00E-99" in rows[1],
    )
    rep.check(
        "edges: zero and negative reals",
        "0.00E+00" in rows[2] and "-4.500E+03" in rows[3],
    )
    rep.check(
        "edges: |x|>=1E100 and |x|<1E-99 widen to three-digit exponents",
        "1.23E+100" in rows[4] and "1.20E-101" in rows[5],
    )
    rep.check(
        "edges: NaN / +Inf / -Inf print as text, right-aligned",
        rows[6].endswith("NaN")
        and rows[7].endswith("Inf")
        and rows[8].endswith("-Inf"),
    )

    rows = blocks["intlog"][1:]
    rep.check(
        "intlog: 10-digit integer fills I10",
        rows[0].endswith("1234567890"),
    )
    rep.check("intlog: negative integer", rows[1].endswith("       -42"))
    rep.check(
        "intlog: logicals print T/F right-aligned",
        rows[2].endswith("T") and rows[3].endswith("F"),
    )
    rep.check("intlog: string field kept", rows[4].endswith("= (BA computed)"))

    rep.check(
        "empty: header prints alone",
        blocks["empty"] == ["EMPTY:"],
    )

    rows = blocks["longlab"][1:]
    c = eq_col(rows)
    rep.check(
        "longlab: >30-char label sets the pad, '=' still aligned",
        c == 4 + 43 + 2,
    )

    rows = blocks["wide"][1:]
    rep.check(
        "wide: capped pad keeps short rows at <=100 cols",
        len(rows[0]) == 100 and len(rows[1]) == 100 and eq_col(rows[:2]) is not None,
    )
    rep.check(
        "wide: over-long label wraps its value to a continuation row",
        len(rows) == 4
        and rows[2].endswith(" =")
        and rows[3].startswith("        ")
        and rows[3].endswith("         3"),
    )

    rows = blocks["compact"]
    rep.check(
        "compact: block collapses to a single legacy line",
        len(rows) == 1
        and rows[0].startswith(" ES_ALI:")
        and "freqs unconverged=98" in rows[0]
        and "max final dJ/J=1.96E-03" in rows[0],
    )

    # Temporary regression hook: CMFGEN_COMPACT_KV is read once at
    # the first ilf_kvb call.  =1 collapses every block; =0 leaves
    # the module value alone so an explicit set still works.
    env1 = dict(os.environ, CMFGEN_COMPACT_KV="1")
    r1 = subprocess.run(
        [str(binary)], capture_output=True, text=True, env=env1
    )
    rep.check(
        "env: CMFGEN_COMPACT_KV=1 collapses every block",
        " ES_ALI: freqs unconverged=98" in r1.stdout
        and len(cases(r1.stdout)["spec"]) == 1,
    )
    env0 = dict(os.environ, CMFGEN_COMPACT_KV="0")
    r0 = subprocess.run(
        [str(binary)], capture_output=True, text=True, env=env0
    )
    b0 = cases(r0.stdout)
    rep.check(
        "env: CMFGEN_COMPACT_KV=0 keeps block form, explicit set still wins",
        len(b0["spec"]) == 6 and len(b0["compact"]) == 1,
    )

    sys.exit(rep.summary())


if __name__ == "__main__":
    main()
