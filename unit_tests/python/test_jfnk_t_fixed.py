#!/usr/bin/env python3
"""Verify that the JFNK temperature lock honors the VADAT FIX_T policy."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, Report, compile_driver, run


def main():
    report = Report("JFNK fixed-temperature policy")
    binary = compile_driver(
        "drv_jfnk_t_fixed",
        [
            TREE / "new_main" / "mod_subs" / "jfnk_t_fixed_mode_v1.f",
            F90DIR / "drv_jfnk_t_fixed.f90",
        ],
    )
    result = run(binary)
    report.check(
        "production fixed-T predicate passes all four mode combinations",
        result.returncode == 0,
        result.stdout.strip(),
    )
    return report.summary()


if __name__ == "__main__":
    sys.exit(main())
