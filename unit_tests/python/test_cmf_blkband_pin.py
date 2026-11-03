#!/usr/bin/env python3
"""Test the production CMF_BLKBAND null-equation completion."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, compile_driver, run


def main():
    binary = compile_driver(
        "drv_cmf_blkband_pin",
        [
            TREE / "tools" / "set_kind_module.f",
            TREE / "new_main" / "subs" / "cmf_blkband_pin_null_v1.f",
            F90DIR / "drv_cmf_blkband_pin.f90",
        ],
    )
    result = run(binary)
    print(result.stdout, end="")
    if "CMF_BLKBAND_PIN_NULL_V1: PASS" not in result.stdout:
        raise SystemExit("CMF_BLKBAND_PIN_NULL_V1 regression failed")


if __name__ == "__main__":
    main()
