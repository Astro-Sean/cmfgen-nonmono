"""Phase 2.1A manufactured validity tests for JRES_TRIG.

Builds the driver against the production jres_validity.o and
asserts the trigger classification for each manufactured field:
finite-converged, interior/boundary NaN, interior Inf, over-cap,
FEDD>2, and finite-but-unconverged (the stock keep-case).
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import TREE, compile_and_link_prod

# (trigger, BAD): BAD is the first non-finite index for trig 1-2
# and the argmax for trig 3 (and harmlessly also when trig=0).
EXPECT = {
    "converged": (0, 1),
    "nan_interior": (1, 5),
    "nan_first": (1, 1),
    "inf_interior": (2, 7),
    "overcap_int": (3, 6),
    "fedd_over2": (3, 3),
    "unconverged": (0, 4),
}


def main():
    drv = compile_and_link_prod(
        "jres_trig",
        [
            Path(__file__).resolve().parents[1] / "f90" / "drv_jres_trig.f",
            TREE / "obs" / "jres_validity.f",
        ],
    )
    print(f"built {drv}")
    out = subprocess.run([str(drv)], capture_output=True, text=True, check=True)
    got = {}
    for line in out.stdout.splitlines():
        parts = line.split()
        got[parts[0]] = (int(parts[1]), int(parts[2]))
    npass = 0
    for name, want in EXPECT.items():
        assert got[name][:2] == want, (
            f"{name}: got {got[name]} want {want}"
        )
        print(f"[PASS] {name}: trig={got[name][0]} bad={got[name][1]}")
        npass += 1
    print(f"{npass}/{npass} jres_trig cases pass")


if __name__ == "__main__":
    main()
