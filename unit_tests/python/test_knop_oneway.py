"""Check production sweep dispatch without replacing the physical gradient."""

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, WORKDIR, Report, compile_and_link_prod, run
from test_knop_ray import HELPERS, parse


def main():
    rep = Report("knop_oneway")
    exe = compile_and_link_prod(
        "knray_oneway",
        [TREE / name for name in HELPERS] + [F90DIR / "drv_knop_ray.f90"],
    )
    for name, expected in (("gposo", 0), ("gzzo", 0), ("gnego", 1), ("ord", 1)):
        result = run(exe, name, cwd=WORKDIR, check=False)
        dispatch = re.search(r"DISPATCH ordered=(\d+) frequency=(\d+)", result.stdout)
        rep.check(
            f"{name} selects expected sweep",
            result.returncode == 0
            and dispatch is not None
            and tuple(map(int, dispatch.groups())) == (expected, 1 - expected),
        )
    for requested, reference in (("gposo", "gpos"), ("gzzo", "gzzf")):
        outputs = [
            run(exe, name, cwd=WORKDIR).stdout for name in (requested, reference)
        ]
        parsed = [parse(output) for output in outputs]
        fields = [
            np.array([value for key, value in store.items() if isinstance(key[0], int)])
            for _, store in parsed
        ]
        rep.check(
            f"{requested} preserves frequency-major field",
            fields[0].shape == fields[1].shape
            and np.isfinite(fields[0]).all()
            and np.allclose(fields[0], fields[1], rtol=1e-13, atol=0),
        )
        rep.check(
            f"{requested} retains converged exit",
            parsed[0][0].get("NXIT") == parsed[1][0].get("NXIT") == 1,
        )
    return rep.summary()


if __name__ == "__main__":
    raise SystemExit(main())
