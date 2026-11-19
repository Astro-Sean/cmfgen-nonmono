"""Independent KNOP P-endpoint (j=1) slope-reconstruction checks.

Builds unit_tests/f90/drv_knop_endpoint.f, a standalone driver that
reconstructs the endpoint source derivative on manufactured S(tau)
cases by two paths: the production secant-extrapolation+limiter
formula (subject under test) and an independent Lagrange-basis
derivative (reference). Asserts exactness where the mathematics
requires it and bounded/order-consistent error elsewhere.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import compile_and_link_prod

DRIVER = Path(__file__).resolve().parents[1] / "f90" / "drv_knop_endpoint.f"


def rows(case):
    drv = compile_and_link_prod("knop_endpoint", [DRIVER])
    out = subprocess.run(
        [str(drv)], capture_output=True, text=True, check=True
    )
    data = {}
    for line in out.stdout.splitlines():
        tok = line.split()
        if len(tok) >= 13 and tok[0] in "ABCDE":
            c, lv = tok[0], int(tok[1])
            data.setdefault(c, {})[lv] = {
                "h": float(tok[2]),
                "ana": float(tok[3]),
                "prod": float(tok[4]),
                "err_p": float(tok[5]),
                "rel_p": float(tok[6]),
                "ref": float(tok[7]),
                "err_r": float(tok[8]),
                "dupd": float(tok[12]),
            }
    return data


def main():
    data = rows(None)
    npass = 0

    for lv, r in data["A"].items():
        assert r["err_p"] <= 1e-10, f"A lv{lv} prod err {r['err_p']}"
        assert r["err_r"] <= 1e-10, f"A lv{lv} ref err {r['err_r']}"
    npass += 1
    print("[PASS] A constant: prod and ref exact")

    for lv, r in data["B"].items():
        assert r["err_p"] <= 1e-10, f"B lv{lv} prod err {r['err_p']}"
        assert r["err_r"] <= 1e-10, f"B lv{lv} ref err {r['err_r']}"
    npass += 1
    print("[PASS] B linear: prod and ref exact")

    for lv, r in data["C"].items():
        assert r["err_p"] <= 1e-10, f"C lv{lv} prod err {r['err_p']}"
        assert r["err_r"] <= 1e-10, f"C lv{lv} ref err {r['err_r']}"
    npass += 1
    print("[PASS] C quadratic: secant extrapolation is node-exact")

    import math

    es = [data["D"][lv]["err_p"] for lv in sorted(data["D"])]
    hs = [data["D"][lv]["h"] for lv in sorted(data["D"])]
    assert all(e < 1e-3 for e in es), f"D errors {es}"
    for i in range(1, len(es)):
        order = math.log(es[i - 1] / es[i]) / math.log(hs[i - 1] / hs[i])
        assert 1.9 < order < 2.1, f"D lv{i + 1} order {order}"
    npass += 1
    print("[PASS] D smooth: O(h^2) error, second-order refinement")

    es = [data["E"][lv]["err_p"] for lv in sorted(data["E"])]
    assert all(e < 1e-2 for e in es), f"E errors {es}"
    assert es == sorted(es, reverse=True), f"E not shrinking {es}"
    for lv, r in data["E"].items():
        assert r["dupd"] == r["err_p"], f"E lv{lv} dupd mismatch"
    npass += 1
    print("[PASS] E extremum: bounded, shrinking under refinement")

    for c, lvs in data.items():
        for lv, r in lvs.items():
            for k in ("prod", "ref", "dupd"):
                v = r[k]
                assert v == v and abs(v) != float("inf"), (
                    f"{c} lv{lv} {k} non-finite"
                )
    npass += 1
    print("[PASS] all cases finite")

    print(f"{npass}/{npass} knop_endpoint checks pass")


if __name__ == "__main__":
    main()
