"""Independent KNOP P-endpoint (j=1) slope-reconstruction checks.

Builds unit_tests/f90/drv_knop_endpoint.f, a standalone driver that
reconstructs the endpoint source derivative on manufactured S(tau)
cases by three paths: the production secant-extrapolation+limiter
formula (subject under test), a three-point Lagrange-basis
derivative (independent reference), and a four-point Lagrange
derivative (higher-order reference). Asserts exactness where the
mathematics requires it and bounded/order-consistent error
elsewhere.

DUPD is a unit-representative formal-update discrepancy
(A4REP=1): it scales the measured DS error by a unit coefficient
to show the magnitude a unit-weight production update would see.
It is not the actual PK_A4-weighted production transfer error.
"""

import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import compile_and_link_prod

DRIVER = Path(__file__).resolve().parents[1] / "f90" / "drv_knop_endpoint.f"


def rows():
    drv = compile_and_link_prod("knop_endpoint", [DRIVER])
    out = subprocess.run(
        [str(drv)], capture_output=True, text=True, check=True
    )
    data = {}
    for line in out.stdout.splitlines():
        tok = line.split()
        if len(tok) >= 12 and tok[0] in "ABCDE":
            c, lv = tok[0], int(tok[1])
            data.setdefault(c, {})[lv] = {
                "h": float(tok[2]),
                "ana": float(tok[3]),
                "prod": float(tok[4]),
                "err_p": float(tok[5]),
                "rel_p": float(tok[6]),
                "ref": float(tok[7]),
                "err_r": float(tok[8]),
                "err_3": float(tok[10]),
                "dupd": float(tok[11]),
            }
    return data


def order(errs, hs, i):
    return math.log(errs[i - 1] / errs[i]) / math.log(hs[i - 1] / hs[i])


def seq(data, case, key):
    lvs = sorted(data[case])
    return [data[case][lv][key] for lv in lvs], [
        data[case][lv]["h"] for lv in lvs
    ]


def main():
    data = rows()
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
        assert r["err_3"] <= 1e-10, f"C lv{lv} ref3 err {r['err_3']}"
    npass += 1
    print("[PASS] C quadratic: all reconstructions node-exact")

    ep, hs = seq(data, "D", "err_p")
    er, _ = seq(data, "D", "err_r")
    e3, _ = seq(data, "D", "err_3")
    assert all(e < 1e-3 for e in ep), f"D prod errors {ep}"
    assert all(e < 1e-3 for e in er), f"D ref errors {er}"
    assert all(e == e and abs(e) != float("inf") for e in e3), (
        f"D ref3 non-finite {e3}"
    )
    op = [order(ep, hs, i) for i in range(1, len(ep))]
    orf = [order(er, hs, i) for i in range(1, len(er))]
    assert all(1.9 < o < 2.1 for o in op), f"D prod orders {op}"
    assert all(1.9 < o < 2.1 for o in orf), f"D ref orders {orf}"
    npass += 1
    print(
        f"[PASS] D smooth: prod order {min(op):.2f}-{max(op):.2f},"
        f" ref order {min(orf):.2f}-{max(orf):.2f}, ref3 finite"
    )

    ee, he = seq(data, "E", "err_p")
    assert all(e < 1e-2 for e in ee), f"E errors {ee}"
    assert all(
        ee[i] < ee[i - 1] for i in range(1, len(ee))
    ), f"E not monotonically shrinking {ee}"
    oe = [order(ee, he, i) for i in range(1, len(ee))]
    assert all(1.0 < o < 5.0 for o in oe), f"E orders {oe}"
    for lv, r in data["E"].items():
        assert r["dupd"] == r["err_p"], f"E lv{lv} dupd mismatch"
    npass += 1
    print(
        f"[PASS] E extremum: bounded, monotone decrease,"
        f" observed order {min(oe):.2f}-{max(oe):.2f} (limiter"
        " region, no fixed order demanded)"
    )

    for c, lvs in data.items():
        for lv, r in lvs.items():
            for k in ("prod", "ref", "err_3", "dupd"):
                v = r[k]
                assert v == v and abs(v) != float("inf"), (
                    f"{c} lv{lv} {k} non-finite"
                )
    npass += 1
    print("[PASS] all cases finite")

    print(f"{npass}/{npass} knop_endpoint checks pass")


if __name__ == "__main__":
    main()
