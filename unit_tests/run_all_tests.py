#!/usr/bin/env python3
"""Run every python/test_*.py in sequence and report the aggregate.

Each test is an independent process that compiles the exact
production Fortran sources from the tree and drives them through
small CLI drivers.  Exit status is the number of failed test files.
"""

import subprocess
import sys
import time
from pathlib import Path

UT = Path(__file__).resolve().parent
TESTS = sorted((UT / "python").glob("test_*.py"))


def main():
    results = []
    for t in TESTS:
        t0 = time.time()
        p = subprocess.run([sys.executable, str(t)], capture_output=True, text=True)
        dt = time.time() - t0
        # Surface the test's own report block (everything from the
        # banner onward) plus any traceback tail on hard errors.
        out = p.stdout
        banner = out.find("=" * 60)
        body = out[banner:] if banner >= 0 else out
        print(body, end="")
        if p.returncode != 0 and p.stderr:
            print(p.stderr[-2000:])
        results.append((t.name, p.returncode, dt))
        print(f"  -> {t.name}: rc={p.returncode} ({dt:.0f}s)")

    print(f"\n{'=' * 60}\nSUITE SUMMARY\n{'=' * 60}")
    nfail = 0
    for name, rc, dt in results:
        tag = "PASS" if rc == 0 else "FAIL"
        if rc != 0:
            nfail += 1
        print(f"  [{tag}] {name:<32} rc={rc} ({dt:.0f}s)")
    print(f"{'-' * 60}")
    print(f"{len(results) - nfail}/{len(results)} test files passed")
    return nfail


if __name__ == "__main__":
    sys.exit(main())
