"""Prepare or run actual-Fortran diagnostic IO tests in an isolated allocation."""

import argparse
import json
import os
import struct
import subprocess
import sys
from pathlib import Path

from phase227_capture import save_json, sha256, utc_now
from phase227_records import census


def row(step, value=0.5141234567890123, channel=1, ray=81, direction=0):
    bits = struct.pack(">d", value).hex().upper()
    return f"0 {ray} {channel} {direction} {step} {value:.17e} {value:.17e} {bits}\n"


def fixtures():
    header = "KNOMG2 81 0 2 1\n"
    valid = header + row(1) + row(2) + row(1, channel=2)
    cases = {
        "valid": valid,
        "missing_pass": header + row(1) + row(1, channel=2),
        "duplicate": header + row(1) + row(1) + row(1, channel=2),
        "wrong_ray": header + row(1, ray=82) + row(2) + row(1, channel=2),
        "wrong_channel": header + row(1, channel=3) + row(2) + row(1, channel=2),
        "wrong_direction": header + row(1, direction=1) + row(2) + row(1, channel=2),
        "nan": header + row(1, float("nan")) + row(2) + row(1, channel=2),
        "infinity": header + row(1, float("inf")) + row(2) + row(1, channel=2),
        "malformed": header + "not a numeric record\n",
        "roundtrip_failure": valid.replace("5.14123456789012336e-01", "0.514", 2),
        "extra": valid + row(3),
    }
    parts = valid.splitlines()
    tokens = parts[1].split()
    tokens[6] = "0.514"
    parts[1] = " ".join(tokens)
    cases["roundtrip_failure"] = "\n".join(parts) + "\n"
    return cases


def prepare(destination, binary):
    destination.mkdir(parents=True, exist_ok=False)
    paths = {}
    for name, text in fixtures().items():
        path = destination / f"{name}.txt"
        path.write_text(text, encoding="ascii")
        paths[name] = str(path)
    paths["missing_file"] = str(destination / "does_not_exist.txt")
    save_json(
        destination / "test_plan.json",
        {
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "binary": str(binary.resolve()),
            "binary_sha256": sha256(binary),
            "cases": paths,
            "fixture_hashes": {
                name: sha256(path)
                for name, path in paths.items()
                if Path(path).exists()
            },
            "status": "PREPARED_NOT_RUN",
            "execution_requires": "existing isolated allocation",
        },
    )


def run(plan_path):
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("No allocation: tests must not run on a login node")
    plan = json.loads(plan_path.read_text())
    binary = Path(plan["binary"])
    if sha256(binary) != plan["binary_sha256"]:
        raise ValueError("Test executable hash mismatch")
    root = plan_path.parent / "results"
    root.mkdir(exist_ok=False)
    outcomes = []
    for name, path in plan["cases"].items():
        process = subprocess.run(
            [str(binary), "omega", path],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        expected_success = name == "valid"
        passed = (process.returncode == 0) == expected_success
        if not expected_success:
            passed = passed and "[FAIL] file=" in process.stderr
        outcomes.append(
            {
                "name": name,
                "passed": passed,
                "returncode": process.returncode,
                "stdout": process.stdout,
                "stderr": process.stderr,
            }
        )
    output = root / "census.txt"
    subprocess.run([str(binary), "census", str(output)], cwd=root, check=True)
    counts = {"FULL": 1, "BULK": 1, "ITLIM": 0, "STAG": 1, "DIV": 0, "INV": 1}
    census(output, range(1, 5), {0: counts})
    lines = output.read_text().splitlines()
    mutations = {
        "missing_ray": lines[:-1],
        "duplicate_ray": lines + [lines[0]],
        "invalid_class": [lines[0].replace("0 1 0 ", "0 1 9 ", 1)] + lines[1:],
        "malformed_census": ["bad row"] + lines[1:],
    }
    for name, text in mutations.items():
        path = root / f"{name}.txt"
        path.write_text("\n".join(text) + "\n")
        try:
            census(path, range(1, 5), {0: counts})
            rejected = False
        except ValueError:
            rejected = True
        outcomes.append({"name": name, "passed": rejected})
    try:
        census(output, range(1, 5), {0: {**counts, "FULL": 2}})
        rejected = False
    except ValueError:
        rejected = True
    outcomes.append({"name": "aggregate_mismatch", "passed": rejected})
    subprocess.run([str(binary), "roundtrip"], cwd=root, check=True)
    application = (root / "KNAPPLY_L081.txt").read_text().split()
    roundtrip = struct.pack(">d", float(application[6])).hex().upper() == application[7]
    outcomes.append({"name": "writer_binary64_roundtrip", "passed": roundtrip})
    save_json(
        root / "results.json",
        {
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "plan_sha256": sha256(plan_path),
            "binary_sha256": sha256(binary),
            "job_id": os.environ["SLURM_JOB_ID"],
            "results": outcomes,
            "census_valid_and_boundaries": "PASS",
            "status": (
                "PASS" if all(result["passed"] for result in outcomes) else "FAIL"
            ),
        },
    )
    return 0 if all(result["passed"] for result in outcomes) else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    try:
        if args.prepare and args.binary:
            prepare(args.prepare.resolve(), args.binary.resolve())
            sys.stdout.write("[ok] Test inputs prepared; nothing executed\n")
            return 0
        if args.run:
            return run(args.run.resolve())
        parser.error("provide --prepare/--binary or --run")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.stderr.write(f"[FAIL] {error}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
