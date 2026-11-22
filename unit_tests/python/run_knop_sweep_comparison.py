"""Compare existing production sweep orderings without changing ray physics."""

import argparse
import importlib.util
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, compile_and_link_prod
from phase227_capture import sha256
from run_knop_replay import HELPERS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--refresh", type=int, default=250)
    parser.add_argument("--rays", type=int, nargs="+", default=[81, 82, 84])
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    args.destination.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location(
        "parse_v2", TREE / "docs/audits/knop_phase226/parse_v2.py"
    )
    deck_parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(deck_parser)
    exe = compile_and_link_prod(
        f"knray_sweep_{args.destination.name}_o3",
        [TREE / name for name in HELPERS] + [F90DIR / "drv_knop_replay.f90"],
        extra_flags=["-O3"],
    )
    results = {
        "executable_sha256": sha256(exe),
        "omp_threads": os.getenv("OMP_NUM_THREADS"),
        "rays": [],
    }
    for ray in args.rays:
        deck = args.source / f"KNRDUMP_L{ray:03d}_r{args.refresh}.bin"
        records = list(deck_parser.read_records(deck))
        block, end = deck_parser.parse_block(records, 0)
        if end != len(records):
            raise ValueError(f"{deck}: trailing records")
        ni, ncf = block["ni"], block["ncf"]
        entry = {
            "ray": ray,
            "deck_sha256": sha256(deck),
            "ni": ni,
            "ncf": ncf,
            "gam_min": float(block["gam"].min()),
            "gam_max": float(block["gam"].max()),
            "gam_negative_cells": int(np.count_nonzero(block["gam"] < 0)),
            "production_residual": block["res_p"],
            "production_exit": block["ixit_p"],
            "fbits": block["fbits"],
            "runs": [],
        }
        results["rays"].append(entry)
        production_m = block["fld_m"].ravel()
        production_p = block["fld_p_end"].ravel()
        for sweep in ("ordered", "frequency"):
            for repeat in (1, 2):
                directory = args.destination / f"L{ray:03d}_{sweep}_{repeat}"
                directory.mkdir()
                logging.info("ray=%d sweep=%s repeat=%d", ray, sweep, repeat)
                command = [str(exe), str(deck), "sweep", sweep, "p2=0"]
                with (directory / "replay.log").open("w") as stream:
                    process = subprocess.run(
                        command,
                        cwd=directory,
                        stdout=stream,
                        stderr=subprocess.STDOUT,
                        timeout=1800,
                        check=False,
                    )
                output = directory / "KN_REPLAY_END.bin"
                run = {
                    "sweep": sweep,
                    "repeat": repeat,
                    "returncode": process.returncode,
                    "result_lines": [
                        line
                        for line in (directory / "replay.log").read_text().splitlines()
                        if line.startswith("RESULT")
                    ],
                }
                if process.returncode != 0 or not output.is_file():
                    raise RuntimeError(f"{directory}: failed replay")
                fields = list(deck_parser.read_records(output))
                p = np.frombuffer(fields[1], dtype="<f8")
                m = np.frombuffer(fields[2], dtype="<f8")
                if p.size != ni * ncf or m.size != p.size:
                    raise ValueError(f"{output}: field size")
                if not np.isfinite(p).all() or not np.isfinite(m).all():
                    raise ValueError(f"{output}: nonfinite field")
                run.update(
                    field_sha256=sha256(output),
                    p_vs_production_abs_max=float(np.max(np.abs(p - production_p))),
                    m_vs_production_abs_max=float(np.max(np.abs(m - production_m))),
                    im_ni_vs_production_abs_max=float(
                        np.max(
                            np.abs(
                                m.reshape((ni, ncf), order="F")[-1]
                                - production_m.reshape((ni, ncf), order="F")[-1]
                            )
                        )
                    ),
                )
                entry["runs"].append(run)
                (args.destination / "summary.json").write_text(
                    json.dumps(results, indent=2) + "\n"
                )
        del records, block, production_m, production_p
    logging.info("sweep comparison complete: %s", args.destination)


if __name__ == "__main__":
    main()
