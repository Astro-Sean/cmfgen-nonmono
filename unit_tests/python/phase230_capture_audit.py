"""Validate selected production captures before any P2-on replay."""

import argparse
import csv
import importlib.util
import json
import struct
from pathlib import Path

import numpy as np
from phase227_capture import sha256
from phase227_records import normalize_omega


def stage_counts(path, refresh):
    header = struct.Struct("<i64sq5i2d")
    counts = {}
    with path.open("rb") as stream:
        while marker := stream.read(4):
            size = struct.unpack("<i", marker)[0]
            payload = stream.read(size)
            if stream.read(4) != marker or size != header.size:
                raise ValueError(f"{path}: invalid stage header")
            magic, label, rf, ray, channel, step, nr, nc, res, omega = header.unpack(
                payload
            )
            if magic != 20261127:
                raise ValueError(f"{path}: stage magic {magic}")
            marker = stream.read(4)
            size = struct.unpack("<i", marker)[0]
            if size != nr * nc * 8:
                raise ValueError(f"{path}: stage field size")
            stream.seek(size, 1)
            if stream.read(4) != marker:
                raise ValueError(f"{path}: stage field marker")
            name = label.decode("ascii").strip()
            if rf == refresh and name in {"M_SOLVE_END", "FLDP_END_DUMP"}:
                if channel in counts:
                    raise ValueError(f"{path}: duplicate channel end")
                counts[channel] = step
    if set(counts) != {1, 2}:
        raise ValueError(f"{path}: missing solve-end counts {counts}")
    return counts[1], counts[2]


def audit(source, destination, refreshes):
    repo = Path(__file__).resolve().parents[2]
    parser_path = repo / "docs/audits/knop_phase226/parse_v2.py"
    spec = importlib.util.spec_from_file_location("parse_v2", parser_path)
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    destination.mkdir(parents=True, exist_ok=False)
    results = []
    for refresh in refreshes:
        for ray in (81, 82, 84):
            path = source / f"KNRDUMP_L{ray:03d}_r{refresh}.bin"
            entry = {"ray": ray, "refresh": refresh, "deck": str(path)}
            try:
                records = list(parser.read_records(path))
                block, end = parser.parse_block(records, 0)
                if end != len(records):
                    raise ValueError("expected exactly one complete deck block")
                if block["ls"] != ray:
                    raise ValueError("deck ray mismatch")
                for key, value in block.items():
                    if isinstance(value, np.ndarray) and not np.isfinite(value).all():
                        raise ValueError(f"nonfinite deck data: {key}")
                if block["fbits"] & 2:
                    raise ValueError("capture unexpectedly has P2 enabled")
                counts = stage_counts(source / f"KNSTAGE_L{ray:03d}.bin", refresh)
                if counts[0] != block["np_p"]:
                    raise ValueError("stage/deck P pass count mismatch")
                schedule = destination / f"omega_L{ray:03d}_r{refresh}.txt"
                normalize_omega(
                    source / f"KNAPPLY_L{ray:03d}.txt",
                    schedule,
                    ray,
                    refresh,
                    counts,
                )
                entry.update(
                    status="VALIDATED_SELECTED_INPUT",
                    records=len(records),
                    sha256=sha256(path),
                    np_p=counts[0],
                    np_m=counts[1],
                    res_p=block["res_p"],
                    ixit_p=block["ixit_p"],
                    ni=block["ni"],
                    ncf=block["ncf"],
                    gam_min=float(block["gam"].min()),
                    fbits=block["fbits"],
                    schedule=str(schedule),
                    schedule_sha256=sha256(schedule),
                    apply_sha256=sha256(source / f"KNAPPLY_L{ray:03d}.txt"),
                )
            except (OSError, ValueError, AssertionError, IndexError) as error:
                entry.update(status="BLOCKED", error=str(error))
            results.append(entry)
    (destination / "inventory.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="ascii"
    )
    with (destination / "census_selected.tsv").open("w", encoding="ascii") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["refresh", "ray", "class", "exit", "nbad", "residual"])
        for refresh in refreshes:
            for line in (
                (source / f"KNRAY_CENSUS_r{refresh}.txt").read_text().splitlines()
            ):
                fields = line.split()
                writer.writerow([fields[i] for i in (0, 1, 2, 4, 5, 6)])
    return int(any(entry["status"] == "BLOCKED" for entry in results))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--refreshes", type=int, nargs="+", default=[1, 10, 250])
    args = parser.parse_args()
    return audit(args.source.resolve(), args.destination.resolve(), args.refreshes)


if __name__ == "__main__":
    raise SystemExit(main())
