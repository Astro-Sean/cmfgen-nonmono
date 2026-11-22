"""Check stage-count validation without interpreting a transfer operator."""

import struct
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import Report
from phase230_capture_audit import stage_counts


def record(payload):
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


def event(label, channel, step, magic=20261127):
    header = struct.pack(
        "<i64sq5i2d",
        magic,
        label.encode("ascii").ljust(64, b" "),
        250,
        81,
        channel,
        step,
        1,
        1,
        0.01,
        1.0,
    )
    return record(header) + record(struct.pack("<d", 1.0))


def main():
    rep = Report("phase230_capture_audit")
    m = event("M_SOLVE_END", 2, 4)
    p = event("FLDP_END_DUMP", 1, 18)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "stages.bin"
        path.write_bytes(m + p)
        rep.check(
            "independent P/M solve-end counts", stage_counts(path, 250) == (18, 4)
        )
        for name, data in (
            ("missing M end rejected", p),
            ("duplicate end rejected", m + p + p),
            ("bad stage magic rejected", event("M_SOLVE_END", 2, 4, 0) + p),
            ("broken field marker rejected", (m + p)[:-4] + b"xxxx"),
        ):
            path.write_bytes(data)
            try:
                stage_counts(path, 250)
            except ValueError:
                rep.check(name, True)
            else:
                rep.check(name, False)
    return rep.summary()


if __name__ == "__main__":
    raise SystemExit(main())
