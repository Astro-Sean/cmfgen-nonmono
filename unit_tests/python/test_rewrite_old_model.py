#!/usr/bin/env python3
"""REWRITE_OLD_MODEL_DATA conversion pass (misc/rewrite_old_model_data.f).

Alteration under test: the program now ends with STOP rather than
RETURN so the shell reliably sees a normal exit after the copy.

The production executable (rebuilt post-alteration) is run against a
fabricated OLD_MODEL_DATA in the 8-byte sequential format:

  DATE; ND,NUM_SPECIES; SN_AGE; 7 x VEC(ND);
  per species: BEG_ID,END_ID,NAME then per ion:
  ID,NXZV_F,ZXV; XZV_F(NXZV_F,ND); VEC(ND)

The output must hold the same records widened to 16-byte reals.  A
NUM_SPECIES=0 deck checks the header path; a one-species deck checks
the ion blocks are copied through.
"""
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

EXE = TREE / "exe" / "rewrite_old_model_data.exe"
ND = 4


def rec(fmt, *vals):
    data = struct.pack(fmt, *vals)
    return struct.pack("<i", len(data)) + data + struct.pack("<i", len(data))


def write_fixture(path, n_species):
    out = [
        rec("<11s", b"24-FEB-2004 "),
        rec("<2i", ND, n_species),
        rec("<d", 12.5),
    ]
    for k in range(7):
        out.append(rec("<%dd" % ND, *[10.0 * k + i for i in range(ND)]))
    if n_species:
        out.append(rec("<2i10s", 1, 1, b"HI        "))
        out.append(rec("<2id", 1, 2, 1.0))
        out.append(rec("<%dd" % (2 * ND), *[1.0 + i for i in range(2 * ND)]))
        out.append(rec("<%dd" % ND, *[7.0 + i for i in range(ND)]))
    path.write_bytes(b"".join(out))


def record_sizes(path):
    raw = path.read_bytes()
    sizes = []
    pos = 0
    while pos < len(raw):
        n = struct.unpack_from("<i", raw, pos)[0]
        sizes.append(n)
        pos += 4 + n + 4
    return sizes, pos


def rec_values(path, idx):
    raw = path.read_bytes()
    pos = 0
    for _ in range(idx):
        n = struct.unpack_from("<i", raw, pos)[0]
        pos += 4 + n + 4
    n = struct.unpack_from("<i", raw, pos)[0]
    return struct.unpack_from("<%dd" % (n // 8), raw, pos + 4)


def run_case(n_species):
    with tempfile.TemporaryDirectory() as d:
        write_fixture(Path(d, "OLD_MODEL_DATA"), n_species)
        r = subprocess.run([str(EXE)], cwd=d, capture_output=True,
                           text=True, timeout=60)
        out_file = Path(d, "NEW_OLD_MODEL_DATA")
        exists = out_file.exists()
        sizes, end = record_sizes(out_file) if exists else ([], -1)
        total = len(out_file.read_bytes()) if exists else 0
        vals = rec_values(out_file, 2) if exists else ()
        return r, exists, sizes, end, total, vals


def main():
    rep = Report("REWRITE_OLD_MODEL_DATA pass")
    rep.check("production exe present", EXE.exists())

    r, exists, sizes, end, total, vals = run_case(0)
    rep.check("no-species deck: clean exit", r.returncode == 0)
    rep.check("no-species deck: output written", exists)
    rep.check("no-species deck: all records consumed", end == total)
    rep.check("no-species deck: header records kept",
              sizes[:2] == [11, 8])
    rep.check("no-species deck: LDP-sized records written",
              sizes[2:] == [8] + [8 * ND] * 7)
    rep.check("no-species deck: SN_AGE copied exactly",
              len(vals) == 1 and abs(vals[0] - 12.5) < 1e-12)

    r, exists, sizes, end, total, vals = run_case(1)
    rep.check("one-species deck: clean exit", r.returncode == 0)
    rep.check("one-species deck: output written", exists)
    # header(3) + 7 vecs + species hdr + ion hdr + XZV_F + vec
    rep.check("one-species deck: species block copied",
              len(sizes) == 14 and sizes[10] == 18
              and sizes[11] == 16 and sizes[12] == 8 * 2 * ND)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
