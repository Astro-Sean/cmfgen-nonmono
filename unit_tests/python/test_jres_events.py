"""Phase 2.2 manufactured tests for the J_RESCUE_EVENTS writer.

Builds drv_jres_events against the production jres_validity.o,
runs it in a scratch directory, and parses the JEV/JTOT/JEND
records for positional-field sanity: one NES_SAN in CONTINUUM
(trig=3, ltc1=T on SIGMA_TRUE=-1.5) and one NRESC in LINE.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import TREE, compile_and_link_prod


def main():
    drv = compile_and_link_prod(
        "jres_events",
        [
            Path(__file__).resolve().parents[1] / "f90" / "drv_jres_events.f",
            TREE / "obs" / "jres_validity.f",
        ],
    )
    print(f"built {drv}")
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([str(drv)], cwd=td, check=True)
        lines = (Path(td) / "J_RESCUE_EVENTS").read_text().splitlines()

    recs = [l.split() for l in lines if l.startswith("JEV")]
    tots = [l.split() for l in lines if l.startswith("JTOT")]
    ends = [l.split() for l in lines if l.startswith("JEND")]
    assert len(recs) == 2, f"expected 2 JEV records, got {len(recs)}"
    assert len(ends) == 1, f"expected 1 JEND record, got {len(ends)}"

    # NES_SAN record: type indx ml fl sect first fg mom L fdnc
    # trig ftrig bad rjmax fedmax rjesmax rcap rcapit snapok
    # djmax dfmax R V SIG SIGT ltc1 luedd
    ev = recs[0]
    assert ev[0] == "JEV" and ev[1] == "NES_SAN"
    assert ev[5] == "CONTINUUM" and ev[6] == "F"
    assert ev[7] == "MOM_J_CMF_V8" and ev[9] == "0" and ev[10] == "F"
    assert ev[11] == "3" and ev[13] == "12"
    assert abs(float(ev[16]) - 2.0e5) < 1.0  # rjesmax
    assert ev[26] == "T" and ev[27] == "T"  # ltc1, luedd
    print(f"[PASS] nesan record: {' '.join(ev[:8])}")

    ev = recs[1]
    assert ev[1] == "NRESC" and ev[5] == "LINE"
    assert ev[7] == "FG_J_CMF_V12" and ev[8] == "LAM_ES_TC"
    assert ev[9] == "14" and ev[10] == "T"  # L, fdnc
    assert ev[19] == "T"  # snapok
    assert abs(float(ev[20]) - 8.0e6) < 1.0  # djmax
    assert ev[26] == "F" and ev[27] == "T"  # ltc1 F (sigt=4), luedd
    print(f"[PASS] nresc record: {' '.join(ev[:8])}")

    # CONTINUUM JTOT flushed on the transition to LINE.
    assert any(t[0] == "JTOT" and t[1] == "CONTINUUM" for t in tots), tots
    assert any(t[0] == "JTOT" and t[1] == "LINE" for t in tots), tots
    e = ends[0]
    # 'key=' and value are separate whitespace tokens
    assert e[1] == "nresc_cont=" and e[2] == "0"
    assert e[3] == "nesan_cont=" and e[4] == "1"
    assert e[5] == "nresc_line=" and e[6] == "1"
    assert e[7] == "nesan_line=" and e[8] == "0"
    print(f"[PASS] totals: {' '.join(ends[0])}")
    print("jres_events: all checks pass")


if __name__ == "__main__":
    main()
