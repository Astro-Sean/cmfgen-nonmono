#!/usr/bin/env python3
"""Phase-5.3K KN_PSWEEP_DIAG offline verification.

Parses the per-pass records the instrumented SWEEP_FLD_ORD writes to
KN_PSWEEP_DIAG.dat and the optional KN_PSWEEP_XY_r*_p*.bin field
snapshots, then recomputes R_MAP and RMS_MAP independently from the
saved X/Y fields using the production normalisation

    SCALE(i) = max(abs(X(i)), 1.0e-4 * VMAXF(L))

with VMAXF(L) = max(max_J |X(J,L)|, 1.0e-25), the same column-peak
floor the march residual uses.  Also provides a manufactured-array
test of the normalisation so the metric itself is verified outside
the model run.

Record layout (list-directed WRITE, one line per captured pass):
    PSD <LS> <CH> <NRF> <NPO_solve_end> <IXIT> <18 real columns>
columns 1..18:
    1  NPO of this pass
    2  OMEGA_RAY at pass entry (pre-adaptation)
    3  OMEGA_RAY at march end (same value; adaptation runs later)
    4  production local residual RES
    5  R_MAP   = max |Y-X|/SCALE
    6  RMS_MAP = sqrt(sum((Y-X)/SCALE)^2 / NVALID)
    7  max absolute update |Y-X|
    8  J of the R_MAP argmax (ray point; 0 if none)
    9  L of the R_MAP argmax (frequency index)
    10 NBAD_CUR (cells above KN_RAY_TOL)
    11 NTOT_CUR (cells above the significance floor)
    12 PSD_NV (finite-delta cells entering the RMS)
    13 PSD_NFLR (cells with |X| < 1e-4*VMAXF)
    14 PSD_NNF (non-finite deltas excluded from RMS)
    15 NCP_P (cap events this pass)
    16 NNG_P (negative-floor events this pass)
    17 GAM sign at argmax (+1/-1, -9 unavailable)
    18 NFM_L(Jmax) (match-freq count at argmax; -9 unavailable)

Snapshot layout (unformatted sequential WRITE, one file per pass):
    rec1: NI, NCF, LS, NPO, INWARD(L4)
    rec2: X(J=1..NI, L=1..NCF)  real*8, Fortran column order
    rec3: Y same shape
"""

import glob
import struct
import sys
from dataclasses import dataclass, field

import numpy as np

REC_COLS = 18


@dataclass
class PsdRecord:
    ls: int = -1
    ch: str = "?"
    nrf: int = -1
    npo_final: int = -1
    ixit: int = -9
    cols: np.ndarray = field(default_factory=lambda: np.zeros(REC_COLS))

    @property
    def npo(self):
        return int(round(self.cols[0]))

    @property
    def omega_in(self):
        return self.cols[1]

    @property
    def res(self):
        return self.cols[3]

    @property
    def r_map(self):
        return self.cols[4]

    @property
    def rms_map(self):
        return self.cols[5]

    @property
    def dabs(self):
        return self.cols[6]

    @property
    def jmx(self):
        return int(round(self.cols[7]))

    @property
    def lmx(self):
        return int(round(self.cols[8]))


def parse_records(path):
    """Parse PSD records; list-directed output wraps each record over
    several lines, so accumulate tokens until REC_COLS floats arrive."""
    recs = []
    cur = None
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) >= 6 and tok[0] == "PSD":
                if cur is not None and len(cur.cols) == REC_COLS:
                    recs.append(cur)
                cur = PsdRecord(
                    ls=int(tok[1]),
                    ch=tok[2],
                    nrf=int(tok[3]),
                    npo_final=int(tok[4]),
                    ixit=int(tok[5]),
                )
                cur.cols = np.array(
                    [float(t.replace("D", "E").replace("d", "e")) for t in tok[6:]]
                )
            elif cur is not None:
                try:
                    extra = np.array(
                        [float(t.replace("D", "E").replace("d", "e")) for t in tok]
                    )
                except ValueError:
                    continue
                cur.cols = np.concatenate([cur.cols, extra])
    if cur is not None and len(cur.cols) == REC_COLS:
        recs.append(cur)
    return [r for r in recs if len(r.cols) == REC_COLS]


def _read_fortran_record(fh):
    (n,) = struct.unpack("<i", fh.read(4))
    payload = fh.read(n)
    (m,) = struct.unpack("<i", fh.read(4))
    assert m == n, "fortran record markers disagree"
    return payload


def read_snapshot(path):
    with open(path, "rb") as fh:
        hdr = _read_fortran_record(fh)
        ni, ncf, ls, npo, inward_i = struct.unpack("<4i?", hdr[:17])
        x = np.frombuffer(_read_fortran_record(fh), dtype="<f8").reshape(
            (ni, ncf), order="F"
        )
        y = np.frombuffer(_read_fortran_record(fh), dtype="<f8").reshape(
            (ni, ncf), order="F"
        )
    return ni, ncf, ls, npo, bool(inward_i), x, y


def map_residuals(x, y, floor_frac=1.0e-4):
    """R_MAP/RMS_MAP with the production SCALE = max(|X|, 1e-4*VMAXF)."""
    vmaxf = np.maximum(np.abs(x).max(axis=0), 1.0e-25)
    scale = np.maximum(np.abs(x), floor_frac * vmaxf[None, :])
    delta = y - x
    finite = np.isfinite(delta)
    nrm = np.zeros_like(delta)
    nrm[finite] = np.abs(delta[finite]) / scale[finite]
    nvalid = int(finite.sum())
    jmx, lmx = np.unravel_index(np.argmax(nrm), nrm.shape)
    return dict(
        r_map=float(nrm.max()),
        rms_map=float(np.sqrt((nrm[finite] ** 2).sum() / max(nvalid, 1))),
        dabs=float(np.abs(delta[finite]).max()) if nvalid else np.nan,
        jmx=int(jmx) + 1,
        lmx=int(lmx) + 1,
        nvalid=nvalid,
        nnf=int((~finite).sum()),
        nflr=int((np.abs(x) < floor_frac * vmaxf[None, :]).sum()),
    )


def check_against_records(rec_path, snap_glob="KN_PSWEEP_XY_r*_p*.bin", rtol=1e-12):
    recs = parse_records(rec_path)
    key = {(r.ls, r.ch, r.npo): r for r in recs}
    results = []
    for path in sorted(glob.glob(snap_glob)):
        ni, ncf, ls, npo, inward, x, y = read_snapshot(path)
        ch = "M" if inward else "P"
        got = map_residuals(x, y)
        rec = key.get((ls, ch, npo))
        row = {"file": path, "ls": ls, "ch": ch, "npo": npo, **got}
        if rec is not None:
            row["fort_rmap"] = rec.r_map
            row["fort_rms"] = rec.rms_map
            row["rmap_ok"] = np.isclose(rec.r_map, got["r_map"], rtol=rtol)
            row["rms_ok"] = np.isclose(rec.rms_map, got["rms_map"], rtol=rtol)
            row["jmx_ok"] = rec.jmx == got["jmx"]
            row["lmx_ok"] = rec.lmx == got["lmx"]
        results.append(row)
    return results


def manufactured_tests():
    """Normalisation checks on synthetic X/Y (Part F3)."""
    rng = np.random.default_rng(7)
    out = []
    x = rng.uniform(1e-3, 10.0, size=(6, 9))
    y = x.copy()
    y[2, 5] += 0.1 * x[2, 5]
    y[0, 0] *= 1.5
    got = map_residuals(x, y)
    scale = np.maximum(np.abs(x), 1e-4 * np.abs(x).max(axis=0)[None, :])
    exp = np.abs(y - x) / scale
    out.append(("r_map matches hand calc", np.isclose(got["r_map"], exp.max())))
    out.append(
        (
            "rms matches hand calc",
            np.isclose(got["rms_map"], np.sqrt((exp**2).mean())),
        )
    )
    out.append(("argmax J", got["jmx"] == 1))
    out.append(("argmax L", got["lmx"] == 1))
    x2 = x.copy()
    y2 = x.copy()
    y2[0, 0] = np.inf
    got2 = map_residuals(x2, y2)
    out.append(("inf delta excluded from RMS", np.isfinite(got2["rms_map"])))
    out.append(("nnf counts the inf", got2["nnf"] == 1))
    x3 = np.full((4, 3), 1e-30)
    y3 = x3.copy()
    got3 = map_residuals(x3, y3)
    out.append(("floor cells: all below 1e-4 peak", got3["nflr"] == 12))
    out.append(("zero delta gives r_map 0", got3["r_map"] == 0.0))
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--mfg":
        for name, ok in manufactured_tests():
            print(f"{'PASS' if ok else 'FAIL'}: {name}")
        sys.exit(0 if all(ok for _, ok in manufactured_tests()) else 1)
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    rows = check_against_records(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "KN_PSWEEP_XY_r*_p*.bin")
    for row in rows:
        flag = "ok" if row.get("rmap_ok") and row.get("rms_ok") else "CHECK"
        print(
            f"{flag} {row['file']} r={row['ls']} ch={row['ch']} p={row['npo']} "
            f"py r_map={row['r_map']:.6e} rms={row['rms_map']:.6e} "
            f"fort={row.get('fort_rmap', float('nan')):.6e}/"
            f"{row.get('fort_rms', float('nan')):.6e} "
            f"jmx_ok={row.get('jmx_ok')} lmx_ok={row.get('lmx_ok')}"
        )
