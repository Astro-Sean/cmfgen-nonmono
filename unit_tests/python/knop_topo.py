#!/usr/bin/env python3
"""Reader for the KNOP_TOPO compact topology checkpoint (P1-T1h).

KNOP_TOPO is written by DUMP_KNOP_TOPO_V1 in
plane/fg_j_cmf_nonmono_v1.f at the serial point just before the
OpenMP ray loop of FG_KNOP_REFRESH_V1.  It is a *topology* record,
not a portable offline operator checkpoint: it carries the
geometry/ray/frequency/dimension/configuration state needed to
rebuild the fixed-opacity transport dependency graph, plus input
store fingerprints.  The ~300 MB opacity packs are deliberately
absent (they set Jacobian weights, not the sparsity pattern).

Format: unformatted sequential, one record per WRITE, framed as
[int32 reclen][payload][int32 reclen] little-endian.  Record order
mirrors the writer exactly; arrays are stored in Fortran
column-major order and are reshaped with order='F'.
"""

import struct
from pathlib import Path

import numpy as np

MAGIC = 20261120
VER = 1
KTOPO = 1

# Integer-array presence bits (APRS_I), in fixed write order.
IARRAYS = [
    ("COMP_DEPTH_IDX", "int"),
    ("J_PNT", "int"),
    ("RAY_PNT", "int"),
    ("NI_RAY", "int"),
    ("MAX_LS", "int"),
    ("NM_RAY", "int"),
]

# Real-array presence bits (APRS_R), in fixed write order.
RARRAYS = [
    ("GAM", "real"),
    ("GAM_MAG", "real"),
    ("dGAMdR", "real"),
    ("R_RAY", "real"),
    ("Z", "real"),
    ("DLN_RAY", "real"),
    ("KN_FREQ", "real"),
    ("KN_DLN_P", "real"),
    ("KN_DLN_N", "real"),
    ("KN_DBB", "real"),
    ("KN_IC", "real"),
    ("KN_R", "real"),
    ("KN_P", "real"),
]

# ICNF slots (writer order).
ICNF_NAMES = [
    "KN_RAY_MAXPASS", "KN_MAX_THREADS", "KN_RAY_MAX_INS",
    "FG_ORD_AA_M", "FG_ORD_P2N", "FG_ORD_P2REL",
    "KN_RAY_HOTPASSES", "KN_HOT_TOTAL_PASSES", "FREQ_IDX_ITER",
    "ND_ADD", "N_TH",
]

# RCF slots (writer order), all LDP reals.
RCF_NAMES = [
    "INEXT_RELAX", "KN_RAY_TOL", "KN_RAY_BULK_TOL", "KN_RAY_OMEGA",
    "KN_RAY_STALL", "KN_RF_TOT_MB", "KN_RAY_PACK_MB",
    "FG_COMMIT_JUMP", "FG_COMMIT_NC", "FG_DIV_HEAL_CAP",
    "NONMONO_INEXT_CAP",
]

# LCFG slots (writer order), packed 0/1 logicals.
LCFG_NAMES = [
    "DO_KNOP_REFRESH", "FREEZE_I_NEXT", "FREEZE_I_PREV",
    "HAS_COMPRESSING", "KN_STORE_READY", "KN_ADAPT_RELAX",
    "FG_FREQ_MATCH", "FG_ORD_SWEEP", "FG_ORD_AITKEN", "FG_ORD_AA",
    "FG_ORD_P2DAMP", "FG_COMMIT_CONV", "FG_SEED_COMMIT",
    "FG_BULK_CONTR", "FG_PLATEAU_DET", "KN_INSERT_SAV",
    "PRINT_KNOP_DEBUG",
]


def _records(path):
    """Yield raw payload bytes for each unformatted-sequential record."""
    data = Path(path).read_bytes()
    off = 0
    n = len(data)
    while off + 8 <= n:
        (rlen,) = struct.unpack_from("<i", data, off)
        if rlen < 0 or off + 8 + rlen > n:
            raise ValueError(
                f"{path}: bad record length {rlen} at offset {off}")
        yield data[off + 4:off + 4 + rlen]
        off += 8 + rlen
    if off != n:
        raise ValueError(f"{path}: trailing {n - off} bytes")


def _i32(payload):
    return np.frombuffer(payload, dtype="<i4").copy()


def _f64(payload):
    return np.frombuffer(payload, dtype="<f8").copy()


def _chars(payload):
    return payload.decode("ascii", "replace").rstrip("\x00 ")


def read_topo(path):
    """Parse KNOP_TOPO into a dict of arrays/config.

    Returns a dict with keys: magic, ver, kind, kind_params,
    dims (dict), icnf (dict), rcf (dict), lcfg (dict), chars (dict),
    arrays (dict of ndarray), fpw (float64[10]).
    """
    recs = list(_records(path))
    pos = 0

    def nxt():
        nonlocal pos
        r = recs[pos]
        pos += 1
        return r

    # 0: magic/version/kind
    hdr = _i32(nxt())
    if hdr.size != 3 or int(hdr[0]) != MAGIC:
        raise ValueError(f"{path}: bad magic {hdr}")
    magic, ver, kind = int(hdr[0]), int(hdr[1]), int(hdr[2])
    if ver != VER or kind != KTOPO:
        raise ValueError(
            f"{path}: unsupported ver={ver} kind={kind}")

    # 1: kind params (LDP kind, default real kind, default int kind)
    kp = _i32(nxt())
    kind_params = {"ldp": int(kp[0]), "real": int(kp[1]),
                   "int": int(kp[2])}

    # 2: dims
    d = _i32(nxt())
    dims = dict(zip(
        ["KN_ND", "KN_NP", "KN_NC", "KN_NCF", "KN_ND_EXT",
         "N_COMP", "NRAY_MAX"], d.astype(int)))

    # 3: integer config
    ic = _i32(nxt())
    icnf = {n: int(v) for n, v in zip(ICNF_NAMES, ic)}

    # 4: real config
    rc = _f64(nxt())
    rcf = {n: float(v) for n, v in zip(RCF_NAMES, rc)}

    # 5: logical config
    lc = _i32(nxt())
    lcfg = {n: bool(v) for n, v in zip(LCFG_NAMES, lc)}

    # 6-8: character config
    chars = {
        "SOLUTION_METHOD": _chars(nxt()),
        "KN_INNER_METH": _chars(nxt()),
        "KN_METHOD_SAV": _chars(nxt()),
    }

    # 9: presence masks
    pm = _i32(nxt())
    aprs_i, aprs_r = int(pm[0]), int(pm[1])

    # array shapes keyed by name -> (dims tuple)
    ND, NP = dims["KN_ND"], dims["KN_NP"]
    NRMAX, NCF = dims["NRAY_MAX"], dims["KN_NCF"]
    shape_of = {
        "COMP_DEPTH_IDX": (ND,),
        "J_PNT": (ND, NP),
        "RAY_PNT": (NRMAX, NP),
        "NI_RAY": (NP,),
        "MAX_LS": (NRMAX,),
        "NM_RAY": (NRMAX,),
        "GAM": (NRMAX, NP),
        "GAM_MAG": (NRMAX, NP),
        "dGAMdR": (NRMAX, NP),
        "R_RAY": (NRMAX, NP),
        "Z": (NRMAX, NP),
        "DLN_RAY": (NRMAX,),
        "KN_FREQ": (NCF,),
        "KN_DLN_P": (NCF,),
        "KN_DLN_N": (NCF,),
        "KN_DBB": (NCF,),
        "KN_IC": (NCF,),
        "KN_R": (ND,),
        "KN_P": (NP,),
    }

    arrays = {}
    for bit, (name, kind_) in enumerate(IARRAYS):
        if aprs_i & (1 << bit):
            raw = _i32(nxt())
            arrays[name] = raw.reshape(shape_of[name], order="F")
    for bit, (name, kind_) in enumerate(RARRAYS):
        if aprs_r & (1 << bit):
            raw = _f64(nxt())
            arrays[name] = raw.reshape(shape_of[name], order="F")

    # fingerprints
    fpw = _f64(nxt())

    # end marker
    end = _i32(nxt())
    if end.size != 1 or int(end[0]) != -MAGIC:
        raise ValueError(f"{path}: missing end marker")
    if pos != len(recs):
        raise ValueError(
            f"{path}: {len(recs) - pos} trailing records")

    return {
        "magic": magic, "ver": ver, "kind": kind,
        "kind_params": kind_params, "dims": dims,
        "icnf": icnf, "rcf": rcf, "lcfg": lcfg, "chars": chars,
        "aprs_i": aprs_i, "aprs_r": aprs_r,
        "arrays": arrays, "fpw": fpw,
    }


def store_fingerprint(topo):
    """Return (P_fp, M_fp) dicts of the 5 fingerprint moments each.

    fpw layout: [P_sum,P_suma,P_sumsq,P_max,P_nnz,
                 M_sum,M_suma,M_sumsq,M_max,M_nnz].
    """
    fpw = topo["fpw"]
    keys = ["sum", "suma", "sumsq", "max", "nnz"]
    return ({k: float(fpw[i]) for i, k in enumerate(keys)},
            {k: float(fpw[5 + i]) for i, k in enumerate(keys)})


if __name__ == "__main__":
    import sys
    t = read_topo(sys.argv[1] if len(sys.argv) > 1 else "KNOP_TOPO")
    print("dims:", t["dims"])
    print("icnf:", t["icnf"])
    print("chars:", t["chars"])
    print("arrays:")
    for k, v in t["arrays"].items():
        print(f"  {k}: shape={v.shape} dtype={v.dtype} "
              f"min={v.min():.4g} max={v.max():.4g}")
    print("fpw:", t["fpw"])
