#!/usr/bin/env python3
"""Independent KNRDUMP v2 deck parser (magic 20261123).

Reads Fortran unformatted sequential records and validates the
record layout against docs/audits/knop_phase226/v2_format_spec.md,
which was derived from the writer source (fg_j_cmf_nonmono_v1.f
lines 11205-11284). Deliberately shares no code with the Fortran
replay driver so a layout bug in one cannot hide in the other.
"""
import hashlib
import json
import struct
import sys

import numpy as np

R8 = np.dtype("<f8")
I4 = np.dtype("<i4")


def read_records(path):
    """Yield all unformatted record payloads as bytes."""
    with open(path, "rb") as fh:
        while True:
            head = fh.read(4)
            if len(head) == 0:
                return
            if len(head) < 4:
                raise ValueError("truncated record marker")
            n = struct.unpack("<i", head)[0]
            if n < 0:
                raise ValueError(f"negative record length {n}")
            body = fh.read(n)
            tail = fh.read(4)
            if len(body) < n or len(tail) < 4:
                raise ValueError(f"truncated record body len={n}")
            if struct.unpack("<i", tail)[0] != n:
                raise ValueError("record marker mismatch")
            yield body


def parse_block(recs, idx):
    """Parse one refresh block starting at recs[idx]."""
    h = np.frombuffer(recs[idx], I4)
    if h[0] != 20261123:
        raise ValueError(f"block {idx}: magic {h[0]} != 20261123")
    (magic, ls, ni, kn_nd, nxd, ncf, kn_nc, n_comp, n_add,
     ncps, imeth, maxpass, cont_np, fbits, kn_np) = h[:15]
    i = idx + 1

    def arr(n, dt=R8):
        nonlocal i
        a = np.frombuffer(recs[i], dt).copy()
        assert a.size == n, f"rec {i}: want {n} got {a.size}"
        i += 1
        return a

    tol = arr(13)
    inner_meth = recs[i].decode(errors="replace").strip()
    i += 1
    blk = {"ls": int(ls), "ni": int(ni), "kn_nd": int(kn_nd),
           "nxd": int(nxd), "ncf": int(ncf), "kn_nc": int(kn_nc),
           "n_comp": int(n_comp), "fbits": int(fbits),
           "imeth": int(imeth), "inner_meth": inner_meth,
           "kn_ray_tol": float(tol[0])}
    blk["r_ext"] = arr(nxd)
    blk["log_r_ext"] = arr(nxd)
    blk["kn_r"] = arr(kn_nd)
    blk["kn_p"] = arr(1)[0]
    blk["z"] = arr(ni)
    blk["r_ray"] = arr(ni)
    blk["ray_pnt"] = arr(ni, I4)
    blk["j_pnt"] = arr(kn_nd, I4)
    blk["comp_depth_idx"] = arr(kn_nd, I4)
    blk["gam"] = arr(ni)
    blk["gam_mag"] = arr(ni)
    blk["dgamdr"] = arr(ni)
    blk["chi_coef"] = arr(nxd * 4 * ncf).reshape(nxd, 4, ncf)
    blk["eta_coef"] = arr(nxd * 4 * ncf).reshape(nxd, 4, ncf)
    if fbits & (1 << 9):
        blk["chi_ext"] = arr(nxd * ncf, np.dtype("<f4")).reshape(
            nxd, ncf)
        blk["eta_ext"] = arr(nxd * ncf, np.dtype("<f4")).reshape(
            nxd, ncf)
    blk["kn_freq"] = arr(ncf)
    blk["dln_p"] = arr(ncf)
    blk["dln_n"] = arr(ncf)
    blk["kn_ic"] = arr(ncf)
    blk["kn_dbb"] = arr(ncf)
    blk["fld_m"] = arr(ni * ncf).reshape(ni, ncf)
    blk["i_m_next"] = arr(n_comp * (ncf - 1)).reshape(
        n_comp, ncf - 1)
    blk["i_p_next"] = arr(n_comp * (ncf - 1)).reshape(
        n_comp, ncf - 1)
    if fbits & (1 << 14):
        blk["jmp_sib_m"] = arr(ncf)
        blk["jmp_sib_p"] = arr(ncf)
    blk["fld_p_seed"] = arr(ni * ncf).reshape(ni, ncf)

    st = np.frombuffer(
        recs[i],
        np.dtype([("res", "<f8"), ("np", "<i4"), ("ix", "<i4"),
                  ("om", "<f8")]))
    assert st.itemsize == len(recs[i]), "status rec size"
    blk["res_p"], blk["np_p"], blk["ixit_p"], blk["omega"] = (
        float(st["res"][0]), int(st["np"][0]), int(st["ix"][0]),
        float(st["om"][0]))
    i += 1
    blk["fld_p_end"] = arr(ni * ncf).reshape(ni, ncf)
    return blk, i


def main(path):
    recs = list(read_records(path))
    out = {"file": path,
           "sha256": hashlib.sha256(open(path, "rb").read())
           .hexdigest(),
           "n_records": len(recs), "blocks": []}
    idx = 0
    while idx < len(recs):
        blk, idx = parse_block(recs, idx)
        out["blocks"].append({
            "ls": blk["ls"], "ni": blk["ni"], "ncf": blk["ncf"],
            "fbits": blk["fbits"], "imeth": blk["imeth"],
            "inner_meth": blk["inner_meth"],
            "res_p": float(blk["res_p"]), "np_p": int(blk["np_p"]),
            "ixit_p": int(blk["ixit_p"]),
            "has_ext": "chi_ext" in blk, "has_sib": "jmp_sib_m" in blk,
        })
    return out


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1]), indent=1))
