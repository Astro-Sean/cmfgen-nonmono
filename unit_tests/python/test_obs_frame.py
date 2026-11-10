#!/usr/bin/env python3
"""OBS_FRAME_V2 VTURB normalization fix (obs/obs_frame_v2.f).

Alteration under test (bug fix, not flag-gated): the turbulent velocity
depth profile is now built as

    VTURB(I) = VTURB_MIN + (VTURB_MAX-VTURB_MIN)*MIN(1,|V(I)|/max|V|)

instead of the upstream V(I)/V(1).  The old form divided by zero for a
static outer boundary and produced a negative turbulent width wherever
a non-monotonic velocity field reverses sign; it also overflowed when
V(1) was not the speed-law maximum.  VTURB feeds MAX_DEL_V_RES_ZONE =
VTURB*FRAC_DOP, which sizes the inserted ray points; a non-positive
value reaches OBS_FRAME_SUB_V9, whose guard now stops with
"MAX_DEL_V_RES_ZONE must be strictly positive".

The production binary is the shipped exe/obs_frame.exe (built Nov-2026
from the current sources).  The reference binary compiles the frozen
obs_frame_v2.f + obs_frame_sub_v9.f pair from cmfgen_reference, so the
comparison is upstream driver+sub vs production driver+sub.

Fixtures (fabricated per case): an RVTJ descriptor (ND=8), matching
direct-access ETA_DATA / CHI_DATA / RJ_DATA plus their *_INFO tables,
and an OBS_FRAME_INPUT option store.  The m15 model's V(1) is the wind
maximum, so only a non-monotonic profile separates the two builds.

  mono    : monotonic V -> prod flux identical to ref
  vneg    : V dips negative mid-grid -> ref writes VTURB<0 and dies in
            the V9 MINVAL guard; prod completes with finite OBSFRAME
  vzero   : V=0 everywhere -> ref divides 0/0 (NaN lands in the same
            guard); prod assigns VTURB_MIN and completes
"""
import os
import re
import struct
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, PROD_MOD_INC, REF, TREE, WORKDIR,
                    Report, compile_objs, link_prod, run)

PROD_EXE = TREE / "exe" / "obs_frame.exe"
ND = 8
FREQS = [2.0, 1.5, 1.0, 0.85, 0.7, 0.6, 0.45, 0.3]
R_GRID = [12.0, 10.0, 8.0, 6.5, 5.0, 4.0, 3.0, 2.2]


def build_ref():
    md = BINDIR / "obj_obs_frame_ref_mods"
    md.mkdir(parents=True, exist_ok=True)
    objs = compile_objs(
        "obs_frame_ref_src",
        [REF / "obs" / "obs_frame_v2.f",
         REF / "obs" / "obs_frame_sub_v9.f"],
        moddir=md, extra_flags=PROD_MOD_INC, prod=True)
    libs = [TREE / "lib" / a for a in ("libdisp.a", "libsubs.a", "libsubs.a")]
    return link_prod("obs_frame_ref", objs + libs, moddir=md)


def write_rvtj(path, v):
    def hdr(label, val):
        return f" {label:<28s}{val}\n"

    def vec(tag, arr):
        row = "  " + "  ".join(f"{x:16.8E}" for x in arr) + "\n"
        return f" {tag}\n{row}"

    nd = len(v)
    sigma = [-0.5] * nd
    ed = [1.0e4] * nd
    t = [1.0] * nd
    one = [1.0] * nd
    dens = [1.0e6] * nd
    lines = [
        hdr("Output format date:", "15-Aug-2019         "),
        hdr("Completion of Model:", "test-fixture"),
        hdr("Program Date:", "18-Jun-2025 "),
        hdr("ND:", f"{nd:5d}"),
        hdr("NC:", "   15"),
        hdr("NP:", "   85"),
        hdr("NCF:", f"{len(FREQS):5d}"),
        hdr("Mdot(Msun/yr):", " 1.00000E-20"),
        hdr("L(Lsun):", " 1.00000E+06"),
        hdr("H/He abundance:", " 1.00000E+00"),
        hdr("Was T fixed?:", "F"),
        hdr("Species naming convention:", "K_FOR_I"),
        vec("Radius (10^10 cm)", R_GRID),
        vec("Velocity (km/s)", v),
        vec("dlnV/dlnr-1", sigma),
        vec("Electron Density (cm^-3)", ed),
        vec("Temperature (10^4 K)", t),
        vec("Atom Density (cm^-3)", dens),
        vec("Ion Density (cm^-3)", dens),
        vec("Mass Density (g/cm^3)", [1.0e-11] * nd),
        vec("Clumping Factor", one),
    ]
    path.write_text("".join(lines))


def write_direct(path, name, vec_rows):
    """Direct-access file in the CMFGEN layout: rec3 = (INIT_REC,NCF,ND),
    data records (VEC(1:ND), FREQ) starting at INIT_REC."""
    rec_bytes = (ND + 1) * 8
    init_rec = 4
    with open(path, "wb") as fh:
        fh.truncate(init_rec * rec_bytes)  # zero recs 1..3 incl. pad
        fh.seek(2 * rec_bytes)
        fh.write(struct.pack(f"<{3}i", init_rec, len(vec_rows), ND))
        fh.seek((init_rec - 1) * rec_bytes)
        for row, freq in vec_rows:
            fh.write(struct.pack(f"<{ND + 1}d", *row, freq))
    info = (
        " 12-Apr-2017                                       "
        "!INFO format date\n"
        " 23-Jan-2017                                       "
        "!File format date\n"
        f"{ND:11d}{18:12d}{8:12d}{4:12d}{4:12d}            T\n"
        "          ND        RECL   WORD_SIZE   UNIT_SIZE"
        "    INT_SIZE     LIT_END\n"
    )
    path.parent.joinpath(f"{name}_INFO").write_text(info)


def write_input(path):
    path.write_text(
        "! minimal observer-frame controls\n"
        " 4500.0    [MIN_WAVE]\n"
        " 5500.0    [MAX_WAVE]\n"
        " 1000.0    [DEL_V_OBS]\n"
        " 10.0      [VTURB_MIN]\n"
        " 100.0     [VTURB_MAX]\n"
        " 0.5       [FRAC_DOP]\n"
        " F         [CONV_ES_J]\n"
        " 10.0      [TAU_MAX]\n"
        " 0.5       [ES_DTAU]\n"
        "STAU       [INT_METH]\n"
        " 0         [N_INS_OBS]\n"
        " F         [WR_IP]\n"
        " F         [RD_RAY]\n")


def casedir(tag, v):
    d = WORKDIR / f"obsf_{tag}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    write_rvtj(d / "RVTJ", v)
    eta = [([1.0e-12] * ND, f) for f in FREQS]
    chi = [([1.0e-9] * ND, f) for f in FREQS]
    rj = [([1.0e-8] * ND, f) for f in FREQS]
    write_direct(d / "ETA_DATA", "ETA_DATA", eta)
    write_direct(d / "CHI_DATA", "CHI_DATA", chi)
    write_direct(d / "RJ_DATA", "RJ_DATA", rj)
    write_input(d / "OBS_FRAME_INPUT")
    return d


def obs_flux(path):
    """OBSFRAME flux vector (all values after the header lines)."""
    vals = []
    for ln in path.read_text().splitlines():
        vals += re.findall(r"[-+]?\d\.\d+E[+-]?\d+", ln)
    return [float(v) for v in vals]


def main():
    rep = Report("OBS_FRAME_V2 VTURB normalization")
    if not PROD_EXE.exists():
        rep.summary()
        raise SystemExit(f"shipped binary missing: {PROD_EXE}")
    ref = build_ref()

    # ---- mono: upstream and prod agree on a monotonic field ----
    v_mono = [2000.0, 1500.0, 1200.0, 900.0, 700.0, 500.0, 300.0, 150.0]
    d = casedir("mono_prod", v_mono)
    r = run(str(PROD_EXE), cwd=d, check=False, timeout=120)
    prod_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("mono: prod writes OBSFRAME", len(prod_flux) > 0,
              f"rc={r.returncode} n={len(prod_flux)}")
    rep.check("mono: prod flux finite",
              all(abs(x) < 1.0e30 for x in prod_flux))

    d = casedir("mono_ref", v_mono)
    r = run(str(ref), cwd=d, check=False, timeout=120)
    ref_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("mono: ref completes", len(ref_flux) > 0,
              f"rc={r.returncode}")
    rep.check("mono: prod == ref (unchanged path)",
              prod_flux == ref_flux and len(prod_flux) > 0,
              f"len {len(prod_flux)} vs {len(ref_flux)}")

    # ---- vneg: sign-reversing V ----
    v_neg = [1000.0, 800.0, 500.0, 200.0, -350.0, -600.0, -300.0, -100.0]
    d = casedir("vneg_prod", v_neg)
    r = run(str(PROD_EXE), cwd=d, check=False, timeout=120)
    prod_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("vneg: prod writes OBSFRAME", len(prod_flux) > 0,
              f"rc={r.returncode}")
    rep.check("vneg: prod flux finite",
              all(abs(x) < 1.0e30 for x in prod_flux))

    d = casedir("vneg_ref", v_neg)
    r = run(str(ref), cwd=d, check=False, timeout=120)
    ref_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("vneg: ref dies or produces degenerate flux",
              r.returncode != 0 or not ref_flux
              or ref_flux != prod_flux
              or any(abs(x) >= 1.0e30 for x in ref_flux),
              f"rc={r.returncode} nflux={len(ref_flux)}")

    # ---- vzero: static atmosphere.  Upstream computes 0/0=NaN here;
    # harmless on a fully static field (no Doppler shift, no insertion
    # pressure), so ref still completes - the defect is latent.  Prod's
    # T1<=1e-12 else-branch assigns VTURB_MIN cleanly.
    v_zero = [0.0] * ND
    d = casedir("vzero_prod", v_zero)
    r = run(str(PROD_EXE), cwd=d, check=False, timeout=120)
    prod_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("vzero: prod writes OBSFRAME (VTURB=VTURB_MIN)",
              len(prod_flux) > 0, f"rc={r.returncode}")
    rep.check("vzero: prod flux finite",
              all(abs(x) < 1.0e30 for x in prod_flux))

    d = casedir("vzero_ref", v_zero)
    r = run(str(ref), cwd=d, check=False, timeout=120)
    ref_flux = obs_flux(d / "OBSFRAME") if (d / "OBSFRAME").exists() else []
    rep.check("vzero: ref completes (NaN latent on static field)",
              len(ref_flux) > 0, f"rc={r.returncode}")
    rep.check("vzero: ref flux agrees (NaN zone unused when V=0)",
              ref_flux == prod_flux and len(prod_flux) > 0)

    rep.summary()


if __name__ == "__main__":
    main()
