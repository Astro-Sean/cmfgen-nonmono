#!/usr/bin/env python3
"""Regenerate CMFGEN POP<species> files from CUR_MODEL_DATA.

CUR_MODEL_DATA is written by WRITE_SEQ_TIME_FILE_V1 at the end of a
CMFGEN run and carries the converged per-ion level populations.  The
OBS executable (cmf_flux) needs the ASCII POP<species> files; when a
directory refresh loses them this script reconstructs them.

POP layout (verified against op_spec_file_v2.f / rd_pop_dim.f /
rite_asc.f and a surviving POPHYD):

    Output format date:          27-JAN-1992
    Completion of Model:         <20-char time>
    ND:                          <nd>
    <SPEC>/He abundance:         <number fraction>
    <nd floats>                  POP_SPECIES (species density per depth)
    per ion (ascending stage):
      Number of <ION_ID> levels: <nlev>
      <ION_ID> Oscillator date:  <date>
      <nlev*nd floats>           XzV_F, Fortran (nlev,nd) order
      <nd floats>                DXzV (superlevel population)

Usage: cur_model_data_to_pop.py <cmfgen_dir> [out_dir]
"""

import sys
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

GEN_ION_ID = [
    "I", "2", "III", "IV", "V", "SIX", "SEV", "VIII", "IX", "X",
    "XI", "XII", "XIII", "XIV", "XV", "XSIX", "XSEV", "X8", "X9", "XX",
]

# (SPECIES name in CUR_MODEL_DATA) -> (SPECIES_ABR for ION_ID, AT_MASS)
# Mirrors the table in obs/cmf_flux_v5.f.
SPECIES_INFO = {
    "HYD": ("H", 1.0), "HE": ("He", 4.0), "CARB": ("C", 12.0),
    "NIT": ("N", 14.0), "OXY": ("O", 16.0), "FLU": ("F", 19.0),
    "NEON": ("Ne", 20.2), "SOD": ("Na", 23.0), "MAG": ("Mg", 24.3),
    "ALUM": ("Al", 27.0), "SIL": ("Sk", 28.1), "PHOS": ("P", 31.0),
    "SUL": ("S", 32.1), "CHL": ("Cl", 35.5), "ARG": ("Ar", 39.9),
    "POT": ("K", 39.1), "CAL": ("Ca", 40.1), "SCAN": ("Sc", 44.96),
    "TIT": ("Tk", 47.88), "VAN": ("V", 50.94), "CHRO": ("Cr", 52.0),
    "MANG": ("Mn", 54.9), "IRON": ("Fe", 55.8), "COB": ("Co", 58.9),
    "NICK": ("Nk", 58.7), "COPP": ("Cu", 63.6), "ZINC": ("Zn", 65.4),
}


def read_cur_model_data(path):
    """Parse the unformatted sequential CUR_MODEL_DATA file."""
    f = FortranFile(str(path))
    fmt_date = f.read_record("S11").tobytes().decode().strip()
    nd, num_species = f.read_ints(np.int32)
    sn_age = f.read_reals(np.float64)[0]
    vectors = {}
    for name in ("R", "V", "SIGMA", "T", "ED", "POP_ATOM", "DENSITY"):
        vectors[name] = f.read_reals(np.float64)
        assert vectors[name].size == nd, (name, vectors[name].size, nd)
    species = []
    for _ in range(num_species):
        rec = f.read_record(np.int32, np.int32, "S10")
        beg, endm1 = int(rec[0][0]), int(rec[1][0])
        spec = rec[2].tobytes().decode().strip()
        ions = []
        for _ in range(beg, endm1 + 1):
            rec = f.read_record(np.int32, np.int32, np.float64)
            iid, nlev, zx = int(rec[0][0]), int(rec[1][0]), float(rec[2][0])
            xzv = f.read_reals(np.float64)
            assert xzv.size == nlev * nd, (spec, iid, xzv.size, nlev * nd)
            dxz = f.read_reals(np.float64)
            ions.append({"id": iid, "nlev": nlev, "zx": zx,
                         "xzv": xzv, "dxz": dxz})
        species.append({"name": spec, "beg": beg, "endm1": endm1,
                        "ions": ions})
    f.close()
    return fmt_date, nd, sn_age, vectors, species


def read_mass_fractions(model_path):
    """Negative entries in the [SPEC/X] MODEL lines are mass fractions."""
    mfrac = {}
    for line in Path(model_path).read_text(errors="replace").splitlines():
        if "/X]" not in line:
            continue
        try:
            val = float(line.split()[0])
        except ValueError:
            continue
        key = line.split("[")[1].split("/")[0].strip().upper()
        mfrac[key] = abs(val)
    return mfrac


def number_fractions(mfrac):
    """Convert mass fractions to number fractions (SET_ABUND_CLUMP)."""
    t2 = 0.0
    nfrac = {}
    for key, mf in mfrac.items():
        mass = SPECIES_INFO.get(key, (key, 1.0))[1]
        nfrac[key] = mf / mass
        t2 += mf / mass
    if t2 > 0:
        for key in nfrac:
            nfrac[key] /= t2
    return nfrac


def write_e167(vals):
    """8 values per line, 1PE16.7, matching '(1X,1P8E16.7)'."""
    vals = np.asarray(vals).ravel()
    return [
        " " + "".join(f"{v:16.7E}" for v in vals[i : i + 8])
        for i in range(0, vals.size, 8)
    ]


def gen_pop_file(spec, nd, abund, time_str):
    name = spec["name"]
    abr = SPECIES_INFO.get(name, (name[:2], 1.0))[0]
    pop_species = np.zeros(nd)
    for ion in spec["ions"]:
        # XzV_F is stored Fortran (nlev, nd) column-major: each group of
        # nlev consecutive values is one depth.
        pop_species += ion["xzv"].reshape(nd, ion["nlev"]).sum(axis=1)
    if spec["ions"]:
        pop_species += spec["ions"][-1]["dxz"]

    lines = [
        f" {'Output format date:':<28s}27-JAN-1992",
        f" {'Completion of Model:':<28s}{time_str:<20s}",
        f" {'ND:':<28s}{nd:d}",
        f" {name + '/He abundance:':<28s}{abund:12.5E}",
    ]
    lines += write_e167(pop_species)
    for ion in spec["ions"]:
        z = int(round(ion["zx"]))
        ion_id = abr + GEN_ION_ID[z - 1]
        lines.append(
            f" {'Number of ' + ion_id + ' levels:':<28s}{ion['nlev']:d}"
        )
        lines.append(f" {ion_id + ' Oscillator date:':<28s}recovered")
        # Raw record order is already (nlev, nd) Fortran order = what
        # RITE_ASC writes.
        lines += write_e167(ion["xzv"])
        lines += write_e167(ion["dxz"])
    return "\n".join(lines) + "\n", pop_species


def main():
    d = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else d
    out_dir.mkdir(parents=True, exist_ok=True)
    fmt_date, nd, sn_age, vectors, species = read_cur_model_data(
        d / "CUR_MODEL_DATA"
    )
    print(f"format={fmt_date} ND={nd} species={len(species)} age={sn_age}")
    nfrac = number_fractions(read_mass_fractions(d / "MODEL"))
    time_str = "recovered"
    rvtj = d / "RVTJ"
    if rvtj.exists():
        for line in rvtj.read_text(errors="replace").splitlines()[:5]:
            if "Completion of Model:" in line:
                time_str = line.split(":", 1)[1].strip()[:20]
    for spec in species:
        abund = nfrac.get(spec["name"].upper(), 0.0)
        text, pop_species = gen_pop_file(spec, nd, abund, time_str)
        if pop_species[-1] == 0.0:
            continue
        out = out_dir / ("POP" + spec["name"].strip())
        out.write_text(text)
        print(
            f"  wrote {out.name} "
            f"({len(spec['ions'])} ions, abund={abund:.5E})"
        )


if __name__ == "__main__":
    main()
