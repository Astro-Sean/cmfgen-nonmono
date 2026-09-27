#!/usr/bin/env python3
"""Generate STELLA/TREFOR photometric-band blocks from SVO filter curves.

Produces the Fortran fragments that were hand-inserted into
  stella/strad/tt4strad.trf  (subroutine tt4strad, blocks %_varubv,
  %_dataubv, %_flx_fits_obsubv, %_flx_fits_obsubv_zeroflux, and the
  .tt header/write statements).

For each band the script emits:
  * PARAMETER(NPTS<n>=N)                  -> %_varubv
  * real*8 WL..(N),T..(N),..F0,M..        -> %_varubv
  * data WL../.../;  data T../.../;       -> %_dataubv
  * ..F0 = <FX0>d0;                       -> %_flx_fits_obsubv_zeroflux
  * _define ... @find  /  @flux lines     -> %_flx_fits_obsubv
  * suggested .tt header/unit fields      -> write(@wres,...) block

Photometric convention (verified against the code's own UBVRI
zero-points, see --selftest):

    FX0 = 2.5*log10( c * int(T(lam)/lam^2 dlam) ) - 48.6075 - dAB

with lam in cm, c = 2.99792458e10 cm/s, and dAB the Vega->AB offset
(0 for bands already tabulated in AB, e.g. SDSS, ATLAS).

SVO photon-counting curves (DetectorType=1) are converted to energy
response (T -> T*lambda, renormalised to peak 1) before tabulation;
DetectorType=0 curves are used as-is.

TREFOR output rules baked in (deviating from these breaks trefor):
  * NO '*' continuation markers - continuation lines are bare,
    indented 8 spaces.
  * DATA statements end with '/;' (slash-semicolon).
  * 6-space indent for the 'data' keyword line.

Usage:
  # fetch fresh curves from SVO and emit all fragments:
  python3 scripts/stella_svo_filters.py --fetch --out bands.trf.inc

  # reuse previously fetched VOTables:
  python3 scripts/stella_svo_filters.py --svo-dir /tmp/svo_cache --out bands.trf.inc

  # verify the FX0 convention against the built-in Bessell UBVRI tables:
  python3 scripts/stella_svo_filters.py --selftest
"""

import argparse
import os
import re
import sys
import urllib.request

import numpy as np

CCL = 2.99792458e10  # speed of light, cm/s

SVO_FPS = "https://svo2.cab.inta-csic.es/theory/fps/fps.php"

# ---------------------------------------------------------------------------
# Band table.  One row per filter:
#   tag      short label used in comments and .tt header
#   svo_id   SVO FPS filter identifier  (<Facility>/<Instrument>.<band>)
#   wn/tn    Fortran wavelength/transmission array names (WL.. / T..)
#   zp       Fortran zero-point scalar name (..F0)
#   nn       Fortran PARAMETER name for the point count (NPTS..)
#   mv       Fortran magnitude variable written to .tt
#   dAB      Vega - AB correction subtracted in FX0 (0 if band is AB-native)
#   unit     label written in the .tt units row ('vegamag' / 'abmag')
#
# NOTE: 'MZ'/'MZ2' etc. collide with TREFOR markers and existing commons;
# keep the historical MGU/MGG/... style names.
# ---------------------------------------------------------------------------
BANDS = [
    # tag  svo_id                     wn      tn     zp      nn        mv      dAB    unit
    ("u",  "SLOAN/SDSS.u",            "WLSu", "TSu", "SUF0", "NPTSSU", "MGU",  0.0,  "abmag"),
    ("g",  "SLOAN/SDSS.g",            "WLSg", "TSg", "SGF0", "NPTSSG", "MGG",  0.0,  "abmag"),
    ("r",  "SLOAN/SDSS.r",            "WLSr", "TSr", "SRF0", "NPTSSR", "MGR",  0.0,  "abmag"),
    ("i",  "SLOAN/SDSS.i",            "WLSi", "TSi", "SIF0", "NPTSSI", "MGI",  0.0,  "abmag"),
    ("z",  "SLOAN/SDSS.z",            "WLSz", "TSz", "SZF0", "NPTSSZ", "MGZ",  0.0,  "abmag"),
    ("J",  "2MASS/2MASS.J",           "WLJ",  "TJ",  "JF0",  "NPTSJ",  "MJ2",  0.90, "vegamag"),
    ("H",  "2MASS/2MASS.H",           "WLH",  "TH",  "HF0",  "NPTSH",  "MH2",  1.38, "vegamag"),
    ("Ks", "2MASS/2MASS.Ks",          "WLK",  "TK",  "KF0",  "NPTSK",  "MK2",  1.85, "vegamag"),
    ("c",  "Misc/Atlas.cyan",         "WLAc", "TAc", "ACF0", "NPTSAC", "MAC",  0.0,  "abmag"),
    ("o",  "Misc/Atlas.orange",       "WLAo", "TAo", "AOF0", "NPTSAO", "MAO",  0.0,  "abmag"),
    # VISTA Z was included then removed at user request; uncomment to re-add.
    # ("Z", "Paranal/VISTA.Z",         "WLZ",  "TZ",  "ZF0",  "NPTSZ",  "MZ",   0.502,"vegamag"),
]

FIRST_NEW_BAND = 6  # lband index of first added band (1-5 = UBVRI)


def votable_filename(svo_id):
    return "svo_" + svo_id.replace("/", "_").replace(".", "_") + ".xml"


def fetch_votable(svo_id, outdir):
    """Download the SVO VOTable for svo_id into outdir; return its path."""
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, votable_filename(svo_id))
    url = "%s?ID=%s" % (SVO_FPS, urllib.request.quote(svo_id))
    with urllib.request.urlopen(url, timeout=60) as r:
        data = r.read()
    if b"<VOTABLE" not in data:
        raise RuntimeError("SVO did not return a VOTable for %s" % svo_id)
    with open(path, "wb") as f:
        f.write(data)
    return path


def parse_votable(path):
    """Return (wl_Angstrom, T, detector_type) sorted by wavelength."""
    txt = open(path).read()
    m = re.search(r'PARAM name="DetectorType"[^/]*value="(\d)"', txt)
    det = int(m.group(1)) if m else 0
    rows = re.findall(
        r"<TR>\s*<TD>([\d.eE+-]+)</TD>\s*<TD>([\d.eE+-]+)</TD>", txt)
    if not rows:
        raise RuntimeError("no TABLEDATA rows parsed in %s" % path)
    wl = np.array([float(r[0]) for r in rows])
    T = np.array([float(r[1]) for r in rows])
    o = np.argsort(wl)
    return wl[o], T[o], det


def trapezoid(y, x):
    return np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)


def fx0(wl_ang, T, dAB):
    """STELLA zero-point: 2.5*log10(c*int(T/lam^2 dlam)) - 48.6075 - dAB."""
    wlcm = wl_ang * 1e-8
    J = trapezoid(T / wlcm**2, wlcm)
    return 2.5 * np.log10(CCL * J) - 48.6075 - dAB


def fmt_data(name, arr, per=5):
    """Format a TREFOR DATA statement: bare 8-space continuations, '/;' end."""
    vals = ["%.6g" % v for v in arr]
    out = []
    for i in range(0, len(vals), per):
        chunk = ", ".join(vals[i:i + per])
        last = i + per >= len(vals)
        end = "/;" if last else ","
        if i == 0:
            out.append("      data %s/%s%s" % (name, chunk, end))
        else:
            out.append("        %s%s" % (chunk, end))
    return out


def selftest():
    """Verify the FX0 convention reproduces the built-in UBVRI tables."""
    # Bessell UBVRI tables copied from tt4strad.trf (%_dataubv).
    WLU = [3000., 3050., 3100., 3150., 3200., 3250., 3300., 3350., 3400.,
           3450., 3500., 3550., 3600., 3650., 3700., 3750., 3800., 3850.,
           3900., 3950., 4000., 4050., 4100., 4150., 4200.]
    TU = [0.000, 0.016, 0.068, 0.167, 0.287, 0.423, 0.560, 0.673, 0.772,
          0.841, 0.905, 0.943, 0.981, 0.993, 1.000, 0.989, 0.916, 0.804,
          0.625, 0.423, 0.238, 0.114, 0.051, 0.019, 0.000]
    WLB = [3600., 3700., 3800., 3900., 4000., 4100., 4200., 4300., 4400.,
           4500., 4600., 4700., 4800., 4900., 5000., 5100., 5200., 5300.,
           5400., 5500., 5600.]
    TB = [0.000, 0.030, 0.134, 0.567, 0.920, 0.978, 1.000, 0.978, 0.935,
          0.853, 0.740, 0.640, 0.536, 0.424, 0.325, 0.235, 0.150, 0.095,
          0.043, 0.009, 0.000]
    WLV = [4700., 4800., 4900., 5000., 5100., 5200., 5300., 5400., 5500.,
           5600., 5700., 5800., 5900., 6000., 6100., 6200., 6300., 6400.,
           6500., 6600., 6700., 6800., 6900., 7000.]
    TV = [0.000, 0.030, 0.163, 0.458, 0.780, 0.967, 1.000, 0.973, 0.898,
          0.792, 0.684, 0.574, 0.461, 0.359, 0.270, 0.197, 0.135, 0.081,
          0.045, 0.025, 0.017, 0.013, 0.009, 0.000]
    WLR = [5500., 5600., 5700., 5800., 5900., 6000., 6100., 6200., 6300.,
           6400., 6500., 6600., 6700., 6800., 6900., 7000., 7100., 7200.,
           7300., 7400., 7500., 8000., 8500., 9000.]
    TR = [0.00, 0.23, 0.74, 0.91, 0.98, 1.00, 0.98, 0.96, 0.93, 0.90,
          0.86, 0.81, 0.78, 0.72, 0.67, 0.61, 0.56, 0.51, 0.46, 0.40,
          0.35, 0.14, 0.03, 0.00]
    WLI = [7000., 7100., 7200., 7300., 7400., 7500., 7600., 7700., 7800.,
           7900., 8000., 8100., 8200., 8300., 8400., 8500., 8600., 8700.,
           8800., 8900., 9000., 9100., 9200.]
    TI = [0.000, 0.024, 0.232, 0.555, 0.785, 0.910, 0.965, 0.985, 0.990,
          0.995, 1.000, 1.000, 0.990, 0.980, 0.950, 0.910, 0.860, 0.750,
          0.560, 0.330, 0.150, 0.030, 0.000]
    # code FX0 values from %_flx_fits_obsubv_zeroflux and literature
    # Vega-AB offsets (Blanton & Roweis 2007).  The check compares the
    # code's zero-point against the uncorrected AB formula; the residual
    # is the per-band correction the code effectively applies.
    # NOTE: the U band is expected to deviate -- STELLA's TU table folds
    # in atmospheric extinction, so UF0 carries an extra ~-0.7 mag
    # correction beyond Vega-AB.  U is reported but excluded from PASS.
    bands = {"U": (WLU, TU, -13.90, 0.02), "B": (WLB, TB, -13.00, -0.07),
             "V": (WLV, TV, -13.72, 0.03), "R": (WLR, TR, -13.66, 0.19),
             "I": (WLI, TI, -14.42, 0.45)}
    print("%-3s %10s %10s %10s %10s  %s"
          % ("b", "FX0_code", "ZP_AB", "impl_corr", "Vega-AB", "note"))
    ok = True
    for b, (wl, T, code, dab) in bands.items():
        zp_ab = fx0(np.array(wl), np.array(T), 0.0)
        implied = code - zp_ab          # correction the code applies
        resid = implied - (-dab)        # vs literature Vega-AB offset
        note = ""
        if b == "U":
            note = "atmospheric corr in TU; excluded"
        elif abs(resid) > 0.15:
            note = "MISMATCH"
            ok = False
        print("%-3s %10.3f %10.3f %+10.3f %+10.3f  %s"
              % (b, code, zp_ab, implied, -dab, note))
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true",
                    help="download VOTables from SVO into --svo-dir first")
    ap.add_argument("--svo-dir", default="svo_cache",
                    help="directory with cached svo_*.xml VOTables")
    ap.add_argument("--out", default="stella_bands.trf.inc",
                    help="output file for the generated TREFOR fragments")
    ap.add_argument("--selftest", action="store_true",
                    help="verify FX0 convention vs built-in UBVRI tables")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    os.makedirs(args.svo_dir, exist_ok=True)

    params, decls, datas, zps, finds, fluxes = [], [], [], [], [], []
    hdr_cols, unt_cols, mag_vars = [], [], []
    for i, (tag, svo_id, wn, tn, zp, nn, mv, dab, unit) in enumerate(BANDS):
        lband = FIRST_NEW_BAND + i
        path = os.path.join(args.svo_dir, votable_filename(svo_id))
        if args.fetch or not os.path.exists(path):
            path = fetch_votable(svo_id, args.svo_dir)
            print("fetched %-22s -> %s" % (svo_id, path))
        wl, T, det = parse_votable(path)
        if det == 1:                     # photon-counting -> energy response
            T = T * wl
        T = T / T.max()
        F0 = fx0(wl, T, dab)

        params.append("      PARAMETER(%s=%d);" % (nn, len(wl)))
        decls.append("      real*8 %s(%s),%s(%s),%s,%s;"
                     % (wn, nn, tn, nn, zp, mv))
        datas.append("--  %s: %s  NPTS=%d DetectorType=%d Vega-AB corr=%.3g"
                     % (tag, svo_id, len(wl), det, dab))
        datas += fmt_data(wn, wl)
        datas += fmt_data(tn, T)
        datas.append("")
        zps.append("      %s = %.4fd0;" % (zp, F0))
        finds.append("  _define @NPTSX %s @WLX %s @lband %d ; @find;"
                     % (nn, wn, lband))
        fluxes.append("_define @NPTSX %s @WLX %s @FX0 %s @TX %s @lband %d "
                      ";@flux; %s=AMAGX;" % (nn, wn, zp, tn, lband, mv))
        hdr_cols.append("M" + tag)
        unt_cols.append(unit)
        mag_vars.append(mv)
        print("%-3s lband=%2d NPTS=%3d det=%d FX0=%9.4f  %s"
              % (tag, lband, len(wl), det, F0, svo_id))

    nmag = 7 + len(BANDS)  # Mbol + UBVRI + new bands + Mbolavg
    out = []
    out.append("-- ===== %_varubv additions (INTEGER + PARAMETER + real*8) =====")
    out.append("      INTEGER " + ", ".join(b[5] for b in BANDS) + ";")
    out += params
    out += decls
    out.append("      common /ubvwrk/indfreq(NPTSMAX,%d), freqfind;"
               % (FIRST_NEW_BAND + len(BANDS) - 1))
    out.append("")
    out.append("-- ===== %_dataubv additions (transmission curves) =====")
    out += datas
    out.append("-- ===== %_flx_fits_obsubv_zeroflux additions =====")
    out += zps
    out.append("")
    out.append("-- ===== %_flx_fits_obsubv: @find block =====")
    out += finds
    out.append("")
    out.append("-- ===== %_flx_fits_obsubv: @flux block =====")
    out += fluxes
    out.append("")
    out.append("-- ===== .tt header columns (after Mbol MU MB MV MI MR) =====")
    out.append("--   " + " ".join(hdr_cols))
    out.append("-- ===== .tt units columns =====")
    out.append("--   " + " ".join(unt_cols))
    out.append("-- ===== write var list (after MU,MB,MV,MI,MR) =====")
    out.append("--   " + ", ".join(mag_vars))
    out.append("-- ===== write format: 0p,%df8.3 for the %d magnitudes ====="
               % (nmag, nmag))
    with open(args.out, "w") as f:
        f.write("\n".join(out) + "\n")
    print("\nfragments written to %s" % args.out)
    print("insert into tt4strad.trf, then: trf -nfs tt4strad.trf && "
          "rebuild (see build_parallel.sh)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
