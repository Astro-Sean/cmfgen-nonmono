#!/usr/bin/env python3
"""DO_VAR_CONT nonmonotonic-branch coverage
(new_main/mod_subs/do_var_cont.f).

Alteration under test: under DO_NONMONO_SOB/DO_NONMONO_KNOP the
continuum dJ/dCHI, dJ/dETA variation diagonals come from the
approximate-Lambda operators stored by the forward solver
(FG_J_CMF_NONMONO_V1, MOM_J_CMF_NONMONO_V1 or FG_J_CMF_V13) --
upstream VAR_MOM_J_CMF_V12 linearizes the monotonic V11 system and
returned NaNs on reversal profiles.  The driver injects the stored
diagonals through the same SAVE variables the solvers write, runs
the production do_var_cont.o, and reads TX/VJ back out.

Cases:
  nostore : KNOP + LAM_ES, no store -> ALO_OK=F, TX/VJ stay zero
  fg      : KNOP + LAM_ES + FG store   -> TX3=DJDCHI, TX4=DJDETA,
            TX5=DJDCHI+DJDETA*RJ_ES, VJ diag propagated
  mom     : KNOP + LAM_ES=F + mom store -> mom diagonal used
  momquar : mom store + MOM_CELL_QUAR  -> masked rows zeroed
  fg13    : SOB + LAM_ES + FG13 store  -> V13 diagonal used
  sobfgk  : SOB + LAM_ES + SOB_FG_KNOP + FG store -> FG diagonal
  sobmom  : SOB + LAM_ES=F + mom store -> moment diagonal used
  esali   : FG store + DO_ES_ALI_BA    -> 1/(1-DJDETA*chi_scat)
            amplification floored by ES_ALI_BA_MIN
  tri     : SOB + DO_FG13_ALO + TRI store -> K+-1 off-diagonals
            filled from GET_FG13_ALO_TRI
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (BINDIR, F90DIR, PROD_MOD_INC, TREE, Report,
                    compile_objs, link_prod, prod_obj, run)

ND = 4
DJDC_FG = (-1.0e-2, -2.0e-2, -3.0e-2, -4.0e-2)
DJDE_FG = (0.4, 0.5, 0.6, 0.7)
DJDC_MOM = (-5.0e-2, -6.0e-2, -7.0e-2, -8.0e-2)
DJDE_MOM = (0.9, 0.8, 0.7, 0.6)
DJDC_13 = (-1.1e-2, -2.2e-2, -3.3e-2, -4.4e-2)
DJDE_13 = (0.41, 0.51, 0.61, 0.71)
DJDE_M = (0.0, 0.12, 0.13, 0.14)
DJDC_M = (0.0, -1.1e-3, -1.2e-3, -1.3e-3)
DJDE_P = (0.21, 0.22, 0.23, 0.0)
DJDC_P = (-2.0e-3, -2.1e-3, -2.2e-3, 0.0)
RJ_ES = 0.8
CHI_SCAT = 0.3
VCHI = 0.05
VETA = 0.02


def build():
    md = BINDIR / "obj_dvc"
    objs = [prod_obj(TREE / "new_main" / "mod_subs" / "do_var_cont.f")]
    objs += compile_objs(
        "dvc_drv", [F90DIR / "drv_dvc.f90"],
        moddir=md, extra_flags=PROD_MOD_INC)
    libs = [TREE / "lib" / a for a in
            ("libsubs.a", "libmod.a", "libsubs.a")]
    return link_prod("dvc", objs + libs, moddir=md)


def row(out, key, n):
    m = re.search(rf"^{key}\s+((?:-?\d\.\d+E[+-]\d+\s*)+)$",
                  out, re.M)
    assert m, f"{key} row missing"
    vals = [float(x) for x in m.group(1).split()]
    assert len(vals) == n, (key, vals)
    return vals


def close(a, b, tol=1.0e-5):
    return abs(a - b) <= tol * max(1.0, abs(b))


def diag_expect(rep, out, djdc, djde, mask=None, denb=None,
                tag=""):
    mask = mask or [False] * ND
    denb = denb or [1.0] * ND
    tx3 = row(out, "TX3", ND)
    tx4 = row(out, "TX4", ND)
    tx5 = row(out, "TX5", ND)
    vjd = row(out, "VJD", ND)
    for k in range(ND):
        if mask[k]:
            rep.check(
                f"{tag} row {k} quarantine-zeroed",
                all(abs(arr[k]) < 1.0e-12
                    for arr in (tx3, tx4, tx5, vjd)))
        else:
            rep.check(
                f"{tag} TX3[{k}]==DJDCHI*DENB",
                close(tx3[k], djdc[k] * denb[k]),
                f"{tx3[k]} vs {djdc[k] * denb[k]}")
            rep.check(
                f"{tag} TX4[{k}]==DJDETA*DENB",
                close(tx4[k], djde[k] * denb[k]),
                f"{tx4[k]} vs {djde[k] * denb[k]}")
            rep.check(
                f"{tag} TX5[{k}]==DJDCHI+DJDETA*RJ_ES",
                close(tx5[k],
                      denb[k] * (djdc[k] + djde[k] * RJ_ES)),
                f"{tx5[k]}")
            rep.check(
                f"{tag} VJD[{k}]==VCHI*TX3+VETA*TX4+ESEC*TX5/ED",
                close(vjd[k],
                      denb[k] * (VCHI * djdc[k] + VETA * djde[k]
                                 + 0.5 * (djdc[k] + RJ_ES * djde[k]))),
                f"{vjd[k]}")


def main():
    rep = Report("dvc")
    exe = build()

    # no stored diagonal -> ALO_OK=F, TX/VJ untouched
    out = run(exe, "nostore").stdout
    rep.check("nostore completes", "DONE_" in out)
    rep.check("nostore skips ALO banner (ALO_OK=F)",
              " ALO " not in out and "ALO (" not in out)
    for key in ("TX3", "TX4", "TX5", "VJD"):
        rep.check(
            f"nostore {key} all zero (ALO_OK=F)",
            all(abs(v) < 1.0e-12 for v in row(out, key, ND)))

    # KNOP + LAM_ES -> FG_J_CMF_NONMONO_V1 diagonal
    out = run(exe, "fg").stdout
    rep.check("fg completes", "DONE_" in out)
    rep.check("fg banner says KNOP", "KNOP" in out)
    diag_expect(rep, out, DJDC_FG, DJDE_FG, tag="fg")

    # KNOP + LAM_ES=F -> MOM_J_CMF_NONMONO_V1 diagonal
    out = run(exe, "mom").stdout
    rep.check("mom completes", "DONE_" in out)
    diag_expect(rep, out, DJDC_MOM, DJDE_MOM, tag="mom")

    # quarantine mask zeroes frozen rows
    out = run(exe, "momquar").stdout
    rep.check("momquar completes", "DONE_" in out)
    diag_expect(rep, out, DJDC_MOM, DJDE_MOM,
                mask=[False, True, False, True], tag="momquar")

    # SOB + LAM_ES + DO_FG13_ALO -> FG_J_CMF_V13 diagonal
    out = run(exe, "fg13").stdout
    rep.check("fg13 completes", "DONE_" in out)
    rep.check("fg13 banner says SOB", "SOB" in out)
    diag_expect(rep, out, DJDC_13, DJDE_13, tag="fg13")

    # SOB + DO_SOB_FG_KNOP -> FG nonmono diagonal
    out = run(exe, "sobfgk").stdout
    rep.check("sobfgk completes", "DONE_" in out)
    diag_expect(rep, out, DJDC_FG, DJDE_FG, tag="sobfgk")

    # SOB + LAM_ES=F -> moment diagonal
    out = run(exe, "sobmom").stdout
    rep.check("sobmom completes", "DONE_" in out)
    diag_expect(rep, out, DJDC_MOM, DJDE_MOM, tag="sobmom")

    # DO_ES_ALI_BA: DENB = 1/max(1-DJDETA*CHI_SCAT, ES_ALI_BA_MIN)
    out = run(exe, "esali").stdout
    rep.check("esali completes", "DONE_" in out)
    denb = [1.0 / max(1.0 - d * CHI_SCAT, 0.05) for d in DJDE_FG]
    diag_expect(rep, out, DJDC_FG, DJDE_FG, denb=denb, tag="esali")

    # DO_FG13_ALO_TRI: TX(K,K+-1,1:2) = TRA/TRC * DENB(K), DENB=1
    out = run(exe, "tri").stdout
    rep.check("tri completes", "DONE_" in out)
    diag_expect(rep, out, DJDC_13, DJDE_13, tag="tri")
    tx3m = row(out, "TX3M", ND - 1)
    tx3p = row(out, "TX3P", ND - 1)
    tx4m = row(out, "TX4M", ND)
    tx4p = row(out, "TX4P", ND)
    for k in range(1, ND):
        rep.check(
            f"tri TX3[K,K-1]==DJDCHI_M at K={k}",
            close(tx3m[k - 1], DJDC_M[k]),
            f"{tx3m[k-1]} vs {DJDC_M[k]}")
        rep.check(
            f"tri TX4[K,K-1]==DJDETA_M at K={k}",
            close(tx4m[k - 1], DJDE_M[k]),
            f"{tx4m[k-1]} vs {DJDE_M[k]}")
    for k in range(ND - 1):
        rep.check(
            f"tri TX3[K,K+1]==DJDCHI_P at K={k}",
            close(tx3p[k], DJDC_P[k]),
            f"{tx3p[k]} vs {DJDC_P[k]}")
        rep.check(
            f"tri TX4[K,K+1]==DJDETA_P at K={k}",
            close(tx4p[k + 1], DJDE_P[k]),
            f"{tx4p[k+1]} vs {DJDE_P[k]}")

    rep.summary()


if __name__ == "__main__":
    main()
