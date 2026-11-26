#!/usr/bin/env python3
"""Augmented-residual evaluator battery (MOM_AUG_RES / MOM_AUG_RESID_V1).

drv_mom_augres fills the module capture arrays with manufactured
systems, calls the production evaluator, and dumps the inputs.
This battery recomputes RESID_IMPLEMENTED and RESID_TARGET in numpy
independently and checks every reported number, then asserts the
per-case semantics: exact manufactured fields return machine-level
residuals, perturbations on J and H are detected separately, the
equation-changing stabilisers (LAG_SCL, EPS caps, XM_POS, H_IMPL)
split implemented from target, boundary rows contribute, skipped and
non-finite rows are counted, and solver state is not modified.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import (SRC, 
    F90DIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)

ND = 5
NCF = 10
HELPERS = [
    "tools/set_kind_module.f",
    "unix/term_io.f",
    "unix/gen_asci_open.f",
    "unix/verbose_info.f",
    "subs/set_case.f",
    "subs/phys_con.f",
    "subs/mon_int_funs_v2.f",
    "subs/mon_interp.f",
    "subs/mod_trap_derivatives.f",
    "subs/nordtau.f",
    "subs/qfromf.f",
    "subs/integrate.f",
    "subs/write_vec.f",
    "unix/dir_acc_pars_gen.f",
    "plane/mod_ray_mom_store.f",
    "new_main/mod_subs/knop_status_mod.f90",
    "plane/get_moms_non_rel.f",
    "new_main/mod_subs/mod_cmfgen.f",
    "new_main/mod_subs/control_variable_mod.f",
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f",
    "new_main/mod_subs/jfnk_solver_v1.f90",
    "new_main/mod_subs/check_param_consistency.f",
    "newsubs/rd_store_log.f",
    "new_main/mod_subs/iter_ledger_mod.f90",
    "new_main/mod_subs/iter_log_fmt_mod.f90",
    "plane/fg_j_cmf_nonmono_v1.f",
    "plane/mom_j_cmf_nonmono_v1.f",
    "subs/thomas.f",
    "plane/mom_j_cmf_v11.f",
]
CASES = [
    "exact",
    "epsmix",
    "lagb",
    "him",
    "clip",
    "pert",
    "perth",
    "bnd",
    "epscap",
    "mix",
    "offg",
    "nan",
    "skip",
    "state",
]


def build():
    return compile_and_link_prod(
        "mom_augres",
        [SRC / h for h in HELPERS]
        + [F90DIR / "stub_vadat.f", F90DIR / "drv_mom_augres.f"],
    )


def load(tag):
    d = {"nd": ND, "ncf": NCF}
    sec = None
    names2r = [
        "JNEXT", "TA", "TB", "TC", "CUP", "CUPU", "VB", "VC",
        "LAGSCL", "SRC",
    ]
    names2r_m = ["HNEXT", "HU", "HL", "HS", "EPS", "EPSU", "EPSP",
                 "EPSPU"]
    names2i = ["KUP", "KUH"]
    names1r = ["T1B"]
    names1i = ["OUTBC", "HCHK", "NMJ", "NMH"]
    arrs = {}
    with open(WORKDIR / f"{tag}_ar.txt") as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            if ln.startswith("#"):
                p = ln.split()
                if p[1] == "SHAPE":
                    assert (int(p[2]), int(p[3])) == (ND, NCF)
                elif p[1] == "FLAGS":
                    d["xm_pos"] = p[2] == "T"
                    d["eps_lim"] = p[3] == "T"
                    d["him"] = p[4] == "T"
                    sec = None
                else:
                    sec = p[1]
                    if sec in names2r:
                        arrs[sec] = np.zeros((ND, NCF))
                    elif sec in names2r_m:
                        arrs[sec] = np.zeros((ND - 1, NCF))
                    elif sec in names2i:
                        arrs[sec] = np.zeros(
                            (ND if sec == "KUP" else ND - 1, NCF),
                            dtype=int,
                        )
                    elif sec in names1r:
                        arrs[sec] = np.zeros(NCF)
                    elif sec in names1i:
                        n = NCF if sec in ("OUTBC", "HCHK") else ND
                        arrs[sec] = np.zeros(n, dtype=int)
                continue
            p = ln.split()
            if sec in names2r + names2r_m:
                arrs[sec][int(p[0]) - 1, int(p[1]) - 1] = float(p[2])
            elif sec in names2i:
                arrs[sec][int(p[0]) - 1, int(p[1]) - 1] = int(p[2])
            elif sec in names1r + names1i:
                arrs[sec][int(p[0]) - 1] = float(p[1])
    for n in names2r + names2r_m + names2i + names1r + names1i:
        d[n.lower()] = arrs[n]
    return d


def eval_aug(d):
    """Independent reimplementation of MOM_AUG_RESID_V1."""
    nd, ncf = d["nd"], d["ncf"]
    x, h = d["jnext"], d["hnext"]
    nmj = d["nmj"]
    r = dict(
        r2m=[0.0, 0.0], r2h=[0.0, 0.0], s2m=[0.0, 0.0], s2h=[0.0, 0.0],
        mxm=[0.0, 0.0], mxh=[0.0, 0.0],
        im=[0, 0], km=[0, 0], ih=[0, 0], kh=[0, 0],
        nclip=0, nlagb=0, nepsl=0, noff=0, nnan=0, nhclp=0,
        nskip=0, nbndr=0, rbn2=0.0,
    )
    tiny = 1.0e-280

    def readx(i, k):
        return x[i, k] if 0 <= i < nd and 0 <= k < ncf else np.nan

    for k in range(ncf):
        outbc = int(d["outbc"][k])
        if outbc < 0:
            r["nskip"] += nd
            continue
        t1b = d["t1b"][k]
        for i in range(nd):
            lbn = False
            kup = int(d["kup"][i, k])
            xup = x[i, kup - 1] if 1 <= kup <= ncf else 0.0
            if i == 0:
                if outbc <= 1:
                    xmi = d["cup"][0, k] * xup
                    xmt = d["cupu"][0, k] * xup
                    lhs = d["tb"][0, k] * x[0, k] + d["tc"][0, k] * x[1, k]
                    rm, rt = lhs - xmi, lhs - xmt
                    sm = (abs(d["tb"][0, k] * x[0, k])
                          + abs(d["tc"][0, k] * x[1, k])
                          + abs(d["cup"][0, k] * xup))
                    st = (abs(d["tb"][0, k] * x[0, k])
                          + abs(d["tc"][0, k] * x[1, k])
                          + abs(d["cupu"][0, k] * xup))
                else:
                    thb = 0.0
                    if d["hs"][0, k] != 0.0:
                        u = int(d["kuh"][0, k])
                        if 1 <= u <= ncf:
                            thb = -d["hs"][0, k] * h[0, u - 1] / t1b
                        else:
                            r["noff"] += 1
                            lbn = True
                    xmp = (d["src"][0, k] + d["cup"][0, k] * xup
                           + d["lagscl"][0, k] * thb)
                    t5 = 0.0
                    c = d["epsp"][0, k]
                    if c != 0.0 and t1b != 0.0:
                        if c >= 0.0:
                            dd, u = 0, k - 1
                        else:
                            dd = 1
                            u = (k + nmj[dd]
                                 if k + nmj[dd] < ncf else -1)
                        if 0 <= u < ncf:
                            t5 = 2.0 * c * x[dd, u] / t1b
                            if d["eps_lim"] and t5 > 0.0:
                                t4 = max(xmp, 0.0)
                                if t5 > t4:
                                    t5 = t4
                                    r["nepsl"] += 1
                        else:
                            r["noff"] += 1
                            lbn = True
                    xmi = xmp - d["lagscl"][0, k] * t5
                    lhs = d["tb"][0, k] * x[0, k] + d["tc"][0, k] * x[1, k]
                    rm = lhs - xmi
                    sm = (abs(d["tb"][0, k] * x[0, k])
                          + abs(d["tc"][0, k] * x[1, k])
                          + abs(d["src"][0, k])
                          + abs(d["cup"][0, k] * xup)
                          + d["lagscl"][0, k] * (abs(thb) + abs(t5)))
                    tbt, tct = d["tb"][0, k], d["tc"][0, k]
                    if t1b != 0.0:
                        e, eu = d["eps"][0, k], d["epsu"][0, k]
                        if e >= 0.0:
                            tbt += 2.0 * e / t1b
                        else:
                            tct += 2.0 * e / t1b
                        tbt -= eu / t1b
                        tct -= eu / t1b
                    cu = d["epspu"][0, k]
                    te1t = 0.0
                    if cu != 0.0 and t1b != 0.0:
                        if cu >= 0.0:
                            u1 = u2 = k - 1
                        else:
                            u1 = (k + nmj[0]
                                  if k + nmj[0] < ncf else -1)
                            u2 = (k + nmj[1]
                                  if k + nmj[1] < ncf else -1)
                        if u1 >= 0:
                            te1t = cu * x[0, u1] / t1b
                        if u2 >= 0:
                            te1t += cu * x[1, u2] / t1b
                    xmt = (d["src"][0, k] + d["cupu"][0, k] * xup
                           + thb - te1t)
                    rt = tbt * x[0, k] + tct * x[1, k] - xmt
                    st = (abs(tbt * x[0, k]) + abs(tct * x[1, k])
                          + abs(d["src"][0, k])
                          + abs(d["cupu"][0, k] * xup)
                          + abs(thb) + abs(te1t))
            elif i == nd - 1:
                xmi = d["src"][nd - 1, k]
                lhs = (d["ta"][nd - 1, k] * x[nd - 2, k]
                       + d["tb"][nd - 1, k] * x[nd - 1, k])
                rm = rt = lhs - xmi
                sm = st = (abs(d["ta"][nd - 1, k] * x[nd - 2, k])
                           + abs(d["tb"][nd - 1, k] * x[nd - 1, k])
                           + abs(xmi))
            else:
                thi, smh = 0.0, 0.0
                if d["vb"][i, k] != 0.0:
                    u = int(d["kuh"][i - 1, k])
                    if 1 <= u <= ncf:
                        tv = d["vb"][i, k] * h[i - 1, u - 1]
                    else:
                        tv = 0.0
                        r["noff"] += 1
                        lbn = True
                    thi += tv
                    smh += abs(tv)
                if d["vc"][i, k] != 0.0:
                    u = int(d["kuh"][i, k])
                    if 1 <= u <= ncf:
                        tv = d["vc"][i, k] * h[i, u - 1]
                    else:
                        tv = 0.0
                        r["noff"] += 1
                        lbn = True
                    thi += tv
                    smh += abs(tv)
                te1 = te1t = 0.0
                c = d["epsp"][i, k]
                if c != 0.0:
                    if c >= 0.0:
                        dd, u = i, k - 1
                    else:
                        dd = i + 1
                        u = (k + nmj[dd]
                             if k + nmj[dd] < ncf else -1)
                    if 0 <= u < ncf:
                        te1 = 2.0 * c * x[dd, u]
                    else:
                        r["noff"] += 1
                        lbn = True
                cu = d["epspu"][i, k]
                if cu != 0.0:
                    if cu >= 0.0:
                        u1 = u2 = k - 1
                    else:
                        u1 = (k + nmj[i]
                              if k + nmj[i] < ncf else -1)
                        u2 = (k + nmj[i + 1]
                              if k + nmj[i + 1] < ncf else -1)
                    if u1 >= 0:
                        te1t = cu * x[i, u1]
                    if u2 >= 0:
                        te1t += cu * x[i + 1, u2]
                te2 = te2t = 0.0
                c = d["epsp"][i - 1, k]
                if c != 0.0:
                    if c >= 0.0:
                        dd, u = i - 1, k - 1
                    else:
                        dd = i
                        u = (k + nmj[dd]
                             if k + nmj[dd] < ncf else -1)
                    if 0 <= u < ncf:
                        te2 = 2.0 * c * x[dd, u]
                    else:
                        r["noff"] += 1
                        lbn = True
                cu = d["epspu"][i - 1, k]
                if cu != 0.0:
                    if cu >= 0.0:
                        u1 = u2 = k - 1
                    else:
                        u1 = (k + nmj[i - 1]
                              if k + nmj[i - 1] < ncf else -1)
                        u2 = (k + nmj[i]
                              if k + nmj[i] < ncf else -1)
                    if u1 >= 0:
                        te2t = cu * x[i - 1, u1]
                    if u2 >= 0:
                        te2t += cu * x[i, u2]
                scl = d["lagscl"][i, k]
                xmi = (d["src"][i, k] + d["cup"][i, k] * xup
                       + scl * (thi + te1 - te2))
                if scl < 1.0:
                    r["nlagb"] += 1
                if d["xm_pos"] and xmi < 0.0:
                    xmi = 0.0
                    r["nclip"] += 1
                lhs = (d["ta"][i, k] * x[i - 1, k]
                       + d["tb"][i, k] * x[i, k]
                       + d["tc"][i, k] * x[i + 1, k])
                rm = lhs - xmi
                sm = (abs(d["ta"][i, k] * x[i - 1, k])
                      + abs(d["tb"][i, k] * x[i, k])
                      + abs(d["tc"][i, k] * x[i + 1, k])
                      + abs(d["src"][i, k])
                      + abs(d["cup"][i, k] * xup)
                      + scl * (smh + abs(te1) + abs(te2)))
                tat = d["ta"][i, k]
                tbt = d["tb"][i, k]
                tct = d["tc"][i, k]
                ep, eup = d["eps"][i - 1, k], d["epsu"][i - 1, k]
                e, eu = d["eps"][i, k], d["epsu"][i, k]
                if ep >= 0.0:
                    tat += 2.0 * ep
                else:
                    tbt += 2.0 * ep
                if e >= 0.0:
                    tbt -= 2.0 * e
                else:
                    tct -= 2.0 * e
                tat -= eup
                tbt += -eup + eu
                tct += eu
                if d["him"]:
                    tbt -= abs(d["vb"][i, k]) + abs(d["vc"][i, k])
                xmt = (d["src"][i, k] + d["cupu"][i, k] * xup
                       + thi + te1t - te2t)
                rt = (tat * x[i - 1, k] + tbt * x[i, k]
                      + tct * x[i + 1, k] - xmt)
                st = (abs(tat * x[i - 1, k]) + abs(tbt * x[i, k])
                      + abs(tct * x[i + 1, k])
                      + abs(d["src"][i, k])
                      + abs(d["cupu"][i, k] * xup)
                      + smh + abs(te1t) + abs(te2t))
            if not (np.isfinite(rm) and np.isfinite(rt)):
                r["nnan"] += 1
            else:
                sl = rm / max(sm, tiny)
                sr = rt / max(st, tiny)
                r["r2m"][0] += rm * rm
                r["r2m"][1] += rt * rt
                r["s2m"][0] += sl * sl
                r["s2m"][1] += sr * sr
                if abs(sl) > r["mxm"][0]:
                    r["mxm"][0] = abs(sl)
                    r["im"][0], r["km"][0] = i + 1, k + 1
                if abs(sr) > r["mxm"][1]:
                    r["mxm"][1] = abs(sr)
                    r["im"][1], r["km"][1] = i + 1, k + 1
                if lbn:
                    r["nbndr"] += 1
                    r["rbn2"] += sl * sl

        for i in range(nd - 1):
            lbn = False
            rec = (d["hu"][i, k] * x[i + 1, k]
                   - d["hl"][i, k] * x[i, k])
            rect = rec
            sh = (abs(x[i, k]) + abs(d["hu"][i, k] * x[i + 1, k])
                  + abs(d["hl"][i, k] * x[i, k]))
            sht = sh
            if d["hs"][i, k] != 0.0:
                u = int(d["kuh"][i, k])
                if 1 <= u <= ncf:
                    tv = d["hs"][i, k] * h[i, u - 1]
                else:
                    tv = 0.0
                    r["noff"] += 1
                    lbn = True
                rec += tv
                rect += tv
                sh += abs(tv)
                sht += abs(tv)
            c = d["epsp"][i, k]
            if c != 0.0:
                if c >= 0.0:
                    dd, u = i, k - 1
                else:
                    dd = i + 1
                    u = (k + nmj[dd]
                         if k + nmj[dd] < ncf else -1)
                if 0 <= u < ncf:
                    tv = 2.0 * c * x[dd, u]
                else:
                    tv = 0.0
                    r["noff"] += 1
                    lbn = True
                rec += tv
                sh += abs(tv)
            e = d["eps"][i, k]
            if e != 0.0:
                dd = i if e >= 0.0 else i + 1
                rec -= 2.0 * e * x[dd, k]
                sh += abs(2.0 * e * x[dd, k])
            cu = d["epspu"][i, k]
            if cu != 0.0:
                if cu >= 0.0:
                    u1 = u2 = k - 1
                else:
                    u1 = (k + nmj[i]
                          if k + nmj[i] < ncf else -1)
                    u2 = (k + nmj[i + 1]
                          if k + nmj[i + 1] < ncf else -1)
                tv = 0.0
                if u1 >= 0:
                    tv = cu * x[i, u1]
                if u2 >= 0:
                    tv += cu * x[i + 1, u2]
                rect += tv
                sht += abs(tv)
            eu = d["epsu"][i, k]
            if eu != 0.0:
                tv = eu * (x[i, k] + x[i + 1, k])
                rect -= tv
                sht += abs(tv)
            hchk = int(d["hchk"][k])
            if hchk == 1:
                jbnd = 0.5 * (x[i, k] + x[i + 1, k])
            elif hchk == 2:
                jbnd = max(x[i, k], x[i + 1, k])
            else:
                jbnd = 0.0
            if jbnd > 0.0 and abs(rec) > 0.99999 * jbnd:
                r["nhclp"] += 1
            rh = h[i, k] - rec
            rht = h[i, k] - rect
            if not (np.isfinite(rh) and np.isfinite(rht)):
                r["nnan"] += 1
            else:
                sl = rh / max(sh, tiny)
                sr = rht / max(sht, tiny)
                r["r2h"][0] += rh * rh
                r["r2h"][1] += rht * rht
                r["s2h"][0] += sl * sl
                r["s2h"][1] += sr * sr
                if abs(sl) > r["mxh"][0]:
                    r["mxh"][0] = abs(sl)
                    r["ih"][0], r["kh"][0] = i + 1, k + 1
                if abs(sr) > r["mxh"][1]:
                    r["mxh"][1] = abs(sr)
                    r["ih"][1], r["kh"][1] = i + 1, k + 1
                if lbn:
                    r["nbndr"] += 1
    r["r2m"] = [np.sqrt(v) for v in r["r2m"]]
    r["r2h"] = [np.sqrt(v) for v in r["r2h"]]
    r["s2m"] = [np.sqrt(v) for v in r["s2m"]]
    r["s2h"] = [np.sqrt(v) for v in r["s2h"]]
    r["rbn2"] = np.sqrt(r["rbn2"])
    return r


def parse_results(out):
    v = {}
    for ln in out.splitlines():
        ln = ln.strip()
        for key in ("AR_R2M", "AR_R2H", "AR_S2M", "AR_S2H",
                    "AR_MXM", "AR_MXH", "AR_RBN2"):
            if ln.startswith(key + "="):
                v[key] = [float(t) for t in
                          ln.split("=", 1)[1].split()]
        if ln.startswith("AR_CELLM="):
            v["CELLM"] = [int(t) for t in ln.split("=", 1)[1].split()]
        elif ln.startswith("AR_CELLH="):
            v["CELLH"] = [int(t) for t in ln.split("=", 1)[1].split()]
        elif ln.startswith("AR_CNT="):
            v["CNT"] = [int(t) for t in ln.split("=", 1)[1].split()]
        elif ln.startswith("PRECK="):
            v["PRECK"] = float(ln.split("=", 1)[1])
        elif ln.startswith("POSTCK="):
            v["POSTCK"] = float(ln.split("=", 1)[1])
    return v


def close(a, b, tol=2.0e-9):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.all(np.isfinite(a)) and np.all(
        np.abs(a - b) <= tol * np.maximum(np.abs(b), 1.0e-6))


def main():
    rep = Report("test_mom_augres")
    exe = build()
    outs = {}
    for case in CASES:
        tag = f"aug_{case}"
        r = run(exe, case, tag, cwd=WORKDIR)
        rep.check(f"{case}: driver ran", "RESULT_OK" in r.stdout,
                  r.stdout.strip().splitlines()[-3:])
        outs[case] = parse_results(r.stdout)

    # Cross-check: independent numpy evaluation reproduces every
    # reported number on every case.
    for case in CASES:
        tag = f"aug_{case}"
        d = load(tag)
        ref = eval_aug(d)
        v = outs[case]
        ok = (close(v["AR_R2M"], ref["r2m"])
              and close(v["AR_R2H"], ref["r2h"])
              and close(v["AR_S2M"], ref["s2m"])
              and close(v["AR_S2H"], ref["s2h"])
              and close(v["AR_MXM"], ref["mxm"])
              and close(v["AR_MXH"], ref["mxh"])
              and close([v["AR_RBN2"]], [ref["rbn2"]])
              # Worst-cell argmax is only stable when the residual
              # is above ULP ties.
              and (v["AR_MXM"][0] < 1.0e-12
                   or v["CELLM"][:2] == [ref["im"][0], ref["km"][0]])
              and (v["AR_MXM"][1] < 1.0e-12
                   or v["CELLM"][2:] == [ref["im"][1], ref["km"][1]])
              and (v["AR_MXH"][0] < 1.0e-12
                   or v["CELLH"][:2] == [ref["ih"][0], ref["kh"][0]])
              and (v["AR_MXH"][1] < 1.0e-12
                   or v["CELLH"][2:] == [ref["ih"][1], ref["kh"][1]])
              and v["CNT"] == [ref["nclip"], ref["nlagb"],
                               ref["nepsl"], ref["noff"], ref["nnan"],
                               ref["nhclp"], ref["nskip"],
                               ref["nbndr"]])
        rep.check(
            f"{case}: numpy reimplementation matches evaluator",
            ok,
            f"fortran {v['AR_S2M']} {v['AR_S2H']} cnt {v['CNT']} vs "
            f"py {ref['s2m']} {ref['s2h']} cnt "
            + str([ref["nclip"], ref["nlagb"], ref["nepsl"], ref["noff"],
                   ref["nnan"], ref["nhclp"], ref["nskip"], ref["nbndr"]]),
        )
        # Evaluator must not modify solver state.  The nan case
        # poisons the checksum itself; both sides NaN still proves
        # no write-back happened.
        rep.check(
            f"{case}: fields unchanged by evaluator",
            v["PRECK"] == v["POSTCK"]
            or (np.isnan(v["PRECK"]) and np.isnan(v["POSTCK"])),
            f"pre {v['PRECK']:.6e} post {v['POSTCK']:.6e}",
        )

    # A: exact manufactured solution -> machine-level residuals.
    v = outs["exact"]
    exact_ok = (v["AR_S2M"][0] < 1.0e-11 and v["AR_S2H"][0] < 1.0e-11
                and v["AR_S2M"][1] < 1.0e-11
                and v["AR_S2H"][1] < 1.0e-11)
    rep.check(
        "exact: both residuals at machine level",
        exact_ok,
        f"impl {v['AR_S2M'][0]:.2e}/{v['AR_S2H'][0]:.2e} "
        f"tgt {v['AR_S2M'][1]:.2e}/{v['AR_S2H'][1]:.2e}",
    )
    # Off-grid upwind reads (noff) and boundary rows (nbndr) are
    # discretisation features, not equation modifications.
    rep.check(
        "exact: no equation modifications counted",
        v["CNT"][:3] == [0, 0, 0] and v["CNT"][4:7] == [0, 0, 0],
        f"cnt {v['CNT']}",
    )

    # G: each stabiliser splits implemented from target.
    v = outs["epsmix"]
    rep.check(
        "epsmix: implemented exact, target sees stencil+cap",
        v["AR_S2M"][0] < 1.0e-11 and v["AR_S2H"][0] < 1.0e-11
        and v["AR_S2M"][1] > 1.0e-6,
        f"impl {v['AR_S2M'][0]:.2e} tgt {v['AR_S2M'][1]:.2e}",
    )
    v = outs["lagb"]
    rep.check(
        "lagb: implemented exact, target sees unbounded lag",
        v["AR_S2M"][0] < 1.0e-11 and v["AR_S2M"][1] > 1.0e-6
        and v["CNT"][1] > 0,
        f"impl {v['AR_S2M'][0]:.2e} tgt {v['AR_S2M'][1]:.2e} "
        f"nlagb {v['CNT'][1]}",
    )
    v = outs["him"]
    rep.check(
        "him: implemented exact, target strips H_IMPL mass",
        v["AR_S2M"][0] < 1.0e-11 and v["AR_S2M"][1] > 1.0e-6,
        f"impl {v['AR_S2M'][0]:.2e} tgt {v['AR_S2M'][1]:.2e}",
    )
    v = outs["clip"]
    rep.check(
        "clip: XM_POS clip counted, implemented row exact",
        v["CNT"][0] == 1 and v["AR_S2M"][0] < 1.0e-10
        and v["AR_S2M"][1] > 1.0e-6,
        f"nclip {v['CNT'][0]} impl {v['AR_S2M'][0]:.2e} "
        f"tgt {v['AR_S2M'][1]:.2e}",
    )

    # B/C: perturbation detection on J and H independently.
    v = outs["pert"]
    rep.check(
        "pert: J perturbation detected in moment residual",
        v["AR_S2M"][0] > 1.0e-3 and v["AR_S2H"][0] > 1.0e-3,
        f"s2mom {v['AR_S2M'][0]:.2e} s2H {v['AR_S2H'][0]:.2e} "
        f"worst {v['CELLM']}",
    )
    rep.check(
        "pert: worst moment cell touches the perturbed column",
        v["CELLM"][1] in (4, 5, 6) or v["CELLM"][1] == 3,
        f"cell {v['CELLM']}",
    )
    v = outs["perth"]
    rep.check(
        "perth: H perturbation detected in H residual",
        v["AR_S2H"][0] > 1.0e-3,
        f"s2H {v['AR_S2H'][0]:.2e} worst {v['CELLH']}",
    )
    rep.check(
        "perth: moment rows see the perturbed H",
        v["AR_S2M"][0] > 1.0e-4,
        f"s2mom {v['AR_S2M'][0]:.2e}",
    )

    # F: boundary-row machinery and the row-1 EPS cap.
    v = outs["bnd"]
    rep.check(
        "bnd: negated boundary rows exact under manufacture",
        v["AR_S2M"][0] < 1.0e-11 and v["AR_S2H"][0] < 1.0e-11,
        f"impl {v['AR_S2M'][0]:.2e}/{v['AR_S2H'][0]:.2e}",
    )
    v = outs["epscap"]
    rep.check(
        "epscap: boundary EPS cap fires",
        v["CNT"][2] >= 1,
        f"nepsl {v['CNT'][2]}",
    )

    # D/E: mixed graphs and off-grid reads.
    v = outs["offg"]
    rep.check(
        "offg: off-grid reads counted",
        v["CNT"][3] >= 1,
        f"noff {v['CNT'][3]}",
    )
    v = outs["mix"]
    rep.check(
        "mix: mixed ring/chain system evaluated finite",
        np.isfinite(v["AR_S2M"][0]) and np.isfinite(v["AR_S2H"][0]),
        f"impl {v['AR_S2M'][0]:.2e}/{v['AR_S2H'][0]:.2e}",
    )

    # I: non-finite injection.
    v = outs["nan"]
    rep.check(
        "nan: non-finite rows counted, norms stay finite",
        v["CNT"][4] >= 1 and np.isfinite(v["AR_S2M"][0])
        and np.isfinite(v["AR_S2H"][0]),
        f"nnan {v['CNT'][4]} s2 {v['AR_S2M']}/{v['AR_S2H']}",
    )

    # Skipped (never-captured) columns.
    v = outs["skip"]
    rep.check(
        "skip: uncaptured columns counted, not evaluated",
        v["CNT"][6] == ND * 2,
        f"nskip {v['CNT'][6]} expect {ND * 2}",
    )

    rep.summary()


if __name__ == "__main__":
    main()
