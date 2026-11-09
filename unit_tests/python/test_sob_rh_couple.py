#!/usr/bin/env python3
"""Unit tests for the Rybicki-Hummer non-local Sobolev coupling,
subs/sob_rh_couple.f (new file, DO_NONLOC_SOB path).

SOB_RH_GEOM builds a cached table of coupling depth pairs (I,J) with
per-pair angular factors T1=|1+mu^2*SIGMA(I)|, T1P, and pair weight
PW=W_J*R_J*V_J*Upsilon.  The driver USEs SOB_RH_MOD so the dumped
tables are the production ones; Python then verifies:

  * the accepted pair set against an independent O(ND^2) enumeration
    of the documented acceptance rules;
  * pair symmetry (I,J) accepted iff (J,I) accepted;
  * the mu^2 identity 1-mu^2 = (rj/ri)^2 (1-mu'^2);
  * no pairs for homologous V=aR (D2=0 exactly);
  * V<=0 rows never appear as target or partner;
  * DRC equals CHIL*dJ/BB_COR recomputed from the dumped tables;
  * NEG_OPACITY and S_L<=0 suppress contributions; BB_COR=0 is safe;
  * continuum chord extinction (DO_SOB_RH_EXTINCT): analytic uniform-
    CHI attenuation, opaque-medium suppression, EXTINCT=F recovery,
    and shell-localized crossing pairs.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (
    F90DIR,
    PLOTDIR,
    TREE,
    WORKDIR,
    Report,
    compile_and_link_prod,
    run,
)

HELPERS = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "exponx.f",
    TREE / "new_main" / "mod_subs" / "control_variable_mod.f",
    TREE / "subs" / "sob_rh_couple.f",
]

ND, NSIM = 30, 2
FL = 1.0


def build():
    return compile_and_link_prod(
        "sobrh", HELPERS + [F90DIR / "drv_sob_rh.f90"]
    )


def write_case(
    name, r, v, sig, chil, etal, bbc, znet, neg, fl=FL,
    chi=None, extinct=True, extra_geom=None,
):
    path = WORKDIR / name
    chi = np.zeros(len(r)) if chi is None else chi
    with open(path, "w") as f:
        f.write(f"{ND} {NSIM} {fl:.6e}\n")
        f.write(" ".join("T" if x else "F" for x in neg) + "\n")
        f.write(("T" if extinct else "F") + "\n")
        f.write(" ".join(f"{x:.6e}" for x in chi) + "\n")
        for i in range(ND):
            f.write(f"{r[i]:.17e} {v[i]:.17e} {sig[i]:.17e}\n")
        for l in range(NSIM):
            for i in range(ND):
                f.write(
                    f"{chil[i, l]:.17e} {etal[i, l]:.17e} "
                    f"{bbc[i, l]:.17e} {znet[i, l]:.17e}\n"
                )
        # Each extra (r, v, sig) block re-runs GEOM+COUPLE in-process
        # (driver phase >= 2), exercising the cache invalidation.
        for r2, v2, sig2 in extra_geom or []:
            for i in range(ND):
                f.write(f"{r2[i]:.17e} {v2[i]:.17e} {sig2[i]:.17e}\n")
    return path


def parse(out):
    """Split the driver dump into per-phase blocks.

    The SOB_RH: rebuild banner is written inside SOB_RH_GEOM, before
    the driver's PHASE line, so banners are carried into the block they
    precede. A single-phase run returns one (pairs, drc, meta) tuple;
    multi-phase runs return a list of them.
    """
    blocks, pairs, drc, meta, pend = [], [], {}, {}, 0
    for ln in out.splitlines():
        t = ln.split()
        if not t:
            continue
        if t[0] == "PHASE":
            if pairs or drc or meta:
                blocks.append((pairs, drc, meta))
            pairs, drc, meta = [], {}, {"banners": pend}
            pend = 0
        elif t[0] == "SOB_RH:":
            pend += 1
        elif t[0] == "TOTPAIRS":
            meta["tot"] = int(t[1])
        elif t[0] == "PAIR":
            pairs.append((int(t[1]), int(t[3]), float(t[4]), float(t[5]), float(t[6])))
        elif t[0] == "NCALL_ND":
            meta["ncall"] = int(t[1])
        elif t[0] == "MAXDJZ":
            meta["maxdj"] = float(t[1])
            meta["maxdz"] = float(t[2])
        elif t[0] == "DRC":
            drc[int(t[1]), int(t[2])] = float(t[3])
        elif t[0] == "ATT":
            meta.setdefault("att", {})[int(t[1]), int(t[2])] = float(t[3])
    blocks.append((pairs, drc, meta))
    return blocks[0] if len(blocks) == 1 else blocks


def expected_pairs(r, v, sig):
    """Independent enumeration of the documented acceptance rules."""
    nd = len(r)
    wj = np.abs(np.gradient(r))
    wj[0] = abs(r[1] - r[0])
    wj[-1] = abs(r[-1] - r[-2])
    out = {}
    for i in range(nd):
        for j in range(nd):
            if i == j or v[i] <= 0.0 or v[j] <= 0.0:
                continue
            d1 = r[j] ** 2 - r[i] ** 2
            d2 = r[j] ** 2 * v[i] ** 2 - r[i] ** 2 * v[j] ** 2
            if d1 * d2 <= 0.0:
                continue
            mu2 = d1 * v[j] ** 2 / d2
            mup2 = d1 * v[i] ** 2 / d2
            if mu2 <= 0.0 or mu2 > 1.0 or mup2 <= 0.0 or mup2 > 1.0:
                continue
            t1 = max(abs(1.0 + mu2 * sig[i]), 1e-30)
            t1p = max(abs(1.0 + mup2 * sig[j]), 1e-30)
            pw = wj[j] * r[j] * v[j] / np.sqrt(d1 * d2)
            out.setdefault(i, []).append((j, t1, t1p, pw))
    return out


def exponx(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    s = np.abs(x) < 1.0e-3
    out[s] = 1.0 - x[s] * (0.5 - x[s] / 6.0 * (1.0 - x[s] / 4.0))
    m = (~s) & (x < 40.0)
    out[m] = (1.0 - np.exp(-x[m])) / x[m]
    out[~s & ~m] = 1.0 / x[~s & ~m]
    return out


def chord_tau(r, v, i, j, chi):
    """Replicate the production trapezoid: continuum tau along the
    chord between resonance points i and j, summed shell by shell in
    the path coordinate s=sqrt(r^2-p^2)."""
    d1 = r[j] ** 2 - r[i] ** 2
    d2 = r[j] ** 2 * v[i] ** 2 - r[i] ** 2 * v[j] ** 2
    mu2 = d1 * v[j] ** 2 / d2
    p = r[i] * np.sqrt(max(1.0 - mu2, 0.0))
    tau = 0.0
    for k in range(min(i, j), max(i, j)):
        sa = np.sqrt(max((r[k] - p) * (r[k] + p), 0.0))
        sb = np.sqrt(max((r[k + 1] - p) * (r[k + 1] + p), 0.0))
        tau += 0.5 * (chi[k] + chi[k + 1]) * abs(sa - sb)
    return tau


def drc_expected(pairs, r, v, chil, etal, bbc, neg, fl, chi=None, extinct=True):
    nd, nsim = chil.shape
    cfac = 3.0e-10 / fl
    drc = np.zeros_like(chil)
    for l in range(nsim):
        gamh = chil[:, l] * cfac * r / np.maximum(np.abs(v), 1.0)
        sl = np.where((chil[:, l] > 0) & (etal[:, l] > 0), etal[:, l] / chil[:, l], 0.0)
        for i, plist in pairs.items():
            if neg[i] or gamh[i] <= 0.0:
                continue
            dj = 0.0
            for j, t1, t1p, pw in plist:
                if neg[j] or sl[j] <= 0.0:
                    continue
                att = 1.0
                if extinct and chi is not None:
                    att = np.exp(-min(chord_tau(r, v, i, j, chi), 700.0))
                ti = gamh[i] / t1
                tj = gamh[j] / t1p
                dj += pw * gamh[j] * sl[j] * att * exponx(ti) * exponx(tj)
            dj *= 0.5
            if dj <= 0.0 or bbc[i, l] == 0.0:
                continue
            drc[i, l] += dj * chil[i, l] / bbc[i, l]
    return drc


def base_inputs(nd=ND, nsim=NSIM):
    rng = np.random.default_rng(7)
    chil = rng.uniform(1.0e4, 1.0e6, (nd, nsim))
    etal = chil * rng.uniform(0.5, 2.0, (nd, nsim))
    bbc = np.ones((nd, nsim))
    znet = np.zeros((nd, nsim))
    neg = np.zeros(nd, bool)
    return chil, etal, bbc, znet, neg


def main():
    rep = Report("SOB_RH non-local coupling")
    binary = build()
    r = np.geomspace(1.0e5, 1.0e3, ND)

    # -- T1: homologous V=aR -> D2=0 exactly, zero pairs ----------
    v = r / r[0] * 8000.0
    sig = np.zeros(ND)
    chil, etal, bbc, znet, neg = base_inputs()
    path = write_case("rh_homo.txt", r, v, sig, chil, etal, bbc, znet, neg)
    pairs, drc, meta = parse(run(binary, path).stdout)
    rep.check(
        "homologous V: TOTPAIRS=0, all DRC=0",
        meta.get("tot") == 0 and all(x == 0 for x in drc.values()),
    )

    # -- T2: non-monotonic dip -> pairs exist and match the -------
    # independent enumeration, with the mu identity holding.
    v = np.linspace(8000.0, 1000.0, ND)
    v[12:18] = [3000.0, 1500.0, 1200.0, 2000.0, 3500.0, 4500.0]
    sig = np.full(ND, -0.5)
    path = write_case("rh_dip.txt", r, v, sig, chil, etal, bbc, znet, neg)
    pairs, drc, meta = parse(run(binary, path).stdout)
    exp_pairs = expected_pairs(r, v, sig)
    got = {}
    for i, j, t1, t1p, pw in pairs:
        got.setdefault(i - 1, []).append((j - 1, t1, t1p, pw))
    set_got = {(i, p[0]) for i, pl in got.items() for p in pl}
    set_exp = {(i, p[0]) for i, pl in exp_pairs.items() for p in pl}
    rep.check(
        "dip: pair set matches independent enumeration",
        set_got == set_exp and meta.get("tot") == len(set_exp),
        f"got={len(set_got)} exp={len(set_exp)}",
    )

    sym = all((j, i) in set_got for (i, j) in set_got)
    rep.check("dip: pair symmetry (i,j)<=>(j,i)", sym)

    ident_ok = True
    for i, pl in got.items():
        for j, t1, t1p, pw in pl:
            # recover mu^2 from T1=|1+mu2*sig|, sig=-0.5
            # use identity instead: 1-mu2 = (rj/ri)^2 (1-mup2)
            # verify via PW ratio instead: T1/T1P encodes nothing
            # directly; recompute mu2 independently.
            d1 = r[j] ** 2 - r[i] ** 2
            d2 = r[j] ** 2 * v[i] ** 2 - r[i] ** 2 * v[j] ** 2
            mu2 = d1 * v[j] ** 2 / d2
            mup2 = d1 * v[i] ** 2 / d2
            if abs(t1 - abs(1.0 + mu2 * sig[i])) > 1e-12:
                ident_ok = False
            if abs(t1p - abs(1.0 + mup2 * sig[j])) > 1e-12:
                ident_ok = False
            # 1-mu2 = (rj/ri)^2 (1-mup2) follows from D1/D2 algebra
            ident = (1.0 - mu2) - (r[j] / r[i]) ** 2 * (1.0 - mup2)
            if abs(ident) > 1e-12:
                ident_ok = False
    rep.check("dip: T1/T1P tables + mu^2 identity correct", ident_ok)

    # -- T2b: W profile (two reversals) -> denser multi-resonance --
    # structure; pair set + DRC still match the enumeration.
    v = np.linspace(8000.0, 1000.0, ND)
    v[8:12] = [2500.0, 1500.0, 3000.0, 4500.0]
    v[18:22] = [2500.0, 1200.0, 2800.0, 4200.0]
    sig = np.full(ND, -0.5)
    path = write_case("rh_w.txt", r, v, sig, chil, etal, bbc, znet, neg)
    pairs, drc, meta = parse(run(binary, path).stdout)
    exp_pairs = expected_pairs(r, v, sig)
    got = {}
    for i, j, t1, t1p, pw in pairs:
        got.setdefault(i - 1, []).append((j - 1, t1, t1p, pw))
    set_got = {(i, p[0]) for i, pl in got.items() for p in pl}
    set_exp = {(i, p[0]) for i, pl in exp_pairs.items() for p in pl}
    drc_exp = drc_expected(exp_pairs, r, v, chil, etal, bbc, neg, FL)
    drc_arr = np.zeros((ND, NSIM))
    for (l, i), val in drc.items():
        drc_arr[i - 1, l - 1] = val
    rep.check(
        "W profile: pair set + DRC match enumeration",
        set_got == set_exp
        and meta.get("tot") == len(set_exp)
        and np.allclose(drc_arr, drc_exp, rtol=1e-9, atol=1e-30),
        f"tot={meta.get('tot')} exp={len(set_exp)} "
        f"max|dDRC|={np.abs(drc_arr - drc_exp).max():.2e}",
    )

    # -- T3: DRC equals the formula recomputed from the tables ----
    dexp = drc_expected(exp_pairs, r, v, chil, etal, bbc, neg, FL)
    dgot = np.zeros_like(chil)
    for (l, i), x in drc.items():
        dgot[i - 1, l - 1] = x
    err = np.abs(dgot - dexp).max() / max(1.0, np.abs(dexp).max())
    rep.check(
        "dip: DRC == recomputed rate term",
        err < 1e-12,
        f"relerr={err:.2e} max={np.abs(dexp).max():.3e}",
    )

    # -- T4: V<=0 rows dropped ------------------------------------
    v2 = v.copy()
    v2[15] = -100.0  # infalling cell
    v2[16] = 0.0  # stagnation
    path = write_case("rh_stag.txt", r, v2, sig, chil, etal, bbc, znet, neg)
    pairs2, _, meta2 = parse(run(binary, path).stdout)
    set2 = {(i - 1, j - 1) for i, j, t1, t1p, pw in pairs2}
    rows_viol = {15, 16} & (set(i for i, j in set2) | set(j for i, j in set2))
    exp2 = expected_pairs(r, v2, sig)
    set_exp2 = {(i, p[0]) for i, pl in exp2.items() for p in pl}
    rep.check("V<=0 rows never couple; set matches", not rows_viol and set2 == set_exp2)

    # -- T5: NEG_OPACITY suppresses target and partner ------------
    neg2 = np.zeros(ND, bool)
    neg2[13] = True
    path = write_case("rh_neg.txt", r, v, sig, chil, etal, bbc, znet, neg2)
    pairs3, drc3, _ = parse(run(binary, path).stdout)
    dexp3 = drc_expected(exp_pairs, r, v, chil, etal, bbc, neg2, FL)
    dgot3 = np.zeros_like(chil)
    for (l, i), x in drc3.items():
        dgot3[i - 1, l - 1] = x
    err = np.abs(dgot3 - dexp3).max() / max(1.0, np.abs(dexp3).max())
    rep.check(
        "NEG_OPACITY row: DRC matches suppressed model",
        err < 1e-12 and dgot3[13].max() == 0.0,
        f"relerr={err:.2e}",
    )

    # -- T6: BB_COR=0 and ETAL=0 are safe -------------------------
    bbc0 = bbc.copy()
    bbc0[:, :] = 0.0
    path = write_case("rh_bbc0.txt", r, v, sig, chil, etal, bbc0, znet, neg)
    _, drc4, _ = parse(run(binary, path).stdout)
    rep.check(
        "BB_COR=0 everywhere: DRC all zero, no NaN",
        all(x == 0.0 for x in drc4.values()),
    )
    etal0 = etal.copy()
    etal0[:, :] = 0.0
    path = write_case("rh_etal0.txt", r, v, sig, chil, etal0, bbc, znet, neg)
    _, drc5, _ = parse(run(binary, path).stdout)
    rep.check("ETAL=0 (S_L=0): DRC all zero", all(x == 0.0 for x in drc5.values()))

    # -- T7: cache invalidation on the SIGMA input (BUG-002) --------
    # Phase 2 reuses R,V but passes sigma=-0.9: the pair selection is
    # unchanged (it depends on R,V only) but T1=|1+mu^2*sigma| must be
    # rebuilt. The old SUM(V)/SUM(R) fingerprint never looked at SIGMA,
    # so phase 2 kept the phase-1 tables.
    # Phase 3 repeats the phase-2 inputs: it must be a cache hit (no
    # rebuild banner) with bit-identical tables.
    sig2 = np.full(ND, -0.9)
    path = write_case(
        "rh_cache.txt",
        r, v, sig, chil, etal, bbc, znet, neg,
        extra_geom=[(r, v, sig2), (r, v, sig2)],
    )
    (pairs_a, _, meta_a), (pairs_b, _, meta_b), (pairs_c, _, meta_c) = parse(
        run(binary, path).stdout
    )
    got_b = {}
    for i, j, t1, t1p, pw in pairs_b:
        got_b.setdefault(i - 1, []).append((j - 1, t1, t1p, pw))
    set_b = {(i, p[0]) for i, pl in got_b.items() for p in pl}
    set_a = {(i - 1, j - 1) for i, j, t1, t1p, pw in pairs_a}
    rep.check(
        "cache: new SIGMA keeps pair selection, triggers a rebuild",
        set_b == set_a and meta_b.get("banners", 0) >= 1,
        f"sel_same={set_b == set_a} banners={meta_b.get('banners')}",
    )
    t_ok = bool(pairs_b)
    for i, j, t1, t1p, pw in pairs_b:
        ii, jj = i - 1, j - 1
        d1 = r[jj] ** 2 - r[ii] ** 2
        d2 = r[jj] ** 2 * v[ii] ** 2 - r[ii] ** 2 * v[jj] ** 2
        mu2 = d1 * v[jj] ** 2 / d2
        mup2 = d1 * v[ii] ** 2 / d2
        if abs(t1 - abs(1.0 + mu2 * sig2[ii])) > 1e-12:
            t_ok = False
        if abs(t1p - abs(1.0 + mup2 * sig2[jj])) > 1e-12:
            t_ok = False
    rep.check(
        "cache: rebuilt T1/T1P use the new SIGMA, not the stale tables",
        t_ok,
    )
    rep.check(
        "cache: identical inputs reuse tables (no rebuild, same values)",
        meta_c.get("banners", 0) == 0 and pairs_c == pairs_b,
        f"banners={meta_c.get('banners')}",
    )

    # -- T8: uniform CHI -> analytic chord attenuation -------------
    # With CHI=c0 the trapezoid telescopes to c0*|s_i - s_j| on the
    # (same-leg) admitted pairs, so ATT has a closed form.
    chi_c = np.full(ND, 1.0e-4)
    path = write_case(
        "rh_att_uniform.txt", r, v, sig, chil, etal, bbc, znet, neg, chi=chi_c
    )
    pairs_a, drc_a, meta_a = parse(run(binary, path).stdout)
    att_map = meta_a.get("att", {})
    att_ok = bool(att_map)
    for i, j, t1, t1p, pw in pairs_a:
        ii, jj = i - 1, j - 1
        k = [n for n, (pi, pj, *_x) in enumerate(
            [p for p in pairs_a if p[0] == i]) if pj == j][0] + 1
        d1 = r[jj] ** 2 - r[ii] ** 2
        d2 = r[jj] ** 2 * v[ii] ** 2 - r[ii] ** 2 * v[jj] ** 2
        mu2 = d1 * v[jj] ** 2 / d2
        p = r[ii] * np.sqrt(1.0 - mu2)
        s_i = np.sqrt((r[ii] - p) * (r[ii] + p))
        s_j = np.sqrt((r[jj] - p) * (r[jj] + p))
        att_exp = np.exp(-min(1.0e-4 * abs(s_i - s_j), 700.0))
        if abs(att_map[i, k] - att_exp) > 1e-9:
            att_ok = False
    dexp_a = drc_expected(
        exp_pairs, r, v, chil, etal, bbc, neg, FL, chi=chi_c
    )
    dgot_a = np.zeros_like(chil)
    for (l, i), x in drc_a.items():
        dgot_a[i - 1, l - 1] = x
    err_a = np.abs(dgot_a - dexp_a).max() / max(1.0, np.abs(dexp_a).max())
    rep.check(
        "uniform CHI: ATT=exp(-c0*chord) analytic; DRC consistent",
        att_ok and err_a < 1e-12,
        f"att_ok={att_ok} relerr={err_a:.2e}",
    )

    # -- T9: opaque CHI kills the remote term; EXTINCT=F restores --
    chi_big = np.full(ND, 1.0e-1)
    path = write_case(
        "rh_att_big.txt", r, v, sig, chil, etal, bbc, znet, neg, chi=chi_big
    )
    _, drc_big, _ = parse(run(binary, path).stdout)
    path = write_case(
        "rh_att_off.txt", r, v, sig, chil, etal, bbc, znet, neg,
        chi=chi_big, extinct=False,
    )
    _, drc_off, _ = parse(run(binary, path).stdout)
    dexp_off = drc_expected(exp_pairs, r, v, chil, etal, bbc, neg, FL)
    dgot_off = np.zeros_like(chil)
    for (l, i), x in drc_off.items():
        dgot_off[i - 1, l - 1] = x
    err_off = np.abs(dgot_off - dexp_off).max() / max(1.0, np.abs(dexp_off).max())
    rep.check(
        "opaque CHI: DRC~0; EXTINCT=F recovers unobstructed kernel",
        max(abs(x) for x in drc_big.values()) < 1e-30 and err_off < 1e-12,
        f"max|DRC|={max(abs(x) for x in drc_big.values()):.2e} "
        f"relerr_off={err_off:.2e}",
    )

    # -- T10: opaque shell extinguishes only crossing pairs --------
    chi_sh = np.zeros(ND)
    chi_sh[10] = 1.0e2  # one opaque cell in the middle of the grid
    chi_sh[11] = 1.0e2
    path = write_case(
        "rh_att_shell.txt", r, v, sig, chil, etal, bbc, znet, neg, chi=chi_sh
    )
    pairs_s, _, meta_s = parse(run(binary, path).stdout)
    att_map = meta_s.get("att", {})
    sh_ok = bool(att_map)
    for n, (i, j, t1, t1p, pw) in enumerate(pairs_s):
        ii, jj = i - 1, j - 1
        k = [m for m, (pi, pj, *_x) in enumerate(
            [p for p in pairs_s if p[0] == i]) if pj == j][0] + 1
        # the trapezoid nodes include both endpoint cells, so an
        # opaque cell attenuates whenever it lies in [min, max]
        crosses = min(ii, jj) <= 11 and max(ii, jj) >= 10
        got = att_map[i, k]
        if crosses and got > 1e-6:
            sh_ok = False
        if not crosses and abs(got - 1.0) > 1e-12:
            sh_ok = False
    rep.check("opaque shell: only chord-crossing pairs attenuated", sh_ok)

    # -- Plot: pair map + DRC profile ------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    m = np.zeros((ND, ND))
    for i, j in set_got:
        m[i, j] = 1
    ax[0].imshow(m, origin="lower", cmap="Blues", extent=[0, ND - 1, 0, ND - 1])
    ax[0].set_xlabel("partner depth J")
    ax[0].set_ylabel("target depth I")
    ax[0].set_title("accepted RH78 pairs (dip profile)")
    ax[1].semilogy(r, np.abs(v) / 1e3, "k-", label="|V|")
    for l in range(NSIM):
        ax[1].semilogy(
            r, dexp[:, l] / dexp.max() * 10, "--", label=f"DRC line {l} (scaled)"
        )
    ax[1].set_xscale("log")
    ax[1].set_xlabel("R")
    ax[1].legend()
    fig.tight_layout()
    out = PLOTDIR / "sob_rh_couple.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  plot -> {out}")

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
