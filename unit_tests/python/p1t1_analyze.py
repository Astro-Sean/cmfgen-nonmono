#!/usr/bin/env python3
"""P1-T1 field-sensitivity analysis (absolute-gain form).

Re-marches the production KNOP short-characteristics formal solve from the
KN_DBG_* coefficient dump of the p1t1 manufactured case, validates against
the production I_M/I_P (IMP lines), then finite-differences the lagged
I_NEXT channel cell-by-cell.

Lagged input per cell:
  compressing (Q<0): I_lag = INEXT_V  (store read; KN_DBG_INEXTM/P)
  expanding  (Q>0): I_lag = I_*_PREV (previous-frequency field)
SP - SR = |coupling| * I_lag  with coupling = |Q| (compressing) or Q (expanding).
Perturbing I_lag by an absolute delta gives SP' = SP + |Q|*delta; all other
state held fixed.  The map gain g = dI_formal/dI_NEXT is the contraction
relevant quantity; K_local = (I_lag/I_formal)*g degenerates to 0 when the
current lagged input vanishes.
"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/p1t1.out")
NDA, NCFA = 8, 6
NI = NDA


def fnum(x):
    if "E" not in x and "e" not in x:
        x = x[:-4] + "E" + x[-4:]
    return float(x)


CHIH = np.full((NDA + 1, NCFA + 1), np.nan)
Q = np.full((NDA + 1, NCFA + 1), np.nan)
DTAU = np.full((NDA + 1, NCFA + 1), np.nan)
SPM = np.full((NDA + 1, NCFA + 1), np.nan)
SPP = np.full((NDA + 1, NCFA + 1), np.nan)
SR = np.full((NDA + 1, NCFA + 1), np.nan)
IMPR = np.full((NDA + 1, NCFA + 1), np.nan)
IPPR = np.full((NDA + 1, NCFA + 1), np.nan)
INEXTM = np.full((NDA + 1, NCFA + 1), np.nan)
INEXTP = np.full((NDA + 1, NCFA + 1), np.nan)
IM = np.full((NDA + 1, NCFA + 1), np.nan)
IP = np.full((NDA + 1, NCFA + 1), np.nan)
A = np.full((5, NDA + 1, NCFA + 1), np.nan)
GAMV = {}

for ln in OUT.read_text().splitlines():
    t = ln.split()
    if not t:
        continue
    if t[0] == "DBS":
        i, l = int(t[1]), int(t[2])
        v = [fnum(x) for x in t[3:14]]
        (CHIH[i, l], Q[i, l], _, DTAU[i, l], SPM[i, l], SPP[i, l],
         _, _, IMPR[i, l], IPPR[i, l], SR[i, l]) = v
    elif t[0] == "DBSA":
        i, l = int(t[1]), int(t[2])
        for k in range(5):
            A[k, i, l] = fnum(t[3 + k])
    elif t[0] == "DBSN":
        i, l = int(t[1]), int(t[2])
        INEXTM[i, l], INEXTP[i, l] = fnum(t[3]), fnum(t[4])
    elif t[0] == "IMP":
        l, i = int(t[1]), int(t[2])
        IM[i, l], IP[i, l] = fnum(t[3]), fnum(t[4])
    elif t[0] == "GAMV":
        GAMV[int(t[1])] = fnum(t[2])


def sgn1(x):
    return np.where(x >= 0.0, 1.0, -1.0)


def build_ds(SP, dtau):
    S = np.zeros(NI + 1)
    for i in range(1, NI):
        S[i] = (SP[i + 1] - SP[i]) / dtau[i] if dtau[i] > 1e-30 else 0.0
    dS = np.zeros(NI + 1)
    td = dtau[1] + dtau[2]
    dS[1] = 0.0 if abs(td) < 1e-30 else S[1] + (S[1] - S[2]) * dtau[1] / td
    for i in range(2, NI):
        td = dtau[i - 1] + dtau[i]
        dS[i] = 0.0 if abs(td) < 1e-30 else \
            (S[i - 1] * dtau[i] + S[i] * dtau[i - 1]) / td
    if NI >= 3:
        td = dtau[NI - 2] + dtau[NI - 1]
        dS[NI] = 0.0 if abs(td) < 1e-30 else \
            S[NI - 1] + (S[NI - 1] - S[NI - 2]) * dtau[NI - 1] / td
    else:
        dS[NI] = S[NI - 1]
    dS[1] = (sgn1(S[1]) + sgn1(dS[1])) * min(abs(S[1]), 0.5 * abs(dS[1]))
    for i in range(2, NI):
        dS[i] = (sgn1(S[i - 1]) + sgn1(S[i])) * \
            min(abs(S[i - 1]), abs(S[i]), 0.5 * abs(dS[i]))
    dS[NI] = (sgn1(S[NI - 1]) + sgn1(dS[NI])) * \
        min(abs(S[NI - 1]), 0.5 * abs(dS[NI]))
    return S, dS


def floor_cap(x, cap=1e20):
    x = np.where((x < 0.0) | ~np.isfinite(x), 1e-30, x)
    return np.minimum(x, cap)


def march_inward(SP, a, dtau):
    _, dS = build_ds(SP, dtau)
    A0, A1, A2, A3, A4 = a
    I = np.zeros(NI + 1)
    for i in range(1, NI):
        I[i + 1] = I[i] * A0[i] + (SP[i] * A1[i] + SP[i + 1] * A2[i]
                                  + dS[i] * A3[i] + dS[i + 1] * A4[i])
    return floor_cap(I)


def march_outward(SP, a, dtau, icore):
    _, dS = build_ds(SP, dtau)
    A0, A1, A2, A3, A4 = a
    I = np.zeros(NI + 1)
    I[NI] = icore
    for i in range(NI - 1, 0, -1):
        I[i] = I[i + 1] * A0[i] + (SP[i + 1] * A1[i] + SP[i] * A2[i]
                                  - dS[i + 1] * A3[i] - dS[i] * A4[i])
    return floor_cap(I)


# ---- validation -------------------------------------------------------
maxrel = {"M": 0.0, "P": 0.0}
for L in range(1, NCFA + 1):
    dt = DTAU[:, L]
    a = [A[k, :, L] for k in range(5)]
    im = march_inward(SPM[:, L], a, dt)
    ip = march_outward(SPP[:, L], a, dt, im[NI])
    for i in range(1, NI + 1):
        if np.isfinite(IM[i, L]) and IM[i, L] > 0:
            maxrel["M"] = max(maxrel["M"], abs(im[i] - IM[i, L]) / IM[i, L])
        if np.isfinite(IP[i, L]) and IP[i, L] > 0:
            maxrel["P"] = max(maxrel["P"], abs(ip[i] - IP[i, L]) / IP[i, L])
print(f"VALIDATION max rel err: I_M={maxrel['M']:.3e}  I_P={maxrel['P']:.3e}")

# ---- finite difference, absolute gain ---------------------------------
EPS = [1e-3, 1e-4]
rows = []
for L in range(1, NCFA + 1):
    dt = DTAU[:, L]
    a = [A[k, :, L] for k in range(5)]
    im0 = march_inward(SPM[:, L], a, dt)          # I_M(NI) -> I_CORE
    for direction, spkey, lagkey, base in (
            ("M", SPM, INEXTM if False else IMPR, im0),
            ("P", SPP, IPPR, None)):
        sp = spkey[:, L].copy()
        sr = SR[:, L]
        q = Q[:, L]
        # outward base needs I_CORE = I_M(NI)
        if direction == "M":
            bse = march_inward(sp, a, dt)
            icore = bse[NI]
        else:
            icore = im0[NI]
            bse = march_outward(sp, a, dt, icore)
        for i in range(1, NI + 1):
            qi = q[i]
            if not np.isfinite(qi) or abs(qi) < 1e-9:
                continue
            coup = abs(qi)
            # lagged input: compressing -> store read; expanding -> I_prev
            if qi < 0:
                ilag = INEXTM[i, L] if direction == "M" else INEXTP[i, L]
                branch_fired = np.isfinite(ilag) and ilag >= 0.0
            else:
                ilag = IMPR[i, L] if direction == "M" else IPPR[i, L]
                branch_fired = np.isfinite(ilag) and ilag > 0.0
            i_formal = bse[i]
            if i_formal <= 0 or not np.isfinite(i_formal):
                continue
            scale = max(abs(ilag) if np.isfinite(ilag) else 0.0,
                        abs(i_formal), 1e-30)
            ders = []
            for eps in EPS:
                d = eps * scale
                spp_ = sp.copy(); spm_ = sp.copy()
                spp_[i] = sp[i] + coup * d
                spm_[i] = sp[i] - coup * d
                if direction == "M":
                    bp = march_inward(spp_, a, dt)
                    bm = march_inward(spm_, a, dt)
                else:
                    bp = march_outward(spp_, a, dt, icore)
                    bm = march_outward(spm_, a, dt, icore)
                ders.append((bp[i] - bm[i]) / (2.0 * d))
            g = 0.5 * (ders[0] + ders[1])
            rel = abs(ders[0] - ders[1]) / max(abs(ders[0]), abs(ders[1]),
                                              1e-300)
            kloc = (ilag / i_formal) * g if (np.isfinite(ilag)
                    and ilag != 0.0) else np.nan
            rows.append(dict(L=L, dir=direction, I=i, Q=qi, absQ=abs(qi),
                             branch=int(branch_fired), ilag=ilag,
                             iformal=i_formal, gain=g, kloc=kloc,
                             epscons=rel, dtau=dt[i], chih=CHIH[i, L],
                             gam=GAMV.get(i, np.nan)))

Path("/tmp/p1t1_fd.json").write_text(json.dumps(rows, indent=1))
print(f"FD cells: {len(rows)}")
for lo, hi, tag in ((0, 0.05, "|Q|<0.05"), (0.05, 0.5, "0.05-0.5"),
                    (0.5, 0.9, "0.5-0.9"), (0.9, 2.0, "|Q|>0.9")):
    g = [r for r in rows if lo <= r["absQ"] < hi]
    if not g:
        continue
    gg = np.array([abs(r["gain"]) for r in g])
    print(f"{tag:10s} n={len(g):3d}  |gain| med={np.median(gg):.3e} "
          f"max={gg.max():.3e}  epscons_max={max(r['epscons'] for r in g):.2e}")
print("\n-- compressing (Q<0) cells --")
for r in rows:
    if r["Q"] < 0:
        print(f"  {r['dir']} L{r['L']} I{r['I']} Q={r['Q']:.3f} "
              f"branch={r['branch']} ilag={r['ilag']:.3e} "
              f"if={r['iformal']:.3e} gain={r['gain']:.3e} "
              f"K={r['kloc'] if np.isfinite(r['kloc']) else float('nan'):.3e}")
