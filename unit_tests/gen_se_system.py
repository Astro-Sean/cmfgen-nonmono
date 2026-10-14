#!/usr/bin/env python3
"""Generate surrogate statistical-equilibrium systems for JFNK unit tests.

The surrogate reproduces the frozen-J residual form used by
JFNK_RESIDUAL_V1 (jfnk_cmfgen_interface_v1.f):

    F_i(x) = sum_j A_ij * n0_j * exp(x_j)  +  sum_j C_ij * x_j  +  r_i

where x_j = ln(n_j / n0_j), A is the collisional rate matrix (the only
nonlinearity is exp(x), matching the log-unknown map), C is the
CMAT_RAD linear response (block-banded across depths), and r holds the
frozen radiative/source remainder.

Regimes reproduced, per depth block (NT=37 rows):

    rows  0-23  normal level rows; loss rates lam_i ~ 1e-2..1e9 s^-1,
                gains from 1-3 donor levels (stiff spectrum like SE)
    rows 24-29  rydberg chain: each row = q * previous row + eps_vec,
                eps/|row| ~ 1e-9 -> near-proportional rows, i.e. near-
                null singular directions like high-n HYD levels where
                detailed balance nearly equates adjacent rate equations
    rows 30-32  floor rows: n0 = 1e-30, weak sink -> rowmax << blockmax
                (BA row is null -> JFNK_ROW_REG candidate); equilibrium
                sits above the floor so x* ~ +46 (demand e46)
    row  33     ionization/conservation row: F = sum n_group - N_ion
    row  34     const row (hard spec only): A and C rows identically
                zero but r != 0 -> structurally unsolvable, like a BA
                row zeroed by the NaN path with a finite STEQ residual
    row  35     T equation - fixed unknown (masked, F forced to 0)
    row  36     Ne equation - fixed unknown (masked)

Outputs per spec directory (plain text, list-directed readable):
    params.txt      nt nd n
    n0.txt          n0(1..n)
    acoll.txt       A^T written row-major (Fortran reads A col-major)
    cmat.txt        C^T likewise
    rvec.txt        r
    xstar.txt       x* (F(x*) = 0 except const rows)
    mask.txt        1 = fixed row (T/Ne), 0 = free
    escale.txt      row scale = max_j |J0_ij|  (DGEEQU-style, for the
                    equation-scaled variant)
    meta.json       row classes, seeds, provenance
"""

import json
import sys
from pathlib import Path

import numpy as np

NT = 37
N_NORMAL = 24
RYD = list(range(24, 30))
FLOOR = list(range(30, 33))
ION_ROW = 33
CONST_ROW = 34
T_ROW = 35
NE_ROW = 36

X_FLOOR_DEMAND = np.log(1e-10 / 1e-30)  # ~46.05


def build_spec(nd, seed, hard, demand_scale=1.0, ryd_eps=1e-9,
               radpin_frac=0.0, radpin_jfac=1e-6, radpin_fmag=1e2,
               tstar=0.0, lam_max=9.0, n0_hi=6.0, ryd_hi=9.0,
               flip=False, macro=False):
    rng = np.random.default_rng(seed)
    n = NT * nd
    a = np.zeros((n, n))
    c = np.zeros((n, n))
    n0 = np.zeros(n)
    mask = np.zeros(n, dtype=int)
    row_class = [""] * n

    lam = 10.0 ** rng.uniform(-2.0, lam_max, N_NORMAL)
    lam[0] = 10.0 ** lam_max  # designated stiff row (HYD lev-23 analogue)

    for k in range(nd):
        base = k * NT

        # --- normal level rows -------------------------------------
        for l in range(N_NORMAL):
            i = base + l
            a[i, i] = -lam[l]
            donors = rng.choice(N_NORMAL, size=3, replace=False)
            for j in donors:
                if j != l:
                    a[i, base + j] += 10.0 ** rng.uniform(-3.0, 6.0)
            n0[i] = 10.0 ** rng.uniform(-8.0, n0_hi)
            row_class[i] = "lev"
        n0[base + 0] = min(1e4, 10.0 ** n0_hi)

        # --- rydberg chain: near-proportional rows ------------------
        prev = None
        for l in RYD:
            i = base + l
            if prev is None:
                a[i, i] = -10.0 ** rng.uniform(ryd_hi - 3.0, ryd_hi)
                for j in rng.choice(N_NORMAL, 3, replace=False):
                    a[i, base + j] += 10.0 ** rng.uniform(
                        0.0, min(5.0, ryd_hi - 1.0))
            else:
                q = rng.uniform(0.7, 1.4)
                a[i, :] = q * a[prev, :]
                a[i, :] += ryd_eps * np.abs(a[prev, :]).max() * \
                    rng.uniform(-1.0, 1.0, n)
            n0[i] = 10.0 ** rng.uniform(-14.0, -2.0)
            row_class[i] = "ryd"
            prev = i

        # --- floor rows ---------------------------------------------
        for l in FLOOR:
            i = base + l
            a[i, i] = -1e-2
            a[i, base + N_NORMAL - 1] = 1e-7
            n0[i] = 1e-30
            row_class[i] = "floor"

        # --- ionization/conservation row ----------------------------
        i = base + ION_ROW
        group = list(range(base, base + 8))
        a[i, group] = 1.0
        n0[i] = 1.0
        row_class[i] = "ion"

        # --- const / normal spare row -------------------------------
        i = base + CONST_ROW
        n0[i] = 1.0
        row_class[i] = "const" if hard else "lev"
        if not hard:
            a[i, i] = -10.0 ** rng.uniform(-2.0, 2.0)
            a[i, base + 1] = 1e-3

        # --- fixed T / Ne rows --------------------------------------
        for l in (T_ROW, NE_ROW):
            i = base + l
            a[i, i] = -1e6
            a[i, base + 0] = 1e2
            n0[i] = 1.0
            mask[i] = 1
            row_class[i] = "T" if l == T_ROW else "Ne"

        # --- CMAT_RAD: diagonal blocks + nearest-depth band ---------
        for l in list(range(N_NORMAL)) + RYD:
            i = base + l
            c[i, i] = -10.0 ** rng.uniform(-1.0, 4.0)
        for l in FLOOR:
            c[base + l, base + l] = -1e-8
        if k < nd - 1:
            nxt = (k + 1) * NT
            for l in range(0, N_NORMAL, 3):
                c[base + l, nxt + l] = 10.0 ** rng.uniform(-2.0, 1.0)
                c[nxt + l, base + l] = 10.0 ** rng.uniform(-2.0, 1.0)

    # --- exact solution x* and constant vector r --------------------
    xstar = np.zeros(n)
    for k in range(nd):
        base = k * NT
        xstar[base : base + N_NORMAL] = (
            rng.normal(0.0, 0.5 * np.log(1e2), N_NORMAL) * demand_scale
        )
        xstar[base + RYD[0] : base + RYD[-1] + 1] = (
            rng.normal(0.0, 0.5 * np.log(1e4), len(RYD)) * demand_scale
        )
        xstar[base + FLOOR[0] : base + FLOOR[-1] + 1] = X_FLOOR_DEMAND
        xstar[base + 3] = -20.0  # overpopulated level
        # tstar != 0 is the STELLA shock: the equilibrium T row sits
        # at x*_T = tstar so F_T(0) ~ A_TT*n0_T*(1 - e^tstar) ~ 2e7,
        # a massive energy residual on a masked-but-live T row.
        xstar[base + T_ROW] = tstar
        xstar[base + NE_ROW] = 0.0
    xstar = np.clip(xstar, -40.0, 46.5)
    nstar = n0 * np.exp(np.clip(xstar, -80, 80))
    r = -a @ nstar - c @ xstar

    if hard:
        # structurally-null rows with finite residual (unsolvable)
        for k in range(nd):
            i = k * NT + CONST_ROW
            a[i, :] = 0.0
            c[i, :] = 0.0
            r[i] = 1e6 * (1.0 + k)

    if radpin_frac > 0.0:
        # Radiation-pinned rows (production signature): the row residual
        # is dominated by the frozen source r_i while the Jacobian row
        # retains only a weak population coupling, so F is nearly
        # orthogonal to range(J) and ||F|| is unreachable by any step.
        # Mirrors the SOB frozen-J residual where radiative terms pin
        # the balance on most rows.
        n_pin = int(round(radpin_frac * N_NORMAL))
        for k in range(nd):
            base = k * NT
            pin = rng.choice(N_NORMAL, size=n_pin, replace=False)
            for l in pin:
                i = base + l
                jrow = a[i, :] * n0[None, :] + c[i, :]
                a[i, :] *= radpin_jfac
                c[i, :] *= radpin_jfac
                fscale = radpin_fmag * (np.abs(r[i]) if r[i] != 0.0
                                        else np.abs(jrow).max())
                r[i] = fscale * rng.choice([-1.0, 1.0])
                row_class[i] = "pin"

    j0 = a * n0[None, :] + c  # exact x-space Jacobian at baseline
    escale = np.abs(j0).max(axis=1)
    escale[escale == 0.0] = 1.0

    flip_spec = None
    if flip:
        # REPLACE analogue: row irow carries its normal rate equation
        # while x[icol] <= thr, and a conservation form
        # sum(w n) - tgt above it.  thr sits just above x[icol]=0 so
        # a single FD step on icol crosses the flip; tgt is offset
        # from the baseline-consistent value by a fraction of the
        # row scale so the two branches disagree by O(F_row) at the
        # crossing, like a conservation-vs-rate residual swap.
        irow, icol = 7, 5
        grp = np.array([5, 6, 7])
        w = np.ones(len(grp))
        gap = 1e-3 * escale[irow]
        flip_spec = dict(
            irow=irow, icol=icol, thr=1e-9,
            tgt=float((w * n0[grp]).sum() + gap),
            grp=grp.tolist(), w=w.tolist(),
        )

    factor_spec = None
    if macro:
        grp = np.array([5, 6, 7])
        target = float(n0[grp].sum())
        factor_spec = dict(
            alpha=0.5 / target, target=target,
            grp=(grp + 1).tolist(),
        )

    return dict(
        a=a, c=c, r=r, n0=n0, xstar=xstar, mask=mask,
        row_class=row_class, nt=NT, nd=nd, n=n,
        seed=seed, hard=hard, escale=escale, flip=flip_spec,
        factor=factor_spec,
    )


def write_spec(spec, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n, nt, nd = spec["n"], spec["nt"], spec["nd"]
    np.savetxt(outdir / "params.txt", np.array([[nt, nd, n]], dtype=int),
               fmt="%d")
    np.savetxt(outdir / "n0.txt", spec["n0"], fmt="%.16e")
    np.savetxt(outdir / "acoll.txt", spec["a"].T, fmt="%.16e")
    np.savetxt(outdir / "cmat.txt", spec["c"].T, fmt="%.16e")
    np.savetxt(outdir / "rvec.txt", spec["r"], fmt="%.16e")
    np.savetxt(outdir / "xstar.txt", spec["xstar"], fmt="%.16e")
    np.savetxt(outdir / "mask.txt", spec["mask"], fmt="%d")
    np.savetxt(outdir / "escale.txt", spec["escale"], fmt="%.16e")
    if spec.get("flip") is not None:
        fl = spec["flip"]
        with open(outdir / "flip.txt", "w") as fh:
            fh.write(
                f"{fl['irow'] + 1} {fl['icol'] + 1} "
                f"{fl['thr']:.16e} {fl['tgt']:.16e}\n"
            )
            fh.write(f"{len(fl['grp'])}\n")
            for j, w in zip(fl["grp"], fl["w"]):
                fh.write(f"{j + 1} {w:.16e}\n")
    if spec.get("factor") is not None:
        fac = spec["factor"]
        with open(outdir / "factor.txt", "w") as fh:
            fh.write(
                f"{fac['alpha']:.16e} {fac['target']:.16e} "
                f"{len(fac['grp'])}\n"
            )
            fh.write(" ".join(str(j) for j in fac["grp"]) + "\n")
    meta = dict(
        nt=nt, nd=nd, n=n, seed=spec["seed"], hard=spec["hard"],
        row_class=spec["row_class"],
        rows=dict(
            normal=[0, N_NORMAL - 1], rydberg=RYD, floor=FLOOR,
            ion=ION_ROW, const=CONST_ROW, T=T_ROW, Ne=NE_ROW,
        ),
    )
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2))


def main():
    outroot = Path(__file__).resolve().parent / "systems"
    build = {"mild": dict(nd=4, seed=11, hard=False, demand_scale=0.3,
                          ryd_eps=1e-4),
             "std": dict(nd=4, seed=7, hard=False, demand_scale=1.0,
                         ryd_eps=1e-6),
             "hard": dict(nd=4, seed=3, hard=True, demand_scale=1.0,
                          ryd_eps=1e-9),
             # production-like: ~2/3 of rows radiation-pinned, so most
             # of ||F|| is left-null and no population step reaches it
             "radpin": dict(nd=4, seed=5, hard=True, demand_scale=1.0,
                            ryd_eps=1e-9, radpin_frac=0.66,
                            radpin_jfac=1e-6, radpin_fmag=1e2),
             # STELLA-shock mock: tame population spectrum (converges
             # in a few Newton its, so RE is the dominant unsolved
             # part like production's coll-dominated residual) plus an
             # energy residual ~2e7 needing x*_T = 3 (a factor-e3 T
             # correction); the masked T row is where FIX_T_AUTO
             # damping must keep the step bounded.
             "shock": dict(nd=4, seed=13, hard=False, demand_scale=0.3,
                           ryd_eps=1e-1, tstar=3.0, lam_max=2.0,
                           n0_hi=2.0, ryd_hi=3.0),
             # 'std' matrices plus a structure-flip row (REPLACE
             # analogue) for the equation-freeze regression test
             "flip": dict(nd=4, seed=7, hard=False, demand_scale=1.0,
                          ryd_eps=1e-6, flip=True),
             # Factor-dependent residual for the lagged implicit-factor
             # test; the baseline BA omits the aggregate-factor derivative.
             "lag": dict(nd=4, seed=7, hard=False, demand_scale=1.0,
                         ryd_eps=1e-6, macro=True)}
    for name, kw in build.items():
        write_spec(build_spec(**kw), outroot / name)
        print("wrote", outroot / name)


if __name__ == "__main__":
    sys.exit(main())
