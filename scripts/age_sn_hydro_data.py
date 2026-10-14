#!/usr/bin/env python3
"""Age a CMFGEN SN_HYDRO_DATA snapshot to a later epoch.

Maps an existing converted snapshot onto an older-but-larger ejecta by
ballistic expansion (R -> R + V*dt, the same mapping
set_rv_hydro_model_v3.f applies when SN_AGE > header age), diluting
rho/ne/na by the Lagrangian volume factor and cooling T as T ~ 1/t
(radiation-dominated homologous expansion).  The velocity is
monotonised with prepare_VADAT.monotonize_velocity BEFORE the mapping
so the file that reaches VEL_LAW=11 is already monotone - editing
RDINR alone does not propagate to the solver.

The Ni56 -> Co56 -> Fe56 chain is re-computed at the target age from
the recovered initial Ni56 (NICK_t = Ni56_0 * exp(-t/8.8d)).  Stable
iron is preserved by subtracting the decay-produced Fe56 from the file's
IRON block.

Usage:
    python3 scripts/age_sn_hydro_data.py SRC_SN_HYDRO_DATA OUT_SN_HYDRO_DATA \
        --target-days 350 --lum-lsun 2.6e8 [--t-floor-k 2500]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(SCRIPTS_DIR / "grid_mapping"))

from prepare_VADAT import monotonize_velocity, parse_sn_hydro
from sn_hydro_io import read_sn_hydro_data, validate_sn_hydro_data

_PROJ_ROOT = SCRIPTS_DIR.parent
if str(_PROJ_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJ_ROOT))

T_NI56_D = 8.8    # Ni56 e-folding decay time (days), matches grid_build.py
T_CO56_D = 111.3  # Co56 e-folding decay time (days)
SEC_PER_DAY = 86400.0
KAPPA_ES = 0.34   # e-scattering floor used for the header Rosseland tau


def _block(label: str, values: np.ndarray, per_line: int = 6) -> str:
    out = [f"{label}\n"]
    for i in range(0, len(values), per_line):
        chunk = values[i : i + per_line]
        out.append("".join(f"    {v:22.15E}" for v in chunk) + "\n")
    out.append("\n")
    return "".join(out)


def _iso_block(label: str, values: np.ndarray, per_line: int = 8) -> str:
    out = [f"{label}\n"]
    for i in range(0, len(values), per_line):
        chunk = values[i : i + per_line]
        out.append("".join(f"    {v:14.7E}" for v in chunk) + "\n")
    out.append("\n")
    return "".join(out)


def decay_chain(ni56_0: np.ndarray, t_days: float):
    """Ni56 -> Co56 -> Fe56 fractions at t_days from initial Ni56."""
    ni = ni56_0 * np.exp(-t_days / T_NI56_D)
    co = (
        ni56_0
        * (T_CO56_D / (T_CO56_D - T_NI56_D))
        * (np.exp(-t_days / T_CO56_D) - np.exp(-t_days / T_NI56_D))
    )
    fe = ni56_0 - ni - co
    return ni, co, fe


def ni56_initial_from_swd(
    swd_path: Path, xni_path: Path, epoch_days: float, r_cm: np.ndarray
):
    """Recover the explosion-time Ni56 mass fraction on the radius grid r_cm.

    The .swd composition columns carry no Ni56 at late epochs (the column
    is the surviving fraction, floored to zero by the writer); the binary
    .xni holds the explosion-time Ni56 per STELLA zone, written
    inner->outer in zone order.  Zone radii at the source epoch come from
    the .swd 'km' zone index, and the profile is interpolated onto the
    source-epoch radii r_cm (the Lagrangian mapping preserves zone mass
    order, so radius interpolation is exact on a shared grid).
    """
    from prepare_SN_HYDRO_DATA.utils.swd_io import read_swd, read_xni

    epochs = read_swd(str(swd_path), target_days=epoch_days)
    t_sel = next(t for t, s in epochs.items() if s is not None)
    snap = epochs[t_sel]
    nzon = len(snap["r"])
    ni0_zone = read_xni(str(xni_path), nzon)  # outer->inner after reader
    r_swd = snap["r"]                        # outer->inner
    order = np.argsort(r_swd)                # -> inner->outer
    ni0_grid = np.interp(
        r_cm, r_swd[order], ni0_zone[order], left=0.0, right=0.0
    )
    print(
        f"[ni56] .xni from t_sel={t_sel:.3f} d, nzon={nzon}, "
        f"max Ni56_0={ni0_grid.max():.3e}"
    )
    return np.clip(ni0_grid, 0.0, 1.0)


def age_sn_hydro(
    src: Path,
    out: Path,
    target_days: float,
    lum_lsun: float,
    t_floor_k: float,
    swd_path: Path = None,
    xni_path: Path = None,
) -> None:
    hydro = parse_sn_hydro(src)          # inner->outer order
    raw = read_sn_hydro_data(src)        # outer->inner order; isotopes
    t_old = hydro.params["time_days"]
    if target_days <= t_old:
        raise ValueError(
            f"target_days={target_days} must exceed source age {t_old:.3f} d"
        )
    n = int(hydro.params["ND"])
    print(f"[age] {src} t={t_old:.3f} d -> {target_days:.3f} d (ND={n})")

    hydro = monotonize_velocity(hydro)   # edits vel_cms only
    v = hydro.vel_cms                    # cm/s, inner->outer
    r = hydro.radius_cm.copy()           # cm, inner->outer
    rho = hydro.rho_cgs.copy()
    if np.any(v <= 0):
        raise ValueError("non-positive velocity after monotonization")
    if np.any(np.diff(v) <= 0):
        raise ValueError("velocity still not strictly increasing outward")

    # Ballistic expansion of every Lagrangian grid point.
    dt_s = (target_days - t_old) * SEC_PER_DAY
    r_new = r + v * dt_s
    drdr = np.gradient(r_new, r)
    if np.any(drdr <= 0):
        raise ValueError("ballistic map r->r+v*dt is not monotone")
    rho_new = rho * (r / r_new) ** 2 / drdr
    dilution = rho_new / rho
    print(
        f"[age] r: {r[0]:.3e} -> {r_new[0]:.3e} cm (inner), "
        f"{r[-1]:.3e} -> {r_new[-1]:.3e} cm (outer)"
    )
    print(f"[age] rho dilution: {dilution.min():.3e} .. {dilution.max():.3e}")

    # Radiation-dominated homologous cooling; floor keeps a sane start T
    # for the FIX_T=F bootstrap (the point source sets the real profile).
    t_scale = t_old / target_days
    t_k = None
    if hydro.temp_104k is not None:
        t_k = np.maximum(hydro.temp_104k * 1e4 * t_scale, t_floor_k)
        print(
            f"[age] T: {hydro.temp_104k.min()*1e4:.0f}-{hydro.temp_104k.max()*1e4:.0f} K "
            f"-> {t_k.min():.0f}-{t_k.max():.0f} K (T~1/t, floor {t_floor_k:.0f} K)"
        )

    ed_new = hydro.ed_cgs * dilution if hydro.ed_cgs is not None else None
    na_new = hydro.na_cgs * dilution if hydro.na_cgs is not None else None
    kappa = hydro.kappa_cgs
    if kappa is not None:
        kappa = np.maximum(kappa, KAPPA_ES)

    sigma_new = np.gradient(np.log(v), np.log(r_new)) - 1.0
    print(f"[age] sigma: {sigma_new.min():.4f} .. {sigma_new.max():.4f}")

    # --- Ni56 decay bookkeeping -------------------------------------------
    # The file's NICK elemental block is the surviving Ni56; recover the
    # explosion-time Ni56 (Ni56_0 = Ni56_t * exp(+t/8.8d)), then decay to
    # the target age.  Fe56 produced by decay is removed from IRON before
    # re-adding the new value so stable progenitor Fe is preserved.
    mf = {k: v.copy() for k, v in hydro.mass_fractions.items()}
    isotopes = {}
    if swd_path is not None and xni_path is not None:
        ni56_0 = ni56_initial_from_swd(swd_path, xni_path, t_old, r)
    else:
        nick = mf.get("NICK", np.zeros(n))
        cob = mf.get("COB", np.zeros(n))
        ni56_0 = nick * np.exp(t_old / T_NI56_D)
        # Where surviving Ni56 is too small to back-extrapolate reliably,
        # recover the initial abundance from the Co56 reservoir instead:
        # Co56_t = Ni56_0 * (tCo/(tCo-tNi)) * (e^-t/tCo - e^-t/tNi).
        co_fac = (T_CO56_D / (T_CO56_D - T_NI56_D)) * (
            np.exp(-t_old / T_CO56_D) - np.exp(-t_old / T_NI56_D)
        )
        ni_lo = (nick < 1e-8) | ~np.isfinite(ni56_0)
        if np.any(ni_lo):
            ni_from_co = np.where(co_fac > 0, cob / co_fac, 0.0)
            ni56_0 = np.where(ni_lo, ni_from_co, ni56_0)
        ni56_0 = np.clip(
            np.nan_to_num(ni56_0, nan=0.0, posinf=0.0, neginf=0.0),
            0.0,
            1.0,
        )
    ni_t, co_t, fe56_t = decay_chain(ni56_0, target_days)
    _, _, fe56_old = decay_chain(ni56_0, t_old)
    iron_old = mf.get("IRON", np.zeros(n))
    stable_fe = np.clip(iron_old - fe56_old, 0.0, 1.0)
    mf["NICK"] = ni_t
    mf["COB"] = co_t
    mf["IRON"] = stable_fe + fe56_t
    if raw.isotopes:
        isotopes[("NICK", 56)] = ni_t
        isotopes[("COB", 56)] = co_t
        isotopes[("IRON", 56)] = mf["IRON"]
    print(
        f"[age] Ni56/Co56/Fe56 at {target_days:.0f} d "
        f"(max): {ni_t.max():.3e} / {co_t.max():.3e} / {fe56_t.max():.3e}"
    )

    # Renormalise so sum(mass fractions) = 1 at every depth (decay is
    # mass-conserving, this only fixes roundoff and the stable-Fe clip).
    mf_sum = np.zeros(n)
    for arr in mf.values():
        mf_sum += arr
    mf_sum = np.where(mf_sum < 1e-30, 1.0, mf_sum)
    for k in mf:
        mf[k] = mf[k] / mf_sum
    for k in isotopes:
        isotopes[k] = isotopes[k] / mf_sum

    # Header optical depths on the new structure (integrated inward from
    # the outer edge; r increases inner->outer so trapezoid over r_new
    # gives the full column).
    chi_col = None
    if kappa is not None:
        chi_col = float(np.trapezoid(kappa * rho_new, x=r_new))
    tau_es = float(np.trapezoid(np.full(n, KAPPA_ES) * rho_new, x=r_new))
    print(f"[age] tau_es={tau_es:.3e}  tau_chi={chi_col if chi_col else float('nan'):.3e}")

    # --- write (outer -> inner order, matching rd_sn_data.f) ----------------
    flip = lambda a: a[::-1].copy() if a is not None else None
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write(f"Number of data points:{n:>10d}\n")
        f.write(f"Number of mass fractions:{len(mf):>10d}\n")
        f.write(f"Number of isotopes:{len(isotopes):>10d}\n")
        f.write(f"Time(days) since explosion:{target_days:>15.7f}\n")
        f.write(f"Energy Lum [lsun]:{lum_lsun:>19.6E}\n")
        f.write(f"Inner boundary radius [cm]:{r_new[0]:>19.6E}\n")
        if t_k is not None:
            f.write(f"Inner boundary temperature [K]:{t_k[0]:>19.6E}\n")
        f.write(f"Inner boundary Rosseland tau:{tau_es:>19.6E}\n")
        if chi_col is not None:
            f.write(f"Inner boundary chi-column tau:{chi_col:>19.6E}\n")
        f.write("\n")
        f.write(_block("Radius grid (10^10cm)", flip(r_new / 1e10)))
        f.write(_block("Velocity (km/s)", flip(v / 1e5)))
        f.write(_block("Sigma (dlnV/dlnr-1)", flip(sigma_new)))
        if t_k is not None:
            f.write(_block("Temperature (10^4 K)", flip(t_k / 1e4)))
        f.write(_block("Density (gm/cm^3)", flip(rho_new)))
        if na_new is not None:
            f.write(_block("Atom density (/cm^3)", flip(na_new)))
        if ed_new is not None:
            f.write(_block("Electron density (/cm^3)", flip(ed_new)))
        if kappa is not None:
            f.write(
                _block(
                    "Rosseland mean opacty (10^{-10} cm^{-1})",
                    flip(kappa * rho_new * 1e10),
                )
            )
            f.write(_block("Kappa (cm^2/gm)", flip(kappa)))
        for elem, frac in mf.items():
            f.write(_block(f"{elem} mass fraction", flip(frac)))
        for (iso_name, baryon), arr in isotopes.items():
            f.write(_iso_block(f"{iso_name} {baryon} mass fraction", flip(arr)))
    print(f"[age] wrote {out}")

    rep = validate_sn_hydro_data(out)
    for m in rep.errors:
        print(f"  VALIDATE ERROR: {m}")
    for m in rep.warnings:
        print(f"  VALIDATE WARN:  {m}")
    if not rep.passed:
        sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("src", type=Path, help="source SN_HYDRO_DATA")
    ap.add_argument("out", type=Path, help="output SN_HYDRO_DATA")
    ap.add_argument("--target-days", type=float, required=True)
    ap.add_argument(
        "--lum-lsun",
        type=float,
        default=0.0,
        help="Energy Lum [lsun] header (point-source luminosity)",
    )
    ap.add_argument(
        "--t-floor-k",
        type=float,
        default=2500.0,
        help="floor on the adiabatically cooled temperature [K]",
    )
    ap.add_argument(
        "--swd",
        type=Path,
        default=None,
        help="STELLA .swd of the source model (zone radii for the .xni map)",
    )
    ap.add_argument(
        "--xni",
        type=Path,
        default=None,
        help="STELLA .xni with the explosion-time Ni56 per zone",
    )
    args = ap.parse_args()
    age_sn_hydro(
        args.src,
        args.out,
        args.target_days,
        args.lum_lsun,
        args.t_floor_k,
        swd_path=args.swd,
        xni_path=args.xni,
    )


if __name__ == "__main__":
    main()
