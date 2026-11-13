#!/usr/bin/env python3
"""Focused checks for OBS_FRAME_SUB_V9.

The tests run the exact production observer-frame routine on a compact
model.  A zero-emissivity optically thin envelope isolates the diffusion
inner-boundary seed.  A reversing velocity field exercises the enlarged
ray-storage bound used for non-monotonic models.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import F90DIR, PLOTDIR, TREE, WORKDIR, Report, compile_driver, run

SOURCES = [
    TREE / "tools" / "set_kind_module.f",
    TREE / "subs" / "phys_con.f",
    TREE / "subs" / "mon_interp_fast.f",
    TREE / "subs" / "open_dir_acc_v1.f",
    TREE / "subs" / "write_direct_info_v3.f",
    TREE / "tools" / "chk_little_endian.f",
    TREE / "tools" / "get_lu.f",
    TREE / "unix" / "term_io.f",
    TREE / "unix" / "dir_acc_pars_gen.f",
    TREE / "unix" / "tune.f",
    TREE / "obs" / "obs_frame_sub_v9.f",
    F90DIR / "drv_obs_frame_ib.f",
]


def main():
    rep = Report("OBS frame diffusion boundary and non-mono rays")
    out_dir = WORKDIR / "obs_frame_ib"
    out_dir.mkdir(parents=True, exist_ok=True)
    binary = compile_driver("drv_obs_frame_ib", SOURCES)
    run(binary, out_dir, cwd=out_dir)

    no_ib = np.loadtxt(out_dir / "flux_no_ib.txt")
    ib = np.loadtxt(out_dir / "flux_ib.txt")
    ib2 = np.loadtxt(out_dir / "flux_ib2.txt")
    nonmono = np.loadtxt(out_dir / "flux_nonmono.txt")
    ib_etaz = np.loadtxt(out_dir / "flux_ib_etaz.txt")
    ib2_etaz = np.loadtxt(out_dir / "flux_ib2_etaz.txt")
    thick_stau = np.loadtxt(out_dir / "flux_thick_stau.txt")
    thick_etaz = np.loadtxt(out_dir / "flux_thick_etaz.txt")
    thick_etaz_hi = np.loadtxt(out_dir / "flux_thick_etaz_hi.txt")
    nonmono_etaz = np.loadtxt(out_dir / "flux_nonmono_etaz.txt")
    negchi_stau = np.loadtxt(out_dir / "flux_negchi_stau.txt")

    rep.check(
        "stock hollow-core path gives zero flux for zero emissivity",
        np.all(np.isfinite(no_ib)) and np.allclose(no_ib, 0.0, atol=1e-20),
        f"flux={no_ib}",
    )
    rep.check(
        "diffusion boundary injects finite positive flux",
        np.all(np.isfinite(ib)) and np.all(ib > 0.0),
        f"flux={ib}",
    )
    rep.check(
        "I(mu)=J+3muH seed scales linearly with boundary J and H",
        np.allclose(ib2, 2.0 * ib, rtol=1.0e-10, atol=0.0),
        f"ratio={ib2 / ib}",
    )
    rep.check(
        "reversing velocity model completes with finite flux",
        np.all(np.isfinite(nonmono)) and np.all(nonmono >= 0.0),
        f"flux={nonmono}",
    )
    rep.check(
        "ETAZ injects the same diffusion boundary seed as STAU",
        np.all(ib_etaz > 0.0) and np.allclose(ib_etaz, ib, rtol=1.0e-6),
        f"stau={ib} etaz={ib_etaz}",
    )
    rep.check(
        "ETAZ boundary seed scales linearly with boundary moments",
        np.allclose(ib2_etaz, 2.0 * ib_etaz, rtol=1.0e-10, atol=0.0),
        f"ratio={ib2_etaz / ib_etaz}",
    )
    rep.check(
        "STAU and ETAZ agree within quadrature error on thick model",
        np.allclose(thick_etaz, thick_stau, rtol=0.15),
        f"stau={thick_stau} etaz={thick_etaz} "
        f"rel={np.abs(thick_etaz / thick_stau - 1.0)}",
    )
    rep.check(
        "ETAZ approaches the STAU result at higher ray resolution",
        np.all(np.abs(thick_etaz_hi - thick_stau) < np.abs(thick_etaz - thick_stau)),
        f"lo_rel={np.abs(thick_etaz / thick_stau - 1.0)} "
        f"hi_rel={np.abs(thick_etaz_hi / thick_stau - 1.0)}",
    )
    rep.check(
        "reversing velocity model completes with finite flux under ETAZ",
        np.all(np.isfinite(nonmono_etaz)) and np.all(nonmono_etaz >= 0.0),
        f"flux={nonmono_etaz}",
    )
    rep.check(
        "negative-opacity cell forces STAU-to-ETAZ fallback cleanly",
        np.all(np.isfinite(negchi_stau)),
        f"flux={negchi_stau}",
    )

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    x = np.arange(no_ib.size)
    ax.plot(x, no_ib, "o-", label="hollow-core stock path")
    ax.plot(x, ib, "s-", label="J+3muH boundary (STAU)")
    ax.plot(x, ib_etaz, "v--", label="J+3muH boundary (ETAZ)")
    ax.plot(x, ib2, "d--", label="2 x boundary moments")
    ax.plot(x, thick_stau, "s:", label="thick model (STAU)")
    ax.plot(x, thick_etaz, "v:", label="thick model (ETAZ)")
    ax.plot(x, nonmono, "^-", label="non-monotonic ray")
    ax.set_xlabel("observer frequency index")
    ax.set_ylabel("formal flux")
    ax.legend()
    fig.savefig(PLOTDIR / "obs_frame_ib.png", dpi=160)
    plt.close(fig)

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
