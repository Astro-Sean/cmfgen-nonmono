#!/usr/bin/env python3
"""Unit tests for ADD_SOB_LINE_PROF (embedded in obs/cmf_flux_sub_v5.f).

Covers the guards added 16-Oct-2026:
  * the stored Doppler core integrates to CHIL*CLUMP over frequency
    (cell-width normalization),
  * a line whose rest frequency falls outside the CMF NU grid is
    skipped with a named diagnostic (not silently dropped),
  * a line masing at >10% of depths triggers the per-line warning and
    is counted in N_DROP / DROP_SUM,
  * a non-monotonic NU aborts before any deposition.

compiles the exact production file:
  obs/cmf_flux_sub_v5.f
with unit_tests/f90/drv_sob_prof.f90, using the prebuilt module files
in mod/ and the production archives in lib/ (the routine's host file
references the full cmf_flux call tree, so the whole library set must
resolve at link time even though the driver only calls one entry
point).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import F90DIR, TREE, Report, compile_driver, run

TREE_SRC = [
    # Module sources whose data cmf_flux_sub_v5 references; their
    # .mod files land in the isolated bin/mod, not the build's mod/.
    TREE / "tools" / "set_kind_module.f",
    TREE / "dust" / "dust_mod.f",
    TREE / "new_main" / "mod_subs" / "mod_lev_dis_blk.f",
    TREE / "new_main" / "mod_subs" / "eddfac_rec_defs_mod.f",
    TREE / "obs" / "cmf_flux_cntrl_var_mod.f",
    TREE / "obs" / "mod_cmf_obs.f",
    TREE / "obs" / "mod_freq_obs.f",
    TREE / "obs" / "cmf_flux_sub_v5.f",
    # Siblings of the host file that its unexercised paths still
    # reference at link time.
    TREE / "obs" / "rd_cmf_flux_controls.f",
    TREE / "obs" / "ins_line_v5.f",
    TREE / "obs" / "ins_line_obs_v4.f",
    TREE / "obs" / "rd_cont_j_obs.f",
    TREE / "obs" / "revise_obs_p.f",
    TREE / "obs" / "obs_frame_sub_v9.f",
    TREE / "obs" / "comp_ew_abs_ew.f",
    TREE / "obs" / "check_cmf_flux_param_consis.f",
    F90DIR / "drv_sob_prof.f90",
    # Archives as positional link inputs (must follow the objects).
    TREE / "lib" / "libsubs.a",
    TREE / "lib" / "libdisp.a",
    TREE / "lib" / "libmod.a",
    TREE / "lib" / "libpgplt.a",
    TREE / "lib" / "libpgplot_stubs.a",
    TREE / "lib" / "libmy_lapack.a",
    TREE / "lib" / "libmy_blas.a",
    TREE / "lib" / "libspec.a",
]


def get_bin():
    return compile_driver(
        "drv_sob_prof",
        TREE_SRC,
        extra_flags=[f"-I{TREE}/mod"],
    )


def field(out, tag, idx):
    for ln in out.splitlines():
        t = ln.split()
        if t and t[0] == tag:
            return float(t[idx])
    return None


def main():
    rep = Report("ADD_SOB_LINE_PROF (cmf_flux SOB injection)")
    bin_ = get_bin()

    # -- T1: deposition normalization.  int chi(nu) dnu must equal
    # CHIL*CLUMP*1e-15 within a few percent (Gaussian truncated at
    # +/-3 Doppler widths, so ~0.1% is lost by design).
    p = run(bin_, "mono")
    integ = field(p.stdout, "CHI_INTEG_EXPECT", 1)
    expect = field(p.stdout, "CHI_INTEG_EXPECT", 2)
    rep.check(
        "T1 stored profile integrates to CHIL*CLUMP",
        integ is not None and abs(integ / expect - 1.0) < 0.02,
        f"integ={integ:.4e} expect={expect:.4e}",
    )
    rep.check(
        "T1b no cells dropped on a clean line",
        field(p.stdout, "N_DROP", 1) == 0.0,
    )

    # -- T2: off-grid line is skipped and named in the log.
    p = run(bin_, "offgrid")
    rep.check(
        "T2 off-grid SOB line named and skipped",
        "outside CMF grid" in p.stdout and field(p.stdout, "CHI_SUM", 1) == 0.0,
    )

    # -- T3: masing line drops >10% of depths -> per-line warning,
    #    counters consistent.
    p = run(bin_, "maser")
    rep.check(
        "T3 masing line warns and counts drops",
        "masing/neg-opacity" in p.stdout and field(p.stdout, "N_DROP", 1) == 8.0,
        f"N_DROP={field(p.stdout, 'N_DROP', 1)}",
    )

    # -- T4: non-monotonic NU aborts before deposition (A2 guard).
    # Production code uses bare STOP, which exits rc=0; assert on the
    # message plus the absence of the post-call output line.
    p = run(bin_, "nonmono", check=False)
    rep.check(
        "T4 non-monotonic NU aborts",
        "not strictly" in p.stdout and "N_DROP" not in p.stdout,
        f"rc={p.returncode}",
    )

    return rep.summary()


if __name__ == "__main__":
    sys.exit(main())
