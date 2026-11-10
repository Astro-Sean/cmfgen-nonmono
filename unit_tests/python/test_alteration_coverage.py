#!/usr/bin/env python3
"""Alteration-coverage manifest for the CMFGEN sources.

Every build input that differs from the pristine reference tree
(../cmfgen_reference/cur_cmf) must be registered in COVERAGE with the
tests that exercise it, a source-pin check, or an explicit exemption.
A new alteration that lands without registration fails this test, so
the suite enforces that changes to the CMFGEN sources stay covered.

Coverage kinds:

  exec    - named tests compile and run code from this file
            (production objects via unit_tests/common.py)
  pin     - named tests assert on this file's source text only
  gap     - shipped in an executable but no test coverage yet;
            reported loudly, does not fail (registered debt)
  exempt  - not compiled into a shipped executable, a plotting or
            standalone utility outside the cmfgen_dev.exe physics
            path, or verified by an external gate (./build.sh)

The manifest is intentionally explicit: adding a source change without
adding or updating a test means adding a COVERAGE entry here, which is
the review checkpoint for 'is this alteration legitimate'.
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import TREE, Report

REF = TREE.parent / "cmfgen_reference" / "cur_cmf"

# Only build inputs count as alterations: Fortran sources, INCLUDE
# files, and Makefiles.  Everything else (.gitignore, docs, run
# products, test trees) is out of scope for this manifest.
SOURCE_RE = re.compile(r"\.(f|f90|F|F90|inc|INC)$|(^|/)Makefile[^/]*$")

# diff -rq is pruned to the same noise the manifest ignores.
DIFF_EXCLUDES = [
    "--exclude=.git",
    "--exclude=*.o", "--exclude=*.a", "--exclude=*.mod",
    "--exclude=core*", "--exclude=*.log", "--exclude=*.out",
    "--exclude=*.err", "--exclude=*.png", "--exclude=*.csv",
    "--exclude=unit_tests", "--exclude=tests", "--exclude=docs",
    "--exclude=scripts", "--exclude=cmfgen_obsidian",
    "--exclude=.devin", "--exclude=mod", "--exclude=lib",
    "--exclude=exe*", "--exclude=bin", "--exclude=blas",
    "--exclude=dust", "--exclude=jfnk_solver_files",
    "--exclude=knop_solver_files",
]

# Directories that never hold shipped build inputs.
OUT_OF_SCOPE = (
    "unit_tests/", "tests/", "docs/", "scripts/", "cmfgen_obsidian/",
    ".devin/", ".git/", "mod/", "mod_dbg/", "lib/", "exe", "bin/",
    "blas/", "dust/", "jfnk_solver_files/", "knop_solver_files/",
)

UT = "unit_tests/python/"
TC = "tests/jfnk_core/"
TF = "tests/f90/"
UF = "unit_tests/f90/"

# path -> (kind, payload).  payload is a list of repo-relative test
# files for exec/pin, a reason string for gap/exempt.
COVERAGE = {
    # ---------------- cmfgen_dev.exe link path: executed ----------------
    "new_main/mod_subs/check_param_consistency.f": ("exec", [
        UT + "test_mom_solve.py", UT + "test_vadat_keys.py",
        UT + "test_lam_handoff.py"]),
    "new_main/mod_subs/control_variable_mod.f": ("exec", [
        UT + "test_vadat_keys.py", UT + "test_sigma_true.py",
        UT + "test_knop_ray.py"]),
    "new_main/mod_subs/iter_ledger_mod.f90": ("exec", [
        UT + "test_iter_ledger.py"]),
    "new_main/mod_subs/jfnk_cmfgen_interface_v1.f": ("exec", [
        UT + "test_jfnk_opint.py", TC + "test_schur.f90",
        TC + "test_precond_blocks.f90"]),
    "new_main/mod_subs/jfnk_cmfgen_int_mod.f": ("exec", [
        UT + "test_jfnk_opint.py"]),
    "new_main/mod_subs/jfnk_op_probe_v1.f90": ("exec", [
        UT + "test_jfnk_opint.py"]),
    "new_main/mod_subs/jfnk_solve_corr_v1.f": ("exec", [
        UT + "test_jfnk_opint.py", UT + "test_jfnk_regressions.py"]),
    "new_main/mod_subs/jfnk_solver_v1.f90": ("exec", [
        UT + "test_jfnk_solver.py", UT + "test_jfnk_opint.py",
        TC + "test_newton_core.f90", TC + "test_nm_window.f90",
        TC + "test_pc_quality.f90", TC + "test_schur.f90",
        TC + "test_gmres_edge.f90"]),
    "new_main/mod_subs/knop_status_mod.f90": ("exec", [
        UT + "test_knop_status_text.py", UT + "test_knop_ray.py"]),
    "new_main/mod_subs/lam_handoff_mod.f90": ("exec", [
        UT + "test_lam_handoff.py"]),
    "new_main/mod_subs/mod_cmfgen.f": ("exec", [
        UT + "test_sigma_true.py", UT + "test_knop_ray.py"]),
    "new_main/mod_subs/mod_lev_dis_blk.f": ("exec", [
        UT + "test_sup_to_full.py", UT + "test_sob_prof.py"]),
    "new_main/subs/cmf_blkband_pin_null_v1.f": ("exec", [
        UT + "test_cmf_blkband_pin.py"]),
    "new_main/subs/cmf_blkband_v3.f": ("exec", [
        UT + "test_mat5pen.py"]),
    "new_main/subs/do_ng_accel_band_v4.f": ("exec", [
        UT + "test_ng_abstain.py"]),
    "new_main/subs/get_edep_shock_power.f": ("exec", [
        UT + "test_edep_shock.py"]),
    "new_main/subs/recompute_sigma_true.f": ("exec", [
        UT + "test_sigma_true.py", UT + "test_velocity_profiles.py"]),
    "new_main/subs/rv_sn_model_snhyd.f": ("exec", [
        UT + "test_rv_snhyd.py"]),
    "new_main/subs/rv_sn_model_sniin.f": ("exec", [
        UT + "test_rv_sniin.py"]),
    "new_main/subs/set_rv_hydro_model_v3.f": ("exec", [
        UT + "test_rv_hydro_rdin.py"]),
    "newsubs/rd_rv_file_v2.f": ("exec", [
        UT + "test_rd_rv_file.py"]),
    "obs/cmf_flux_cntrl_var_mod.f": ("exec", [
        UT + "test_sob_prof.py", UT + "test_obs_nonmono_dispatch.py"]),
    "obs/cmf_flux_sub_v5.f": ("exec", [
        UT + "test_sob_prof.py", UT + "test_obs_nonmono_dispatch.py"]),
    "obs/mod_cmf_obs.f": ("exec", [
        UT + "test_sob_prof.py"]),
    "obs/obs_frame_sub_v9.f": ("exec", [
        UT + "test_sob_prof.py", UT + "test_obs_nonmono_dispatch.py"]),
    "obs/rd_cmf_flux_controls.f": ("exec", [
        UT + "test_sob_prof.py", UT + "test_obs_nonmono_dispatch.py"]),
    "plane/fg_j_cmf_nonmono_v1.f": ("exec", [
        UT + "test_knop_ray.py", UT + "test_mom_solve.py",
        UT + "test_fg_gs.py"]),
    "plane/get_moms_non_rel.f": ("exec", [
        UT + "test_knop_ray.py", UT + "test_mom_solve.py"]),
    "plane/mod_ray_mom_store.f": ("exec", [
        UT + "test_knop_ray.py", UT + "test_mom_solve.py"]),
    "plane/mom_j_cmf_nonmono_v1.f": ("exec", [
        UT + "test_mom_solve.py", UT + "test_mom_sweep_table.py",
        UT + "test_mom_gs.py"]),
    "pgplt/do_vec_op.f": ("exec", [
        UT + "test_sp_mon_interp.py"]),
    "pgplt/mon_interp_sp.f": ("exec", [
        UT + "test_sp_mon_interp.py"]),
    "subs/comp_obs_v2.f": ("exec", [
        UT + "test_comp_obs.py", UT + "test_comp_obs_retain.py"]),
    "subs/comp_obs_v3.f": ("exec", [
        UT + "test_comp_obs_v3_retain.py"]),
    "subs/lin_interp.f": ("exec", [
        UT + "test_interp_guards.py"]),
    "subs/mon_interp.f": ("exec", [
        UT + "test_interp_guards.py"]),
    "subs/mon_interp_fast.f": ("exec", [
        UT + "test_interp_guards.py"]),
    "subs/mon_interp_fast_v2.f": ("exec", [
        UT + "test_interp_guards.py"]),
    "subs/mon_int_funs_v2.f": ("exec", [
        UT + "test_interp_guards.py"]),
    "subs/sobew.f": ("exec", [
        UT + "test_sob_tau_floor.py"]),
    "subs/sobew_grad_v2.f": ("exec", [
        UT + "test_sob_tau_floor.py"]),
    "subs/sobew_grad_v3.f": ("exec", [
        UT + "test_sob_tau_floor.py"]),
    "subs/sobjbar_sim.f": ("exec", [
        UT + "test_sob_tau_floor.py"]),
    "subs/sob_rh_couple.f": ("exec", [
        UT + "test_sob_rh_couple.py"]),
    "subs/sup_to_full_v3b.f": ("exec", [
        UT + "test_sup_to_full.py"]),

    # ------------- link path: source-text pins only -------------------
    "new_main/cmfgen_sub.f": ("pin", [
        UT + "test_iter_conv_gate.py", UT + "test_jfnk_regressions.py",
        UT + "test_knop_status_text.py", UT + "test_lam_handoff.py"]),
    "new_main/mod_subs/comp_j_blank.f": ("pin", [
        UT + "test_knop_status_text.py"]),
    "new_main/mod_subs/rd_control_variables.f": ("pin", [
        UT + "test_vadat_keys.py", UT + "test_iter_conv_gate.py",
        UT + "test_gs_drift_diag.py", UT + "test_pop_audit.py",
        UT + "test_jfnk_regressions.py", UT + "test_knop_lagged_field.py",
        UT + "test_knop_status_text.py", UT + "test_lam_handoff.py"]),
    "new_main/mod_subs/solve_for_pops.f": ("pin", [
        UT + "test_iter_conv_gate.py", UT + "test_iter_ledger.py",
        UT + "test_gs_drift_diag.py", UT + "test_pop_audit.py",
        UT + "test_jfnk_regressions.py", UT + "test_knop_lagged_field.py",
        UT + "test_knop_status_text.py", UT + "test_lam_handoff.py"]),
    "new_main/subs/solveba_v13.f": ("pin", [
        UT + "test_knop_lagged_field.py"]),
    "obs/COMP_JCONT_V4.INC": ("pin", [
        UT + "test_obs_nonmono_dispatch.py"]),

    # ------------- link path: shipped, no test coverage yet ------------
    "new_main/cmfgen.f": ("gap",
        "main program; SIGMA_TRUE allocation runs only in full model runs"),
    "new_main/mod_subs/adjust_r_grid_v4.f": ("gap",
        "grid revision under nonmono flags; needs model state"),
    "new_main/mod_subs/comp_grey_v2.f": ("gap",
        "grey-temperature update; untested"),
    "new_main/mod_subs/comp_grey_v4.f": ("gap",
        "grey-temperature update; untested"),
    "new_main/mod_subs/do_var_cont.f": ("gap",
        "SOB/KNOP continuum-J source selection; needs model state"),
    "new_main/mod_subs/line_mod.f": ("exec",
        [UT + "test_sub_sob_line.py"]),
    "new_main/mod_subs/rd_cont_j.f": ("exec",
        [UT + "test_rd_cont_j.py"]),
    "new_main/mod_subs/set_abund_clump.f": ("exec",
        [UT + "test_vturb.py"]),
    "new_main/mod_subs/sub_sob_line_v3.f": ("exec",
        [UT + "test_sub_sob_line.py"]),
    "new_main/mod_subs/sup_to_full_v4.f": ("exec",
        [UT + "test_sup_to_full_v4.py"]),
    "new_main/subs/auto_add_ion.f": ("gap", "ion auto-add; untested"),
    "new_main/subs/check_ions_present.f": ("exec",
        [UT + "test_check_ions.py"]),
    "new_main/subs/do_tau_regrid_v2.f": ("exec",
        [UT + "test_tau_regrid.py"]),
    "new_main/subs/do_wind_vel_v1.f": ("exec",
        [UT + "test_do_wind_vel.py"]),
    "new_main/subs/eval_adiabatic_v3.f": ("gap", "adiabatic eval; untested"),
    "new_main/subs/fiddle_pop_corrections_v2.f": ("gap",
        "population-correction fiddling; untested"),
    "new_main/subs/fixpop_in_ba_v3.f": ("gap",
        "BA population fix; untested"),
    "new_main/subs/generate_full_matrix_v3.f": ("gap",
        "matrix generation; untested"),
    "new_main/subs/get_pops_at_prev_time_step_v5.f": ("exec",
        [UT + "test_get_pops_v5.py"]),
    "new_main/subs/open_rw_eddfactor.f": ("exec",
        [UT + "test_eddfactor.py"]),
    "new_main/subs/rd_sn_data.f": ("exec",
        [UT + "test_rd_sn_data.py"]),
    "new_main/subs/read_bcd_mat.f": ("exec",
        [UT + "test_bcd_mat.py"]),
    "new_main/subs/steq_advec_v4.f": ("gap",
        "statistical-equilibrium advection; untested"),
    "obs/obs_frame_v2.f": ("gap",
        "observer-frame transform V2; untested"),
    "opac/par_opac_sub.f": ("exec", [UT + "test_par_opac.py"]),
    "pgplt/set_upper_axis.f": ("exec",
        [UT + "test_set_upper_axis.py"]),
    "plane/add_dust_eta.f": ("exec", [UT + "test_add_dust.py"]),
    "plane/characteristics_v2.f": ("exec",
        [UT + "test_characteristics_v2.py"]),
    "plane/characteristics_v2_gam.f": ("exec",
        [UT + "test_characteristics_v2.py"]),
    "plane/cmf_formal_rel_v4.f": ("gap",
        "relativistic formal solver; untested"),
    "plane/cmf_formal_rel_v4_gam.f": ("gap",
        "relativistic formal solver GAM; untested"),
    "plane/fg_j_cmf_v12.f": ("gap", "V12 formal J; untested"),
    "plane/fg_j_cmf_v13.f": ("gap",
        "V13 formal J + DO_FG13_ALO diagonal store; needs model state"),
    "plane/get_ibound.f": ("exec",
        [UT + "test_get_ibound.py"]),
    "plane/get_jh_at_prev_time_step.f": ("exec",
        [UT + "test_jh_prev.py"]),
    "subs/gff.f": ("exec", [UT + "test_gff.py"]),
    "subs/ins_line_obs_v5.f": ("exec",
        [UT + "test_ins_line.py"]),
    "subs/ins_line_v6.f": ("exec",
        [UT + "test_ins_line.py"]),
    "subs/scr_read_v2.f": ("exec",
        [UT + "test_scr_read.py"]),

    # ------------- shipped, outside cmfgen_dev.exe physics path --------
    "disp/maingen.f": ("exempt",
        "standalone dispgen visualisation executable"),
    "disp/subs/gam_abs.f": ("exempt",
        "standalone dispgen visualisation executable"),
    "disp/subs/log_mon_interp.f": ("exempt",
        "standalone dispgen visualisation executable"),
    "disp/subs/tst_clump_law.f": ("exempt",
        "standalone dispgen visualisation executable"),
    "misc/plt_ion_mf.f": ("exempt", "standalone plotting tool"),
    "misc/plt_net.f": ("exempt", "standalone plotting tool"),
    "misc/plt_rates.f": ("exempt", "standalone plotting tool"),
    "misc/plt_steq.f": ("exempt", "standalone plotting tool"),
    "spec_plt/Makefile": ("exempt",
        "build system for the spectroscopy plotting library"),
    "spec_plt/plt_dfr.f": ("exempt", "standalone plotting executable"),
    "spec_plt/plt_jh.f": ("exempt", "standalone plotting executable"),
    "spec_plt/plt_jh_cur.f": ("exempt", "standalone plotting executable"),
    "spec_plt/subs/compute_tau_radial.f": ("exempt",
        "spectroscopy/plotting library (libspec.a), not cmfgen_dev.exe"),
    "spec_plt/subs/do_sn_shell_shift.f": ("exempt",
        "spectroscopy/plotting library (libspec.a), not cmfgen_dev.exe"),

    # ------------- shipped tools that rewrite model data ----------------
    "misc/adjust_atm_r_grid.f": ("gap",
        "rewrites atmosphere grid; untested"),
    "misc/ave_dc.f": ("gap",
        "averages model data on a V coordinate; V-guard untested"),
    "misc/extrap_scr.f": ("gap", "scratch-file extrapolation tool; untested"),
    "misc/merge_scr.f": ("gap", "scratch-file merge tool; untested"),
    "misc/rev_hydro_turb.f": ("gap",
        "rewrites hydro turbulence data; untested"),
    "misc/rev_rdinr.f": ("gap",
        "rewrites inner radii; untested"),
    "misc/rev_rvsig.f": ("gap",
        "rewrites velocity/sigma data; untested"),
    "misc/rewrite_old_model_data.f": ("gap",
        "rewrites stored model data; untested"),

    # ------------- compiled but not reachable from shipped executables --
    "new_main/subs/get_pops_at_prev_time_step_v2.f": ("exempt",
        "superseded by V5; absent from every Makefile, not compiled"),
    "new_main/subs/get_pops_at_prev_time_step_v3.f": ("exempt",
        "superseded by V5; absent from every Makefile, not compiled"),
    "new_main/subs/get_pops_at_prev_time_step_v4.f": ("exempt",
        "superseded by V5; absent from every Makefile, not compiled"),
    "subs/sobew_grad.f": ("exempt",
        "legacy routine never called (V2/V3 supersede it); "
        "tau-floor guard mirrors the tested V2/V3 path"),

    # ------------- build system (verified by ./build.sh) -----------------
    "Makefile": ("exempt", "build system; make all must stay clean"),
    "Makefile_definitions": ("exempt",
        "build system; make all must stay clean"),
    "misc/Makefile": ("exempt", "build system; make all must stay clean"),
    "new_main/Makefile": ("exempt", "build system; make all must stay clean"),
    "new_main/mod_subs/Makefile": ("exempt",
        "build system; make all must stay clean"),
    "new_main/subs/Makefile": ("exempt",
        "build system; make all must stay clean"),
    "obs/Makefile": ("exempt", "build system; make all must stay clean"),
    "plane/Makefile": ("exempt", "build system; make all must stay clean"),
    "subs/Makefile": ("exempt", "build system; make all must stay clean"),
}


def altered_files():
    """Repo-relative paths of build inputs differing from the reference."""
    cmd = ["diff", "-rq"] + DIFF_EXCLUDES + [str(REF), "."]
    out = subprocess.run(cmd, capture_output=True, text=True,
                         cwd=TREE).stdout
    found = set()
    for ln in out.splitlines():
        if ln.startswith("Files ") and ln.endswith(" differ"):
            # 'Files <ref>/<a> and ./<b> differ'
            b = ln.split(" and ./", 1)[1][:-len(" differ")]
            found.add(b)
        elif ln.startswith("Only in ./"):
            rest = ln[len("Only in ./"):]
            d, _, name = rest.partition(": ")
            found.add(f"{d}/{name}" if d else name)
        elif ln.startswith("Only in ") and "cmfgen_reference" in ln:
            rest = ln.split("cur_cmf/", 1)[1]
            d, _, name = rest.partition(": ")
            found.add(f"{d}/{name}" if d else name)
    return {p for p in found
            if SOURCE_RE.search(p) and not p.startswith(OUT_OF_SCOPE)}


def main():
    rep = Report("alteration coverage manifest")
    altered = altered_files()

    unregistered = sorted(altered - set(COVERAGE))
    rep.check("every alteration is registered in COVERAGE",
              not unregistered,
              detail="; ".join(unregistered) if unregistered
              else f"{len(altered)} files registered")

    stale = sorted(set(COVERAGE) - altered)
    rep.check("no stale COVERAGE entries", not stale,
              detail="; ".join(stale) if stale else "manifest is current")

    missing_tests = []
    for path, (kind, payload) in COVERAGE.items():
        if kind in ("exec", "pin"):
            for t in payload:
                if not (TREE / t).exists():
                    missing_tests.append(f"{path} -> {t}")
    rep.check("registered covering tests exist", not missing_tests,
              detail="; ".join(missing_tests) if missing_tests
              else "all covering tests present")

    counts = {}
    for kind, _ in COVERAGE.values():
        counts[kind] = counts.get(kind, 0) + 1
    gaps = sorted(p for p, (k, _) in COVERAGE.items() if k == "gap")
    print(f"\ncoverage: {counts.get('exec', 0)} exec, "
          f"{counts.get('pin', 0)} pin, {counts.get('gap', 0)} gap, "
          f"{counts.get('exempt', 0)} exempt "
          f"({len(altered)} altered files)")
    if gaps:
        print("known gaps (shipped, untested):")
        for g in gaps:
            print(f"  {g} -- {COVERAGE[g][1]}")

    rep.summary()
    sys.exit(0)


if __name__ == "__main__":
    main()
