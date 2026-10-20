#!/usr/bin/env python3
import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict

try:
    from log_format import Tee, print_header, print_subheader
except ImportError:
    try:
        from scripts.log_format import Tee, print_header, print_subheader
    except ImportError:
        # Fallback: keep the script usable if log_format is not importable.
        def print_header(title: str = "") -> None:
            print(f"\n=== {title} ===")

        def print_subheader(title: str = "") -> None:
            print(f"\n--- {title} ---")

        class Tee:
            def __init__(self, stdout, f):
                self._stdout = stdout
                self._file = f

            def write(self, s):
                self._stdout.write(s)
                self._file.write(s)
                self._file.flush()

            def flush(self):
                self._stdout.flush()
                self._file.flush()

# IN_ITS iteration controls for CMFGEN. D2_EQ_D1: when T, depth-2 populations are replaced
# by depth-1 in non-LAMBDA iterations. Default F (independent depth 2); use --d2-eq-d1 to enable T.
# NOTE: REVISE_R_GRID is a VADAT keyword (REV_RGRID), not an IN_ITS keyword.
KEY_ORDER = [
    "NUM_ITS",
    "DO_LAM_IT",
    "DO_LAM_AUTO",
    "DO_GT_AUTO",
    "DO_T_AUTO",
    "D2_EQ_D1",
]

def parse_in_its(path: Path) -> Dict[str, str]:
    """Parse a simple IN_ITS, returning {KEY: VALUE} for known keys."""
    cfg: Dict[str, str] = {}
    if not path.is_file():
        return cfg
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or not line.startswith("[") or "]" not in line:
            continue
        key = line[1:line.index("]")].strip()
        rest = line[line.index("]") + 1 :].strip()
        if not rest:
            continue
        value = rest.split()[0]
        cfg[key] = value
    return cfg


def check_vadat_fix_t(directory: Path) -> bool:
    """
    Check if VADAT in the given directory has FIX_T=F (T allowed to vary).

    When FIX_T=F, T is a variable in the BA matrix and DO_T_AUTO must be T
    so that T is updated between iterations. When FIX_T=T (default), T is
    fixed to SN_HYDRO_DATA and DO_T_AUTO=F is correct.

    Note: FIX_T_AUTO=T is NOT sufficient to release T. With TAU_SCL_T=0
    (our default), VARFIXT has no runtime effect in solve_for_pops.f.
    FIX_T_AUTO=T only signals DO_LAM_AUTO=T (exit Lambda mode for BA
    convergence). T release requires FIX_T=F explicitly.

    Returns True if FIX_T=F (T should vary), False otherwise.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return False
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[FIX_T]" in line:
            # Format: "F [FIX_T]" or "T [FIX_T]"
            value = line.split("[")[0].strip()
            if value.upper() == "F":
                return True
    return False


def check_vadat_depth_independent(directory: Path) -> bool:
    """
    Check if VADAT in the given directory uses DIAG (depth-independent solver).

    Returns True if SOL_METH=DIAG (or N_TYPE=G_ONLY), False if TRIDIAG.
    Used to determine DO_LAM_IT: continuous LAMBDA (T) for DIAG, alternating
    (F) for TRIDIAG. DO_LAM_IT=T also avoids the BAMATPNT read bug that
    occurs with DO_LAM_IT=F when no BA file exists yet on the first iteration.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return True  # default: depth-independent
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[SOL_METH]" in line:
            value = line.split("[")[0].strip().upper()
            return value == "DIAG"
        if "[N_TYPE]" in line:
            value = line.split("[")[0].strip().upper()
            if value == "G_ONLY":
                return True
            elif value == "N_ON_J":
                return False
    return True  # default: depth-independent


def check_vadat_fix_t_auto(directory: Path) -> bool:
    """
    Check if VADAT in the given directory has FIX_T_AUTO=T.

    FIX_T_AUTO=T signals that CMFGEN should automatically manage when to
    fix/unfix T. This is used as a signal to enable DO_LAM_AUTO=T in IN_ITS,
    which allows the model to exit LAMBDA mode when MAXCH < 50%. In non-LAMBDA
    mode, T becomes a variable in the BA matrix (when FIX_T=F) and is solved
    self-consistently with the radiation field.

    This is needed because during LAMBDA iterations, CMFGEN hardcodes
    FIXED_T=.TRUE. (cmfgen_sub.f:611), so FIX_T=F alone has no effect.
    The model must exit LAMBDA mode (via DO_LAM_AUTO=T) for T to be solved.

    Safe only when LSTAR != 0 (DO_LAM_AUTO=T causes DTDR NaN when LSTAR=0).
    The caller is responsible for ensuring LSTAR is non-zero.

    Returns True if FIX_T_AUTO=T, False otherwise.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return False
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[FIX_T_AUTO]" in line:
            value = line.split("[")[0].strip()
            if value.upper() == "T":
                return True
    return False


def check_vadat_use_fixed_j(directory: Path) -> bool:
    """
    Check if VADAT in the given directory has USE_FIXED_J=T.

    USE_FIXED_J=T holds J at B(T) (LTE) during LAMBDA iterations, breaking
    the J-population feedback loop that causes trace species explosions.
    J is released when MAXCH < 50% via DO_LAM_AUTO=T (cmfgen_sub.f:4727).

    When USE_FIXED_J=T, prepare_IN_INTS.py sets DO_LAM_AUTO=T in IN_ITS
    to enable the automatic release of J. T remains fixed (DO_T_AUTO=F)
    regardless of USE_FIXED_J — T solving is controlled by FIX_T/DO_T_AUTO.

    Returns True if USE_FIXED_J=T, False otherwise.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return False
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[USE_FIXED_J]" in line:
            value = line.split("[")[0].strip()
            if value.upper() == "T":
                return True
    return False


def check_vadat_do_lam_auto(directory: Path) -> bool | None:
    """
    Check if VADAT in the given directory explicitly sets DO_LAM_AUTO.

    DO_LAM_AUTO is an IN_ITS keyword, but prepare_VADAT.py writes it to the
    VADAT file as a transport mechanism (update_in_its_for_profile reads it
    back). This function lets prepare_IN_INTS.py respect an explicit
    DO_LAM_AUTO setting from the VADAT file (e.g. set via vadat_overrides
    in the YAML config), rather than always deriving it from FIX_T_AUTO
    and USE_FIXED_J.

    Returns True if DO_LAM_AUTO=T, False if DO_LAM_AUTO=F, None if not found.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return None
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[DO_LAM_AUTO]" in line:
            value = line.split("[")[0].strip().upper()
            if value == "T":
                return True
            elif value == "F":
                return False
    return None


def check_vadat_do_lam_it(directory: Path) -> bool | None:
    """
    Check if VADAT in the given directory explicitly sets DO_LAM_IT.

    DO_LAM_IT is an IN_ITS keyword, but prepare_VADAT.py writes it to the
    VADAT file as a transport mechanism (update_in_its_for_profile reads it
    back). This function lets prepare_IN_INTS.py respect an explicit
    DO_LAM_IT setting from the VADAT file (e.g. set via vadat_overrides
    in the YAML config), rather than always defaulting to T.

    Returns True if DO_LAM_IT=T, False if DO_LAM_IT=F, None if not found.
    """
    vadat_path = directory / "VADAT"
    if not vadat_path.is_file():
        return None
    for line in vadat_path.read_text().splitlines():
        line = line.strip()
        if "[DO_LAM_IT]" in line:
            value = line.split("[")[0].strip().upper()
            if value == "T":
                return True
            elif value == "F":
                return False
    return None


def format_in_its(cfg: Dict[str, str]) -> str:
    """Format an IN_ITS dictionary back to text."""
    lines = []
    for key in KEY_ORDER:
        if key in cfg:
            lines.append(f"{cfg[key]} [{key}] ")
    lines.append("")
    return "\n".join(lines)


def write_init_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    d2_eq_d1: bool = False,
    gt_auto: str | None = None,
):
    """
    Write IN_ITS for the *initial* Lambda-iteration run,
    e.g. first ./batch.sh to get a stable model.
    """
    cfg: Dict[str, str] = {}

    # Check if VADAT has FIX_T=F (T allowed to vary)
    t_varies = check_vadat_fix_t(directory)
    # Check if VADAT has FIX_T_AUTO=T (signals DO_LAM_AUTO=T for T solving)
    do_lam_auto = check_vadat_fix_t_auto(directory)
    # Check if VADAT has USE_FIXED_J=T (signals DO_LAM_AUTO=T for J release)
    use_fixed_j = check_vadat_use_fixed_j(directory)

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 3)
    # Lambda iteration with USE_FIXED_J=T (holds J at LTE to break oscillation).
    # DO_LAM_IT: respect explicit VADAT setting; default T (Lambda cold-start).
    # DO_LAM_AUTO: F by default (pure Lambda); T when FIX_T_AUTO=T in VADAT
    # (allows exit from Lambda mode when MAXCH < 50% so T can be solved).
    vadat_lam_it = check_vadat_do_lam_it(directory)
    cfg["DO_LAM_IT"] = "T" if vadat_lam_it is None else ("T" if vadat_lam_it else "F")
    cfg["DO_LAM_AUTO"] = "T" if (do_lam_auto or use_fixed_j) else "F"
    # Temperature control for static non-monotonic snapshots:
    # DO_T_AUTO and DO_GT_AUTO are T only when FIX_T=F (T is a variable in
    # the BA matrix). FIX_T_AUTO=T alone does NOT release T — it only signals
    # DO_LAM_AUTO=T (exit Lambda mode for BA convergence). With TAU_SCL_T=0
    # (our default), VARFIXT has no runtime effect in solve_for_pops.f, so
    # keeping T fixed during BA iterations is the physically consistent choice
    # for a static hydro snapshot.
    cfg["DO_T_AUTO"] = "T" if t_varies else "F"
    cfg["DO_GT_AUTO"] = (gt_auto or "T") if t_varies else "F"  # grey-T iteration only when T varies; gt_auto cannot unpin FIX_T
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    t_status = "T varies (FIX_T=F)" if t_varies else "T fixed (FIX_T=T, DO_T_AUTO=F)"
    lam_status = "DO_LAM_AUTO=T (FIX_T_AUTO=T)" if do_lam_auto else ("DO_LAM_AUTO=T (USE_FIXED_J=T)" if use_fixed_j else "DO_LAM_AUTO=F")
    print(f"Wrote initialization IN_ITS (Lambda only, {t_status}; {lam_status}; {d2m}) in {directory}")


def write_stable_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    d2_eq_d1: bool = False,
    gt_auto: str | None = None,
):
    """
    Write IN_ITS for the *stable* Lambda-only run between init and final.
    This phase runs pure Lambda iterations to create a more stable system
    before transitioning to full Non-Lambda iterations.
    """
    cfg: Dict[str, str] = {}

    # Check if VADAT has FIX_T=F (T allowed to vary)
    t_varies = check_vadat_fix_t(directory)
    # Check if VADAT has FIX_T_AUTO=T (signals DO_LAM_AUTO=T for T solving)
    do_lam_auto = check_vadat_fix_t_auto(directory)
    # Check if VADAT has USE_FIXED_J=T (signals DO_LAM_AUTO=T for J release)
    use_fixed_j = check_vadat_use_fixed_j(directory)

    # Default to more iterations than init for Lambda settling
    cfg["NUM_ITS"] = str(num_its if num_its is not None else 20)
    # Pure Lambda-iteration phase - no auto-switch by default (stable profile).
    # DO_LAM_IT: respect explicit VADAT setting; default T (pure Lambda).
    # DO_LAM_AUTO=T when FIX_T_AUTO=T in VADAT (allows exit from Lambda for T solving).
    vadat_lam_it = check_vadat_do_lam_it(directory)
    cfg["DO_LAM_IT"] = "T" if vadat_lam_it is None else ("T" if vadat_lam_it else "F")
    cfg["DO_LAM_AUTO"] = "T" if (do_lam_auto or use_fixed_j) else "F"
    # Temperature control for static non-monotonic snapshots:
    # DO_T_AUTO and DO_GT_AUTO are T only when FIX_T=F (T is a variable in
    # the BA matrix). FIX_T_AUTO=T alone does NOT release T — it only signals
    # DO_LAM_AUTO=T (exit Lambda mode for BA convergence). With TAU_SCL_T=0
    # (our default), VARFIXT has no runtime effect in solve_for_pops.f, so
    # keeping T fixed during BA iterations is the physically consistent choice
    # for a static hydro snapshot.
    cfg["DO_T_AUTO"] = "T" if t_varies else "F"
    cfg["DO_GT_AUTO"] = (gt_auto or "T") if t_varies else "F"  # grey-T iteration only when T varies; gt_auto cannot unpin FIX_T
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    t_status = "T varies (FIX_T=F)" if t_varies else "T fixed (FIX_T=T, DO_T_AUTO=F)"
    lam_status = "DO_LAM_AUTO=T (FIX_T_AUTO=T)" if do_lam_auto else ("DO_LAM_AUTO=T (USE_FIXED_J=T)" if use_fixed_j else "DO_LAM_AUTO=F")
    print(f"Wrote stable IN_ITS (pure Lambda only, {t_status}; {lam_status}; {d2m}) in {directory}")


def write_final_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    d2_eq_d1: bool = False,
    gt_auto: str | None = None,
):
    """
    Write IN_ITS for the *final* run.
    Lambda bootstrap then auto-switch to full iterations (MAXCH < 50%) for quadratic convergence.
    """
    cfg: Dict[str, str] = {}

    # Check if VADAT has FIX_T=F (T allowed to vary)
    t_varies = check_vadat_fix_t(directory)
    # Check if VADAT has FIX_T_AUTO=T (signals DO_LAM_AUTO=T for T solving)
    do_lam_auto = check_vadat_fix_t_auto(directory)
    # Check if VADAT has USE_FIXED_J=T (signals DO_LAM_AUTO=T for J release)
    use_fixed_j = check_vadat_use_fixed_j(directory)

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 10000)
    # DO_LAM_IT: respect explicit VADAT setting (from vadat_overrides in YAML);
    #   otherwise default to T (start in LAMBDA mode, stable cold-start from LTE).
    # DO_LAM_AUTO: respect explicit VADAT setting (from vadat_overrides in YAML);
    #   otherwise T when FIX_T_AUTO=T or USE_FIXED_J=T signals it, F for pure LAMBDA.
    vadat_lam_auto = check_vadat_do_lam_auto(directory)
    vadat_lam_it = check_vadat_do_lam_it(directory)
    if vadat_lam_it is not None:
        cfg["DO_LAM_IT"] = "T" if vadat_lam_it else "F"
    else:
        cfg["DO_LAM_IT"] = "T"
    if vadat_lam_auto is not None:
        cfg["DO_LAM_AUTO"] = "T" if vadat_lam_auto else "F"
    else:
        cfg["DO_LAM_AUTO"] = "T" if (do_lam_auto or use_fixed_j) else "F"
    # Temperature control for static non-monotonic snapshots:
    # DO_T_AUTO and DO_GT_AUTO are T only when FIX_T=F (T is a variable in
    # the BA matrix). FIX_T_AUTO=T alone does NOT release T — it only signals
    # DO_LAM_AUTO=T (exit Lambda mode for BA convergence). With TAU_SCL_T=0
    # (our default), VARFIXT has no runtime effect in solve_for_pops.f, so
    # keeping T fixed during BA iterations is the physically consistent choice
    # for a static hydro snapshot.
    cfg["DO_T_AUTO"] = "T" if t_varies else "F"
    cfg["DO_GT_AUTO"] = (gt_auto or "T") if t_varies else "F"  # grey-T iteration only when T varies; gt_auto cannot unpin FIX_T
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    t_status = "T varies (FIX_T=F)" if t_varies else "T fixed (FIX_T=T, DO_T_AUTO=F)"
    lam_status = "DO_LAM_AUTO=T (FIX_T_AUTO=T)" if do_lam_auto else ("DO_LAM_AUTO=T (USE_FIXED_J=T)" if use_fixed_j else "DO_LAM_AUTO=F (pure LAMBDA)")
    if vadat_lam_auto is not None:
        lam_status = f"DO_LAM_AUTO={'T' if vadat_lam_auto else 'F'} (from VADAT vadat_overrides)"
    print(f"Wrote final IN_ITS (pure LAMBDA, {t_status}; {lam_status}; {d2m}) in {directory}")


def write_single_final_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    d2_eq_d1: bool = False,
    gt_auto: str | None = None,
):
    """
    Write IN_ITS for the *single-final* LTE cold-start run.
    Lambda bootstrap then auto-switch to full iterations (MAXCH < 50%) for quadratic convergence.
    """
    cfg: Dict[str, str] = {}

    # Check if VADAT has FIX_T=F (T allowed to vary)
    t_varies = check_vadat_fix_t(directory)
    # Check if VADAT has FIX_T_AUTO=T (signals DO_LAM_AUTO=T for T solving)
    do_lam_auto = check_vadat_fix_t_auto(directory)
    # Check if VADAT has USE_FIXED_J=T (signals DO_LAM_AUTO=T for J release)
    use_fixed_j = check_vadat_use_fixed_j(directory)

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 10000)
    # DO_LAM_IT: respect explicit VADAT setting (from vadat_overrides in YAML);
    #   otherwise default to T (start in LAMBDA mode, stable cold-start from LTE).
    # DO_LAM_AUTO: respect explicit VADAT setting (from vadat_overrides in YAML);
    #   otherwise T when FIX_T_AUTO=T or USE_FIXED_J=T signals it, F for pure LAMBDA.
    vadat_lam_auto = check_vadat_do_lam_auto(directory)
    vadat_lam_it = check_vadat_do_lam_it(directory)
    if vadat_lam_it is not None:
        cfg["DO_LAM_IT"] = "T" if vadat_lam_it else "F"
    else:
        cfg["DO_LAM_IT"] = "T"
    if vadat_lam_auto is not None:
        cfg["DO_LAM_AUTO"] = "T" if vadat_lam_auto else "F"
    else:
        cfg["DO_LAM_AUTO"] = "T" if (do_lam_auto or use_fixed_j) else "F"
    # Temperature control for static non-monotonic snapshots:
    # DO_T_AUTO=T only when FIX_T=F (T is a variable in the BA matrix).
    # FIX_T_AUTO=T alone does NOT release T — it only signals DO_LAM_AUTO=T
    # (exit Lambda mode for BA convergence). With TAU_SCL_T=0 (our default),
    # VARFIXT has no runtime effect in solve_for_pops.f, so keeping T fixed
    # during BA iterations is the physically consistent choice for a static
    # hydro snapshot.
    cfg["DO_T_AUTO"] = "T" if t_varies else "F"
    # DO_GT_AUTO: only relevant when T is being released (FIX_T=F).
    # When T is fixed to the hydro snapshot, grey-T iteration is unnecessary
    # and can destabilise the non-monotonic formal solver.
    cfg["DO_GT_AUTO"] = (gt_auto or "T") if t_varies else "F"  # grey-T iteration only when T varies; gt_auto cannot unpin FIX_T
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    t_status = "T varies (FIX_T=F)" if t_varies else "T fixed (FIX_T=T, DO_T_AUTO=F)"
    lam_status = "DO_LAM_AUTO=T (FIX_T_AUTO=T)" if do_lam_auto else ("DO_LAM_AUTO=T (USE_FIXED_J=T)" if use_fixed_j else "DO_LAM_AUTO=F (pure LAMBDA)")
    if vadat_lam_auto is not None:
        lam_status = f"DO_LAM_AUTO={'T' if vadat_lam_auto else 'F'} (from VADAT vadat_overrides)"
    print(f"Wrote single-final IN_ITS (pure LAMBDA, {t_status}; {lam_status}; {d2m}) in {directory}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create CMFGEN IN_ITS for either the Lambda-init run "
            "or the final full-iteration run."
        )
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Target model directory (default: current directory).",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--final",
        action="store_true",
        help="Write IN_ITS for the *final* run (warm restart; Lambda bootstrap, auto-switch to full iterations at MAXCH < 50%).",
    )
    group.add_argument(
        "--single-final",
        action="store_true",
        help="Write IN_ITS for the *single-final* LTE cold-start run (no prior run; DO_LAM_AUTO=F)."
    )
    parser.add_argument(
        "--num-its",
        type=int,
        default=None,
        help="Override NUM_ITS (default: 10000 for final and single-final).",
    )
    d2_group = parser.add_mutually_exclusive_group()
    d2_group.add_argument(
        "--d2-eq-d1",
        action="store_true",
        help="Set D2_EQ_D1=T (tie depth-2 populations to depth-1). Default is F (independent depth-2 populations).",
    )
    d2_group.add_argument(
        "--no-d2-eq-d1",
        action="store_true",
        help="Deprecated no-op; D2_EQ_D1=F is now the default. Kept for backward compatibility.",
    )
    parser.add_argument(
        "--gt-auto",
        choices=["T", "F"],
        default=None,
        help="Override DO_GT_AUTO (grey-T update). Default: T when VADAT has FIX_T=F, else F.",
    )
    parser.add_argument(
        "--no-outgen",
        action="store_true",
        help="Do not append output to OUTGEN (default: append to OUTGEN).",
    )

    args = parser.parse_args()
    directory = Path(args.directory)
    directory.mkdir(parents=True, exist_ok=True)

    profile = "single_final" if args.single_final else "final"

    # Tee stdout to OUTGEN so the IN_ITS settings land in the same log as
    # CMFGEN (same mechanism as prepare_VADAT.py).
    outgen_file = None
    old_stdout = sys.stdout
    if not args.no_outgen:
        try:
            outgen_file = open(directory / "OUTGEN", "a")
            sys.stdout = Tee(old_stdout, outgen_file)
        except OSError:
            outgen_file = None

    try:
        print_header(f"prepare_IN_INTS.py [{profile}]")
        print(f"  started:   {datetime.now().isoformat()}")
        print(f"  directory: {directory.resolve()}")

        # Use existing IN_ITS as a template if present
        template_cfg = parse_in_its(directory / "IN_ITS")

        # Temperature control: DO_T_AUTO is set based on FIX_T in VADAT.
        # When FIX_T=F (set via vadat_overrides), DO_T_AUTO=T allows T to vary.
        # When FIX_T=T (default for static snapshots), DO_T_AUTO=F keeps T fixed
        # to SN_HYDRO_DATA. FIX_T_AUTO=T alone does NOT release T — it only
        # signals DO_LAM_AUTO=T (exit Lambda mode for BA convergence).
        d2 = bool(args.d2_eq_d1)
        if args.single_final:
            write_single_final_in_its(
                directory, template_cfg, args.num_its, d2_eq_d1=d2, gt_auto=args.gt_auto
            )
        else:  # args.final
            write_final_in_its(
                directory, template_cfg, args.num_its, d2_eq_d1=d2, gt_auto=args.gt_auto
            )

        # Echo the written file so OUTGEN records the actual IN_ITS keys.
        print_subheader("IN_ITS written")
        in_its_path = directory / "IN_ITS"
        if in_its_path.is_file():
            for line in in_its_path.read_text().splitlines():
                if line.strip():
                    print(f"    {line.strip()}")
        else:
            print("    (IN_ITS not found)")
    finally:
        if outgen_file is not None:
            sys.stdout = old_stdout
            outgen_file.close()


if __name__ == "__main__":
    main()
