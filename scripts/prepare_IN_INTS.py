#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict

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
    hold_t: bool = False,
    d2_eq_d1: bool = True,
):
    """
    Write IN_ITS for the *initial* Lambda-iteration run,
    e.g. first ./batch.sh to get a stable model.
    """
    cfg: Dict[str, str] = {}

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 3)
    # Pure Lambda-iteration phase to settle LTE-like start before any full linearisation.
    cfg["DO_LAM_IT"] = "T"
    cfg["DO_LAM_AUTO"] = "F"
    # Temperature control: disable auto-updates since FIX_T=T fixes T to SN_HYDRO_DATA
    cfg["DO_T_AUTO"] = "F"
    cfg["DO_GT_AUTO"] = "T"
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    mode = "T fixed" if hold_t else "T varies"
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    print(f"Wrote initialization IN_ITS (Lambda only, {mode}; {d2m}) in {directory}")


def write_stable_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    hold_t: bool = False,
    d2_eq_d1: bool = True,
):
    """
    Write IN_ITS for the *stable* Lambda-only run between init and final.
    This phase runs pure Lambda iterations to create a more stable system
    before transitioning to full Non-Lambda iterations.
    """
    cfg: Dict[str, str] = {}

    # Default to more iterations than init for Lambda settling
    cfg["NUM_ITS"] = str(num_its if num_its is not None else 20)
    # Pure Lambda-iteration phase - stay in Lambda, no auto-switch to Non-Lambda
    cfg["DO_LAM_IT"] = "T"
    cfg["DO_LAM_AUTO"] = "F"
    # Temperature control: disable auto-updates since FIX_T=T fixes T to SN_HYDRO_DATA
    cfg["DO_T_AUTO"] = "F"
    cfg["DO_GT_AUTO"] = "T"
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    mode = "T fixed" if hold_t else "T varies"
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    print(f"Wrote stable IN_ITS (pure Lambda only, {mode}; {d2m}) in {directory}")


def write_final_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    hold_t: bool = False,
    d2_eq_d1: bool = True,
):
    """
    Write IN_ITS for the *final* run.
    """
    cfg: Dict[str, str] = {}

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 500)
    # Full-iteration phase: enable auto-toggle from LAMBDA to Non-Lambda iterations
    cfg["DO_LAM_IT"] = "T"
    cfg["DO_LAM_AUTO"] = "T"
    # Temperature control: disable auto-updates since FIX_T=T fixes T to SN_HYDRO_DATA
    cfg["DO_T_AUTO"] = "F"
    cfg["DO_GT_AUTO"] = "T"
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    mode = "T fixed" if hold_t else "T varies"
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    print(f"Wrote final IN_ITS ({mode}; {d2m}) in {directory}")


def write_single_final_in_its(
    directory: Path,
    template: Dict[str, str],
    num_its: int | None,
    hold_t: bool = False,
    d2_eq_d1: bool = True,
):
    """
    Write IN_ITS for the *single-final* LTE cold-start run.
    Starts from LTE populations and auto-switches Lambda->Non-Lambda.
    """
    cfg: Dict[str, str] = {}

    cfg["NUM_ITS"] = str(num_its if num_its is not None else 100000)
    cfg["DO_LAM_IT"] = "T"
    cfg["DO_LAM_AUTO"] = "T"
    # Temperature control: disable auto-updates since FIX_T=T fixes T to SN_HYDRO_DATA
    cfg["DO_T_AUTO"] = "F"
    cfg["DO_GT_AUTO"] = "T"
    cfg["D2_EQ_D1"] = "T" if d2_eq_d1 else "F"

    (directory / "IN_ITS").write_text(format_in_its(cfg))
    mode = "T fixed" if hold_t else "T varies"
    d2m = "D2_EQ_D1=T" if d2_eq_d1 else "D2_EQ_D1=F"
    print(f"Wrote single-final IN_ITS (LTE cold start, Lambda->Non-Lambda auto, {mode}; {d2m}) in {directory}")


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
        "--init",
        action="store_true",
        help="Write IN_ITS for the *initial* Lambda-iteration run.",
    )
    group.add_argument(
        "--stable",
        action="store_true",
        help="Write IN_ITS for the *stable* pure Lambda run (between init and final).",
    )
    group.add_argument(
        "--final",
        action="store_true",
        help="Write IN_ITS for the *final* non-Lambda run.",
    )
    group.add_argument(
        "--single-final",
        action="store_true",
        help="Write IN_ITS for the *single-final* LTE cold-start run (no prior INIT; DO_LAM_AUTO=T).",
    )
    parser.add_argument(
        "--num-its",
        type=int,
        default=None,
        help="Override NUM_ITS (default: 3 for init, 20 for stable, 500 for final).",
    )
    parser.add_argument(
        "--d2-eq-d1",
        action="store_true",
        help="Set D2_EQ_D1=T (tie depth-2 populations to depth-1). Default is F.",
    )

    args = parser.parse_args()
    directory = Path(args.directory)
    directory.mkdir(parents=True, exist_ok=True)

    # Use existing IN_ITS as a template if present
    template_cfg = parse_in_its(directory / "IN_ITS")

    # Temperature is always allowed to vary (hold_t=False)
    d2 = bool(args.d2_eq_d1)
    if args.init:
        write_init_in_its(directory, template_cfg, args.num_its, hold_t=False, d2_eq_d1=d2)
    elif args.stable:
        write_stable_in_its(directory, template_cfg, args.num_its, hold_t=False, d2_eq_d1=d2)
    elif args.single_final:
        write_single_final_in_its(directory, template_cfg, args.num_its, hold_t=False, d2_eq_d1=d2)
    else:  # args.final
        write_final_in_its(directory, template_cfg, args.num_its, hold_t=False, d2_eq_d1=d2)


if __name__ == "__main__":
    main()
