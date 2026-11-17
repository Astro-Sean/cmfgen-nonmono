"""P8.3 KNOP spectrum-gate tests (obs/knop_spec_gate.f).

Builds a driver against the production knop_spec_gate.o and runs it
from staged model/obs fixture dirs.  Covers the acceptance matrix:
SOB (gate not consulted), missing/malformed GS_STATE, immature,
mature-stalled, converged, and manifest-flagged incomplete runs.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import TREE, compile_and_link_prod

GS_CONVERGED = """\
n_frozen=0
iter=700 frozen= F
mode=FULL
pops_held= F
mom_refresh_conv= F
lum_ratio=  1.0001E+00
n_mature=3
"""
GS_IMMATURE = """\
n_frozen=12
iter=157 frozen= T
mode=LAMBDA_GS
pops_held= T
mom_refresh_conv= F
lum_ratio=  3.7561E-01
n_mature=0
"""
GS_STALL_NOMATURE = GS_IMMATURE.replace("n_frozen=12", "n_frozen=172")
GS_STALL_MATURE = GS_STALL_NOMATURE.replace("n_mature=0", "n_mature=4")
GS_MALFORMED = "garbage without keys\n"
MANIFEST_KILL = '{"segments": [], "x": 1}\n{"exit_code": 152, "termination": "signal_kill:SIGXCPU"}\n'
MANIFEST_OK = '{"exit_code": 0, "termination": "normal_exit"}\n'

EXPECT = {
    "converged": (1, "KNOP_CONVERGED", "T"),
    "immature": (3, "KNOP_IMMATURE", "F"),
    "stall_nomature": (3, "KNOP_IMMATURE", "F"),
    "stall_mature": (2, "KNOP_MATURE_STALLED", "F"),
    "missing": (4, "KNOP_EXECUTION_INCOMPLETE", "F"),
    "malformed": (4, "KNOP_EXECUTION_INCOMPLETE", "F"),
    "killed": (4, "KNOP_EXECUTION_INCOMPLETE", "F"),
    "marker": (3, "KNOP_IMMATURE", "F"),
}


def run_case(drv, model_dir, gs_text=None, marker=False, manifest=None):
    obs = model_dir / "obs"
    obs.mkdir(parents=True, exist_ok=True)
    if gs_text is not None:
        (model_dir / "GS_STATE").write_text(gs_text)
    if marker:
        (model_dir / ".CMFGEN_KNOP_STALLED").touch()
    if manifest is not None:
        (model_dir / "manifest.json").write_text(manifest)
    out = subprocess.run(
        [str(drv)], cwd=obs, capture_output=True, text=True, check=True
    )
    row = [l for l in out.stdout.splitlines() if l.startswith("RESULT")][0]
    _, state, cls, allow = row.split()
    status = (obs / "KNOP_SPEC_STATUS").read_text()
    return int(state), cls, allow, status


def main():
    drv = compile_and_link_prod(
        "knop_spec_gate",
        [
            Path(__file__).resolve().parents[1] / "f90" / "drv_knop_spec_gate.f",
            TREE / "obs" / "knop_spec_gate.f",
            TREE / "obs" / "cmf_flux_cntrl_var_mod.f",
        ],
    )
    print(f"built {drv}")
    npass = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cases = {
            "converged": dict(gs_text=GS_CONVERGED),
            "immature": dict(gs_text=GS_IMMATURE),
            "stall_nomature": dict(gs_text=GS_STALL_NOMATURE),
            "stall_mature": dict(gs_text=GS_STALL_MATURE),
            "missing": dict(),
            "malformed": dict(gs_text=GS_MALFORMED),
            "killed": dict(gs_text=GS_CONVERGED, manifest=MANIFEST_KILL),
            "marker": dict(gs_text=GS_IMMATURE, marker=True),
        }
        for name, kw in cases.items():
            st, cls, allow, status = run_case(drv, td / name, **kw)
            est, ecls, eallow = EXPECT[name]
            assert (st, cls, allow) == (est, ecls, eallow), (
                f"{name}: got ({st},{cls},{allow}) "
                f"want ({est},{ecls},{eallow})"
            )
            assert "class=" + ecls in status, f"{name}: bad status file"
            print(f"[PASS] {name}: {cls} state={st} allow={allow}")
            npass += 1
        # converged model plus a clean manifest stays allowed
        st, cls, allow, _ = run_case(
            drv, td / "ok", gs_text=GS_CONVERGED, manifest=MANIFEST_OK
        )
        assert (st, cls, allow) == (1, "KNOP_CONVERGED", "T")
        print(f"[PASS] manifest_ok: {cls} state={st} allow={allow}")
        npass += 1
    print(f"{npass}/{npass} knop_spec_gate cases pass")


if __name__ == "__main__":
    main()
