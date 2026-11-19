"""P8.3 KNOP spectrum-gate tests (obs/knop_spec_gate.f).

Builds a driver against the production knop_spec_gate.o and runs it
from staged model/obs fixture dirs.  Covers the acceptance matrix:
provenance-bound convergence, open OBS segments (null exit_code),
stale GS_STATE after warm restart, bound vs stale stall markers,
multi-segment FINAL manifests, missing/malformed evidence, and the
non-KNOP path.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import TREE, compile_and_link_prod

JOB1 = "99111"
JOB2 = "99222"

GS_CONVERGED = f"""\
n_frozen=0
run_id={JOB1}
iter=700 frozen= F
mode=FULL
pops_held= F
mom_refresh_conv= F
lum_ratio=  1.0001E+00
n_mature=3
"""
GS_IMMATURE = f"""\
n_frozen=12
run_id={JOB1}
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
GS_NOID = GS_CONVERGED.replace(f"run_id={JOB1}\n", "")


def manifest(*segs):
    lines = ["{", '  "schema": 1,', '  "runs": [']
    body = []
    for s in segs:
        seg = ["    {"]
        seg.append(f'      "label": "{s[0]}",')
        seg.append(f'      "slurm_job_id": "{s[1]}",')
        seg.append(f'      "exit_code": {s[2]},')
        seg.append(f'      "termination": {s[3]},')
        seg.append('      "ended_at": null')
        seg.append("    }")
        body.append("\n".join(seg))
    lines.append(",\n".join(body))
    lines += ["  ]", "}"]
    return "\n".join(lines) + "\n"


NORMAL = ('"normal_exit"',)
FINAL_OK = ("FINAL", JOB1, 0, '"normal_exit"')
FINAL_OK2 = ("FINAL", JOB2, 0, '"normal_exit"')
FINAL_KILL = ("FINAL", JOB1, 152, '"signal_kill:SIGXCPU"')
FINAL_STALL = ("FINAL", JOB1, 7, '"knop_stall"')
FINAL_OPEN = ("FINAL", JOB2, "null", "null")
OBS_OPEN = ("OBS", JOB2, "null", "null")

EXPECT = {
    # A. healthy GS + normal FINAL, bound run_id -> converged
    "A_bound": ((1, "KNOP_CONVERGED", "T"),
                dict(gs_text=GS_CONVERGED,
                     manifest=manifest(FINAL_OK))),
    # B. healthy GS + abnormal FINAL -> incomplete
    "B_abnormal": ((4, "KNOP_EXECUTION_INCOMPLETE", "F"),
                   dict(gs_text=GS_CONVERGED,
                        manifest=manifest(FINAL_KILL))),
    # C. healthy stale GS (old run_id) + new open FINAL -> incomplete
    "C_stale_open": ((4, "KNOP_EXECUTION_INCOMPLETE", "F"),
                     dict(gs_text=GS_CONVERGED,
                          manifest=manifest(FINAL_OPEN))),
    # D. multi-segment: FINAL ok + trailing open OBS -> converged
    "D_obs_open": ((1, "KNOP_CONVERGED", "T"),
                   dict(gs_text=GS_CONVERGED,
                        manifest=manifest(FINAL_OK, OBS_OPEN))),
    # E. multi-segment: failed FINAL then open OBS -> incomplete
    "E_final_fail": ((4, "KNOP_EXECUTION_INCOMPLETE", "F"),
                     dict(gs_text=GS_CONVERGED,
                          manifest=manifest(FINAL_KILL, OBS_OPEN))),
    # F. missing manifest -> unverified, allowed but marked
    "F_no_manifest": ((5, "KNOP_PROVENANCE_UNVERIFIED", "T"),
                      dict(gs_text=GS_CONVERGED)),
    # G. malformed manifest (no FINAL segment) -> unverified
    "G_no_final": ((5, "KNOP_PROVENANCE_UNVERIFIED", "T"),
                   dict(gs_text=GS_CONVERGED,
                        manifest=manifest(OBS_OPEN))),
    # H. missing GS_STATE -> incomplete
    "H_no_gs": ((4, "KNOP_EXECUTION_INCOMPLETE", "F"), dict()),
    # I. stale stall marker (old job id) + verified new run -> converged
    "I_stale_marker": ((1, "KNOP_CONVERGED", "T"),
                       dict(gs_text=GS_CONVERGED,
                            marker_text=JOB2 + "\n",
                            manifest=manifest(FINAL_OK))),
    # J. stalled FINAL seg superseded by a normal FINAL seg -> converged
    # (GS bound to the surviving segment's job)
    "J_stall_then_ok": ((1, "KNOP_CONVERGED", "T"),
                        dict(gs_text=GS_CONVERGED.replace(JOB1, JOB2),
                             manifest=manifest(FINAL_STALL,
                                               FINAL_OK2))),
    # K. unverified renders anyway (allow=T) - diagnostic opt-out tier
    "K_noid": ((5, "KNOP_PROVENANCE_UNVERIFIED", "T"),
               dict(gs_text=GS_NOID, manifest=manifest(FINAL_OK))),
    # run_id mismatch: healthy GS but bound to a different job -> 5
    "M_mismatch": ((5, "KNOP_PROVENANCE_UNVERIFIED", "T"),
                   dict(gs_text=GS_CONVERGED.replace(JOB1, JOB2),
                        manifest=manifest(FINAL_OK))),
    # bound marker for THIS run stalls it
    "marker_bound": ((3, "KNOP_IMMATURE", "F"),
                     dict(gs_text=GS_IMMATURE,
                          marker_text=JOB1 + "\n",
                          manifest=manifest(FINAL_OK))),
    # legacy empty marker + no run_id anywhere -> still applies
    # (converged-shaped state + stall evidence -> mature stall, blocked)
    "marker_legacy": ((2, "KNOP_MATURE_STALLED", "F"),
                      dict(gs_text=GS_NOID, marker_text="")),
    # stall field evidence
    "stall_mature": ((2, "KNOP_MATURE_STALLED", "F"),
                     dict(gs_text=GS_STALL_MATURE,
                          manifest=manifest(FINAL_OK))),
    "stall_nomature": ((3, "KNOP_IMMATURE", "F"),
                       dict(gs_text=GS_STALL_NOMATURE,
                            manifest=manifest(FINAL_OK))),
    "immature": ((3, "KNOP_IMMATURE", "F"),
                 dict(gs_text=GS_IMMATURE,
                      manifest=manifest(FINAL_OK))),
    "malformed": ((4, "KNOP_EXECUTION_INCOMPLETE", "F"),
                  dict(gs_text=GS_MALFORMED)),
}


def run_case(drv, model_dir, sob=False, gs_text=None,
             marker_text=None, manifest=None):
    obs = model_dir / "obs"
    obs.mkdir(parents=True, exist_ok=True)
    if gs_text is not None:
        (model_dir / "GS_STATE").write_text(gs_text)
    if marker_text is not None:
        (model_dir / ".CMFGEN_KNOP_STALLED").write_text(marker_text)
    if manifest is not None:
        (model_dir / "manifest.json").write_text(manifest)
    cmd = [str(drv)] + (["sob"] if sob else [])
    out = subprocess.run(
        cmd, cwd=obs, capture_output=True, text=True, check=True
    )
    row = [l for l in out.stdout.splitlines() if l.startswith("RESULT")][0]
    _, state, cls, allow = row.split()
    sf = obs / "KNOP_SPEC_STATUS"
    status = sf.read_text() if sf.exists() else ""
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
        for name, (want, kw) in EXPECT.items():
            st, cls, allow, status = run_case(drv, td / name, **kw)
            assert (st, cls, allow) == want, (
                f"{name}: got ({st},{cls},{allow}) want {want}"
            )
            assert "class=" + want[1] in status, f"{name}: bad status file"
            print(f"[PASS] {name}: {cls} state={st} allow={allow}")
            npass += 1
        # L. non-KNOP path: gate does not engage at all
        st, cls, allow, status = run_case(
            drv, td / "sob", sob=True, gs_text=GS_STALL_MATURE,
            marker_text=JOB1 + "\n", manifest=manifest(FINAL_STALL),
        )
        assert (st, cls, allow) == (0, "NOT_APPLICABLE", "T"), (
            f"sob: got ({st},{cls},{allow})"
        )
        print(f"[PASS] sob: {cls} state={st} allow={allow}")
        npass += 1
    print(f"{npass}/{npass} knop_spec_gate cases pass")


if __name__ == "__main__":
    main()
