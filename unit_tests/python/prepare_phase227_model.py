"""Prepare separate capture candidates; unreviewed input closure blocks execution."""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from phase227_capture import command, prepare, save_json, sha256, utc_now, validate

INPUT_NAMES = {
    "ADJUST_R_DEFAULTS",
    "CFDAT_IN",
    "CUR_MODEL_DATA",
    "EDDFACTOR_INFO",
    "FIN_CAL_GRID",
    "GAMMAS_IN",
    "IN_ITS",
    "MODEL",
    "MODEL_SPEC",
    "NUC_DECAY_DATA",
    "RDINR",
    "RVTJ",
    "SN_HYDRO_DATA",
    "SPECIES_MASSES",
    "T_IN",
    "VADAT",
    "VADAT_INIT",
    "arnaud_rothenflug.dat",
    "ca2_col_data_filtered",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source-model", type=Path, required=True)
    parser.add_argument("--build-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    repo, model, destination = (
        args.repo.resolve(),
        args.source_model.resolve(),
        args.destination.resolve(),
    )
    destination.mkdir(parents=True, exist_ok=False)
    inputs, excluded, missing = (
        {},
        {},
        sorted(name for name in INPUT_NAMES if not (model / name).is_file()),
    )
    inventory = []
    for path in sorted(model.iterdir()):
        record = {"path": str(path), "symlink": path.is_symlink()}
        if path.is_symlink():
            record["link_text"] = os.readlink(path)
            record["resolved_path"] = str(path.resolve())
        if path.name in INPUT_NAMES or (
            path.is_symlink() and not path.resolve().is_relative_to(model)
        ):
            if not path.is_file():
                raise ValueError(f"Unresolved candidate input: {path}")
            inputs[path.name] = str(path)
            record["selected"] = True
        else:
            excluded[path.name] = (
                "not in reviewed static-input candidates; restart/output role unresolved"
            )
            record["selected"] = False
        inventory.append(record)
    save_json(
        destination / "model_inventory.json",
        {
            "source_model": str(model),
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "entries": inventory,
            "excluded": excluded,
            "missing_optional_candidates": missing,
            "input_closure": "NOT_CERTIFIED; cold versus restart seed and scratch inputs require review",
        },
    )
    manifests = {}
    for variant in ("normal", "debug"):
        build_root = args.build_root.resolve() / variant
        build = json.loads((build_root / "build_manifest.json").read_text())
        if build["status"] != "BUILT_NOT_RUNTIME_TESTED":
            raise ValueError(f"Incomplete build: {build_root}")
        if any(
            sha256(path) != digest for path, digest in build["linked_inputs"].items()
        ):
            raise ValueError(f"Linked inputs changed after build: {variant}")
        if any(sha256(path) != digest for path, digest in build["executables"].items()):
            raise ValueError(f"Executable changed after build: {variant}")
        if any(sha256(path) != digest for path, digest in build["sources"].items()):
            raise ValueError(f"Source/build hash mismatch: {variant}")
        libraries = {}
        for path in re.findall(r"(/[\w/+.=-]+)", build["dynamic_libraries"]):
            library = Path(path).resolve(strict=True)
            libraries[library.name] = {"path": str(library), "sha256": sha256(library)}
        provenance = {
            "build_manifest.json": str(build_root / "build_manifest.json"),
            "build.log": str(build_root / "build.log"),
            "model_inventory.json": str(destination / "model_inventory.json"),
            "phase227_capture.py": str(repo / "unit_tests/python/phase227_capture.py"),
        }
        provenance.update(
            {
                "source/" + str(Path(path).relative_to(repo)): path
                for path in build["sources"]
            }
        )
        root = destination / variant
        spec = {
            "capture_id": f"{destination.parent.name}_{destination.name}_{variant}",
            "intended_job_type": variant,
            "source_commit": build["source_commit"],
            "source_hashes": build["sources"],
            "source_dirty": True,
            "compiler": build["compiler"],
            "compiler_flags": build["flags"],
            "linked_libraries": libraries,
            "linked_static_inputs": build["linked_inputs"],
            "executable": str(build_root / "cmfgen_diagnostic.exe"),
            "inputs": inputs,
            "environment": {
                "OMP_NUM_THREADS": "8",
                "MKL_NUM_THREADS": "8",
                "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
                "LD_PRELOAD": "",
                "SCRATCH_DIR": str(root / "scratch"),
                "TMPDIR": str(root / "scratch"),
            },
            "diagnostic_configuration": {
                "dump_rays": [81, 82, 84],
                "unchanged_VADAT": True,
                "probe_ray": 83,
                "input_inventory_complete": False,
            },
            "provenance_files": provenance,
            "limitations": [
                "input role inventory is not a complete read-set trace",
                "debug checks cover the new solver/driver, not all archived objects",
                "allocated IO and replay acceptance tests not run",
            ],
        }
        inherited_keys = [
            "OMP_STACKSIZE",
            "OMP_PROC_BIND",
            "OMP_PLACES",
            "OMP_DYNAMIC",
            "OMP_SCHEDULE",
            "OMP_THREAD_LIMIT",
            "KMP_STACKSIZE",
            "KMP_AFFINITY",
            "KMP_HW_SUBSET",
            "MKL_DYNAMIC",
            "MKL_CBWR",
            "CMF_KN_RF_MB",
        ]
        # SLURM_JOB_ID is required by the runner itself, PWD/USER only feed
        # a log line in wr_pwd; none belong in the pinned environment.
        spec["environment"].update(
            {key: os.environ[key] for key in inherited_keys if key in os.environ}
        )
        spec["environment_absent"] = [
            key for key in inherited_keys if key not in os.environ
        ]
        spec["runtime_dirs"] = ["data"]
        spec["scratch_files"] = [
            "POINT1",
            "POINT2",
            "SCRTEMP",
            "STEQ_VALS",
            "J_COMP",
            "JEW",
            "EWDATA",
            "CFDAT_OUT",
            "CONT_FREQ",
            "LINEHEAT",
            "NETRATE",
            "TOTRATE",
            "TRANS_INFO",
            "ADIABAT_CHK",
        ]
        spec_path = destination / f"{variant}_spec.json"
        save_json(spec_path, spec)
        prepare(argparse.Namespace(spec=spec_path, destination=root))
        manifest, _ = validate(root / "manifest.json")
        manifests[variant] = manifest
        launcher = root / "run_prepared.sh"
        with launcher.open("x", encoding="ascii") as stream:
            stream.write("#!/bin/sh\nset -eu\n")
            for key in spec["environment_absent"]:
                stream.write(f"unset {key}\n")
            for key, value in spec["environment"].items():
                stream.write(f"export {key}={shlex.quote(value)}\n")
            runtime_command = [
                sys.executable,
                str(root / "provenance/phase227_capture.py"),
                "runtime",
                str(root / "manifest.json"),
                "--execute",
            ]
            stream.write(f"exec {shlex.join(runtime_command)}\n")
        launcher.chmod(0o555)
        slurm = root / "run_prepared.slurm"
        with slurm.open("x", encoding="ascii") as stream:
            stream.write("#!/bin/sh\n")
            stream.write(f"#SBATCH --job-name=phase227b_{variant}\n")
            stream.write("#SBATCH --partition=interactive\n#SBATCH --qos=interactive\n")
            stream.write(
                "#SBATCH --cpus-per-task=8\n#SBATCH --mem=29G\n#SBATCH --time=01:59:00\n"
            )
            stream.write(f"#SBATCH --chdir={root / 'work'}\n")
            stream.write(f"#SBATCH --output={root / 'output/slurm-%j.out'}\n")
            stream.write(f"#SBATCH --error={root / 'output/slurm-%j.err'}\n")
            stream.write(f"set -eu\nexec /bin/sh {shlex.quote(str(launcher))}\n")
        slurm.chmod(0o444)
    if {
        name: entry["sha256"] for name, entry in manifests["normal"]["inputs"].items()
    } != {
        name: entry["sha256"] for name, entry in manifests["debug"]["inputs"].items()
    }:
        raise ValueError("Normal/debug input snapshots differ")
    save_json(
        destination / "preparation_result.json",
        {
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "source_commit": command(["git", "rev-parse", "HEAD"], repo).strip(),
            "status": "PREPARED_NOT_VALIDATED",
            "submission_allowed": False,
            "input_files_per_capture": len(inputs),
            "input_snapshots_identical": True,
            "manifests": {
                name: {
                    "path": str(destination / name / "manifest.json"),
                    "sha256": sha256(destination / name / "manifest.json"),
                }
                for name in manifests
            },
        },
    )
    sys.stdout.write(
        "[ok] Isolated candidate snapshots prepared; input-closure review blocks execution\n"
    )


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.stderr.write(f"[FAIL] {error}\n")
        raise SystemExit(1) from error
