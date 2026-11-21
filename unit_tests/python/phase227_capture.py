"""Prepare hash-sealed diagnostic inputs without submitting scheduler jobs."""

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command(args, cwd=None):
    return subprocess.check_output(args, cwd=cwd, text=True)


def save_json(path, data):
    with Path(path).open("x", encoding="ascii") as stream:
        json.dump(data, stream, indent=2, ensure_ascii=True)
        stream.write("\n")


def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def index_artifacts(args):
    repo = args.repo.resolve(strict=True)
    paths = [path for path in args.external_root.resolve().rglob("*") if path.is_file()]
    paths += [repo / name for name in args.files]
    save_json(
        args.destination,
        {
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "source_commit": command(["git", "rev-parse", "HEAD"], repo).strip(),
            "python": sys.version,
            "hostname": socket.gethostname(),
            "machine": list(os.uname()),
            "files": {str(path): sha256(path) for path in sorted(set(paths))},
            "source_locator_note": "Routine anchors refer to the hashed source, not the pre-2.27B line numbering",
            "scheduler": {
                job: command(["scontrol", "show", "job", job])
                for job in ("13464217", "13464311")
            },
        },
    )


def preserve(args):
    repo = args.repo.resolve(strict=True)
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    sources = [
        "plane/fg_j_cmf_nonmono_v1.f",
        "unit_tests/f90/drv_knop_replay.f90",
        "new_main/mod_subs/rd_control_variables.f",
        "new_main/mod_subs/control_variable_mod.f",
        "Makefile_definitions",
    ]
    sources += [
        str(path.relative_to(repo))
        for path in sorted((repo / "docs/audits/knop_phase227").glob("*"))
        if path.is_file() and path.suffix in {".md", ".json", ".csv", ".tsv", ".py"}
    ]
    sources.append("docs/audits/CMFGEN_PHASE227_ADAPTIVE_OMEGA_AND_RAY_CENSUS.md")
    for name in sources:
        target = destination / "source" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(repo / name, target)
    for name, argv in {
        "git_status.txt": ["git", "status", "--short", "--untracked-files=all"],
        "tracked.diff": ["git", "diff", "--binary", "HEAD"],
        "scheduler_13464217.txt": ["scontrol", "show", "job", "13464217"],
        "scheduler_13464311.txt": ["scontrol", "show", "job", "13464311"],
    }.items():
        (destination / name).write_text(command(argv, repo), encoding="ascii")
    binary_paths = [repo / "exe/cmfgen_dev.exe", repo / "exe/cmfgen_dev_v2dump.exe"]
    save_json(
        destination / "preservation.json",
        {
            "timestamp": utc_now(),
            "generation_command": sys.argv,
            "source_commit": command(["git", "rev-parse", "HEAD"], repo).strip(),
            "sources": {name: sha256(repo / name) for name in sources},
            "executables": {str(path): sha256(path) for path in binary_paths},
            "untracked_source_note": "Replay driver copied in full; tracked.diff excludes untracked files",
        },
    )


def copy_sealed(source, target):
    source = source.resolve(strict=True)
    if not source.is_file():
        raise ValueError(f"Input is not a regular file: {source}")
    if target.exists():
        raise ValueError(f"Refusing existing destination: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    before = sha256(source)
    shutil.copyfile(source, target)
    if sha256(source) != before or sha256(target) != before:
        raise ValueError(f"Input changed during snapshot: {source}")
    target.chmod(0o444)
    return {"path": str(target), "sha256": before, "original_path": str(source)}


def prepare(args):
    spec = json.loads(args.spec.read_text(encoding="ascii"))
    required = {
        "capture_id",
        "intended_job_type",
        "source_commit",
        "source_hashes",
        "compiler",
        "compiler_flags",
        "linked_libraries",
        "inputs",
        "environment",
        "diagnostic_configuration",
        "executable",
        "provenance_files",
    }
    if required - spec.keys():
        raise ValueError(
            f"Missing specification keys: {sorted(required - spec.keys())}"
        )
    root = args.destination.resolve()
    root.mkdir(parents=True, exist_ok=False)
    records = {}
    for name, original in spec["inputs"].items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe input name: {name}")
        records[name] = copy_sealed(Path(original), root / "inputs" / relative)
    executable = copy_sealed(
        Path(spec["executable"]), root / "bin/cmfgen_diagnostic.exe"
    )
    Path(executable["path"]).chmod(0o555)
    for directory in ("scratch", "work", "output", "build"):
        (root / directory).mkdir()
    provenance = {}
    for name, source in spec["provenance_files"].items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe provenance name: {name}")
        provenance[name] = copy_sealed(Path(source), root / "provenance" / relative)
    manifest = dict(spec)
    manifest["provenance_files"] = provenance
    manifest.update(
        {
            "schema_version": 1,
            "generated_utc": utc_now(),
            "generation_command": sys.argv,
            "capture_root": str(root),
            "executable": executable,
            "inputs": records,
            "scratch_path": str(root / "scratch"),
            "working_path": str(root / "work"),
            "output_path": str(root / "output"),
            "omega_path": str(root / "work"),
            "census_path": str(root / "work"),
            "submission_allowed": False,
            "readiness": "PREPARED_NOT_VALIDATED",
            "immutability": "read-only copies plus SHA256; not filesystem-enforced immutable",
        }
    )
    save_json(root / "manifest.json", manifest)
    (root / "manifest.json").chmod(0o444)
    save_json(root / "manifest_seal.json", {"sha256": sha256(root / "manifest.json")})
    (root / "manifest_seal.json").chmod(0o444)
    for name in ("inputs", "provenance", "bin"):
        directory = root / name
        for path in directory.rglob("*"):
            if path.is_dir():
                path.chmod(0o555)
        directory.chmod(0o555)


def validate(manifest_path):
    path = manifest_path.resolve(strict=True)
    root = path.parent
    seal = json.loads((root / "manifest_seal.json").read_text(encoding="ascii"))
    if sha256(path) != seal["sha256"]:
        raise ValueError(f"Manifest seal mismatch: {path}")
    manifest = json.loads(path.read_text(encoding="ascii"))
    records = {
        "executable": manifest["executable"],
        **{f"input:{key}": value for key, value in manifest["inputs"].items()},
        **{
            f"provenance:{key}": value
            for key, value in manifest["provenance_files"].items()
        },
    }
    actual = {}
    for name, record in records.items():
        target = Path(record["path"])
        if not target.is_absolute() or not target.resolve().is_relative_to(root):
            raise ValueError(f"Unisolated path: {target}")
        actual[name] = sha256(target)
        if actual[name] != record["sha256"]:
            raise ValueError(f"Hash mismatch: {target}")
    return manifest, actual


def runtime(args):
    started = utc_now()
    output = args.manifest.resolve().parent / "output" / "runtime_manifest.json"
    result = {
        "start_utc": started,
        "hostname": socket.gethostname(),
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "environment": {},
        "generation_command": sys.argv,
        "status": "INVALID",
    }
    process_status = 1
    try:
        manifest, actual = validate(args.manifest)
        result["actual_hashes"] = actual
        for name, record in manifest["linked_libraries"].items():
            if sha256(record["path"]) != record["sha256"]:
                raise ValueError(f"Dynamic library mismatch: {name}")
        if manifest["readiness"] != "READY_FOR_CAPTURE":
            raise ValueError(
                "Capture has not passed input-closure and runtime acceptance review"
            )
        environment = manifest["environment"]
        result["environment"] = {key: os.environ.get(key) for key in environment}
        for key, expected in environment.items():
            if os.environ.get(key) != expected:
                raise ValueError(f"Environment mismatch: {key}")
        for key in manifest.get("environment_absent", []):
            if key in os.environ:
                raise ValueError(f"Unexpected inherited environment: {key}")
        if not os.environ.get("SLURM_JOB_ID"):
            raise ValueError("Runtime execution requires an allocated scheduler job")
        if not args.execute:
            raise ValueError(
                "Execution requires explicit --execute; preparation never runs a model"
            )
        work = Path(manifest["working_path"])
        for name, record in manifest["inputs"].items():
            target = work / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() or target.is_symlink():
                raise ValueError(f"Working directory is not fresh: {target}")
            shutil.copyfile(record["path"], target)
        result["working_input_hashes"] = {
            name: sha256(work / name) for name in manifest["inputs"]
        }
        for name, record in manifest["inputs"].items():
            if result["working_input_hashes"][name] != record["sha256"]:
                raise ValueError(f"Working input mismatch: {name}")
        for name in manifest.get("runtime_dirs", []):
            if Path(name).name != name:
                raise ValueError(f"Unsafe runtime dir: {name}")
            (work / name).mkdir()
        for name in manifest.get("scratch_files", []):
            if Path(name).name != name:
                raise ValueError(f"Unsafe scratch name: {name}")
            (work / name).symlink_to(Path(manifest["scratch_path"]) / name)
        with (output.parent / "stdout.log").open("x") as stdout:
            with (output.parent / "stderr.log").open("x") as stderr:
                process = subprocess.run(
                    [manifest["executable"]["path"]],
                    cwd=work,
                    stdout=stdout,
                    stderr=stderr,
                    check=False,
                )
        process_status = process.returncode
        validate(args.manifest)
        result["status"] = "PROCESS_EXITED" if process_status == 0 else "FAILED"
        result["returncode"] = process_status
        result["scientific_validity"] = "NOT_ASSESSED"
    except (OSError, ValueError, KeyError) as error:
        process_status = 1
        result["error"] = str(error)
    finally:
        result["end_utc"] = utc_now()
        save_json(output, result)
    return process_status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    snapshot = actions.add_parser("preserve")
    snapshot.add_argument("--repo", type=Path, required=True)
    snapshot.add_argument("--destination", type=Path, required=True)
    preparation = actions.add_parser("prepare")
    preparation.add_argument("--spec", type=Path, required=True)
    preparation.add_argument("--destination", type=Path, required=True)
    check = actions.add_parser("validate")
    check.add_argument("manifest", type=Path)
    run = actions.add_parser("runtime")
    run.add_argument("manifest", type=Path)
    run.add_argument("--execute", action="store_true")
    index = actions.add_parser("index")
    index.add_argument("--repo", type=Path, required=True)
    index.add_argument("--destination", type=Path, required=True)
    index.add_argument("--external-root", type=Path, required=True)
    index.add_argument("files", nargs="+")
    args = parser.parse_args()
    try:
        if args.action == "preserve":
            preserve(args)
        elif args.action == "prepare":
            prepare(args)
        elif args.action == "validate":
            validate(args.manifest)
        elif args.action == "index":
            index_artifacts(args)
        else:
            return runtime(args)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.stderr.write(f"[FAIL] {error}\n")
        return 1
    sys.stdout.write(f"[ok] {args.action} completed\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
