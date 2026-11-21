"""Build Phase 2.27 diagnostics outside installed module/library/executable paths."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from phase227_capture import command, save_json, sha256, utc_now


def build(repo, destination, variant):
    destination.mkdir(parents=True, exist_ok=False)
    modules = destination / "mod"
    modules.mkdir()
    flags = ["-qopenmp", "-mkl", "-nbs", "-diag-disable=10448", "-traceback"]
    if variant == "normal":
        flags += [
            "-O3",
            "-mavx512f",
            "-mavx512dq",
            "-mavx512bw",
            "-mavx512vl",
            "-mavx512cd",
            "-align",
            "array64byte",
        ]
    else:
        flags += ["-O0", "-g", "-check", "bounds"]
    include = ["-I" + str(repo / "mod"), "-module", str(modules)]
    solver = repo / "plane/fg_j_cmf_nonmono_v1.f"
    driver = repo / "unit_tests/f90/drv_knop_replay.f90"
    harness = repo / "unit_tests/f90/drv_phase227_io.f90"
    objects = [destination / "solver.o", destination / "driver.o"]
    commands = [
        [
            "ifort",
            *flags,
            *include,
            "-fixed",
            "-extend-source",
            "-c",
            str(solver),
            "-o",
            str(objects[0]),
        ],
        [
            "ifort",
            *flags,
            "-I" + str(modules),
            *include,
            "-c",
            str(driver),
            "-o",
            str(objects[1]),
        ],
    ]
    libraries = [str(repo / "lib/libsubs.a"), str(repo / "lib/libmod.a")]
    commands.append(
        [
            "ifort",
            *flags,
            *map(str, objects),
            *libraries,
            "-o",
            str(destination / "knray_replay"),
        ]
    )
    if harness.exists():
        commands.append(
            [
                "ifort",
                *flags,
                "-I" + str(modules),
                *include,
                str(harness),
                str(objects[0]),
                *libraries,
                "-o",
                str(destination / "phase227_io"),
            ]
        )
    commands.append(
        [
            "ifort",
            *flags,
            str(repo / "new_main/cmfgen.o"),
            str(repo / "new_main/cmfgen_sub.o"),
            str(objects[0]),
            str(repo / "lib/libsubs.a"),
            str(repo / "lib/libpgplt.a"),
            str(repo / "lib/libmod.a"),
            str(repo / "lib/libpgplot_stubs.a"),
            "-lgfortran",
            "-o",
            str(destination / "cmfgen_diagnostic.exe"),
        ]
    )
    metadata = {
        "variant": variant,
        "generated_utc": utc_now(),
        "generation_command": sys.argv,
        "commands": commands,
        "source_commit": command(["git", "rev-parse", "HEAD"], repo).strip(),
        "compiler": command(["ifort", "--version"]),
        "flags": flags,
        "fma_policy": "compiler default; no additional override",
        "sources": {
            str(p): sha256(p)
            for p in [solver, driver] + ([harness] if harness.exists() else [])
        },
        "linked_inputs": {
            str(p): sha256(p)
            for p in [repo / "new_main/cmfgen.o", repo / "new_main/cmfgen_sub.o"]
            + list((repo / "lib").glob("*.a"))
        },
        "scope": "new solver and drivers only; remaining objects are installed production objects",
        "status": "BUILDING",
    }
    try:
        with (destination / "build.log").open("x") as log:
            for argv in commands:
                log.write(json.dumps(argv) + "\n")
                log.flush()
                subprocess.run(argv, cwd=modules, stdout=log, stderr=log, check=True)
        metadata["executables"] = {
            str(p): sha256(p)
            for p in destination.iterdir()
            if p.name in {"knray_replay", "phase227_io", "cmfgen_diagnostic.exe"}
        }
        metadata["dynamic_libraries"] = command(
            ["ldd", str(destination / "cmfgen_diagnostic.exe")]
        )
        metadata["status"] = "BUILT_NOT_RUNTIME_TESTED"
    except subprocess.CalledProcessError:
        metadata["status"] = "FAILED"
        raise
    finally:
        save_json(destination / "build_manifest.json", metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--variant", choices=("normal", "debug", "all"), default="all")
    args = parser.parse_args()
    variants = ("normal", "debug") if args.variant == "all" else (args.variant,)
    for variant in variants:
        sys.stdout.write(f"[info] Building isolated {variant} diagnostics\n")
        sys.stdout.flush()
        build(args.repo.resolve(), args.destination.resolve() / variant, variant)
    sys.stdout.write("[ok] Builds complete; no executable was installed or executed\n")


if __name__ == "__main__":
    main()
