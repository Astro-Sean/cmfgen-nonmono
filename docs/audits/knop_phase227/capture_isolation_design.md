# Phase 2.27B isolated capture design

Status: PREPARATION ONLY. No new scheduler submission is authorized. Jobs
13464217 and 13464311 remain held; their scripts and executable references
are unchanged. Releasing them would still invoke the old shared-directory
capture, NOT the candidates described here.

## Locations and preservation

Dedicated external root:

    /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604

- preservation/: first attempt; retained after scontrol rejected a comma-
  separated job argument. This is NOT the complete preservation receipt.
- preservation_v2/: complete initial source snapshot, full tracked binary
  diff, untracked-file status, source and executable hashes, and separate
  scontrol records for both held jobs. The untracked replay source is
  copied in full because git diff cannot represent it against HEAD.
- before_doc_updates/: a second, separately timestamped preservation,
  including pre-existing Phase 2.27 documents before their updates.
- build_v1/: first isolated link, retained with the census ledger-kind
  warning. Not an acceptance build.
- build_v2/: census kind corrected; a pre-existing Makefile-derived debug
  flag was rejected by ifort and ignored. Retained, not the final candidate.
- build_v3/: retained prior candidates; ledger refresh identity repeats
  across great iterations, so they are not current capture candidates.
- build_v4/: final candidate builds, separate normal and debug module and
  object directories. No installed library, module or executable is used
  as an output destination.
- captures_v4/normal/ and captures_v4/debug/: separately generated candidate
  directories. No shared mutable files or cross-capture links are allowed.
- io_tests_normal/ and io_tests_debug/: prepared fixture plans only.

Known pre-existing related diagnostic locations, all left untouched:

- /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/m15_2p82d_KNOP_diagdump_i/cmfgen/2.82days
  and its SCRATCH subdirectory, including the three large v2 raw decks.
- /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/docs/audits/knop_phase226/replay_out
- /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/docs/audits/knop_phase227/omsweep
- /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk/docs/audits/knop_phase227/replay_o3

The candidate preparation receipt records 251 copied input files per
capture with identical normal/debug input hashes. Manifest SHA256:
normal fe9cc2c53fecf26f4bb65c388c876d3b7af0fcea2e45b30d8f1c45d0ce203529;
debug f1e41e22cb83ecfc408a8ddea1d6b10c7aaf54f9303059088ba3ab6f8730a9b3.

Within EACH captures_v4/<variant>/ directory:

    bin/cmfgen_diagnostic.exe   read-only executable copy
    inputs/                    independent read-only regular-file snapshot
    provenance/                copied sources, build log, build manifest,
                               original model inventory and frozen runner
    manifest.json              pre-submission manifest, never job-modified
    manifest_seal.json         SHA256 of pre-submission manifest
    scratch/                   private mutable scratch
    work/                      private mutable model working files
    output/                    stdout/stderr and separate runtime manifest
    build/                     reserved capture-specific build metadata
    run_prepared.sh             pinned runner and explicit environment
    run_prepared.slurm          absolute paths; template only, NOT submitted

Raw KNAPPLY, KNOMG, KNSTAGE, KNCELL and KNRAY_CENSUS files remain in their
capture's work directory, which is their diagnostic output location in the
manifest. They are never moved over a different run's output.

## Snapshot guarantees and limits

The preparation helper copies bytes rather than hard-linking data. Atomic
input symlinks are dereferenced into private regular files. Source and
copy are hashed before/after copying; copies are mode 0444, executable
0555, and snapshot directories 0555. Normal/debug input hashes must agree
exactly. No existing destination is reused or overwritten.

These are read-only, hash-sealed snapshots, NOT filesystem-enforced WORM
storage. The owner can chmod and change them. The runner rejects seal,
executable, input, copied-provenance and dynamic-library hash mismatches.
Static archive and main-object hashes are recorded at build time and
rechecked during preparation. Compiled static code is embedded in each
private executable. Dynamic libraries remain external dependencies and
are hash-validated, not copied into each capture.

Source identity is HEAD e28833fe9d1286c9b92f909f559d92d28f4f7cd1 PLUS the
actual dirty source hashes. HEAD alone does not identify these binaries.
The final build is incremental and isolated: the revised solver object is
linked before the unchanged archives, together with existing main
objects. It is not a clean source rebuild of every archived routine.
Debug means O0/bounds checks in the revised solver/driver, NOT an entirely
debug-compiled CMFGEN. The same revised solver object supplies the full
candidate executable, replay driver and Fortran IO harness per variant.

## Manifests and runtime validation

capture_manifest_schema.json describes the pre-submission document.
Each contains source hashes, compiler identity and exact commands/flags,
linked static-object hashes, dynamic-library paths/hashes, input identity,
absolute paths, diagnostic configuration and generation command/time.
The frozen runtime runner creates output/runtime_manifest.json separately
with job ID, node, start/end time, relevant environment, actual hashes,
working-input hashes and process status. A successful process exit is NOT
scientific acceptance. Invalid provenance is status INVALID and nonzero.
No existing runtime manifest or log is overwritten.

The prepared candidates are deliberately PREPARED_NOT_VALIDATED and
submission_allowed=false. The runtime runner rejects this readiness state
even if invoked manually inside an allocation. A later approved capture
requires a NEW reviewed manifest/snapshot; neither job nor operator should
edit a pre-submission manifest while it is running.

The template requests interactive/interactive, 8 CPUs, 29G, 01:59:00.
It launches only an absolute private runner/executable. OMP and MKL thread
counts are eight, matching the intended original launcher/batch setup.
Loader/thread variables and their absence are pinned explicitly. Scratch
links, if later enabled after acceptance, point only to that capture's
absolute private scratch directory. No helper submits jobs or writes to
an existing model directory.

## Commands and remaining input-closure blocker

Build generation command:

    python3 unit_tests/python/build_phase227.py --repo /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk --destination /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604/build_v4

Candidate-generation command:

    python3 unit_tests/python/prepare_phase227_model.py --repo /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk --source-model /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/m15_2p82d_KNOP_diagdump_i/cmfgen/2.82days --build-root /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604/build_v4 --destination /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604/captures_v4

Exact receipts and output hashes are in preparation_result.json and the
artifact index. Candidates copy the existing prepared VADAT and static
inputs without running prepare_VADAT or modifying the originals. The
old directory contains generated MODEL/RVTJ, population/checkpoint data,
mutable IN_ITS and scratch outputs from an earlier execution. The helper
records every top-level entry, selected candidate input, symlink target,
exclusion and missing optional input. This is NOT a complete read-set
trace. In particular, cold-start versus warm-restart scratch/population
roles and effective IN_ITS stop policy remain unreviewed. The candidates
must not be submitted as a production-faithful restart until this is
resolved. No delayed background edit of IN_ITS is included.

## Submission prerequisites

1. Pass allocated actual-writer/reader tests and control regressions.
2. Certify the complete input/restart read set and intended stop policy.
3. Review diagnostic stage gaps and missing census fields.
4. Verify field-stage agreement, actual applied schedules, and provenance.
5. Obtain separate approval for any scheduler action. The two old jobs
   remain held regardless of candidate preparation progress.
