# Phase 2.27B build and runtime status

Overall: PASS WITH LIMITATIONS for isolated compilation/linking; BLOCKED
for numerical execution. No jobs submitted, released, canceled or replaced.
Both original captures remain held as JobHeldUser. No idle diagnostic
allocation was found; production and unrelated allocations were not used.
Project rules prohibit these diagnostics on login nodes, and this task
prohibits new submissions. No runtime test was attempted on the login node.

## Built, not installed

Dedicated root:

    /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604

Final build command:

    python3 unit_tests/python/build_phase227.py --repo /viper/ptmp1/sbrennan/SN2024pba/cmfgen_jfnk --destination /viper/ptmp1/sbrennan/SN2024pba/cmfgen_validation/phase227b_20261121_0604/build_v4

Normal is O3 with the production AVX-512/alignment flags. Debug is O0,
-g, -check bounds. Both link the revised solver object explicitly before
unchanged production archives. The full candidate executable, replay and IO
harness use that same revised object within each variant. Other CMFGEN
objects remain production-built; this is not a whole-program debug build.
Commands, compiler identity, source/static-input hashes and ldd output are
in each build_manifest.json; complete compiler/linker output is build.log.

Earlier build_v1 is retained with the diagnosed ledger-integer argument
mismatch. KN_CENSUS_ROW now accepts the actual 64-bit ledger kind. build_v2
is retained with an unsupported Makefile-derived -fpe-summary=none option;
that option was removed in build_v4. The final build has no source type
warnings or ignored options. Remaining output: the inherited -mkl deprecation
remark and the pre-existing executable-stack requirement in cmfgen_sub.o.
No compiler/linker security setting was weakened to hide that warning.
build_v3 is retained but superseded: the original ledger refresh counter
resets each great iteration. build_v4 uses a process-local diagnostic serial
for all related output identity, without changing physical solver counters.
The full build.sh was not run because it publishes installed executables.

Final source SHA256:

- plane/fg_j_cmf_nonmono_v1.f:
  5cbe33221f4d65e03e99140e63181591cf5d3567edd2599d61299613abcdbe14
- unit_tests/f90/drv_knop_replay.f90:
  dd4be2f10b1478b135ba0094ebe741d6817bcc7c2afe18da0280b6be731d666f
- unit_tests/f90/drv_phase227_io.f90:
  22bc29bf5daf0a3f7c9ffbc89041139472e6242cf18be4ec2f7288e46eab3cfb

| Binary | SHA256 |
|---|---|
| normal/cmfgen_diagnostic.exe | e39a61a7d1485519d4ed5d59b6f7ee70ff5faf6f1fbefb7ef75fd468ae5e5461 |
| debug/cmfgen_diagnostic.exe | 04b8261d2b3c2f90a7e87359b84a5598228f4b13a6cdfa7542608cb7e9e651af |
| normal/knray_replay | 23a3fc6f9d6146107a43a52d85b85f9a295bb5b9530c29485a07361fb5fc4b04 |
| debug/knray_replay | b18bf7207fdb845969e22b00d704c767a4e17c62a4f58e7fccfdb778e15298aa |
| normal/phase227_io | 49091db86f8fc5dc889d3e11c0a3c002967e8adf4afe0c1f80c4e9d76cafdecc |
| debug/phase227_io | 251dbfa982b560160cc0f1bf67ba1320d817cfccf0ee95ffc4234732cbc151d0 |

## Prepared actual-Fortran tests

The harness calls KN_REPLAY_LOAD, KN_CENSUS_ROW and KN_OMEGA_APPLIED from
the actual revised solver object, not a Python imitation of the writer.
The Python runner independently validates authoritative census codes,
ray completeness, duplication, finite fields and aggregate counts; it
never invents a second convergence classifier.

Omega cases: valid complete two-channel history; absent file; missing pass;
duplicate pass; wrong ray/channel/direction; NaN; infinity; malformed numeric
record; mismatched decimal/IEEE binary64 encoding; extra ambiguous record.
Census cases: complete four-ray synthetic format fixture; all four requested
classes; missing ray; duplicate ray; invalid class; malformed row; integer
and floating-point format boundaries; exact aggregate reconciliation and
mismatch rejection. This is a format/parser test, NOT an 84-ray model census
or a test of the physical correctness of a classification.

Use io_tests_v4_normal/test_plan.json and io_tests_v4_debug/test_plan.json.
Earlier io_tests_normal/debug plans are retained but superseded: inspection
found their nominal NaN case contained a malformed numeric suffix. The v2
plan uses an actual NaN token and bit encoding. No version was executed.
Fixture preparation records the binary and fixture hashes and exact command.
The runner refuses numerical execution without SLURM_JOB_ID. The existence
of an allocation still requires permission to use it; no production job
may be borrowed for this test.

## Checks and unavailable results

- Fortran compile/link: PASS WITH LIMITATIONS above, both variants.
- Python syntax, black (24.4.2), isort (5.13.2): checked separately from
  numerical runtime tests; receipts/commands indexed with the artifacts.
- ruff and fprettify were not found in the current environment. No dependency
  or project security policy was changed to obtain them.
- git diff --check: PASS at the recorded verification point.
- Actual Fortran IO execution: NOT RUN / BLOCKED on an isolated allocation.
- Ray 84 control: NOT RUN. Historical 3.487919e-8 is a reference, not a new
  measurement. Planned printed-value tolerance: absolute 5e-15, zero relative
  tolerance; retain full-precision fields for the subsequent comparison.
- O0/O3 bytewise end fields: NOT RUN. The new builds vary solver optimization,
  unlike the historical driver-only comparison, and may legitimately expose
  differences that must be investigated rather than silently tolerated.
- v2 framing/content runtime regression: NOT RUN; deck reader unchanged.
- Rays 81/82 identical-stage trajectories and M/P percentile metrics:
  AWAITING CAPTURE and allocated replay execution.

No runtime result is inferred from compilation, fixture generation or
successful metadata-hash validation.
