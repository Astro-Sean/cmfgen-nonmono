# REPOSITORY_HYGIENE_REPORT.md (Proxy Populations)

**Date:** 2026-08-21
**Scope:** Repository hygiene for the proxy populations validation work.

---

## 1. New Files Created This Session

### Source code and tests
| File | Location | Purpose |
|---|---|---|
| `AGENTS_PROXY_POPULATIONS.md` | project root | Proxy populations project guidance |
| `test_independent_in_parser.py` | tests/ | Independent _IN file parser/validator |
| `build_manifests.py` | tests/ | Manifest builder for atomic/model/grid |
| `compare_out_files.py` | scripts/ | OUT file structural comparison |

### Generated outputs (runs/)
| File | Location | Purpose |
|---|---|---|
| `SN_12C_0.5d/` | runs/proxy_validation/ | Proxy _IN files + diagnostics |
| `cmfgen_startup_test/` | runs/proxy_validation/ | CMFGEN startup test dir |
| `reference/` | runs/proxy_validation/ | Reference manifests |

### Manifests and reports (runs/proxy_validation/)
| File | Purpose |
|---|---|
| `SN_12C_0.5d/ATOMIC_DATA_MANIFEST.json` | Atomic data checksums (65 files) |
| `SN_12C_0.5d/MODEL_ATOM_MANIFEST.json` | Model atom (10 species, 80 NF) |
| `SN_12C_0.5d/INPUT_GRID_MANIFEST.json` | Input grid (ND=400, outer-to-inner) |
| `SN_12C_0.5d/independent_parse_report.json` | Independent parser results |
| `OUT_FILE_INVENTORY.csv` | Proxy + reference file inventory |
| `OUT_FILE_COMPARISON.csv` | Shared species structural comparison |
| `OUT_FILE_COMPARISON.md` | Comparison report |

### Deliverable documents (docs/)
| File | Purpose |
|---|---|
| `PROXY_POPULATIONS_AUDIT.md` | Consolidated audit |
| `PROXY_FORMAT_SPECIFICATION.md` | _IN format spec |
| `PROXY_ROUNDTRIP_REPORT.md` | Round-trip validation |
| `REFERENCE_COMPARABILITY_REPORT.md` | Reference comparability |
| `POPULATION_NORMALIZATION_REPORT.md` | Normalization validation |
| `FRESH_MODEL_CONVERGENCE_REPORT.md` | CMFGEN startup test |

### Reference data (data/reference/)
| Directory | Purpose |
|---|---|
| `M15-He-O-Ca-evo-new-4-all/33/` | Downloaded reference files (34 files + checksums) |

---

## 2. Separation of Concerns (VERIFIED)

| Category | Location | Status |
|---|---|---|
| Raw atomic data | /u/sbrennan/atomic/ (symlinked) | OK (external) |
| CMFGEN source | cmfgen_source/ | OK (immutable) |
| CMFGEN patches | cmfgen/patches/ | OK (empty — no patches needed) |
| MODEL_SPEC/VADAT templates | cmfgen_base/ | OK |
| Raw hydrodynamic input | SN_12C/cmfgen/0.5days/SN_HYDRO_DATA | OK |
| Generated proxy _IN files | runs/proxy_validation/SN_12C_0.5d/ | OK (outside source) |
| Generated _OUT files | runs/proxy_validation/cmfgen_startup_test/ | OK (outside source) |
| Reference files | data/reference/M15-He-O-Ca-evo-new-4-all/33/ | OK (outside source) |
| Comparison products | runs/proxy_validation/ | OK |
| Logs | runs/proxy_validation/cmfgen_startup_test/batch.log | OK |
| Manifests | runs/proxy_validation/SN_12C_0.5d/ | OK |

No generated files were placed in source directories.

---

## 3. Untracked Files (VERIFIED)

All proxy-related files are currently untracked in git:
- `scripts/proxy_populations.py` (v1)
- `scripts/proxy_populations_v2.py` (v2)
- `scripts/test_proxy_populations_v2.py`
- `scripts/compare_proxy_v1_v2.py`
- `scripts/plot_proxy_population_comparison.py`
- `scripts/proxy_audit_report.md`
- `configs/proxy_grey_example.sh`
- `configs/proxy_lte_example.sh`

**Recommendation:** Commit the proxy v2 script, tests, and deliverables once the
orchestrator integration is complete. The v1 script should be archived.

---

## 4. Duplicate/Stale Files (VERIFIED)

| File | Status | Action |
|---|---|---|
| `scripts/proxy_populations.py` (v1) | Superseded by v2 | Archive to scripts/deprecated/ |
| `scripts/proxy_audit_report.md` | Superseded by docs/PROXY_POPULATIONS_AUDIT.md | Keep as historical reference |
| `docs/proxy_populations_audit.md` | Prior detailed audit | Keep (cited by new audit) |
| `docs/proxy_v1_v2_comparison.md` | Prior v1/v2 comparison | Keep (historical) |
| `docs/proxy_validation_roundtrip_report.md` | Prior roundtrip report | Keep (cited by new report) |
| `docs/proxy_unresolved_assumptions.md` | Prior assumptions | Keep (updated by new audit) |
| `docs/proxy_changelog.md` | Prior changelog | Keep (historical) |

No stale binaries or broken symlinks found in the proxy validation directories.

---

## 5. Provenance (VERIFIED)

### Reference files
- URL: https://ttt.astro.su.se/~maer0651/M15-He-O-Ca-evo-new-4-all/33/
- Retrieval date: 2026-08-21
- 34 files downloaded with SHA-256 checksums
- Stored in data/reference/M15-He-O-Ca-evo-new-4-all/33/

### Atomic data
- 65 files symlinked to /u/sbrennan/atomic/ (19apr23 versions for most species)
- H/He use 5dec96 versions
- All checksums in ATOMIC_DATA_MANIFEST.json

### CMFGEN binary
- Path: cmfgen_source/exe/cmfgen_dev.exe
- Build date: 2026-08-21
- HDKT: 4.79924307 (computed from physical constants)

---

## 6. Hygiene Assessment

| Check | Status |
|---|---|
| No generated files in source directories | PASS |
| Reference files in data/reference/ | PASS |
| Generated outputs in runs/ | PASS |
| Manifests with checksums | PASS |
| No duplicate active configurations | PASS |
| No stale binaries | PASS |
| No broken symlinks in proxy dirs | PASS |
| Provenance recorded | PASS |
| Untracked files documented | PASS (recommend commit) |

**Hygiene status: PASS**
