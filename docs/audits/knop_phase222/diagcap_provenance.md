# Part A - diagcap job 13429272 provenance (snapshot 2026-11-20 22:45 CET)

State: PENDING (QOSGrpCpuLimit) - never started; no per-pass data
exists. Preserved, not modified. The run dir's lone OUTGEN is a
staging error log, not solver output.

## Identity (manifest.json + on-disk hashes, verified this session)

| item | value |
|---|---|
| run dir | cmfgen_validation/m15_2p82d_KNOP_diagcap |
| deck | clone of B1_PREFLAG + KN_PSWEEP_DIAG=T, KN_PSWEEP_LS=-1 |
| cmfgen_dev.exe | pinned exe_frozen/93364cace2fabdd9_484f34b19ee9e205 (nested doubled path resolves; sha256 93364cac e2fabdd9...) |
| cmf_flux.exe | sha256 484f34b1 9ee9e205... |
| source commit | 217049ee (pre-flag era, same binary as A1/A2/B1) |
| VADAT_INIT | sha256 f238af97 9b9613d3... (identical to A1/B1/B2) |
| SN_HYDRO_DATA | sha256 347eda19 ce4f645c... (identical) |
| SCRTEMP | sha256 5704dc48 6d6f5a16... (iso_M it157 donor, identical to B1/B2) |
| FGINEXT_GS | sha256 48077bb3 65c20760... (identical) |
| KN_MAT_FREEZE | T (inherited - frozen matter, same limitation as A/B arms) |
| KN_SEED_ARM | 3 (arm A/A control: retained-store seed on probe ray 83) |
| KN_WATCH | ray 83, K=2, slot 1436 (probe telemetry, same as A1) |
| KN_PSWEEP_DIAG | T; KN_PSWEEP_LS=-1 (arms ALL rays) |

## What it will measure when it runs

KN_PSWEEP_DIAG (phase-5.3K instrumentation, confirmed present in the
93364cac binary via strings) records per ray, per completed pass:
NPO, OMEGA_RAY used, production RES, R_MAP (complete-map defect
max|Y-X|/SCALE), RMS_MAP, argmax ray-point J and frequency L,
NBAD/NTOT census, floor/non-finite counts, cap/floor activations
(NCP_P/NNG_P), sign(GAM(argmax)) and NFM_L(argmax) - the donor slot
offset at the worst cell. Written to KN_PSWEEP_DIAG.dat.

This answers per-ray: per-pass residual history, terminal residual,
argmax (J,L), GAM sign and donor offset at the argmax cell - the
residual-attribution record the A1 telemetry lacks.

Not recorded even by diagcap: per-cell CHI/CHI_HAT, donor indices
JPR/IKC/WINX values, per-pass field dumps (KN_PSWEEP_SNAP off,
KN_PSWEEP_SNHI=0). The working tree's KN_RDUMP_LS/LS2 replay-deck
diagnostic (phase-2.22, serializes ray geometry/GAM/CHI/ETA cubics/
store columns/seeded field for an offline driver) is NOT in the
93364cac binary - it exists only in the current live build
(d50dcacd...). A future replay run needs that binary, not diagcap's.

## Termination

Not applicable - job has not run.
