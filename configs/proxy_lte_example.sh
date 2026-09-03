#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# proxy_lte_example.sh -- Conservative LTE baseline (Mode A) configuration.
#
# This is the DEFAULT and the reference against which every other mode is
# compared.  It is the safest starting point for CMFGEN:
#
#   * Local physical gas temperature is written to the _IN T column (never
#     a synthetic T_ion).  Verified: the species _IN T column is only an
#     interpolation coordinate for REGRID_LOG_DC_V1; the physical T used by
#     LTEPOP_WLD_V2 comes from T_IN (set_new_model_estimates.f).
#   * DC = 1.0 everywhere (full LTE), written LINEAR because a modern format
#     date => TAKE_LOGS=.TRUE. in regrid_log_dc_v1.f:127-134.
#   * DI = next-higher ion-stage density from a log-space Saha chain using the
#     validated GAMMAS_IN electron density (set_new_model_estimates.f:365-375).
#   * No grey radiation correction, no median smoothing, no abundance fallback.
#   * Numerical floor (1e-30) is a safeguard only; it is reported and excluded
#     from abundance / charge closure.
#   * Round-trip is exact (target == CMFGEN LTE population, DC=1).
#
# Output is written to a NEW directory; the source model is never overwritten
# unless --force-overwrite-source is given explicitly.
#
# Required inputs in <model_dir>:
#   SN_HYDRO_DATA, VADAT (or VADAT_INIT), GAMMAS_IN, MODEL_SPEC, *_F_OSCDAT,
#   *_F_TO_S (only when NF != NS for a species).
# ---------------------------------------------------------------------------
set -euo pipefail

MODEL_DIR="${1:?usage: proxy_lte_example.sh <model_dir> [output_dir]}"
OUTPUT_DIR="${2:-${MODEL_DIR}/proxy_lte_out}"

python scripts/proxy_populations_v2.py \
    "$MODEL_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --strategy LTE \
    --write-levels NF \
    --min-ion-frac 1e-30 \
    --t-exc-min 0.3

# Exit status: 0 = all species round-trip OK; 2 = round-trip failed for some
# species (do NOT use the output for a production CMFGEN run in that case).
#
# Inspect, in $OUTPUT_DIR:
#   proxy_validation.txt  -- per-species DI/DC ranges, omitted-stage risk,
#                            recomputed N_e vs input, abundance residuals.
#   proxy_audit.json      -- machine-readable diagnostics.
#   proxy_strategy.json   -- selected strategy + abundance provenance.
#
# If proxy_validation.txt reports "highest modelled stage carries >50%" for a
# species, the highest included ion stage is NOT a safe sink: add the next
# stage to MODEL_SPEC or accept that CMFGEN's first iteration will move a lot
# of population out of that stage.
