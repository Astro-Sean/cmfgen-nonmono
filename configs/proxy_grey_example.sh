#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# proxy_grey_example.sh -- Grey dilute-LTE DIAGNOSTIC mode (Mode B).
#
# This is an OPT-IN diagnostic strategy, NOT the default.  Use it only when you
# have evidence that the LTE baseline produces large first-iteration CMFGEN
# corrections in the optically-thin outer ejecta and want to test whether a
# dilute-LTE starting point reduces them.
#
# What it does:
#   * Computes a grey optical depth, geometric dilution W, T_rad and T_ion.
#   * Uses T_ion + W for the Saha *diagnostic* ionization balance.
#   * By DEFAULT still writes T_local (not T_ion) to the _IN T column, because
#     writing a synthetic T_ion is only safe once the round-trip proves the
#     species _IN T column cannot alter the physical model (audit B5).  Pass
#     --write-t-ion to override AFTER you have run the round-trip test.
#   * DC = target / LTE (target = Boltzmann-partitioned current-ion population),
#     written LINEAR.  Abundance closure is enforced and DC is recomputed after
#     closure, independently of the numerical floor.
#   * Reports the ratio of dilute-LTE to local-LTE populations per species/depth
#     in proxy_validation.txt, and auto-falls back to Mode A behaviour (DC=1,
#     T_local) wherever the grey correction exceeds the configured limits.
#   * Grey corrections go continuously to zero in the optically thick region
#     (W -> 1, T_ion -> T_local), so the deep ejecta is identical to Mode A.
#
# When NOT to use it:
#   * If you cannot verify the CMFGEN conventions (the selector refuses to run
#     production output in that case -- use Mode A).
#   * As a replacement for CMFGEN's NLTE solution -- it is a preconditioner.
#   * If any species has an omitted-stage fraction above your threshold: the
#     selector warns and excludes that species from automatic NLTE enhancement.
#
# Damping is used ONLY in the grey iteration and convergence of both N_e and
# opacity is reported.  The final N_e is checked for consistency with the
# chosen population set, not merely with the input hydro N_e.
# ---------------------------------------------------------------------------
set -euo pipefail

MODEL_DIR="${1:?usage: proxy_grey_example.sh <model_dir> [output_dir]}"
OUTPUT_DIR="${2:-${MODEL_DIR}/proxy_grey_out}"

python scripts/proxy_populations_v2.py \
    "$MODEL_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --strategy GREY \
    --write-levels NF \
    --min-ion-frac 1e-30 \
    --t-exc-min 0.3 \
    --tau-trans 1.0
    # Add --write-t-ion ONLY after a round-trip test confirms the species _IN
    # T column does not alter the physical model (audit B5 / unresolved #1).

# Compare $OUTPUT_DIR/proxy_validation.txt against the Mode A run to confirm
# the grey corrections reduce the predicted first-iteration corrections without
# violating abundance closure or charge consistency.  If they do not, use the
# LTE baseline (proxy_lte_example.sh) instead.
