# Phase 2.26 level-1 reproducibility package

- v2_format_spec.md / parse_v2.py: deck format + independent parser
- deck_validation.json: block census, fbits, dtypes, per-ray status
- provenance.txt: exe/source/archive hashes, git state, compiler
- deck_sha256.txt: SHA256 of the three v2 decks (40 GB total,
  retained in the run dir - level 2)
- discrepancy_reconciliation.md: the A-E replay-fidelity table
- refresh_sequences_v2.json: all 38 per-refresh residuals/ixit/np
- topdefects_*.tsv, summary_*.json: corrected per-cell tables
- manifest_sha256.txt: hashes of the replay binaries + logs
- knop_phase226_full_archive.tar.gz (parent dir): replay binaries
  (~470 MB -> 134 MB gz), verified extraction 2026-11-21
