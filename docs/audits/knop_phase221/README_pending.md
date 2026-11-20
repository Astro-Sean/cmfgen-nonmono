# knop_phase221 artifacts - pending items

Live analyses awaiting job completion:

- A1 13414391 / A2 13414398: pre-flag replicate runs (exe 93364cac),
  started 2026-11-20 10:52, RUNNING at snapshot time. Extract at
  completion: great iterations, release/freeze counts, census history,
  GS_MAX_ALL, GS_RF_DCHG, MAXCH, termination.
- B1 13427111 / B2 13427115: controlled pre/post pair (93364cac vs
  58a699d0), identical VADAT_INIT/SN_HYDRO_DATA/donor checkpoint.
- OBS chain: A-arms 13414392/13414399, B-arms 13427112/13427116,
  standalone p221 arms 13427285 (diag on) / 13427286 (diag off).
- Byte-compare OBSFLUX/obs_fin_10/OUT_FLUX between 13427285/13427286
  after completion.

## Update 20-Nov-2026 - MAT_FREEZE finding + diagcap arm

A1/A2/B1/B2 all carry KN_MAT_FREEZE=T (inherited from the iso_M
donor): populations are held by design. The 50% convergence plateau
and held=T transcript are NOT a convergence stall. Field-level
finding: refresh residual plateaus ~6.5E-2 over 122 refreshes;
offender set {22,45,50} + probe 83 stable since refresh 2; 2 capped
commits/refresh. diagcap arm 13429272 (KN_PSWEEP all-rays) submitted
to identify the capped rays. See diagcap_arm.txt and the audit Part B.
