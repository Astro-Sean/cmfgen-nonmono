# SN_B_X_grid model study rules

Rules for adding, running, and evaluating models in this grid. Read
this plus ../MODELS_STUDY.md before submitting anything.

## Hard constraints

- **SN_A ejecta stays fully intact.** Never set `ejecta_r_inner_cm` or
  `ejecta_v_max_kms`. Wind/CSM may only occupy shell surface to the
  SN_A inner edge (~7e15 cm).
- **No 56Ni.** Keep `nickel56.enabled: false` (zero .xni). The user
  explicitly excluded radioactive nickel.
- **Intact models end in `i`** (shell scans append `_sh*` after the i;
  the plot filter accepts `endswith('i') or 'i_sh' in name`).
- Reference config: `X_V10000_w2m16i.yaml` (current best,
  score 639 / early 169 / tail 471 under the r-anchored metric).

## Workflow for a new model

1. Copy a reference YAML; change `model_name`, `description`,
   `model_dir`, and the physics keys only. Validate YAML parses.
2. Always set `cd_max_dlogrho: 60.0` and `cd_transition_zones: 5`
   when the outer CSM is dense - otherwise the CD ramp erases the
   dense segment (hardcoded |dlog rho/dlog r| < 16 otherwise).
3. Submit via `python3 scripts/run_model_grid.py <yaml...>` from the
   project root. Never run STELLA by hand.
4. Verify `SN_B_X_grid/<name>/stella/input_profile`: print
   w = 4 pi r^2 rho at 3/4.5/5.5/6.5/7.0e15 cm and confirm the
   requested column actually lands (blends and the outer_w
   normalization at 1e15 cm both rescale things).
5. Register the model in `scripts/plot_x_grid_vs_obs.py` MODEL_STYLE +
   MODEL_ORDER, then rerun the plot once the .tt file lands.

## Config gotchas

- `csm_outer_w_gcm` normalizes at 1e15 cm, not at the split radius.
  For outer_n != 2, rescale: w(split) = outer_w * (r_split/1e15)^(2-n).
- `csm_outer_n: 0` is silently swallowed (`or -1.0` in
  run_model_grid.py ~line 1749). Use n >= 0.1.
- `csm_split_blend_dex` (default 0.15) controls the tanh blend width
  of the density step; wider values spread the luminosity onset and
  reduce the +80 d bump but also lower the peak column.
- `velocity_max_kms` must be >= `shell_velocity_kms`.
- `wind_mdot_msun_yr` only applies to beta-wind mode; power-law CSM
  uses `csm_w_gcm`.
- CMFGEN temperature release: the runner pins T to hydro
  (FIX_T=T, FIX_T_AUTO/DO_T_AUTO/DO_GT_AUTO=F) and rejects
  vadat_overrides asking otherwise.  `cmfgen.solve_temperature: true`
  is the opt-out.  To solve T only in the optically thin region:
  FIX_T=F + FIX_INB_T=T + FIX_X_DPTH=<N inner depths pinned> +
  DO_JFNK=F (JFNK ignores the FIX_X_DPTH mask).  FIX_T_AUTO is NOT
  this domain -- VARFIXT pins T in the thin region and frees it
  deep.  Working example: configs/dessart_m15_2p82d_SOB_l14_freet.yaml;
  recipe is also in the prepare_VADAT.py header.

## Evaluation

- Metric: `best_phase_shift` in scripts/plot_dessart_grid_vs_obs.py.
  r-band foot gate at -19 d (crossing allowed -19..-9 d), o-band peak
  gate 0.4 mag, early (<= +80 d) weight 5.0, tail (> +80 d) weight 0.5
  on r/g only, colors weight 0. Missing coverage counts as dark.
- Check the score decomposition (early/tail) - a low total score with
  a bad early term is not a good model.
- Compare absolute r/g magnitudes +80-250 d directly, not just scores.
- The rebrightening at ~+330 d should come from SN_A's overdensity at
  1.5e16 cm; if it appears at ~+200 d the shock is too fast (needs
  ~5300 km/s average, not ~7900).

## What NOT to do again (established by tests)

- Flat high-Mdot winds (>=0.12) broaden/delay the peak.
- Wind v_inf = 100 km/s wrecks the rise (wind09li/12li early ~1800-4500).
- Homologous eruption cocoon (t_cs ~180 d) failed outright (w2m15i).
- Sharp w steps produce a luminosity step = +80 d bump.
- Dense segments that die by ~5e15 (old CD taper) cannot hold the tail.
