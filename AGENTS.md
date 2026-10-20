# SN_B_X_grid model study rules

Rules for adding, running, and evaluating models in this grid. Read
this plus ../MODELS_STUDY.md before submitting anything.

## Hard constraints

- **SN_A ejecta stays fully intact.** Never set `ejecta_r_inner_cm` or
  `ejecta_v_max_kms`. Wind/CSM may only occupy shell surface to the
  SN_A inner edge (~7e15 cm).
- **No 56Ni.** Keep `nickel56.enabled: false` (zero .xni). The user
  explicitly excluded radioactive nickel.
- **Intact models end in `i`** (variant suffixes `_sh*`, `_tcs` may
  follow the i; the plot filter accepts `endswith('i') or 'i_sh' or
  'i_tcs' in name`).
- Reference config: `X_V6000M1p5_w2m54a390i.yaml` (1.5 Msun shell @
  6000 km/s, inner w=4e17 -> 2.5e17 @4.5e15, SN_A ejecta at
  source_age 390 d).  The lighter shell gives the best g-band offset
  in the grid; the 2 Msun parent `X_V6000M2_w2m54a390i.yaml` holds
  the best total score (156.2) but runs ~0.3 mag bright in g.

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
- `ejecta_T_floor_K` is a clamp, not a target: set it BELOW the
  source ejecta's minimum T (~1100 K at ~319 d in Y_C05) or every
  zone gets flattened to the floor. Use 500.0, not 2000.0.
- `wind_mdot_msun_yr` only applies to beta-wind mode; power-law CSM
  uses `csm_w_gcm`.
- `ce_v_ramp_zones` (try 10) ramps the outer CSM velocity up to the
  ejecta inner-edge speed; without it a 20000-d `t_cs` leaves a
  ~37 vs ~2230 km/s jump at the interface and a rarefaction opens
  at t=0. `csm_age_days: 0` instead auto-sets t_cs = source age
  (319 d), making the whole CSM co-moving - velocity-continuous but
  it cuts the relative shock speed and dims the early lightcurve.
- The ejecta monotonic-velocity cleanup in
  `prepare_SN_B_model.py::clean_combined_outer_boundary` once used
  `v[i] = v[i-1] + mean(dv)`: a single source dip propagated a
  ~67 km/s/zone staircase across ~90 zones and inflated the shelf
  to ~6200 vs the true ~5500 km/s. Now clamps to `v[i-1]`. Never
  reintroduce a global-slope step there; verify u(1.59e16) ~ 5530
  km/s vs the source snapshot after every build.
- CMFGEN temperature policy (three modes, default = thin release):
  T solves in the optically thin region by default -- FIX_T=F +
  FIX_INB_T=T + FIX_X_DPTH auto-filled from the tau=2/3 photosphere.
  `cmfgen.fix_temperature: true` restores the old all-depths hydro
  pin (required for A/B comparisons against older pinned runs);
  `cmfgen.solve_temperature: true` means fully free T (no stamps).
  For thin runs also set DO_JFNK=F (JFNK ignores the FIX_X_DPTH
  mask).  FIX_T_AUTO is NOT this domain -- VARFIXT pins T in the
  thin region and frees it deep.  Working example:
  configs/dessart_m15_2p82d_SOB_l14_freet.yaml; recipe is also in
  the prepare_VADAT.py header.

## Evaluation

- Metric: `best_phase_shift` in scripts/plot_dessart_grid_vs_obs.py.
  r-band foot gate at -19 d (crossing allowed -19..-9 d), o-band peak
  gate 0.4 mag, early (<= +80 d) weight 5.0, tail (> +80 d) weight 0.5
  on r/g only, colors weight 0. Missing coverage counts as dark.
- Check the score decomposition (early/tail) - a low total score with
  a bad early term is not a good model.
- Compare absolute r/g magnitudes +80-250 d directly, not just scores.
- The +330 d feature in the built models is the shock crossing the
  CSM/ejecta interface and entering SN_A's inner density shelf
  (w ~ 5e16 -> 3e17 at 6.5e15-1.6e16), NOT the overdensity peak at
  1.5e16 - that recedes at ~5500-6000 km/s and outruns the swept-up
  shell (~1900-5500 km/s). To delay the feature toward +330 d,
  raise the outer CSM column (`csm_outer_w_gcm` 2.5e17 -> 4e17/6e17
  adds ~0.1/0.25 Msun) rather than slowing the shell, which wrecks
  the early light. A ~250-270 d dip is the shock leaving the dense
  CSM; more outer column fills it inside the 245-325 d data gap.

## What NOT to do again (established by tests)

- Flat high-Mdot winds (>=0.12) broaden/delay the peak.
- Wind v_inf = 100 km/s wrecks the rise (wind09li/12li early ~1800-4500).
- Homologous eruption cocoon (t_cs ~180 d) failed outright (w2m15i).
- Sharp w steps produce a luminosity step = +80 d bump.
- Dense segments that die by ~5e15 (old CD taper) cannot hold the tail.
