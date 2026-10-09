#!/usr/bin/env python3
"""Generate the SN_B "Dessart 2010jl analog" grid.

Dessart, Audit & Hillier (2015, MNRAS 449, 4304; docs/1503.05463v1.txt)
model X: a 9.8 Msun / 1 Bethe ejecta shell (flat core below V0=Vmax/3,
rho ~ v^-8 envelope; Vmax=9608 km/s) colliding with a steady
100 km/s wind of Mdot=0.1 Msun/yr that is truncated at 1e16 cm
(Mdot=1e-3 beyond).  That light curve matches the 2024pba morphology
but is ~1 mag too bright and decays too slowly.

This grid keeps the Dessart shell/wind construction (v1
prepare_SN_B_model.py path: shell + steady 1/r^2 wind + SN_A ejecta
appended at its homologous position) and scans the CSM loading
(w = Mdot/v_w) around the 2024pba-appropriate value
(~w=1.5e17 g/cm -> Mdot ~ 0.02-0.03 at v_w=100 km/s), plus the
cavity width (t_start sets the shell/ejecta contact radius).

Parent: SN_B_double_peak_grid/SN_B_T6k_M7.0_V10000_ts15/config.yaml
(SN_A_g_T20k_Ni06 ejecta at 400 d, Ni56 off, STELLA only, 100 d).

Grid points
-----------
DX_M7_ts15_W1e-2            Mdot=1e-2   thin-Dessart bridge point
DX_M7_ts15_W3e-2            Mdot=3e-2   2024pba luminosity scale
DX_M7_ts15_W1e-1            Mdot=1e-1   Dessart X density
DX_M7_ts8_W3e-2             t_start=8   wider cavity (~3x the wind mass)
DX_M9p8_V9608_ts11p6_W1e-1  exact model-X shell (9.8 Msun, 9608 km/s,
                            11.6 d) + SN_A ejecta appended
DX_M7_ts15_Wdm100           density-matched wind (Mdot derived from
                            ejecta inner edge -> continuous join)
DX_M7_ts8_W3e-2_dtr         truncated wind: dense to r_lead~3e15 cm
                            then Mdot=1e-3 background to the ejecta
                            (Dessart's two-region decline driver)

DX2_* block (t_start=2 d, SN_A ejecta fixed at 330 d, run to 170 d
observer time): Mdot ladder 1e-2..1e-1 testing whether a heavier
inner-wind column sustains the post-peak decline (linear-in-mag)
through continued shell-wind interaction rather than ejecta
envelope loading.  DX2_M9p8 adds the model-X shell mass.

Each config is the parent with only the keys listed above modified.
Existing config.yaml files are not overwritten; re-running only writes
missing grid points.
"""
import math
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
from orchestration.config import (
    OBSERVED_SN_A_AGE_DAYS,
    SN_A_LAUNCH_MJD,
    SN_B_DESSART_TIMING_FAMILY_START,
    SN_B_T0_MJD,
)

BASE_DIR = Path("/viper/ptmp1/sbrennan/SN2024pba")
GRID_DIR = BASE_DIR / "SN_B_dessart_grid"
PARENT_CONFIG = (
    BASE_DIR / "SN_B_double_peak_grid" / "SN_B_T6k_M7.0_V10000_ts15"
    / "config.yaml"
)

# (name, shell_mass, shell_v, t_start_d, wind_mdot, extras dict)
GRID_POINTS = [
    ("DX_M7_ts15_W1e-2", 7.0, 10000, 15.0, 1e-2, {}),
    ("DX_M7_ts15_W3e-2", 7.0, 10000, 15.0, 3e-2, {}),
    ("DX_M7_ts15_W1e-1", 7.0, 10000, 15.0, 1e-1, {}),
    ("DX_M7_ts8_W3e-2", 7.0, 10000, 8.0, 3e-2, {}),
    ("DX_M9p8_V9608_ts11p6_W1e-1", 9.8, 9608, 11.6, 1e-1, {}),
    ("DX_M7_ts15_Wdm100", 7.0, 10000, 15.0, 5e-2,
     {"density_matched_wind": True, "allow_steep_interface": True}),
    ("DX_M7_ts8_W3e-2_dtr", 7.0, 10000, 8.0, 3e-2,
     {"wind_duration_days": 3472.0, "wind_mdot_bg_msun_yr": 1e-3}),
    # --- DX2 decay-phase block: t_start=2 d, SN_A fixed at 330 d ---
    # SN_A ejecta always at source_age=330 d (fixed two-event
    # separation).  The sustained post-peak shock power should come
    # from the inner shell<->wind interaction, so the scanned axis is
    # the wind loading Mdot (column in the shell->ejecta cavity), not
    # the ejecta density.  target_days=170 so the decay is visible.
    ("DX2_M7_ts2_W1e-2_sa330", 7.0, 10000, 2.0, 1e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX2_M7_ts2_W3e-2_sa330", 7.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX2_M7_ts2_W6e-2_sa330", 7.0, 10000, 2.0, 6e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX2_M7_ts2_W1e-1_sa330", 7.0, 10000, 2.0, 1e-1,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX2_M9p8_ts2_W3e-2_sa330", 9.8, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    # --- DX3 smeared-decline block: shallow power-law CSM ---
    # v2 builder (powerlaw_csm=True): rho ~ r^-n homologous CSM,
    # density-matched to the ejecta inner edge, SN_A still at 330 d.
    # For n<2 the swept column 4*pi*r^2*rho ~ r^(2-n) grows outward,
    # so the shock keeps finding optically thick material and the
    # peak blends into the tail.  wind_mdot is ignored by v2
    # (parsed placeholder only).
    ("DX3_M7_ts2_pl10_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0}),
    ("DX3_M7_ts2_pl15_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.5}),
    ("DX3_M7_ts2_pl18_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.8}),
    # Explicit normalisation: w = 4 pi r^2 rho at 1e15 cm fixed at the
    # 2024pba luminosity scale (w ~ 1e17 -> L ~ 0.5 w v_s^3 ~ few e42
    # at v_s ~ 5e3 km/s).  csm_density_match=False keeps rho15 fixed;
    # the ejecta interface keeps a density step absorbed by the CD ramp.
    ("DX3_M7_ts2_pl15_w1e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.5,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "cd_transition_zones": 15}),
    ("DX3_M7_ts2_pl10_w1e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "cd_transition_zones": 15}),
    ("DX3_M7_ts2_pl18_w1e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.8,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "cd_transition_zones": 15}),
    # --- DX4 shell-mass ladder: fixed V=10000 (KE ~ M) at the best
    # wind column; M3_V15000 is the fixed-KE control (0.5*3*1500^2 ~
    # 6.8e50 ~ same as M7@10000).  Lighter shells decelerate faster in
    # the same CSM, so the rise/decline clock shortens.
    ("DX4_M5_ts2_W3e-2_sa330", 5.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX4_M4_ts2_W3e-2_sa330", 4.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX4_M3_ts2_W3e-2_sa330", 3.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX4_M2_ts2_W3e-2_sa330", 2.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    ("DX4_M3_V15000_ts2_W3e-2_sa330", 3.0, 15000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 170.0}),
    # --- DX5 delayed-onset block: thin inner wind + dense outer wind ---
    # wind_mdot fills shell edge -> r_lead; wind_mdot_bg fills r_lead ->
    # ejecta.  bg > mdot inverts the usual ordering: a thin recent wind
    # (faint early phase, possibly powering the -15.5 plateau) then a
    # dense older outflow the shock runs into ~15-20 d before peak.
    # r_lead = v_wind * wind_duration_days.
    # NB: model names must stay <=29 chars - STELLA reads the strad.1
    # name line into a fixed 160-char buffer (5*N+14 <= 160).
    ("DX5_M7_i1e3_o3e2_rl2_sa330", 7.0, 10000, 2.0, 1e-3,
     {"source_age_days": 330.0, "target_days": 170.0,
      "wind_duration_days": 2315.0, "wind_mdot_bg_msun_yr": 3e-2}),
    ("DX5_M7_i1e3_o6e2_rl2_sa330", 7.0, 10000, 2.0, 1e-3,
     {"source_age_days": 330.0, "target_days": 170.0,
      "wind_duration_days": 2315.0, "wind_mdot_bg_msun_yr": 6e-2}),
    ("DX5_M7_i3e4_o3e2_rl1_sa330", 7.0, 10000, 2.0, 3e-4,
     {"source_age_days": 330.0, "target_days": 170.0,
      "wind_duration_days": 1157.0, "wind_mdot_bg_msun_yr": 3e-2}),
    # --- DX6 outer-CSM tail block ---
    # A tenuous wind OUTSIDE the SN_A ejecta edge: the ejecta's leading
    # edge (~1e9 cm/s) shocks it continuously, giving a nearly constant
    # L ~ 0.5*w*v^3 ~ few e41 erg/s - the observed -16 tail floor that
    # persists to +245d.  Physically SN_A's own pre-explosion wind.
    # Outer edge at 8e16 so the ejecta front never exits the wind in
    # the 300 d run.  Base model = DX4_M5 (best rise/peak score).
    # r_outer=1.2e17 (0.6 dex span): the inner-edge density ramp needs
    # ~dex/(10*dlogr) zones, so a wider CSM span is required for the
    # thinner winds before n_ramp hits the n_zones//2 cap.
    ("DX6_M5_o1e4_sa330", 5.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 300.0,
      "outer_csm_mdot_msun_yr": 1e-4, "outer_csm_velocity_kms": 850.0,
      "outer_csm_r_outer_cm": 1.2e17}),
    ("DX6_M5_o3e4_sa330", 5.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 300.0,
      "outer_csm_mdot_msun_yr": 3e-4, "outer_csm_velocity_kms": 850.0,
      "outer_csm_r_outer_cm": 1.2e17}),
    ("DX6_M5_o1e3_sa330", 5.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 300.0,
      "outer_csm_mdot_msun_yr": 1e-3, "outer_csm_velocity_kms": 850.0,
      "outer_csm_r_outer_cm": 1.2e17}),
    ("DX6_M5_o3e3_sa330", 5.0, 10000, 2.0, 3e-2,
     {"source_age_days": 330.0, "target_days": 300.0,
      "outer_csm_mdot_msun_yr": 3e-3, "outer_csm_velocity_kms": 850.0,
      "outer_csm_r_outer_cm": 1.2e17}),
    # --- DX7 peak-width/tail-level block ---
    # DX3 pl10 w1e17 wins the tail-inclusive score: n=1.0 gives the
    # right ~1 mag/100d decline but the peak is a bit narrow and the
    # tail sits ~0.5-1 mag faint.  L ~ 0.5 w v_s^3 scales with the
    # column, so a w ladder at fixed n=1.0 brightens the tail and
    # lengthens the shock-crossing time (broader peak) without
    # changing the decline slope.  pl09 brackets the rate: n=0.9 grows
    # the column slightly faster at large r.
    ("DX7_M7_ts2_pl10_w15e16_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX7_M7_ts2_pl10_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX7_M7_ts2_pl10_w3e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX7_M7_ts2_pl09_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # --- w ladder extension ---
    # The g-r colour plot shows the densest column (w3e17) is the
    # reddest model and closest to the observed +0.2..+0.3 tail colour
    # - denser CSM reprocesses more of the blue continuum.  Extend the
    # w ladder upward; if these overshoot the tail luminosity the
    # colour-versus-flux tension is itself diagnostic.
    ("DX7_M7_ts2_pl09_w4e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX7_M7_ts2_pl10_w4e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX7_M7_ts2_pl10_w5e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 5e17,
      "cd_transition_zones": 15}),
    # --- DX9 density-law scan ---
    # DX7 winner (pl09/pl10, w=2e17) with only the CSM slope varied:
    # rho ~ r^-n for n = 3, 2, 0.  n=1 is covered by the identical
    # DX7_M7_ts2_pl10_w2e17_sa330 run, n=0.9 by the pl09 run.
    # n=3 drops the swept column outward (shock runs out of fuel),
    # n=0 keeps it growing as r^2 (thick outer reservoir).
    ("DX9_M7_ts2_pl30_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 3.0,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX9_M7_ts2_pl20_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 2.0,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # n=0 (constant rho) ends at rho15=1.6e-14 vs ejecta 8e-17 - a
    # 2.3 dex step that needs ~40 CD ramp zones to keep the per-zone
    # density gradient under the |dlog rho/dlog r|=16 solver limit.
    ("DX9_M7_ts2_pl00_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.0,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 40}),
    # Intermediate slopes: the DX7 n=0.9 winner declines slightly too
    # steep, n=0 probably too flat - bracket with 0.5/0.75.
    ("DX9_M7_ts2_pl075_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX9_M7_ts2_pl05_w2e17_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # --- DX10 compact inner CEE block ---
    # Goal: brighten the early (pre-peak) bump to match the observed
    # r ~ -15.5 plateau at -45..-19 d.  The DX7 bump is breakout
    # cooling: ~5e40 erg/s bolometric but only R ~ -9.4 - hot and
    # compact, so its flux is in the UV.  A dense, slow shell (a ~3
    # week old pre-collapse ejection: cee_age=20 d -> u ~ 1000..3500
    # km/s) is shocked hard within the first days; its optical depth
    # releases the stored shock energy over ~weeks at lower T_eff,
    # which is what the plateau needs.  Outer law stays pl09 w=2e17
    # (the DX7 winner) so only the early phase changes.
    # Grid constraint: the shell-CSM buffer (sw_transition_dex=0.5)
    # overwrites rho out to ~5.5e14 cm, so the inner CEE law is only
    # sampled if r_split > ~5.5e14 - cee_r_outer=6e14 puts a flat
    # dense shell across 5.5-6e14 with the buffer ramping up into it.
    # Scan: CEE mass 0.3 / 0.6 / 1.2 Msun (rho ~ 1.1/2.3/4.5e-12).
    ("DX10_M7_ce6_m03_pl09_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cee_age_days": 20.0, "cee_r_outer_cm": 6e14,
      "cee_mass_msun": 0.3, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX10_M7_ce6_m06_pl09_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cee_age_days": 20.0, "cee_r_outer_cm": 6e14,
      "cee_mass_msun": 0.6, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX10_M7_ce6_m12_pl09_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cee_age_days": 20.0, "cee_r_outer_cm": 6e14,
      "cee_mass_msun": 1.2, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    # --- DX8 two-region CSM block ---
    # Structure: shell -> compact CEE shell (60 d old at SN_B launch,
    # i.e. a terminal pre-collapse ejection) -> fast wind bridging to
    # the SN_A ejecta -> SN_A ejecta.  The CEE is hit almost
    # immediately (compact, dense), powering the plateau/rise and
    # broadening the peak; the outer wind column carries the tail.
    # The wind must be fast (~2e3 km/s): anything sitting between a
    # 60 d-old shell and the 330 d-old ejecta must outrun SN_A's
    # slowest ejecta (~2350 km/s inner edge) to avoid engulfment.
    # v_inf=2000, beta=1, R*=2e14 gives v(1e15)=1600 km/s ~ the CEE
    # edge speed (u=1929 km/s for t_cee=60 d), so the junction is
    # near velocity-continuous.  Mdot 0.15/0.3 = w~5e16/1e17 at
    # v~1.7e3 km/s, bracketing the DX3 winner's column.
    # The _pl10 control keeps the winning n=1 envelope outside the
    # CEE: isolates what the inner CEE adds without changing the tail
    # law (velocity drops 1929->349 km/s at the split, blended).
    # Names <=29 chars for the STELLA strad.1 name buffer.
    ("DX8_M7_cee60_m05_wd15_sa330", 7.0, 10000, 2.0, 0.15,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_density_match": False,
      "beta_wind": True, "wind_velocity_kms": 2000.0,
      "wind_beta": 1.0, "wind_rstar_cm": 2e14,
      "cee_age_days": 60.0, "cee_r_outer_cm": 1e15,
      "cee_mass_msun": 0.5, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX8_M7_cee60_m05_wd30_sa330", 7.0, 10000, 2.0, 0.3,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_density_match": False,
      "beta_wind": True, "wind_velocity_kms": 2000.0,
      "wind_beta": 1.0, "wind_rstar_cm": 2e14,
      "cee_age_days": 60.0, "cee_r_outer_cm": 1e15,
      "cee_mass_msun": 0.5, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX8_M7_cee60_m10_wd30_sa330", 7.0, 10000, 2.0, 0.3,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_density_match": False,
      "beta_wind": True, "wind_velocity_kms": 2000.0,
      "wind_beta": 1.0, "wind_rstar_cm": 2e14,
      "cee_age_days": 60.0, "cee_r_outer_cm": 1e15,
      "cee_mass_msun": 1.0, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX8_M7_cee60_m05_pl10_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "cee_age_days": 60.0, "cee_r_outer_cm": 1e15,
      "cee_mass_msun": 0.5, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    # --- DX11 shell mass/velocity colour scan ---
    # All DX7/DX9/DX10 models hold the shell at M7, V=10000, yet every
    # model sits ~0.2-0.3 mag too blue in g-r while w3e17 is reddest.
    # The shell controls two colour levers the CSM column does not:
    # a denser shell buries the photosphere deeper (larger radius,
    # cooler, more reprocessing -> redder), and shell velocity sets the
    # shock speed at the ejecta contact (hotter shock -> bluer).
    # Outer law fixed at the winner (pl09, w=2e17, sa330) so any g-r
    # shift is attributable to the shell alone.  Watch for degeneracy:
    # M and V both scale the shock power (KE ~ M v^2), so a redder
    # model that overshoots the peak is not a win.
    ("DX11_M4_ts2_pl09_sa330", 4.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX11_M10_ts2_pl09_sa330", 10.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX11_M12_ts2_pl09_sa330", 12.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX11_M7_V8000_ts2_pl09_sa330", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX11_M7_V14000_ts2_pl09_sa330", 7.0, 14000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # --- DX12 source-age scan ---
    # DX7 winner (M7, V10000, pl09, w=2e17) with only source_age_days
    # varied 300..400 in 25 d steps around sa330.  sa sets the SN_A
    # ejecta expansion epoch AND t_cs (CSM homologous time): at fixed
    # w(1e15) the same column is stretched over a different radius
    # range (u = r/t_cs), so the density near the shell, the shock
    # crossing time and the interface radius all shift together.
    ("DX12_M7_ts2_pl09_w2e17_sa300", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 300.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX12_M7_ts2_pl09_w2e17_sa325", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 325.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX12_M7_ts2_pl09_w2e17_sa350", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX12_M7_ts2_pl09_w2e17_sa375", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX12_M7_ts2_pl09_w2e17_sa400", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # --- DX13 hybrid: compact inner ejection + pl075 envelope ---
    # Combines the two scoring winners in one shell -> envelope ->
    # SN_A system: DX10's dense ~20-day-old inner ejection (ce6,
    # r_outer=6e14, which holds the early light and posts the best
    # tail chi2) with DX9's pl075 outer law (n=0.75, the slope the
    # decline actually wants).  Kept inside the wind/envelope
    # picture: the ce6 shell is a terminal pre-collapse outflow, not
    # a separate CSM component.
    ("DX13_M7_ce6_m03_pl075_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cee_age_days": 20.0, "cee_r_outer_cm": 6e14,
      "cee_mass_msun": 0.3, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    ("DX13_M7_ce6_m06_pl075_sa330", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cee_age_days": 20.0, "cee_r_outer_cm": 6e14,
      "cee_mass_msun": 0.6, "cee_powerlaw_n": 0.0,
      "cee_blend_dex": 0.2, "cd_transition_zones": 15}),
    # --- DX14 combined-winner optimisation round ---
    # All scans at M7, ts2, pl075 (n=0.75 is the new slope winner) and
    # V=8000 km/s (DX11: best rise chi2 of the grid, slower shock =
    # cooler early light).  Varied: source age 330/350/375 (DX12 put
    # the optimum at 350-375) and column w = 2/2.5/3e17 (w3e17 was the
    # reddest lightcurve-family model; the slower shell may afford
    # more column before overshooting).  Names imply M7, ts2, pl075.
    ("DX14_V8000_sa330", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX14_V8000_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX14_V8000_sa375", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX14_V8000_w25e17_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2.5e17,
      "cd_transition_zones": 15}),
    ("DX14_V8000_w3e17_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX14_V8000_w3e17_sa375", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    # --- DX15 physical reddening round ---
    # All models are ~0.2 mag too blue in g-r at fixed MW-only
    # extinction (E(B-V)=0.15, host excluded by choice).  Push the
    # levers the scans identified without the luminosity overshoot
    # that killed w4e17 at V10000: even slower shell (V7000, cooler
    # shock), heavier shell (M10, deeper photosphere), thick column
    # on the slow shell, and the pl05 envelope (best g-r among the
    # lightcurve-matching models).  pl075/sa350/w2e17 baseline.
    ("DX15_V7000_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX15_V7000_w3e17_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX15_M10_V8000_sa350", 10.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX15_M10_V8000_w3e17_sa350", 10.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX15_V8000_pl05_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX15_V8000_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    # --- DX16 tail+colour corner: sa400 x flat n x thick column ---
    # Residual analysis of the finished grid: all models are too
    # bright at 0-60 d and too faint at 90-130 d, and g-r sits ~0.15
    # mag too blue nearly constant in phase.  sa375 halves the late
    # fade vs sa330 - the shock exhausts the CSM at the SN_A
    # interface (6.7e15 cm at sa330, reached ~+80 obs d).  Colour
    # residuals shrink with outer column (pl05 -0.10, w3e17 -0.12
    # vs pl075 -0.17).  Combine: sa400 pushes the interface to
    # ~8e15 cm, n=0.5-0.75 brackets the tail slope, w=3e17 adds the
    # column the colour wants - all on the V8000 shell (best rise),
    # plus one V10000 control to isolate the velocity contribution.
    ("DX16_V8000_pl05_w3e17_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 80}),
    ("DX16_V8000_pl06_w3e17_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.6,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 80}),
    ("DX16_V8000_pl075_w3e17_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 80}),
    ("DX16_pl05_w3e17_sa400", 7.0, 10000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "cd_transition_zones": 80}),
    # --- DX17 dim-the-colour-matchers round ---
    # The w4e17/w5e17 columns reproduce the observed colour decay
    # (chi2_col 24-31) but at V10000 they run ~0.5-1 mag too bright.
    # L ~ w v^3, so V8000 cuts ~half the luminosity while keeping
    # the column (and hence the reprocessing-driven red colour):
    # w4e17@V8000 ~ w2e17@V10000 in power with ~2x the column.
    # V7000 brackets the dimmer side.
    ("DX17_V8000_pl09_w4e17_sa330", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX17_V8000_pl10_w4e17_sa330", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 1.0,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX17_V7000_pl09_w4e17_sa330", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.9,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX17_V8000_pl05_w4e17_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 4e17,
      "cd_transition_zones": 80}),
    # --- DX18 outer-reservoir round ---
    # Residual structure of DX14_V8000_sa330 (best overall): ~0 mag
    # through +65..+90 d then +0.3..+0.5 mag faint at +95..+140 d -
    # the shock reaches the tenuous SN_A ejecta.  sa375 sustains the
    # tail but is ~0.25 mag too bright at all earlier phases because
    # the whole column grows.  The new csm_r_split_cm option decouples
    # the two: inner w=1.5e17 (thinner -> dimmer mid) then the density
    # jumps UP to a thicker outer segment at r_split=4.5e15 cm (~65 d
    # at 8000 km/s, where the sa330 residual crosses zero) - a second,
    # later eruptive shell.  sa330/375/400 brackets move the SN_A
    # interface while the reservoir supplies the late luminosity.
    ("DX18_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15}),
    ("DX18_sa330_i15_w30", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX18_sa375_i15_w30", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX18_sa375_i15_w35", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3.5e17,
      "cd_transition_zones": 15}),
    ("DX18_sa400_i15_w30", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX18_sa375_i15_w30_n05", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3e17,
      "csm_outer_n": 0.5,
      "cd_transition_zones": 15}),
    # --- DX19 CSM-velocity (t_cs) round ---
    # The CSM velocity u = r/t_cs is auto-pinned to the SN_A ejecta
    # age for interface continuity, but it is a physical lever of its
    # own: L ~ w*v_rel^3 and T_shock ~ v_rel^2 with v_rel = v_shell -
    # u_csm.  A FASTER CSM (smaller t_cs) dims AND reddens at fixed
    # column - the direction the w3e17 colour matchers need (they hold
    # chi2_col ~ 26 but sit ~0.4 mag too bright in the mid-phase; the
    # required dimming is v_rel x 0.88, i.e. u(4.5e15) 1400 -> ~2200
    # km/s, t_cs ~ 230 d).  A t_cs=150/200/250 scan brackets it.  The
    # reverse test on sa330 (t_cs=450, slower CSM -> brighter)
    # probes whether extra interaction power fills its +95..140 d
    # tail deficit without growing the column.
    ("DX19_w3e17_sa350_tcs250", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "csm_age_days": 250.0,
      "cd_transition_zones": 15}),
    ("DX19_w3e17_sa350_tcs200", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "csm_age_days": 200.0,
      "cd_transition_zones": 15}),
    ("DX19_w3e17_sa350_tcs150", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 3e17,
      "csm_age_days": 150.0,
      "cd_transition_zones": 15}),
    ("DX19_sa330_tcs450", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 330.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "csm_age_days": 450.0,
      "cd_transition_zones": 15}),
    # --- DX20 dim-the-colour-matchers via split CSM ---
    # The best colour matchers (DX14_w3e17_sa350, col chi2 26) sit
    # ~0.4-0.7 mag too bright because the SAME thick column that makes
    # the red colour is swept during the peak.  Decouple them with
    # csm_r_split_cm = 3e15 cm (~43 d at V8000, where the overshoot is
    # worst): a thin inner segment dims the peak/mid luminosity, then
    # the density jumps up to the colour-producing column so tau_es
    # and the late reprocessing are preserved.
    ("DX20_i20_w30_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "csm_r_split_cm": 3e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX20_i15_w30_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 3e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    ("DX20_i15_w40_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "csm_r_split_cm": 3e15, "csm_outer_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    ("DX20_i10_w40_sa350", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "csm_r_split_cm": 3e15, "csm_outer_w_gcm": 4e17,
      "cd_transition_zones": 15}),
    # --- DX21 smooth-CSM round (no piecewise split) ---
    # Single smooth rho ~ r^-0.75 envelope throughout, per the
    # structure requirement.  Within smooth laws the late-time deficit
    # (+0.3..+0.5 mag at +95..140 d) is the shock reaching the SN_A
    # interface: sa330 crosses at ~+90 d, sa375 at ~+105 d (its tail
    # residual is only +0.14 at +125).  The sa375 penalty is a uniform
    # ~0.25 mag mid-phase overshoot from the larger column - so push
    # sa to 375-450 AND drop w to 1.5-1.75e17, which should cancel the
    # mid excess while keeping the envelope extended to the end of the
    # observed window.
    ("DX21_sa375_w15e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 375.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX21_sa400_w15e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX21_sa400_w175e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "cd_transition_zones": 15}),
    ("DX21_sa450_w15e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 450.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX21_sa450_w10e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 450.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1e17,
      "cd_transition_zones": 15}),
    # --- DX22 smooth-CSM round 2 (no piecewise split) ---
    # Residual map of the best smooth models: DX15_V7000_sa350 is
    # ~0.15 mag bright mid-phase then +0.28 faint at +125 (interface
    # reached ~+100); DX21_sa400_w15e17 is flat to +110 then +0.22 at
    # +125; DX21_sa450_w10e17 overshot faint overall.  Optimum: push
    # sa to 400-450 on the V7000 shell where the slower shock reaches
    # the interface later anyway, and fill the w gap between 1.0e17
    # (too faint) and 1.5e17 on V8000 sa425-450.
    ("DX22_V7000_sa400", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2.0e17,
      "cd_transition_zones": 15}),
    ("DX22_V7000_sa400_w175e17", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "cd_transition_zones": 15}),
    ("DX22_V7000_sa450_w15e17", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 450.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX22_V6500_sa400", 7.0, 6500, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2.0e17,
      "cd_transition_zones": 15}),
    ("DX22_sa425_w15e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 425.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX22_sa450_w125e17", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 450.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.25e17,
      "cd_transition_zones": 15}),
    # --- DX23 decay-phase round (submitted 2026-10-04) ---
    # Per-epoch residuals of DX15_V7000_sa350 at its best shift: ~0.1-
    # 0.35 mag too bright from the peak through +75 d, on the data at
    # +75-90 d, then drifting faint to +0.3 at +125.  The initial
    # interpretation was that the shock reached the SN_A ejecta edge near
    # +95.  Later STELLA profiles showed the compressed front still inside
    # the gap CSM through about +125.  This round tests whether a shallower
    # or denser outer CSM sustains power, and thins the inner column to
    # limit the mid-phase excess, on the V7000 backbone.
    ("DX23_V7000_w175e17_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "cd_transition_zones": 15}),
    ("DX23_V7000_i175_w30_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "csm_r_split_cm": 4.5e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    # Reservoir placement followed the initial estimate of where the
    # residual began. DX15 profile inspection later showed the compressed
    # front still inside the gap CSM through about obs +125.
    ("DX23_V7000_i175_w30_rs60_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "csm_r_split_cm": 6e15, "csm_outer_w_gcm": 3e17,
      "cd_transition_zones": 15}),
    # Sensitivity test: ejecta_r_inner_cm=9e15 removes the tenuous inner
    # SN_A ejecta zones (rho ~ 5e-17 g/cc, <<0.1 Msun) so the smooth
    # pl075 law extends farther out. This modifies the physical ejecta
    # profile and is not a baseline model. Later DX15 profile inspection
    # showed the compressed front still inside the gap CSM through about
    # obs +125, so an ejecta-edge crossing is not established.
    ("DX23_V7000_w175e17_rei9e15", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "ejecta_r_inner_cm": 9e15,
      "cd_transition_zones": 15}),
    # Same extension on the best V8000 model (DX21_sa400_w15e17) as an
    # A/B test of shell velocity degeneracy; sa400 also pushes the
    # physical interface out so the remaining gap to 9e15 is smaller.
    ("DX23_V8000_w15e17_rei9e15_sa400", 7.0, 8000, 2.0, 0.0,
     {"source_age_days": 400.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "ejecta_r_inner_cm": 9e15,
      "cd_transition_zones": 15}),
    # Radioactive-tail falsification test: the observed +90..130 d
    # flattening at Mr ~ -16.3 matches 56Co decay from ~0.03 Msun of
    # Ni56 in the SN_B shell (L ~ 1.3-1.5e41 erg/s over the window).
    # Same DX15 backbone with Ni56 enabled; if the tail is radioactive
    # this lands on it directly rather than via extra CSM column.
    ("DX23_V7000_sa350_Ni003", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 170.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.75,
      "csm_density_match": False, "csm_w_gcm": 2e17,
      "cd_transition_zones": 15,
      "stella_nickel56": {"enabled": True, "total_mass_msun": 0.03,
                          "distribution_mass_msun": 1.0,
                          "distribution_type": "uniform",
                          "inner_zone_fraction": 0.5}}),
    ("DX24_V7000_n05_w15_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX24_V7000_n05_w175_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "cd_transition_zones": 15}),
    ("DX24_V7000_n065_w15_sa350", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 350.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.65,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX25_V7000_n05_w15_sa334", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 334.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
    ("DX25_V7000_n05_w175_sa334", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 334.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.5,
      "csm_density_match": False, "csm_w_gcm": 1.75e17,
      "cd_transition_zones": 15}),
    ("DX25_V7000_n065_w15_sa334", 7.0, 7000, 2.0, 0.0,
     {"source_age_days": 334.0, "target_days": 200.0,
      "powerlaw_csm": True, "csm_powerlaw_n": 0.65,
      "csm_density_match": False, "csm_w_gcm": 1.5e17,
      "cd_transition_zones": 15}),
]

TARGET_DAYS = 100.0


def make_observational_timing(model_name, source_age, extras):
    policy = extras.pop("source_age_policy", None)
    reason = extras.pop("source_age_reason", None)
    family = re.match(r"DX(\d+)_", model_name)
    family_number = int(family.group(1)) if family else None
    is_observed_age = math.isclose(
        source_age, OBSERVED_SN_A_AGE_DAYS, abs_tol=1e-6
    )

    if policy is None:
        if is_observed_age:
            policy = "observed"
        elif family_number is not None and family_number >= SN_B_DESSART_TIMING_FAMILY_START:
            raise ValueError(
                f"{model_name} uses source_age_days={source_age:g}, not the "
                f"observed {OBSERVED_SN_A_AGE_DAYS:g} d. Add an explicit "
                "source_age_policy='sensitivity' and source_age_reason only "
                "for an age-sensitivity model."
            )
        else:
            policy = "sensitivity"
            reason = reason or "Legacy source-age grid point; not observation-anchored."

    if policy == "observed" and not is_observed_age:
        raise ValueError(
            f"{model_name} marks its source age observed, but "
            f"source_age_days={source_age:g}; expected "
            f"{OBSERVED_SN_A_AGE_DAYS:g}."
        )
    if policy == "sensitivity" and not str(reason or "").strip():
        raise ValueError(
            f"{model_name} source-age sensitivity requires source_age_reason."
        )
    if policy not in ("observed", "sensitivity"):
        raise ValueError(
            f"{model_name} source_age_policy must be observed or sensitivity."
        )

    timing = {
        "sn_a_launch_mjd": SN_A_LAUNCH_MJD,
        "sn_b_t0_mjd": SN_B_T0_MJD,
        "observed_source_age_days": OBSERVED_SN_A_AGE_DAYS,
        "configured_source_age_days": source_age,
        "source_age_policy": policy,
    }
    if policy == "sensitivity":
        timing["sensitivity_reason"] = reason
    return timing


def main():
    with open(PARENT_CONFIG) as f:
        parent = yaml.safe_load(f)

    GRID_DIR.mkdir(parents=True, exist_ok=True)
    configs = []
    for name, mass, vel, ts, mdot, extras in GRID_POINTS:
        model_dir = GRID_DIR / name
        config_path = model_dir / "config.yaml"
        if config_path.exists():
            print(f"  Exists, skipping: {config_path}")
            continue
        model_dir.mkdir(parents=True, exist_ok=True)

        extras = dict(extras)
        target_days = extras.pop("target_days", TARGET_DAYS)
        nickel56 = extras.pop("stella_nickel56", None)
        source_age = float(
            extras.pop("source_age_days", OBSERVED_SN_A_AGE_DAYS)
        )
        profile_age = source_age + ts
        timing = make_observational_timing(name, source_age, extras)
        timing["profile_age_days"] = profile_age

        cfg = dict(parent)
        cfg["model_name"] = name
        if extras.get("powerlaw_csm", False):
            csm_setup = (
                f"homologous CSM rho~r^-{extras['csm_powerlaw_n']:g}, "
                f"w={extras['csm_w_gcm']:.3g} g/cm"
            )
        else:
            csm_setup = f"steady wind Mdot={mdot} Msun/yr at v=100 km/s"
        cfg["description"] = (
            f"SN_B interaction model: H/He shell M={mass} Msun, "
            f"v={vel} km/s, solar comp; {csm_setup} from the shell to "
            f"SN_A age {source_age:.0f} d at SN_B t=0 and "
            f"{profile_age:.0f} d at STELLA TcurA={ts:g} d; "
            f"Ni56 disabled. STELLA to {target_days:.0f} d."
        )
        cfg["observational_timing"] = timing

        cfg["run_cmfgen"] = False
        cfg["run_sirocco"] = False
        cfg["run_sedona"] = False
        cfg["run_superlite"] = False
        cfg["skip_stella"] = False
        cfg["reuse_stella"] = False
        cfg["submission"] = "normal"
        cfg["debug"] = False

        x = dict(parent["x_model_from_sn"])
        x["shell_mass_msun"] = mass
        x["shell_velocity_kms"] = vel
        x["t_start_days"] = ts
        x["wind_mdot_msun_yr"] = mdot
        x["wind_velocity_kms"] = 100.0
        x["wind_T_K"] = 2000.0
        # Clear any inherited two-region wind settings unless requested
        x.pop("wind_duration_days", None)
        x.pop("wind_mdot_bg_msun_yr", None)
        x.pop("density_matched_wind", None)
        x.update(extras)
        x["source_age_days"] = source_age
        cfg["x_model_from_sn"] = x

        st = dict(parent["stella"])
        st["target_days"] = target_days
        if nickel56 is not None:
            st["nickel56"] = nickel56
        cfg["stella"] = st

        out = dict(parent["output"])
        out["model_dir"] = f"./SN_B_dessart_grid/{name}"
        cfg["output"] = out

        header = (
            f"# SN_B interaction grid model: {name}\n"
            f"# Shell {mass} Msun @ {vel} km/s; {csm_setup}\n"
            f"# SN_A MJD {SN_A_LAUNCH_MJD:g}; SN_B t=0 MJD {SN_B_T0_MJD:g}\n"
            f"# SN_A age {source_age:.0f} d at SN_B t=0; "
            f"{profile_age:.0f} d at STELLA TcurA={ts:g} d; "
            f"age policy={timing['source_age_policy']}\n"
            f"# STELLA {target_days:.0f} d; STELLA only\n\n"
        )
        with open(config_path, "w") as f:
            f.write(header)
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)
        configs.append(str(config_path))
        print(f"  Written: {config_path}")

    print(f"\nGenerated {len(configs)} configs in {GRID_DIR}/")
    print("\nSubmit all with:")
    for c in configs:
        print(f"  python scripts/run_model_grid.py {c}")


if __name__ == "__main__":
    main()
