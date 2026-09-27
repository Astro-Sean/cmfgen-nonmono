-- SEDONA parameter file template
-- Settings are auto-configured by run_model_grid.py based on YAML config.
--   use_homologous_expansion=true  -> hydro_module="homologous", fuzz lines, NLTE off
--   use_homologous_expansion=false -> hydro_module="none", static snapshot, NLTE per config
--
-- Line opacity for non-monotonic velocity (SN 2024pba interaction models):
--   Direct Voigt (bound_bound=1) is the non-monotonic-safe path. Opacity is
--   evaluated per-cell in the comoving frame; a photon can resonate with the
--   same line in multiple cells at different radii (multiple resonance points).
--   Sobolev (line_expansion=1) uses tau = n_l * sigma * f * lambda * t_exp,
--   which assumes homologous expansion (t_exp = 1/(dv/dr) = const). This is
--   WRONG for non-monotonic velocity (shocks, reversals) because the Sobolev
--   length is not t_exp * v but 1/|dv/dr_sob|, which varies with radius and
--   ray angle. See AGENTS.md §40 and docs/non_monotonic_solver.md.
--   VERIFIED (2026-09-10): line_expansion_opacity in AtomicSpecies_opacities.cpp
--   uses `time` = t_now_ (global expansion time), NOT local dv/dr.
--   VERIFIED (2026-09-10): bound_bound_opacity evaluates Voigt profile per-cell
--   on the frequency grid; no velocity-gradient assumption.
-- DDMC ON (transport_use_ddmc=2): handles photon trapping in optically thick
--   zones where direct Voigt produces high peak optical depths. Without DDMC,
--   direct Voigt in static snapshots causes severe photon trapping (tau_line
--   >> tau_es). DDMC (Fleck & Canfield 1984) converts MC transport to diffusion
--   in zones with tau > threshold, enabling efficient escape.
-- LTE: opacity_use_nlte=0, opacity_atoms_in_nlte={} (Saha-Boltzmann populations)

grid_type    = "grid_1D_sphere"
model_file   = "@model_name@.mod"
hydro_module = "none"

-- SEDONA home directory
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"
data_fuzzline_file = sedona_home.."/data/kurucz_cd23_fuzz.hdf5"  -- Iron-group line blanketing

-- Transport properties
transport_nu_grid  = {1.0e14, 2.0e16, 0.005, 1}  -- Transport grid: q=0.005 (~1060 bins, finer to reduce zero-bin fraction)
transport_radiative_equilibrium  = 0  -- Hold T fixed from hydro (paper: no temperature update)
transport_boundary_in_reflect   = 1  -- Reflect photons at inner boundary (needed for core source)
transport_steady_iterate         = 25  -- Will be overridden by config if specified

-- DDMC (Discrete Diffusion Monte Carlo) for optically thick zones.
-- Direct Voigt line opacity produces high peak optical depths in static
-- snapshots. DDMC=2 (Fleck & Canfield 1984) converts MC to diffusion in
-- zones with tau > threshold, enabling efficient photon escape.
-- Without DDMC, direct Voigt causes severe photon trapping (tau_line >> tau_es).
transport_use_ddmc               = 2  -- 2=DDMC (Fleck & Canfield 1984)
transport_ddmc_tau_threshold     = 10.0  -- DDMC active in zones with tau > 10

-- Output spectrum frequency grid
-- Narrowed to optical range (2.5e14-1.2e15 Hz = 2500-12000 A) for maximum resolution.
-- With q=0.0005: ~3100 bins, ~3 A/bin, ~150 km/s velocity resolution.
-- The transport_nu_grid below still covers the full range for correct UV photon transport.
spectrum_nu_grid   = {2.5e14, 1.2e15, 0.0005, 1}  -- Spectrum grid: optical range, q=0.0005 (~3100 bins, ~3 A/bin)

-- Time of spectrum calculation (days since energy injection)
tstep_time_start = @snapshot_days@ * 3600.0 * 24.0

-- Radioactive particle emission (static snapshot)
particles_n_emit_radioactive = 1e4
particles_n_emit_thermal     = 5e8  -- Thermal photons from all zones (high count for fine spectrum grid + low escape fraction)
particles_max_total = 2e+09  -- Override default max for high-photon runs
particles_last_iter_pump     = 20  -- Average radiation field over last N iterations (more iters = more photons in spectrum)

-- Inner boundary source from hydro model innermost region
-- For STELLA models, this is set from the .swd snapshot data
-- For CMFGEN models, this is set from SN_HYDRO_DATA
-- Core photon source: emits blackbody photons from core_radius at core_temperature
-- In addition to thermal emission from all zones (particles_n_emit_thermal)
core_n_emit      = 0    -- Number of core photons per iteration (0 = no core source)
core_radius      = @r_inner@      -- cm, innermost radius from hydro model
core_temperature = @temperature@  -- K, innermost temperature from hydro model
core_luminosity  = @luminosity@   -- erg/s, L=4piR^2sigmaT^4 from innermost region
core_fix_luminosity = 1  -- Rescale spectrum by 1/escape_fraction to correct for absorbed photons

-- Opacity settings: auto-configured by run_model_grid.py
--   use_homologous_expansion=true  -> fuzz_exp=1, line_exp=1, bound_bound=0 (Sobolev OK for homologous)
--   use_homologous_expansion=false -> fuzz_exp=0, line_exp=0, bound_bound=1 (direct Voigt for non-monotonic)
-- Direct Voigt (bound_bound=1) is NON-MONOTONIC-SAFE: opacity is per-cell in
--   comoving frame, handles multiple resonance points in velocity reversals.
-- Sobolev (line_expansion=1) assumes homologous expansion (tau = n*sigma*f*lambda*t_exp).
--   WRONG for non-monotonic velocity: Sobolev length should be 1/|dv/dr_sob|, not t_exp*v.
-- Fuzz expansion provides iron-group line blanketing (uses 1/t_exp Sobolev approx).
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 0  -- Sobolev OFF (assumes homologous expansion; wrong for non-monotonic velocity)
opacity_bound_bound          = 1  -- Direct Voigt ON (non-monotonic-safe: per-cell opacity, multiple resonance points)
opacity_bound_free           = 1
opacity_free_free            = 1
opacity_epsilon              = 1
opacity_use_nlte             = 0  -- LTE: Saha-Boltzmann populations (set by run_model_grid.py from yaml use_nlte)
opacity_nlte_dens_min        = 1e-15  -- min density (g/cm^3) for NLTE; below this, fall back to LTE (avoids singular rate matrix in low-density outer zones)
opacity_atoms_in_nlte        = {}  -- Empty: NLTE fully disabled (LTE mode). GasState forces use_nlte=1 if list non-empty.

-- Line Doppler width: turbulent/macro velocity broadening for direct Voigt lines.
-- Thermal Doppler width (~1-13 km/s) gives ~10^3x higher peak optical depth than Sobolev,
-- causing severe line trapping in static snapshots. A turbulent velocity of 200 km/s
-- reduces peak tau by ~40x (H) to ~300x (Fe) while preserving distinct P-Cygni features.
-- Units: cm/s (200 km/s = 2.0e7 cm/s). 0 = use thermal Doppler width.
line_velocity_width          = 2.0e7  -- 200 km/s turbulent broadening (override thermal Doppler width)

-- Output files
output_write_radiation = 1
output_file = "spectrum.h5"
