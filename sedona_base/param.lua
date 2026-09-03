-- SEDONA parameter file template
-- Settings are auto-configured by run_model_grid.py based on YAML config.
--   use_homologous_expansion=true  -> hydro_module="homologous", fuzz lines, NLTE off
--   use_homologous_expansion=false -> hydro_module="none", static snapshot, NLTE per config
--
-- Line opacity: Sobolev ON with local |dv/dr| (non-monotonic velocity fix in source)
-- DDMC OFF: causes photon accumulation with many low-tau zones (standard MC is stable)
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
transport_steady_iterate         = 25  -- Will be overridden by config if specified

-- DDMC disabled: causes photon accumulation with many low-tau zones (standard MC is stable)
transport_use_ddmc               = 0  -- 0=off (standard MC: stable escape fraction, no accumulation)
transport_ddmc_tau_threshold     = 1000.0  -- Effectively unused with DDMC off

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
--   use_homologous_expansion=true  -> fuzz_exp=1, line_exp=1, use_nlte=0
--   use_homologous_expansion=false -> fuzz_exp=0, line_exp=1, bound_bound=0, use_nlte per config
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 1  -- Sobolev ON with local |dv/dr| (non-monotonic velocity safe, source code fix)
opacity_bound_bound          = 0  -- Direct Voigt OFF (causes photon trapping in optically thick ejecta)
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
