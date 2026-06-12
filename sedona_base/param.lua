-- SEDONA parameter file template
-- NLTE (all 8 species) + full opacities
-- This is a template - model_file and tstep_time_start should be updated per model

grid_type    = "grid_1D_sphere"
model_file   = "@model_name@.mod"
hydro_module = "none"

-- SEDONA home directory
sedona_home        = os.getenv('SEDONA_HOME')
defaults_file      = sedona_home.."/defaults/sedona_defaults.lua"
data_atomic_file   = sedona_home.."/data/cmfgen_levelcap100.hdf5"

-- Transport properties: improved convergence
transport_nu_grid  = {0.8e14, 1.0e16, 0.001, 1}
transport_radiative_equilibrium  = 1
transport_steady_iterate         = 5   -- 5 iterations, NLTE is slow

-- Output spectrum frequency grid
spectrum_nu_grid   = {0.8e14, 1.0e16, 0.002, 1}

-- Time of spectrum calculation (days since energy injection)
-- Updated by run_model_grid from SN_HYDRO_DATA Time(days) since explosion
tstep_time_start = @snapshot_days@ * 3600.0 * 24.0

-- Radioactive particle emission (static snapshot)
particles_n_emit_radioactive = 1e4
particles_last_iter_pump     = 10

-- Inner boundary source from tau-trimmed SN_HYDRO_DATA innermost region
core_n_emit      = 2e5
core_radius      = @r_inner@      -- cm, innermost radius from hydro model
core_temperature = @temperature@  -- K, innermost temperature from hydro model
core_luminosity  = @luminosity@   -- erg/s, L=4piR^2sigmaT^4 from innermost region

-- Opacity settings: FULL OPACITIES + NLTE (all species)
opacity_grey_opacity         = 0
opacity_electron_scattering  = 1
opacity_fuzz_expansion       = 0
opacity_line_expansion       = 1
opacity_bound_bound          = 1
opacity_bound_free           = 1
opacity_free_free            = 1
opacity_epsilon              = 1
opacity_use_nlte             = 1 
opacity_atoms_in_nlte        = {1, 2, 6, 8, 14, 26, 27, 28}  -- species in X1-test-debug.mod (H, He, C, O, Si, Fe, Co, Ni)
opacity_use_collisions_nlte  = 1

-- Output files
output_write_radiation = 1
output_file = "spectrum.h5"
