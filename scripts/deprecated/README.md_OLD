# Deprecated Scripts

These scripts reference HERACLES, a hydrodynamics code that does not fit the
science case for this project (interacting Type IIn supernovae with CSM
interaction).

**Deprecated:** 2026-07-09

**Reason:** HERACLES lacks the radiation-hydrodynamics capabilities needed to
self-consistently model shock breakout and CSM interaction in Type IIn SNe.
The project now uses the STELLA -> CMFGEN pipeline exclusively.

**Active alternatives:**

| Deprecated script | Active replacement |
|---|---|
| `convert_heracles_to_cmfgen.py` | `../convert_stella_swd_to_cmfgen.py` |
| `convert_heracles_to_cmfgen_completed.py` | `../convert_stella_swd_to_cmfgen.py` |
| `convert_heracles_to_cmfgen_wrapper.py` | `../run_model_grid.py` (with `hydro_backend: stella`) |
| `generate_heracles_namelist.py` | Not needed (STELLA uses `stella_job_template.slurm`) |
| `find_closest_snapshot.py` | STELLA `.swd` files selected via `--target-days` in the converter |

These files are retained for reference only and should not be used for new runs.
