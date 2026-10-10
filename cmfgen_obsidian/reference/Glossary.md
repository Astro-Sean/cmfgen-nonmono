# Glossary

Terms and symbols as used in this tree and the Knop/CMFGEN literature.

## Symbols

- `sigma` ($\sigma$) = $d\ln V/d\ln R - 1$. Homology: 0. Compression
  zone: < -1. `SIGMA` is the clamped infrastructure copy;
  `SIGMA_TRUE` is unclamped.
- `GAM` ($\Gamma$) = $(V/R)(1+\sigma\mu^2)/c$ - the velocity-coupling
  coefficient of the comoving-frame equation; its sign picks the
  frequency upwind direction.
- `chi_hat` ($\hat\chi$) = $\chi + |\Gamma|/d\ln\nu$ - Knop's
  generalized opacity; positive in compression zones where plain
  `chi + GAM/dln nu` can go negative.
- `Q` = $\Gamma/(d\ln\nu\,\hat\chi)$ - the coupling weight;
  `|Q|/(chi+|Q|) -> 1` marks marginally contractive links.
- `I_NEXT` - lagged intensity store at the not-yet-computed upwind
  frequency (compressing cells only). `J_NEXT_FRQ`/`H_NEXT_FRQ` are
  the moment-solver analogues.
- `STEQ` - statistical-equilibrium residual vector (and matrix
  workspace).
- `MAXCH` - max fractional population change. Two flavours matter:
  **demanded** (raw correction the solve wants) vs **applied**
  (post-cap, post-rollback). The `1e7` sentinel = saturated demand.
- `BA` - the block-diagonal-per-depth approximate Jacobian
  ("the BA matrix"); also the name for a full/Newton iteration.
- `ND`, `NT`, `NCF` - depth points, unknowns per depth, frequency
  count. m15: `N = NT*ND ~ 1.3e5`, NCF ~ 4-9e4.

## Solver vocabulary

- **Lambda iteration**: fixed-point SE update against the current J;
  cheap, bounded by `MAX_LAM_COR`.
- **BA / full step**: the linearized solve `A x = F` via the BA
  matrix.
- **resacc** (`DO_BA_RES_ACC`): residual acceptance - take the BA
  step, re-evaluate the SE residual, roll back on growth.
- **GS freeze** (`GS_FREEZE`): the controller state where the lagged
  field is not yet a fixed point; forces bounded Lambda updates.
- **Mature field**: the lagged store's drift below the tolerance/bulk
  gates for `GS_REL_STABLE` consecutive sweeps.
- **MOM_TRIAL**: bounded resacc-gated BA probe allowed after
  `MOM_TRIAL_MIN` consecutive marginal-frozen iterations.
- **Refresh**: `FG_KNOP_REFRESH_V1`/`MOM_KNOP_REFRESH_V1` - re-solve
  the lagged coupling to its fixed point after each sweep.
- **Ordered march** (`KN_RAY_ORDERED`/`SWEEP_FLD_ORD`): march along
  the characteristic; frequency inside each row in the row's upwind
  direction.
- **C-corner**: Knop's spatially-downwind x wavelength-upwind corner -
  the only term that must stay lagged in the ordered march.
- **Frequency matching** (`DO_NONMONO_FREQ_MATCH`): remap the
  frequency grid so the upwind sample lands on the characteristic
  foot (Ceulemans 2025).
- **JFNK**: Jacobian-free Newton-Krylov; GMRES + finite-difference
  Jacobian-vector products; here in log-x unknowns with PTC.
- **Frozen vs coupled residual**: frozen uses stored radiative terms +
  `CMAT_RAD` matvec (cheap); coupled re-runs the radiation sweep per
  eval (expensive, sees dJ/dn).
- **PTC**: pseudo-transient continuation - `G(x) = F(x) + x*s/dtau`.
- **EW forcing**: Eisenstat-Walker adaptive GMRES tolerance (floored
  at 3e-4).
- **F(0) test**: the residual identity check before a Newton solve.
- **ALO / Lambda***: approximate Lambda operator; here diagonal
  `dJ/dchi, dJ/deta` from the formal solver (`TX(K,K,1/2)`).
- **NG acceleration**: CMFGEN's Ng outer-sequence acceleration
  (`do_ng_accel_band_v4`); acts on the outer fixed-point map, gated
  off during freeze/warmup.

## Process / infrastructure

- **Great iteration**: one pass of the full loop (transfer sweep + SE
  update).
- **Requeue**: SLURM segment resubmission; process-local state (the
  lagged store) is lost unless checkpointed.
- **Run dir**: `cmfgen_validation*/...` under `/ptmp`; config in
  `SN2024pba/configs/*.yaml`; jobs via `run_model_grid.py`.
- **OBS stage**: `cmf_flux.exe` + `obs_frame.exe` producing
  `obs_fin_10`/`OBSFLUX`; requires `USE_LAM_ES=T` heritage.

## File tags

- `#status/verified` - demonstrated in a model run or bitwise audit.
- `#status/partial` - validated in part; limits documented.
- `#status/unverified` - code reviewed/compiled only.
