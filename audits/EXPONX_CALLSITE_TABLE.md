# EXPONX / d_EXPONX_dX call-site table

Source of truth: `cmfgen_jfnk` @ commit `a029fd3709078a60cff04bc0102cebad52f28a48`
(main, origin/main). Audited model was built from `9c5d5a25` — same guard logic
verified present at that commit's line numbers via run output (`NEG_OPAC` log).

Definitions: `GAMH = CHIL * 3e-10 * R / V / FL` (angle-free Sobolev tau,
CHIL in CMFGEN opacity units, FL in 1e15 Hz). `GAM = GAMH / T1` where
`T1 = MAX(|1+sigma*mu^2|, -0.2*(1+sigma), 1e-4)` for sigma<0 (DAH15 floor),
`T1 = |1+sigma*mu^2|` otherwise.

## Production call sites of EXPONX

| file:line | routine | argument | can be <0? | active in DO_NONMONO_SOB | active in stock path | neg-CHIL handling |
|---|---|---|---|---|---|---|
| subs/sobjbar_sim.f:181 | SOBJBAR_SIM | GAM=GAMH_sum/T1 | YES in principle | yes (rates+BETA) | yes | upstream crush on summed CHIL (see guards) |
| subs/sobjbar_sim.f:182 | SOBJBAR_SIM | d_EXPONX_dX(GAM) | same arg | yes | yes | same |
| subs/sobew.f:188 | SOBEW | GAM=GAMH_sum/T1 | YES in principle | yes (EW diag, SIGMA_TRUE) | yes (SIGMA) | fed same sanitized CHIL |
| subs/sobew_grad.f:181 | SOBEW_GRAD | GAM | YES | not called in prod | legacy | same CHIL |
| subs/sobew_grad_v2.f:232 | SOBEW_GRAD_V2 | GAM | YES | obs/disp only | obs/disp | same CHIL |
| subs/sobew_grad_v3.f:235 | SOBEW_GRAD_V3 | GAM | YES | TEST-ONLY (drv_sob_floor) | no | same CHIL |
| subs/sob_rh_couple.f:317,320 | SOB_RH_COUPLE | EXPONX(TI)*EXPONX(TJ), per-line GAMH_L | NO - guarded | yes when DO_SOB_RH | n/a | explicit `GAMH_L<=0` and `NEG_OPACITY` skips (lines ~309-315) |
| disp/subs/sobjbar.f:127,133 | SOBJBAR (disp) | GAM, signed T1 denom | YES + no T1 floor | display tool only | display | none in routine |

Callers: SOBJBAR_SIM <- SUB_SOB_LINE{,_V2,_V3} <- cmfgen_sub.f:4186-4201 and
JFNK sweep replay ~7294-7350. SOBEW <- SUB_SOB_LINE_V3:194-198 and
cmfgen_sub.f:2269 (dielectronic EW diag). SOB_RH_COUPLE <-
sub_sob_line_v3.f:172,176.

## Negative-opacity machinery

| site | action |
|---|---|
| new_main/cmfgen_sub.f:4101-4108 | CHIL(I)=sum CHIL_MAT over SIM group (SOBOLEV or NUM_SIM_LINES=1) |
| new_main/cmfgen_sub.f:4113-4130 | CHECK_LINE_OPAC+SOBOLEV: if CHIL*3e-10/FL*R/V < -0.5 -> NEG_OPACITY(I)=T, log to NEG_OPAC (last iter), CHIL(I)=1.0 |
| new_main/cmfgen_sub.f:4137-4141 | crush: CHIL_MAT(I, all lines)=1/NUM_SIM_LINES at flagged cells |
| new_main/cmfgen_sub.f:4148-4160 | ELSE IF(CHECK_LINE_OPAC) non-Sobolev branch: line-center opacity vs continuum check, CHIL=1e-4*ESEC/PROF; CHIL_MAT NOT crushed (CMF ignores CHIL_MAT) |
| new_main/cmfgen_sub.f:2543-2580 | NEG_OPAC_OPT: ESEC_CHK floors total CHI at 0.1*CHI_SCAT; SRCE_CHK alt; sets NEG_OPAC_FAC=0 |
| new_main/cmfgen_sub.f:6302-6360 | identical block inside JFNK_RAD_SWEEP_V1 (+ KN_NO_GAIN floor for DO_NONMONO_KNOP) |
| new_main/cmfgen_sub.f:7294-7306 | JFNK coupled-SOB replay duplicates the crush verbatim |
| new_main/mod_subs/sub_sob_line_v3.f:210,297 | NEG_OPACITY(I) -> VB_SIM=0, VB_2=0 (Jacobian term excised) |
| new_main/mod_subs/do_var_cont.f:644-872 | NEG_OPAC_FAC multiplies line-opacity terms in variation -> 0 at flagged cells |
| obs/cmf_flux_sub_v5.f:2880-2900 | identical CHECK_LINE_OPAC crush in observed-frame code |
| obs/cmf_flux_sub_v5.f:3129+ ADD_SOB_LINE_PROF | NEG_OPACITY depths dropped from observed-spectrum line profiles |
| lte_hydro/lte_sub.f:893 | separate LTE hydro code, same NEG_OPAC logging pattern |
| rd_control_variables.f:815 | CHECK_LINE_OPAC default = TRUE (VADAT key CHK_L_POS) |
| VADAT:108 | NEG_OPAC_OPT=ESEC_CHK (model config) |

## Reachability verdict

Per-line CHIL_MAT<0 NEVER reaches EXPONX individually — GAM is built from
the group-summed CHIL only; CHIL_MAT enters only the apportionment
(ST_ON_SL, VB, VC, BETAC) at sobjbar_sim.f:204-218.

Summed CHIL<0 reaches EXPONX only in the window
-0.5 < angle-free tau < 0 (below -0.5 the crush fires). With T1>=1e-4 the
resulting GAM can in principle reach -5000 — an overflow-capable latent
window — but the run's own NEG_OPAC log shows all flagged cells at
TAU(Sob) in [-0.99, -0.5], and the model completed, so realized GAM stayed
positive at flagged cells and bounded elsewhere.
