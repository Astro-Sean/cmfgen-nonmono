# BLANK-render continuum rescue (J_RESCUE)

This note explains what makes an OBS BLANK render fail on dense
interacting-SN models, and how the J_RESCUE machinery in
`obs/COMP_JCONT_V4.INC` fixes it.  The full iteration-level audit is
`docs/audits/CMF_FLUX_P8_2_ES_CONVERGENCE_AUDIT.md` section 7; this
document is the short version aimed at users and reviewers.

## 1. The failure

On `mc170_H10` (the SN_A candidate), the BLANK observed-frame render
(`GLOBAL_LINE=BLANK`, `NUM_ES=2`, `COH_ES=F`, `N_TYPE=N_ON_J`)
completed the continuum march but `obs_fin_10` carried an infrared
intensity tail rising to ~1e61 Jy - twelve orders above the
physical scale.  The SOB render of the same model (`obs_sob_es1`)
logged 153 convergence warnings but returned a bounded spectrum
(max ~2.1e6 Jy), and the KNOP render (`obs_knop`) converged
everywhere (max ~1.6e6 Jy).  Only the BLANK render's output was
poisoned.

The failure is localised to a band around FL ~ 0.031-0.053
(nu ~ 3.1-5.3e13 Hz, lambda ~ 5.6-9.6 um).  Below and above the
band the same solve converges on its own.

## 2. Root cause: an amplitude-blind fixed point

In the IR band the atmosphere is nearly pure electron scattering:
the divergence dump at FL=0.0532 shows `ESEC~2e-5` dominating the
absorptive opacity `CHI~2e-6`, with true emissivity `ETA~1e-21`.

The formal-solve source is `TA = ETA + CHI_SCAT*RJ`.  With
`ETA~0`, `TA` is nearly homogeneous in J: the transfer equation
`J = Lambda*TA` becomes `J = A.J + b.ETA ~ A.J`.  Any amplitude is
a near-fixed point - the iteration cannot distinguish `J=1e-7`
from `J=1e8` because the equations are almost scale-free there.

The loop's only convergence test is on the Eddington factor
`FEDD = K/J`, which is a ratio.  If J and K drift upward together,
`FEDD` converges while the amplitude runs away.  The stored J then
ratchets upward ~1.14x per frequency: each frequency stores
slightly-larger-but-still-"converged" garbage, `JNU_PREV` (the
warm-start field for the next frequency) absorbs it, and the next
iterate grows from there.  EDDFACTOR records in the band show Jmax
ramping smoothly from ~1e-4 to ~3e8 across adjacent frequencies.

Every iteration-level discriminator fails by construction here:

- `FEDD` convergence: amplitude-blind, as above.
- `FEDD > 2` validity: a smooth multiplicative ratchet keeps f in
  bounds (f is a ratio of two drifting fields).
- `L > 20` exit status: resc5 showed the poisoned band frequencies
  exit *converged* - the residual `T1` is small because a wrong-
  amplitude fixed point is still a fixed point.
- Residual size at exit: `T1 <= 0.33` across the whole band, no
  separation between sane and garbage fields.

The KNOP moment solver escapes this because it couples all
frequencies: the neighbouring optically-thick frequencies anchor
the amplitude, and its stored field is flat J ~ 1e-7-1e-5 across
the whole band - the reference answer.

## 3. Why the simple fixes did not work

- `DAMP_FJ=0.3` under-relaxation: a fixed point is still a fixed
  point after damping.  Tested (obs_blank_damp): no improvement;
  damping is off in the final config (`DAMP_FJ=1.0`).
- Luminosity-derived `DTDR` (the nonmono cool-inner-zone boundary
  fix, extended to `J_RESCUE`): `DBB` went from -2.0e-12 to
  +1.5e-11 and the runaway was identical (obs_blank_resc3).  The
  boundary sign is not the driver.
- Convergence-gated restore (any `L>20` exit restores): froze J
  at the first frequency - stiff frequencies routinely exit L>20
  with sane fields (stock semantics keep them; sob_es1 has 153
  such exits and a good spectrum).  Reverted.
- `RESC_CAP = MAX(1e4*prev, 1e-3)`: the constant is a *floor* -
  once `JNU_PREV` inflated past ~1e-7, the `1e4*prev` term
  recalibrated the cap upward and the ratchet walked underneath.

## 4. The fix

`J_RESCUE` (now default T) adds a field-validity check at the
exit of each frequency's J/f iteration, plus an in-loop runaway
detector.  Two caps split the roles:

- `RESC_CAP` (stored-field validity + `RJ_ES` sanitize):
  `MIN(1e4 * MAXVAL(ABS(JNU_PREV)), RESC_JMAX)`.  The relative
  term catches single-step explosions and garbage in dark
  regions; `RESC_JMAX` (new key, default `1e-3` code units) is an
  absolute ceiling that no converged field on this model class
  approaches - legitimate J peaks at 4.8e-5 across the whole grid
  (KNOP and SOB agree), and the largest legitimate adjacent-
  frequency Jmax ratio measured is ~44x, so `1e4x` has ~200x
  margin and cannot false-fire on a real transition.

- `RESC_CAP_IT` (in-loop early exit only): the previous loose
  bound `MIN(1e10, MAX(1e4*prev, 1e-3))`.  Iterates may
  transiently overshoot a tight cap, so the early exit stays on
  runaway scale - it trims wasted iterations, it does not decide
  validity.

At loop exit, if the stored field is non-finite or exceeds
`RESC_CAP` at any depth, the code restores the full previous-
frequency field - `JNU_PREV`, `FEDD_PREV`, `GEDD_PREV`,
`N_ON_J_PREV`, `RSQHNU_PREV`, `HBC_PREV`, `NBC_PREV` - warns, and
continues the march.  Mild `L>20` exits with bounded fields keep
their fields (stock semantics), so the warm-start chain is only
broken where the field is actually invalid.

Under incoherent electron scattering (`COH_ES=F`), `RJ_ES` is
read from `ES_J_CONV` for the emissivity term.  If that field is
non-finite or exceeds `RESC_CAP` it is replaced by the local
coherent `RJ` - incoherent falls back to coherent at that
frequency only, with a throttled warning.

`J_RESCUE` also joins the `DO_NONMONO_*` gate for the luminosity-
derived `DTDR` at the inner boundary and the `DBB>=0` floor - a
model that needs rescue may also carry an inverted inner-zone
temperature gradient, and the floor is the documented physical
limit (no diffusion flux through the boundary).

The first rescued frequency dumps a per-depth diagnostic block
(`R, V, ED, CHI, ETA, ESEC, CHI_SCAT, ZETA, RJ, FEDD, JNU_PREV`)
so the origin of each failure band can be inspected.

## 5. Controls and defaults

`CMF_FLUX_PARAM` keys (all optional since 14-Nov-2026):

| key | default | role |
|---|---|---|
| `J_RESCUE` | `T` | master switch; explicit `F` restores stock abort-on-divergence behaviour |
| `RESC_JMAX` | `1e-3` | absolute J ceiling in code units; raise for models whose legitimate J exceeds ~5e-5 |
| `DAMP_FJ` | `1.0` | under-relaxation for the f/J update; kept off - it does not fix the eigenproblem |
| `ALLOW_NONMONO_OBS` | `F` | required for BLANK renders on non-monotonic RVTJ models (pipeline stamps it) |

The grid pipeline (`scripts/prepare_CMF_PARAM.py`) stamps
`J_RESCUE=T`, `DAMP_FJ=1.0`, `ALLOW_NONMONO_OBS=T` into every
non-monotonic BLANK and KNOP render control file; explicit
`cmf_flux_param_overrides` still win.  `RESC_JMAX` rides the code
default.

## 6. Validation

`obs_blank_resc6` (job 13158209) on mc170_H10, compared against
the fully-converged KNOP render and the SOB render:

| render | max I (Jy) | points >1e10 | verdict |
|---|---|---|---|
| `obs_knop` (reference) | 1.6e6 | 0 | converged everywhere |
| `obs_sob_es1` | 2.1e6 | 0 | converged, 153 warnings |
| `obs_blank_resc5` | 5.9e18 | 4109 | poisoned (pre-fix) |
| `obs_blank_resc6` | 2.3e6 | 0 | bounded |

Ratio to the converged KNOP spectrum:

- FL 0.5-5 (optical/UV): **1.00** - the rescue is invisible in
  the region the photometry sees.
- FL < 0.031 (deep IR, self-converged): median 1.36.
- FL 0.031-0.053 (restored band): median ~4.3, max ~72 - the
  flat last-converged field sits slightly above KNOP's local
  level.  Bounded, logged, and ~10 orders below the previous
  garbage; the band sits at lambda >~5.6 um, outside every
  photometric band used in the light-curve analysis.
- 57 mild `L>20` warnings kept (stock semantics, matching
  sob_es1's pattern); zero `RJ_ES` sanitize events.

Regression on the flipped default: P8.1 battery 116/116, P8.2
battery 20/20.

## 7. What the rescue is not

A rescued frequency is a flagged fallback - the previous
frequency's converged field written into a band the stock
iteration cannot solve - not an independently converged
solution.  The band edge and every rescue are logged in OUTGEN
(`J_RESCUE:` lines) and counted in `NRESC`.  Where the model
needs the IR band itself solved, the answer is the KNOP moment
solver (`obs_line_treatment='knop'` or `DO_NONMONO_KNOP`), which
converges these frequencies cleanly; the rescue exists so that
the BLANK render always returns a usable, bounded spectrum
instead of aborting or emitting a poisoned tail.

## 8. Limitations

- All J_RESCUE machinery lives in the `ND` (non-`INC_GRID`)
  branch of `COMP_JCONT_V4.INC` - the path every production
  render on these models takes.  `INC_GRID` renders keep stock
  behaviour under `J_RESCUE=T`.
- At the first frequency there is no previous field to restore;
  a divergence there still stops the run (honest behaviour -
  there is nothing valid to fall back to).
- `RESC_JMAX` is a code-units scale tuned on this model class
  (legitimate J <= ~5e-5).  A model family with systematically
  larger J should set `RESC_JMAX` explicitly rather than rely on
  the 1e-3 default.
