#!/usr/bin/env python3
"""
proxy_populations_v2.py

Generate CMFGEN ``*_IN`` departure-coefficient files that are a better
initial guess than ``LTE_EST`` for interacting supernovae, where
``T_gas``, ``T_rad`` and the dilution factor ``W`` are strongly decoupled.

Strategy hierarchy (preferred first)
------------------------------------
1. ``CONTINUATION`` — regrid converged DCs from a neighbouring CMFGEN
   model (requires ``--prev-model``).  Closest to how CMFGEN restarts
   actually work.
2. ``GREY`` + ``--stella-file`` — STELLA's own T_rad / tau / energy
   density replace the grey field; W_eff is then inferred from
   J/B(T_rad) when an independent T_rad + U_rad/L are available.
3. ``GREY`` — luminosity-anchored T_rad + ML93 + dilute-LTE.  The
   cold-start default; a radiation-field proxy, not a physical solution.
4. ``LTE`` — only for deliberately testing LTE convergence.

The default ``--tau-mode plane`` computes the honest grey continuum
optical depth.  ``--tau-mode nonmonotonic`` multiplies tau by a fixed
trapping factor inside velocity reversals — a numerical heuristic, not
a non-monotonic radiative-transfer treatment; velocity reversal changes
the CMF line-transfer problem, it does not multiply tau_continuum.

Physics (GREY / nebular, the production default)
------------------------------------------------
Ionisation (Mazzali & Lucy 1993; TARDIS ``PhiSahaNebular``)::

    Φ = W [δ ζ + W (1-ζ)] (T_e / T_R)^{1/2} Φ_LTE(T_R)

with ``T_e = T_gas``, ``T_R`` from the grey/STELLA radiation field, and
``W`` from geometric dilution plus a τ-blend.  Optically thick depths
(τ ≳ τ_trans) are mixed back to LTE(T_gas).

Excitation (Abbott & Lucy 1985; TARDIS dilute-LTE)::

    ground + metastable (E < 2 eV): Boltzmann at T_R
    other levels:                  × W

Departure coefficients written for CMFGEN::

    DC = n_dilute-LTE(T_R, W, n_ion) / n_LTE(T_gas, n_e^STELLA, DI)

The denominator is always LTE at the hydrodynamic gas temperature,
because CMFGEN multiplies DC by its own LTE populations at ``T_IN``.
Writing ``T_ion`` into the ``*_IN`` T column while computing LTE at
``T_ion`` is wrong and is no longer done.

``t_exc_min`` floors *excitation* temperature only.
``t_ion_floor`` is a numerical floor on ionisation temperatures.
Do not use a 10 000 K floor on ionisation — that over-ionises cool
outer ejecta.

Optional modes (all off by default)
-----------------------------------
``--excitation-mode nlte_se``: replace the single-T_exc dilute-LTE
excitation with a per-depth statistical-equilibrium solve over the
bound-bound transitions in ``*_F_OSCDAT`` (A, f, i-j).  Radiative rates
use a dilute blackbody J_nu = W B_nu(T_rad); collisional rates use a
van-Regemorter form at T_gas.  Each level then carries its own
thermalisation balance epsilon_u = C_u/(C_u + R_u): dense depths
approach LTE(T_gas), cold dilute depths approach dilute-LTE(T_rad, W).
Levels without transitions keep the dilute-LTE estimate.

``--dc-log-cap C``: symmetric log cap |log10 DC| <= C on the written
departure coefficients — limits the step away from LTE rather than the
absolute DC (preferred over ``--dc-max`` for cold models).

``--stella-file``: overlay the grey radiation field with STELLA's own
T_rad / tau / t_ion / luminosity / energy-density columns interpolated
onto the hydro grid.  When T_rad comes from a STELLA column *and* U_rad
or L is also supplied, W is inferred self-consistently as J/B(T_rad)
rather than taken from the geometric grey dilution.

``--strategy CONTINUATION --prev-model DIR``: regrid converged
departure coefficients from a previous CMFGEN model directory onto the
new hydro grid; species without a previous ``*_IN`` fall back to the
proxy build.

Saha inversion of STELLA n_e is retained as a *diagnostic overlay*
only.  It must not drive the ionisation temperature: the proxy models
a subset of species, so the inverted T is unphysical.
"""

# The implementation was split into the prepare_PROXY_POPULATIONS
# package (prepare_PROXY_POPULATIONS/utils/, one module per pipeline
# stage). This file is the stable entry point: it re-exports the full
# public API so all existing callers keep working unchanged:
#   import proxy_populations_v2                    (sys.path=scripts/)
#   python scripts/proxy_populations_v2.py ...     (CLI, submit_cmfgen.sh)
#   scripts_dir / "proxy_populations_v2.py"        (run_model_grid.py)

import sys
from pathlib import Path

# Keep the original module-level imports visible on this module too:
# callers sometimes reach e.g. proxy_populations_v2.np in tests.
import argparse  # noqa: F401,E402
import dataclasses  # noqa: F401,E402
import json  # noqa: F401,E402
import math  # noqa: F401,E402
import os  # noqa: F401,E402
import re  # noqa: F401,E402
import shutil  # noqa: F401,E402
from collections import defaultdict  # noqa: F401,E402
from dataclasses import dataclass, field  # noqa: F401,E402
from datetime import datetime  # noqa: F401,E402
from typing import Any, Dict, List, Optional, Sequence, Tuple  # noqa: F401,E402
import numpy as np  # noqa: F401,E402

_PKG_ROOT = Path(__file__).resolve().parent.parent
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))

from prepare_PROXY_POPULATIONS.utils import *  # noqa: F401,F403,E402
from prepare_PROXY_POPULATIONS.utils import __all__ as _utils_all  # noqa: E402
from prepare_PROXY_POPULATIONS.utils.driver import main  # noqa: F401,E402

__all__ = list(_utils_all) + ["main"]

if __name__ == "__main__":
    sys.exit(main())
