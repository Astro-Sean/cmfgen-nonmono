#!/usr/bin/env python3
"""test_iter_ledger -- T-06 per-iteration work ledger (audit F-22/F-24).

Compiles the production ITER_LEDGER_MOD (new_main/mod_subs/
iter_ledger_mod.f90) with a small driver and asserts on the emitted
ITER_LEDGER CSV:

  * header written exactly once (file-new only)
  * one row per LEDG_WRITE_ITER call
  * every counter/timer/state column reconciles with the values the
    caller supplied (the "ledger totals reconcile with OUTGEN
    counters" acceptance: the row mirrors the same state ITER_SUM
    prints, plus the per-iteration work counters)
  * n_dem counts '*demote' tokens inside ITER_REASON
  * jfnk_ev/jfnk_jvp are per-iteration deltas vs the reset snapshot
  * rss_mb is a real positive MB value on Linux (or -1 elsewhere)
  * a call before LEDG_ITER_RESET leaves no file (stock path silent)

Run: python3 unit_tests/python/test_iter_ledger.py
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import common

REP = common.Report("test_iter_ledger (T-06 ITER_LEDGER)")

WORK = common.WORKDIR / "iter_ledger"
WORK.mkdir(parents=True, exist_ok=True)


def run_driver(mode):
    binpath = common.compile_and_link_prod(
        "drv_iter_ledger",
        [
            common.TREE / "tools" / "set_kind_module.f",
            common.TREE / "new_main" / "mod_subs" / "iter_ledger_mod.f90",
            common.F90DIR / "drv_iter_ledger.f90",
        ],
    )
    ledger = WORK / "ITER_LEDGER"
    if ledger.exists():
        ledger.unlink()
    r = common.run(binpath, mode, cwd=WORK)
    return r, ledger


def parse_ledger(path):
    text = path.read_text().splitlines()
    header = [l for l in text if l.startswith("# ")]
    rows = [l for l in text if l and not l.startswith("#")]
    cols = header[0][2:].split(",") if header else []
    parsed = [dict(zip(cols, row.split(","))) for row in rows]
    return text, cols, parsed


# --------------------------------------------------------------------
# 'on' leg: armed iteration, two writes, distinct caller state
# --------------------------------------------------------------------
r, ledger = run_driver("on")
REP.check("driver ran (on)", r.returncode == 0)
REP.check("ITER_LEDGER created", ledger.exists())

text, cols, rows = parse_ledger(ledger)
EXPECTED_COLS = (
    "iter,mode,next,ba,frozen,held,reason,"
    "n_fg,n_mom,n_rf,np_rf,ba_build,ba_reuse,"
    "jfnk,jfnk_cpl,jfnk_ev,jfnk_jvp,n_rej,n_dem,n_frz,"
    "n_frozen,n_rel,n_raw_rel,"
    "maxch,maxch_app,maxatt,se_res,se_res2,trust,"
    "t_iter,t_rad,t_rf,t_pops,rss_mb"
).split(",")
REP.check("header written once", len([l for l in text if l.startswith("#")]) == 1)
REP.check("header columns", cols == EXPECTED_COLS,
          detail=f"got {cols}")
REP.check("two rows", len(rows) == 2, detail=f"{len(rows)} rows")

if len(rows) == 2:
    a, b = rows
    # Row 1: caller-supplied state and counters must land verbatim.
    REP.check("row1 identity", a["iter"] == "12"
              and a["mode"] == "LAMBDA_TFIX" and a["next"] == "FULL"
              and a["ba"] == "LAMBDA")
    REP.check("row1 flags", a["frozen"] == "1" and a["held"] == "1")
    REP.check("row1 reason verbatim",
              a["reason"] == "gs_immature_demote+hold_pops")
    REP.check("row1 solve counters", a["n_fg"] == "1300"
              and a["n_mom"] == "1300")
    REP.check("row1 refresh", a["n_rf"] == "1" and a["np_rf"] == "7")
    REP.check("row1 ba build/reuse",
              a["ba_build"] == "1" and a["ba_reuse"] == "0")
    REP.check("row1 jfnk idle", a["jfnk"] == "0" and a["jfnk_cpl"] == "0")
    REP.check("row1 jfnk eval deltas",
              a["jfnk_ev"] == "12" and a["jfnk_jvp"] == "12")
    REP.check("row1 demote counted", a["n_dem"] == "1")
    REP.check("row1 freeze edge", a["n_frz"] == "1")
    REP.check("row1 cumulative state", a["n_frozen"] == "3"
              and a["n_rel"] == "2" and a["n_raw_rel"] == "1")
    REP.check("row1 demand/residual",
              abs(float(a["maxch"]) - 4.5e2) < 1e-3
              and abs(float(a["maxch_app"]) - 1.1e2) < 1e-3
              and abs(float(a["se_res"]) - 3.2e4) < 1e-3
              and abs(float(a["trust"]) - 0.5) < 1e-3)
    # t_iter = wt1 - ledg_t0 = 300 - 100 = 200; t_pops = 300 - 200 = 100
    REP.check("row1 stage times",
              abs(float(a["t_iter"]) - 200.0) < 0.5
              and abs(float(a["t_rad"]) - 55.25) < 0.5
              and abs(float(a["t_rf"]) - 12.5) < 0.5
              and abs(float(a["t_pops"]) - 100.0) < 0.5)

    # Row 2: second write reuses the same reset -> eval deltas are
    # still measured against the iteration-start snapshot, not the
    # previous row.
    REP.check("row2 identity", b["iter"] == "13" and b["mode"] == "FULL")
    REP.check("row2 ba reuse", b["ba_build"] == "0" and b["ba_reuse"] == "1")
    REP.check("row2 jfnk fired",
              b["jfnk"] == "1" and b["jfnk_cpl"] == "1")
    REP.check("row2 jfnk eval deltas vs reset",
              b["jfnk_ev"] == "20" and b["jfnk_jvp"] == "18")
    REP.check("row2 no demote token", b["n_dem"] == "0")
    REP.check("row2 reject counted", b["n_rej"] == "1")
    REP.check("row2 t_iter", abs(float(b["t_iter"]) - 260.0) < 0.5)

    # RSS: positive MB on Linux, or -1 where /proc is unavailable.
    rss = float(a["rss_mb"])
    REP.check("rss present or sentinel", rss > 0.0 or rss == -1.0,
              detail=f"rss_mb={rss}")

# --------------------------------------------------------------------
# 'off' leg: no reset armed -> the internal LEDG_ACTIVE guard leaves
# no file at all (the stock-path acceptance: zero ledger output)
# --------------------------------------------------------------------
r, ledger = run_driver("off")
REP.check("driver ran (off)", r.returncode == 0)
REP.check("no ITER_LEDGER without reset", not ledger.exists())

sys.exit(REP.summary())
