#!/usr/bin/env python3
"""Shared log-banner formatting for the modelling pipeline.

run_model_grid.py, prepare_VADAT.py and prepare_IN_INTS.py use these
helpers so that output written to submit.log / batch.log / OUTGEN has a
consistent, scannable structure:

    ============================ SECTION TITLE =============================
    --------------------------- sub-section -------------------------------

Plain ASCII only: OUTGEN is a fixed-format Fortran log.
"""

WIDTH = 78


def _banner(title: str, fill: str) -> str:
    """Center `title` in a WIDTH-wide line padded with `fill` characters."""
    text = f" {title.strip()} "
    pad = max(4, WIDTH - len(text))
    left = pad // 2
    return fill * left + text + fill * (pad - left)


def header(title: str) -> str:
    """Top-level section banner: ======== TITLE ========."""
    return _banner(title, "=")


def subheader(title: str) -> str:
    """Sub-section banner: -------- TITLE --------."""
    return _banner(title, "-")


def print_header(title: str = "") -> None:
    """Print a blank line followed by a top-level banner."""
    print()
    print(header(title))


def print_subheader(title: str = "") -> None:
    """Print a blank line followed by a sub-section banner."""
    print()
    print(subheader(title))


class Tee:
    """Write to stdout and a file simultaneously (like tee(1)).

    Used so that preparation-script output appears on the console/batch.log
    and is also appended to OUTGEN, keeping the whole run log in one place.
    """

    def __init__(self, stdout, f):
        self._stdout = stdout
        self._file = f

    def write(self, s):
        self._stdout.write(s)
        self._file.write(s)
        self._file.flush()

    def flush(self):
        self._stdout.flush()
        self._file.flush()
