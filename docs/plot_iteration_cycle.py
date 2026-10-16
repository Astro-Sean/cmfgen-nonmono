#!/usr/bin/env python3
"""
Great-iteration cycle of this CMFGEN fork, drawn as a flowchart with
the non-monotonic additions marked against the unchanged steps.

One pass through the diagram is one great iteration of the main
loop.  Boxes are coloured by provenance: grey steps are code this
fork does not touch, amber steps are existing routines carrying new
flag-gated options, and blue steps are new code.  Every blue/amber
feature is off by default.

Box widths are measured from the rendered text, so no label can
overflow its box.  The canvas is widened as needed to fit the legend
row under the title.

Output: iteration_cycle.png (next to this script).
"""


from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------
# Box styles by provenance.
# ---------------------------------------------------------------
STYLE = {
    "stock": dict(fc="#f2f2f2", ec="#777777"),   # unchanged CMFGEN
    "mod":   dict(fc="#fdf1dc", ec="#c9821e"),   # existing + new options
    "new":   dict(fc="#e2eefb", ec="#1f6fb4"),   # new code in this fork
}

FS_TITLE = 11
FS_BODY = 9.2
LH_TITLE = 0.34        # axis units per title line
LH_BODY = 0.265        # axis units per body line
PAD_V = 0.16
BOX_X = 0.75           # left edge of the main column
TXT_PAD_L = 0.28       # text indent inside a box
TXT_PAD_R = 0.25       # clearance between text and box edge
RMARGIN = 0.3          # right edge margin
GAP = 0.45             # vertical gap left for arrows between boxes
DPI = 160
FIG_W_IN = 13.6
AX_W_FRAC = 0.775      # default subplots_adjust: right 0.9 - left 0.125


BOXES = [
    dict(
        kind="stock",
        title="Opacities",
        lines=[
            "how much light each layer absorbs and emits",
        ]),
    dict(
        kind="new",
        title="Radiation field",
        lines=[
            "intensity at each depth and frequency;",
            "in compression zones the coupling direction",
            "flips, so the missing neighbour is taken",
            "from the previous sweep (lagged store)",
        ]),
    dict(
        kind="stock",
        title="Rate equations",
        lines=[
            "residuals and Jacobian terms for the update:",
            "levels, electron density, temperature",
        ]),
    dict(
        kind="new",
        title="Settle the stored field",
        lines=[
            "repeat the radiation solve until the stored",
            "intensities stop drifting between rounds",
        ]),
    dict(
        kind="mod",
        title="Update the populations",
        lines=[
            "the mode controller picks one per call:",
            "  Lambda step   - safe fixed-point update",
            "                  (warm-up, drifting field)",
            "  Newton step   - full Jacobian per depth",
            "  Newton-Krylov - no Jacobian built (opt-in)",
        ]),
    dict(
        kind="mod",
        title="Acceleration (scheduled)",
        lines=[
            "extrapolate the last few iterates",
            "average levels that keep oscillating",
            "reuse the Jacobian while demand is small",
        ]),
    dict(
        kind="mod",
        title="Apply and checkpoint",
        lines=[
            "cap the step; roll back if the residual grows",
            "write restart files",
        ]),
]

DONE_TITLE = "Done"
DONE_LINES = [
    "populations have stopped changing;",
    "a separate tool computes the spectrum",
]

DEC_TEXT = "converged?"

LEGEND = [("new", "new code"),
          ("mod", "existing code, new switches"),
          ("stock", "unchanged CMFGEN")]


# ---------------------------------------------------------------
# Measure rendered text widths (pixels) on a throwaway figure, then
# convert to axis units.  This guarantees the boxes are wide enough.
# ---------------------------------------------------------------
def measure_px(items):
    """items: list of (text, fontsize, weight). Returns {text: px}."""
    figm = plt.figure(figsize=(8, 3), dpi=DPI)
    axm = figm.add_axes([0, 0, 1, 1])
    artists = [axm.text(0, 0, t, fontsize=fs, fontweight=w)
               for t, fs, w in items]
    figm.canvas.draw()
    out = {a.get_text(): a.get_window_extent().width for a in artists}
    plt.close(figm)
    return out


def box_height(lines):
    return PAD_V * 2 + LH_TITLE + LH_BODY * len(lines)


def draw_box(ax, xc, y_top, title, lines, kind, width):
    """Draw a rounded box hanging down from y_top; return its centre y."""
    h = box_height(lines)
    st = STYLE[kind]
    ax.add_patch(FancyBboxPatch(
        (xc - width / 2, y_top - h), width, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        fc=st["fc"], ec=st["ec"], lw=1.6, zorder=2))
    ax.text(xc, y_top - PAD_V - LH_TITLE / 2, title,
            ha="center", va="center", fontsize=FS_TITLE,
            fontweight="bold", zorder=3)
    y = y_top - PAD_V - LH_TITLE - LH_BODY / 2
    for ln in lines:
        ax.text(xc - width / 2 + TXT_PAD_L, y, ln,
                ha="left", va="center", fontsize=FS_BODY, zorder=3)
        y -= LH_BODY
    return y_top - h, h


def arrow(ax, x, y0, y1, label=None, lx=None, colour="#444444"):
    ax.add_patch(FancyArrowPatch(
        (x, y0), (x, y1), arrowstyle="-|>", mutation_scale=16,
        lw=1.4, color=colour, zorder=1))
    if label:
        ax.text(lx if lx is not None else x + 0.12, (y0 + y1) / 2,
                label, fontsize=8.6, color=colour, va="center")


# ---------------------------------------------------------------
# Measure every text element so columns can be sized to fit.
# ---------------------------------------------------------------
main_items = []
for b in BOXES:
    main_items.append((b["title"], FS_TITLE, "bold"))
    main_items += [(ln, FS_BODY, "normal") for ln in b["lines"]]
main_items += [(DONE_TITLE, FS_TITLE, "bold")]
main_items += [(ln, FS_BODY, "normal") for ln in DONE_LINES]
main_items += [(DEC_TEXT, FS_TITLE, "bold")]
leg_items = [(lab, 8.8, "normal") for _, lab in LEGEND]

px_main = measure_px(main_items)
px_leg = measure_px(leg_items)

# Convert pixel widths to axis units: with a fixed figure width, each
# axis unit is aw_px / xlim_span pixels wide.  The x-axis span must fit
# both the box column and the horizontal legend row.
aw_px = FIG_W_IN * DPI * AX_W_FRAC
wbox_px = max(px_main.values())
leg_text_px = sum(px_leg[lab] for _, lab in LEGEND)
leg_fixed = 3 * 0.46 + 2 * 0.5 + 0.6   # swatches, gaps, side margins

S_box = (BOX_X + TXT_PAD_L + TXT_PAD_R + RMARGIN) / (1.0 - wbox_px / aw_px)
S_leg = leg_fixed / (1.0 - leg_text_px / aw_px)
S = max(S_box, S_leg)

BOX_W = wbox_px * S / aw_px + TXT_PAD_L + TXT_PAD_R
XC = S / 2                 # centre the column on the canvas
XLIM = S

# ---------------------------------------------------------------
# Layout: stack the main cycle down the left column.
# ---------------------------------------------------------------
total_h = sum(box_height(b["lines"]) for b in BOXES)
dec_h, done_h = 0.85, box_height(DONE_LINES)
# +1.0 headroom for the legend row under the title
fig_h = total_h + GAP * (len(BOXES) + 1) + dec_h + done_h + 1.9

fig, ax = plt.subplots(figsize=(FIG_W_IN, fig_h * 0.58))
ax.set_xlim(0, XLIM)
ax.set_ylim(0, fig_h)
ax.axis("off")

y = fig_h - 0.55
ax.text(XLIM / 2, y + 0.25,
        "The iteration loop for non-monotonic models",
        ha="center", fontsize=13, fontweight="bold")

# Horizontal legend row under the title.
leg_y = y - 0.15
units_per_px = XLIM / aw_px
leg_w = sum(0.34 + 0.12 + px_leg[lab] * units_per_px
            for _, lab in LEGEND) + 0.5 * (len(LEGEND) - 1)
lx = XLIM / 2 - leg_w / 2
for kind, lab in LEGEND:
    st = STYLE[kind]
    ax.add_patch(FancyBboxPatch(
        (lx, leg_y - 0.14), 0.34, 0.28,
        boxstyle="round,pad=0.01,rounding_size=0.05",
        fc=st["fc"], ec=st["ec"], lw=1.4))
    lx += 0.46
    ax.text(lx, leg_y, lab, fontsize=8.8, va="center")
    lx += px_leg[lab] * units_per_px + 0.5
ax.text(XLIM / 2, leg_y - 0.4,
        "new features are off by default",
        fontsize=8.4, ha="center", va="top", style="italic",
        color="#444444")
y -= 1.0

bounds = []
for b in BOXES:
    yb, h = draw_box(ax, XC, y, b["title"], b["lines"], b["kind"], BOX_W)
    bounds.append((yb, h))
    if len(bounds) > 1:
        arrow(ax, XC, bounds[-2][0], y)
    y = yb - GAP

# Decision box.
dec_top = y
dec_w = max(0.5 * BOX_W, px_main[DEC_TEXT] * S / aw_px + 0.7)
ax.add_patch(FancyBboxPatch(
    (XC - dec_w / 2, dec_top - dec_h), dec_w, dec_h,
    boxstyle="round,pad=0.02,rounding_size=0.12",
    fc="#f2f2f2", ec="#777777", lw=1.6, zorder=2))
ax.text(XC, dec_top - dec_h / 2, DEC_TEXT,
        ha="center", va="center", fontsize=FS_TITLE, fontweight="bold")
arrow(ax, XC, bounds[-1][0], dec_top)

# Done box below the decision.
done_top = dec_top - dec_h - GAP
yb, _ = draw_box(ax, XC, done_top, DONE_TITLE, DONE_LINES,
                 "stock", BOX_W)
arrow(ax, XC, dec_top - dec_h, done_top, label="yes", lx=XC + 0.15)

# Loop-back arrow: decision "no" -> left margin -> back to box 1.
x_loop = 0.28
box1_mid = bounds[0][0] - bounds[0][1] / 2
ax.add_patch(FancyArrowPatch(
    (XC - dec_w / 2, dec_top - dec_h / 2), (x_loop, dec_top - dec_h / 2),
    arrowstyle="-", lw=1.4, color="#444444"))
ax.add_patch(FancyArrowPatch(
    (x_loop, dec_top - dec_h / 2), (x_loop, box1_mid),
    arrowstyle="-", lw=1.4, color="#444444"))
ax.add_patch(FancyArrowPatch(
    (x_loop, box1_mid), (XC - BOX_W / 2, box1_mid),
    arrowstyle="-|>", mutation_scale=16, lw=1.4, color="#444444"))
ax.text(XC - dec_w / 2 - 0.15, dec_top - dec_h / 2 + 0.16,
        "no: repeat", fontsize=8.6, color="#444444", ha="right")

out = Path(__file__).resolve().parent / "iteration_cycle.png"
plt.savefig(out, dpi=DPI, bbox_inches="tight")
print(f"wrote {out}  (BOX_W={BOX_W:.2f}, XLIM={XLIM:.2f})")
