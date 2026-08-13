#!/usr/bin/env python3
"""Cartoon documenting the ISI experimental design.

Two panels. Left: a flower is pollinated with its own pollen (S = selfed) and
produces no seeds because self-incompatibility rejects the self pollen.
Right: the same flower is pollinated with pollen from a different plant
(O = outcrossed) and produces seeds. Companion illustration for
analyze_ISI_breeding_system.py; use as an introduction slide before
showing the ISI density figure.

Flowers are drawn with 4 white petals arranged as a cross (Brassicaceae
convention — Lepidium papilliferum has a cruciform corolla).

Output:
    figures/ISI_experimental_design_cartoon.pdf|.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Rectangle

DEFAULT_FIGURES = "figures"

# ------------------------------------------------------- palette
COL_PETAL = "#FFFFFF"        # LEPA petals are white
COL_PETAL_EDGE = "#B8B8B8"   # grey outline so white shows against white bg
COL_CENTER = "#F4C430"       # golden yellow (anthers + stigma)
COL_ANTHER_DOT = "#6B4E1E"
COL_STEM = "#4CAF50"
COL_LEAF = "#66BB6A"
COL_POD = "#F5EFE0"          # translucent silique wall (cream)
COL_POD_EDGE = "#7D5A2B"
COL_SEED = "#A0522D"         # sienna — realistic LEPA seed colour
COL_SEED_EDGE = "#5D3A1A"
COL_NO = "#B0B0B0"
COL_SI = "#b2182b"           # ISI SI band red — SI outcome only
COL_SC = "#1b7837"           # ISI SC band green — SC outcome only
COL_NEUTRAL = "#2C5F8D"      # medium blue — distinct from red/green; outcross control
COL_HEADER = "#2c3e50"       # dark slate — panel headers and treatment badges
COL_ARROW_SELF = "#5D4E7A"   # purple
COL_ARROW_CROSS = COL_NEUTRAL

FONT_HEADER = 15
FONT_LABEL = 12
FONT_MICRO = 10

# ------------------------------------------------------- primitives
def draw_flower(ax, cx, cy, radius=0.55, stem_bottom=None):
    """4-petal white Brassicaceae flower with optional stem."""
    if stem_bottom is not None:
        ax.plot([cx, cx], [stem_bottom, cy - radius * 0.55],
                color=COL_STEM, linewidth=4.5, solid_capstyle="round", zorder=1)
        leaf = Ellipse((cx - radius * 0.75, (cy + stem_bottom) / 2 + 0.05),
                       radius * 1.0, radius * 0.32,
                       angle=30, facecolor=COL_LEAF, edgecolor=COL_STEM,
                       linewidth=1, zorder=2)
        ax.add_patch(leaf)
    for angle in (0, 90, 180, 270):                     # cruciform 4 petals
        rad = np.deg2rad(angle)
        px = cx + radius * 0.55 * np.cos(rad)
        py = cy + radius * 0.55 * np.sin(rad)
        petal = Ellipse((px, py), radius * 0.85, radius * 0.55,
                        angle=angle, facecolor=COL_PETAL,
                        edgecolor=COL_PETAL_EDGE, linewidth=1.5, zorder=3)
        ax.add_patch(petal)
    center = Circle((cx, cy), radius * 0.24, facecolor=COL_CENTER,
                    edgecolor="#B8930A", linewidth=1.2, zorder=4)
    ax.add_patch(center)
    for a in (45, 135, 225, 315):              # 4 stamens (Lepidium reduction)
        r = np.deg2rad(a)
        dx = cx + radius * 0.13 * np.cos(r)
        dy = cy + radius * 0.13 * np.sin(r)
        ax.add_patch(Circle((dx, dy), radius * 0.05,
                            facecolor=COL_ANTHER_DOT, edgecolor="none", zorder=5))

def draw_silique_with_seeds(ax, cx, cy, size=0.5):
    """Silique (LEPA fruit) shown as translucent pod with visible oval seeds
    arranged in two rows separated by a septum — anatomically correct for
    Brassicaceae siliques."""
    pod = Ellipse((cx, cy), size * 1.3, size * 2.1,
                  facecolor=COL_POD, edgecolor=COL_POD_EDGE,
                  linewidth=2, zorder=2)
    ax.add_patch(pod)
    ax.plot([cx, cx], [cy - size * 0.90, cy + size * 0.90],
            color=COL_POD_EDGE, linewidth=1, linestyle=(0, (2, 2)), zorder=2)
    for py in (0.65, 0.30, -0.05, -0.40, -0.72):
        for dx in (-0.22, 0.22):
            seed = Ellipse((cx + dx * size, cy + py * size),
                           size * 0.28, size * 0.20,
                           facecolor=COL_SEED, edgecolor=COL_SEED_EDGE,
                           linewidth=1, zorder=3)
            ax.add_patch(seed)

def draw_empty_silique(ax, cx, cy, size=0.5):
    """Dashed empty outline of what the silique would look like, with a red ×."""
    pod = Ellipse((cx, cy), size * 1.3, size * 2.1,
                  facecolor="none", edgecolor=COL_NO, linewidth=2,
                  linestyle=(0, (5, 3)), zorder=2)
    ax.add_patch(pod)
    off_x, off_y = size * 0.55, size * 0.9
    ax.plot([cx - off_x, cx + off_x], [cy - off_y, cy + off_y],
            color=COL_SI, linewidth=4, solid_capstyle="round", zorder=3)
    ax.plot([cx - off_x, cx + off_x], [cy + off_y, cy - off_y],
            color=COL_SI, linewidth=4, solid_capstyle="round", zorder=3)

def draw_self_loop(ax, cx, cy, radius):
    """Loop arrow bulging UP: pollen leaves anther, arcs over, returns to stigma.
    Endpoints sit at anther level; peak clears the label above."""
    end_y = cy + radius * 0.4                 # endpoints at anther/stigma level
    start = (cx + radius * 0.55, end_y)
    end   = (cx - radius * 0.55, end_y)
    arr = FancyArrowPatch(start, end, connectionstyle="arc3,rad=1.3",
                          arrowstyle="->", mutation_scale=22,
                          color=COL_ARROW_SELF, linewidth=2.2, zorder=6)
    ax.add_patch(arr)
    ax.text(cx, cy + radius * 1.4, "self pollen",
            ha="center", va="bottom", fontsize=FONT_MICRO,
            color=COL_ARROW_SELF, fontstyle="italic", fontweight="bold")

def draw_cross_arrow(ax, donor, recip, radius):
    """Arrow showing pollen going from donor flower to recipient flower."""
    dcx, dcy = donor
    rcx, rcy = recip
    arr = FancyArrowPatch(
        (dcx + radius * 0.35, dcy + radius * 0.4),
        (rcx - radius * 0.35, rcy + radius * 0.4),
        connectionstyle="arc3,rad=-0.45",
        arrowstyle="->", mutation_scale=20,
        color=COL_ARROW_CROSS, linewidth=2.0, zorder=6,
    )
    ax.add_patch(arr)
    ax.text((dcx + rcx) / 2, max(dcy, rcy) + radius * 1.7, "cross pollen",
            ha="center", va="bottom", fontsize=FONT_MICRO,
            color=COL_ARROW_CROSS, fontstyle="italic")

def draw_down_arrow(ax, cx, y_top, y_bot):
    arr = FancyArrowPatch((cx, y_top), (cx, y_bot), arrowstyle="->",
                          mutation_scale=17, color="#666666", linewidth=1.6)
    ax.add_patch(arr)

def draw_treatment_badge(ax, cx, cy, letter, fill_color):
    """Small coloured circle badge with S or O letter — matches ISI figure notation."""
    ax.add_patch(Circle((cx, cy), 0.18, facecolor=fill_color,
                        edgecolor="black", linewidth=1.5, zorder=8))
    ax.text(cx, cy, letter, ha="center", va="center", fontsize=13,
            fontweight="bold", color="white", zorder=9)

# ------------------------------------------------------- panels
def build_panel_selfing(ax):
    """Self-pollination panel: one flower, TWO possible outcomes (SI vs SC)
    depending on the plant's breeding system."""
    ax.set_xlim(0, 5); ax.set_ylim(0, 6); ax.set_aspect("equal"); ax.axis("off")
    # No overall panel tint — outcome boxes carry the ISI colours

    ax.text(2.5, 5.75, "S  —  Self-pollinated", ha="center", va="top",
            fontsize=FONT_HEADER, fontweight="bold", color=COL_HEADER)

    r = 0.42
    fx, fy = 2.5, 4.6
    draw_flower(ax, fx, fy, radius=r, stem_bottom=3.85)
    draw_self_loop(ax, fx, fy, r)
    draw_treatment_badge(ax, fx + r * 1.30, fy + r * 0.95, "S", COL_HEADER)
    ax.text(fx, 3.55, "Plant A", ha="center", va="top",
            fontsize=FONT_MICRO, fontweight="bold", color="#444444")
    ax.text(fx, 3.32, "(same plant  —  source and recipient)",
            ha="center", va="top",
            fontsize=FONT_MICRO - 1, color="#666666", fontstyle="italic")

    # Branching arrow: single stem down, forks to two outcomes
    ax.plot([fx, fx], [3.10, 2.85], color="#666666", linewidth=1.6, zorder=1)
    for target_x in (1.20, 3.80):
        arr = FancyArrowPatch((fx, 2.85), (target_x, 2.35),
                              arrowstyle="->", mutation_scale=16,
                              color="#666666", linewidth=1.6, zorder=1)
        ax.add_patch(arr)

    # Left outcome — SI (red tint box)
    ax.add_patch(Rectangle((-0.05, 0.10), 2.50, 2.05, facecolor=COL_SI, alpha=0.10,
                           edgecolor="none", zorder=-1))
    draw_empty_silique(ax, 1.20, 1.35, size=0.36)
    ax.text(1.20, 0.55, "SI", ha="center", va="bottom",
            fontsize=FONT_LABEL + 1, fontweight="bold", color=COL_SI)
    ax.text(1.20, 0.22, "self-incompatible  →  no seeds",
            ha="center", va="bottom", fontsize=FONT_MICRO - 1,
            color="#555555", fontstyle="italic")

    # Right outcome — SC (green tint box)
    ax.add_patch(Rectangle((2.55, 0.10), 2.50, 2.05, facecolor=COL_SC, alpha=0.10,
                           edgecolor="none", zorder=-1))
    draw_silique_with_seeds(ax, 3.80, 1.35, size=0.36)
    ax.text(3.80, 0.55, "SC", ha="center", va="bottom",
            fontsize=FONT_LABEL + 1, fontweight="bold", color=COL_SC)
    ax.text(3.80, 0.22, "self-compatible  →  seeds",
            ha="center", va="bottom", fontsize=FONT_MICRO - 1,
            color="#555555", fontstyle="italic")

def build_panel_outcross(ax):
    """Cross-pollination panel: control condition, produces seeds regardless
    of breeding system. Colour is neutral (blue-slate) — deliberately not
    reusing the ISI SC/SI colours."""
    ax.set_xlim(0, 5); ax.set_ylim(0, 6); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 5, 6, facecolor=COL_NEUTRAL, alpha=0.10,
                           edgecolor="none", zorder=-1))

    ax.text(2.5, 5.75, "O  —  Cross-pollinated", ha="center", va="top",
            fontsize=FONT_HEADER, fontweight="bold", color=COL_HEADER)

    r = 0.42
    donor = (1.35, 4.6); recip = (3.65, 4.6)
    draw_flower(ax, *donor, radius=r, stem_bottom=3.85)
    draw_flower(ax, *recip, radius=r, stem_bottom=3.85)
    draw_cross_arrow(ax, donor, recip, r)
    draw_treatment_badge(ax, recip[0] + r * 1.30, recip[1] + r * 0.95, "O", COL_HEADER)

    ax.text(donor[0], 3.55, "Plant A", ha="center", va="top",
            fontsize=FONT_MICRO, fontweight="bold", color="#444444")
    ax.text(donor[0], 3.32, "(pollen donor)", ha="center", va="top",
            fontsize=FONT_MICRO - 1, color="#666666", fontstyle="italic")
    ax.text(recip[0], 3.55, "Plant B", ha="center", va="top",
            fontsize=FONT_MICRO, fontweight="bold", color="#444444")
    ax.text(recip[0], 3.32, "(recipient)", ha="center", va="top",
            fontsize=FONT_MICRO - 1, color="#666666", fontstyle="italic")

    draw_down_arrow(ax, recip[0], 3.05, 2.55)
    draw_silique_with_seeds(ax, recip[0], 1.55, size=0.42)

    ax.text(2.5, 0.55, "Seeds produced  —  control",
            ha="center", va="bottom", fontsize=FONT_LABEL, fontweight="bold",
            color=COL_NEUTRAL)
    ax.text(2.5, 0.22, "outcrossing succeeds for both SI and SC plants",
            ha="center", va="bottom", fontsize=FONT_MICRO - 1,
            color="#666666", fontstyle="italic")

# ------------------------------------------------------- main
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures-dir", type=Path, default=Path(DEFAULT_FIGURES))
    args = parser.parse_args()
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax_s, ax_o) = plt.subplots(1, 2, figsize=(12, 6),
                                     gridspec_kw={"width_ratios": [1, 1]})
    build_panel_selfing(ax_s)
    build_panel_outcross(ax_o)
    fig.suptitle(
        "The ISI experiment — comparing selfed vs cross-pollinated fruit set",
        fontsize=15, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_pdf = args.figures_dir / "ISI_experimental_design_cartoon.pdf"
    out_png = args.figures_dir / "ISI_experimental_design_cartoon.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"[cartoon] {out_pdf}")
    print(f"[cartoon] {out_png}")

if __name__ == "__main__":
    main()
