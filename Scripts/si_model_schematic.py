#!/usr/bin/env python3
"""Pedagogical schematic of sporophytic SI with Class I / Class II
dominance in tetraploid LEPA (Phase 5 § A.6).

Three panels:
  A. Dominance within one plant — Case A (Class I dominant) vs Case B
     (all Class II co-dominant), showing which alleles get expressed
     on pollen + stigma.
  B. Between-plant recognition — worked example: mother M vs three
     candidate fathers F1 / F2 / F3, with the verdict per cross.
  C. Verdict card — the three cross types (Class I × Class I,
     Class I × Class II, Class II × Class II) and the compatibility
     rule per type.

Output
------
    figures/Phase5/step30_A_si_model_schematic.png/pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

OUT_DIR = Path("figures/Phase5")

# Okabe-Ito palette
CLASS_I     = "#D55E00"   # vermillion — Class I expressed
CLASS_II    = "#6B6B6B"   # dark gray — Class II expressed
SILENT      = "#E5E5E5"   # very faded — silent alleles
COMPATIBLE  = "#009E73"   # green
REJECTED    = "#CC3333"   # red
INK         = "#222222"
BG_MOTHER   = "#FBEEDF"   # light peach mother panel background
BG_FATHER   = "#F1F1F1"   # light gray father panel background


def allele_tile(ax, x, y, label, class_type, silent=False,
                w=0.10, h=0.09, fontsize=8.5):
    """Draw a single allele tile at (x, y) as (left, bottom).
    class_type: 'I' or 'II'; silent=True dims it (Case-A Class-II slots).
    """
    if silent:
        fc, ec, txt = SILENT, "#CCCCCC", "#999999"
    elif class_type == "I":
        fc, ec, txt = CLASS_I, CLASS_I, "white"
    else:
        fc, ec, txt = CLASS_II, CLASS_II, "white"
    ax.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.01",
        facecolor=fc, edgecolor=ec, linewidth=0.7))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=fontsize, color=txt, weight="bold")


def draw_arrow_down(ax, x, y_start, y_end, colour="#666", lw=1.4):
    ax.add_patch(FancyArrowPatch((x, y_start), (x, y_end),
                                 arrowstyle="->", mutation_scale=16,
                                 color=colour, lw=lw,
                                 shrinkA=0, shrinkB=0))


def draw_genotype_row(ax, alleles, classes, x_start, y, tile_w=0.10,
                       tile_h=0.09, gap=0.02, silent_mask=None,
                       fontsize=8.5):
    """Draw a row of allele tiles for a plant's genotype."""
    if silent_mask is None:
        silent_mask = [False] * len(alleles)
    for i, (a, c, sil) in enumerate(zip(alleles, classes, silent_mask)):
        allele_tile(ax, x_start + i * (tile_w + gap), y, a, c,
                    silent=sil, w=tile_w, h=tile_h, fontsize=fontsize)


# ---------------------------------------------------------------------------
# Panel A — Dominance within one plant
# ---------------------------------------------------------------------------
def panel_A(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.96,
            "A. Dominance within one plant — which alleles show up on "
            "pollen and stigma?",
            ha="center", va="top", fontsize=13, weight="bold", color=INK)

    # ---- Case A (left) ----
    axA_x0 = 0.05
    ax.text(axA_x0 + 0.20, 0.85,
            "Case A — plant carries ≥ 1 Class I allele",
            ha="center", fontsize=11, weight="bold", color=INK)
    # Genotype row
    ax.text(axA_x0 - 0.02, 0.71, "Genotype (4 SRK\ncopies, tetraploid):",
            ha="right", va="center", fontsize=9, color=INK)
    draw_genotype_row(ax,
        alleles=["FG001", "FG002", "FG010", "FG015"],
        classes=["I", "I", "II", "II"],
        x_start=axA_x0, y=0.68)
    # Arrow
    draw_arrow_down(ax, axA_x0 + 0.20, 0.66, 0.53)
    # Expressed row
    ax.text(axA_x0 - 0.02, 0.44, "Expressed on\npollen + stigma:",
            ha="right", va="center", fontsize=9, color=INK)
    draw_genotype_row(ax,
        alleles=["FG001", "FG002", "FG010", "FG015"],
        classes=["I", "I", "II", "II"],
        x_start=axA_x0, y=0.41,
        silent_mask=[False, False, True, True])
    ax.text(axA_x0 + 0.20, 0.26,
            "Class I dominant → only FG001 & FG002 "
            "participate in SI.\nFG010 and FG015 are present in the "
            "genome but silent for SI.",
            ha="center", va="center", fontsize=9, color=INK, style="italic")

    # ---- Case B (right) — no left labels (redundant with Case A) ----
    axB_x0 = 0.62
    tile_w_B = 0.085; gap_B = 0.012   # slightly tighter to fit in [0, 1]
    ax.text(axB_x0 + 0.18, 0.85,
            "Case B — plant has only Class II alleles",
            ha="center", fontsize=11, weight="bold", color=INK)
    draw_genotype_row(ax,
        alleles=["FG010", "FG015", "FG018", "FG024"],
        classes=["II", "II", "II", "II"],
        x_start=axB_x0, y=0.68, tile_w=tile_w_B, gap=gap_B, fontsize=8)
    draw_arrow_down(ax, axB_x0 + 0.18, 0.66, 0.53)
    draw_genotype_row(ax,
        alleles=["FG010", "FG015", "FG018", "FG024"],
        classes=["II", "II", "II", "II"],
        x_start=axB_x0, y=0.41, tile_w=tile_w_B, gap=gap_B, fontsize=8)
    ax.text(axB_x0 + 0.18, 0.26,
            "No Class I present → all four Class II alleles are\n"
            "expressed co-dominantly and participate in SI.",
            ha="center", va="center", fontsize=9, color=INK, style="italic")

    # Legend at the bottom
    lg_y = 0.08
    for i, (colour, label) in enumerate([
        (CLASS_I,  "Class I expressed"),
        (CLASS_II, "Class II expressed"),
        (SILENT,   "silent (present in genome, not expressed for SI)"),
    ]):
        cx = 0.10 + i * 0.30
        allele_tile(ax, cx, lg_y, "", "I" if colour == CLASS_I else "II",
                    silent=(colour == SILENT), w=0.04, h=0.05)
        ax.text(cx + 0.05, lg_y + 0.025, label,
                ha="left", va="center", fontsize=9, color=INK)


# ---------------------------------------------------------------------------
# Panel B — Worked example: M vs F1 / F2 / F3
# ---------------------------------------------------------------------------
def panel_B(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.96,
            "B. Between-plant recognition — worked example, one mother "
            "vs three candidate fathers",
            ha="center", va="top", fontsize=13, weight="bold", color=INK)

    # ---- Mother (top-center) ----
    mom_x0 = 0.35
    mom_tiles_y = 0.72
    # Background panel — encompass title (top) + tiles + expressed set
    ax.add_patch(mpatches.FancyBboxPatch(
        (mom_x0 - 0.06, 0.66), 0.62, 0.22,
        boxstyle="round,pad=0.005,rounding_size=0.015",
        facecolor=BG_MOTHER, edgecolor="#F0B85F", linewidth=1.0))
    # Title above tiles
    ax.text(mom_x0 + 0.20, 0.855,
            "MOTHER  M",
            ha="center", va="center", fontsize=11, weight="bold", color=INK)
    ax.text(mom_x0 - 0.07, mom_tiles_y + 0.045,
            "Genotype:", ha="right", va="center",
            fontsize=9, color=INK)
    draw_genotype_row(ax,
        alleles=["FG001", "FG002", "FG024", "FG031"],
        classes=["I", "I", "II", "II"],
        x_start=mom_x0, y=mom_tiles_y)
    ax.text(mom_x0 + 0.20, 0.685,
            "Expressed set (Case A → Class I only):  {FG001, FG002}",
            ha="center", va="center",
            fontsize=9, color=CLASS_I, weight="bold")

    # ---- Fathers ----
    father_data = [
        {
            "name": "F1", "x0": 0.02,
            "alleles": ["FG001", "FG007", "FG024", "FG032"],
            "classes": ["I", "II", "II", "II"],
            "expr":    "{FG001}",
            "expr_col": CLASS_I,
            "case":    "Case A (Class I present)",
            "verdict": "REJECTED",
            "reason":  "FG001 is in M's\nexpressed set",
            "colour":  REJECTED,
        },
        {
            "name": "F2", "x0": 0.35,
            "alleles": ["FG015", "FG018", "FG024", "FG032"],
            "classes": ["II", "II", "II", "II"],
            "expr":    "{FG015, FG018, FG024, FG032}",
            "expr_col": CLASS_II,
            "case":    "Case B (all Class II)",
            "verdict": "COMPATIBLE",
            "reason":  "no overlap with\nM's Class I set",
            "colour":  COMPATIBLE,
        },
        {
            "name": "F3", "x0": 0.68,
            "alleles": ["FG002", "FG003", "FG010", "FG018"],
            "classes": ["I", "I", "II", "II"],
            "expr":    "{FG002, FG003}",
            "expr_col": CLASS_I,
            "case":    "Case A (Class I present)",
            "verdict": "REJECTED",
            "reason":  "FG002 is in M's\nexpressed set",
            "colour":  REJECTED,
        },
    ]

    for f in father_data:
        x0 = f["x0"]
        # Panel background
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0 - 0.005, 0.02), 0.32, 0.55,
            boxstyle="round,pad=0.005,rounding_size=0.01",
            facecolor=BG_FATHER, edgecolor="#BBBBBB", linewidth=0.7))
        # Arrow from mother down to this father, with verdict alongside
        arrow_start_x = mom_x0 + 0.20
        arrow_start_y = 0.665           # bottom of mother panel
        arrow_end_x = x0 + 0.155
        arrow_end_y = 0.57
        ax.add_patch(FancyArrowPatch(
            (arrow_start_x, arrow_start_y), (arrow_end_x, arrow_end_y),
            arrowstyle="->", mutation_scale=18,
            color=f["colour"], lw=1.6, shrinkA=1, shrinkB=1))
        # Father header
        ax.text(x0 + 0.155, 0.51,
                f'FATHER  {f["name"]}',
                ha="center", fontsize=10, weight="bold", color=INK)
        # Genotype
        ax.text(x0 + 0.155, 0.46,
                "Genotype:", ha="center", fontsize=8, color=INK)
        draw_genotype_row(ax, alleles=f["alleles"], classes=f["classes"],
                          x_start=x0 + 0.01, y=0.37,
                          tile_w=0.072, gap=0.005, fontsize=7.5)
        # Case
        ax.text(x0 + 0.155, 0.31,
                f["case"], ha="center", fontsize=8, style="italic",
                color=INK)
        # Expressed set
        ax.text(x0 + 0.155, 0.26,
                "Expressed:", ha="center", fontsize=8, color=INK)
        ax.text(x0 + 0.155, 0.22,
                f["expr"], ha="center", fontsize=8.5,
                color=f["expr_col"], weight="bold")
        # Verdict badge
        badge_col = f["colour"]
        ax.add_patch(mpatches.FancyBboxPatch(
            (x0 + 0.045, 0.09), 0.22, 0.08,
            boxstyle="round,pad=0.005,rounding_size=0.015",
            facecolor=badge_col, edgecolor=badge_col, linewidth=0))
        verdict_symbol = "✓" if f["verdict"] == "COMPATIBLE" else "✗"
        ax.text(x0 + 0.155, 0.13,
                f"{verdict_symbol}  {f['verdict']}",
                ha="center", va="center",
                fontsize=10, color="white", weight="bold")
        # Reason
        ax.text(x0 + 0.155, 0.045,
                f["reason"], ha="center", va="center",
                fontsize=8, color=INK, style="italic")


# ---------------------------------------------------------------------------
# Panel C — Verdict card: the three cross types
# ---------------------------------------------------------------------------
def panel_C(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.96,
            "C. Compatibility rule by cross type",
            ha="center", va="top", fontsize=13, weight="bold", color=INK)

    row_specs = [
        {
            "y_center": 0.72,
            "title":    "Class I × Class I",
            "mom":      (["FG001", "FG002", "FG010", "FG015"],
                          ["I", "I", "II", "II"],
                          [False, False, True, True]),
            "mom_expr": "{FG001, FG002}",
            "dad":      (["FG003", "FG004", "FG010", "FG015"],
                          ["I", "I", "II", "II"],
                          [False, False, True, True]),
            "dad_expr": "{FG003, FG004}",
            "rule":     "Compatible if their **expressed Class I alleles differ**.\n"
                        "Class II alleles on both sides are silent — they can be\n"
                        "shared without affecting SI. Even a Class II match is\n"
                        "irrelevant here.",
            "verdict":  "COMPATIBLE",
            "verdict_col": COMPATIBLE,
            "verdict_reason": "no Class I overlap",
        },
        {
            "y_center": 0.44,
            "title":    "Class I × Class II  (between-class)",
            "mom":      (["FG001", "FG002", "FG010", "FG015"],
                          ["I", "I", "II", "II"],
                          [False, False, True, True]),
            "mom_expr": "{FG001, FG002}",
            "dad":      (["FG010", "FG015", "FG018", "FG024"],
                          ["II", "II", "II", "II"],
                          [False, False, False, False]),
            "dad_expr": "{FG010, FG015, FG018, FG024}",
            "rule":     "**Always compatible.** One side expresses only Class I,\n"
                        "the other side expresses only Class II — the two\n"
                        "expressed sets belong to disjoint classes by construction,\n"
                        "so no allele-level check is needed.",
            "verdict":  "ALWAYS COMPATIBLE",
            "verdict_col": COMPATIBLE,
            "verdict_reason": "disjoint classes",
        },
        {
            "y_center": 0.16,
            "title":    "Class II × Class II",
            "mom":      (["FG010", "FG015", "FG018", "FG024"],
                          ["II", "II", "II", "II"],
                          [False, False, False, False]),
            "mom_expr": "{FG010, FG015, FG018, FG024}",
            "dad":      (["FG010", "FG020", "FG028", "FG031"],
                          ["II", "II", "II", "II"],
                          [False, False, False, False]),
            "dad_expr": "{FG010, FG020, FG028, FG031}",
            "rule":     "Compatible if their four expressed Class II alleles do\n"
                        "**not share any allele**. All four are expressed on both\n"
                        "sides, so the whole genotype matters.",
            "verdict":  "REJECTED",
            "verdict_col": REJECTED,
            "verdict_reason": "shared FG010",
        },
    ]

    # Small tile geometry for Panel C so labels don't collide
    TW, TH, TG = 0.048, 0.045, 0.004
    genotype_width = 4 * TW + 3 * TG   # ~0.204

    for r in row_specs:
        yc = r["y_center"]
        # Row background (taller — 0.26)
        ax.add_patch(mpatches.Rectangle(
            (0.01, yc - 0.13), 0.98, 0.26,
            facecolor="#FAFAFA", edgecolor="#DDDDDD", linewidth=0.7))
        # Row title (above content, inside row bg)
        ax.text(0.02, yc + 0.11,
                r["title"], ha="left", va="center",
                fontsize=11, weight="bold", color=INK)

        # Column 1 — mother / × / father schematic ----------------------------
        mom_x = 0.03
        dad_x = mom_x + genotype_width + 0.05     # cross symbol in the gap
        cross_x = mom_x + genotype_width + 0.025

        # MOTHER label above the tile row
        ax.text(mom_x + genotype_width / 2, yc + 0.055,
                "MOTHER", ha="center", va="center",
                fontsize=8, color=INK, weight="bold")
        # "expressed:" prefix at the far-left of the expressed-set row
        ax.text(mom_x - 0.005, yc - 0.05, "expressed:",
                ha="right", va="center", fontsize=7.5,
                color="#555", style="italic")
        draw_genotype_row(ax, alleles=r["mom"][0], classes=r["mom"][1],
                          silent_mask=r["mom"][2],
                          x_start=mom_x, y=yc,
                          tile_w=TW, tile_h=TH, gap=TG, fontsize=7)
        # Expressed set label under the tiles (short — no "expressed:" prefix)
        ax.text(mom_x + genotype_width / 2, yc - 0.05,
                r["mom_expr"],
                ha="center", va="center", fontsize=7.5,
                color=CLASS_I if any(c == "I" and not s for c, s
                                     in zip(r["mom"][1], r["mom"][2]))
                                     else CLASS_II,
                weight="bold")

        # Cross symbol between mother and father
        ax.text(cross_x, yc + TH / 2, "×", ha="center", va="center",
                fontsize=20, color=INK, weight="bold")

        # FATHER label above the tile row
        ax.text(dad_x + genotype_width / 2, yc + 0.055,
                "FATHER", ha="center", va="center",
                fontsize=8, color=INK, weight="bold")
        draw_genotype_row(ax, alleles=r["dad"][0], classes=r["dad"][1],
                          silent_mask=r["dad"][2],
                          x_start=dad_x, y=yc,
                          tile_w=TW, tile_h=TH, gap=TG, fontsize=7)
        ax.text(dad_x + genotype_width / 2, yc - 0.05,
                r["dad_expr"],
                ha="center", va="center", fontsize=7.5,
                color=CLASS_I if any(c == "I" and not s for c, s
                                     in zip(r["dad"][1], r["dad"][2]))
                                     else CLASS_II,
                weight="bold")

        # Column 2 — rule text -----------------------------------------------
        rule_x = dad_x + genotype_width + 0.035
        ax.text(rule_x, yc + 0.075, r["rule"].replace("**", ""),
                ha="left", va="top", fontsize=9, color=INK, wrap=True)

        # Column 3 — verdict badge -------------------------------------------
        bx, by, bw, bh = 0.86, yc - 0.045, 0.13, 0.09
        ax.add_patch(mpatches.FancyBboxPatch(
            (bx, by), bw, bh,
            boxstyle="round,pad=0.005,rounding_size=0.012",
            facecolor=r["verdict_col"], edgecolor=r["verdict_col"]))
        verdict_symbol = "✓" if r["verdict"].startswith(
            ("COMP", "ALWAYS")) else "✗"
        ax.text(bx + bw / 2, by + bh / 2,
                f"{verdict_symbol}  {r['verdict']}",
                ha="center", va="center",
                fontsize=8, color="white", weight="bold")
        # Reason above the badge
        ax.text(bx + bw / 2, by + bh + 0.02,
                r["verdict_reason"], ha="center", va="center",
                fontsize=8, style="italic", color=INK)


# ---------------------------------------------------------------------------
# Assemble figure
# ---------------------------------------------------------------------------
def main() -> None:
    fig = plt.figure(figsize=(13.5, 15.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.35, 1.35],
                          hspace=0.10)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])

    panel_A(axA)
    panel_B(axB)
    panel_C(axC)

    fig.suptitle(
        "Sporophytic self-incompatibility with Class I / Class II dominance "
        "in tetraploid LEPA  (Phase 5 § A.6)",
        fontsize=14, weight="bold", color=INK, y=0.995,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "step30_A_si_model_schematic.png",
                dpi=200, bbox_inches="tight")
    fig.savefig(OUT_DIR / "step30_A_si_model_schematic.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"[si_model_schematic] Figure in {OUT_DIR}/")


if __name__ == "__main__":
    main()
