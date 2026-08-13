#!/usr/bin/env python3
"""ISI (Index of Self-Incompatibility) analysis for Lepidium papilliferum from
the Billinge & Robertson outcrossing dataset. Bootstraps

    ISI = 1 - mean(S fruit set) / mean(O fruit set)

where S = selfed and O = outcrossed. Runs at both species-wide and per-EO
levels and classifies breeding system:

    ISI  < 0.2         -> SC       (self-compatible)
    0.2 <= ISI < 0.8   -> Partial  (partially self-compatible)
    ISI >= 0.8         -> SI       (self-incompatible)

Selfing (S) and outcrossing (O) samples are resampled independently with
replacement. Zero values in percent fruit set are replaced with a small
epsilon (1e-6) to prevent divide-by-zero in bootstrap resamples that draw
all-zero O. Bootstrapped ISI is then clipped at 0 (negative ISI means
selfed > outcrossed, biologically uninformative for classification).

By default, analysis is restricted to sites (Site Acronym) that carry BOTH
S and O observations, so mean(S) and mean(O) are computed on the same set
of sites (apples-to-apples). Pass --include-unmatched to pool across all
sites (adds S-only sites and inflates ISI upward when those sites are low).

Raw XLS uses "AS" for outcrossed observations; the loader renames to "O".

Outputs:
    Tables/ISI_breeding_system_summary.tsv    per-level summary + bootstrap CI
                                              (species-wide and per-EO rows kept for reference)
    figures/ISI_breeding_system.pdf|.png      species-wide bootstrap density
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

try:
    from draw_ISI_cartoon import (
        draw_silique_with_seeds, draw_empty_silique,
        draw_flower, draw_self_loop,
    )
    from matplotlib.patches import FancyArrowPatch
    _CARTOON_GLYPHS_AVAILABLE = True
except ImportError:
    _CARTOON_GLYPHS_AVAILABLE = False

# ---------------------------------------------------------------- constants
DEFAULT_INPUT = ("/Users/sven/Documents/Current_projects/SRK_bioinformatics/"
                 "Billinge and Robertson outcrossing data_Percent Fruit Set.xls")
DEFAULT_TABLES = "Tables"
DEFAULT_FIGURES = "figures"

TREATMENT_S = "S"          # selfing
TREATMENT_O = "O"          # outcrossing (raw XLS code "AS" is renamed on load)
RAW_TREATMENT_OUTCROSS = "AS"

BAND_SC = (0.0, 0.2)
BAND_PARTIAL = (0.2, 0.8)
BAND_SI = (0.8, 1.0)

COL_SC = "#1b7837"        # matches SI-status green
COL_PARTIAL = "#e08214"   # matches SI-status amber
COL_SI = "#b2182b"        # matches SI-status red

DEFAULT_B = 10_000
DEFAULT_MIN_N = 3         # per-EO minimum for BOTH n_S and n_O
DEFAULT_SEED = 20260624

EPSILON = 1e-6            # replaces 0 in percent fruit set to avoid /0 in bootstrap

SPECIES_LABEL = "Lepidium papilliferum"

# ---------------------------------------------------------------- I/O
def load_data(path: Path) -> pd.DataFrame:
    warnings.filterwarnings("ignore")
    df = pd.read_excel(path)
    df["EO Name"] = df["EO Name"].ffill()
    df = df[df["treatment"].isin([TREATMENT_S, RAW_TREATMENT_OUTCROSS])].copy()
    df = df.rename(columns={
        "EO Name": "EO_name",
        "Site Acronym": "Site",
        "slick spot": "SlickSpot",
        "treatment": "Treatment",
        "percent fruit set": "PercentFruitSet",
    })
    df["Treatment"] = df["Treatment"].replace({RAW_TREATMENT_OUTCROSS: TREATMENT_O})
    df["PercentFruitSet"] = df["PercentFruitSet"].astype(float)
    df.loc[df["PercentFruitSet"] == 0, "PercentFruitSet"] = EPSILON
    return df

def restrict_to_matched_sites(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Keep only sites (Site Acronym) present in BOTH S and O."""
    per_site = df.groupby("Site")["Treatment"].nunique()
    matched = sorted(per_site[per_site >= 2].index.tolist())
    unmatched = sorted(set(df["Site"].unique()) - set(matched))
    return df[df["Site"].isin(matched)].copy(), matched, unmatched

# ---------------------------------------------------------------- ISI + bootstrap
def isi_point(selfed: np.ndarray, outcross: np.ndarray) -> float:
    m_as = float(np.mean(outcross))
    if m_as == 0:
        return float("nan")
    return max(0.0, 1.0 - float(np.mean(selfed)) / m_as)

def bootstrap_isi(selfed: np.ndarray, outcross: np.ndarray, B: int,
                  rng: np.random.Generator) -> np.ndarray:
    s = np.asarray(selfed, dtype=float)
    a = np.asarray(outcross, dtype=float)
    s_boot = rng.choice(s, size=(B, s.size), replace=True).mean(axis=1)
    a_boot = rng.choice(a, size=(B, a.size), replace=True).mean(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        isi = 1.0 - s_boot / a_boot
    isi = np.where(a_boot == 0, np.nan, isi)
    return np.maximum(0.0, isi)

def classify(isi: float) -> str:
    if not np.isfinite(isi):
        return "Undefined"
    if isi < BAND_SC[1]:
        return "Self-compatible"
    if isi >= BAND_SI[0]:
        return "Self-incompatible"
    return "Partially self-incompatible"

def summarise(boot: np.ndarray, point: float) -> dict:
    finite = boot[np.isfinite(boot)]
    if finite.size == 0:
        return dict(ISI_point=point, ISI_median=np.nan, CI_lo=np.nan, CI_hi=np.nan,
                    Prop_SC=np.nan, Prop_Partial=np.nan, Prop_SI=np.nan,
                    Classification="Undefined", N_boot_finite=0)
    median = float(np.median(finite))
    ci_lo, ci_hi = np.percentile(finite, [2.5, 97.5])
    return dict(
        ISI_point=point,
        ISI_median=median,
        CI_lo=float(ci_lo),
        CI_hi=float(ci_hi),
        Prop_SC=float(np.mean(finite < BAND_SC[1])),
        Prop_Partial=float(np.mean((finite >= BAND_PARTIAL[0]) & (finite < BAND_PARTIAL[1]))),
        Prop_SI=float(np.mean(finite >= BAND_SI[0])),
        Classification=classify(median),
        N_boot_finite=int(finite.size),
    )

# ---------------------------------------------------------------- analysis
def analyse(df: pd.DataFrame, B: int, min_n: int, seed: int
            ) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    boots: dict[str, np.ndarray] = {}

    # Species-wide (pool all sites)
    s_all = df.loc[df["Treatment"] == TREATMENT_S, "PercentFruitSet"].to_numpy()
    a_all = df.loc[df["Treatment"] == TREATMENT_O, "PercentFruitSet"].to_numpy()
    b = bootstrap_isi(s_all, a_all, B, rng)
    rows.append({"Level": "Species", "Group": SPECIES_LABEL,
                 "N_S": int(s_all.size), "N_O": int(a_all.size),
                 **summarise(b, isi_point(s_all, a_all))})
    boots[f"Species::{SPECIES_LABEL}"] = b

    # Per-EO (only if both treatments meet min_n)
    for eo, sub in df.groupby("EO_name"):
        s = sub.loc[sub["Treatment"] == TREATMENT_S, "PercentFruitSet"].to_numpy()
        a = sub.loc[sub["Treatment"] == TREATMENT_O, "PercentFruitSet"].to_numpy()
        if s.size < min_n or a.size < min_n:
            continue
        b_eo = bootstrap_isi(s, a, B, rng)
        rows.append({"Level": "EO", "Group": eo, "N_S": int(s.size), "N_O": int(a.size),
                     **summarise(b_eo, isi_point(s, a))})
        boots[f"EO::{eo}"] = b_eo

    return pd.DataFrame(rows), boots

# ---------------------------------------------------------------- plotting
def shade_bands(ax) -> None:
    ax.axvspan(*BAND_SC, color=COL_SC, alpha=0.10, zorder=0)
    ax.axvspan(*BAND_PARTIAL, color=COL_PARTIAL, alpha=0.10, zorder=0)
    ax.axvspan(*BAND_SI, color=COL_SI, alpha=0.10, zorder=0)
    for x in (BAND_SC[1], BAND_SI[0]):
        ax.axvline(x, color="grey", linestyle=":", linewidth=0.8, zorder=1)
    trans = ax.get_xaxis_transform()          # x = data, y = axes fraction
    label_kw = dict(ha="center", va="top", fontsize=14,
                    fontweight="bold", transform=trans, zorder=1)
    ax.text(0.10, 0.97, "SC",         color=COL_SC,      **label_kw)
    ax.text(0.50, 0.97, "Partial-SI", color=COL_PARTIAL, **label_kw)
    ax.text(0.90, 0.97, "SI",         color=COL_SI,      **label_kw)

def plot_species(ax, boot: np.ndarray, row: dict, blank: bool = False,
                 line_only: bool = False) -> None:
    finite = boot[np.isfinite(boot)]
    kde = gaussian_kde(finite, bw_method=0.15)
    xs = np.linspace(-0.05, 1.05, 400)
    ys = kde(xs)                             # computed in every mode for consistent y-axis
    shade_bands(ax)
    med, lo, hi = row["ISI_median"], row["CI_lo"], row["CI_hi"]
    if line_only:
        # Minimalist: LEPA bootstrap distribution as a thick black line + median
        ax.plot(xs, ys, color="black", linewidth=3.0, zorder=5)
        ax.axvline(med, color="black", linewidth=2.0, zorder=5)
        # No ax title (talk figure — kept clean for slide use)
    elif not blank:
        ax.fill_between(xs, ys, color="#4a4a4a", alpha=0.55)
        ax.plot(xs, ys, color="black", linewidth=1.2)
        ax.axvline(med, color="black", linewidth=1.5)
        y_mark = ys.max() * 0.08
        ax.plot([lo, hi], [y_mark, y_mark], color="black", linewidth=3, solid_capstyle="butt")
    # No ax title on any variant (talk-ready)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, ys.max() * 1.45)
    ax.set_xlabel("Index of Self-incompatibility (ISI)")
    ax.set_ylabel("Bootstrap density")

def _add_band_annotation(ax, kind: str) -> None:
    """Overlay a small flower → silique glyph inside a classification band (blank
    figure only). kind='SC' → flower + silique-with-seeds in the green band;
    kind='SI' → flower + empty-silique-with-red-× in the red band. No-op if the
    cartoon glyphs are unavailable."""
    if not _CARTOON_GLYPHS_AVAILABLE:
        return
    if kind == "SC":
        inset = ax.inset_axes([0.02, 0.04, 0.185, 0.85])
        subtitle = "self  →  seeds"
        draw_result_glyph = draw_silique_with_seeds
    else:  # SI
        inset = ax.inset_axes([0.795, 0.04, 0.185, 0.85])
        subtitle = "self  →  no seeds"
        draw_result_glyph = draw_empty_silique
    inset.set_xlim(0, 3)
    inset.set_ylim(0, 6.5)
    inset.set_aspect("equal")
    inset.axis("off")
    inset.patch.set_alpha(0)
    # Flower + self-pollen loop at top
    flower_r = 0.55
    draw_flower(inset, 1.5, 4.75, radius=flower_r, stem_bottom=3.55)
    draw_self_loop(inset, 1.5, 4.75, flower_r)
    # Down arrow between flower and silique
    arr = FancyArrowPatch((1.5, 3.30), (1.5, 2.60),
                          arrowstyle="->", mutation_scale=17,
                          color="#666666", linewidth=1.8, zorder=6)
    inset.add_patch(arr)
    # Silique below
    draw_result_glyph(inset, 1.5, 1.60, size=0.50)
    # Subtitle at the bottom (top band label already carries "SC"/"SI")
    inset.text(1.5, 0.35, subtitle, ha="center", va="bottom",
               fontsize=10, fontstyle="italic", color="#555555")

def plot_figure(summary: pd.DataFrame, boots: dict[str, np.ndarray],
                out_pdf: Path, out_png: Path, blank: bool = False,
                annotate_sc: bool = False, annotate_si: bool = False,
                line_only: bool = False) -> None:
    species_row = summary[summary["Level"] == "Species"].iloc[0].to_dict()
    species_boot = boots[f"Species::{species_row['Group']}"]

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    plot_species(ax, species_boot, species_row, blank=blank, line_only=line_only)
    if annotate_sc:
        _add_band_annotation(ax, "SC")
    if annotate_si:
        _add_band_annotation(ax, "SI")

    suptitle = "Breeding system classification via ISI  —  Billinge & Robertson outcrossing data"
    if blank and not line_only:
        suptitle += "  [blank]"
    fig.suptitle(suptitle, fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

# ---------------------------------------------------------------- main
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, default=Path(DEFAULT_INPUT),
                   help="Path to the Billinge & Robertson xls file")
    p.add_argument("--tables-dir", type=Path, default=Path(DEFAULT_TABLES))
    p.add_argument("--figures-dir", type=Path, default=Path(DEFAULT_FIGURES))
    p.add_argument("-B", "--bootstrap-samples", type=int, default=DEFAULT_B)
    p.add_argument("--min-n", type=int, default=DEFAULT_MIN_N,
                   help="Min n_S AND n_O per EO for per-EO bootstrap")
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--include-unmatched", action="store_true",
                   help="Include S observations from sites without matched O "
                        "(default: restrict to sites with both treatments)")
    args = p.parse_args()

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    df = load_data(args.input)
    if not args.include_unmatched:
        df, matched, unmatched = restrict_to_matched_sites(df)
        print(f"[ISI] Restricted to {len(matched)} matched sites: {matched}")
        if unmatched:
            print(f"[ISI] Dropped S-only sites: {unmatched}")
    else:
        print("[ISI] Pooling across all sites (--include-unmatched)")
    summary, boots = analyse(df, args.bootstrap_samples, args.min_n, args.seed)

    for col in ("ISI_point", "ISI_median", "CI_lo", "CI_hi",
                "Prop_SC", "Prop_Partial", "Prop_SI"):
        summary[col] = summary[col].astype(float).round(4)

    out_tsv = args.tables_dir / "ISI_breeding_system_summary.tsv"
    out_pdf = args.figures_dir / "ISI_breeding_system.pdf"
    out_png = args.figures_dir / "ISI_breeding_system.png"
    out_pdf_blank = args.figures_dir / "ISI_breeding_system_blank.pdf"
    out_png_blank = args.figures_dir / "ISI_breeding_system_blank.png"
    out_pdf_blank_sc = args.figures_dir / "ISI_breeding_system_blank_SC.pdf"
    out_png_blank_sc = args.figures_dir / "ISI_breeding_system_blank_SC.png"
    out_pdf_blank_scsi = args.figures_dir / "ISI_breeding_system_blank_SC_SI.pdf"
    out_png_blank_scsi = args.figures_dir / "ISI_breeding_system_blank_SC_SI.png"
    out_pdf_annotated_data = args.figures_dir / "ISI_breeding_system_annotated_data.pdf"
    out_png_annotated_data = args.figures_dir / "ISI_breeding_system_annotated_data.png"

    summary.to_csv(out_tsv, sep="\t", index=False, encoding="utf-8")
    plot_figure(summary, boots, out_pdf, out_png, blank=False)
    plot_figure(summary, boots, out_pdf_blank, out_png_blank, blank=True)
    plot_figure(summary, boots, out_pdf_blank_sc, out_png_blank_sc,
                blank=True, annotate_sc=True)
    plot_figure(summary, boots, out_pdf_blank_scsi, out_png_blank_scsi,
                blank=True, annotate_sc=True, annotate_si=True)
    plot_figure(summary, boots, out_pdf_annotated_data, out_png_annotated_data,
                blank=True, annotate_sc=True, annotate_si=True, line_only=True)

    n_eo = int((summary["Level"] == "EO").sum())
    print(f"[ISI] Bootstrap B = {args.bootstrap_samples}, seed = {args.seed}")
    print(f"[ISI] Species-wide row + {n_eo} EO row(s) written to {out_tsv}")
    print(summary.to_string(index=False))
    print(f"[ISI] Figure: {out_pdf}")
    print(f"[ISI] Figure: {out_png}")
    print(f"[ISI] Figure (blank for talks): {out_pdf_blank}")
    print(f"[ISI] Figure (blank for talks): {out_png_blank}")
    print(f"[ISI] Figure (blank + SC glyph): {out_pdf_blank_sc}")
    print(f"[ISI] Figure (blank + SC glyph): {out_png_blank_sc}")
    print(f"[ISI] Figure (blank + SC + SI glyphs): {out_pdf_blank_scsi}")
    print(f"[ISI] Figure (blank + SC + SI glyphs): {out_png_blank_scsi}")
    print(f"[ISI] Figure (annotated + line data):  {out_pdf_annotated_data}")
    print(f"[ISI] Figure (annotated + line data):  {out_png_annotated_data}")

if __name__ == "__main__":
    main()
