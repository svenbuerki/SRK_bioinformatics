#!/usr/bin/env python3
"""Step 28c — Within-location pollen connectivity.

Reads the event-level spatial frame written by Step 28
(`step28_events_spatial_neighborhood.tsv`) and computes, for each LEPA
location, how much of its adult mating pool is actually connected via
pollen flow at three flight radii (10 m, 25 m, 50 m).

Method
------
For each location we build a graph on its events: two events are
connected by an edge if their haversine distance is <= R metres. We
then find connected components on that graph and report:

    n_events                 total events at the location (2025 field only)
    total_n_fertile          sum of N_fertile across those events
    n_components_R           number of connected components at radius R
    largest_component_share_R
                             fraction of the location's adults that sit in
                             the largest connected component
    connected_share_R        fraction of adults that sit in a component
                             containing more than one event (i.e. exchange
                             pollen with at least one other event)

Interpretation
--------------
`connected_share_R = 1.0` → the whole location behaves as one mating unit
at that pollen-flight radius. `connected_share_R = 0.0` → every event is
an isolated island; the "location" is really several disconnected
populations. Values in between reveal partial connectivity — some events
form mating clusters, others sit alone.

Outputs
-------
Tables/Phase5/step29_location_connectivity.tsv
    One row per location with the metrics above at R = 10 / 25 / 50 m.
figures/Phase5/step29_location_connectivity.png
    Per-location bar chart of `connected_share_10m`, panelled by BL
    (using the project-wide BL palette / order).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from step28_seed_sampling_per_mother import haversine_meters
from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl

DEFAULT_TABLES = Path("Tables/Phase5")
DEFAULT_FIGURES = Path("figures/Phase5")
INPUT_SPATIAL = DEFAULT_TABLES / "step28_events_spatial_neighborhood.tsv"
INPUT_LOCATIONS = DEFAULT_TABLES / "step29_sampling_per_location.tsv"
RADII_M = (10.0, 25.0, 50.0)


def _connected_components(adj: np.ndarray) -> list[list[int]]:
    """Simple BFS on a boolean adjacency matrix. Returns a list of lists
    of node indices, one list per connected component."""
    n = adj.shape[0]
    seen = np.zeros(n, bool)
    comps: list[list[int]] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [start]
        comp: list[int] = []
        while stack:
            v = stack.pop()
            if seen[v]:
                continue
            seen[v] = True
            comp.append(v)
            for u in np.where(adj[v] & ~seen)[0]:
                stack.append(int(u))
        comps.append(comp)
    return comps


def build_connectivity(spatial: pd.DataFrame,
                       radii_m: tuple[float, ...] = RADII_M) -> pd.DataFrame:
    rows = []
    for loc_id, sub in spatial.groupby("locationID"):
        n_ev = len(sub)
        n_fert = sub["n_fertile"].astype(int).values
        total_fert = int(n_fert.sum())
        lat = sub["lat"].values; lon = sub["lon"].values
        d = haversine_meters(
            lat[:, None], lon[:, None], lat[None, :], lon[None, :],
        )
        row = {
            "locationID":       int(loc_id),
            "n_events":         n_ev,
            "total_n_fertile":  total_fert,
        }
        for R in radii_m:
            R_i = int(round(R))
            adj = (d > 0) & (d <= R)
            comps = _connected_components(adj)
            n_components = len(comps)
            # size (in adults) of each component
            sizes = np.array([int(n_fert[c].sum()) for c in comps])
            largest = int(sizes.max()) if len(sizes) else 0
            # "connected" adults = in a component of >1 event
            multi_event = np.array([len(c) for c in comps]) > 1
            connected_adults = int(sizes[multi_event].sum())
            row[f"n_components_{R_i}m"]         = n_components
            row[f"largest_component_share_{R_i}m"] = (
                largest / total_fert if total_fert > 0 else 0.0
            )
            row[f"connected_share_{R_i}m"] = (
                connected_adults / total_fert if total_fert > 0 else 0.0
            )
        rows.append(row)
    return pd.DataFrame(rows)


def plot_connectivity_by_bl(df: pd.DataFrame, out_png: Path, out_pdf: Path,
                            radius_m: int = 25):
    """Per-location bars of `connected_share_{R}m`, panelled by BL."""
    df = df.copy()
    df["locationCode"] = df["locationCode"]
    df["BL"] = locationCode_to_bl(df["locationCode"]).values
    df["BL"] = df["BL"].fillna("Unassigned")

    bls = [b for b in BL_ORDER if b in df["BL"].values]
    if (df["BL"] == "Unassigned").any():
        bls.append("Unassigned")
    palette = {**BL_COLORS, "Unassigned": "#8a8a8a"}
    metric = f"connected_share_{radius_m}m"

    heights = [max(int((df["BL"] == b).sum()), 1) for b in bls]
    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(9.0, max(6.0, 0.26 * sum(heights) + 1.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub = df[df["BL"] == bl].sort_values(
            metric, ascending=True,
        ).reset_index(drop=True)
        y = np.arange(len(sub))
        colour = palette[bl]
        ax.barh(y, sub[metric], color=colour,
                edgecolor="white", height=0.8)
        # threshold markers: 50 % and 90 %
        for x, lab, col in [(0.5, " 50 %", "#e08214"),
                             (0.9, " 90 %", "#1b7837")]:
            ax.axvline(x, color=col, ls=":", lw=1.0, alpha=0.6)
        labels = [f"{code}   (events = {int(e)}, adults = {int(f)})"
                  for code, e, f in zip(sub["locationCode"],
                                        sub["n_events"],
                                        sub["total_n_fertile"])]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.7, len(sub) - 0.3)
        ax.text(1.01, 0.5, bl, transform=ax.transAxes,
                fontsize=13, fontweight="bold", color=colour,
                va="center", ha="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel(
        f"Fraction of adults connected to another event by pollen flow  "
        f"(within {radius_m} m)",
        fontsize=11,
    )
    fig.suptitle(
        f"Within-location pollen connectivity across LEPA locations  —  "
        f"pollen-flight radius {radius_m} m",
        fontsize=13, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 0.94, 0.97])
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_radius_sensitivity(df: pd.DataFrame, out_png: Path, out_pdf: Path):
    """Radius-sensitivity view: aggregate fraction of locations reaching
    50 % and 90 % within-location connectivity at each pollen-flight
    radius. Shows how much the fragmentation story shifts with the
    assumed foraging distance."""
    thresholds = [0.5, 0.9]
    cols = [c for c in df.columns
            if c.startswith("connected_share_") and c.endswith("m")]
    radii = sorted(
        int(c.replace("connected_share_", "").replace("m", ""))
        for c in cols
    )
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    x = np.arange(len(radii))
    width = 0.35
    for i, thr in enumerate(thresholds):
        frac = [(df[f"connected_share_{R}m"] >= thr).mean() for R in radii]
        offset = (i - 0.5) * width
        colour = "#e08214" if thr == 0.5 else "#1b7837"
        ax.bar(x + offset, frac, width=width, color=colour,
               edgecolor="white",
               label=f"≥ {int(thr*100)} % adults connected")
        for xi, f in zip(x + offset, frac):
            ax.text(xi, f + 0.015, f"{f:.0%}",
                    ha="center", va="bottom", fontsize=9, color=colour)

    ax.set_xticks(x); ax.set_xticklabels([f"{R} m" for R in radii],
                                          fontsize=11)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Assumed pollen-flight radius", fontsize=11)
    ax.set_ylabel("Fraction of LEPA locations", fontsize=11)
    ax.set_title(
        "Within-location pollen connectivity — radius sensitivity",
        fontsize=12,
    )
    ax.legend(loc="upper left", fontsize=10, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200); fig.savefig(out_pdf); plt.close(fig)


def main() -> None:
    if not INPUT_SPATIAL.exists():
        raise SystemExit(
            f"Missing {INPUT_SPATIAL} — run step28 first.")
    spatial = pd.read_csv(INPUT_SPATIAL, sep="\t", encoding="utf-8-sig")

    # Optional location codes (join for display purposes)
    if INPUT_LOCATIONS.exists():
        loc_meta = pd.read_csv(INPUT_LOCATIONS, sep="\t",
                                encoding="utf-8-sig")
        loc_code_map = (loc_meta.drop_duplicates("locationID")
                        .set_index("locationID")["locationCode"].to_dict())
    else:
        loc_code_map = {}

    df = build_connectivity(spatial)
    df["locationCode"] = df["locationID"].map(loc_code_map).fillna(
        df["locationID"].astype(str))
    # tidy column order
    df = df[
        ["locationID", "locationCode", "n_events", "total_n_fertile"]
        + [c for c in df.columns
           if c.startswith(("n_components_", "largest_component_share_",
                            "connected_share_"))]
    ]
    out_tsv = DEFAULT_TABLES / "step29_location_connectivity.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[step29] Wrote {out_tsv}")

    # High-level report
    for R in RADII_M:
        R_i = int(round(R))
        n_fully = int((df[f"connected_share_{R_i}m"] >= 0.9).sum())
        n_frag  = int((df[f"connected_share_{R_i}m"] <  0.5).sum())
        print(f"[step29] R = {R_i} m: "
              f"{n_fully}/{len(df)} locations ≥ 90 % connected · "
              f"{n_frag}/{len(df)} locations < 50 % connected.")

    DEFAULT_FIGURES.mkdir(parents=True, exist_ok=True)
    # Primary figure at the small-bee 10 m assumption, and sensitivity
    # figures at 25 m and 50 m for reviewers who want to see how the
    # picture changes with a longer pollen-flight assumption.
    for R in RADII_M:
        R_i = int(round(R))
        # 25 m is the primary radius (unlabelled), 10 m and 50 m are sensitivity.
        suffix = "" if R_i == 25 else f"_{R_i}m"
        plot_connectivity_by_bl(
            df,
            out_png=DEFAULT_FIGURES / f"step29_location_connectivity{suffix}.png",
            out_pdf=DEFAULT_FIGURES / f"step29_location_connectivity{suffix}.pdf",
            radius_m=R_i,
        )
    plot_radius_sensitivity(
        df,
        out_png=DEFAULT_FIGURES / "step29_location_connectivity_radius_sensitivity.png",
        out_pdf=DEFAULT_FIGURES / "step29_location_connectivity_radius_sensitivity.pdf",
    )
    print(f"[step29] Figures in {DEFAULT_FIGURES}/")


if __name__ == "__main__":
    main()
