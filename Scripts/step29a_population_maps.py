"""Step 29a-maps — population and across-year occupancy maps.

Two figure sets:

1. Overview map (all 44 populations, Snake River Plain).
   One dot per population, size ∝ total n_fertile across 2025 + 2026,
   colour-coded by occupancy pattern (both_years / 2025_only /
   2026_only), panelled by Bottleneck Lineage.

2. Per-population zoom panels for the selected candidates.
   For each candidate population: 2025 events (blue), 2026 events
   (vermillion), both-years slickspots (purple ring), 50 m deme
   partition (per-year polygons via convex hull), 500 m population
   outline. One figure per candidate; a 2x2 panel compares
   candidates.

Reads Tables/Phase5/step29a_* and (if present) the step30g
candidate selection; writes figures/Phase5/step29a_*.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from scipy.spatial import ConvexHull

from srk_bl_constants import BL_COLORS, BL_ORDER, locationCode_to_bl, make_location_label

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

OCC_COLORS = {
    "both_years": "#5e3c99",   # purple
    "2025_only":  "#0072B2",   # Okabe-Ito blue
    "2026_only":  "#D55E00",   # Okabe-Ito vermillion
}


# ---------------------------------------------------------------------------
# 1. Overview map
# ---------------------------------------------------------------------------
def plot_overview(pop_summary: pd.DataFrame,
                   slickspots: pd.DataFrame,
                   out_png: Path, out_pdf: Path) -> None:
    pop_centroids = (slickspots.groupby("populationID", as_index=False)
                                .agg(lat=("lat", "mean"),
                                     lon=("lon", "mean")))
    df = pop_summary.merge(pop_centroids, on="populationID")

    # BL from the first locationCode of each population (modal)
    def _modal_bl(codes: str) -> str:
        codes_list = [c.strip() for c in codes.split(",") if c.strip()]
        if not codes_list:
            return "Unassigned"
        bls = locationCode_to_bl(pd.Series(codes_list)).dropna()
        return bls.mode().iat[0] if len(bls) else "Unassigned"

    df["BL"] = df["locationCodes"].astype(str).apply(_modal_bl)

    fig, ax = plt.subplots(figsize=(10, 6.5))
    size_min, size_max = 25, 450
    n_max = max(df["n_fertile_total"].max(), 1)
    for occ, color in OCC_COLORS.items():
        sub = df[df["occupancy"] == occ]
        if sub.empty:
            continue
        sizes = size_min + (size_max - size_min) * np.sqrt(
            sub["n_fertile_total"] / n_max
        )
        ax.scatter(sub["lon"], sub["lat"], s=sizes, c=color,
                   alpha=0.75, edgecolor="white", linewidth=0.8,
                   label=f"{occ} ({len(sub)})")

    # Light grey dots for BL framing
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(
        "Phase 5 populations across the Snake River Plain — "
        "2025 + 2026 above-ground occupancy\n"
        "One dot per population; size ∝ total n_fertile across both "
        "years; colour = occupancy pattern.",
        fontsize=11,
    )
    ax.legend(loc="upper right", fontsize=9, frameon=True,
              title="Occupancy")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step29a-maps] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# 2. Per-population zoom panels
# ---------------------------------------------------------------------------
def _plot_one_population(ax: plt.Axes,
                          events_pop: pd.DataFrame,
                          slickspots_pop: pd.DataFrame,
                          title: str) -> None:
    # Group events by year for colour
    for yr, color in [(2025, OCC_COLORS["2025_only"]),
                       (2026, OCC_COLORS["2026_only"])]:
        sub = events_pop[events_pop["event_year"] == yr]
        if sub.empty:
            continue
        sizes = 25 + 60 * np.sqrt(np.clip(sub["n_fertile"], 1, None))
        ax.scatter(sub["lon"], sub["lat"], s=sizes, c=color, alpha=0.65,
                   edgecolor="white", linewidth=0.6,
                   label=f"{yr} events (n = {len(sub)}, "
                         f"N_fertile = {int(sub['n_fertile'].sum())})")

    # Both-years slickspots: ring
    both = slickspots_pop[slickspots_pop["both_years"]]
    if len(both):
        ax.scatter(both["lon"], both["lat"], s=180, facecolor="none",
                   edgecolor=OCC_COLORS["both_years"], linewidth=2.0,
                   label=f"both-year slickspot (n = {len(both)})")

    # Convex hull of all event coords = rough population footprint
    pts = events_pop[["lat", "lon"]].to_numpy()
    if len(pts) >= 3:
        try:
            hull = ConvexHull(pts[:, [1, 0]])  # lon,lat
            poly = pts[hull.vertices][:, [1, 0]]
            poly = np.vstack([poly, poly[:1]])
            ax.plot(poly[:, 0], poly[:, 1], ":",
                    color="#666666", linewidth=1.2, alpha=0.7,
                    label="event hull")
        except Exception:
            pass

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=7, frameon=True)


def plot_candidate_panels(candidates: pd.DataFrame,
                            events_crosswalk: pd.DataFrame,
                            slickspots: pd.DataFrame,
                            out_png: Path, out_pdf: Path,
                            suptitle: str | None = None) -> None:
    pops = candidates["populationID"].tolist()
    n = len(pops)
    if n == 0:
        print("[step29a-maps] No candidates to plot.")
        return

    cols = min(n, 2)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(6.5 * cols, 5.5 * rows),
                              squeeze=False)
    for i, pop_id in enumerate(pops):
        r, c = divmod(i, cols)
        ax = axes[r, c]
        events_pop = events_crosswalk[
            events_crosswalk["populationID"] == pop_id
        ]
        slick_pop = slickspots[slickspots["populationID"] == pop_id]
        info = candidates[candidates["populationID"] == pop_id].iloc[0]
        role = info.get("role", "")
        # Project-standard labels: {locationCode}_{locationID} per
        # merged location, comma-separated.
        pairs = (events_pop[["locationCode", "locationID"]]
                    .drop_duplicates().sort_values("locationID"))
        labels = []
        for _, p in pairs.iterrows():
            code = str(p["locationCode"]).strip()
            lid = int(p["locationID"])
            if code and code != "nan":
                labels.append(make_location_label(code, lid))
            else:
                labels.append(f"locID_{lid}")
        location_label = ", ".join(labels)
        title = (f"populationID = {pop_id}  {role}\n"
                 f"locations: {location_label}  |  "
                 f"N_fert 2025/2026 = {info['n_fertile_2025']}/{info['n_fertile_2026']}  "
                 f"|  slickspots {info['n_slickspots']} "
                 f"(both-years {info['n_slickspots_both_years']})")
        _plot_one_population(ax, events_pop, slick_pop, title)
    # hide any extra axes
    for i in range(n, rows * cols):
        r, c = divmod(i, cols)
        axes[r, c].axis("off")
    fig.suptitle(
        suptitle or
        "Per-population zoom panels — 2025 vs 2026 event occupancy",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step29a-maps] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# 3. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    pop_summary = pd.read_csv(TABLES / "step29a_population_summary.tsv",
                               sep="\t", encoding="utf-8-sig")
    slickspots  = pd.read_csv(TABLES / "step29a_slickspot_summary.tsv",
                               sep="\t", encoding="utf-8-sig")
    crosswalk   = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                               sep="\t", encoding="utf-8-sig")

    # ---- Overview ----
    plot_overview(
        pop_summary, slickspots,
        FIGURES / "step29a_populations_overview.png",
        FIGURES / "step29a_populations_overview.pdf",
    )

    # ---- Candidate pair (stable LARGE + SMALL) ----
    pair_path = TABLES / "step30g_stable_candidate_pair.tsv"
    if pair_path.exists():
        pair = pd.read_csv(pair_path, sep="\t", encoding="utf-8-sig")
        cand_ids = pair["populationID"].astype(int).tolist()
        roles = {
            int(r["populationID"]): f"({r['role']} — stable in both years)"
            for _, r in pair.iterrows()
        }
        cand_df = pop_summary[pop_summary["populationID"].isin(cand_ids)].copy()
        cand_df["role"] = cand_df["populationID"].map(roles)
        cand_df["sort_order"] = cand_df["populationID"].map(
            {pid: i for i, pid in enumerate(cand_ids)}
        )
        cand_df = cand_df.sort_values("sort_order").reset_index(drop=True)
        plot_candidate_panels(
            cand_df, crosswalk, slickspots,
            FIGURES / "step29a_candidate_populations_zoom.png",
            FIGURES / "step29a_candidate_populations_zoom.pdf",
            suptitle=(
                "Stable candidate pair — 2025 vs 2026 event occupancy "
                "at the LARGE and SMALL stable populations"
            ),
        )

    # ---- Crash populations (clearest within-pipeline H2 candidates) ----
    crash_path = TABLES / "step30g_crash_candidates.tsv"
    if crash_path.exists():
        crash = pd.read_csv(crash_path, sep="\t", encoding="utf-8-sig")
        if len(crash):
            crash_ids = crash["populationID"].astype(int).tolist()
            crash_df = pop_summary[pop_summary["populationID"].isin(crash_ids)].copy()
            crash_roles = {}
            for _, r in crash.iterrows():
                pid = int(r["populationID"])
                crash_roles[pid] = (
                    f"(crash — 2026/2025 = {r['crash_ratio']:.2f}, "
                    f"Δdiv = {r['delta_diversity']:+.1f})"
                )
            crash_df["role"] = crash_df["populationID"].map(crash_roles)
            crash_df["sort_order"] = crash_df["populationID"].map(
                {pid: i for i, pid in enumerate(crash_ids)}
            )
            crash_df = crash_df.sort_values("sort_order").reset_index(drop=True)
            plot_candidate_panels(
                crash_df, crosswalk, slickspots,
                FIGURES / "step29a_crash_populations_zoom.png",
                FIGURES / "step29a_crash_populations_zoom.pdf",
                suptitle=(
                    "Crash populations — N_fert 2026 << 2025 with "
                    "negative Δdiversity (strongest within-pipeline H2 candidates)"
                ),
            )


if __name__ == "__main__":
    main()
