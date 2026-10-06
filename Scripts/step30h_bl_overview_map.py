"""Step 30h / Phase V — BL overview map on the Snake River Plain.

Phase 5 § A.4.5. Simple geographic map: one dot per population,
coloured by new BL, sized by `n_fert_total`, with the new population
IDs labeled. Thin BL convex-hull outlines reinforce the cluster
structure.

Outputs
-------
figures/Phase5/step30h_overview_by_BL.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

# Okabe-Ito palette, extended for 5 BLs
BL_COLORS = {
    "BL1": "#0072B2",   # blue
    "BL2": "#009E73",   # teal green
    "BL3": "#E69F00",   # orange
    "BL4": "#D55E00",   # vermillion
    "BL5": "#CC79A7",   # pink
}


def main() -> None:
    pts = pd.read_csv(TABLES / "step30h_cluster_assignments.tsv",
                       sep="\t", encoding="utf-8-sig")
    bl_def = pd.read_csv(TABLES / "step30h_bl_definition.tsv",
                           sep="\t", encoding="utf-8-sig")
    # Phase III remapped cluster_assignments: it already carries the
    # new populationID + BL columns.

    fig, ax = plt.subplots(figsize=(11.5, 7.5))
    size_min, size_max = 30, 500
    n_max = max(pts["n_fert_total"].max(), 1)
    sizes = size_min + (size_max - size_min) * np.sqrt(
        pts["n_fert_total"] / n_max
    )

    # BL hulls
    for bl, grp in pts.groupby("BL"):
        color = BL_COLORS.get(bl, "#777777")
        if len(grp) >= 3:
            hull_pts = grp[["lon", "lat"]].to_numpy()
            try:
                hull = ConvexHull(hull_pts)
                poly = hull_pts[hull.vertices]
                poly = np.vstack([poly, poly[:1]])
                ax.fill(poly[:, 0], poly[:, 1], color=color,
                        alpha=0.08, zorder=1)
                ax.plot(poly[:, 0], poly[:, 1], color=color,
                        linewidth=1.3, linestyle="--", alpha=0.75,
                        zorder=2)
            except Exception:
                pass

    # Population dots coloured by BL
    for bl, grp in pts.groupby("BL"):
        color = BL_COLORS.get(bl, "#777777")
        s = sizes.loc[grp.index].to_numpy()
        row = bl_def[bl_def["BL"] == bl].iloc[0]
        ax.scatter(grp["lon"], grp["lat"], s=s, c=color,
                   alpha=0.85, edgecolor="white", linewidth=0.7,
                   zorder=3,
                   label=f"{bl}  (n = {int(row['n_populations'])}, "
                         f"{row['convex_hull_area_km2']:.0f} km²)")

    # Label populations (only when there's room — skip the densest BL5)
    for _, r in pts.iterrows():
        ax.annotate(f"P{int(r['populationID'])}",
                    (r["lon"], r["lat"]),
                    xytext=(5, 3), textcoords="offset points",
                    fontsize=7, color="#333", alpha=0.85)

    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(
        "Phase 5 populations across the Snake River Plain — "
        "new BL framework (Ward's D2, k = 5, silhouette = 0.732)\n"
        "Dot size ∝ total N_fertile across 2025 + 2026; "
        "dashed hulls = BL convex outlines.",
        fontsize=11,
    )
    ax.legend(loc="upper right", fontsize=9, frameon=True, title="BL")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_png = FIGURES / "step30h_overview_by_BL.png"
    out_pdf = FIGURES / "step30h_overview_by_BL.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h-V] Wrote {out_png.name} + .pdf")


if __name__ == "__main__":
    main()
