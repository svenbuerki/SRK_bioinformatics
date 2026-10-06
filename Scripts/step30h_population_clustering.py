"""Step 30h — hierarchical geographic clustering of the 44 Phase 5
populations.

Phase 5 § A.4.5 (new). Replicates the content of the external
LEPA_EO_spatial_clustering repo (Ward's D2 at the EO level) but now
at the POPULATION level, using the 2 years of pooled 2025+2026 data
that produced the 44 populations in step29a. The purpose is to
re-check whether the Bottleneck-Lineage (BL) grouping used
throughout Phase 5 still holds, and to let the data pick the BL
cardinality k rather than forcing k = 5.

Approach
--------
1. Compute per-population centroid (lat, lon) from the
   step29a_population_crosswalk.tsv — a size-weighted mean of event
   coordinates across both years so that the centroid tracks the
   bulk of the above-ground sample.
2. Build the pairwise haversine distance matrix (metres).
3. Ward's D2 linkage on that distance matrix (standard choice for
   geographic hierarchical clustering of this scale, matches the
   external repo).
4. Silhouette for k = 2 … 10; k_optimal = argmax silhouette.
5. Dendrogram with the k_optimal cut highlighted; silhouette curve
   alongside.

Outputs
-------
Tables/Phase5/step30h_population_distances.tsv
    long-form distance matrix (populationA, populationB,
    distance_m, pair label).
Tables/Phase5/step30h_cluster_assignments.tsv
    populationID → cluster label at k_optimal, k=5 (for comparison
    with the old BL cardinality), and a tagged `cluster_k_optimal`
    column that becomes the input to Phase II / III.
Tables/Phase5/step30h_silhouette_curve.tsv
    k vs silhouette coefficient; the whole curve for the report.
figures/Phase5/step30h_dendrogram.png / .pdf
figures/Phase5/step30h_silhouette_curve.png / .pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

from step28_seed_sampling_per_mother import haversine_meters

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

K_RANGE = range(2, 11)


# ---------------------------------------------------------------------------
# 1. Centroids
# ---------------------------------------------------------------------------
def load_population_centroids() -> pd.DataFrame:
    """Per-population centroid = size-weighted mean of event (lat, lon)
    across both years, so the centroid sits where the plants actually
    are rather than at the midpoint of a sparsely-sampled hull."""
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                      sep="\t", encoding="utf-8-sig")
    # weight by n_fertile (fall back to 1 if NaN or 0)
    cw["w"] = cw["n_fertile"].fillna(0).clip(lower=0)
    cw.loc[cw["w"] == 0, "w"] = 1
    rows = []
    for pop_id, g in cw.groupby("populationID"):
        rows.append({
            "populationID": int(pop_id),
            "lat":          float(np.average(g["lat"], weights=g["w"])),
            "lon":          float(np.average(g["lon"], weights=g["w"])),
            "n_fert_total": int(g["n_fertile"].sum()),
            "n_events":     int(g["eventID"].nunique()),
        })
    return pd.DataFrame(rows).sort_values("populationID").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 2. Distances + linkage
# ---------------------------------------------------------------------------
def pairwise_distance_matrix(cents: pd.DataFrame) -> np.ndarray:
    lat = cents["lat"].to_numpy()
    lon = cents["lon"].to_numpy()
    return haversine_meters(
        lat[:, None], lon[:, None],
        lat[None, :], lon[None, :],
    )


def long_form_distances(cents: pd.DataFrame,
                          D: np.ndarray) -> pd.DataFrame:
    rows = []
    ids = cents["populationID"].to_numpy()
    n = len(ids)
    for i in range(n):
        for j in range(i + 1, n):
            rows.append({
                "populationA": int(ids[i]),
                "populationB": int(ids[j]),
                "distance_m":  float(D[i, j]),
            })
    df = pd.DataFrame(rows).sort_values("distance_m").reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 3. Silhouette for k selection
# ---------------------------------------------------------------------------
def silhouette_curve(D: np.ndarray, Z: np.ndarray,
                       k_range) -> pd.DataFrame:
    rows = []
    for k in k_range:
        labels = fcluster(Z, t=k, criterion="maxclust")
        if len(set(labels)) < 2:
            s = np.nan
        else:
            s = silhouette_score(D, labels, metric="precomputed")
        rows.append({"k": int(k), "silhouette": float(s),
                      "n_clusters_actual": int(len(set(labels)))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
def plot_dendrogram(Z: np.ndarray, cents: pd.DataFrame,
                      k_optimal: int,
                      out_png: Path, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 7))
    labels = [f"P{int(r):>2}" for r in cents["populationID"]]
    dendrogram(
        Z,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=Z[-k_optimal + 1, 2] if k_optimal > 1 else 0,
        ax=ax,
    )
    # Cut line
    if k_optimal > 1:
        cut = Z[-k_optimal + 1, 2]
        ax.axhline(cut, color="#444", linestyle="--", linewidth=1.2,
                   label=f"k = {k_optimal} cut (silhouette optimum)")
        ax.legend(loc="upper right", fontsize=10)
    ax.set_ylabel("Haversine linkage distance (m, Ward's D2)")
    ax.set_title(
        f"Population hierarchical clustering — Ward's D2 on pairwise "
        f"haversine distances\n(44 populations, pooled 2025 + 2026 data; "
        f"silhouette-optimal k = {k_optimal})",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h] Wrote {out_png.name} + .pdf")


def plot_silhouette_curve(sil: pd.DataFrame, k_optimal: int,
                            out_png: Path, out_pdf: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(sil["k"], sil["silhouette"], "-o",
            color="#5e3c99", linewidth=2, markersize=7)
    ax.axvline(k_optimal, color="#D55E00", linestyle="--", linewidth=1.5,
               label=f"k_optimal = {k_optimal}")
    ax.axvline(5, color="#aaaaaa", linestyle=":", linewidth=1.3,
               label="k = 5 (old BL framework)")
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Mean silhouette coefficient")
    ax.set_title(
        "Silhouette-based k selection for Ward's D2 population "
        f"clustering (n = 44)\nOptimal k = {k_optimal} "
        f"(silhouette = {sil.loc[sil['k'] == k_optimal, 'silhouette'].iat[0]:.3f})",
        fontsize=11,
    )
    ax.set_xticks(list(sil["k"]))
    ax.legend(loc="best", fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30h] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# 5. Driver
# ---------------------------------------------------------------------------
def main() -> None:
    cents = load_population_centroids()
    print(f"[step30h] Loaded {len(cents)} population centroids")

    D = pairwise_distance_matrix(cents)
    long_D = long_form_distances(cents, D)
    long_D.to_csv(TABLES / "step30h_population_distances.tsv",
                   sep="\t", index=False)
    print(f"[step30h] Wrote step30h_population_distances.tsv "
          f"(nearest pair {long_D.iloc[0]['distance_m']:.0f} m; "
          f"farthest pair {long_D.iloc[-1]['distance_m'] / 1000:.1f} km)")

    # Linkage expects a condensed distance matrix (upper triangle flat)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="ward")

    sil = silhouette_curve(D, Z, K_RANGE)
    sil.to_csv(TABLES / "step30h_silhouette_curve.tsv",
                sep="\t", index=False)
    k_optimal = int(sil.loc[sil["silhouette"].idxmax(), "k"])
    print(f"[step30h] Silhouette curve:")
    print(sil.to_string(index=False))
    print(f"[step30h] k_optimal = {k_optimal} "
          f"(silhouette = {sil.loc[sil['k'] == k_optimal, 'silhouette'].iat[0]:.3f})")

    # Cluster assignments: always expose both an explicit
    # `cluster_k_optimal` column (downstream scripts key on it) AND
    # the `cluster_k5` column for continuity with the old EO BL
    # framework. When k_optimal happens to be 5 the two are
    # identical, which is itself the message.
    labels_opt = fcluster(Z, t=k_optimal, criterion="maxclust")
    labels_5   = fcluster(Z, t=5,         criterion="maxclust")
    assigns = cents.copy()
    assigns["cluster_k_optimal"] = labels_opt
    assigns["cluster_k5"]        = labels_5
    assigns["k_optimal"]         = k_optimal
    assigns.to_csv(TABLES / "step30h_cluster_assignments.tsv",
                    sep="\t", index=False)
    print(f"[step30h] Wrote step30h_cluster_assignments.tsv")

    # Figures
    plot_dendrogram(
        Z, cents, k_optimal,
        FIGURES / "step30h_dendrogram.png",
        FIGURES / "step30h_dendrogram.pdf",
    )
    plot_silhouette_curve(
        sil, k_optimal,
        FIGURES / "step30h_silhouette_curve.png",
        FIGURES / "step30h_silhouette_curve.pdf",
    )

    # Headline: cluster sizes at both k levels
    print()
    print(f"[step30h] Cluster sizes at k = {k_optimal} "
          f"(silhouette-optimal):")
    print(pd.Series(labels_opt).value_counts().sort_index().to_string())
    print()
    print(f"[step30h] Cluster sizes at k = 5 (old BL cardinality):")
    print(pd.Series(labels_5).value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()
