"""Step 30f — SRK diversity hypothesis test via radius sweep.

Phase 5 § C.0.c. Companion to § C.0.b (step30e), which decomposes the
pollen-compatibility prediction. This script addresses the parallel
question for SRK **allele diversity**: at EO67, EO70, EO76 (the three
clean-overlap locations), the 50 m prediction over-estimates observed
SRK allele counts by ~20 alleles (predicted 11 / 28 / 32 vs observed
7 / 6 / 9). Two competing hypotheses:

  H1 — operational deme too wide at 50 m. Realised gene flow is
       tighter; a smaller radius (10-25 m) partitions each location
       into more and smaller demes, each drifting independently, so
       the per-location union of per-deme Fg sets is smaller.
       Diagnostic signature: predicted-vs-observed curves cross a
       smaller radius than 50 m.

  H2 — per-deme drift history beyond the species-wide prior P1.
       Decades of local drift have pushed each deme's SRK pool below
       what P1 predicts even at the correct deme size. Diagnostic
       signature: predicted-vs-observed gap persists at every radius.

Approach
--------
For each (EO, radius) combination:
  1. Build the connectivity graph on the 2025 LEPA events belonging
     to that location, with an edge wherever two events sit within
     `radius_m` metres (haversine). Connected components = demes.
  2. For each deme c: draw `PLOIDY * component_N_fertile_c` alleles
     from P1 (point-estimate frequency vector), union the resulting
     Fg sets across demes to get the per-replicate location pool
     size. Repeat `N_REPLICATES` times for a 95 % CI.
  3. Compare predicted mean + CI against the observed count
     (`obs_distinct_Fgs` in step30_B_partC_clean_overlap_per_location.tsv).

Reuses `step28.load_all_events + haversine_meters`, `step28.PLOIDY`,
`step30.build_p1_prior`.

Outputs
-------
Tables/Phase5/step30f_srk_diversity_radius_sweep.tsv
Figures/Phase5/step30f_srk_diversity_radius_sweep.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from step28_seed_sampling_per_mother import (
    DEFAULT_DB, load_all_events, haversine_meters, PLOIDY,
)
from step30_srk_diversity_prediction_vs_observed import build_p1_prior

TABLES = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

RADII_M = [10, 25, 50, 75, 100, 150]
N_REPLICATES = 2000
RNG_SEED = 2030

# Clean-overlap locations with meaningful fragmentation across the
# radius sweep (deme count varies 10 -> 150 m). EO67 has only 2
# events > 150 m apart, so its deme partition is INVARIANT at every
# tested radius -> no information for the H1-vs-H2 test; it is
# included in the TSV for the record but omitted from the figure.
# locationIDs taken from Phase 2 Locations table and confirmed in
# step30_B_partC_clean_overlap_per_location.tsv.
TARGETS = [
    ("EO67", 39, "#8a8a8a", False),  # grey; invariant deme count, TSV only
    ("EO70", 26, "#D55E00", True),   # vermillion; plotted
    ("EO76",  2, "#009E73", True),   # teal green; plotted
]


def components_at_radius(sub: pd.DataFrame, radius_m: float) -> list[list[int]]:
    """Return a list of positional-index lists, one per connected
    component in the <= radius_m adjacency graph. BFS."""
    n = len(sub)
    if n == 0:
        return []
    lat = sub["lat"].values
    lon = sub["lon"].values
    d = haversine_meters(
        lat[:, None], lon[:, None],
        lat[None, :], lon[None, :],
    )
    adj = (d > 0) & (d <= radius_m)
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


def simulate_location_pool(comp_sizes: list[int],
                            p1_vec: np.ndarray,
                            rng: np.random.Generator,
                            n_reps: int) -> np.ndarray:
    """Per replicate: for each deme, draw PLOIDY * size alleles from P1
    and record its set of present Fgs; location pool = union across
    demes. Returns an array of n_reps distinct-Fg counts."""
    K = len(p1_vec)
    counts = np.zeros(n_reps, dtype=int)
    for i in range(n_reps):
        present = np.zeros(K, bool)
        for size in comp_sizes:
            n_alleles = PLOIDY * int(size)
            if n_alleles > 0:
                draws = rng.choice(K, size=n_alleles, replace=True, p=p1_vec)
                present[draws] = True
        counts[i] = int(present.sum())
    return counts


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)

    prior = build_p1_prior(TABLES / "step26i_L1_carrier_inventory.tsv")
    p1_vec = prior["f_mean"].values
    K_fg = len(p1_vec)
    print(f"[step30f] P1 vector: K_fg = {K_fg}, "
          f"FG001 frequency = {p1_vec[0]:.3f}")

    events = load_all_events(DEFAULT_DB, year=2025)
    print(f"[step30f] Loaded {len(events)} LEPA events from 2025")

    obs = pd.read_csv(TABLES / "step30_B_partC_clean_overlap_per_location.tsv",
                      sep="\t", encoding="utf-8-sig")
    obs_by_id = obs.set_index("locationID")["obs_distinct_Fgs"].to_dict()

    rows: list[dict] = []
    for eo_code, loc_id, _color, _plot in TARGETS:
        sub = events[events["locationID"] == loc_id].copy()
        print(f"\n[step30f] {eo_code} (locationID={loc_id}) — "
              f"{len(sub)} events, N_fertile total = {int(sub['n_fertile'].sum())}")
        for r in RADII_M:
            comps = components_at_radius(sub, r)
            comp_sizes = [int(sub.iloc[idx]["n_fertile"].sum()) for idx in comps]
            counts = simulate_location_pool(comp_sizes, p1_vec, rng, N_REPLICATES)
            obs_count = int(obs_by_id[loc_id])
            rows.append({
                "EO":                  eo_code,
                "locationID":          loc_id,
                "radius_m":            r,
                "n_components":        len(comps),
                "n_fert_total":        int(sum(comp_sizes)),
                "component_sizes":     ",".join(str(s) for s in
                                                 sorted(comp_sizes, reverse=True)),
                "largest_comp_N":      int(max(comp_sizes)) if comp_sizes else 0,
                "pred_mean":           float(counts.mean()),
                "pred_lo95":           float(np.quantile(counts, 0.025)),
                "pred_hi95":           float(np.quantile(counts, 0.975)),
                "observed":            obs_count,
                "gap_pred_minus_obs":  float(counts.mean() - obs_count),
                "obs_inside_CI":       bool(
                    np.quantile(counts, 0.025) <= obs_count
                    <= np.quantile(counts, 0.975)
                ),
            })
            print(f"  r = {r:>3} m  |  "
                  f"{len(comps):>2} demes  "
                  f"(largest {max(comp_sizes) if comp_sizes else 0:>4}, "
                  f"smallest {min(comp_sizes) if comp_sizes else 0:>3})  |  "
                  f"predicted {counts.mean():5.1f} "
                  f"[{np.quantile(counts, 0.025):4.1f}, "
                  f"{np.quantile(counts, 0.975):4.1f}]  "
                  f"vs observed {obs_count}  "
                  f"{'(inside CI)' if np.quantile(counts, 0.025) <= obs_count <= np.quantile(counts, 0.975) else ''}")

    df = pd.DataFrame(rows)
    out_tsv = TABLES / "step30f_srk_diversity_radius_sweep.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[step30f] Wrote {out_tsv}")

    # ---- Figure — gap-focused single panel ----
    # y-axis = predicted - observed (distinct SRK alleles). Zero line =
    # perfect match. Observed values appear only as right-side labels
    # (not horizontal lines, because they are single-measurement adult
    # counts, not an outcome of the radius sweep).
    fig, ax = plt.subplots(figsize=(9.0, 5.2))

    plotted_any = False
    for eo_code, loc_id, color, plot_it in TARGETS:
        if not plot_it:
            continue
        sub = df[df["EO"] == eo_code].sort_values("radius_m")
        obs_val = int(sub["observed"].iloc[0])
        gap_mean = sub["pred_mean"].values - obs_val
        gap_lo = sub["pred_lo95"].values - obs_val
        gap_hi = sub["pred_hi95"].values - obs_val

        ax.plot(sub["radius_m"], gap_mean, "-o",
                color=color, linewidth=2.4, markersize=7,
                label=f"{eo_code}  (observed = {obs_val} alleles)")
        ax.fill_between(sub["radius_m"], gap_lo, gap_hi,
                         color=color, alpha=0.18)

        # Right-side label: "+22" style gap at the largest radius
        x_right = sub["radius_m"].iloc[-1]
        y_right = gap_mean[-1]
        ax.annotate(f"gap: +{y_right:.0f}", xy=(x_right, y_right),
                    xytext=(10, 0), textcoords="offset points",
                    color=color, fontsize=10, fontweight="bold",
                    va="center", ha="left")

        plotted_any = True

    ax.axhline(0, color="#444", linewidth=1.4, zorder=0)
    ax.text(RADII_M[0] * 0.95, 0.5, "perfect match",
            fontsize=8.5, color="#444", va="bottom", ha="left")

    # 50 m marker
    ax.axvline(50, color="#777", linestyle=":", linewidth=1.3, alpha=0.8)
    y_top = ax.get_ylim()[1]
    ax.text(50, y_top * 0.97,
            "50 m\n(current\noperational\ndeme)",
            fontsize=8.5, color="#333", ha="center", va="top")

    ax.set_xscale("log")
    ax.set_xticks(RADII_M)
    ax.set_xticklabels([str(r) for r in RADII_M])
    ax.set_xlabel("Pollinator radius used to build the deme partition  (m, log scale)")
    ax.set_ylabel("Predicted − observed distinct SRK alleles\n(the diversity gap at each radius)")
    ax.set_title(
        "SRK diversity gap vs pollinator radius — EO70 and EO76\n"
        "If gene flow were tighter than 50 m (H1), the gap should\n"
        "shrink as radius decreases. It does not.",
        fontsize=11,
    )
    ax.legend(loc="upper right", fontsize=9, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"step30f_srk_diversity_radius_sweep.{ext}",
                    dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30f] Wrote {FIGURES}/step30f_srk_diversity_radius_sweep.{{png,pdf}}  "
          f"(EO67 omitted from plot: deme count invariant across sweep)")


if __name__ == "__main__":
    main()
