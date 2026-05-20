"""srk_bl_constants — single source of truth for BL ordering and colours.

Imported by every SRK pipeline script that orders or colours bottleneck
lineages so that ordering and palette stay consistent across the project.

Both orderings are DERIVED FROM CSVs mirrored from the
LEPA_EO_spatial_clustering project, so the entire pipeline reorders
automatically when sampling grows — only the CSVs in ``Tables/`` need
updating.

Exports
-------
BL_ORDER : list[str]
    BLs sorted high-to-low by within-BL connectivity, computed at import
    time from ``Tables/EO_BL_summary.csv`` using the formula
        connectivity = n_locations - n_groups        (primary, desc)
        total_pop_size                               (tie-break, desc)
    Same formula used to order the panels in the BL drift figure
    produced by the spatial-clustering project.

BL_COLORS : dict[str, str]
    Set1 palette mapped to BL by cluster-index, matching the rendered
    figures in LEPA_EO_spatial_clustering (dendrogram, BL geographic
    context map, drift panel, fragmentation network). Verified against
    live ``bl_strip_cols`` on 2026-05-18.

get_eo_order_within_bl(eo_codes=None, bl_summary_path=...)
    EO codes ordered by (BL_ORDER, ascending mean Drift_index within BL).
    Reads ``Tables/EO_group_BL_summary.csv``. EO name is the tie-break so
    the order is deterministic and matches the R helper.
"""

from pathlib import Path
import csv

_BL_SUMMARY_FILE     = Path("Tables/EO_BL_summary.csv")
_EO_GROUP_BL_SUMMARY = Path("Tables/EO_group_BL_summary.csv")

# BL -> hex colour (cluster-index Set1 mapping)
BL_COLORS = {
    "BL1": "#984EA3",  # purple
    "BL2": "#377EB8",  # blue
    "BL3": "#E41A1C",  # red
    "BL4": "#FF7F00",  # orange
    "BL5": "#4DAF4A",  # green
}


def _derive_bl_order(path=_BL_SUMMARY_FILE):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"BL summary file not found: {path} — mirror it from "
            f"LEPA_EO_spatial_clustering/data/EO_BL_summary.csv"
        )
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            rows.append((
                r["BL"],
                int(r["n_locations"]) - int(r["n_groups"]),  # connectivity
                int(r["total_pop_size"]),
            ))
    rows.sort(key=lambda x: (-x[1], -x[2], x[0]))
    return [bl for bl, _conn, _pop in rows]


BL_ORDER = _derive_bl_order()

# Numerical (alphanumeric) BL order — used for SCATTER-PLOT LEGENDS only
# (TP1, TP2, GFS scatters, allele-accumulation curve legend, etc.) where
# the x/y axes are not BL and the legend just needs an easy-to-read order.
# Pass to plotting code (e.g. matplotlib legend handles, ggplot breaks=...)
# wherever the BL appears in a colour legend.  Categorical axes / facets /
# tables should still use BL_ORDER (connectivity-driven).
BL_ORDER_NUMERIC = sorted(BL_ORDER)


def get_eo_order_within_bl(eo_codes=None, bl_summary_path=_EO_GROUP_BL_SUMMARY):
    """EO codes ordered by (BL_ORDER, ascending mean Drift_index within BL).

    Parameters
    ----------
    eo_codes : iterable[str] or None
        If provided, restrict the returned list to these EOs (preserving the
        connectivity order). Any extras not in the summary file are appended
        lexically at the end so downstream code is robust to germplasm
        sub-codes etc.
    bl_summary_path : str or pathlib.Path
        Path to EO_group_BL_summary.csv (mirrored under Tables/ in this repo).

    Returns
    -------
    list[str]
    """
    path = Path(bl_summary_path)
    if not path.is_file():
        raise FileNotFoundError(f"EO group BL summary file not found: {path}")

    import re
    sums = {}
    counts = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            # Composite EO entries like "EO118; EO76" represent one geographic
            # group spanning two EOs — split so downstream code sees both
            # codes (each inherits the group's mean Drift_index).
            eos = [e.strip() for e in re.split(r"[;,]", row["EO"]) if e.strip()]
            di = float(row["Drift_index"])
            for eo in eos:
                key = (eo, row["BL"])
                sums[key] = sums.get(key, 0.0) + di
                counts[key] = counts.get(key, 0) + 1

    means = {key: sums[key] / counts[key] for key in sums}
    bl_rank = {bl: i for i, bl in enumerate(BL_ORDER)}
    # Secondary sort on EO name keeps ties deterministic and identical to R.
    ordered = sorted(
        means.keys(),
        key=lambda k: (bl_rank.get(k[1], len(BL_ORDER)), means[k], k[0]),
    )
    ordered_eos = [eo for eo, _bl in ordered]

    if eo_codes is not None:
        eo_set = set(eo_codes)
        kept = [eo for eo in ordered_eos if eo in eo_set]
        extras = sorted(eo_set - set(kept))
        ordered_eos = kept + extras
    return ordered_eos
