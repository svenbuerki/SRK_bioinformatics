"""bl_constants — BL (Bottleneck Lineage) ordering, colours, and EO lookup.

Mirrored from ~/Repos/SRK_bioinformatics/Scripts/srk_bl_constants.py
(svenbuerki/SRK_bioinformatics) and SRK_BL_integration.py.

Single source of truth for any code in the polyploid model that orders,
colours, or stratifies output by Bottleneck Lineage. Locked Set1 palette
is shared across both projects so figures stay visually consistent.

Resync discipline
-----------------
The two BL summary CSVs (``EO_BL_summary.csv`` and
``EO_group_BL_summary.csv``) live in ``data/salleles/`` and are copied from
the upstream ``tables/`` directory. When upstream refreshes them, copy the
new files over and re-run notebook 00. The BL_ORDER constant is derived
at import time from those CSVs so the ordering tracks the current data.
"""

from __future__ import annotations

import csv
import os
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Paths (auto-discovered from notebook, repo-root, or script contexts)
# ---------------------------------------------------------------------------

def _find_data_file(filename: str) -> str:
    """Locate a data file in ``data/salleles/`` from several typical CWDs.

    Order of search (first hit wins):
      1. notebook-relative: ``../data/salleles/<filename>``
      2. repo-root-relative: ``data/salleles/<filename>``
      3. module-relative: ``<src>/../data/salleles/<filename>``
    Returns the first path that exists; otherwise returns the notebook-
    relative default so the caller gets a useful error message.
    """
    notebook_rel = os.path.join("..", "data", "salleles", filename)
    repo_rel     = os.path.join("data", "salleles", filename)
    module_rel   = os.path.join(os.path.dirname(__file__), "..", "data", "salleles", filename)
    for candidate in (notebook_rel, repo_rel, module_rel):
        if os.path.isfile(candidate):
            return candidate
    return notebook_rel


_DEFAULT_BL_SUMMARY = _find_data_file("EO_BL_summary.csv")
_DEFAULT_EO_GROUP   = _find_data_file("EO_group_BL_summary.csv")


# ---------------------------------------------------------------------------
# Locked palette — RColorBrewer Set1 mapping to BLs by cluster-index.
# Shared with svenbuerki/SRK_bioinformatics and
# svenbuerki/LEPA_EO_spatial_clustering. Do not substitute.
# ---------------------------------------------------------------------------

BL_COLORS: dict[str, str] = {
    "BL1": "#984EA3",  # purple
    "BL2": "#377EB8",  # blue
    "BL3": "#E41A1C",  # red
    "BL4": "#FF7F00",  # orange
    "BL5": "#4DAF4A",  # green
}


# ---------------------------------------------------------------------------
# Manual overrides for Pop codes that do not match an EO label directly.
# Mirrored from SRK_BL_integration.py. Add entries here when more germplasm
# sub-codes are resolved upstream.
# ---------------------------------------------------------------------------

POP_TO_EO_OVERRIDE: dict[str, str] = {
    "15": "EO18",
}


def derive_bl_order(path: str = _DEFAULT_BL_SUMMARY) -> list[str]:
    """Return BLs sorted by total habitat area (desc), within-BL connectivity
    (desc), then BL name (asc). Area is the Ne proxy; connectivity is the
    secondary stratification for BLs of similar size.

    Same sort key as ``srk_bl_constants.py`` upstream and the panel order in
    the BL drift figure produced by LEPA_EO_spatial_clustering.
    """
    path = os.fspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"BL summary file not found: {path} -- copy "
            f"~/Repos/SRK_bioinformatics/tables/EO_BL_summary.csv into "
            f"data/salleles/"
        )
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append((
                r["BL"],
                float(r["total_area_ha"]),
                int(r["n_locations"]) - int(r["n_groups"]),
            ))
    rows.sort(key=lambda x: (-x[1], -x[2], x[0]))
    return [bl for bl, _, _ in rows]


def numeric_order() -> list[str]:
    """Alphanumeric BL1..BL5 for scatter-plot legends."""
    return ["BL1", "BL2", "BL3", "BL4", "BL5"]


def load_eo_to_bl(path: str = _DEFAULT_EO_GROUP) -> dict[str, dict]:
    """Return ``{EO: {BL, Group, Drift_index}}`` keyed by individual EO.

    Composite EO entries like ``"EO118; EO76"`` are split so each EO can be
    looked up independently. Within a single EO, multiple group rows are
    collapsed: ``Group`` becomes a comma-separated list of group IDs and
    ``Drift_index`` is the mean across that EO's groups.
    """
    path = os.fspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"EO group BL summary file not found: {path} -- copy "
            f"~/Repos/SRK_bioinformatics/tables/EO_group_BL_summary.csv "
            f"into data/salleles/"
        )

    rows = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append(r)

    # Split composite EOs; each EO inherits the group's BL and Drift_index.
    by_eo: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        for eo in re.split(r"[;,]", r["EO"]):
            eo = eo.strip()
            if eo:
                by_eo[eo].append(r)

    out: dict[str, dict] = {}
    for eo, group_rows in by_eo.items():
        bls = sorted({g["BL"] for g in group_rows})
        if len(bls) > 1:
            import sys
            sys.stderr.write(f"WARNING: EO {eo} spans multiple BLs {bls}; using first.\n")
        groups = sorted({g["Group"] for g in group_rows}, key=lambda x: int(x))
        dis = [float(g["Drift_index"]) for g in group_rows if g.get("Drift_index")]
        out[eo] = {
            "BL": bls[0],
            "Group": ", ".join(groups),
            "Drift_index": round(sum(dis) / len(dis), 4) if dis else None,
        }
    return out


def normalise_pop_to_eo(pop: str, known_eos: set[str]) -> tuple[str | None, str]:
    """Map a raw ``Pop`` value (from sampling_metadata.csv) to an EO label.

    Returns ``(eo_label_or_None, status)`` where status is one of:
      - ``"Assigned"`` — Pop matched directly to a numeric EO (e.g. ``27 -> EO27``)
      - ``"Inferred"`` — matched via override dict or prefix split (e.g. ``26-3 -> EO26``)
      - ``"Unassigned"`` — no EO match (typically germplasm sub-codes)
    """
    pop = (pop or "").strip()
    if not pop:
        return None, "Unassigned"

    if pop in POP_TO_EO_OVERRIDE:
        eo = POP_TO_EO_OVERRIDE[pop]
        return (eo, "Inferred") if eo in known_eos else (None, "Unassigned")

    def numeric_to_eo(s: str) -> str | None:
        if not s.isdigit():
            return None
        candidate = f"EO{int(s):02d}"
        return candidate if candidate in known_eos else None

    eo = numeric_to_eo(pop)
    if eo:
        return eo, "Assigned"

    if "-" in pop:
        prefix = pop.split("-", 1)[0]
        eo = numeric_to_eo(prefix)
        if eo:
            return eo, "Inferred"

    return None, "Unassigned"


def get_eo_order_within_bl(
    eo_codes: list[str] | None = None,
    bl_summary_path: str = _DEFAULT_EO_GROUP,
) -> list[str]:
    """EO codes ordered by (BL_ORDER, ascending mean Drift_index within BL).

    EO name is the deterministic tie-break so the order matches the R helper
    in upstream ``srk_bl_constants.R``.

    Parameters
    ----------
    eo_codes : list[str] or None
        If provided, restrict the returned list to these EOs (preserving the
        connectivity order). EOs not in the summary file are appended at the
        end sorted lexically.
    bl_summary_path : str
        Path to ``EO_group_BL_summary.csv``.
    """
    eo_to_bl = load_eo_to_bl(bl_summary_path)
    bl_order = derive_bl_order()
    bl_rank = {bl: i for i, bl in enumerate(bl_order)}

    ordered = sorted(
        eo_to_bl.keys(),
        key=lambda e: (
            bl_rank.get(eo_to_bl[e]["BL"], len(bl_order)),
            eo_to_bl[e]["Drift_index"] if eo_to_bl[e]["Drift_index"] is not None else float("inf"),
            e,
        ),
    )

    if eo_codes is not None:
        eo_set = set(eo_codes)
        kept = [e for e in ordered if e in eo_set]
        extras = sorted(eo_set - set(kept))
        ordered = kept + extras
    return ordered


# Eager-derived constants for convenience. Re-import this module if the
# underlying CSVs change.
try:
    BL_ORDER: list[str] = derive_bl_order()
except FileNotFoundError:
    # Allow import in test contexts where CSVs are not present; callers can
    # still use derive_bl_order(path=...) with an explicit path.
    BL_ORDER = []

BL_ORDER_NUMERIC: list[str] = numeric_order()
