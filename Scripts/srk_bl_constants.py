"""srk_bl_constants.py — Python mirror of srk_bl_constants.R.

Single source of truth for the Bottleneck-Lineage (BL) ordering, colour
palette, and locationCode → BL mapping used across the SRK pipeline.
Reads the same `Tables/EO_BL_summary.csv` and `Tables/EO_group_BL_summary.csv`
that the R constants file consumes, so ordering and colours stay in lock-step
with the LEPA_EO_spatial_clustering figures.

Public API
----------
BL_COLORS               dict {BLn: hex colour} — Set1 cluster-index mapping.
BL_ORDER                list of BLn strings, area-then-connectivity DESC.
load_eo_to_bl(...)      map {EOcode: BL} from EO_group_BL_summary.csv
locationCode_to_bl(...) map {locationCode: BL} — locationCode → base EO → BL
                        (locationCodes like "EO18-7" map to "EO18" → its BL).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
DEFAULT_BL_SUMMARY_CSV       = Path("Tables/EO_BL_summary.csv")
DEFAULT_EO_GROUP_BL_CSV      = Path("Tables/EO_group_BL_summary.csv")

# BL → hex colour (Set1 palette, cluster-index mapping).
BL_COLORS: dict[str, str] = {
    "BL1": "#984EA3",   # purple
    "BL2": "#377EB8",   # blue
    "BL3": "#E41A1C",   # red
    "BL4": "#FF7F00",   # orange
    "BL5": "#4DAF4A",   # green
}


def load_bl_order(path: Path = DEFAULT_BL_SUMMARY_CSV) -> list[str]:
    """BL_ORDER: area DESC, connectivity (n_locations − n_groups) DESC,
    BL name ASC. Same key as srk_bl_constants.R."""
    df = pd.read_csv(path)
    df["connectivity"] = df["n_locations"] - df["n_groups"]
    df = df.sort_values(
        ["total_area_ha", "connectivity", "BL"],
        ascending=[False, False, True],
    )
    return df["BL"].tolist()


try:
    BL_ORDER = load_bl_order()
except FileNotFoundError:
    BL_ORDER = ["BL4", "BL5", "BL3", "BL1", "BL2"]   # documented default


# ---------------------------------------------------------------------------
def load_eo_to_bl(path: Path = DEFAULT_EO_GROUP_BL_CSV) -> dict[str, str]:
    """Map each EO code to its BL. Composite EO entries like 'EO118; EO76'
    are split so both codes are returned. If an EO spans multiple BLs
    (uncommon), the first one encountered wins."""
    df = pd.read_csv(path)
    parts = df["EO"].astype(str).str.split(r"[;,]\s*")
    df = df.assign(EO=parts).explode("EO").copy()
    df["EO"] = df["EO"].str.strip()
    eo_to_bl: dict[str, str] = {}
    for _, row in df.iterrows():
        eo = row["EO"]
        if eo and eo not in eo_to_bl:
            eo_to_bl[eo] = row["BL"]
    return eo_to_bl


# Match "EO" + one-or-more digits at the start. Any trailing suffix
# (a dash + subunit, or letters like "RT" for reintroduction / reserve
# treatment) is stripped so subunits inherit their parent EO's BL.
_EO_BASE = re.compile(r"^(EO)(\d+)")


def base_eo(location_code: str) -> str:
    """Return the base EO code matching the LEPA_EO_spatial_clustering
    naming convention (zero-padded to two digits for EO0-EO9).

    Examples::

        'EO18-7'  → 'EO18'
        'EO25-A'  → 'EO25'
        'EO27RT'  → 'EO27'    (reintroduction / reserve-treatment subunit)
        'EO8'     → 'EO08'    (zero-padded to match the BL summary CSV)
        'EO118'   → 'EO118'   (three-digit codes unchanged)

    Returns an empty string when the code does not start with 'EO'.
    """
    m = _EO_BASE.match(str(location_code))
    if not m:
        return ""
    n = int(m.group(2))
    # Zero-pad only the one-digit EOs; the CSV keeps single-digit EOs as
    # 'EO01'..'EO09'. Two- and three-digit EOs are unchanged.
    return f"EO{n:02d}" if n < 10 else f"EO{n}"


def locationCode_to_bl(location_codes,
                       eo_to_bl: dict[str, str] | None = None
                       ) -> pd.Series:
    """Vectorised map locationCode → BL using base_eo() + eo_to_bl."""
    if eo_to_bl is None:
        eo_to_bl = load_eo_to_bl()
    codes = pd.Series(location_codes, dtype=str).reset_index(drop=True)
    bases = codes.map(base_eo)
    return bases.map(eo_to_bl)
