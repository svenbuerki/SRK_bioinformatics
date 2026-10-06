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
make_location_label(...) str — `{locationCode}_{locationID}` single-site
                        display label. Project-wide convention adopted
                        2026-10-03: every per-location row across every
                        figure uses this label so sites are never
                        confused when a locationCode is shared between
                        multiple locationIDs (e.g. EO8 → EO8_27 /
                        EO8_28 / EO8_29).
location_label_series(...) pd.Series — vectorised form of the above for
                        DataFrames.
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


def make_location_label(location_code: str, location_id: int | str) -> str:
    """Project-wide display label for a single Phase 5 location.

    Returns `{location_code}_{location_id}`.

    Rationale (adopted 2026-10-03). The DB has multiple distinct
    physical sites that share one `locationCode` (EO8 covers four
    locationIDs; EO18-7, EO26-3, EO27-1 each cover four; EO27 covers
    three; EO26 covers two). Those sites are > 500 m apart and the
    Phase 5 pipeline treats each `locationID` as a separate location
    (correct), but the raw `locationCode` is ambiguous for display.
    Joining code + locationID with `_` gives every row across every
    figure (step29d, step30 Figures 3 / 4, step30b, step30c, step30d,
    step30e) a unique, biologically interpretable label.
    """
    return f"{location_code}_{int(location_id)}"


def location_label_series(df: pd.DataFrame,
                           code_col: str = "locationCode",
                           id_col:   str = "locationID") -> pd.Series:
    """Vectorised `make_location_label` for a DataFrame.

    Returns a Series of `{code}_{id}` strings aligned with `df`'s index.
    """
    return (df[code_col].astype(str).str.strip()
            + "_"
            + df[id_col].astype(int).astype(str))


# ---------------------------------------------------------------------------
# NEW BL framework (2026-10-06, Phase 5 § A.4.5) — populationID-based
# ---------------------------------------------------------------------------
# From 2026-10-06 the Bottleneck Lineage (BL) is a POPULATION-level
# property, not an EO property. The population-level BL framework was
# built in step30h_* on the 44 populations defined by step29a (500 m
# connectivity on pooled 2025+2026 events + historical centroids).
#
# Transition status (2026-10-06): both frameworks coexist in the code
# base.
#   - locationCode_to_bl()  — OLD EO-based BL (via Tables/EO_*_summary.csv)
#   - populationID_to_bl()  — NEW population-based BL (via
#                             Tables/Phase5/step30h_populationID_crosswalk.tsv)
# New Phase 5 scripts (step30h_*, step30g_* post-remap) use the new
# BL. Legacy per-locationCode figures (step30, step30c-e) still use
# the EO-based mapping; migrate those scripts individually as needed.
#
# BL_COLORS stay the same 5-colour Set1 palette (BL names BL1..BL5
# are reused with new meanings). The new BL_ORDER is read from
# Tables/Phase5/step30h_bl_definition.tsv (area DESC → connectivity
# DESC rule) and will typically read BL1 → BL5 directly.

DEFAULT_POPULATION_BL_TSV = Path(
    "Tables/Phase5/step30h_bl_definition.tsv"
)
DEFAULT_POPULATION_CROSSWALK_TSV = Path(
    "Tables/Phase5/step30h_populationID_crosswalk.tsv"
)


def load_population_bl_order(
        path: Path = DEFAULT_POPULATION_BL_TSV) -> list[str]:
    """BL_ORDER for the NEW population-based BL framework.

    Reads step30h_bl_definition.tsv; falls back to the hard-coded
    BL1..BL5 if the file does not exist (so this import does not
    break scripts that run before Phase III).
    """
    if not path.exists():
        return ["BL1", "BL2", "BL3", "BL4", "BL5"]
    df = pd.read_csv(path, sep="\t", encoding="utf-8-sig")
    return df.sort_values("BL_rank")["BL"].tolist()


def populationID_to_bl(
        population_ids,
        crosswalk_path: Path = DEFAULT_POPULATION_CROSSWALK_TSV
) -> pd.Series:
    """Vectorised NEW BL lookup: populationID → BL. Reads the Phase III
    crosswalk."""
    if not crosswalk_path.exists():
        return pd.Series([pd.NA] * len(list(population_ids)),
                          dtype="object")
    cw = pd.read_csv(crosswalk_path, sep="\t", encoding="utf-8-sig")
    m = cw.set_index("populationID_new")["BL"].to_dict()
    return pd.Series(population_ids, dtype="Int64").map(m)
