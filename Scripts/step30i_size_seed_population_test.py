"""Step 30i — size → seed-yield test at the Phase 5 **population** level.

Adapted from the Field-work protocol's `size_seed_model.py` / `size_matched_
comparison.py` / `per_eo_trends.py` family (ported from EO level to the
Phase 5 populationID level, 2026-10-07).

Goal
----
Behavioural validation of the Phase A per-population pollen-compatibility
prediction (Figure 11 of the compact doc). If a mother plant's observed
seed yield is systematically **less than what her plant size predicts**,
and if that shortfall is concentrated in populations whose Phase A
prediction says they should be mate-limited, we have a direct seed-set
signal of the mate-limitation pathway — **without needing SRK genotypes
yet**. Each plant is its own control via the species-wide size → yield
allometry.

Key change from the Field-work protocol version
-----------------------------------------------
The EO-level version applied a `min-n = 15 mother plants` threshold that
discarded small populations. **Small populations are the ones we most
want to see** (user instruction 2026-10-07): they are where drift is
expected to be strongest and where the Phase A P_compat prediction is
most informative. This Phase 5 version **keeps every population with
≥ 2 plants** (minimum for a CI) and **does not discard any**.

Pipeline
--------
STAGE 1  Predictor selection — 10-fold cross-validated RMSE on log10
         seed yield among candidate allometries:
             height / crown / area (ellipse) / crown + height
         Winner = lowest mean CV-RMSE.
STAGE 2  Expectation model — refit the winner on all plants (both
         years pooled); fitted value = expected log10(yield) for a
         plant of that size. Residual = observed − expected.
STAGE 3  Below-expectation detection per (populationID, year):
             - mean residual = log10 fold-of-expectation (0 = on curve)
             - 95 % CI from the per-population residual spread
             - one-sample t-test (where n ≥ 3) or Wilcoxon (n ≥ 6)
             - Benjamini–Hochberg FDR across all (populationID × year)
               combinations with n ≥ 3
             - Flag "below expectation" (seed shortfall) if FDR-sig and
               mean_resid < 0
STAGE 4  Outputs:
             tables/Phase5/step30i_size_seed_population_strata.tsv
             figures/Phase5/step30i_size_seed_population.png / .pdf

Conventions
-----------
BL colour palette + BL ordering (BL1 → BL5, area DESC → connectivity
DESC) match the rest of Phase 5 — [[feedback_bl_order_area_first]] and
the new population-level BL framework (§ Bottleneck Lineages of the
compact doc). Figures follow the project convention: **2025 = open
circle, 2026 = filled circle**, same BL colour, per-population row.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from sklearn.model_selection import KFold
from statsmodels.stats.multitest import fdrcorrection

from srk_bl_constants import BL_COLORS

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

DB_PATH = Path(
    "/Users/sven/Documents/Current_projects/LEPA_fieldwork_protocol/"
    "SQL_DB/LEPA_SQL.db"
)

NEW_BL_ORDER = ["BL1", "BL2", "BL3", "BL4", "BL5"]
CANDIDATES = {
    "height":       "logy ~ logh",
    "crown":        "logy ~ logc",
    "area":         "logy ~ logarea",
    "crown+height": "logy ~ logc + logh",
}
N_FOLDS = 10
CV_SEED = 1


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def _year_from_date(s: pd.Series) -> pd.Series:
    """occurrenceDate is `MM-DD-YYYY` for 2025; the 2026 'MM-DD' rows
    without a year (field-note data issue, events 660-724) are counted
    as 2026. Returns Int64 year or NaN."""
    s = s.astype(str)
    years = pd.to_numeric(s.str[-4:], errors="coerce")
    short_mask = (s.str.len() == 5) & years.isna()
    years = years.where(~short_mask, 2026)
    return years.astype("Int64")


def load_occurrences() -> pd.DataFrame:
    """Pull per-occurrence crown + height + seed yield + eventID."""
    con = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql_query(
            """
            SELECT occurrenceID,
                   eventID,
                   occurrenceDate,
                   occurrenceHeight    AS h,
                   occurrenceCrownSize AS c,
                   seedQuantityTotal   AS y
            FROM vOccurrenceTraits
            WHERE occurrenceHeight    > 0
              AND occurrenceCrownSize > 0
              AND seedQuantityTotal   > 0
            """,
            con,
        )
    finally:
        con.close()
    df["year"] = _year_from_date(df["occurrenceDate"])
    df = df.dropna(subset=["year"]).copy()
    df["year"] = df["year"].astype(int)
    df["area"] = np.pi / 4.0 * df["c"] * df["h"]
    for col in ("h", "c", "area", "y"):
        df[f"log{col[0] if col != 'area' else 'area'}"] = np.log10(df[col])
    df["logy"] = np.log10(df["y"])
    return df


def load_population_crosswalk() -> pd.DataFrame:
    """populationID per (event_year, eventID). BL is pulled from the
    classified TSV instead of the crosswalk because the crosswalk's
    BL column was NOT updated when the Ward cluster → BL renumbering
    (step30h_define_bl.py) happened — see
    `step30g_populations_classified.tsv` for the authoritative BL."""
    cw = pd.read_csv(TABLES / "step29a_population_crosswalk.tsv",
                     sep="\t", encoding="utf-8-sig")
    return cw[["event_year", "eventID", "populationID"]].drop_duplicates()


def load_population_meta() -> pd.DataFrame:
    """populationID → BL + population_label (authoritative source)."""
    cl = pd.read_csv(TABLES / "step30g_populations_classified.tsv",
                     sep="\t", encoding="utf-8-sig")
    return (cl[["populationID", "BL", "population_label"]]
              .drop_duplicates())


# ---------------------------------------------------------------------------
# Predictor selection
# ---------------------------------------------------------------------------
def cv_rmse(df: pd.DataFrame, formula: str,
            folds: int, seed: int) -> tuple[float, float]:
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    err = []
    for tr, te in kf.split(df):
        m = smf.ols(formula, df.iloc[tr]).fit()
        pred = m.predict(df.iloc[te])
        err.append(np.sqrt(np.mean((df.iloc[te].logy.values - pred.values) ** 2)))
    return float(np.mean(err)), float(np.std(err))


def pick_best_predictor(df: pd.DataFrame) -> tuple[str, str, pd.DataFrame]:
    rows = []
    for name, f in CANDIDATES.items():
        m = smf.ols(f, df).fit()
        rmse, sd = cv_rmse(df, f, N_FOLDS, CV_SEED)
        rows.append((name, rmse, sd, m.rsquared, m.aic))
    rows.sort(key=lambda r: r[1])
    table = pd.DataFrame(rows, columns=["predictor", "cv_rmse", "cv_sd",
                                          "r2_in_sample", "aic"])
    best = rows[0][0]
    return best, CANDIDATES[best], table


# ---------------------------------------------------------------------------
# Per-population test
# ---------------------------------------------------------------------------
def per_population_year_test(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (pid, yr), g in df.groupby(["populationID", "year"]):
        n = len(g)
        r = g["resid"].to_numpy()
        bl = str(g["BL"].iloc[0])
        label = str(g["population_label"].iloc[0])
        if n >= 3:
            mean = float(r.mean())
            sd = float(r.std(ddof=1))
            ci = float(stats.t.ppf(0.975, n - 1) * sd / np.sqrt(n))
            t, p_t = stats.ttest_1samp(r, 0.0)
            p_t = float(p_t)
        else:
            mean = float(r.mean()) if n >= 1 else np.nan
            ci   = np.nan
            t    = np.nan
            p_t  = np.nan
        if n >= 6:
            try:
                _, p_w = stats.wilcoxon(r)
                p_w = float(p_w)
            except ValueError:
                p_w = np.nan
        else:
            p_w = np.nan
        out.append(dict(populationID=pid, year=int(yr), BL=bl,
                        population_label=label, n_plants=int(n),
                        mean_resid=mean, ci95=ci, t_stat=float(t)
                        if not np.isnan(t) else np.nan,
                        p_ttest=p_t, p_wilcoxon=p_w))
    res = pd.DataFrame(out)
    # BH-FDR across all (populationID, year) rows with n ≥ 3
    mask = res["n_plants"] >= 3
    res["q_ttest"] = np.nan
    if mask.sum() >= 2:
        q = fdrcorrection(res.loc[mask, "p_ttest"].fillna(1.0).to_numpy())[1]
        res.loc[mask, "q_ttest"] = q
    res["fold_of_expectation"] = 10 ** res["mean_resid"]
    res["fold_lo95"] = 10 ** (res["mean_resid"] - res["ci95"])
    res["fold_hi95"] = 10 ** (res["mean_resid"] + res["ci95"])
    res["flag_below_expectation"] = ((res["q_ttest"].fillna(1.0) < 0.05)
                                       & (res["mean_resid"] < 0))
    res["flag_above_expectation"] = ((res["q_ttest"].fillna(1.0) < 0.05)
                                       & (res["mean_resid"] > 0))

    # Conservation-priority bucket per (populationID, year). The test
    # separates populations into **buffered** (reproduction is near
    # or above the species expectation — mate availability holds) vs
    # **mate-limited** (reproduction is below size expectation —
    # candidate for prioritization). Four tiers so the "which
    # populations to intervene on" question has a one-shot answer.
    def _bucket(r):
        if r["n_plants"] < 3:
            return "insufficient_data"
        if r["flag_below_expectation"]:
            return "mate_limited"
        if r["flag_above_expectation"]:
            return "buffered_surplus"
        if r["fold_of_expectation"] < 0.8:
            # Trend below but not FDR-sig — watch
            return "candidate_mate_limited"
        return "buffered"

    res["conservation_bucket"] = res.apply(_bucket, axis=1)
    return res


def summarise_per_population(res: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-year rows into one row per populationID with the
    worst-case conservation bucket across years (conservative — if
    EITHER year is flagged mate-limited the population is flagged).
    Useful as the one-shot 'which populations to prioritize' table."""
    PRIORITY = {
        "mate_limited":            0,  # highest priority
        "candidate_mate_limited":  1,
        "insufficient_data":       2,
        "buffered":                3,
        "buffered_surplus":        4,  # lowest priority
    }
    rows = []
    for pid, g in res.groupby("populationID"):
        # Worst bucket across years (lowest PRIORITY score)
        buckets = g["conservation_bucket"].tolist()
        worst = min(buckets, key=lambda b: PRIORITY.get(b, 99))
        years_flagged_below = sorted(
            int(y) for y in g.loc[g["flag_below_expectation"], "year"]
        )
        years_flagged_above = sorted(
            int(y) for y in g.loc[g["flag_above_expectation"], "year"]
        )
        bl = str(g["BL"].iloc[0])
        label = str(g["population_label"].iloc[0])
        folds_by_year = {int(y): float(f) for y, f in
                           zip(g["year"], g["fold_of_expectation"])}
        rows.append(dict(
            populationID=int(pid),
            BL=bl,
            population_label=label,
            conservation_priority=worst,
            priority_rank=PRIORITY.get(worst, 99),
            years_observed=sorted(int(y) for y in g["year"]),
            years_sig_below=years_flagged_below,
            years_sig_above=years_flagged_above,
            fold_2025=folds_by_year.get(2025, np.nan),
            fold_2026=folds_by_year.get(2026, np.nan),
            n_plants_2025=int(g.loc[g["year"] == 2025, "n_plants"].sum())
                           if (g["year"] == 2025).any() else 0,
            n_plants_2026=int(g.loc[g["year"] == 2026, "n_plants"].sum())
                           if (g["year"] == 2026).any() else 0,
        ))
    out = pd.DataFrame(rows).sort_values(["priority_rank", "populationID"])
    return out


# ---------------------------------------------------------------------------
# Figure A — species-wide calibration scatter
# ---------------------------------------------------------------------------
def plot_species_calibration(df: pd.DataFrame, gm,
                               best: str, cv_table: pd.DataFrame,
                               out_png: Path, out_pdf: Path) -> None:
    """Panel A of the original Field-work-protocol figure, ported to
    the Phase 5 population frame: observed seed yield vs the
    species-wide allometric expectation. One dot per plant, coloured
    by its population's BL. **This is the 'across-species' half of
    the test** — the species-wide curve that every per-population
    row in Figure 17b is measured against."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.0, 7.0),
                                     gridspec_kw={"width_ratios": [3.0, 2.0]})

    bls_present = [b for b in NEW_BL_ORDER if (df["BL"] == b).any()]
    expected = 10 ** df["expected"].to_numpy()
    observed = df["y"].to_numpy()

    # Left panel: species calibration scatter
    for bl in bls_present:
        g = df[df["BL"] == bl]
        col = BL_COLORS.get(bl, "#777777")
        axL.scatter(10 ** g["expected"], g["y"],
                     s=16, alpha=0.55, color=col,
                     edgecolor="white", linewidth=0.3,
                     label=f"{bl}  (n = {len(g)} plants)")
    lim_lo = max(min(observed.min(), expected.min()) * 0.6, 0.3)
    lim_hi = max(observed.max(), expected.max()) * 1.4

    # Species curve = 1:1 — black dashed, prominent
    axL.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
              color="black", linestyle="--", lw=1.6, alpha=0.95,
              label="y = x · on the species curve  (fold = 1)",
              zorder=1.5)
    # ½× — SHORTFALL side (observed is HALF of expected): red,
    # dash-dot. Sits BELOW the 1:1 line.
    axL.plot([lim_lo, lim_hi], [lim_lo * 0.5, lim_hi * 0.5],
              color="#b2182b", linestyle=(0, (5, 2, 1, 2)), lw=1.3,
              alpha=0.85,
              label="y = x / 2 · ½ × expectation  (shortfall band)",
              zorder=1.4)
    # 2× — SURPLUS side (observed is DOUBLE of expected): green,
    # densely dotted. Sits ABOVE the 1:1 line.
    axL.plot([lim_lo, lim_hi], [lim_lo * 2.0, lim_hi * 2.0],
              color="#1b7837", linestyle=(0, (1, 2)), lw=1.3,
              alpha=0.85,
              label="y = 2 x · 2 × expectation  (surplus band)",
              zorder=1.4)

    # End-of-line labels at the right edge so the three lines are
    # self-explanatory even without the legend.
    for yfrac, text, col in [
        (lim_hi,       "y = x",      "black"),
        (lim_hi * 0.5, "y = x / 2",  "#b2182b"),
        (lim_hi * 2.0, "y = 2x",     "#1b7837"),
    ]:
        if yfrac <= lim_hi * 1.05 and yfrac >= lim_lo:
            axL.text(lim_hi * 0.96, yfrac, "  " + text,
                      color=col, fontsize=9, fontweight="bold",
                      ha="left", va="center",
                      bbox=dict(facecolor="white", edgecolor="none",
                                 alpha=0.85, pad=1.0))
    axL.set_xscale("log"); axL.set_yscale("log")
    axL.set_xlim(lim_lo, lim_hi); axL.set_ylim(lim_lo, lim_hi)
    axL.set_xlabel("Expected seed yield per plant  (from the species-wide "
                    f"{best} model, log scale)", fontsize=10)
    axL.set_ylabel("Observed seed yield per plant  (log scale)", fontsize=10)
    axL.legend(loc="lower right", fontsize=8, frameon=True, ncol=1)
    axL.spines["top"].set_visible(False); axL.spines["right"].set_visible(False)
    axL.grid(True, which="both", alpha=0.3)
    # Panel label as a bold text in the upper-left corner of the
    # data area (NOT via set_title, which collides with the suptitle
    # at the top of the figure).
    axL.text(0.02, 0.985,
              "A. Species-wide calibration scatter",
              transform=axL.transAxes, fontsize=11.5, va="top", ha="left",
              fontweight="bold")
    axL.text(0.02, 0.945,
              f"R² = {gm.rsquared:.2f}  ·  one dot per plant  "
              f"(n = {len(df)} plants, {df['populationID'].nunique()} "
              f"populations, {df['year'].nunique()} years)",
              transform=axL.transAxes, fontsize=8.5, va="top", ha="left",
              color="#444")

    # Right panel: two clearly separated blocks on a plain
    # background. Panel letter + heading sit INSIDE the data area
    # (not via set_title, which was colliding with the suptitle).
    axR.axis("off")

    # Panel label at the very top
    axR.text(0.02, 0.985,
              "B. Model build  —  how the species-wide expectation is set",
              transform=axR.transAxes, fontsize=11.5, va="top", ha="left",
              fontweight="bold")

    # --- TOP block: STAGE 1 predictor-selection CV table -----------
    axR.text(0.02, 0.90,
              "STAGE 1  ·  predictor selection (10-fold cross-validation)",
              transform=axR.transAxes, fontsize=10.5, va="top", ha="left",
              fontweight="bold")
    axR.text(0.02, 0.855,
              "Lower CV-RMSE = better out-of-sample fit  →  wins Stage 2.",
              transform=axR.transAxes, fontsize=8.5, va="top", ha="left",
              color="#555")
    cv_hdr = f"{'predictor':<13}{'CV-RMSE':>9}{'±sd':>7}{'R²_in':>7}{'AIC':>9}"
    cv_lines = [cv_hdr]
    for _, r in cv_table.iterrows():
        marker = "  ★" if r["predictor"] == best else "   "
        cv_lines.append(
            f"{r['predictor']:<13}"
            f"{r['cv_rmse']:>9.3f}"
            f"{r['cv_sd']:>7.3f}"
            f"{r['r2_in_sample']:>7.3f}"
            f"{r['aic']:>9.1f}"
            f"{marker}"
        )
    axR.text(0.02, 0.80, "\n".join(cv_lines), transform=axR.transAxes,
              fontsize=9.5, va="top", ha="left", family="monospace")
    axR.text(0.02, 0.53,
              f"→  Winner: {best}  (★)",
              transform=axR.transAxes, fontsize=10, va="top", ha="left",
              color="#1b7837", fontweight="bold")

    # Thin horizontal separator between the two blocks
    axR.axhline(0.47, color="#aaa", lw=0.7, alpha=0.8)

    # --- BOTTOM block: STAGE 2 refitted expectation model ----------
    axR.text(0.02, 0.42,
              "STAGE 2  ·  refitted expectation model (all plants pooled)",
              transform=axR.transAxes, fontsize=10.5, va="top", ha="left",
              fontweight="bold")
    axR.text(0.02, 0.375,
              "Fitted value = expected log₁₀(yield) for a plant of that "
              "size.\nResidual = observed − expected.",
              transform=axR.transAxes, fontsize=8.5, va="top", ha="left",
              color="#555")

    coef_bits = ",   ".join(
        f"{k} = {v:+.3f}" for k, v in gm.params.to_dict().items()
    )
    sd_log = float(np.sqrt(gm.scale))
    stats_lines = [
        f"Formula       :  log₁₀(yield) ~ {CANDIDATES[best].split('~',1)[1].strip()}",
        f"Coefficients  :  {coef_bits}",
        f"R² (in sample):  {gm.rsquared:.3f}",
        f"Residual SD   :  {sd_log:.3f} log₁₀   (≈ × / ÷ {10**sd_log:.1f})",
    ]
    axR.text(0.02, 0.27, "\n".join(stats_lines), transform=axR.transAxes,
              fontsize=9.5, va="top", ha="left", family="monospace")
    # Footer
    axR.text(0.02, 0.03,
              "Every row of Figure 17b is measured against this expectation.",
              transform=axR.transAxes, fontsize=8.5, va="bottom", ha="left",
              color="#333", style="italic")

    fig.suptitle(
        "Phase 5 § C.0.d · Species-wide size → seed-yield expectation  "
        "(across-species build; the per-population test is Figure 17b)",
        fontsize=12.5, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30i] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# Figure B — forest plot per BL, 2025 open vs 2026 filled per population
# ---------------------------------------------------------------------------
def plot_forest(res: pd.DataFrame, df: pd.DataFrame,
                 out_png: Path, out_pdf: Path) -> None:
    # Ordering: within each BL, by populationID (same as other per-pop figures)
    bls = [b for b in NEW_BL_ORDER if (res["BL"] == b).any()]
    heights = []
    for b in bls:
        n_pops = res.loc[res["BL"] == b, "populationID"].nunique()
        heights.append(max(int(n_pops), 1))

    fig, axes = plt.subplots(
        len(bls), 1,
        figsize=(10.5, max(6.5, 0.33 * sum(heights) + 2.5)),
        gridspec_kw={"height_ratios": heights},
        sharex=True,
    )
    if len(bls) == 1:
        axes = [axes]

    for ax, bl in zip(axes, bls):
        sub = (res[res["BL"] == bl]
                 .sort_values(["populationID", "year"])
                 .reset_index(drop=True))
        pops = sorted(sub["populationID"].unique())
        y_pos = {pid: i for i, pid in enumerate(pops)}
        colour = BL_COLORS.get(bl, "#777777")

        # Reference line at fold = 1 (on the species curve)
        ax.axvline(1.0, color="#444", ls="--", lw=1.0, alpha=0.9)
        # 0.5× and 2× light reference bands
        ax.axvline(0.5, color="#888", ls=":", lw=0.7, alpha=0.6)
        ax.axvline(2.0, color="#888", ls=":", lw=0.7, alpha=0.6)

        for _, row in sub.iterrows():
            pid = int(row["populationID"])
            yr  = int(row["year"])
            yi  = y_pos[pid] + (-0.14 if yr == 2025 else 0.14)
            fold = row["fold_of_expectation"]
            lo   = row["fold_lo95"]
            hi   = row["fold_hi95"]
            # Error bar if CI exists
            if not np.isnan(lo) and not np.isnan(hi):
                ax.plot([lo, hi], [yi, yi], color=colour, lw=1.3,
                        alpha=0.55, zorder=1)
            # Marker — open 2025 / filled 2026
            n = int(row["n_plants"])
            # Dot size scaled by n_plants (sqrt), floor at 30
            s = 30 + 4.5 * np.sqrt(max(n, 1))
            mfc = ("white" if yr == 2025 else colour)
            mew = (1.4 if yr == 2025 else 0.6)
            # Halo ring BEHIND the dot ties the FDR-flag to the
            # specific (population, year) data point and keeps the
            # open/filled fill visible. Red = below expectation
            # (seed shortfall → mate limitation surfacing);
            # green = above expectation (seed surplus → population
            # is mate-replete). Both halos are predictions tested
            # against observed seed yield (user feedback 2026-10-08).
            if bool(row["flag_below_expectation"]):
                ax.scatter([fold], [yi], s=s * 2.3,
                            facecolor="none", edgecolor="#b2182b",
                            linewidth=2.0, zorder=2)
            elif bool(row["flag_above_expectation"]):
                ax.scatter([fold], [yi], s=s * 2.3,
                            facecolor="none", edgecolor="#1b7837",
                            linewidth=2.0, zorder=2)
            ax.scatter([fold], [yi], s=s, facecolor=mfc,
                        edgecolor=colour, linewidth=mew, zorder=3)

        # Left margin: P{N} identifier
        ax.set_yticks(list(y_pos.values()))
        ax.set_yticklabels([f"P{pid:>2}" for pid in pops], fontsize=9)

        # Finalise ax limits + invert BEFORE building the twinx right
        # margin, so the right-margin ylim is in sync with the plotted
        # points. (If ax2 is set up first and ax.invert_yaxis() is
        # called afterwards, ax2's labels end up at the wrong vertical
        # positions — bug fix 2026-10-08.)
        ax.set_xscale("log")
        ax.set_xlim(0.1, 10.0)
        ax.set_ylim(-0.7, len(pops) - 0.3)
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="x", which="both", alpha=0.3)

        # Right margin: per-year n_plants annotation, aligned to the
        # same vertical positions as the plotted markers.
        right = []
        for pid in pops:
            parts = []
            for yr in (2025, 2026):
                r = sub[(sub["populationID"] == pid) & (sub["year"] == yr)]
                if not r.empty:
                    parts.append(f"{yr} n={int(r['n_plants'].iat[0])}")
            right.append("  ·  ".join(parts))
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(list(y_pos.values()))
        ax2.set_yticklabels(right, fontsize=7.5, color="#444")
        ax2.tick_params(axis="y", length=0, pad=2)
        for s in ("top", "right", "left"):
            ax2.spines[s].set_visible(False)
        # BL badge
        ax.text(0.008, 0.97, bl,
                transform=ax.transAxes,
                fontsize=12, fontweight="bold", color=colour,
                va="top", ha="left",
                bbox=dict(facecolor="white", edgecolor=colour,
                           boxstyle="round,pad=0.25", alpha=0.9,
                           linewidth=1.0))

    axes[-1].set_xlabel(
        "Fold of expectation  =  observed yield ÷ expected from size "
        "(log scale)   ·   < 1 = seed shortfall  ·  > 1 = seed surplus"
        "\nred halo = FDR < 0.05 BELOW   ·   green halo = FDR < 0.05 "
        "ABOVE   (both for that specific population × year)",
        fontsize=10,
    )

    # Figure-level legend above subplots
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="white", markeredgecolor="#444",
                markeredgewidth=1.4, markersize=8, label="2025 (open)"),
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="#444", markeredgecolor="white",
                markeredgewidth=0.5, markersize=8, label="2026 (filled)"),
        Line2D([0], [0], color="#444", ls="--", lw=1.0,
                label="on species curve  (fold = 1)"),
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="none", markeredgecolor="#b2182b",
                markeredgewidth=2.0, markersize=13,
                label="red halo = FDR < 0.05 BELOW  (seed shortfall)"),
        Line2D([0], [0], marker="o", linestyle="",
                markerfacecolor="none", markeredgecolor="#1b7837",
                markeredgewidth=2.0, markersize=13,
                label="green halo = FDR < 0.05 ABOVE  (seed surplus)"),
    ]
    fig.legend(handles=legend_handles,
                loc="upper center", bbox_to_anchor=(0.5, 0.965),
                fontsize=9, frameon=True, ncol=4)

    fig.suptitle(
        "Phase 5 — per-population observed seed yield vs size expectation  "
        f"(all populations with ≥ 2 plants, n = {res['populationID'].nunique()} "
        f"populations across {df['year'].nunique()} years)",
        fontsize=12, y=0.998,
    )
    fig.tight_layout(rect=[0, 0, 0.90, 0.93])
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30i] Wrote {out_png.name} + .pdf")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    print("[step30i] Loading occurrences from LEPA DB …")
    occ = load_occurrences()
    print(f"[step30i]   n = {len(occ)} occurrences with (height, crown, seed yield)")

    cw = load_population_crosswalk()
    merged = occ.merge(cw, left_on=["year", "eventID"],
                        right_on=["event_year", "eventID"], how="inner")
    print(f"[step30i]   n = {len(merged)} occurrences joinable to a Phase 5 populationID")

    meta = load_population_meta()
    merged = merged.merge(meta, on="populationID", how="left")
    # Populations present in the crosswalk but missing from the
    # classified TSV would land with NaN BL — flag + drop so stats are
    # not polluted. Should be zero when both outputs are in sync.
    n_missing_bl = int(merged["BL"].isna().sum())
    if n_missing_bl:
        print(f"[step30i]   WARNING: {n_missing_bl} occurrences have no BL "
              f"assignment (populationID not in classified TSV) — dropped.")
        merged = merged.dropna(subset=["BL"]).copy()

    # ---- STAGE 1 — predictor selection ----
    print("[step30i] STAGE 1 — predictor selection (10-fold CV)")
    best, best_formula, cv_table = pick_best_predictor(merged)
    print(cv_table.to_string(index=False))
    print(f"[step30i]   -> selected predictor: {best}  ({best_formula})")

    # ---- STAGE 2 — global expectation model ----
    gm = smf.ols(best_formula, merged).fit()
    merged["expected"] = gm.fittedvalues
    merged["resid"]    = merged["logy"] - merged["expected"]
    print(f"[step30i] STAGE 2 — global expectation model: "
          f"R² = {gm.rsquared:.3f}, "
          f"residual SD = {np.sqrt(gm.scale):.3f} log10 "
          f"(≈ × / ÷ {10**np.sqrt(gm.scale):.1f})")

    # ---- STAGE 3 — per-population test ----
    res = per_population_year_test(merged)
    # BL ordering + populationID within-BL ordering for display
    res["_bl_rank"] = res["BL"].map(
        {b: i for i, b in enumerate(NEW_BL_ORDER)}).fillna(99).astype(int)
    res = res.sort_values(["_bl_rank", "populationID", "year"])
    print(f"[step30i] STAGE 3 — per-population test "
          f"({res['populationID'].nunique()} populations across both years, "
          f"min-n = 2)")
    print(res[["populationID", "year", "BL", "n_plants",
                "fold_of_expectation", "fold_lo95", "fold_hi95",
                "p_ttest", "q_ttest", "flag_below_expectation"]]
          .to_string(index=False))

    flagged_below = res[res["flag_below_expectation"]]
    flagged_above = res[res["flag_above_expectation"]]
    print(f"[step30i] Below-expectation (FDR < 0.05 AND mean_resid < 0): "
          f"{len(flagged_below)} (population, year) rows")
    for _, r in flagged_below.iterrows():
        print(f"           P{int(r['populationID']):>2} {r['BL']} "
              f"{int(r['year'])}  fold = {r['fold_of_expectation']:.2f}  "
              f"n = {int(r['n_plants'])}  "
              f"({r['population_label']})")
    print(f"[step30i] Above-expectation (FDR < 0.05 AND mean_resid > 0): "
          f"{len(flagged_above)} (population, year) rows")

    # ---- STAGE 4 — outputs ----
    out_tsv = TABLES / "step30i_size_seed_population_strata.tsv"
    res.drop(columns=["_bl_rank"]).to_csv(out_tsv, sep="\t", index=False)
    print(f"[step30i] Strata table → {out_tsv}")

    # Keep the predictor-selection CV table for the methods section
    cv_table.to_csv(TABLES / "step30i_size_seed_cv_rmse.tsv",
                     sep="\t", index=False)

    # Per-population aggregated conservation priority (worst-case
    # bucket across years → one row per populationID). This is the
    # one-shot "which populations to prioritize" answer.
    prio = summarise_per_population(res)
    out_prio = TABLES / "step30i_conservation_priority.tsv"
    prio.to_csv(out_prio, sep="\t", index=False)
    print(f"[step30i] Conservation priority (per population) → {out_prio}")
    print(f"[step30i] Priority counts:")
    for bucket in ("mate_limited", "candidate_mate_limited",
                     "insufficient_data", "buffered", "buffered_surplus"):
        n = int((prio["conservation_priority"] == bucket).sum())
        if n:
            print(f"           {bucket:<28} {n}")

    FIGURES.mkdir(parents=True, exist_ok=True)
    # Figure 17a — species-wide calibration (across-species evidence)
    plot_species_calibration(
        merged, gm, best, cv_table,
        FIGURES / "step30i_species_calibration.png",
        FIGURES / "step30i_species_calibration.pdf",
    )
    # Figure 17b — per-population forest (within-population test)
    plot_forest(res, merged,
                FIGURES / "step30i_size_seed_population.png",
                FIGURES / "step30i_size_seed_population.pdf")


if __name__ == "__main__":
    main()
