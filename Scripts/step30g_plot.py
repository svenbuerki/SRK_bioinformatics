"""Step 30g companion — across-year prediction scatter.

Reads step30g_across_year_comparison.tsv and plots:

  Panel A — SRK diversity 2025 (x) vs 2026 (y) per population
  Panel B — pollen compatibility 2025 (x) vs 2026 (y) per population

Dots sized by total N_fertile; coloured by stable_across_years flag
(both CIs overlap). 1:1 equality diagonal. Candidate large + small
annotated with populationID.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TABLES  = Path("Tables/Phase5")
FIGURES = Path("figures/Phase5")

STABLE_COLOR   = "#5e3c99"  # purple — same as both-years occupancy
UNSTABLE_COLOR = "#D55E00"  # vermillion


def main() -> None:
    comp = pd.read_csv(TABLES / "step30g_across_year_comparison.tsv",
                        sep="\t", encoding="utf-8-sig")
    stable = comp[comp["stable_across_years"]].sort_values(
        "N_fert_total", ascending=False
    ).reset_index(drop=True)
    large_id = int(stable.iloc[0]["populationID"]) if len(stable) else None
    small_id = int(stable.iloc[-1]["populationID"]) if len(stable) >= 2 else None

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    size_min, size_max = 25, 420
    n_max = max(comp["N_fert_total"].max(), 1)
    sizes = size_min + (size_max - size_min) * np.sqrt(
        comp["N_fert_total"] / n_max
    )

    for i, (ax, metric, lo_col, hi_col, label) in enumerate([
        (axes[0], "diversity", "diversity_CI_overlap",   None,
            "Predicted distinct SRK alleles"),
        (axes[1], "pcompat",   "pcompat_CI_overlap",     None,
            "Predicted pollen compatibility"),
    ]):
        x = comp[f"pred_{metric}_2025"]
        y = comp[f"pred_{metric}_2026"]
        ok = comp[lo_col]

        ax.scatter(x[ok], y[ok], s=sizes[ok], c=STABLE_COLOR, alpha=0.7,
                   edgecolor="white", linewidth=0.6,
                   label=f"stable ({int(ok.sum())})")
        ax.scatter(x[~ok], y[~ok], s=sizes[~ok], c=UNSTABLE_COLOR, alpha=0.7,
                   edgecolor="white", linewidth=0.6,
                   label=f"CI mismatch ({int((~ok).sum())})")

        # 1:1 diagonal
        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        pad = (hi - lo) * 0.08 + 0.001
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                "-", color="#444", linewidth=1.1, alpha=0.6)

        # Annotate candidates
        if large_id is not None:
            r = comp[comp["populationID"] == large_id].iloc[0]
            ax.annotate(f"P{large_id} (LARGE)",
                        xy=(r[f"pred_{metric}_2025"], r[f"pred_{metric}_2026"]),
                        xytext=(8, -4), textcoords="offset points",
                        fontsize=9, color=STABLE_COLOR, fontweight="bold")
        if small_id is not None and small_id != large_id:
            r = comp[comp["populationID"] == small_id].iloc[0]
            ax.annotate(f"P{small_id} (SMALL)",
                        xy=(r[f"pred_{metric}_2025"], r[f"pred_{metric}_2026"]),
                        xytext=(8, -4), textcoords="offset points",
                        fontsize=9, color=STABLE_COLOR, fontweight="bold")

        ax.set_xlabel(f"{label} — 2025")
        ax.set_ylabel(f"{label} — 2026")
        ax.set_title(f"Panel {'A' if i == 0 else 'B'} — "
                      f"{label.lower()} across years",
                      fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="lower right", fontsize=9, frameon=True)

    fig.suptitle(
        "Phase 5 across-year prediction comparison per population "
        f"(both-year populations only; n = {len(comp)})\n"
        "Dot size ∝ total N_fertile. Purple = both 95 % CIs overlap "
        "(stable). Vermillion = CI mismatch (prediction shifted).",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    out_png = FIGURES / "step30g_across_year_scatter.png"
    out_pdf = FIGURES / "step30g_across_year_scatter.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[step30g-plot] Wrote {out_png.name} + .pdf")


if __name__ == "__main__":
    main()
