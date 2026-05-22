#!/usr/bin/env Rscript
# =============================================================================
# SRK_allele_eulerr_BLs.R — Step 18 companion
# Area-proportional 5-set Euler diagram of S-allele sets per Bottleneck Lineage
# =============================================================================
#
# Companion to SRK_allele_sharing_EOs.py (the BL UpSet plot). UpSet is the
# information-dense view; eulerr is the stakeholder-friendly view that shows
# the same data with area-proportional ellipses. Regenerate alongside the
# UpSet plot whenever new samples are added.
#
# Inputs:
#   Tables/SRK_individual_allele_genotypes.tsv  (Step 11)
#   SRK_individual_BL_assignments.tsv           (Step 13)
#
# Outputs:
#   figures/SRK_allele_eulerr_BLs.png
#   figures/SRK_allele_eulerr_BLs.pdf
#
# Note: with 5 sets, eulerr fits ellipses to the actual region counts via
# nonlinear optimisation; small "stress" / "diagError" values are reported.
# If the geometry cannot fit perfectly, the package issues a warning but
# still produces the best-fit diagram. The total allele count printed in
# the subtitle is authoritative.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(eulerr)
})

source("srk_bl_constants.R")

GENOTYPES <- "Tables/SRK_individual_allele_genotypes.tsv"
BL_FILE   <- "SRK_individual_BL_assignments.tsv"

# =============================================================================
# 1. LOAD + BUILD BL-LEVEL ALLELE SETS
# =============================================================================
cat("Loading allele matrix and BL assignments ...\n")

geno <- read_tsv(GENOTYPES, show_col_types = FALSE)
bl   <- read_tsv(BL_FILE,   show_col_types = FALSE) %>%
  filter(!is.na(BL), BL != "", BL != "Unassigned") %>%
  select(Individual, BL)

allele_cols <- setdiff(colnames(geno), "Individual")

# Long format: one row per (individual, allele present)
long <- geno %>%
  pivot_longer(all_of(allele_cols), names_to = "allele", values_to = "count") %>%
  filter(count > 0) %>%
  inner_join(bl, by = "Individual")

# Per-BL allele set
bl_sets <- split(long$allele, long$BL) |> lapply(unique)
bl_sets <- bl_sets[intersect(names(BL_COLORS), names(bl_sets))]

cat("\nAllele counts per BL:\n")
for (b in names(bl_sets)) {
  cat(sprintf("  %s: %d alleles (N individuals = %d)\n",
              b, length(bl_sets[[b]]),
              sum(bl$BL == b)))
}
total_alleles <- length(unique(unlist(bl_sets)))
cat(sprintf("\nTotal distinct alleles across BLs: %d\n", total_alleles))

# =============================================================================
# 2. FIT EULER + RENDER
# =============================================================================
cat("\nFitting Euler diagram ...\n")
fit <- euler(bl_sets, shape = "ellipse")

cat(sprintf("  stress = %.4f, diagError = %.4f\n",
            fit$stress, max(fit$diagError)))

# Match palette to set names; eulerr applies them in the order of the
# input list, so reorder to match BL1..BL5 for stable colour mapping.
fill_cols <- BL_COLORS[names(bl_sets)]

subtitle <- sprintf(
  "Per-BL allele counts: %s. Total distinct alleles across all BLs: %d.",
  paste(sprintf("%s = %d", names(bl_sets), lengths(bl_sets)), collapse = "; "),
  total_alleles
)

p <- plot(
  fit,
  fills = list(fill = fill_cols, alpha = 0.55),
  edges = list(col = "grey15", lwd = 1.4),
  labels = list(font = 2, cex = 1.1),
  quantities = list(type = "counts", cex = 1.0, font = 2),
  main = list(
    label = "S-allele sets per Bottleneck Lineage (area-proportional Euler)",
    cex = 1.15,
    font = 2
  )
)

# eulerr's plot does not natively render a subtitle, so add via grid.
grid::grid.text(
  subtitle,
  x = 0.5, y = 0.05,
  gp = grid::gpar(fontsize = 9, fontface = "italic", col = "grey25")
)

# =============================================================================
# 3. SAVE
# =============================================================================
dir.create("figures", showWarnings = FALSE)

# PNG
png("figures/SRK_allele_eulerr_BLs.png",
    width = 10, height = 8, units = "in", res = 200, bg = "white")
print(p)
grid::grid.text(subtitle, x = 0.5, y = 0.04,
                gp = grid::gpar(fontsize = 9, fontface = "italic", col = "grey25"))
dev.off()
cat("  Written: figures/SRK_allele_eulerr_BLs.png\n")

# PDF
pdf("figures/SRK_allele_eulerr_BLs.pdf", width = 10, height = 8)
print(p)
grid::grid.text(subtitle, x = 0.5, y = 0.04,
                gp = grid::gpar(fontsize = 9, fontface = "italic", col = "grey25"))
dev.off()
cat("  Written: figures/SRK_allele_eulerr_BLs.pdf\n")

cat("\nStep 18 (eulerr companion) complete.\n")
