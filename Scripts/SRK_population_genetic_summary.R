############################################################
# SRK population genetic summary
# Step 14 of the Canu_amplicon pipeline
#
# Produces population genetic statistics at TWO levels:
#   1. Element Occurrence (EO) — the management unit
#   2. Bottleneck Lineage (BL) — the evolutionary unit (5 independent
#      bottleneck lineages BL1-BL5 from the LEPA spatial-clustering
#      project; see SRK_BL_integration.py / Step 13)
#
# Outputs
#   SRK_population_genetic_summary.tsv      EO-level (with BL column),
#                                           sorted by BL then EO
#   SRK_population_genetic_summary_BL.tsv   BL-level aggregates
#   SRK_population_genetic_summary.pdf      2-page PDF: page 1 EO bars
#                                           sorted by BL, page 2 BL bars
#                                           Bars colored by parent BL
#                                           (RColorBrewer Dark2 palette).
############################################################

cat("Starting SRK population genetic summary analysis (Step 14)\n")

############################
# Shared BL color palette (used here and in Steps 15, 17, 18, 20)
############################
# RColorBrewer Set1, BL1 -> BL5. Matches the locked palette used in the
# sibling LEPA_EO_spatial_clustering project so BL identity is visually
# consistent across all figures in both repositories.
BL_PALETTE <- c(
  BL1 = "#E41A1C",   # red
  BL2 = "#377EB8",   # blue
  BL3 = "#4DAF4A",   # green
  BL4 = "#984EA3",   # purple
  BL5 = "#FF7F00"    # orange
)
BL_UNASSIGNED_COLOR <- "#999999"

############################
# 1. Load files
############################

geno <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)

zyg <- read.table(
  "SRK_individual_zygosity.tsv",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE
)

bl_assign <- read.table(
  "SRK_individual_BL_assignments.tsv",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE
)

cat("Loaded BL assignments for", nrow(bl_assign), "individuals\n")

############################
# 2. Restrict to BL-assigned individuals
############################
# Drops the 10 individuals whose Pop codes (germplasm sub-codes) could
# not be resolved to an EO in the spatial-clustering BL summary.

bl_ok <- bl_assign[bl_assign$BL_status %in% c("Assigned", "Inferred"), ]
cat("BL-assigned individuals (Assigned + Inferred):", nrow(bl_ok), "\n")

geno <- geno[geno$Individual %in% bl_ok$Individual, ]
zyg  <- zyg [zyg$Individual  %in% bl_ok$Individual, ]

############################
# 3. Attach EO and BL labels
############################

geno$EO <- bl_ok$EO[match(geno$Individual, bl_ok$Individual)]
geno$BL <- bl_ok$BL[match(geno$Individual, bl_ok$Individual)]
zyg$EO  <- bl_ok$EO[match(zyg$Individual,  bl_ok$Individual)]
zyg$BL  <- bl_ok$BL[match(zyg$Individual,  bl_ok$Individual)]

# Drift index per EO (mean across that EO's groups; one value per EO)
eo_drift <- tapply(
  as.numeric(bl_ok$Drift_index),
  bl_ok$EO,
  function(x) round(mean(x, na.rm = TRUE), 3)
)

############################
# 4. Standardised genotype processing
############################
# Same logic as the accumulation-curve script — keeps allele-count
# dosage information, robust to character/factor input columns.

standardized_genotype_processing <- function(data, exclude_columns) {
  protein_cols <- setdiff(colnames(data), exclude_columns)
  cat("Total protein columns detected:", length(protein_cols), "\n")

  protein_data <- data[, protein_cols, drop = FALSE]

  for (i in seq_along(protein_cols)) {
    col_data <- protein_data[, i]
    if (is.character(col_data) || is.factor(col_data)) {
      col_data <- as.character(col_data)
      col_data[is.na(col_data) | col_data == "" | col_data == "NA"] <- "0"

      numeric_col <- suppressWarnings(as.numeric(col_data))

      if (any(is.na(numeric_col) & col_data != "0")) {
        numeric_col <- ifelse(col_data == "0" | is.na(col_data), 0, 1)
      }

      protein_data[, i] <- numeric_col
    } else {
      protein_data[, i] <- ifelse(is.na(protein_data[, i]), 0, protein_data[, i])
    }
  }

  mat <- as.matrix(protein_data)
  mat <- apply(mat, 2, as.numeric)
  mat[is.na(mat)] <- 0
  list(matrix = mat, protein_cols = protein_cols)
}

processed   <- standardized_genotype_processing(
  geno, exclude_columns = c("Individual", "EO", "BL")
)
geno_matrix <- processed$matrix
allele_cols <- processed$protein_cols
geno[, allele_cols] <- geno_matrix
cat("Genotype matrix processing complete\n")

############################
# 5. Generic per-group analysis
############################
# Computes the same set of metrics regardless of whether `group_label`
# is an EO ("EO27") or a BL ("BL2"). `g` and `z` must already be
# subsetted to the group's individuals.

analyze_group <- function(group_label, g, z) {
  mat <- as.matrix(g[, allele_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  if (nrow(mat) == 1) {
    allele_counts <- as.numeric(mat[1, ])
    names(allele_counts) <- colnames(mat)
  } else {
    allele_counts <- colSums(mat, na.rm = TRUE)
  }

  present   <- allele_counts > 0
  N_alleles <- sum(present)

  cat("  ", group_label, ": ", N_alleles, " alleles in ",
      nrow(mat), " individuals\n", sep = "")

  if (N_alleles > 0) {
    freqs <- allele_counts[present] / sum(allele_counts[present])
    Ne    <- 1 / sum(freqs^2)
    top   <- head(sort(allele_counts[present], decreasing = TRUE), 3)
    top_str <- paste0(names(top), "(", top, ")", collapse = ", ")
  } else {
    freqs   <- numeric(0)
    Ne      <- NA_real_
    top_str <- "None"
  }

  if (nrow(z) > 0) {
    prop_het  <- mean(z$Zygosity == "Heterozygous", na.rm = TRUE)
    mean_prot <- mean(z$N_distinct_alleles,         na.rm = TRUE)
  } else {
    prop_het  <- NA_real_
    mean_prot <- NA_real_
  }

  if (nrow(mat) == 1) {
    individuals_with_alleles <- if (sum(allele_counts) > 0) 1L else 0L
  } else {
    individuals_with_alleles <- sum(rowSums(mat) > 0)
  }

  data.frame(
    label                       = group_label,
    N_individuals               = nrow(g),
    N_individuals_with_alleles  = individuals_with_alleles,
    N_alleles                   = N_alleles,
    Effective_alleles_Ne        = round(Ne, 3),
    Prop_heterozygous           = round(prop_het, 3),
    Prop_homozygous             = round(1 - prop_het, 3),
    Mean_alleles                = round(mean_prot, 3),
    Top_alleles                 = top_str,
    stringsAsFactors            = FALSE
  )
}

############################
# 6. EO-level analysis (sorted by BL then EO)
############################

eo_to_bl <- unique(bl_ok[, c("EO", "BL")])
eo_to_bl <- eo_to_bl[order(eo_to_bl$BL, eo_to_bl$EO), ]
eos_sorted <- eo_to_bl$EO

cat("\nEO-level analysis (", length(eos_sorted), " EOs, sorted by BL):\n", sep = "")

eo_results <- do.call(rbind, lapply(eos_sorted, function(eo) {
  res <- analyze_group(eo,
                       geno[geno$EO == eo, ],
                       zyg [zyg$EO  == eo, ])
  res$EO           <- eo
  res$BL           <- eo_to_bl$BL[eo_to_bl$EO == eo]
  res$Drift_index  <- eo_drift[eo]
  res
}))

eo_results$label <- NULL
eo_results <- eo_results[, c(
  "BL", "EO", "Drift_index",
  "N_individuals", "N_individuals_with_alleles",
  "N_alleles", "Effective_alleles_Ne",
  "Prop_heterozygous", "Prop_homozygous", "Mean_alleles", "Top_alleles"
)]
# Backwards-compat: Population mirrors EO so legacy downstream scripts keep working
eo_results$Population <- eo_results$EO

write.table(
  eo_results,
  "SRK_population_genetic_summary.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("\nWritten EO summary -> SRK_population_genetic_summary.tsv (",
    nrow(eo_results), " rows)\n", sep = "")

############################
# 7. BL-level analysis
############################

bls_sorted <- sort(unique(bl_ok$BL))

cat("\nBL-level analysis (", length(bls_sorted), " BLs):\n", sep = "")

bl_results <- do.call(rbind, lapply(bls_sorted, function(bl) {
  res <- analyze_group(bl,
                       geno[geno$BL == bl, ],
                       zyg [zyg$BL  == bl, ])
  res$BL           <- bl
  res$N_EOs        <- length(unique(bl_ok$EO[bl_ok$BL == bl]))
  res$Mean_drift   <- round(mean(as.numeric(
    bl_ok$Drift_index[bl_ok$BL == bl]), na.rm = TRUE), 3)
  res
}))

bl_results$label <- NULL
bl_results <- bl_results[, c(
  "BL", "N_EOs", "Mean_drift",
  "N_individuals", "N_individuals_with_alleles",
  "N_alleles", "Effective_alleles_Ne",
  "Prop_heterozygous", "Prop_homozygous", "Mean_alleles", "Top_alleles"
)]

write.table(
  bl_results,
  "SRK_population_genetic_summary_BL.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Written BL summary -> SRK_population_genetic_summary_BL.tsv (",
    nrow(bl_results), " rows)\n", sep = "")

############################
# 8. Console summaries
############################

cat("\n=== EO-LEVEL SUMMARY (sorted by BL) ===\n")
print(eo_results[, c("BL", "EO", "N_individuals", "N_alleles",
                     "Effective_alleles_Ne", "Prop_heterozygous")],
      row.names = FALSE)

cat("\n=== BL-LEVEL SUMMARY ===\n")
print(bl_results[, c("BL", "N_EOs", "N_individuals", "N_alleles",
                     "Effective_alleles_Ne", "Prop_heterozygous")],
      row.names = FALSE)

############################
# 9. Plots — page 1: EO sorted by BL; page 2: BL aggregate
############################

eo_colors <- BL_PALETTE[eo_results$BL]
eo_colors[is.na(eo_colors)] <- BL_UNASSIGNED_COLOR

bl_colors <- BL_PALETTE[bl_results$BL]
bl_colors[is.na(bl_colors)] <- BL_UNASSIGNED_COLOR

bar_panel <- function(values, labels, fill, ylab, main) {
  bp <- barplot(values, names.arg = labels, col = fill,
                ylab = ylab, main = main, las = 2,
                ylim = c(0, max(values, na.rm = TRUE) * 1.15))
  text(x = bp, y = values + max(values, na.rm = TRUE) * 0.03,
       labels = ifelse(is.na(values), "", values),
       cex = 0.7)
  invisible(bp)
}

pdf("SRK_population_genetic_summary.pdf", width = 14, height = 9)

# ---- Page 1: EO bars sorted by BL, colored by BL ----
par(mfrow = c(2, 2), mar = c(6, 4, 3, 1), oma = c(0, 0, 3, 0))
bar_panel(eo_results$Effective_alleles_Ne, eo_results$EO, eo_colors,
          "Effective number of alleles (Ne)", "Effective Allele Diversity")
legend("topleft", legend = names(BL_PALETTE), fill = unname(BL_PALETTE),
       bty = "n", cex = 0.85, ncol = 1, title = "BL", title.font = 2)
bar_panel(eo_results$Prop_heterozygous,    eo_results$EO, eo_colors,
          "Proportion heterozygous",         "Heterozygosity")
bar_panel(eo_results$N_alleles,            eo_results$EO, eo_colors,
          "Number of alleles",               "Allele Richness")
bar_panel(eo_results$N_individuals,        eo_results$EO, eo_colors,
          "Number of individuals",           "Sample Sizes")
mtext("EO-level statistics  bars sorted by BL, colored by parent BL",
      side = 3, outer = TRUE, line = 0.8, cex = 1.05, font = 2)

# ---- Page 2: BL aggregate bars ----
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))
bar_panel(bl_results$Effective_alleles_Ne, bl_results$BL, bl_colors,
          "Effective number of alleles (Ne)", "Effective Allele Diversity")
bar_panel(bl_results$Prop_heterozygous,    bl_results$BL, bl_colors,
          "Proportion heterozygous",         "Heterozygosity")
bar_panel(bl_results$N_alleles,            bl_results$BL, bl_colors,
          "Number of alleles",               "Allele Richness")
bar_panel(bl_results$N_individuals,        bl_results$BL, bl_colors,
          "Number of individuals",           "Sample Sizes")
mtext("BL-level aggregates  five independent bottleneck lineages",
      side = 3, outer = TRUE, line = 0.8, cex = 1.05, font = 2)

dev.off()
cat("\nPlots saved -> SRK_population_genetic_summary.pdf\n")
cat("\nStep 14 complete.\n")
