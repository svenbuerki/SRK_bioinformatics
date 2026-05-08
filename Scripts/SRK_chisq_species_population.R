###############################
# SRK allele frequency chi-square tests
# Step 16 of the Canu_amplicon pipeline
#
# Tests whether observed allele copy-count frequencies deviate from the
# equal-frequency NFDS expectation. Three levels are tested:
#   1. Species (all 272 ingroup individuals)
#   2. Element Occurrence (EO; only EOs with N >= 5 individuals)
#   3. Bottleneck Lineage (BL1-BL5; all five always tested)
#
# The species pool MUST include every ingroup individual we have
# evidence for (including the 10 unresolved germplasm sub-codes that
# are not BL-assigned). The species_optimum baseline read from
# SRK_species_richness_estimates.tsv reflects that inclusive estimate.
###############################

cat("\nStarting SRK allele frequency chi-square analysis (Step 16)\n")

# Locked BL color palette (RColorBrewer Set1; must match Steps 14/15/17/18/20
# and the sibling LEPA_EO_spatial_clustering project).
BL_PALETTE <- c(
  BL1 = "#E41A1C", BL2 = "#377EB8", BL3 = "#4DAF4A",
  BL4 = "#984EA3", BL5 = "#FF7F00"
)
BL_UNASSIGNED_COLOR <- "#999999"

###############################
# 1. Load species richness optimum
###############################

species_optimum <- 50L  # fallback default

est_file <- "SRK_species_richness_estimates.tsv"
if (file.exists(est_file)) {
  est_df <- read.table(est_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
  if (!is.na(est_df$MM_estimate[1])) {
    species_optimum <- as.integer(est_df$MM_estimate[1])
    cat("Species richness optimum (MM):", species_optimum, "alleles\n")
  } else if (!is.na(est_df$Chao1_estimate[1])) {
    species_optimum <- as.integer(est_df$Chao1_estimate[1])
    cat("Species richness optimum (Chao1):", species_optimum, "alleles\n")
  } else if (!is.na(est_df$Consensus_estimate[1])) {
    species_optimum <- as.integer(est_df$Consensus_estimate[1])
    cat("Species richness optimum (consensus):", species_optimum, "alleles\n")
  } else {
    cat("No richness estimate available - using default optimum:",
        species_optimum, "\n")
  }
} else {
  cat("Richness estimates file not found - using default optimum:",
      species_optimum, "\n")
  cat("Run SRK_allele_accumulation_analysis.R first.\n")
}

###############################
# 2. Load genotype data, metadata, BL assignments
###############################

geno_all <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)
meta_df   <- read.csv("sampling_metadata.csv", stringsAsFactors = FALSE)
bl_assign <- read.table(
  "SRK_individual_BL_assignments.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Species-level filter: all ingroup individuals (272)
meta_df  <- meta_df[meta_df$Ingroup == 1, ]
geno_all <- geno_all[geno_all$Individual %in% meta_df$SampleID, ]
cat("Species-level individuals (all ingroup):", nrow(geno_all), "\n")

# EO/BL-level filter: BL-assigned subset (262)
bl_ok <- bl_assign[bl_assign$BL_status %in% c("Assigned", "Inferred"), ]
geno  <- geno_all[geno_all$Individual %in% bl_ok$Individual, ]
geno$EO <- bl_ok$EO[match(geno$Individual, bl_ok$Individual)]
geno$BL <- bl_ok$BL[match(geno$Individual, bl_ok$Individual)]
cat("EO/BL-level individuals (BL-assigned):", nrow(geno), "\n")

allele_cols <- setdiff(colnames(geno_all), "Individual")
geno_all[, allele_cols] <- lapply(geno_all[, allele_cols], as.numeric)
geno    [, allele_cols] <- lapply(geno    [, allele_cols], as.numeric)

cat("Allele bins:", length(allele_cols), "\n")

###############################
# 3. Chi-square helper
###############################

run_chisq <- function(mat, allele_cols) {
  allele_counts <- colSums(mat[, allele_cols, drop = FALSE])
  allele_counts <- allele_counts[allele_counts > 0]

  n_individuals <- nrow(mat)
  n_alleles     <- length(allele_counts)

  if (n_alleles < 2) {
    return(list(
      stats = data.frame(
        N_individuals = n_individuals,
        N_alleles     = n_alleles,
        X2 = NA_real_, df = NA_real_, p_value = NA_real_
      ),
      counts = allele_counts
    ))
  }

  obs      <- as.numeric(allele_counts)
  exp_prob <- rep(1 / length(obs), length(obs))

  chisq <- suppressWarnings(chisq.test(x = obs, p = exp_prob))

  list(
    stats = data.frame(
      N_individuals = n_individuals,
      N_alleles     = n_alleles,
      X2      = as.numeric(chisq$statistic),
      df      = as.numeric(chisq$parameter),
      p_value = chisq$p.value
    ),
    counts = allele_counts
  )
}

###############################
# 4. Run tests at all three levels
###############################

results_list <- list()
plot_list    <- list()
plot_meta    <- list()   # carries level, label, color, BL for plotting

push_result <- function(key, level, label, bl_for_color, res) {
  st <- res$stats
  st$Level      <- level
  st$Population <- label
  st$BL         <- if (level == "EO") bl_for_color else
                   if (level == "BL") label else NA_character_
  results_list[[key]] <<- st
  plot_list   [[key]] <<- res$counts
  plot_meta   [[key]] <<- list(
    level = level, label = label,
    color = if (!is.na(bl_for_color)) BL_PALETTE[bl_for_color] else "grey50"
  )
}

# ---- Species ----
res_sp <- run_chisq(geno_all, allele_cols)
push_result("Species_All", "Species", "All", NA_character_, res_sp)

# ---- EO loop (sorted by BL, then EO; only EOs with N >= 5) ----
eo_to_bl  <- unique(geno[, c("EO", "BL")])
eo_to_bl  <- eo_to_bl[order(eo_to_bl$BL, eo_to_bl$EO), ]
eo_sizes  <- table(geno$EO)
valid_eos <- intersect(eo_to_bl$EO, names(eo_sizes[eo_sizes >= 5]))
cat("EOs analyzed (N >= 5, sorted by BL):",
    paste(valid_eos, collapse = ", "), "\n")

for (eo in valid_eos) {
  bl_for_eo <- eo_to_bl$BL[eo_to_bl$EO == eo]
  df_eo     <- geno[geno$EO == eo, ]
  res       <- run_chisq(df_eo, allele_cols)
  push_result(paste0("EO_", eo), "EO", eo, bl_for_eo, res)
}

# ---- BL loop (all 5 BLs) ----
valid_bls <- sort(unique(geno$BL))
cat("BLs analyzed:", paste(valid_bls, collapse = ", "), "\n")

for (bl in valid_bls) {
  df_bl <- geno[geno$BL == bl, ]
  res   <- run_chisq(df_bl, allele_cols)
  push_result(paste0("BL_", bl), "BL", bl, bl, res)
}

###############################
# 5. Save TSV (Species first, then EO sorted by BL, then BL)
###############################

final_results <- do.call(rbind, results_list)
# Order: Species first, then EO rows (in valid_eos order), then BL rows
ord <- c(
  which(final_results$Level == "Species"),
  match(paste0("EO_", valid_eos), names(results_list)),
  match(paste0("BL_", valid_bls), names(results_list))
)
final_results <- final_results[ord, ]
rownames(final_results) <- NULL

write.table(
  final_results, "SRK_chisq_species_population.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Wrote SRK_chisq_species_population.tsv (",
    nrow(final_results), " rows)\n", sep = "")

print(final_results, row.names = FALSE)

###############################
# 6. Frequency plots
###############################

pdf("SRK_chisq_species_population_frequency_plots.pdf", width = 10, height = 6)

plot_one <- function(meta_entry, allele_counts, n_alleles_plot, p_value) {
  level <- meta_entry$level
  label <- meta_entry$label
  bl_color <- meta_entry$color

  if (length(allele_counts) < 2) {
    plot.new()
    title(main = paste(level, label, sep = " - "),
          sub  = "fewer than 2 alleles - chi-square not computable")
    return(invisible())
  }

  obs_freq <- sort(allele_counts / sum(allele_counts), decreasing = TRUE)
  exp_freq <- rep(1 / length(obs_freq), length(obs_freq))

  p_stars <- if (is.na(p_value)) "" else
             if (p_value < 0.001) "***" else
             if (p_value < 0.01)  "**"  else
             if (p_value < 0.05)  "*"   else "ns"
  p_label <- if (is.na(p_value)) "p = NA" else
             if (p_value < 0.001) paste0("p = ",
                  formatC(p_value, format = "e", digits = 2)) else
             paste0("p = ", round(p_value, 3))
  plot_sub <- paste0(n_alleles_plot, " alleles observed | ",
                     p_label, "  ", p_stars)

  bar_fill <- adjustcolor(bl_color, alpha.f = 0.65)

  bp <- barplot(
    obs_freq, space = 0, border = "white", col = bar_fill,
    names.arg = rep("", length(obs_freq)),
    xlim = c(0, species_optimum),
    ylim = c(0, max(obs_freq) * 1.2),
    main = paste(level, label, sep = " - "),
    sub  = plot_sub,
    ylab = "Allele frequency", xlab = "Allele (ranked)"
  )

  # Missing-to-optimum zone
  lgd <- legend(
    "topright",
    legend = c("Observed", "NFDS expectation",
               "Drift-smoothed", "Missing to optimum"),
    col    = c(bl_color, "blue", "red", "tomato"),
    lwd = c(2, 2, 2, NA), lty = c(1, 2, 1, NA),
    pch = c(NA, NA, NA, 22),
    pt.bg = c(NA, NA, NA, adjustcolor("tomato", alpha.f = 0.3)),
    pt.cex = 1.5, bty = "n", plot = FALSE
  )

  n_obs     <- length(obs_freq)
  n_missing <- max(0L, species_optimum - n_obs)
  if (n_missing > 0) {
    last_x   <- max(bp) + 0.5
    zone_mid <- (last_x + species_optimum) / 2
    rect_top <- lgd$rect$top - lgd$rect$h
    rect(last_x, 0, species_optimum, rect_top,
         col = adjustcolor("tomato", alpha.f = 0.15),
         border = adjustcolor("tomato", alpha.f = 0.5), lty = 2)
    text(x = zone_mid, y = rect_top / 2,
         labels = paste0(n_missing, " alleles\nto optimum"),
         col = "tomato4", cex = 0.85, font = 2)
  }

  abline(h = exp_freq[1], col = "blue", lwd = 2, lty = 2)
  lines(seq_along(obs_freq), obs_freq, lwd = 2, col = bl_color)

  lo <- loess(obs_freq ~ seq_along(obs_freq))
  lines(seq_along(obs_freq), predict(lo), col = "red", lwd = 2)

  legend(
    "topright",
    legend = c("Observed", "NFDS expectation",
               "Drift-smoothed", "Missing to optimum"),
    col = c(bl_color, "blue", "red", "tomato"),
    lwd = c(2, 2, 2, NA), lty = c(1, 2, 1, NA),
    pch = c(NA, NA, NA, 22),
    pt.bg = c(NA, NA, NA, adjustcolor("tomato", alpha.f = 0.3)),
    pt.cex = 1.5, bty = "n"
  )
}

# Plot in the same order as the TSV: Species, EOs (BL-sorted), BLs
plot_keys <- c(
  "Species_All",
  paste0("EO_", valid_eos),
  paste0("BL_", valid_bls)
)
for (key in plot_keys) {
  if (is.null(plot_meta[[key]])) next
  i <- which(rownames(final_results) == key |
             paste0(final_results$Level, "_", final_results$Population) == sub("^Species_", "Species_", key))
  # Resolve N_alleles and p_value from final_results matching this row
  level_lab <- plot_meta[[key]]$level
  label_lab <- plot_meta[[key]]$label
  row_idx <- which(final_results$Level == level_lab &
                   final_results$Population == label_lab)[1]
  plot_one(plot_meta[[key]],
           plot_list[[key]],
           final_results$N_alleles[row_idx],
           final_results$p_value [row_idx])
}

dev.off()
cat("Wrote SRK_chisq_species_population_frequency_plots.pdf (",
    length(plot_keys), " pages)\n", sep = "")

cat("\nStep 16 complete.\n")
