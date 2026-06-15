############################################################
# SRK population genetic summary — NULL-AWARE variant
# Step 14b of the Canu_amplicon pipeline (parallel to Step 14)
#
# Reads the Step 26 augmented tables that carry explicit Allele_NULL
# copies for pSI_high and SC individuals, and recomputes the standard
# population-genetic metrics so the bias introduced by treating pSI
# individuals as fully SI is corrected.
#
# Key differences from the functional-only Step 14:
#   * Reads `SRK_individual_allele_genotypes_with_nulls.tsv`
#     and `SRK_individual_zygosity_with_nulls.tsv`.
#   * EO / BL labels come from the augmented zygosity table
#     (no `SRK_individual_BL_assignments.tsv` join needed — the seven
#     SC individuals added by Step 26 already carry EO/BL).
#   * `N_alleles` reports the number of distinct FUNCTIONAL S-alleles
#     (Allele_NULL is excluded — it is not an S-allele, it is a broken copy).
#   * `Effective_alleles_Ne` is computed on the full count vector INCLUDING
#     Allele_NULL, so that the broken-copy frequency contributes to the
#     allele-frequency denominator (per Sven, 2026-06-11).
#   * `Prop_heterozygous` is recomputed on N_distinct_functional_alleles >= 2.
#     SC individuals (genotype 0000) are reported separately as `N_SC`.
#   * New columns: `N_SC` (count of all-null individuals), `Mean_null_copies`
#     (average Allele_NULL count, 0-4), `Prop_pSI` (fraction with at least
#     one null copy but not all four), `Frac_nonfunctional_alleles` (sum of
#     null copies / total copies — population-level broken-copy frequency).
#
# Outputs
#   SRK_population_genetic_summary_with_nulls.tsv      EO-level (with BL)
#   SRK_population_genetic_summary_BL_with_nulls.tsv   BL-level aggregates
#   SRK_population_genetic_summary_with_nulls.pdf      2-page PDF
############################################################

cat("Starting SRK population genetic summary — NULL-AWARE (Step 14b)\n")

source("srk_bl_constants.R")
BL_PALETTE <- BL_COLORS
BL_UNASSIGNED_COLOR <- "#999999"

############################
# 1. Load files
############################
geno <- read.table("Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv",
                   header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE)
zyg  <- read.table("Tables/Phase4/step23_individual_zygosity_with_nulls.tsv",
                   header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)

cat("Loaded", nrow(geno), "individuals with null-aware genotypes\n")

############################
# 2. Restrict to BL-assigned individuals
############################
# Drop Unassigned BL (consistent with Step 14 which restricts to BL1-BL5).
geno <- geno[geno$BL_inferred %in% BL_ORDER, ]
zyg  <- zyg [zyg$BL_inferred  %in% BL_ORDER, ]
cat("After BL filter:", nrow(geno), "individuals retained\n")

# Drift_index — read from the EO_group_BL_summary the same way
# get_eo_order_within_bl() does internally.
eo_di_src <- read.csv("Tables/EO_group_BL_summary.csv", stringsAsFactors = FALSE)
parts <- strsplit(eo_di_src$EO, "[;,]\\s*")
eo_di_src <- eo_di_src[rep(seq_len(nrow(eo_di_src)), lengths(parts)), ]
eo_di_src$EO <- trimws(unlist(parts))
eo_drift <- tapply(eo_di_src$Drift_index, eo_di_src$EO,
                   function(x) round(mean(x, na.rm = TRUE), 3))

############################
# 3. Allele columns
############################
exclude <- c("Individual", "Genotype_class_flag", "EO_normalised", "BL_inferred",
             "SI_status", "pSI_confidence")
allele_cols      <- setdiff(colnames(geno), exclude)         # 49 functional + Allele_NULL
allele_cols_func <- setdiff(allele_cols, "Allele_NULL")      # 49 functional only

# Ensure numeric
for (a in allele_cols) geno[[a]] <- as.numeric(geno[[a]])
geno[is.na(geno)] <- 0

############################
# 4. Generic per-group analysis
############################
analyze_group <- function(label, g, z) {
  if (nrow(g) == 0) {
    return(data.frame(label = label, N_individuals = 0,
                      N_alleles_functional = 0, Effective_alleles_Ne_with_NULL = NA,
                      Prop_heterozygous = NA, Mean_null_copies = NA,
                      Frac_nonfunctional_alleles = NA, N_SC = 0,
                      N_pSI = 0, N_SI = 0, Prop_pSI = NA,
                      Top_alleles = "None", stringsAsFactors = FALSE))
  }

  mat_full <- as.matrix(g[, allele_cols,      drop = FALSE]); storage.mode(mat_full) <- "numeric"
  mat_func <- as.matrix(g[, allele_cols_func, drop = FALSE]); storage.mode(mat_func) <- "numeric"

  counts_full <- if (nrow(mat_full) == 1) as.numeric(mat_full[1,]) else colSums(mat_full)
  names(counts_full) <- allele_cols
  counts_func <- if (nrow(mat_func) == 1) as.numeric(mat_func[1,]) else colSums(mat_func)
  names(counts_func) <- allele_cols_func

  N_alleles_func <- sum(counts_func > 0)

  # Ne computed on full count vector including Allele_NULL — broken copies
  # are a real population-genetic signal at the locus.
  freqs_full <- counts_full / sum(counts_full)
  Ne <- if (sum(counts_full) > 0) 1 / sum(freqs_full^2) else NA_real_

  # Heterozygosity based on number of distinct functional S-alleles >= 2.
  # SC individuals (0 distinct functional) are NOT counted in either bucket
  # for proportion heterozygous; their effect is captured in N_SC + Mean_null.
  n_func_distinct <- z$N_distinct_functional_alleles
  het_eligible    <- n_func_distinct >= 1
  prop_het        <- if (sum(het_eligible) > 0) mean(n_func_distinct[het_eligible] >= 2) else NA_real_

  mean_null  <- mean(z$N_null_copies, na.rm = TRUE)
  total_cp   <- sum(z$N_total_copies, na.rm = TRUE)
  frac_nf    <- if (total_cp > 0) sum(z$N_null_copies, na.rm = TRUE) / total_cp else NA_real_

  n_sc  <- sum(z$SI_status == "SC", na.rm = TRUE)
  n_psi <- sum(z$SI_status == "pSI" & z$pSI_confidence == "high", na.rm = TRUE)
  n_si  <- sum(z$SI_status == "SI" |
               (z$SI_status == "pSI" & z$pSI_confidence == "low"), na.rm = TRUE)
  prop_psi <- if (nrow(z) > 0) n_psi / nrow(z) else NA_real_

  if (N_alleles_func > 0) {
    top   <- head(sort(counts_func[counts_func > 0], decreasing = TRUE), 3)
    top_str <- paste0(names(top), "(", top, ")", collapse = ", ")
  } else {
    top_str <- "None"
  }

  data.frame(
    label                          = label,
    N_individuals                  = nrow(g),
    N_alleles_functional           = N_alleles_func,
    Effective_alleles_Ne_with_NULL = round(Ne, 3),
    Prop_heterozygous              = round(prop_het, 3),
    Mean_null_copies               = round(mean_null, 3),
    Frac_nonfunctional_alleles     = round(frac_nf, 3),
    N_SI                           = n_si,
    N_pSI                          = n_psi,
    N_SC                           = n_sc,
    Prop_pSI                       = round(prop_psi, 3),
    Top_alleles                    = top_str,
    stringsAsFactors               = FALSE
  )
}

############################
# 5. EO-level analysis
############################
eo_to_bl <- unique(geno[, c("EO_normalised", "BL_inferred")])
eos_sorted <- get_eo_order_within_bl(eo_to_bl$EO_normalised)
eos_sorted <- intersect(eos_sorted, eo_to_bl$EO_normalised)
eo_to_bl <- eo_to_bl[match(eos_sorted, eo_to_bl$EO_normalised), ]

cat("\nEO-level (n=", length(eos_sorted), " EOs):\n", sep = "")
eo_results <- do.call(rbind, lapply(eos_sorted, function(eo) {
  res <- analyze_group(eo,
                       geno[geno$EO_normalised == eo, ],
                       zyg [zyg$EO_normalised  == eo, ])
  res$EO          <- eo
  res$BL          <- eo_to_bl$BL_inferred[eo_to_bl$EO_normalised == eo]
  res$Drift_index <- eo_drift[eo]
  res
}))
eo_results$label <- NULL
eo_results <- eo_results[, c(
  "BL", "EO", "Drift_index",
  "N_individuals", "N_SI", "N_pSI", "N_SC",
  "N_alleles_functional", "Effective_alleles_Ne_with_NULL",
  "Prop_heterozygous", "Mean_null_copies", "Frac_nonfunctional_alleles",
  "Prop_pSI", "Top_alleles"
)]

write.table(eo_results, "Tables/Phase4/step14b_population_genetic_summary_with_nulls.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Written EO summary -> SRK_population_genetic_summary_with_nulls.tsv (",
    nrow(eo_results), " rows)\n", sep = "")

############################
# 6. BL-level analysis
############################
bls_sorted <- intersect(BL_ORDER, unique(geno$BL_inferred))

cat("\nBL-level (n=", length(bls_sorted), " BLs):\n", sep = "")
bl_results <- do.call(rbind, lapply(bls_sorted, function(bl) {
  res <- analyze_group(bl,
                       geno[geno$BL_inferred == bl, ],
                       zyg [zyg$BL_inferred  == bl, ])
  res$BL    <- bl
  res$N_EOs <- length(unique(geno$EO_normalised[geno$BL_inferred == bl]))
  res
}))
bl_results$label <- NULL
bl_results <- bl_results[, c(
  "BL", "N_EOs", "N_individuals", "N_SI", "N_pSI", "N_SC",
  "N_alleles_functional", "Effective_alleles_Ne_with_NULL",
  "Prop_heterozygous", "Mean_null_copies", "Frac_nonfunctional_alleles",
  "Prop_pSI", "Top_alleles"
)]

write.table(bl_results, "Tables/Phase4/step14b_population_genetic_summary_BL_with_nulls.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Written BL summary -> SRK_population_genetic_summary_BL_with_nulls.tsv (",
    nrow(bl_results), " rows)\n", sep = "")

############################
# 7. Console summary
############################
cat("\n=== EO-LEVEL SUMMARY (null-aware) ===\n")
print(eo_results[, c("BL", "EO", "N_individuals", "N_alleles_functional",
                     "Effective_alleles_Ne_with_NULL", "Prop_heterozygous",
                     "Frac_nonfunctional_alleles", "Prop_pSI", "N_SC")],
      row.names = FALSE)

cat("\n=== BL-LEVEL SUMMARY (null-aware) ===\n")
print(bl_results[, c("BL", "N_EOs", "N_individuals", "N_alleles_functional",
                     "Effective_alleles_Ne_with_NULL", "Prop_heterozygous",
                     "Frac_nonfunctional_alleles", "Prop_pSI", "N_SC")],
      row.names = FALSE)

############################
# 8. Plots
############################
eo_colors <- BL_PALETTE[eo_results$BL]; eo_colors[is.na(eo_colors)] <- BL_UNASSIGNED_COLOR
bl_colors <- BL_PALETTE[bl_results$BL]; bl_colors[is.na(bl_colors)] <- BL_UNASSIGNED_COLOR

bar_panel <- function(values, labels, fill, ylab, main) {
  bp <- barplot(values, names.arg = labels, col = fill,
                ylab = ylab, main = main, las = 2,
                ylim = c(0, max(values, na.rm = TRUE) * 1.18))
  text(x = bp, y = values + max(values, na.rm = TRUE) * 0.04,
       labels = ifelse(is.na(values), "", round(values, 2)),
       cex = 0.7)
  invisible(bp)
}

pdf("figures/Phase4/step14b_population_genetic_summary_with_nulls.pdf", width = 14, height = 9)
par(mfrow = c(2, 2), mar = c(6, 4, 3, 1), oma = c(0, 0, 3, 0))
bar_panel(eo_results$Effective_alleles_Ne_with_NULL, eo_results$EO, eo_colors,
          "Effective allele diversity (incl. NULL)", "Ne (with NULL)")
legend("topleft", legend = BL_ORDER_NUMERIC,
       fill   = unname(BL_PALETTE[BL_ORDER_NUMERIC]),
       bty = "n", cex = 0.85, ncol = 1, title = "BL", title.font = 2)
bar_panel(eo_results$Frac_nonfunctional_alleles, eo_results$EO, eo_colors,
          "Frac non-functional copies", "Null-allele Frequency")
bar_panel(eo_results$N_alleles_functional, eo_results$EO, eo_colors,
          "Number of functional S-alleles", "Functional Allele Richness")
bar_panel(eo_results$Prop_pSI, eo_results$EO, eo_colors,
          "Proportion pSI individuals", "Partial-SI Prevalence")
mtext("EO-level null-aware statistics — bars sorted by BL, coloured by parent BL",
      side = 3, outer = TRUE, line = 0.8, cex = 1.05, font = 2)

par(mfrow = c(2, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 3, 0))
bar_panel(bl_results$Effective_alleles_Ne_with_NULL, bl_results$BL, bl_colors,
          "Effective allele diversity (incl. NULL)", "Ne (with NULL)")
bar_panel(bl_results$Frac_nonfunctional_alleles, bl_results$BL, bl_colors,
          "Frac non-functional copies", "Null-allele Frequency")
bar_panel(bl_results$N_alleles_functional, bl_results$BL, bl_colors,
          "Number of functional S-alleles", "Functional Allele Richness")
bar_panel(bl_results$Prop_pSI, bl_results$BL, bl_colors,
          "Proportion pSI individuals", "Partial-SI Prevalence")
mtext("BL-level null-aware aggregates — five independent bottleneck lineages",
      side = 3, outer = TRUE, line = 0.8, cex = 1.05, font = 2)
dev.off()
cat("\nPlots saved -> SRK_population_genetic_summary_with_nulls.pdf\n")
cat("\nStep 14b complete.\n")
