###############################
# SRK protein frequency χ² tests
# Compatible with functional protein pipeline
###############################

cat("\nStarting SRK chi-square analysis\n")

###############################
# Load species richness estimate
# (produced by SRK_allele_accumulation_analysis.R)
# Used as the "optimum" allele count for frequency plots.
###############################

species_optimum <- 50L  # fallback default

est_file <- "SRK_species_richness_estimates.tsv"
if (file.exists(est_file)) {
  est_df <- read.table(est_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Prefer MM: its asymptote directly reflects the observed accumulation curve.
  # Chao1 can be unreliable (very high SE) when many alleles are rare, which
  # is common under NFDS. Fall back to Chao1 only if MM failed, then consensus.
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
    cat("No richness estimate available — using default optimum:", species_optimum, "\n")
  }
} else {
  cat("Richness estimates file not found — using default optimum:", species_optimum, "\n")
  cat("Run SRK_allele_accumulation_analysis.R first to generate estimates.\n")
}

###############################
# Load data
###############################

geno <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_df <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

###############################
# Keep ingroup only
###############################

meta_df <- meta_df[meta_df$Ingroup == 1, ]

geno <- geno[geno$Individual %in% meta_df$SampleID, ]

geno$Population <- meta_df$Pop[match(geno$Individual, meta_df$SampleID)]

allele_cols <- setdiff(colnames(geno), c("Individual", "Population"))

geno[, allele_cols] <- lapply(geno[, allele_cols], as.numeric)

cat("Individuals analyzed:", nrow(geno), "\n")
cat("Allele bins:", length(allele_cols), "\n")

###############################
# Safe chi-square function
###############################

run_chisq <- function(mat, allele_cols) {

  # Sum allele copy counts across individuals; remove absent alleles
  allele_counts <- colSums(mat[, allele_cols, drop = FALSE])
  allele_counts <- allele_counts[allele_counts > 0]

  n_individuals <- nrow(mat)
  n_alleles     <- length(allele_counts)

  if (n_alleles < 2) {

    return(list(

      stats = data.frame(

        N_individuals = n_individuals,
        N_alleles     = n_alleles,
        X2 = NA,
        df = NA,
        p_value = NA
      ),

      counts = allele_counts
    ))
  }

  obs      <- as.numeric(allele_counts)
  exp_prob <- rep(1 / length(obs), length(obs))

  chisq <- suppressWarnings(
    chisq.test(x = obs, p = exp_prob)
  )

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
# Storage
###############################

results_list <- list()
plot_list <- list()

###############################
# Species-level
###############################

res_species <- run_chisq(geno, allele_cols)

stats_species <- res_species$stats
stats_species$Level <- "Species"
stats_species$Population <- "All"

results_list[["Species_All"]] <- stats_species
plot_list[["Species_All"]] <- res_species$counts

###############################
# Population-level
###############################

# One row per individual in geno — table() correctly counts unique individuals
pop_sizes <- table(geno$Population)

valid_pops <- names(pop_sizes[pop_sizes >= 5])

cat("Populations analyzed:", length(valid_pops), "\n")

for (pop in valid_pops) {

  df_pop <- geno[geno$Population == pop, ]

  res <- run_chisq(df_pop, allele_cols)

  stats <- res$stats

  stats$Level <- "Population"
  stats$Population <- pop

  key <- paste0("Population_", pop)

  results_list[[key]] <- stats
  plot_list[[key]] <- res$counts
}

###############################
# Combine results
###############################

final_results <- do.call(rbind, results_list)

final_results <- final_results[
  order(final_results$Level,
        final_results$Population),
]

###############################
# Save TSV
###############################

write.table(

  final_results,

  "SRK_chisq_species_population.tsv",

  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################
# Frequency plots
###############################

pdf(
  "SRK_chisq_species_population_frequency_plots.pdf",
  width = 10,
  height = 6
)

for (i in seq_len(nrow(final_results))) {

  level <- final_results$Level[i]
  pop   <- final_results$Population[i]

  key <- if (level == "Species") {

    "Species_All"

  } else {

    paste0("Population_", pop)

  }

  allele_counts <- plot_list[[key]]

  if (length(allele_counts) < 2) next

  obs_freq <- sort(
    allele_counts / sum(allele_counts),
    decreasing = TRUE
  )

  exp_freq <- rep(
    1 / length(obs_freq),
    length(obs_freq)
  )

  n_alleles_plot <- final_results$N_alleles[i]
  p_val_plot     <- final_results$p_value[i]

  p_stars <- if (is.na(p_val_plot)) {
    ""
  } else if (p_val_plot < 0.001) {
    "***"
  } else if (p_val_plot < 0.01) {
    "**"
  } else if (p_val_plot < 0.05) {
    "*"
  } else {
    "ns"
  }

  p_label <- if (is.na(p_val_plot)) {
    "p = NA"
  } else if (p_val_plot < 0.001) {
    paste0("p = ", formatC(p_val_plot, format = "e", digits = 2))
  } else {
    paste0("p = ", round(p_val_plot, 3))
  }

  plot_sub <- paste0(n_alleles_plot, " alleles observed  |  ", p_label, "  ", p_stars)

  bp <- barplot(

    obs_freq,

    space  = 0,
    border = "white",
    col = "gray70",
    names.arg = rep("", length(obs_freq)),

    xlim = c(0, species_optimum),

    ylim = c(
      0,
      max(obs_freq) * 1.2
    ),

    main = paste(level, pop, sep = " – "),
    sub  = plot_sub,
    ylab = "Allele frequency",
    xlab = "Allele (ranked)"
  )

  # --- Missing alleles zone ---
  # Dry-run the legend first to obtain its bounding box, so the rectangle
  # can be capped exactly at the legend's bottom edge with no overlap.
  lgd <- legend(
    "topright",
    legend = c(
      "Observed",
      "NFDS expectation",
      "Drift-smoothed",
      "Missing to optimum"
    ),
    col    = c("black", "blue", "red", "tomato"),
    lwd    = c(2, 2, 2, NA),
    lty    = c(1, 2, 1, NA),
    pch    = c(NA, NA, NA, 22),
    pt.bg  = c(NA, NA, NA, adjustcolor("tomato", alpha.f = 0.3)),
    pt.cex = 1.5,
    bty    = "n",
    plot   = FALSE
  )

  n_obs     <- length(obs_freq)
  n_missing <- max(0L, species_optimum - n_obs)

  if (n_missing > 0) {

    last_x   <- max(bp) + 0.5
    zone_mid <- (last_x + species_optimum) / 2
    rect_top <- lgd$rect$top - lgd$rect$h  # bottom edge of legend bounding box

    rect(
      last_x, 0, species_optimum, rect_top,
      col    = adjustcolor("tomato", alpha.f = 0.15),
      border = adjustcolor("tomato", alpha.f = 0.5),
      lty    = 2
    )

    text(
      x      = zone_mid,
      y      = rect_top / 2,
      labels = paste0(n_missing, " alleles\nto optimum"),
      col    = "tomato4",
      cex    = 0.85,
      font   = 2
    )

  }

  abline(
    h = exp_freq[1],
    col = "blue",
    lwd = 2,
    lty = 2
  )

  lines(
    seq_along(obs_freq),
    obs_freq,
    lwd = 2
  )

  lo <- loess(obs_freq ~ seq_along(obs_freq))

  lines(
    seq_along(obs_freq),
    predict(lo),
    col = "red",
    lwd = 2
  )

  legend(
    "topright",

    legend = c(
      "Observed",
      "NFDS expectation",
      "Drift-smoothed",
      "Missing to optimum"
    ),

    col = c(
      "black",
      "blue",
      "red",
      "tomato"
    ),

    lwd  = c(2, 2, 2, NA),
    lty  = c(1, 2, 1, NA),
    pch  = c(NA, NA, NA, 22),
    pt.bg  = c(NA, NA, NA, adjustcolor("tomato", alpha.f = 0.3)),
    pt.cex = 1.5,
    bty  = "n"
  )
}

dev.off()

cat("\nChi-square analysis complete\n")
