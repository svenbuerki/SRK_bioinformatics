###############################
# SRK protein frequency χ² tests
# Compatible with functional protein pipeline
###############################

cat("\nStarting SRK chi-square analysis\n")

###############################
# Load data
###############################

allele_df <- read.table(
  "SRK_individual_protein_table.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

meta_df <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

###############################
# Detect allele/protein column automatically
###############################

if ("Protein" %in% colnames(allele_df)) {

  allele_col <- "Protein"
  cat("Detected functional protein column\n")

} else if ("Allele" %in% colnames(allele_df)) {

  allele_col <- "Allele"
  cat("Detected phylogenetic allele column\n")

} else {

  stop("ERROR: No Allele or Protein column found")

}

###############################
# Keep ingroup only
###############################

meta_df <- meta_df[meta_df$Ingroup == 1, ]

allele_df <- allele_df[
  allele_df$Individual %in% meta_df$SampleID,
]

allele_df$Population <- meta_df$Pop[
  match(allele_df$Individual, meta_df$SampleID)
]

cat("Individuals analyzed:",
    length(unique(allele_df$Individual)), "\n")

cat("Unique proteins:",
    length(unique(allele_df[[allele_col]])), "\n")

###############################
# Safe chi-square function
###############################

run_chisq <- function(df, allele_col) {

  allele_counts <- table(df[[allele_col]])

  if (length(allele_counts) < 2) {

    return(list(

      stats = data.frame(

        N_individuals = length(unique(df$Individual)),
        N_alleles = length(allele_counts),
        X2 = NA,
        df = NA,
        p_value = NA
      ),

      counts = allele_counts
    ))
  }

  obs <- as.numeric(allele_counts)
  exp_prob <- rep(1 / length(obs), length(obs))

  chisq <- suppressWarnings(
    chisq.test(x = obs, p = exp_prob)
  )

  list(

    stats = data.frame(

      N_individuals = length(unique(df$Individual)),
      N_alleles = length(obs),
      X2 = as.numeric(chisq$statistic),
      df = as.numeric(chisq$parameter),
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

res_species <- run_chisq(allele_df, allele_col)

stats_species <- res_species$stats
stats_species$Level <- "Species"
stats_species$Population <- "All"

results_list[["Species_All"]] <- stats_species
plot_list[["Species_All"]] <- res_species$counts

###############################
# Population-level
###############################

pop_sizes <- table(allele_df$Population)

valid_pops <- names(pop_sizes[pop_sizes >= 5])

cat("Populations analyzed:", length(valid_pops), "\n")

for (pop in valid_pops) {

  df_pop <- allele_df[
    allele_df$Population == pop,
  ]

  res <- run_chisq(df_pop, allele_col)

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

  barplot(

    obs_freq,

    border = NA,
    col = "gray70",

    ylim = c(
      0,
      max(obs_freq) * 1.2
    ),

    main = paste(level, pop, sep = " – "),
    ylab = "Protein frequency",
    xlab = "Proteins (ranked)"
  )

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
      "Drift-smoothed"
    ),

    col = c(
      "black",
      "blue",
      "red"
    ),

    lwd = 2,
    lty = c(1,2,1),
    bty = "n"
  )
}

dev.off()

cat("\nChi-square analysis complete\n")
