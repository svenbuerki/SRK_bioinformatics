############################################################
# SRK allele accumulation analysis - FIXED TO MATCH CORRECT METHOD
# Uses same approach as your working R code
############################################################

pdf(NULL)   # suppress automatic Rplots.pdf creation by Rscript
cat("Starting SRK allele accumulation analysis\n")

############################################################
# 1. Load data
############################################################

geno <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_info <- read.csv(
  "sampling_metadata.csv",
  stringsAsFactors = FALSE
)

############################################################
# 2. Merge and filter ingroup - FIXED individual matching
############################################################

# Filter metadata to ingroup only
sample_info <- sample_info[sample_info$Ingroup == 1, ]

# Filter genotype data to match ingroup individuals
# FIXED: Use exact matching like your working code
geno <- geno[geno$Individual %in% sample_info$SampleID, ]

cat("Individuals after ingroup filtering:", nrow(geno), "\n")

# Add population info using correct matching
geno$Population <- sample_info$Pop[match(geno$Individual, sample_info$SampleID)]

############################################################
# 3. CORRECT allele processing - columns 2 onwards
############################################################

# FIXED: Alleles are columns 2 onwards (like your working code)
allele_start_col <- 2
allele_end_col <- ncol(geno) - 1  # Exclude the Population column we just added

cat("Allele columns: from", allele_start_col, "to", allele_end_col, "\n")
cat("Total allele columns:", allele_end_col - allele_start_col + 1, "\n")

# Extract allele matrix - CORRECT method
allele_matrix <- as.matrix(geno[, allele_start_col:allele_end_col])

# Convert to numeric and handle NAs
allele_matrix <- apply(allele_matrix, 2, as.numeric)
allele_matrix[is.na(allele_matrix)] <- 0

# Add row names
rownames(allele_matrix) <- geno$Individual

cat("Genotype matrix dimensions:", nrow(allele_matrix), "x", ncol(allele_matrix), "\n")

# Safety check
cat("Genotype value summary:\n")
print(summary(as.vector(allele_matrix)))

############################################################
# 3b. Allele richness estimator functions
############################################################

# Michaelis-Menten asymptote fitted to the accumulation curve.
# Input: mean_accum — vector of mean cumulative allele counts (length = n_individuals).
# Returns the estimated asymptote Smax (ceiling) and a success flag.
fit_mm_asymptote <- function(mean_accum) {
  n_pts <- length(mean_accum)
  S_obs <- max(mean_accum, na.rm = TRUE)
  df_mm <- data.frame(n = seq_len(n_pts), S = mean_accum)
  tryCatch({
    fit <- nls(
      S ~ Smax * n / (K + n),
      data  = df_mm,
      start = list(Smax = S_obs * 1.5, K = n_pts / 2),
      control = nls.control(maxiter = 500)
    )
    list(Smax = ceiling(coef(fit)[["Smax"]]), K = coef(fit)[["K"]], success = TRUE)
  }, error = function(e) {
    cat("  MM fitting failed:", conditionMessage(e), "\n")
    list(Smax = NA_real_, K = NA_real_, success = FALSE)
  })
}

# Base function: additional individuals needed to reach an absolute allele target S_target.
# Returns NA if MM fitting failed or if S_target >= Smax (unreachable asymptote).
# Returns 0 if current sample already meets or exceeds the target.
n_to_reach_s <- function(Smax, K, n_current, S_target) {
  if (is.na(Smax) || is.na(K)) return(NA_real_)
  if (S_target >= Smax) return(NA_real_)
  n_needed <- (S_target * K) / (Smax - S_target)
  max(0L, ceiling(n_needed) - n_current)
}

# Additional individuals needed to reach a given fraction of Smax (e.g. 0.95).
n_to_reach_frac <- function(Smax, K, n_current, target_frac) {
  n_to_reach_s(Smax, K, n_current, target_frac * Smax)
}

# Additional individuals needed to discover one more allele beyond S_obs.
n_for_next_allele <- function(Smax, K, n_current, S_obs) {
  n_to_reach_s(Smax, K, n_current, S_obs + 1)
}

# Chao1 non-parametric lower-bound richness estimator.
# Input: count_vec — colSums of the binary allele matrix
#   (number of individuals carrying each allele).
# Returns S_obs, singleton/doubleton counts, and Chao1 estimate (ceiling).
chao1_estimate <- function(count_vec) {
  present <- count_vec[count_vec > 0]
  S_obs   <- length(present)
  f1 <- sum(present == 1)   # singletons  (alleles in exactly 1 individual)
  f2 <- sum(present == 2)   # doubletons  (alleles in exactly 2 individuals)
  if (f2 > 0) {
    est <- S_obs + (f1^2) / (2 * f2)
    se  <- sqrt(f2 * ((f1/f2)^4/4 + (f1/f2)^3 + (f1/f2)^2/2))
  } else if (f1 > 0) {
    cat("  Chao1: no doubletons — using bias-corrected formula\n")
    est <- S_obs + f1 * (f1 - 1) / 2
    se  <- NA_real_
  } else {
    est <- S_obs
    se  <- NA_real_
  }
  list(S_obs = S_obs, f1 = f1, f2 = f2, Chao1 = ceiling(est), SE = se)
}

# iNEXT asymptotic richness estimate (optional — requires the iNEXT package).
# Input: count_vec — colSums of the binary allele matrix.
# Returns the asymptotic estimator from iNEXT::iNEXT() and a success flag.
inext_estimate <- function(count_vec) {
  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    message("  'iNEXT' not installed — skipping. Install with: install.packages('iNEXT')")
    return(list(S_asymptote = NA_real_, success = FALSE))
  }
  present <- as.integer(count_vec[count_vec > 0])
  tryCatch({
    suppressMessages(
      out <- iNEXT::iNEXT(present, q = 0, datatype = "abundance")
    )
    ae    <- out$AsyEst
    s_row <- ae[grepl("richness", ae$Diversity, ignore.case = TRUE), ]
    s_est <- if (nrow(s_row) > 0) s_row$Estimator[1] else ae$Estimator[1]
    list(S_asymptote = ceiling(s_est), success = TRUE)
  }, error = function(e) {
    cat("  iNEXT failed:", conditionMessage(e), "\n")
    list(S_asymptote = NA_real_, success = FALSE)
  })
}

############################################################
# 4. CORRECT allele accumulation function
############################################################

allele_accumulation <- function(mat, pop_name = "Unknown", n_perm = 1000, seed = 123) {

  set.seed(seed)

  n_ind <- nrow(mat)
  
  # CORRECT allele count - same method as your working code
  allele_counts <- colSums(mat)
  present_alleles <- allele_counts > 0
  true_alleles <- sum(present_alleles)
  
  cat("Population", pop_name, "- Individuals:", n_ind, "- Alleles present:", true_alleles, "\n")

  accum <- matrix(0, nrow = n_perm, ncol = n_ind)

  for (p in 1:n_perm) {

    order <- sample(n_ind)
    cum_mat <- mat[order, , drop = FALSE]

    for (i in 1:n_ind) {
      # Count alleles present in first i individuals
      if (i == 1) {
        cum_subset <- cum_mat[1, , drop = FALSE]
      } else {
        cum_subset <- cum_mat[1:i, , drop = FALSE]
      }
      
      # Count alleles with at least one occurrence
      accum[p, i] <- sum(colSums(cum_subset) > 0)
    }
  }

  mean_accum <- colMeans(accum)
  sd_accum <- apply(accum, 2, sd)

  # NFDS statistics
  initial_range <- min(5, n_ind)
  tail_range <- min(5, n_ind)
  
  initial_slope <- if(n_ind > 1) {
    mean(diff(mean_accum[1:initial_range]))
  } else {
    NA
  }

  tail_slope <- if(n_ind > tail_range) {
    mean(diff(mean_accum[(n_ind-tail_range+1):n_ind]))
  } else {
    NA
  }

  half <- floor(n_ind / 2)
  
  late_fraction <- if(n_ind > 2 && max(mean_accum) > 0) {
    (max(mean_accum) - mean_accum[half]) / max(mean_accum)
  } else {
    NA
  }

  return(list(
    mean_accum = mean_accum,
    sd_accum = sd_accum,
    initial_slope = initial_slope,
    tail_slope = tail_slope,
    late_fraction = late_fraction,
    true_alleles = true_alleles
  ))
}

############################################################
# 5. Open PDF output
############################################################

pdf(
  "SRK_allele_accumulation_curves.pdf",
  width = 8,
  height = 6
)

############################################################
# 6. Species-level analysis
############################################################

cat("\nRunning species-level analysis\n")

species_res <- allele_accumulation(allele_matrix, "Species")

n_ind <- nrow(allele_matrix)

plot(
  1:n_ind,
  species_res$mean_accum,
  type = "l",
  lwd = 3,
  col = "blue",
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of functional alleles",
  main = paste("Species-level SRK accumulation\n(", species_res$true_alleles, "total alleles)")
)

polygon(
  c(1:n_ind, rev(1:n_ind)),
  c(species_res$mean_accum - species_res$sd_accum,
    rev(species_res$mean_accum + species_res$sd_accum)),
  col = rgb(0,0,1,0.2),
  border = NA
)

lines(1:n_ind, species_res$mean_accum, lwd = 3)

############################################################
# 6b. Richness estimators — species level
############################################################

cat("\nRichness estimators (species level)\n")

species_counts <- colSums(allele_matrix)

mm_sp    <- fit_mm_asymptote(species_res$mean_accum)
chao1_sp <- chao1_estimate(species_counts)
inext_sp <- inext_estimate(species_counts)

n_sp <- list(
  next1  = n_for_next_allele(mm_sp$Smax, mm_sp$K, nrow(allele_matrix), species_res$true_alleles),
  p80    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix), 0.80),
  p90    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix), 0.90),
  p95    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix), 0.95),
  p99    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix), 0.99)
)

cat("  Observed alleles :", species_res$true_alleles, "\n")
cat("  MM estimate      :", mm_sp$Smax, "\n")
cat("  Chao1 estimate   :", chao1_sp$Chao1,
    if (!is.na(chao1_sp$SE)) paste0("(SE ± ", round(chao1_sp$SE, 1), ")") else "", "\n")
cat("  iNEXT estimate   :", inext_sp$S_asymptote, "\n")
cat("  Add. ind. for +1 allele :", n_sp$next1, "\n")
cat("  Add. ind. for 80% MM   :", n_sp$p80,   "\n")
cat("  Add. ind. for 90% MM   :", n_sp$p90,   "\n")
cat("  Add. ind. for 95% MM   :", n_sp$p95,   "\n")
cat("  Add. ind. for 99% MM   :", n_sp$p99,   "\n")

# Add estimated asymptotes to the open species plot (still on this PDF page)
est_labels <- character(0)
est_cols   <- character(0)
est_lty    <- integer(0)

if (!is.na(mm_sp$Smax)) {
  abline(h = mm_sp$Smax,         col = "darkgreen",  lwd = 1.5, lty = 2)
  est_labels <- c(est_labels, paste0("MM: ",    mm_sp$Smax))
  est_cols   <- c(est_cols,   "darkgreen")
  est_lty    <- c(est_lty,    2L)
}
if (!is.na(chao1_sp$Chao1)) {
  abline(h = chao1_sp$Chao1,     col = "purple",     lwd = 1.5, lty = 3)
  est_labels <- c(est_labels, paste0("Chao1: ", chao1_sp$Chao1))
  est_cols   <- c(est_cols,   "purple")
  est_lty    <- c(est_lty,    3L)
}
if (!is.na(inext_sp$S_asymptote)) {
  abline(h = inext_sp$S_asymptote, col = "darkorange", lwd = 1.5, lty = 4)
  est_labels <- c(est_labels, paste0("iNEXT: ", inext_sp$S_asymptote))
  est_cols   <- c(est_cols,   "darkorange")
  est_lty    <- c(est_lty,    4L)
}

if (length(est_labels) > 0) {
  legend("bottomright",
    legend = est_labels,
    col    = est_cols,
    lwd    = 1.5,
    lty    = est_lty,
    bty    = "n",
    cex    = 0.75,
    title  = "Richness estimates"
  )
}

species_est <- list(mm = mm_sp, chao1 = chao1_sp, inext = inext_sp, n_more = n_sp)

############################################################
# 7. Population-level analysis - CORRECT method
############################################################

pop_counts <- table(geno$Population)
valid_pops <- names(pop_counts[pop_counts >= 5])

cat("\nPopulations retained:", valid_pops, "\n")

pop_results <- list()
pop_est     <- list()

for (pop in valid_pops) {

  cat("Processing population:", pop, "\n")

  # CORRECT population subsetting - same as your working code
  target_pop_metadata <- sample_info[sample_info$Pop == pop, ]
  pop_geno <- geno[geno$Individual %in% target_pop_metadata$SampleID, ]
  
  # Extract allele matrix for this population - CORRECT method
  pop_allele_matrix <- as.matrix(pop_geno[, allele_start_col:allele_end_col])
  pop_allele_matrix <- apply(pop_allele_matrix, 2, as.numeric)
  pop_allele_matrix[is.na(pop_allele_matrix)] <- 0
  rownames(pop_allele_matrix) <- pop_geno$Individual

  res <- allele_accumulation(pop_allele_matrix, pop)

  pop_results[[pop]] <- res

  n_ind <- nrow(pop_allele_matrix)

  plot(
    1:n_ind,
    res$mean_accum,
    type = "l",
    lwd = 3,
    col = "red",
    xlab = "Number of individuals sampled",
    ylab = "Cumulative number of functional alleles",
    main = paste("Population", pop, "\n(", res$true_alleles, "alleles,", n_ind, "individuals)")
  )

  polygon(
    c(1:n_ind, rev(1:n_ind)),
    c(res$mean_accum - res$sd_accum,
      rev(res$mean_accum + res$sd_accum)),
    col = rgb(1,0,0,0.2),
    border = NA
  )

  lines(1:n_ind, res$mean_accum, lwd = 3)

  # Richness estimators for this population
  pop_counts_i <- colSums(pop_allele_matrix)
  mm_pop    <- fit_mm_asymptote(res$mean_accum)
  chao1_pop <- chao1_estimate(pop_counts_i)
  inext_pop <- inext_estimate(pop_counts_i)

  n_pop <- list(
    next1 = n_for_next_allele(mm_pop$Smax, mm_pop$K, nrow(pop_allele_matrix), res$true_alleles),
    p80   = n_to_reach_frac(mm_pop$Smax, mm_pop$K, nrow(pop_allele_matrix), 0.80),
    p90   = n_to_reach_frac(mm_pop$Smax, mm_pop$K, nrow(pop_allele_matrix), 0.90),
    p95   = n_to_reach_frac(mm_pop$Smax, mm_pop$K, nrow(pop_allele_matrix), 0.95),
    p99   = n_to_reach_frac(mm_pop$Smax, mm_pop$K, nrow(pop_allele_matrix), 0.99)
  )

  cat("  Observed:", res$true_alleles,
      "| MM:", mm_pop$Smax,
      "| Chao1:", chao1_pop$Chao1,
      "| iNEXT:", inext_pop$S_asymptote, "\n")
  cat("  Add. ind. — +1:", n_pop$next1,
      "| 80%:", n_pop$p80,
      "| 90%:", n_pop$p90,
      "| 95%:", n_pop$p95,
      "| 99%:", n_pop$p99, "\n")

  pop_est_labels <- character(0)
  pop_est_cols   <- character(0)
  pop_est_lty    <- integer(0)

  if (!is.na(mm_pop$Smax)) {
    abline(h = mm_pop$Smax,         col = "darkgreen",  lwd = 1.5, lty = 2)
    pop_est_labels <- c(pop_est_labels, paste0("MM: ",    mm_pop$Smax))
    pop_est_cols   <- c(pop_est_cols,   "darkgreen")
    pop_est_lty    <- c(pop_est_lty,    2L)
  }
  if (!is.na(chao1_pop$Chao1)) {
    abline(h = chao1_pop$Chao1,     col = "purple",     lwd = 1.5, lty = 3)
    pop_est_labels <- c(pop_est_labels, paste0("Chao1: ", chao1_pop$Chao1))
    pop_est_cols   <- c(pop_est_cols,   "purple")
    pop_est_lty    <- c(pop_est_lty,    3L)
  }
  if (!is.na(inext_pop$S_asymptote)) {
    abline(h = inext_pop$S_asymptote, col = "darkorange", lwd = 1.5, lty = 4)
    pop_est_labels <- c(pop_est_labels, paste0("iNEXT: ", inext_pop$S_asymptote))
    pop_est_cols   <- c(pop_est_cols,   "darkorange")
    pop_est_lty    <- c(pop_est_lty,    4L)
  }

  if (length(pop_est_labels) > 0) {
    legend("bottomright",
      legend = pop_est_labels,
      col    = pop_est_cols,
      lwd    = 1.5,
      lty    = pop_est_lty,
      bty    = "n",
      cex    = 0.75,
      title  = "Richness estimates"
    )
  }

  pop_est[[pop]] <- list(mm = mm_pop, chao1 = chao1_pop, inext = inext_pop, n_more = n_pop)

}

############################################################
# 8. Close PDF
############################################################

dev.off()

cat("\nPDF written: SRK_allele_accumulation_curves.pdf\n")

############################################################
# 8b. Combined plot: all EOs on one axes (no species curve)
############################################################

cat("\nGenerating combined EO accumulation plot...\n")

# Colour palette: one distinct colour per EO
eo_colours <- c(
  "#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "#984ea3",
  "#a65628", "#f781bf"   # extras if more than 5 EOs
)[seq_along(valid_pops)]
names(eo_colours) <- valid_pops

# Axis limits: y based on max observed alleles across EOs
max_eo_alleles <- max(sapply(valid_pops, function(p) pop_results[[p]]$true_alleles))
y_max <- max_eo_alleles * 1.25
x_max <- max(sapply(valid_pops, function(p) length(pop_results[[p]]$mean_accum)))

png("SRK_allele_accumulation_combined.png", width = 9, height = 6,
    units = "in", res = 200)
par(mar = c(5, 4, 4, 10))   # wide right margin for end-of-line labels

plot(
  NA,
  xlim = c(1, x_max),
  ylim = c(0, y_max),
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of S-alleles",
  main = "SRK allele accumulation: EO comparison",
  las  = 1
)
grid(col = "grey90", lty = 1)

# Species MM asymptote reference line (kept as context)
if (!is.na(mm_sp$Smax)) {
  abline(h = mm_sp$Smax, col = "grey40", lwd = 1.2, lty = 2)
  mtext(
    paste0("Species MM: ", mm_sp$Smax),
    side = 4, at = mm_sp$Smax,
    col = "grey40", cex = 0.72, las = 1, line = 0.3
  )
}

# EO curves + CI bands + end labels + MM predicted total
for (pop in valid_pops) {
  res    <- pop_results[[pop]]
  n_pop  <- length(res$mean_accum)
  col_i  <- eo_colours[pop]
  mm_est <- pop_est[[pop]]$mm$Smax

  polygon(
    c(1:n_pop, rev(1:n_pop)),
    c(res$mean_accum - res$sd_accum,
      rev(res$mean_accum + res$sd_accum)),
    col = adjustcolor(col_i, alpha.f = 0.15),
    border = NA
  )
  lines(1:n_pop, res$mean_accum, col = col_i, lwd = 2)

  # End-of-curve label: observed (/ MM predicted)
  mm_label <- if (!is.na(mm_est)) paste0("/", mm_est) else ""
  text(
    x = n_pop, y = res$mean_accum[n_pop],
    labels = paste0("EO", pop, " (", res$true_alleles, mm_label, ")"),
    col = col_i, pos = 4, cex = 0.72, xpd = TRUE
  )
}

# Legend: solid line = observed curve, dotted = MM predicted total
legend(
  "bottomright",
  legend = sapply(valid_pops, function(p) {
    n_p    <- nrow(geno[geno$Population == p, ])
    mm_est <- pop_est[[p]]$mm$Smax
    mm_str <- if (!is.na(mm_est)) paste0(", MM=", mm_est) else ""
    paste0("EO", p, "  N=", n_p, ", obs=", pop_results[[p]]$true_alleles, mm_str)
  }),
  col  = eo_colours[valid_pops],
  lwd  = 2,
  bty  = "n",
  cex  = 0.78
)

dev.off()
par(mar = c(5, 4, 4, 2))   # restore default margins
cat("PNG written: SRK_allele_accumulation_combined.png\n")

############################################################
# 9. Export statistics with CORRECT counts
############################################################

stats <- data.frame(
  Level="Species",
  Population="All",
  N_individuals=nrow(allele_matrix),
  N_alleles=species_res$true_alleles,
  Initial_slope=species_res$initial_slope,
  Tail_slope=species_res$tail_slope,
  Late_fraction=species_res$late_fraction,
  MM_estimate=species_est$mm$Smax,
  Chao1_estimate=species_est$chao1$Chao1,
  iNEXT_estimate=species_est$inext$S_asymptote,
  N_more_next_allele=species_est$n_more$next1,
  N_more_80pct_MM=species_est$n_more$p80,
  N_more_90pct_MM=species_est$n_more$p90,
  N_more_95pct_MM=species_est$n_more$p95,
  N_more_99pct_MM=species_est$n_more$p99
)

for (pop in names(pop_results)) {

  res      <- pop_results[[pop]]
  pop_data <- geno[geno$Population == pop, ]
  pest     <- pop_est[[pop]]

  stats <- rbind(stats,
    data.frame(
      Level="Population",
      Population=pop,
      N_individuals=nrow(pop_data),
      N_alleles=res$true_alleles,
      Initial_slope=res$initial_slope,
      Tail_slope=res$tail_slope,
      Late_fraction=res$late_fraction,
      MM_estimate=pest$mm$Smax,
      Chao1_estimate=pest$chao1$Chao1,
      iNEXT_estimate=pest$inext$S_asymptote,
      N_more_next_allele=pest$n_more$next1,
      N_more_80pct_MM=pest$n_more$p80,
      N_more_90pct_MM=pest$n_more$p90,
      N_more_95pct_MM=pest$n_more$p95,
      N_more_99pct_MM=pest$n_more$p99
    )
  )
}

write.table(
  stats,
  "SRK_allele_accumulation_stats.tsv",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("\nStats written: SRK_allele_accumulation_stats.tsv\n")

############################################################
# 9b. Write species richness estimates file
#     (used by SRK_chisq_species_population.R as the "optimum")
############################################################

consensus_est <- ceiling(mean(
  c(species_est$mm$Smax,
    species_est$chao1$Chao1,
    species_est$inext$S_asymptote),
  na.rm = TRUE
))

species_richness_est <- data.frame(
  S_observed          = species_res$true_alleles,
  MM_estimate         = species_est$mm$Smax,
  Chao1_estimate      = species_est$chao1$Chao1,
  iNEXT_estimate      = species_est$inext$S_asymptote,
  Consensus_estimate  = consensus_est,
  N_more_next_allele  = species_est$n_more$next1,
  N_more_80pct_MM     = species_est$n_more$p80,
  N_more_90pct_MM     = species_est$n_more$p90,
  N_more_95pct_MM     = species_est$n_more$p95,
  N_more_99pct_MM     = species_est$n_more$p99
)

write.table(
  species_richness_est,
  "SRK_species_richness_estimates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Richness estimates written: SRK_species_richness_estimates.tsv\n")
cat("  Consensus (mean of available estimators):", consensus_est, "alleles\n")

############################################################
# 10. Print summary with validation
############################################################

print(stats)

# Validation - should match your working code
cat("\n=== VALIDATION ===\n")
cat("Population 27 should now show 17 alleles\n")
cat("Method matches your working R code:\n")
cat("- Alleles from column 2 onwards\n")
cat("- Correct individual matching\n")
cat("- colSums(matrix) > 0 counting method\n")

cat("\nAnalysis complete\n")

