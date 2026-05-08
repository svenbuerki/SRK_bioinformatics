############################################################
# SRK allele accumulation analysis - Step 15 of Canu_amplicon pipeline
#
# Three levels of analysis are produced:
#   1. Species-level     accumulation curve + MM/Chao1/iNEXT estimators
#                        Output: species optimum used by TP1 (Step 17).
#   2. Bottleneck Lineage (BL1-BL5) accumulation curves
#                        Tests whether independent demographic lineages are
#                        approaching saturation independently.
#   3. Element Occurrence (EO) accumulation curves
#                        For all EOs with N >= 5 individuals; curves colored
#                        by parent BL on the combined plot.
############################################################

pdf(NULL)   # suppress automatic Rplots.pdf creation by Rscript
cat("Starting SRK allele accumulation analysis (Step 15)\n")

# Shared BL color palette (RColorBrewer Set1, locked at Step 14 - must
# match across all Phase 3/4 plots and matches the sibling
# LEPA_EO_spatial_clustering project).
BL_PALETTE <- c(
  BL1 = "#E41A1C",   # red
  BL2 = "#377EB8",   # blue
  BL3 = "#4DAF4A",   # green
  BL4 = "#984EA3",   # purple
  BL5 = "#FF7F00"    # orange
)
BL_UNASSIGNED_COLOR <- "#999999"

############################################################
# 1. Load data
############################################################

geno_all <- read.table(
  "SRK_individual_allele_genotypes.tsv",
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)

sample_info <- read.csv("sampling_metadata.csv", stringsAsFactors = FALSE)

bl_assign <- read.table(
  "SRK_individual_BL_assignments.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

cat("Loaded BL assignments for", nrow(bl_assign), "individuals\n")

############################################################
# 2. Two filtered datasets
#
#   geno_all  -- ALL ingroup individuals (272). Used for the species-level
#                accumulation curve and richness estimators. The species
#                pool must include every individual we have evidence for,
#                including those whose Pop code could not be resolved to
#                a known EO (10 unresolved germplasm sub-codes carry
#                ~4 alleles not seen in the BL-assigned subset).
#
#   geno      -- BL-assigned subset (262). Used for the EO and BL loops,
#                where geographic provenance is required.
############################################################

# Species-level filter: keep all ingroup individuals
sample_info <- sample_info[sample_info$Ingroup == 1, ]
geno_all    <- geno_all[geno_all$Individual %in% sample_info$SampleID, ]
cat("Species-level individuals (all ingroup):", nrow(geno_all), "\n")

# EO/BL filter: keep only BL-assigned individuals
bl_ok <- bl_assign[bl_assign$BL_status %in% c("Assigned", "Inferred"), ]
geno  <- geno_all[geno_all$Individual %in% bl_ok$Individual, ]
geno$EO <- bl_ok$EO[match(geno$Individual, bl_ok$Individual)]
geno$BL <- bl_ok$BL[match(geno$Individual, bl_ok$Individual)]
cat("EO/BL-level individuals (BL-assigned):", nrow(geno), "\n")

############################################################
# 3. Allele matrices
############################################################

allele_cols <- setdiff(colnames(geno_all), "Individual")
cat("Total allele columns:", length(allele_cols), "\n")

build_matrix <- function(df) {
  m <- as.matrix(df[, allele_cols, drop = FALSE])
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  rownames(m) <- df$Individual
  m
}

allele_matrix_all <- build_matrix(geno_all)   # species-level
allele_matrix     <- build_matrix(geno)       # EO/BL-level

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

species_res <- allele_accumulation(allele_matrix_all, "Species")

n_ind <- nrow(allele_matrix_all)

plot(
  1:n_ind,
  species_res$mean_accum,
  type = "l",
  lwd = 3,
  col = "blue",
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of functional alleles",
  main = paste("Species-level S-allele accumulation\n(", species_res$true_alleles, "total alleles)")
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

species_counts <- colSums(allele_matrix_all)

mm_sp    <- fit_mm_asymptote(species_res$mean_accum)
chao1_sp <- chao1_estimate(species_counts)
inext_sp <- inext_estimate(species_counts)

n_sp <- list(
  next1  = n_for_next_allele(mm_sp$Smax, mm_sp$K, nrow(allele_matrix_all), species_res$true_alleles),
  p80    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix_all), 0.80),
  p90    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix_all), 0.90),
  p95    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix_all), 0.95),
  p99    = n_to_reach_frac(mm_sp$Smax, mm_sp$K, nrow(allele_matrix_all), 0.99)
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
# 6c. Species accumulation as a standalone PNG (figures/)
############################################################

png("figures/SRK_allele_accumulation_species.png",
    width = 9, height = 6, units = "in", res = 200)
par(mar = c(5, 4.5, 4, 2))

# y-axis must accommodate the asymptote lines, not just the observed curve
y_max_species <- max(
  c(species_res$mean_accum + species_res$sd_accum,
    mm_sp$Smax, chao1_sp$Chao1, inext_sp$S_asymptote),
  na.rm = TRUE
) * 1.05

plot(
  1:n_ind, species_res$mean_accum,
  type = "l", lwd = 3, col = "#2166ac",
  ylim = c(0, y_max_species),
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of S-alleles",
  main = paste0("Species-level S-allele accumulation (",
                species_res$true_alleles, " observed in ", n_ind,
                " individuals)"),
  las = 1
)
grid(col = "grey90", lty = 1)
polygon(
  c(1:n_ind, rev(1:n_ind)),
  c(species_res$mean_accum - species_res$sd_accum,
    rev(species_res$mean_accum + species_res$sd_accum)),
  col = adjustcolor("#2166ac", alpha.f = 0.2), border = NA
)
lines(1:n_ind, species_res$mean_accum, lwd = 3, col = "#2166ac")

png_labels <- character(0); png_cols <- character(0); png_lty <- integer(0)
if (!is.na(mm_sp$Smax)) {
  abline(h = mm_sp$Smax, col = "darkgreen", lwd = 1.5, lty = 2)
  png_labels <- c(png_labels, paste0("MM: ", mm_sp$Smax))
  png_cols   <- c(png_cols, "darkgreen"); png_lty <- c(png_lty, 2L)
}
if (!is.na(chao1_sp$Chao1)) {
  abline(h = chao1_sp$Chao1, col = "purple", lwd = 1.5, lty = 3)
  png_labels <- c(png_labels, paste0("Chao1: ", chao1_sp$Chao1))
  png_cols   <- c(png_cols, "purple"); png_lty <- c(png_lty, 3L)
}
if (!is.na(inext_sp$S_asymptote)) {
  abline(h = inext_sp$S_asymptote, col = "darkorange", lwd = 1.5, lty = 4)
  png_labels <- c(png_labels, paste0("iNEXT: ", inext_sp$S_asymptote))
  png_cols   <- c(png_cols, "darkorange"); png_lty <- c(png_lty, 4L)
}
if (length(png_labels) > 0) {
  legend("bottomright", legend = png_labels, col = png_cols,
         lwd = 1.5, lty = png_lty, bty = "n", cex = 0.85,
         title = "Asymptotic richness estimates")
}
dev.off()
par(mar = c(5, 4, 4, 2))
cat("PNG written: figures/SRK_allele_accumulation_species.png\n")

############################################################
# 7. Group-level analysis (EO and BL)
############################################################
# Reusable per-group routine. `level_label` is "EO" or "BL"; it controls
# the curve color and panel title only - the underlying logic is identical.

run_group_accumulation <- function(group_label, member_idx, panel_color,
                                   level_label) {
  group_mat <- allele_matrix[member_idx, , drop = FALSE]
  cat("Processing", level_label, group_label, "(", nrow(group_mat),
      "individuals)\n")

  res    <- allele_accumulation(group_mat, group_label)
  n_ind  <- nrow(group_mat)

  # Per-group page in the main PDF
  plot(
    1:n_ind, res$mean_accum,
    type = "l", lwd = 3, col = panel_color,
    xlab = "Number of individuals sampled",
    ylab = "Cumulative number of functional alleles",
    main = paste0(level_label, " ", group_label,
                  "  (", res$true_alleles, " alleles, ", n_ind, " individuals)")
  )
  polygon(
    c(1:n_ind, rev(1:n_ind)),
    c(res$mean_accum - res$sd_accum,
      rev(res$mean_accum + res$sd_accum)),
    col = adjustcolor(panel_color, alpha.f = 0.2),
    border = NA
  )
  lines(1:n_ind, res$mean_accum, lwd = 3, col = panel_color)

  group_counts <- colSums(group_mat)
  mm_g    <- fit_mm_asymptote(res$mean_accum)
  chao1_g <- chao1_estimate(group_counts)
  inext_g <- inext_estimate(group_counts)

  n_more <- list(
    next1 = n_for_next_allele(mm_g$Smax, mm_g$K, n_ind, res$true_alleles),
    p80   = n_to_reach_frac(mm_g$Smax, mm_g$K, n_ind, 0.80),
    p90   = n_to_reach_frac(mm_g$Smax, mm_g$K, n_ind, 0.90),
    p95   = n_to_reach_frac(mm_g$Smax, mm_g$K, n_ind, 0.95),
    p99   = n_to_reach_frac(mm_g$Smax, mm_g$K, n_ind, 0.99)
  )

  cat("  Observed:", res$true_alleles,
      "| MM:",      mm_g$Smax,
      "| Chao1:",   chao1_g$Chao1,
      "| iNEXT:",   inext_g$S_asymptote, "\n")

  # Asymptote reference lines on the per-group panel
  est_labels <- character(0); est_cols <- character(0); est_lty <- integer(0)
  if (!is.na(mm_g$Smax)) {
    abline(h = mm_g$Smax, col = "darkgreen", lwd = 1.5, lty = 2)
    est_labels <- c(est_labels, paste0("MM: ", mm_g$Smax))
    est_cols   <- c(est_cols,   "darkgreen"); est_lty <- c(est_lty, 2L)
  }
  if (!is.na(chao1_g$Chao1)) {
    abline(h = chao1_g$Chao1, col = "purple", lwd = 1.5, lty = 3)
    est_labels <- c(est_labels, paste0("Chao1: ", chao1_g$Chao1))
    est_cols   <- c(est_cols,   "purple"); est_lty <- c(est_lty, 3L)
  }
  if (!is.na(inext_g$S_asymptote)) {
    abline(h = inext_g$S_asymptote, col = "darkorange", lwd = 1.5, lty = 4)
    est_labels <- c(est_labels, paste0("iNEXT: ", inext_g$S_asymptote))
    est_cols   <- c(est_cols,   "darkorange"); est_lty <- c(est_lty, 4L)
  }
  if (length(est_labels) > 0) {
    legend("bottomright", legend = est_labels, col = est_cols,
           lwd = 1.5, lty = est_lty, bty = "n", cex = 0.75,
           title = "Richness estimates")
  }

  list(
    accum  = res,
    est    = list(mm = mm_g, chao1 = chao1_g, inext = inext_g, n_more = n_more)
  )
}

# ---- 7a. EO-level (only EOs with N >= 5 individuals) ----

eo_counts <- table(geno$EO)
valid_eos <- names(eo_counts[eo_counts >= 5])
# Sort EOs by parent BL then alphabetically within BL
eo_to_bl <- unique(geno[, c("EO", "BL")])
eo_to_bl <- eo_to_bl[order(eo_to_bl$BL, eo_to_bl$EO), ]
valid_eos <- intersect(eo_to_bl$EO, valid_eos)
cat("\nEOs retained (N >= 5, sorted by BL):",
    paste(valid_eos, collapse = ", "), "\n")

eo_results <- list()
eo_est     <- list()
for (eo in valid_eos) {
  bl_for_eo <- eo_to_bl$BL[eo_to_bl$EO == eo]
  eo_color  <- BL_PALETTE[bl_for_eo]
  if (is.na(eo_color)) eo_color <- BL_UNASSIGNED_COLOR
  out <- run_group_accumulation(eo, which(geno$EO == eo), eo_color, "EO")
  eo_results[[eo]] <- out$accum
  eo_est[[eo]]     <- out$est
}

# ---- 7b. BL-level (all 5 BLs) ----

valid_bls <- sort(unique(geno$BL))
cat("\nBLs retained:", paste(valid_bls, collapse = ", "), "\n")

bl_results <- list()
bl_est     <- list()
for (bl in valid_bls) {
  bl_color <- BL_PALETTE[bl]
  if (is.na(bl_color)) bl_color <- BL_UNASSIGNED_COLOR
  out <- run_group_accumulation(bl, which(geno$BL == bl), bl_color, "BL")
  bl_results[[bl]] <- out$accum
  bl_est[[bl]]     <- out$est
}

############################################################
# 8. Close PDF
############################################################

dev.off()

cat("\nPDF written: SRK_allele_accumulation_curves.pdf\n")

############################################################
# 8b. Combined plot: EOs on one axes, colored by parent BL
############################################################

cat("\nGenerating combined EO accumulation plot (colored by BL)...\n")

# Each EO inherits its parent BL's color from the locked palette
eo_bl <- setNames(eo_to_bl$BL, eo_to_bl$EO)
eo_colours <- BL_PALETTE[eo_bl[valid_eos]]
eo_colours[is.na(eo_colours)] <- BL_UNASSIGNED_COLOR
names(eo_colours) <- valid_eos

max_eo_alleles <- max(sapply(valid_eos, function(p) eo_results[[p]]$true_alleles))
y_max <- max_eo_alleles * 1.25
x_max <- max(sapply(valid_eos, function(p) length(eo_results[[p]]$mean_accum)))

png("figures/SRK_allele_accumulation_combined.png", width = 9, height = 6,
    units = "in", res = 200)
par(mar = c(5, 4, 4, 10))

plot(
  NA, xlim = c(1, x_max), ylim = c(0, y_max),
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of S-alleles",
  main = "S-allele accumulation per EO (colored by parent BL)",
  las  = 1
)
grid(col = "grey90", lty = 1)

if (!is.na(mm_sp$Smax)) {
  abline(h = mm_sp$Smax, col = "grey40", lwd = 1.2, lty = 2)
  mtext(paste0("Species MM: ", mm_sp$Smax), side = 4, at = mm_sp$Smax,
        col = "grey40", cex = 0.72, las = 1, line = 0.3)
}

for (eo in valid_eos) {
  res    <- eo_results[[eo]]
  n_pop  <- length(res$mean_accum)
  col_i  <- eo_colours[eo]
  mm_est <- eo_est[[eo]]$mm$Smax

  polygon(
    c(1:n_pop, rev(1:n_pop)),
    c(res$mean_accum - res$sd_accum, rev(res$mean_accum + res$sd_accum)),
    col = adjustcolor(col_i, alpha.f = 0.15), border = NA
  )
  lines(1:n_pop, res$mean_accum, col = col_i, lwd = 2)

  mm_label <- if (!is.na(mm_est)) paste0("/", mm_est) else ""
  text(x = n_pop, y = res$mean_accum[n_pop],
       labels = paste0(eo, " (", res$true_alleles, mm_label, ")"),
       col = col_i, pos = 4, cex = 0.72, xpd = TRUE)
}

# Legend: one row per EO, grouped visually by BL color
legend(
  "bottomright",
  legend = sapply(valid_eos, function(p) {
    n_p    <- sum(geno$EO == p)
    mm_est <- eo_est[[p]]$mm$Smax
    mm_str <- if (!is.na(mm_est)) paste0(", MM=", mm_est) else ""
    paste0(eo_bl[p], "  ", p, "  N=", n_p,
           ", obs=", eo_results[[p]]$true_alleles, mm_str)
  }),
  col = eo_colours[valid_eos], lwd = 2, bty = "n", cex = 0.78
)

dev.off()
par(mar = c(5, 4, 4, 2))
cat("PNG written: figures/SRK_allele_accumulation_combined.png\n")

############################################################
# 8b-ii. Combined plot: BL aggregates on one axes
############################################################

cat("Generating combined BL accumulation plot...\n")

bl_colours <- BL_PALETTE[valid_bls]
bl_colours[is.na(bl_colours)] <- BL_UNASSIGNED_COLOR
names(bl_colours) <- valid_bls

max_bl_alleles <- max(sapply(valid_bls, function(b) bl_results[[b]]$true_alleles))
y_max_bl <- max_bl_alleles * 1.25
x_max_bl <- max(sapply(valid_bls, function(b) length(bl_results[[b]]$mean_accum)))

png("figures/SRK_allele_accumulation_BL_combined.png", width = 9, height = 6,
    units = "in", res = 200)
par(mar = c(5, 4, 4, 10))

plot(
  NA, xlim = c(1, x_max_bl), ylim = c(0, y_max_bl),
  xlab = "Number of individuals sampled",
  ylab = "Cumulative number of S-alleles",
  main = "S-allele accumulation per bottleneck lineage (BL)",
  las  = 1
)
grid(col = "grey90", lty = 1)

if (!is.na(mm_sp$Smax)) {
  abline(h = mm_sp$Smax, col = "grey40", lwd = 1.2, lty = 2)
  mtext(paste0("Species MM: ", mm_sp$Smax), side = 4, at = mm_sp$Smax,
        col = "grey40", cex = 0.72, las = 1, line = 0.3)
}

for (bl in valid_bls) {
  res    <- bl_results[[bl]]
  n_pop  <- length(res$mean_accum)
  col_i  <- bl_colours[bl]
  mm_est <- bl_est[[bl]]$mm$Smax

  polygon(
    c(1:n_pop, rev(1:n_pop)),
    c(res$mean_accum - res$sd_accum, rev(res$mean_accum + res$sd_accum)),
    col = adjustcolor(col_i, alpha.f = 0.15), border = NA
  )
  lines(1:n_pop, res$mean_accum, col = col_i, lwd = 2.5)

  mm_label <- if (!is.na(mm_est)) paste0("/", mm_est) else ""
  text(x = n_pop, y = res$mean_accum[n_pop],
       labels = paste0(bl, " (", res$true_alleles, mm_label, ")"),
       col = col_i, pos = 4, cex = 0.78, xpd = TRUE, font = 2)
}

legend(
  "bottomright",
  legend = sapply(valid_bls, function(b) {
    n_b    <- sum(geno$BL == b)
    mm_est <- bl_est[[b]]$mm$Smax
    mm_str <- if (!is.na(mm_est)) paste0(", MM=", mm_est) else ""
    paste0(b, "  N=", n_b, ", obs=",
           bl_results[[b]]$true_alleles, mm_str)
  }),
  col = bl_colours[valid_bls], lwd = 2.5, bty = "n", cex = 0.85
)

dev.off()
par(mar = c(5, 4, 4, 2))
cat("PNG written: figures/SRK_allele_accumulation_BL_combined.png\n")

############################################################
# 8c. Stacked bar: observed / predicted undetected / lost to drift
############################################################

cat("\nGenerating drift erosion stacked bar chart...\n")

# Species MM from the freshly computed estimate (same run)
species_MM <- mm_sp$Smax

# Per-EO data: observed alleles, EO-level MM predicted total, parent BL
# color (used as a thin baseline strip for visual grouping).
eo_df <- data.frame(
  EO       = valid_eos,
  observed = sapply(valid_eos, function(p) eo_results[[p]]$true_alleles),
  mm_eo    = sapply(valid_eos, function(p) {
    mm <- eo_est[[p]]$mm$Smax
    if (is.na(mm)) eo_results[[p]]$true_alleles else mm
  }),
  bar_color = unname(eo_colours[valid_eos]),
  stringsAsFactors = FALSE
)

# Three stacked components
eo_df$predicted_extra <- pmax(0, eo_df$mm_eo    - eo_df$observed)
eo_df$lost_to_drift   <- pmax(0, species_MM      - eo_df$mm_eo)

# Keep BL-sorted order from valid_eos (do NOT re-sort by observed - we
# want EOs grouped by their parent BL so the BL strip below the bars is
# visually contiguous).

col_observed  <- "#2166ac"   # dark blue  — detected
col_predicted <- "#92c5de"   # light blue — predicted undetected
col_lost      <- "#d73027"   # red        — lost to genetic drift

png("figures/SRK_allele_accumulation_drift_erosion.png",
    width = 7, height = 6, units = "in", res = 200)

par(mar = c(5, 5, 7, 5))   # extra top margin for legend

bar_mat <- rbind(
  eo_df$observed,
  eo_df$predicted_extra,
  eo_df$lost_to_drift
)

bp <- barplot(
  bar_mat,
  beside    = FALSE,
  col       = c(col_observed, col_predicted, col_lost),
  names.arg = eo_df$EO,
  ylim      = c(0, species_MM),
  ylab      = "Number of S-alleles",
  xlab      = "Element Occurrence (EO; sorted by parent BL)",
  main      = "S-allele erosion by genetic drift per EO",
  las       = 1,
  border    = "white"
)

# Thin BL color strip along the x-axis baseline (visual grouping by BL)
strip_y_top <- -species_MM * 0.018
strip_y_bot <- -species_MM * 0.045
strip_w <- diff(bp[1:2]) * 0.85
rect(bp - strip_w / 2, strip_y_bot, bp + strip_w / 2, strip_y_top,
     col = eo_df$bar_color, border = NA, xpd = TRUE)

# Species MM reference line
abline(h = species_MM, col = "grey30", lwd = 1.5, lty = 2)
mtext(paste0("Species MM: ", species_MM),
      side = 4, at = species_MM,
      col = "grey30", cex = 0.72, las = 1, line = 0.3)

# Labels inside bars
for (i in seq_along(bp)) {
  text(bp[i], eo_df$observed[i] / 2,
       labels = eo_df$observed[i],
       col = "white", font = 2, cex = 0.9)
  if (eo_df$predicted_extra[i] >= 2) {
    text(bp[i],
         eo_df$observed[i] + eo_df$predicted_extra[i] / 2,
         labels = paste0("+", eo_df$predicted_extra[i]),
         col = "white", cex = 0.78)
  }
  pct_lost <- round(eo_df$lost_to_drift[i] / species_MM * 100)
  y_red_mid <- eo_df$mm_eo[i] + eo_df$lost_to_drift[i] / 2
  text(bp[i], y_red_mid,
       labels = paste0(pct_lost, "%"),
       col = "white", font = 2, cex = 0.85)
}

# Legend placed in the top margin above the plot area
legend(
  x      = mean(par("usr")[1:2]),
  y      = par("usr")[4] + diff(par("usr")[3:4]) * 0.18,
  legend = c(
    "Observed alleles",
    "Predicted undetected (EO MM \u2212 observed)",
    "Lost to genetic drift (species MM \u2212 EO MM)"
  ),
  fill   = c(col_observed, col_predicted, col_lost),
  border = "white",
  bty    = "n",
  cex    = 0.78,
  xjust  = 0.5,
  yjust  = 1,
  xpd    = TRUE
)

dev.off()
par(mar = c(5, 4, 4, 2))
cat("PNG written: figures/SRK_allele_accumulation_drift_erosion.png\n")

# ---- 8c-ii. BL-level drift erosion ----

bl_df <- data.frame(
  BL       = valid_bls,
  observed = sapply(valid_bls, function(b) bl_results[[b]]$true_alleles),
  mm_bl    = sapply(valid_bls, function(b) {
    mm <- bl_est[[b]]$mm$Smax
    if (is.na(mm)) bl_results[[b]]$true_alleles else mm
  }),
  bar_color = unname(bl_colours[valid_bls]),
  stringsAsFactors = FALSE
)
bl_df$predicted_extra <- pmax(0, bl_df$mm_bl  - bl_df$observed)
bl_df$lost_to_drift   <- pmax(0, species_MM   - bl_df$mm_bl)

png("figures/SRK_allele_accumulation_BL_drift_erosion.png",
    width = 7, height = 6, units = "in", res = 200)
par(mar = c(5, 5, 7, 5))

bar_mat_bl <- rbind(bl_df$observed, bl_df$predicted_extra, bl_df$lost_to_drift)
bp_bl <- barplot(
  bar_mat_bl, beside = FALSE,
  col = c(col_observed, col_predicted, col_lost),
  names.arg = bl_df$BL, ylim = c(0, species_MM),
  ylab = "Number of S-alleles",
  xlab = "Bottleneck Lineage (BL)",
  main = "S-allele erosion by genetic drift per BL",
  las = 1, border = "white"
)
strip_w_bl <- diff(bp_bl[1:2]) * 0.85
rect(bp_bl - strip_w_bl / 2, strip_y_bot, bp_bl + strip_w_bl / 2, strip_y_top,
     col = bl_df$bar_color, border = NA, xpd = TRUE)
abline(h = species_MM, col = "grey30", lwd = 1.5, lty = 2)
mtext(paste0("Species MM: ", species_MM), side = 4, at = species_MM,
      col = "grey30", cex = 0.72, las = 1, line = 0.3)

for (i in seq_along(bp_bl)) {
  text(bp_bl[i], bl_df$observed[i] / 2, labels = bl_df$observed[i],
       col = "white", font = 2, cex = 0.9)
  if (bl_df$predicted_extra[i] >= 2) {
    text(bp_bl[i], bl_df$observed[i] + bl_df$predicted_extra[i] / 2,
         labels = paste0("+", bl_df$predicted_extra[i]),
         col = "white", cex = 0.78)
  }
  pct_lost <- round(bl_df$lost_to_drift[i] / species_MM * 100)
  y_red_mid <- bl_df$mm_bl[i] + bl_df$lost_to_drift[i] / 2
  text(bp_bl[i], y_red_mid, labels = paste0(pct_lost, "%"),
       col = "white", font = 2, cex = 0.85)
}

legend(
  x = mean(par("usr")[1:2]),
  y = par("usr")[4] + diff(par("usr")[3:4]) * 0.18,
  legend = c(
    "Observed alleles",
    "Predicted undetected (BL MM \u2212 observed)",
    "Lost to genetic drift (species MM \u2212 BL MM)"
  ),
  fill = c(col_observed, col_predicted, col_lost),
  border = "white", bty = "n", cex = 0.78,
  xjust = 0.5, yjust = 1, xpd = TRUE
)

dev.off()
par(mar = c(5, 4, 4, 2))
cat("PNG written: figures/SRK_allele_accumulation_BL_drift_erosion.png\n")

############################################################
# 9. Export statistics with CORRECT counts
############################################################

stats <- data.frame(
  Level="Species",
  Population="All",
  N_individuals=nrow(allele_matrix_all),
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

build_group_row <- function(level_label, group_label, n_ind, res, est) {
  data.frame(
    Level              = level_label,
    Population         = group_label,
    N_individuals      = n_ind,
    N_alleles          = res$true_alleles,
    Initial_slope      = res$initial_slope,
    Tail_slope         = res$tail_slope,
    Late_fraction      = res$late_fraction,
    MM_estimate        = est$mm$Smax,
    Chao1_estimate     = est$chao1$Chao1,
    iNEXT_estimate     = est$inext$S_asymptote,
    N_more_next_allele = est$n_more$next1,
    N_more_80pct_MM    = est$n_more$p80,
    N_more_90pct_MM    = est$n_more$p90,
    N_more_95pct_MM    = est$n_more$p95,
    N_more_99pct_MM    = est$n_more$p99,
    stringsAsFactors   = FALSE
  )
}

# EO rows (sorted by parent BL via valid_eos order)
for (eo in valid_eos) {
  n_ind <- sum(geno$EO == eo)
  stats <- rbind(stats,
                 build_group_row("EO", eo, n_ind,
                                 eo_results[[eo]], eo_est[[eo]]))
}

# BL rows
for (bl in valid_bls) {
  n_ind <- sum(geno$BL == bl)
  stats <- rbind(stats,
                 build_group_row("BL", bl, n_ind,
                                 bl_results[[bl]], bl_est[[bl]]))
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

print(stats, row.names = FALSE)

cat("\nStep 15 complete.\n")

