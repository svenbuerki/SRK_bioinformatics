############################################################
# SRK FUNCTIONAL-ALLELE accumulation analysis — parallel to Step 15
#
# Uses the collapsed individual × functional-group matrix produced by
# srk_functional_definitions.py, with alleles grouped by cross-genera
# HV union + polymorphic linker positions (see step26a_functional_
# site_positions.tsv).
#
# Produces species-level and per-BL accumulation curves + MM/Chao1
# estimators, side-by-side with the existing sequence-allele curves
# (Step 15). The two side-by-side reports are the input to the paper /
# NSF proposal comparison.
#
# Outputs:
#   Tables/Phase3/step15_functional_allele_accumulation_stats.tsv
#   Tables/Phase3/step15_functional_vs_sequence_comparison.tsv
#   figures/Phase3/step15_functional_allele_accumulation_species.pdf|.png
#   figures/Phase3/step15_functional_allele_accumulation_BL_combined.pdf|.png
#   figures/Phase3/step15_functional_vs_sequence_comparison.pdf|.png
############################################################

pdf(NULL)
cat("Starting SRK FUNCTIONAL-allele accumulation analysis\n")

suppressMessages({
  source("srk_bl_constants.R")
})
BL_PALETTE <- BL_COLORS

# ----------------------------------------------------------
# Load data
# ----------------------------------------------------------
geno_func <- read.table(
  "Tables/Phase2/step11_functional_allele_genotypes.tsv",
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
geno_seq <- read.table(
  "Tables/Phase2/step11_individual_allele_genotypes.tsv",
  header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
)
sample_info <- read.csv("Tables/sampling_metadata.csv", stringsAsFactors = FALSE)
bl_assign <- read.table(
  "Tables/Phase3/step13_individual_BL_assignments.tsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
cat("Loaded:\n")
cat("  functional matrix:", nrow(geno_func), "ind ×", ncol(geno_func) - 1, "groups\n")
cat("  sequence   matrix:", nrow(geno_seq),  "ind ×", ncol(geno_seq)  - 1, "alleles\n")
cat("  BL assignments  :", nrow(bl_assign), "individuals\n")

# ----------------------------------------------------------
# Filter ingroup & attach EO/BL
# ----------------------------------------------------------
sample_info <- sample_info[sample_info$Ingroup == 1, ]
keep_ingroup <- function(df) df[df$Individual %in% sample_info$SampleID, ]
geno_func_all <- keep_ingroup(geno_func)
geno_seq_all  <- keep_ingroup(geno_seq)

bl_ok <- bl_assign[bl_assign$BL_status %in% c("Assigned", "Inferred"), ]
attach_bl <- function(df) {
  m <- df[df$Individual %in% bl_ok$Individual, ]
  m$EO <- bl_ok$EO[match(m$Individual, bl_ok$Individual)]
  m$BL <- bl_ok$BL[match(m$Individual, bl_ok$Individual)]
  m
}
geno_func <- attach_bl(geno_func_all)
geno_seq  <- attach_bl(geno_seq_all)
cat("Species-level (all ingroup): ", nrow(geno_func_all), "individuals\n")
cat("BL-assigned subset          :", nrow(geno_func),      "individuals\n\n")

# ----------------------------------------------------------
# Accumulation-curve helpers
# ----------------------------------------------------------
build_presence <- function(df, count_cols) {
  m <- as.matrix(df[, count_cols, drop = FALSE])
  storage.mode(m) <- "numeric"
  m[is.na(m)] <- 0
  rownames(m) <- df$Individual
  (m > 0) * 1L
}

# Chao1 estimator on incidence data
chao1_est <- function(pres_mat) {
  freq <- colSums(pres_mat)
  freq <- freq[freq > 0]
  S_obs <- length(freq)
  f1 <- sum(freq == 1); f2 <- sum(freq == 2)
  if (f2 > 0) S_obs + (f1^2) / (2 * f2)
  else        S_obs + (f1 * (f1 - 1)) / 2
}

# Michaelis-Menten estimator via NLS on rarefaction curve
mm_est <- function(pres_mat, n_boot = 1000, seed = 20260818) {
  set.seed(seed)
  n_ind <- nrow(pres_mat)
  if (n_ind < 5) return(NA_real_)
  xs <- unique(round(seq(1, n_ind, length.out = min(30, n_ind))))
  ys <- sapply(xs, function(x) {
    mean(sapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n_ind, x, replace = FALSE)
      sum(colSums(pres_mat[idx, , drop = FALSE]) > 0)
    }))
  })
  fit <- tryCatch(nls(ys ~ Vmax * xs / (K + xs),
                      start = list(Vmax = max(ys) * 1.5, K = median(xs)),
                      control = nls.control(warnOnly = TRUE, maxiter = 200)),
                  error = function(e) NULL)
  if (is.null(fit)) NA_real_ else as.numeric(coef(fit)["Vmax"])
}

# Rarefied accumulation curve: mean richness at each sample size 1..n
accumulation_curve <- function(pres_mat, n_boot = 1000, seed = 20260818) {
  set.seed(seed)
  n_ind <- nrow(pres_mat)
  if (n_ind < 2) return(data.frame(N = integer(0), S_mean = numeric(0)))
  xs <- 1:n_ind
  ys <- sapply(xs, function(x) {
    mean(sapply(seq_len(n_boot), function(b) {
      idx <- sample.int(n_ind, x, replace = FALSE)
      sum(colSums(pres_mat[idx, , drop = FALSE]) > 0)
    }))
  })
  data.frame(N = xs, S_mean = ys)
}

summarize_level <- function(pres_mat, label) {
  data.frame(
    Group     = label,
    N_ind     = nrow(pres_mat),
    S_obs     = sum(colSums(pres_mat) > 0),
    Chao1     = round(chao1_est(pres_mat), 2),
    MM_Vmax   = round(mm_est(pres_mat), 2),
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# Compute stats: species + per-BL, for both functional and sequence
# ----------------------------------------------------------
func_cols <- setdiff(colnames(geno_func_all), c("Individual", "EO", "BL"))
seq_cols  <- setdiff(colnames(geno_seq_all),  c("Individual", "EO", "BL"))

pres_func_species <- build_presence(geno_func_all, func_cols)
pres_seq_species  <- build_presence(geno_seq_all,  seq_cols)

stats_rows <- list(
  cbind(Level = "Species", Definition = "Functional",
        summarize_level(pres_func_species, "L. papilliferum")),
  cbind(Level = "Species", Definition = "Sequence",
        summarize_level(pres_seq_species,  "L. papilliferum"))
)

for (bl in unique(sort(geno_func$BL))) {
  sub_f <- geno_func[geno_func$BL == bl, ]
  sub_s <- geno_seq[geno_seq$BL == bl, ]
  if (nrow(sub_f) < 2) next
  stats_rows[[length(stats_rows) + 1]] <- cbind(Level = "BL", Definition = "Functional",
    summarize_level(build_presence(sub_f, func_cols), bl))
  stats_rows[[length(stats_rows) + 1]] <- cbind(Level = "BL", Definition = "Sequence",
    summarize_level(build_presence(sub_s, seq_cols),  bl))
}
stats <- do.call(rbind, stats_rows)
rownames(stats) <- NULL

# Also add side-by-side comparison (Functional vs Sequence per group)
compare <- reshape(stats, direction = "wide",
                   idvar = c("Level", "Group"),
                   timevar = "Definition",
                   v.names = c("N_ind", "S_obs", "Chao1", "MM_Vmax"))
compare$S_obs_ratio <- round(compare$S_obs.Functional / compare$S_obs.Sequence, 3)

dir.create("Tables/Phase3", showWarnings = FALSE, recursive = TRUE)
write.table(stats,   "Tables/Phase3/step15_functional_allele_accumulation_stats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(compare, "Tables/Phase3/step15_functional_vs_sequence_comparison.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("→ wrote Tables/Phase3/step15_functional_allele_accumulation_stats.tsv\n")
cat("→ wrote Tables/Phase3/step15_functional_vs_sequence_comparison.tsv\n\n")
print(stats)
cat("\n")

# ----------------------------------------------------------
# Species-level curves (functional overlay on sequence)
# ----------------------------------------------------------
dir.create("figures/Phase3", showWarnings = FALSE, recursive = TRUE)

curve_f <- accumulation_curve(pres_func_species)
curve_s <- accumulation_curve(pres_seq_species)

save_species_plot <- function(out_path, is_pdf = TRUE) {
  if (is_pdf) pdf(out_path, width = 8, height = 5.5)
  else        png(out_path, width = 8 * 200, height = 5.5 * 200, res = 200)
  par(mar = c(4.5, 4.5, 3, 2))
  ymax <- max(curve_s$S_mean, curve_f$S_mean) * 1.05
  plot(curve_s$N, curve_s$S_mean, type = "l", lwd = 2.5, col = "#666666",
       xlab = "Individuals sampled", ylab = "Distinct alleles",
       main = sprintf("SRK allele accumulation (species-wide, n = %d)", nrow(pres_func_species)),
       ylim = c(0, ymax), las = 1)
  lines(curve_f$N, curve_f$S_mean, lwd = 2.5, col = "#b2182b")
  abline(h = stats$S_obs[stats$Level == "Species" & stats$Definition == "Sequence"],
         lty = 3, col = "#666666")
  abline(h = stats$S_obs[stats$Level == "Species" & stats$Definition == "Functional"],
         lty = 3, col = "#b2182b")
  legend("bottomright",
         legend = c(
           sprintf("Sequence alleles   (S_obs = %d, MM = %s, Chao1 = %s)",
                   stats$S_obs[stats$Level == "Species" & stats$Definition == "Sequence"],
                   stats$MM_Vmax[stats$Level == "Species" & stats$Definition == "Sequence"],
                   stats$Chao1[stats$Level == "Species" & stats$Definition == "Sequence"]),
           sprintf("Functional groups  (S_obs = %d, MM = %s, Chao1 = %s)",
                   stats$S_obs[stats$Level == "Species" & stats$Definition == "Functional"],
                   stats$MM_Vmax[stats$Level == "Species" & stats$Definition == "Functional"],
                   stats$Chao1[stats$Level == "Species" & stats$Definition == "Functional"])
         ),
         col = c("#666666", "#b2182b"), lwd = 2.5, bty = "n", cex = 0.9)
  dev.off()
}
save_species_plot("figures/Phase3/step15_functional_allele_accumulation_species.pdf", TRUE)
save_species_plot("figures/Phase3/step15_functional_allele_accumulation_species.png", FALSE)
cat("→ wrote species-level accumulation figure\n")

# ----------------------------------------------------------
# Per-BL curves (functional definition), coloured by BL
# ----------------------------------------------------------
bl_levels <- if (exists("BL_ORDER")) BL_ORDER else sort(unique(geno_func$BL))
bl_levels <- intersect(bl_levels, unique(geno_func$BL))

save_bl_plot <- function(out_path, is_pdf = TRUE) {
  # Compute per-BL curves for BOTH definitions
  bl_curves <- function(matrix_df, cols) {
    out <- lapply(bl_levels, function(bl) {
      sub <- matrix_df[matrix_df$BL == bl, ]
      if (nrow(sub) < 2) return(NULL)
      curve <- accumulation_curve(build_presence(sub, cols))
      curve$BL <- bl
      curve
    })
    Filter(Negate(is.null), out)
  }
  curves_seq  <- bl_curves(geno_seq,  seq_cols)
  curves_func <- bl_curves(geno_func, func_cols)

  # Shared axes so the visual drop from sequence to functional is direct
  xmax <- max(c(sapply(curves_seq, function(d) max(d$N)),
                sapply(curves_func, function(d) max(d$N))))
  ymax <- max(c(sapply(curves_seq, function(d) max(d$S_mean)),
                sapply(curves_func, function(d) max(d$S_mean)))) * 1.05

  # Two-panel side-by-side layout, sequence left / functional right
  if (is_pdf) pdf(out_path, width = 12, height = 5.5)
  else        png(out_path, width = 12 * 200, height = 5.5 * 200, res = 200)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5), oma = c(0, 0, 1.5, 0))

  draw_panel <- function(curves, panel_title, y_label, stats_def) {
    plot(NA, xlim = c(1, xmax), ylim = c(0, ymax), las = 1,
         xlab = "Individuals sampled", ylab = y_label,
         main = panel_title)
    for (d in curves) {
      col_bl <- BL_PALETTE[d$BL[1]]
      if (is.na(col_bl)) col_bl <- "#999999"
      lines(d$N, d$S_mean, lwd = 2.4, col = col_bl)
    }
    bl_labels <- sapply(curves, function(d) {
      row <- stats[stats$Level == "BL" & stats$Definition == stats_def &
                   stats$Group == d$BL[1], ]
      sprintf("%s (n = %d, S = %d)", d$BL[1], row$N_ind, row$S_obs)
    })
    legend("bottomright", legend = bl_labels,
           col = sapply(curves, function(d) {
             v <- BL_PALETTE[d$BL[1]]; if (is.na(v)) "#999999" else v
           }),
           lwd = 2.4, bty = "n", cex = 0.85)
  }

  draw_panel(curves_seq,  "Sequence alleles",   "Distinct sequence alleles",  "Sequence")
  draw_panel(curves_func, "Functional groups",  "Distinct functional groups", "Functional")
  mtext("SRK allele accumulation per bottleneck lineage  —  sequence (left) vs functional (right)",
        outer = TRUE, cex = 1.05, font = 2)
  dev.off()
}
save_bl_plot("figures/Phase3/step15_functional_allele_accumulation_BL_combined.pdf", TRUE)
save_bl_plot("figures/Phase3/step15_functional_allele_accumulation_BL_combined.png", FALSE)
cat("→ wrote per-BL accumulation figure\n")

# ----------------------------------------------------------
# Side-by-side comparison bar plot
# ----------------------------------------------------------
save_compare_plot <- function(out_path, is_pdf = TRUE) {
  if (is_pdf) pdf(out_path, width = 8, height = 5.5)
  else        png(out_path, width = 8 * 200, height = 5.5 * 200, res = 200)
  par(mar = c(5, 4.5, 3, 2))
  m <- rbind(Sequence   = compare$S_obs.Sequence,
             Functional = compare$S_obs.Functional)
  colnames(m) <- compare$Group
  cols <- c("#666666", "#b2182b")
  bp <- barplot(m, beside = TRUE, col = cols, border = NA, las = 1,
                ylim = c(0, max(m) * 1.15),
                ylab = "Distinct alleles (S_obs)",
                main = "Sequence-allele vs functional-allele richness (side-by-side)")
  text(x = bp, y = m, labels = m, pos = 3, cex = 0.8)
  legend("topright", legend = c("Sequence alleles", "Functional groups"),
         fill = cols, border = NA, bty = "n")
  dev.off()
}
save_compare_plot("figures/Phase3/step15_functional_vs_sequence_comparison.pdf", TRUE)
save_compare_plot("figures/Phase3/step15_functional_vs_sequence_comparison.png", FALSE)
cat("→ wrote side-by-side comparison figure\n")

cat("\nDone. Side-by-side outputs written for direct sequence-vs-functional comparison.\n")
