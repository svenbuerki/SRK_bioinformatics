#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_tipping_point.R - Step 17 of Canu_amplicon pipeline
# Tipping Point 1 (TP1) - Allele richness and frequency evenness status
# =============================================================================
#
# TP1 is breached when a population has lost so many S-alleles relative to the
# species optimum that inter-population allele transfers are required to
# restore richness. Two criteria are evaluated per group:
#
#   prop_optimum < 0.50   group retains < 50% of the species-level S-allele pool
#   evenness     < 0.80   Ne/N_alleles ratio below 0.80 (allele frequencies
#                         substantially skewed from NFDS equal-frequency
#                         expectation)
#
# Status logic:
#   CRITICAL -- both criteria breached
#   AT RISK  -- either criterion breached
#   OK       -- neither criterion breached
#
# Two analysis levels are tested in parallel and overlaid on the same scatter:
#   - Element Occurrence (EO)        circles, colored by parent BL
#   - Bottleneck Lineage (BL1-BL5)   triangles, colored by BL
# Both share the same Set1 BL palette as the spatial-clustering project.
# =============================================================================
# INPUTS
#   SRK_population_genetic_summary.tsv     -- from Step 14 (EO-level, BL column)
#   SRK_population_genetic_summary_BL.tsv  -- from Step 14 (BL-level aggregates)
#   SRK_species_richness_estimates.tsv     -- from Step 15
#
# OUTPUTS
#   SRK_TP1_summary.tsv                    -- per-EO TP1 status with BL column
#   SRK_TP1_summary_BL.tsv                 -- per-BL TP1 status (NEW)
#   SRK_TP1_tipping_point.pdf              -- scatter (root, for backwards compat)
#   figures/SRK_TP1_tipping_point.png      -- same as PDF
#   figures/SRK_TP1_tipping_point_blank.png -- empty zones for presentations
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) {
  message("Note: install ggrepel for non-overlapping labels.")
}

# =============================================================================
# USER SETTINGS
# =============================================================================

POPGEN_FILE      <- "SRK_population_genetic_summary.tsv"
POPGEN_BL_FILE   <- "SRK_population_genetic_summary_BL.tsv"
RICHNESS_FILE    <- "SRK_species_richness_estimates.tsv"

# Locked Set1 BL palette (matches Steps 14-16 and LEPA_EO_spatial_clustering)
BL_PALETTE <- c(
  BL1 = "#E41A1C", BL2 = "#377EB8", BL3 = "#4DAF4A",
  BL4 = "#984EA3", BL5 = "#FF7F00"
)

# Minimum sample size for an EO to be plotted (BLs are always plotted)
EO_MIN_N <- 5

# TP1 thresholds
TP1_PROP_OPTIMUM <- 0.50
TP1_EVENNESS     <- 0.80

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data (Step 17)...\n")

popgen     <- read_tsv(POPGEN_FILE,    show_col_types = FALSE)
popgen_bl  <- read_tsv(POPGEN_BL_FILE, show_col_types = FALSE)
richness   <- read_tsv(RICHNESS_FILE,  show_col_types = FALSE)

species_optimum <- richness$MM_estimate[1]
cat(sprintf("  Species optimum (MM estimate): %d alleles\n", species_optimum))

# =============================================================================
# 2. COMPUTE TP1 METRICS (EO and BL in parallel)
# =============================================================================

compute_tp1 <- function(df, label_col) {
  df %>%
    mutate(
      prop_optimum = N_alleles / species_optimum,
      evenness     = Effective_alleles_Ne / N_alleles,
      TP1_richness_breach = prop_optimum < TP1_PROP_OPTIMUM,
      TP1_evenness_breach = evenness     < TP1_EVENNESS,
      TP1_status = case_when(
        TP1_richness_breach & TP1_evenness_breach ~ "CRITICAL",
        TP1_richness_breach | TP1_evenness_breach ~ "AT RISK",
        TRUE                                       ~ "OK"
      )
    )
}

# ---- EO level ----
cat("\nComputing TP1 metrics per EO (N >=", EO_MIN_N, ")...\n")

tp1_eo <- popgen %>%
  filter(N_individuals >= EO_MIN_N) %>%
  compute_tp1("EO") %>%
  select(BL, EO, N_individuals, N_alleles, Effective_alleles_Ne,
         prop_optimum, evenness,
         TP1_richness_breach, TP1_evenness_breach, TP1_status) %>%
  arrange(BL, prop_optimum)

write_tsv(tp1_eo, "SRK_TP1_summary.tsv")
cat("  Written: SRK_TP1_summary.tsv\n")

cat("\n-- EO-level TP1 Summary --\n")
print(tp1_eo, n = Inf)

# ---- BL level ----
cat("\nComputing TP1 metrics per BL...\n")

tp1_bl <- popgen_bl %>%
  compute_tp1("BL") %>%
  select(BL, N_EOs, N_individuals, N_alleles, Effective_alleles_Ne,
         prop_optimum, evenness,
         TP1_richness_breach, TP1_evenness_breach, TP1_status) %>%
  arrange(BL)

write_tsv(tp1_bl, "SRK_TP1_summary_BL.tsv")
cat("  Written: SRK_TP1_summary_BL.tsv\n")

cat("\n-- BL-level TP1 Summary --\n")
print(tp1_bl, n = Inf)

# =============================================================================
# 3. PLOT: combined EO + BL scatter
# =============================================================================
cat("\nGenerating combined TP1 scatter plot...\n")

# Combine for plotting
tp1_eo_pl <- tp1_eo %>%
  mutate(label = EO, level = "EO")
tp1_bl_pl <- tp1_bl %>%
  mutate(label = BL, level = "BL", EO = BL) %>%   # EO field reused for fill key
  select(any_of(c("BL", "EO", "N_individuals", "N_alleles",
                  "Effective_alleles_Ne", "prop_optimum", "evenness",
                  "TP1_richness_breach", "TP1_evenness_breach",
                  "TP1_status", "label", "level")))
tp1_combined <- bind_rows(tp1_eo_pl, tp1_bl_pl) %>%
  mutate(BL = factor(BL, levels = names(BL_PALETTE)),
         level = factor(level, levels = c("EO", "BL")))

# Base plot: zone polygons + threshold lines (no points)
p_base <- ggplot(tp1_combined,
                 aes(x = prop_optimum, y = evenness)) +
  # OK zone: top-right
  annotate("rect",
           xmin = TP1_PROP_OPTIMUM, xmax = 1,
           ymin = TP1_EVENNESS,     ymax = 1,
           fill = "#4575b4", alpha = 0.08) +
  annotate("text",
           x = (TP1_PROP_OPTIMUM + 1) / 2,
           y = (TP1_EVENNESS + 1) / 2,
           label = "OK",
           colour = "#4575b4", fontface = "bold", size = 4.5) +
  # AT RISK zone A: top-left (richness breach only)
  annotate("rect",
           xmin = 0,                xmax = TP1_PROP_OPTIMUM,
           ymin = TP1_EVENNESS,     ymax = 1,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = TP1_PROP_OPTIMUM / 2,
           y = (TP1_EVENNESS + 1) / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # AT RISK zone B: bottom-right (evenness breach only)
  annotate("rect",
           xmin = TP1_PROP_OPTIMUM, xmax = 1,
           ymin = 0,                ymax = TP1_EVENNESS,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = (TP1_PROP_OPTIMUM + 1) / 2,
           y = TP1_EVENNESS / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # CRITICAL zone: bottom-left (label nudged toward the corner so it doesn't
  # overlap with the cluster of CRITICAL EO/BL points typically lying near
  # the centre of this zone)
  annotate("rect",
           xmin = 0,                xmax = TP1_PROP_OPTIMUM,
           ymin = 0,                ymax = TP1_EVENNESS,
           fill = "#d73027", alpha = 0.08) +
  annotate("text",
           x = TP1_PROP_OPTIMUM * 0.18,
           y = TP1_EVENNESS    * 0.22,
           label = "CRITICAL",
           colour = "#d73027", fontface = "bold", size = 4) +
  # threshold lines
  geom_vline(xintercept = TP1_PROP_OPTIMUM, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = TP1_EVENNESS, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.05)
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01),
    limits = c(0, 1),
    expand = expansion(mult = 0.05)
  ) +
  labs(
    title    = "Tipping Point 1 - Allele richness and frequency evenness",
    subtitle = paste0(
      "Circles = EO (N >= ", EO_MIN_N, "), triangles = BL aggregates. ",
      "Colors = parent BL (Set1 palette, matching LEPA_EO_spatial_clustering)."
    ),
    x = paste0("Proportion of species optimum retained  (species optimum = ",
               species_optimum, " alleles)"),
    y = "Frequency evenness  (Ne / N alleles)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.box = "vertical")

ggsave("figures/SRK_TP1_tipping_point_blank.png", plot = p_base,
       width = 10, height = 8, dpi = 200)
cat("  Written: figures/SRK_TP1_tipping_point_blank.png\n")

# Full plot: EO circles + BL triangles, both colored by BL
p <- p_base +
  geom_point(data = tp1_combined,
             aes(fill = BL, shape = level),
             size = 5.5, stroke = 0.9, colour = "grey20", alpha = 0.95) +
  scale_fill_manual(values = BL_PALETTE, name = "Bottleneck lineage", drop = FALSE) +
  scale_shape_manual(values = c(EO = 21, BL = 24), name = "Level") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4.5)),
         shape = guide_legend(override.aes = list(size = 4.5,
                                                  fill = "grey80")))

# Labels
if (use_repel) {
  p <- p + ggrepel::geom_text_repel(
    data = tp1_combined,
    aes(label = label, fontface = ifelse(level == "BL", "bold.italic", "bold")),
    size = 4.2, show.legend = FALSE,
    box.padding = 0.55, point.padding = 0.4,
    seed = 42, max.overlaps = 20
  )
} else {
  p <- p + geom_text(
    data = tp1_combined,
    aes(label = label),
    vjust = -1.2, size = 4, fontface = "bold", show.legend = FALSE
  )
}

# Save the combined plot
ggsave("figures/SRK_TP1_tipping_point.png", plot = p,
       width = 11, height = 8, dpi = 200)
cat("  Written: figures/SRK_TP1_tipping_point.png\n")

ggsave("SRK_TP1_tipping_point.pdf", plot = p, width = 11, height = 8)
cat("  Written: SRK_TP1_tipping_point.pdf\n")

# =============================================================================
# 4. CONSOLE SUMMARY
# =============================================================================
cat("\n-- Groups breaching TP1 thresholds --\n")
cat("\nEO level:\n")
tp1_eo %>%
  filter(TP1_status != "OK") %>%
  select(BL, EO, N_alleles, prop_optimum, evenness, TP1_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nBL level:\n")
tp1_bl %>%
  filter(TP1_status != "OK") %>%
  select(BL, N_alleles, prop_optimum, evenness, TP1_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nStep 17 complete.\n")
