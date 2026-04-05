#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_tipping_point.R
# Tipping Point 1 (TP1) — Allele richness and frequency evenness status
# =============================================================================
#
# TP1 is breached when a population has lost so many S-alleles relative to the
# species optimum that inter-population allele transfers are required to restore
# richness.  Two criteria are evaluated per EO:
#
#   prop_optimum < 0.50   EO retains < 50% of the species-level S-allele pool
#   evenness     < 0.80   Ne/N_alleles ratio below 0.80 (allele frequencies
#                         substantially skewed from NFDS equal-frequency
#                         expectation)
#
# Status logic (mirrors TP2):
#   CRITICAL -- both criteria breached
#   AT RISK  -- either criterion breached
#   OK       -- neither criterion breached
#
# =============================================================================
# INPUTS
#   SRK_population_genetic_summary.tsv  -- from Step 13
#   SRK_species_richness_estimates.tsv  -- from Step 14
#
# OUTPUTS
#   SRK_TP1_summary.tsv          -- per-EO richness, evenness, and TP1 status
#   SRK_TP1_tipping_point.pdf    -- scatter plot (prop_optimum x evenness)
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) {
  message("Note: install ggrepel for labelled scatter plot.")
}

# =============================================================================
# USER SETTINGS
# =============================================================================

POPGEN_FILE   <- "SRK_population_genetic_summary.tsv"
RICHNESS_FILE <- "SRK_species_richness_estimates.tsv"

# EOs to include (must match Population values in POPGEN_FILE)
EO_MAP <- c(
  "25" = "EO25",
  "27" = "EO27",
  "67" = "EO67",
  "70" = "EO70",
  "76" = "EO76"
)

EO_FOCUS <- c("EO25", "EO27", "EO67", "EO70", "EO76")

# TP1 thresholds
TP1_PROP_OPTIMUM <- 0.50   # minimum proportion of species optimum retained
TP1_EVENNESS     <- 0.80   # minimum Ne / N_alleles evenness ratio

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data...\n")

popgen   <- read_tsv(POPGEN_FILE,   show_col_types = FALSE)
richness <- read_tsv(RICHNESS_FILE, show_col_types = FALSE)

species_optimum <- richness$MM_estimate[1]
cat(sprintf("  Species optimum (MM estimate): %d alleles\n", species_optimum))

# =============================================================================
# 2. COMPUTE TP1 METRICS
# =============================================================================
cat("Computing TP1 metrics per EO...\n")

tp1 <- popgen %>%
  mutate(Population = as.character(Population)) %>%
  filter(Population %in% names(EO_MAP)) %>%
  mutate(
    EO           = EO_MAP[Population],
    prop_optimum = N_alleles / species_optimum,
    evenness     = Effective_alleles_Ne / N_alleles
  ) %>%
  mutate(
    TP1_richness_breach = prop_optimum < TP1_PROP_OPTIMUM,
    TP1_evenness_breach = evenness     < TP1_EVENNESS,
    TP1_status = case_when(
      TP1_richness_breach & TP1_evenness_breach ~ "CRITICAL",
      TP1_richness_breach | TP1_evenness_breach ~ "AT RISK",
      TRUE                                       ~ "OK"
    )
  ) %>%
  select(EO, N_individuals, N_alleles, Effective_alleles_Ne,
         prop_optimum, evenness,
         TP1_richness_breach, TP1_evenness_breach, TP1_status) %>%
  arrange(prop_optimum)

write_tsv(tp1, "SRK_TP1_summary.tsv")
cat("  Written: SRK_TP1_summary.tsv\n")

cat("\n-- EO-level TP1 Summary --\n")
print(tp1, n = Inf)

# =============================================================================
# 3. PLOT
# =============================================================================
cat("\nGenerating TP1 tipping point plot...\n")

tp1_plot <- tp1 %>%
  mutate(EO = factor(EO, levels = EO_FOCUS))

pdf("SRK_TP1_tipping_point.pdf", width = 10, height = 8)

p <- ggplot(tp1_plot,
            aes(x = prop_optimum, y = evenness,
                colour = TP1_status, label = EO)) +
  geom_point(size = 5.5, alpha = 0.9) +
  geom_vline(xintercept = TP1_PROP_OPTIMUM, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = TP1_EVENNESS, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_colour_manual(
    values = c("CRITICAL" = "#d73027", "AT RISK" = "#fc8d59", "OK" = "#4575b4"),
    name = "TP1 status"
  ) +
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
  annotate("rect",
           xmin = 0,                  xmax = TP1_PROP_OPTIMUM,
           ymin = 0,                  ymax = TP1_EVENNESS,
           fill = "#d73027", alpha = 0.08) +
  annotate("text", x = 0.02, y = 0.03, label = "CRITICAL",
           colour = "#d73027", fontface = "bold", hjust = 0, size = 4) +
  labs(
    title    = "Tipping Point 1 -- EO allele richness and frequency evenness status",
    subtitle = paste0(
      "X axis: proportion of species optimum retained (TP1 threshold = ",
      percent(TP1_PROP_OPTIMUM, accuracy = 1), "). ",
      "Y axis: frequency evenness Ne/N (TP1 threshold = ", TP1_EVENNESS, ")."
    ),
    x = paste0("Proportion of species optimum retained  (species optimum = ",
               species_optimum, " alleles)"),
    y = "Frequency evenness  (Ne / N alleles)"
  ) +
  theme_bw(base_size = 13)

if (use_repel) {
  p <- p + ggrepel::geom_text_repel(
    size = 4.5, fontface = "bold", show.legend = FALSE,
    box.padding = 0.5, seed = 42
  )
} else {
  p <- p + geom_text(vjust = -1.2, size = 4, fontface = "bold",
                     show.legend = FALSE)
}

print(p)
dev.off()
cat("  Written: SRK_TP1_tipping_point.pdf\n")

# =============================================================================
# 4. CONSOLE SUMMARY
# =============================================================================
cat("\n-- EOs breaching TP1 thresholds --\n")
tp1 %>%
  filter(TP1_status != "OK") %>%
  select(EO, N_alleles, prop_optimum, evenness, TP1_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nDone.\n")
