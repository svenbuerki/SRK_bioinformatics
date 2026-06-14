#!/usr/bin/env Rscript
# =============================================================================
# SRK_P_compat_traffic_light_with_nulls.R — NULL-AWARE variant
# Random-mating viability per EO under strict tetraploid sporophytic SI,
# computed on the null-aware genotype matrix from Step 26. Outputs are
# suffixed `_with_nulls`; the functional-only originals stay frozen.
# =============================================================================
#
# Stand-alone presentation figure built on top of the Step 17 metrics table.
# Shows the strict-SI compatible-pair fraction (P_compat at L = 0) for each
# focal EO as a horizontal bar, coloured by a traffic-light scheme:
#
#   green   P_compat >= 0.40  random mating sustainable
#   amber   0.20 <= P_compat < 0.40  random mating struggling
#   red     P_compat < 0.20  random mating effectively failed
#
# A parent-BL colour strip on the left of each bar gives the lineage
# context without competing with the traffic-light message.
#
# Inputs:
#   Tables/SRK_EO_allele_richness.tsv  (from SRK_TP1_compatibility_metrics.py)
#
# Outputs:
#   figures/SRK_P_compat_traffic_light_EO_with_nulls.png
#   figures/SRK_P_compat_traffic_light_EO_with_nulls.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source("srk_bl_constants.R")

# Thresholds
T_AMBER  <- 0.40   # green / amber boundary — random mating threshold
T_RED    <- 0.20   # amber / red boundary — effectively failed

# Traffic-light palette
LIGHT_COLS <- c(
  green = "#4CAF50",
  amber = "#FFB300",
  red   = "#E53935"
)

# =============================================================================
# 1. LOAD
# =============================================================================
cat("Loading SRK_EO_allele_richness_with_nulls.tsv ...\n")
dat <- read_tsv("Tables/Phase5/step17b_EO_allele_richness_with_nulls.tsv",
                show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(
    P_compat = P_compat_L0,
    P_lo     = P_compat_L0_lo,
    P_hi     = P_compat_L0_hi,
    light = case_when(
      P_compat >= T_AMBER ~ "green",
      P_compat >= T_RED   ~ "amber",
      TRUE                ~ "red"
    ),
    light = factor(light, levels = c("red", "amber", "green")),
    label_y = paste0("EO", group, " (", BL, ")")
  ) %>%
  arrange(desc(P_compat)) %>%
  mutate(label_y = factor(label_y, levels = label_y))

cat("\nTraffic-light tally:\n")
print(table(dat$light))

# =============================================================================
# 2. PLOT
# =============================================================================
xmax <- 1.0

# Base plot: zones, threshold lines, EO labels on y, BL strip on left, theme.
# Includes geom_blank() over the bar aesthetics so the y-axis factor levels
# (EO names sorted worst-to-best) are realised even with no bars drawn.
p_base <- ggplot(dat, aes(x = P_compat, y = label_y, fill = light)) +
  # Zone bands as light vertical stripes — at-a-glance reference
  annotate("rect", xmin = 0, xmax = T_RED,
           ymin = -Inf, ymax = Inf,
           fill = LIGHT_COLS["red"], alpha = 0.07) +
  annotate("rect", xmin = T_RED, xmax = T_AMBER,
           ymin = -Inf, ymax = Inf,
           fill = LIGHT_COLS["amber"], alpha = 0.07) +
  annotate("rect", xmin = T_AMBER, xmax = xmax,
           ymin = -Inf, ymax = Inf,
           fill = LIGHT_COLS["green"], alpha = 0.07) +
  # Threshold lines
  geom_vline(xintercept = T_AMBER, linetype = "dashed",
             colour = "grey25", linewidth = 0.6) +
  geom_vline(xintercept = T_RED,   linetype = "dotted",
             colour = "grey45", linewidth = 0.5) +
  # BL colour strip on the left
  geom_point(aes(x = -0.025, colour = BL),
             shape = 15, size = 5, show.legend = FALSE) +
  # Invisible geom_blank() to register the y-axis levels and the fill scale
  # even when no bars are drawn (the blank presentation variant)
  geom_blank(aes(x = P_compat, fill = light)) +
  scale_fill_manual(values = LIGHT_COLS, drop = FALSE,
                    name = "Random-mating status",
                    labels = c(red   = "Failed  (< 0.20)",
                               amber = "Struggling  (0.20 – 0.40)",
                               green = "Sustainable  (>= 0.40)")) +
  scale_colour_manual(values = BL_COLORS, drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.04, xmax),
    breaks = seq(0, 1, by = 0.20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Random-mating viability per Element Occurrence",
    subtitle = paste0(
      "Strict tetraploid sporophytic SI (L = 0). ",
      "Left-edge square = parent BL.\n",
      "Error bars = bootstrap 95% CI (1000 resamples).\n",
      "Dashed line = 0.40 (random-mating threshold); ",
      "dotted line = 0.20 (effectively failed)."
    ),
    x = "Compatible-pair fraction  (fraction of random pairs that produce seed)",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    plot.title         = element_text(face = "bold"),
    axis.text.y        = element_text(face = "bold", size = 12)
  )

# Full plot = base + bars + error bars + numeric labels
p_full <- p_base +
  geom_col(width = 0.65, colour = "grey15", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = P_lo, xmax = P_hi),
                 height = 0.25, colour = "grey20", linewidth = 0.6) +
  geom_text(aes(x = P_hi,
                label = sprintf("%.2f  [%.2f, %.2f]",
                                P_compat, P_lo, P_hi)),
            hjust = -0.05, fontface = "bold", size = 3.8)

# Blank plot for presentations — zones + threshold lines + EO labels + BL
# strip, but no data drawn. The legend is forced via a fully transparent
# geom_col so the audience still sees the red/amber/green key.
p_blank <- p_base +
  geom_col(width = 0.65, alpha = 0, colour = NA) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "grey15")))

# =============================================================================
# 3. SAVE
# =============================================================================
dir.create("figures/Phase5", recursive = TRUE, showWarnings = FALSE)

# Full figure
ggsave("figures/Phase5/step17b_P_compat_traffic_light_with_nulls.png", plot = p_full,
       width = 10, height = 5.5, dpi = 200)
ggsave("figures/Phase5/step17b_P_compat_traffic_light_with_nulls.pdf", plot = p_full,
       width = 10, height = 5.5)

# Blank figure (zones + threshold lines + EO labels + BL strip; no data)
ggsave("figures/Phase5/step17b_P_compat_traffic_light_with_nulls_blank.png", plot = p_blank,
       width = 10, height = 5.5, dpi = 200)
ggsave("figures/Phase5/step17b_P_compat_traffic_light_with_nulls_blank.pdf", plot = p_blank,
       width = 10, height = 5.5)

cat("\nWritten:\n")
cat("  figures/SRK_P_compat_traffic_light_EO_with_nulls.png\n")
cat("  figures/SRK_P_compat_traffic_light_EO_with_nulls.pdf\n")
cat("  figures/SRK_P_compat_traffic_light_EO_with_nulls_blank.png\n")
cat("  figures/SRK_P_compat_traffic_light_EO_with_nulls_blank.pdf\n")
