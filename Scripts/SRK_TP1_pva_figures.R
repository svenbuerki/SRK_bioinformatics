#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_pva_figures.R — figures for Step 19b PVA outputs
# =============================================================================
#
# Renders summary figures for the v1 baseline PVA. Inputs are produced by
# SRK_TP1_pva.py. v1 reports 3 checkpoints (year 10, 20, 50); figures are
# step-style and add a y0 anchor (initial state from the metrics TSV).
#
# Outputs:
#   figures/SRK_TP1_pva_trajectory.{png,pdf}    — P_compat + N + L over time
#   figures/SRK_TP1_pva_collapse_summary.{png,pdf} — P(collapsed) per EO
#   figures/SRK_TP1_pva_allele_heatmap.{png,pdf}   — per-allele loss heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("srk_bl_constants.R")

CHECKPOINTS <- c(10, 20, 50)
P_FLOOR <- 0.05
P_VIABLE <- 0.40

# ----- Read inputs ----------------------------------------------------------
metrics <- read_tsv("Tables/SRK_EO_allele_richness.tsv",
                    show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(EO = paste0("EO", group)) %>%
  select(EO, BL, P_compat_0 = P_compat_L0,
         k_obs_0 = k_observed, N_0 = N, prop_AAAA_0 = prop_AAAA)

# Read from v1.6 prediction outputs (preferred — empirical per-EO σ_env).
# Falls back to v1 baseline outputs if prediction not yet run.
PREDICTION_TSV <- "Tables/SRK_TP1_pva_prediction_summary.tsv"
PREDICTION_ALLELES <- "Tables/SRK_TP1_pva_prediction_alleles.tsv"
BASELINE_TSV <- "Tables/SRK_TP1_pva_checkpoint_summary.tsv"
BASELINE_ALLELES <- "Tables/SRK_TP1_pva_allele_extinction.tsv"
if (file.exists(PREDICTION_TSV)) {
  cat("Using v1.6 prediction outputs.\n")
  summary <- read_tsv(PREDICTION_TSV, show_col_types = FALSE)
  alleles <- read_tsv(PREDICTION_ALLELES, show_col_types = FALSE)
} else {
  cat("Using v1 baseline outputs (run SRK_TP1_pva_prediction.py for empirical-parameter version).\n")
  summary <- read_tsv(BASELINE_TSV, show_col_types = FALSE)
  alleles <- read_tsv(BASELINE_ALLELES, show_col_types = FALSE)
}

# Add baseline (year 0) row per EO from the metrics TSV.
baseline <- metrics %>%
  transmute(EO,
            checkpoint_year = 0L,
            P_collapsed = 0,
            N_adult_p10 = N_0, N_adult_p50 = N_0, N_adult_p90 = N_0,
            k_alleles_p10 = k_obs_0, k_alleles_p50 = k_obs_0,
            k_alleles_p90 = k_obs_0,
            P_compat_p10 = P_compat_0, P_compat_p50 = P_compat_0,
            P_compat_p90 = P_compat_0,
            L_p10 = prop_AAAA_0 / 3.5, L_p50 = prop_AAAA_0 / 3.5,
            L_p90 = prop_AAAA_0 / 3.5)

summary_num <- summary %>%
  mutate(across(-c(EO, checkpoint_year), as.numeric))

traj <- bind_rows(baseline, summary_num) %>%
  arrange(EO, checkpoint_year) %>%
  left_join(metrics %>% select(EO, BL), by = "EO")

# Order EOs worst-to-best by initial P_compat (most-vulnerable on top)
eo_order <- metrics %>%
  arrange(P_compat_0) %>%
  pull(EO)
traj$EO <- factor(traj$EO, levels = eo_order)

# =============================================================================
# Figure A — multi-metric trajectory
# =============================================================================
build_trajectory_panel <- function(metric, ymin, ymax, ylab, hlines = NULL) {
  p10 <- paste0(metric, "_p10")
  p50 <- paste0(metric, "_p50")
  p90 <- paste0(metric, "_p90")
  dat <- traj %>%
    transmute(EO, checkpoint_year, BL,
              p10 = .data[[p10]], p50 = .data[[p50]], p90 = .data[[p90]])
  g <- ggplot(dat, aes(x = checkpoint_year, y = p50, colour = BL, fill = BL)) +
    geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4, shape = 21, colour = "grey20", stroke = 0.3) +
    facet_wrap(~ EO, nrow = 2) +
    scale_colour_manual(values = BL_COLORS, name = "Bottleneck lineage",
                        breaks = BL_ORDER_NUMERIC) +
    scale_fill_manual(values = BL_COLORS, guide = "none",
                       breaks = BL_ORDER_NUMERIC) +
    scale_x_continuous(breaks = c(0, CHECKPOINTS),
                       expand = expansion(mult = 0.02)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(x = "Years from now", y = ylab) +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey95", colour = "grey60"),
          legend.position = "right")
  if (!is.null(hlines)) {
    for (hl in hlines) {
      g <- g + geom_hline(yintercept = hl$y, linetype = hl$lt,
                          colour = hl$col, linewidth = 0.4)
    }
  }
  g
}

version_tag <- if (file.exists(PREDICTION_TSV)) "v1.6 prediction" else "v1 baseline"
p_compat_panel <- build_trajectory_panel(
  "P_compat", 0, 1, "Compatible-pair fraction (P_compat)",
  hlines = list(list(y = P_FLOOR, lt = "solid", col = "#d73027"),
                list(y = P_VIABLE, lt = "dashed", col = "grey25"))
) + labs(
  title = paste0("PVA ", version_tag, " — compatible-pair fraction over 50 years"),
  subtitle = paste0(
    "Median + 80 % envelope across 1000 Monte Carlo replicates per EO. ",
    "Red line = collapse threshold (P_compat = ", P_FLOOR, "); ",
    "dashed line = random-mating viability (0.40)."))

n_panel <- build_trajectory_panel(
  "N_adult", 0, NA_real_, "Adult census N"
) + labs(
  title = paste0("PVA ", version_tag, " — adult census trajectory over 50 years"),
  subtitle = "Median + 80 % envelope per EO; K is anchored on observed census × 1.5 with empirical per-EO σ.")

L_panel <- build_trajectory_panel(
  "L", 0, NA_real_, "SI leakage rate L (prop_AAAA / 3.5)"
) + labs(
  title = paste0("PVA ", version_tag, " — SI leakage L over 50 years"),
  subtitle = "Median + 80 % envelope per EO. Initial L reflects historical selfing; L decays as continued outcrossing washes out AAAA homozygotes.")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/SRK_TP1_pva_trajectory_P_compat.png", p_compat_panel,
       width = 12, height = 7, dpi = 200)
ggsave("figures/SRK_TP1_pva_trajectory_P_compat.pdf", p_compat_panel,
       width = 12, height = 7)
ggsave("figures/SRK_TP1_pva_trajectory_N.png", n_panel,
       width = 12, height = 7, dpi = 200)
ggsave("figures/SRK_TP1_pva_trajectory_N.pdf", n_panel,
       width = 12, height = 7)
ggsave("figures/SRK_TP1_pva_trajectory_L.png", L_panel,
       width = 12, height = 7, dpi = 200)
ggsave("figures/SRK_TP1_pva_trajectory_L.pdf", L_panel,
       width = 12, height = 7)
cat("Wrote trajectory figures (P_compat, N, L)\n")

# =============================================================================
# Figure B — collapse-probability summary
# =============================================================================
collapse_df <- summary_num %>%
  select(EO, checkpoint_year, P_collapsed) %>%
  left_join(metrics %>% select(EO, BL), by = "EO") %>%
  mutate(EO = factor(EO, levels = eo_order),
         checkpoint_year = factor(checkpoint_year, levels = CHECKPOINTS,
                                   labels = paste0("y", CHECKPOINTS)))

collapse_plot <- ggplot(collapse_df,
                         aes(x = P_collapsed, y = EO, fill = BL)) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.7,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * P_collapsed)),
            position = position_dodge2(width = 0.7, reverse = TRUE),
            hjust = -0.15, size = 3.2, colour = "grey15") +
  facet_grid(. ~ checkpoint_year) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(),
                     expand = expansion(mult = c(0.005, 0.18))) +
  scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage",
                    breaks = BL_ORDER_NUMERIC) +
  labs(title = paste0("PVA ", version_tag, " — probability of collapse (P_compat < 0.05)"),
       subtitle = paste0(
         "1000 Monte Carlo replicates per EO. Collapse = compatible-pair ",
         "fraction below ", P_FLOOR, " at the checkpoint year. ",
         "EOs sorted worst-to-best by current P_compat. ",
         "NOTE: see SRK_TP1_pva_v17_lineage_fragility for the maternal-dormancy refinement."),
       x = "Probability of collapse at the checkpoint", y = NULL) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey60"),
        legend.position = "right",
        panel.grid.major.y = element_blank())

ggsave("figures/SRK_TP1_pva_collapse_summary.png", collapse_plot,
       width = 11, height = 6, dpi = 200)
ggsave("figures/SRK_TP1_pva_collapse_summary.pdf", collapse_plot,
       width = 11, height = 6)
cat("Wrote collapse-summary figure\n")

# =============================================================================
# Figure C — per-allele loss heatmap
# =============================================================================
allele_df <- alleles %>%
  mutate(P_lost_by_y50 = as.numeric(P_lost_by_y50),
         EO = factor(EO, levels = eo_order))

# Allele order: sort by ID
allele_order <- sort(unique(allele_df$allele))
allele_df$allele <- factor(allele_df$allele, levels = allele_order)

heatmap_plot <- ggplot(allele_df, aes(x = allele, y = EO, fill = P_lost_by_y50)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  geom_text(aes(label = ifelse(P_lost_by_y50 >= 0.01,
                                sprintf("%.0f%%", 100 * P_lost_by_y50),
                                "")),
            size = 2.6, colour = "grey15") +
  scale_fill_gradient2(low = "#ffffff", mid = "#fee08b", high = "#d73027",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "P(lost by y50)",
                       labels = scales::percent_format()) +
  labs(title = paste0("PVA ", version_tag, " — per-allele loss probability by year 50"),
       subtitle = paste0(
         "Each cell = P(this allele is absent from EO's adults AND all seed cohorts ",
         "at year 50), across 1000 Monte Carlo replicates. ",
         "Empty cells = allele not present in EO initially."),
       x = "S-allele", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        panel.grid = element_blank())

ggsave("figures/SRK_TP1_pva_allele_heatmap.png", heatmap_plot,
       width = 14, height = 5, dpi = 200)
ggsave("figures/SRK_TP1_pva_allele_heatmap.pdf", heatmap_plot,
       width = 14, height = 5)
cat("Wrote allele-loss heatmap\n")

cat("\nAll PVA figures written to figures/SRK_TP1_pva_*.{png,pdf}\n")
