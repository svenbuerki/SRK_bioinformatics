#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_pva_sensitivity_figures.R — v1.5 sensitivity sweep figures
# =============================================================================
#
# Reads v1.5 PVA sensitivity outputs and produces:
#   1. Sensitivity tornado per EO — Δ P(collapse at y50) per scenario vs baseline
#   2. Per-scenario P_compat comparison panel (small multiples)
#   3. Per-EO summary table heatmap
#
# Outputs:
#   figures/SRK_TP1_pva_sensitivity_tornado.{png,pdf}
#   figures/SRK_TP1_pva_sensitivity_p_compat.{png,pdf}
#   figures/SRK_TP1_pva_sensitivity_grid.{png,pdf}
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

source("srk_bl_constants.R")

# ----- Read inputs ----------------------------------------------------------
summ <- read_tsv("Tables/SRK_TP1_pva_sensitivity_summary.tsv",
                  show_col_types = FALSE) %>%
  mutate(across(-c(EO, scenario), as.numeric))
metrics <- read_tsv("Tables/SRK_EO_allele_richness.tsv",
                     show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(EO = paste0("EO", group)) %>%
  select(EO, BL)

# Scenario ordering (worst-case on bottom, baseline as reference)
scenario_levels <- c("baseline", "si_proper", "sc_inbreeding", "catastrophe",
                      "K_decline_1pct", "area_capped", "worst_case")
scenario_labels <- c(
  baseline       = "Baseline (v1)",
  si_proper      = "+ proper SI rejection (NFDS)",
  sc_inbreeding  = "+ SC + inbreeding depression",
  catastrophe    = "+ 5 %/yr recruit catastrophe",
  K_decline_1pct = "+ 1 %/yr K decline",
  area_capped    = "+ habitat-area cap on K",
  worst_case     = "Worst-case combined"
)

summ <- summ %>%
  left_join(metrics, by = "EO") %>%
  mutate(scenario = factor(scenario, levels = scenario_levels),
         scenario_label = scenario_labels[as.character(scenario)])

# Order EOs worst-to-best by initial P_compat (baseline y0)
eo_order <- summ %>%
  filter(scenario == "baseline", checkpoint_year == 10) %>%
  arrange(P_compat_p50) %>%
  pull(EO)
summ$EO <- factor(summ$EO, levels = eo_order)

# =============================================================================
# Figure 1 — Sensitivity tornado at y50
# =============================================================================
tornado <- summ %>%
  filter(checkpoint_year == 50, scenario != "baseline") %>%
  mutate(scenario_label = factor(scenario_label,
                                  levels = rev(scenario_labels[-1])),
         delta_sign = ifelse(delta_P_collapsed_vs_baseline > 0,
                             "increases collapse", "decreases collapse"))

tornado_plot <- ggplot(tornado, aes(x = delta_P_collapsed_vs_baseline,
                                      y = scenario_label,
                                      fill = delta_sign)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.1f%%",
                                 100 * delta_P_collapsed_vs_baseline)),
            hjust = ifelse(tornado$delta_P_collapsed_vs_baseline >= 0,
                            -0.15, 1.15),
            size = 3.2, colour = "grey15") +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey25") +
  facet_wrap(~ EO, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = c("increases collapse" = "#d73027",
                                "decreases collapse" = "#4575b4"),
                    name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0.18)) +
  labs(title = "PVA v1.5 sensitivity sweep — Δ P(collapse) at year 50 vs baseline",
       subtitle = paste0(
         "For each EO, change in probability of collapse (P_compat < 0.05) ",
         "under each alternative scenario, vs v1 baseline. Positive = scenario ",
         "drives more collapse; negative = scenario protects."),
       x = "Δ P(collapsed at y50) vs baseline", y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95", colour = "grey60"),
        legend.position = "bottom",
        panel.grid.major.y = element_blank())

dir.create("figures", showWarnings = FALSE)
ggsave("figures/SRK_TP1_pva_sensitivity_tornado.png", tornado_plot,
       width = 14, height = 7, dpi = 200)
ggsave("figures/SRK_TP1_pva_sensitivity_tornado.pdf", tornado_plot,
       width = 14, height = 7)
cat("Wrote sensitivity tornado\n")

# =============================================================================
# Figure 2 — P_compat at y50 per scenario per EO
# =============================================================================
p_compat_plot <- ggplot(summ %>% filter(checkpoint_year == 50),
                          aes(x = scenario_label, y = P_compat_p50,
                              colour = BL, fill = BL)) +
  geom_linerange(aes(ymin = P_compat_p10, ymax = P_compat_p90),
                  alpha = 0.55, linewidth = 0.6, show.legend = FALSE) +
  geom_point(shape = 21, size = 3.6, colour = "grey20", stroke = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "solid", colour = "#d73027",
             linewidth = 0.45) +
  geom_hline(yintercept = 0.40, linetype = "dashed", colour = "grey25",
             linewidth = 0.4) +
  scale_colour_manual(values = BL_COLORS, guide = "none",
                      breaks = BL_ORDER_NUMERIC) +
  scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage",
                    breaks = BL_ORDER_NUMERIC) +
  facet_wrap(~ EO, nrow = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "PVA v1.5 sensitivity sweep — P_compat at year 50 across scenarios",
       subtitle = paste0(
         "Per EO and scenario: median compatible-pair fraction (point) with ",
         "10–90 % envelope (vertical bar). Red line = collapse threshold ",
         "(0.05); dashed = random-mating viability (0.40)."),
       x = NULL, y = "P_compat (year 50)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
        strip.background = element_rect(fill = "grey95", colour = "grey60"))

ggsave("figures/SRK_TP1_pva_sensitivity_p_compat.png", p_compat_plot,
       width = 14, height = 8, dpi = 200)
ggsave("figures/SRK_TP1_pva_sensitivity_p_compat.pdf", p_compat_plot,
       width = 14, height = 8)
cat("Wrote sensitivity P_compat panel\n")

# =============================================================================
# Figure 3 — full sensitivity grid (P_collapsed colour-encoded heatmap)
# =============================================================================
grid_dat <- summ %>%
  filter(checkpoint_year %in% c(10, 20, 50)) %>%
  mutate(checkpoint_label = factor(paste0("Year ", checkpoint_year),
                                    levels = c("Year 10", "Year 20", "Year 50")))

grid_plot <- ggplot(grid_dat,
                     aes(x = scenario_label, y = EO, fill = P_collapsed)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * P_collapsed)),
            size = 3.2, colour = "grey15") +
  facet_wrap(~ checkpoint_label, nrow = 1) +
  scale_fill_gradient2(low = "#ffffff", mid = "#fee08b", high = "#d73027",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "P(collapsed)",
                       labels = scales::percent_format()) +
  labs(title = "PVA collapse probability — EO × scenario × checkpoint",
       subtitle = paste0(
         "P(P_compat < 0.05 by the checkpoint year) for each EO under each ",
         "scenario. 1000 Monte Carlo replicates per cell."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        panel.grid = element_blank())

ggsave("figures/SRK_TP1_pva_sensitivity_grid.png", grid_plot,
       width = 16, height = 5, dpi = 200)
ggsave("figures/SRK_TP1_pva_sensitivity_grid.pdf", grid_plot,
       width = 16, height = 5)
cat("Wrote sensitivity grid heatmap\n")

cat("\nAll v1.5 sensitivity figures written to figures/SRK_TP1_pva_sensitivity_*\n")
