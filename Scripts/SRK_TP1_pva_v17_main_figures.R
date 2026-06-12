#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_pva_v17_main_figures.R — canonical PVA figures (v1.7 only)
# =============================================================================
#
# v1.7 (maternal-dormancy inheritance) is the canonical PVA going forward.
# This script reads Tables/SRK_TP1_pva_v17_summary.tsv and produces three
# focused figures + the lineage-fragility heatmap (built separately by
# SRK_TP1_pva_v17_figure.R, retained for the maternal-class breakdown).
#
# Inputs:
#   Tables/SRK_TP1_pva_v17_summary.tsv          (v1.7 canonical)
#   Tables/SRK_EO_allele_richness.tsv           (for year-0 baselines)
#
# Outputs:
#   figures/SRK_TP1_pva_v17_population_size.{png,pdf}     — N trajectory
#   figures/SRK_TP1_pva_v17_P_compat.{png,pdf}            — P_compat trajectory
#   figures/SRK_TP1_pva_v17_collapse_probability.{png,pdf} — P(collapse) per EO
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("srk_bl_constants.R")

P_FLOOR  <- 0.05
P_VIABLE <- 0.40
CHECKPOINTS <- c(10, 20, 50)
HEADLINE_EO <- "EO67"   # 2026-05-28: K_area cap dropped; EO70 collapse risk
                        # fell from 14.7 % to 0.6 %.  EO67 (chronically <50
                        # adults, 83 % near-miss N<10) is now the standout.

# ----- Load -----------------------------------------------------------------
v17 <- read_tsv("Tables/SRK_TP1_pva_v17_summary.tsv",
                 show_col_types = FALSE) %>%
  mutate(across(-c(EO, checkpoint_year), as.numeric))

metrics <- read_tsv("Tables/SRK_EO_allele_richness.tsv",
                     show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(EO = paste0("EO", group)) %>%
  select(EO, BL, P_compat_0 = P_compat_L0,
         k_obs_0 = k_observed, N_0 = N)

# Year-0 baseline rows to anchor trajectories.
baseline <- metrics %>%
  transmute(EO,
            checkpoint_year = 0L,
            P_collapsed = 0,
            N_adult_p10 = N_0, N_adult_p50 = N_0, N_adult_p90 = N_0,
            P_compat_p10 = P_compat_0, P_compat_p50 = P_compat_0,
            P_compat_p90 = P_compat_0)

traj <- bind_rows(baseline, v17) %>%
  arrange(EO, checkpoint_year) %>%
  left_join(metrics %>% select(EO, BL), by = "EO")

# Order EOs worst-to-best by initial P_compat (EO70 worst → most visible).
eo_order <- metrics %>% arrange(P_compat_0) %>% pull(EO)
traj$EO <- factor(traj$EO, levels = eo_order)

# Highlight EO70 (headline) with a slightly heavier line + bigger points.
is_headline <- function(EO) EO == HEADLINE_EO

# Reusable trajectory builder.
build_traj <- function(metric, ylab, ymin = 0, ymax = NA_real_,
                        hlines = NULL, log_y = FALSE) {
  p10 <- paste0(metric, "_p10")
  p50 <- paste0(metric, "_p50")
  p90 <- paste0(metric, "_p90")
  dat <- traj %>%
    transmute(EO, checkpoint_year, BL,
              p10 = .data[[p10]], p50 = .data[[p50]], p90 = .data[[p90]],
              headline = is_headline(EO))
  g <- ggplot(dat, aes(x = checkpoint_year, y = p50,
                        colour = BL, fill = BL,
                        size = headline, linewidth = headline)) +
    geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.18, colour = NA) +
    geom_line() +
    geom_point(aes(size = headline), shape = 21, colour = "grey20",
               stroke = 0.3) +
    facet_wrap(~ EO, nrow = 2) +
    scale_colour_manual(values = BL_COLORS, name = "Bottleneck lineage",
                         breaks = BL_ORDER_NUMERIC) +
    scale_fill_manual(values = BL_COLORS, guide = "none",
                       breaks = BL_ORDER_NUMERIC) +
    scale_size_manual(values = c(`TRUE` = 3.4, `FALSE` = 2.2), guide = "none") +
    scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.85),
                            guide = "none") +
    scale_x_continuous(breaks = c(0, CHECKPOINTS),
                        expand = expansion(mult = 0.02))
  if (log_y) {
    g <- g + scale_y_log10()
  } else {
    g <- g + coord_cartesian(ylim = c(ymin, ymax))
  }
  g <- g +
    labs(x = "Years from now", y = ylab) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          strip.background = element_rect(fill = "grey95", colour = "grey60"),
          legend.position = "right",
          plot.caption = element_text(colour = "grey45", size = 8))
  if (!is.null(hlines)) {
    for (hl in hlines) {
      g <- g + geom_hline(yintercept = hl$y, linetype = hl$lt,
                          colour = hl$col, linewidth = 0.4)
    }
  }
  g
}

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Figure 1: Population size N(t)
# -----------------------------------------------------------------------------
N_plot <- build_traj("N_adult", "Adult census N", ymin = 0, log_y = TRUE) +
  labs(title = "PVA v1.7: population size over 50 years — EO67 the only chronically small population",
       subtitle = paste0(
         "Median + 80 % envelope across 1000 Monte Carlo replicates per EO. ",
         "K = observed census × 1.5, with per-EO empirical σ_env. ",
         "Log y-axis. EO67 (N=20 starter) equilibrates around 48 adults — the ",
         "only EO chronically below 100. EO70 stabilises around 534 once K is ",
         "anchored to census; the other four EOs settle at 1000-2700 medians."),
       caption = "Source: Tables/SRK_TP1_pva_v17_summary.tsv / 1000 reps / 50-year horizon")
ggsave("figures/SRK_TP1_pva_v17_population_size.png", N_plot,
       width = 12, height = 6.5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v17_population_size.pdf", N_plot,
       width = 12, height = 6.5)
cat("Wrote v1.7 population-size trajectory figure\n")

# -----------------------------------------------------------------------------
# Figure 2: P_compat(t)
# -----------------------------------------------------------------------------
P_plot <- build_traj("P_compat", "Compatible-pair fraction (P_compat)",
                      ymin = 0, ymax = 1,
                      hlines = list(list(y = P_FLOOR, lt = "solid", col = "#d73027"),
                                    list(y = P_VIABLE, lt = "dashed", col = "grey25"))) +
  labs(title = "PVA v1.7: P_compat rises under NFDS in every EO; all medians stay above the collapse floor",
       subtitle = paste0(
         "Median + 80 % envelope per EO. Red line = collapse threshold (P_compat = 0.05); ",
         "dashed line = random-mating viability (0.40). Proper SI rejection pushes allele ",
         "frequencies toward the NFDS attractor, raising P_compat monotonically. EO70 (6 ",
         "alleles) and EO67 (11) sit below the 0.40 viability line; the four large EOs sit ",
         "comfortably above it. No EO's lower envelope touches the collapse floor."),
       caption = "Source: Tables/SRK_TP1_pva_v17_summary.tsv / 1000 reps / 50-year horizon")
ggsave("figures/SRK_TP1_pva_v17_P_compat.png", P_plot,
       width = 12, height = 6.5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v17_P_compat.pdf", P_plot,
       width = 12, height = 6.5)
cat("Wrote v1.7 P_compat trajectory figure\n")

# -----------------------------------------------------------------------------
# Figure 3: P(collapse) per EO and checkpoint — bar chart
# -----------------------------------------------------------------------------
collapse_df <- v17 %>%
  select(EO, checkpoint_year, P_collapsed) %>%
  left_join(metrics %>% select(EO, BL), by = "EO") %>%
  mutate(EO = factor(EO, levels = eo_order),
         checkpoint_label = factor(paste0("Year ", checkpoint_year),
                                    levels = c("Year 10", "Year 20", "Year 50")),
         is_headline = EO == HEADLINE_EO)

highlight_strip <- data.frame(EO = factor(HEADLINE_EO, levels = eo_order))

collapse_plot <- ggplot(collapse_df, aes(x = P_collapsed, y = EO, fill = BL)) +
  geom_tile(data = highlight_strip,
            aes(x = 0.10, y = EO), width = Inf, height = 1,
            fill = "#fff3b0", colour = NA, inherit.aes = FALSE) +
  geom_col(width = 0.65, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * P_collapsed)),
            hjust = -0.15, size = 3.4, colour = "grey15") +
  facet_grid(. ~ checkpoint_label) +
  scale_x_continuous(limits = c(0, 0.20), labels = scales::percent_format(),
                      expand = expansion(mult = c(0.005, 0.18))) +
  scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage",
                     breaks = BL_ORDER_NUMERIC) +
  labs(title = "PVA v1.7: no EO has meaningful 50-year SI-collapse risk under realistic K",
       subtitle = paste0(
         "Probability that a replicate's P_compat falls below 0.05 by the checkpoint year, ",
         "out of 1000 Monte Carlo replicates. With K = census x 1.5 (no area cap), only ",
         "EO67 (yellow strip) shows a non-zero rate (0.5 % at y50), and it is essentially ",
         "the Monte Carlo noise floor. EO67's real concern is demographic, not genetic: ",
         "83 % of replicates dip below 10 adults at some point (see fragility figure)."),
       x = "Probability of collapse at the checkpoint", y = NULL,
       caption = "Source: Tables/SRK_TP1_pva_v17_summary.tsv / 1000 reps / 50-year horizon") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        legend.position = "right",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(
          face = ifelse(levels(collapse_df$EO) == HEADLINE_EO, "bold", "plain")
        ),
        plot.caption = element_text(colour = "grey45", size = 8))

ggsave("figures/SRK_TP1_pva_v17_collapse_probability.png", collapse_plot,
       width = 13, height = 5.5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v17_collapse_probability.pdf", collapse_plot,
       width = 13, height = 5.5)
cat("Wrote v1.7 collapse-probability figure\n")

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
cat("\nv1.7 headline numbers (year 50):\n")
v17 %>% filter(checkpoint_year == 50) %>%
  transmute(EO,
            N_p10_p50_p90 = sprintf("%.0f / %.0f / %.0f",
                                     N_adult_p10, N_adult_p50, N_adult_p90),
            P_compat_p10_p50_p90 = sprintf("%.2f / %.2f / %.2f",
                                            P_compat_p10, P_compat_p50, P_compat_p90),
            P_collapsed = sprintf("%.1f%%", 100 * P_collapsed),
            P_minN_lt10 = sprintf("%.1f%%", 100 * P_minN_lt10_to_date)) %>%
  as.data.frame() %>% print(row.names = FALSE)
