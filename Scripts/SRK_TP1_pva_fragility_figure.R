#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_pva_fragility_figure.R — demographic fragility readout for v1.6 PVA
# =============================================================================
#
# Companion to SRK_TP1_pva_figures.R. v1.6 reports near-zero SI-collapse
# probability for most EOs because P_compat rises under NFDS selection — but
# the genetic-collapse threshold (P_compat < 0.05) doesn't capture demographic
# fragility. This figure surfaces it: per-EO probability that the adult
# census dips into critically-small territory at the checkpoints, plus the
# lifetime "did N ever crash" tracker.
#
# Inputs:
#   Tables/SRK_TP1_pva_prediction_summary.tsv  (with fragility columns added
#                                               by SRK_TP1_pva_prediction.py)
#
# Outputs:
#   figures/SRK_TP1_pva_fragility.{png,pdf}
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("srk_bl_constants.R")

summ <- read_tsv("Tables/SRK_TP1_pva_prediction_summary.tsv",
                  show_col_types = FALSE)

metrics <- read_tsv("Tables/SRK_EO_allele_richness.tsv",
                     show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(EO = paste0("EO", group)) %>%
  select(EO, BL, P_compat_0 = P_compat_L0)

# Worst-to-best by initial P_compat (same ordering as other PVA figures)
eo_order <- metrics %>% arrange(P_compat_0) %>% pull(EO)

frag <- summ %>%
  select(EO, checkpoint_year, P_collapsed,
         P_N_lt50, P_N_lt10, P_minN_lt10_to_date, P_adults_extinct_ever) %>%
  mutate(across(-c(EO, checkpoint_year), as.numeric)) %>%
  left_join(metrics %>% select(EO, BL), by = "EO") %>%
  mutate(EO = factor(EO, levels = eo_order))

frag_long <- frag %>%
  pivot_longer(
    cols = c(P_collapsed, P_N_lt50, P_N_lt10,
             P_minN_lt10_to_date, P_adults_extinct_ever),
    names_to = "metric", values_to = "prob"
  ) %>%
  mutate(
    metric = factor(metric, levels = c(
      "P_collapsed",
      "P_N_lt50",
      "P_N_lt10",
      "P_minN_lt10_to_date",
      "P_adults_extinct_ever"
    ), labels = c(
      "SI collapse (P_compat < 0.05)",
      "N < 50 at checkpoint (chronic)",
      "N < 10 at checkpoint (acute)",
      "min N < 10 by checkpoint (near-miss)",
      "adults extinct at any point (N = 0)"
    )),
    checkpoint_label = factor(paste0("Year ", checkpoint_year),
                              levels = c("Year 10", "Year 20", "Year 50"))
  )

frag_plot <- ggplot(frag_long,
                     aes(x = metric, y = EO, fill = prob)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  geom_text(aes(label = ifelse(prob >= 0.005,
                                sprintf("%.0f%%", 100 * prob),
                                "0%")),
            size = 3.2, colour = "grey15") +
  facet_wrap(~ checkpoint_label, nrow = 1) +
  scale_fill_gradient2(low = "#ffffff", mid = "#fee08b", high = "#d73027",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Probability",
                       labels = scales::percent_format()) +
  labs(title = "PVA v1.6 (empirical per-EO model): EO67 is critically fragile; EO70 chronically small",
       subtitle = paste0(
         "Five failure modes per EO and checkpoint. SI collapse alone (leftmost column) ",
         "misses demographic fragility: EO67 stays above the SI-collapse threshold but ",
         "76 % of replicates dip below 10 adults at some point in 50 yr (median min N = 6). ",
         "EO70 is chronically small (40 % of replicates < 50 adults at y50) but rarely ",
         "near-extinct. EO18, 25, 27, 76 are demographically robust on every axis. ",
         "1000 MC replicates per cell."),
       x = NULL, y = NULL,
       caption = "Source: Tables/SRK_TP1_pva_prediction_summary.tsv (v1.6) / 1000 reps / 50-year horizon") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 8.5),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        panel.grid = element_blank(),
        plot.caption = element_text(colour = "grey45", size = 8))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/SRK_TP1_pva_v16_demographic_fragility.png", frag_plot,
       width = 14, height = 5.5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v16_demographic_fragility.pdf", frag_plot,
       width = 14, height = 5.5)
# Remove the old unversioned filename to avoid confusion
for (ext in c("png", "pdf")) {
  old <- paste0("figures/SRK_TP1_pva_fragility.", ext)
  if (file.exists(old)) file.remove(old)
}
cat("Wrote v1.6 demographic fragility figure (renamed from SRK_TP1_pva_fragility.*)\n")

# Also print a tidy console summary for the user
frag_console <- frag %>%
  filter(checkpoint_year == 50) %>%
  select(EO, BL, P_collapsed, P_N_lt50, P_N_lt10,
         P_minN_lt10_to_date, P_adults_extinct_ever)
cat("\nFragility at year 50:\n")
print(as.data.frame(frag_console), row.names = FALSE, digits = 3)
