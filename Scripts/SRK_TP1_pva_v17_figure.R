#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_pva_v17_figure.R — v1.7 maternal-dormancy lineage-fragility figure
# =============================================================================
#
# Companion to SRK_TP1_pva_fragility_figure.R (which covers v1.6 demographic
# fragility). This figure adds the v1.7-specific lineage-extinction columns:
# probability that a given maternal dormancy class has zero mothers at the
# checkpoint year.
#
# Reads:
#   Tables/SRK_TP1_pva_v17_summary.tsv         (v1.7 output)
#   Tables/SRK_TP1_pva_prediction_summary.tsv  (v1.6 output, for direct
#                                                P_collapse comparison)
#
# Writes:
#   figures/SRK_TP1_pva_v17_lineage_fragility.{png,pdf}
#   figures/SRK_TP1_pva_v16_vs_v17_collapse.{png,pdf}
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("srk_bl_constants.R")

# v1.7 source TSV. Set env SRK_V17_TSV to override (e.g., the stochastic-init
# backup at SRK_TP1_pva_v17_stoinit_summary.tsv) when the canonical file is
# stale or being regenerated.
.v17_path <- Sys.getenv("SRK_V17_TSV", unset = "Tables/SRK_TP1_pva_v17_summary.tsv")
cat("Reading v1.7 from:", .v17_path, "\n")
v17 <- read_tsv(.v17_path, show_col_types = FALSE)
v16 <- read_tsv("Tables/SRK_TP1_pva_prediction_summary.tsv",
                show_col_types = FALSE)

metrics <- read_tsv("Tables/SRK_EO_allele_richness.tsv",
                     show_col_types = FALSE) %>%
  filter(level == "EO") %>%
  mutate(EO = paste0("EO", group)) %>%
  select(EO, BL, P_compat_0 = P_compat_L0)

eo_order <- metrics %>% arrange(P_compat_0) %>% pull(EO)

# EO dormancy shares (from the script DORMANCY_SPLIT) — used to mask cells
# where a class was never present at init (so "no X mothers" is trivial).
dormancy_share <- tribble(
  ~EO,     ~f_ND, ~f_PD, ~f_PH, ~f_NV,
  "EO18",  0.60,  0.13,  0.27,  0.00,
  "EO25",  0.67,  0.08,  0.23,  0.02,
  "EO27",  0.25,  0.40,  0.20,  0.15,
  "EO67",  0.68,  0.08,  0.22,  0.02,
  "EO70",  0.46,  0.15,  0.36,  0.03,
  "EO76",  0.93,  0.00,  0.07,  0.00
)

# -----------------------------------------------------------------------------
# Figure 1 — lineage fragility heatmap (per EO × class × checkpoint)
# -----------------------------------------------------------------------------
lineage <- v17 %>%
  select(EO, checkpoint_year,
         P_no_ND_mothers, P_no_PD_mothers, P_no_PH_mothers) %>%
  mutate(across(starts_with("P_no_"), as.numeric)) %>%
  pivot_longer(cols = starts_with("P_no_"),
               names_to = "class", values_to = "prob") %>%
  mutate(class = sub("^P_no_(.*)_mothers$", "\\1", class)) %>%
  left_join(dormancy_share %>%
              pivot_longer(cols = starts_with("f_"),
                           names_to = "class", values_to = "init_share") %>%
              mutate(class = sub("^f_", "", class)),
            by = c("EO", "class")) %>%
  mutate(EO = factor(EO, levels = eo_order),
         class = factor(class, levels = c("ND", "PD", "PH"),
                        labels = c("Non-dormant (ND)",
                                   "Physiological dormancy (PD)",
                                   "Physical dormancy (PH)")),
         checkpoint_label = factor(paste0("Year ", checkpoint_year),
                                   levels = c("Year 10", "Year 20", "Year 50")),
         # Mask cells where class was never present at init (init_share = 0).
         displayed = ifelse(is.na(init_share) | init_share == 0, NA, prob))

lineage_plot <- ggplot(lineage,
                       aes(x = class, y = EO, fill = displayed)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.na(displayed), "n/a",
                               ifelse(displayed >= 0.005,
                                      sprintf("%.0f%%", 100 * displayed),
                                      "0%"))),
            size = 3.2, colour = "grey15") +
  facet_wrap(~ checkpoint_label, nrow = 1) +
  scale_fill_gradient2(low = "#ffffff", mid = "#fee08b", high = "#d73027",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "P(no mothers\nof this class)",
                       labels = scales::percent_format(),
                       na.value = "grey90") +
  labs(title = "PVA v1.7 (maternal-dormancy model): EO67 loses both dormant lineages in most replicates",
       subtitle = paste0(
         "Probability that a given dormancy class has zero mother plants at the checkpoint, ",
         "1000 MC replicates. EO67 is the headline: ~50 % of replicates lose the PD lineage and ",
         "~87 % lose PH by year 50 (population is too small to buffer). EO70 shows modest 9-10 % ",
         "lineage loss that propagates to genetic collapse. 'n/a' = class never present at this EO ",
         "initially. Lineage loss means the corresponding seed-bank compartment can no longer refill."),
       x = "Maternal dormancy class", y = NULL,
       caption = "Source: SRK_TP1_pva_v17_maternal_dormancy.py / 1000 reps / 50-year horizon") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.caption = element_text(colour = "grey45", size = 8))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/SRK_TP1_pva_v17_lineage_fragility.png", lineage_plot,
       width = 13, height = 5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v17_lineage_fragility.pdf", lineage_plot,
       width = 13, height = 5)
cat("Wrote v1.7 lineage fragility figure\n")

# -----------------------------------------------------------------------------
# Figure 2 — v1.6 vs v1.7 collapse-probability comparison
# -----------------------------------------------------------------------------
comp <- bind_rows(
  v16 %>% select(EO, checkpoint_year, P_collapsed) %>%
    mutate(version = "v1.6 (fixed dormancy split)"),
  v17 %>% select(EO, checkpoint_year, P_collapsed) %>%
    mutate(version = "v1.7 (maternal-inherited dormancy)")
) %>%
  mutate(P_collapsed = as.numeric(P_collapsed),
         EO = factor(EO, levels = eo_order),
         checkpoint_label = factor(paste0("Year ", checkpoint_year),
                                   levels = c("Year 10", "Year 20", "Year 50")))

# Highlight strip for EO70 row (the headline).
highlight_strip <- data.frame(EO = factor("EO70", levels = eo_order))

comp_plot <- ggplot(comp, aes(x = P_collapsed, y = EO, fill = version)) +
  geom_tile(data = highlight_strip,
            aes(x = 0.10, y = EO), width = Inf, height = 1,
            fill = "#fff3b0", colour = NA, inherit.aes = FALSE) +
  geom_col(position = position_dodge2(reverse = TRUE), width = 0.7,
           colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * P_collapsed)),
            position = position_dodge2(width = 0.7, reverse = TRUE),
            hjust = -0.15, size = 3.2, colour = "grey15") +
  facet_grid(. ~ checkpoint_label) +
  scale_x_continuous(limits = c(0, NA), labels = scales::percent_format(),
                     expand = expansion(mult = c(0.005, 0.22))) +
  scale_fill_manual(values = c("v1.6 (fixed dormancy split)" = "#9ecae1",
                                "v1.7 (maternal-inherited dormancy)" = "#d73027"),
                    name = NULL) +
  labs(title = "v1.7 reveals 14x collapse risk for EO70 that v1.6 missed",
       subtitle = paste0(
         "Probability of genetic collapse (P_compat < 0.05) at each checkpoint, ",
         "comparing the two PVA models. v1.6 (light blue, baseline) holds dormancy ",
         "as a fixed environmental EO parameter. v1.7 (red) treats dormancy as a ",
         "heritable maternal trait, so the seed-bank refill rate is contingent on ",
         "which mothers survive each year. The yellow strip marks EO70: in v1.7, ",
         "stochastic loss of the PD/PH dormant-lineage mothers drains the seed-bank ",
         "buffer over ~20 years, raising 50-yr collapse risk from 1.1 % to 16 %. ",
         "Other EOs are unaffected because their populations are large enough to ",
         "preserve all maternal lineages, or because P_compat is buffered by ",
         "higher allele richness (EO67)."),
       x = "Probability of collapse (P_compat < 0.05)", y = NULL,
       caption = paste0("Source: Tables/SRK_TP1_pva_prediction_summary.tsv (v1.6) + ",
                        "Tables/SRK_TP1_pva_v17_summary.tsv (v1.7) / 1000 reps each")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "grey95", colour = "grey60"),
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(
          face = ifelse(levels(comp$EO) == "EO70", "bold", "plain")
        ),
        plot.caption = element_text(colour = "grey45", size = 8))

ggsave("figures/SRK_TP1_pva_v16_vs_v17_collapse.png", comp_plot,
       width = 13, height = 5.5, dpi = 200)
ggsave("figures/SRK_TP1_pva_v16_vs_v17_collapse.pdf", comp_plot,
       width = 13, height = 5.5)
cat("Wrote v1.6 vs v1.7 collapse comparison figure\n")

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------
cat("\nv1.7 lineage status at year 50:\n")
v17 %>% filter(checkpoint_year == 50) %>%
  select(EO, P_collapsed,
         P_no_ND_mothers, P_no_PD_mothers, P_no_PH_mothers,
         n_ND_mothers_p50, n_PD_mothers_p50, n_PH_mothers_p50) %>%
  as.data.frame() %>% print(row.names = FALSE)
