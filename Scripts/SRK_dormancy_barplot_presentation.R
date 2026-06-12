#!/usr/bin/env Rscript
# =============================================================================
# SRK_dormancy_barplot_presentation.R — slide-ready dormancy stacked bars
# =============================================================================
#
# Six EOs (27, 70, 18, 25, 67, 76) shown as stacked composition of:
#   * non-dormant (green)
#   * dormant     (red; = PD + PH merged from PVA DORMANCY_SPLIT)
#   * non-viable  (black)
#
# Source values are the visual estimates already used in the PVA simulator
# (SRK_TP1_pva_v17_maternal_dormancy.py :: DORMANCY_SPLIT). EO70 — for which
# no germination assay exists — is overridden to a 50/50 non-dormant/dormant
# placeholder per the presentation brief.
#
# Bars are ordered by non-dormant fraction ascending (worst → best).
#
# Output:
#   figures/presentation/SRK_dormancy_barplot.{png,pdf}   9 x 7 in
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

PRES_DIR    <- "figures/presentation"
PRES_WIDTH  <- 9
PRES_HEIGHT <- 7

COL_ND <- "#90EE90"   # non-dormant — R "lightgreen", matches original figure
COL_D  <- "#f57c00"   # dormant     — orange (PD + PH merged); avoids drift-red clash
COL_NV <- "#000000"   # non-viable  — black

# (non_dormant, dormant, non_viable) — dormant = PD + PH from DORMANCY_SPLIT.
# EO70 is the user-supplied placeholder (no assay data).
raw <- tribble(
  ~EO,    ~non_dormant, ~dormant, ~non_viable,
  "EO27", 0.25,         0.60,     0.15,
  "EO70", 0.50,         0.50,     0.00,
  "EO18", 0.60,         0.40,     0.00,
  "EO25", 0.67,         0.31,     0.02,
  "EO67", 0.68,         0.30,     0.02,
  "EO76", 0.93,         0.07,     0.00
)

eo_order <- raw %>% arrange(non_dormant) %>% pull(EO)

long <- raw %>%
  pivot_longer(c(non_dormant, dormant, non_viable),
               names_to = "category", values_to = "fraction") %>%
  mutate(
    EO       = factor(EO, levels = eo_order),
    category = factor(category,
                      levels = c("non_viable", "dormant", "non_dormant"),
                      labels = c("Non-viable", "Dormant", "Non-dormant"))
  )

g <- ggplot(long, aes(x = EO, y = fraction, fill = category)) +
  geom_col(position = "stack", width = 0.78, colour = "grey20",
           linewidth = 0.35) +
  scale_fill_manual(values = c("Non-dormant" = COL_ND,
                                "Dormant"     = COL_D,
                                "Non-viable"  = COL_NV),
                    breaks = c("Non-dormant", "Dormant", "Non-viable"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Soil Seed Bank Proportion (%)") +
  theme_classic(base_size = 16) +
  theme(
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    legend.key.size   = unit(0.8, "cm"),
    axis.text.x       = element_text(size = 16, face = "bold"),
    axis.text.y       = element_text(size = 14),
    axis.title.y      = element_text(size = 15, margin = margin(r = 8)),
    axis.line         = element_line(linewidth = 0.5),
    axis.ticks        = element_line(linewidth = 0.5),
    plot.margin       = margin(12, 14, 10, 10)
  )

dir.create(PRES_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(PRES_DIR, "SRK_dormancy_barplot.png"), g,
       width = PRES_WIDTH, height = PRES_HEIGHT, dpi = 300)
ggsave(file.path(PRES_DIR, "SRK_dormancy_barplot.pdf"), g,
       width = PRES_WIDTH, height = PRES_HEIGHT)
cat("Wrote dormancy presentation barplot (PNG + PDF) to ", PRES_DIR, "\n", sep = "")
