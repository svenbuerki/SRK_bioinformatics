#!/usr/bin/env Rscript
# =============================================================================
# SRK_depletion_ranking_presentation.R — Step 17 PRESENTATION-ONLY variant
# =============================================================================
#
# Identical analysis to SRK_depletion_ranking.R but renders SLIDE-READY
# figures with:
#   * no title, no subtitle, no caption
#   * no legend
#   * identical canvas dimensions across all three variants so the figures
#     drop into a slide deck without resizing artefacts
#
# Quadrant colour code (urgency-graded; aligned with the breeding-strategies
# vocabulary). Kept in sync with the canonical script:
#   HEALTHY                                  blue   (safe)
#   INFORMED BREEDING (frequency skew)       orange (caution)
#   INFORMED BREEDING + INJECTION preventive light red (warning)
#   INFORMED BREEDING + INJECTION urgent     dark red (emergency)
#
# Inputs:
#   Tables/SRK_EO_allele_richness.tsv
#
# Outputs (3 variants x 2 panels = 6 figures, all 9 x 8 in):
#   figures/presentation/SRK_depletion_ranking_{observed,predicted}_{blank,EOs,all}.{png,pdf}
#
# NOT referenced from any online doc — these are presentation-only.
# Canonical, fully labelled version lives in SRK_depletion_ranking.R.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) message("Note: install ggrepel for non-overlapping labels.")

source("srk_bl_constants.R")

DI_THRESHOLD    <- 0.50
PCOMP_THRESHOLD <- 0.40

# Slide canvas
PRES_WIDTH  <- 9
PRES_HEIGHT <- 8
PRES_DIR    <- "figures/presentation"

cat("Loading SRK_EO_allele_richness.tsv ...\n")
dat <- read_tsv("Tables/SRK_EO_allele_richness.tsv", show_col_types = FALSE) %>%
  mutate(
    level = factor(level, levels = c("EO", "BL")),
    BL    = factor(BL,    levels = BL_ORDER_NUMERIC),
    label = ifelse(level == "EO", paste0("EO", group), group)
  )

PANEL_NAMES <- c("HEALTHY", "BREEDING", "PREVENT", "URGENT")

build_panel_pres <- function(di_col, show = "all", panels = PANEL_NAMES) {
  stopifnot(show %in% c("none", "EO", "all"))
  stopifnot(all(panels %in% PANEL_NAMES))
  pdat <- dat %>%
    filter(!is.na(.data[[di_col]])) %>%
    mutate(DI = .data[[di_col]],
           P_compat = P_compat_L0,
           P_lo = P_compat_L0_lo,
           P_hi = P_compat_L0_hi)
  pdat_show <- switch(show,
    "none" = pdat[0, ],
    "EO"   = pdat %>% filter(level == "EO"),
    "all"  = pdat
  )

  g <- ggplot(pdat, aes(x = DI, y = P_compat))
  if ("HEALTHY" %in% panels) {
    g <- g +
      annotate("rect", xmin = 0, xmax = DI_THRESHOLD,
               ymin = PCOMP_THRESHOLD, ymax = 1,
               fill = "#4575b4", alpha = 0.08) +
      annotate("text", x = 0.015, y = 0.985,
               label = "HEALTHY\n(random mating)",
               colour = "#4575b4", fontface = "bold", size = 4,
               hjust = 0, vjust = 1, lineheight = 0.95)
  }
  if ("PREVENT" %in% panels) {
    g <- g +
      annotate("rect", xmin = DI_THRESHOLD, xmax = 1,
               ymin = PCOMP_THRESHOLD, ymax = 1,
               fill = "#fc9272", alpha = 0.14) +
      annotate("text", x = 0.985, y = 0.985,
               label = "INFORMED BREEDING\n+ ALLELE INJECTION\n(preventive)",
               colour = "#a50f15", fontface = "bold", size = 3.8,
               hjust = 1, vjust = 1, lineheight = 0.95)
  }
  if ("BREEDING" %in% panels) {
    g <- g +
      annotate("rect", xmin = 0, xmax = DI_THRESHOLD,
               ymin = 0, ymax = PCOMP_THRESHOLD,
               fill = "#f16913", alpha = 0.20) +
      annotate("text", x = 0.015, y = 0.015,
               label = "INFORMED BREEDING\n(frequency skew)",
               colour = "#8c2d04", fontface = "bold", size = 3.8,
               hjust = 0, vjust = 0, lineheight = 0.95)
  }
  if ("URGENT" %in% panels) {
    g <- g +
      annotate("rect", xmin = DI_THRESHOLD, xmax = 1,
               ymin = 0, ymax = PCOMP_THRESHOLD,
               fill = "#cb181d", alpha = 0.20) +
      annotate("text", x = 0.985, y = PCOMP_THRESHOLD - 0.015,
               label = "INFORMED BREEDING\n+ ALLELE INJECTION\n(urgent)",
               colour = "#67000d", fontface = "bold", size = 3.8,
               hjust = 1, vjust = 1, lineheight = 0.95)
  }
  g <- g +
    geom_vline(xintercept = DI_THRESHOLD, linetype = "dashed",
               colour = "grey25", linewidth = 0.6) +
    geom_hline(yintercept = PCOMP_THRESHOLD, linetype = "dashed",
               colour = "grey25", linewidth = 0.6) +
    geom_linerange(data = pdat_show,
                   aes(ymin = P_lo, ymax = P_hi, colour = BL),
                   linewidth = 0.6, alpha = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = BL_COLORS, guide = "none") +
    geom_point(data = pdat_show,
               aes(fill = BL, shape = level),
               size = 5.5, stroke = 0.9, colour = "grey20", alpha = 0.95) +
    scale_fill_manual(values = BL_COLORS, guide = "none",
                      breaks = BL_ORDER_NUMERIC, drop = FALSE) +
    scale_shape_manual(values = c(EO = 21, BL = 24), guide = "none",
                       drop = FALSE) +
    {if (use_repel && nrow(pdat_show) > 0)
      ggrepel::geom_text_repel(
        data = pdat_show,
        aes(label = label,
            fontface = ifelse(level == "BL", "bold.italic", "bold")),
        size = 4.2, show.legend = FALSE,
        box.padding = 0.55, point.padding = 0.4,
        seed = 42, max.overlaps = 20)
     else if (nrow(pdat_show) > 0)
      geom_text(data = pdat_show, aes(label = label), vjust = -1.2, size = 4,
                fontface = "bold", show.legend = FALSE)
     else NULL} +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    labs(
      x = "Depletion Index (1 = no alleles, 0 = at species equilibrium)",
      y = "Compatible-pair fraction (fraction of random pairs that produce seed)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption = element_blank())
}

dir.create(PRES_DIR, recursive = TRUE, showWarnings = FALSE)

render_set_pres <- function(di_col, base_stem) {
  variants <- list(blank = "none", EOs = "EO", all = "all")
  for (label in names(variants)) {
    p <- build_panel_pres(di_col, show = variants[[label]])
    out_stem <- file.path(PRES_DIR, paste0(base_stem, "_", label))
    ggsave(paste0(out_stem, ".png"), plot = p,
           width = PRES_WIDTH, height = PRES_HEIGHT, dpi = 200)
    ggsave(paste0(out_stem, ".pdf"), plot = p,
           width = PRES_WIDTH, height = PRES_HEIGHT)
    cat(sprintf("  Written: %s.{png,pdf}\n", out_stem))
  }
}

render_set_pres("DI_observed",  "SRK_depletion_ranking_observed")
render_set_pres("DI_predicted", "SRK_depletion_ranking_predicted")

# -----------------------------------------------------------------------------
# Animation frames (predicted only) — progressive panel reveal for slide builds.
# Each frame shows the axes + selected coloured quadrants, no data points.
#   _blank_healthy           = top-left blue panel only
#   _blank_healthy_breeding  = top-left blue + bottom-left orange
# Combine with the existing _blank (all four panels) and _EOs / _all (with
# points) to build a 5-step animation in the presentation.
# -----------------------------------------------------------------------------
animation_frames <- list(
  blank_healthy          = c("HEALTHY"),
  blank_healthy_breeding = c("HEALTHY", "BREEDING")
)
for (label in names(animation_frames)) {
  p <- build_panel_pres("DI_predicted", show = "none",
                         panels = animation_frames[[label]])
  out_stem <- file.path(PRES_DIR,
                        paste0("SRK_depletion_ranking_predicted_", label))
  ggsave(paste0(out_stem, ".png"), plot = p,
         width = PRES_WIDTH, height = PRES_HEIGHT, dpi = 200)
  ggsave(paste0(out_stem, ".pdf"), plot = p,
         width = PRES_WIDTH, height = PRES_HEIGHT)
  cat(sprintf("  Written: %s.{png,pdf}\n", out_stem))
}

cat(sprintf("\nPresentation figures written to %s/ (%d x %d in).\n",
            PRES_DIR, PRES_WIDTH, PRES_HEIGHT))
