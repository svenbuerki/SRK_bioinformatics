#!/usr/bin/env Rscript
# =============================================================================
# SRK_depletion_ranking.R — Step 17 companion
# Compatible-pair fraction x Depletion Index ranking, with EO + BL panels
# =============================================================================
#
# Builds on the Step 17 metrics TSV. Plots each group on a 2D axis of:
#   x = Depletion Index (DI) — distance from species k-equilibrium
#   y = Compatible-pair fraction (P_compat, strict SI)
#
# Two figures are produced — one using DI_observed (k_rarefied30 / k_species)
# and one using DI_predicted (MM-extrapolated k / k_species). The observed
# version is the conservative read; the predicted version shows the upper
# bound on achievable richness without further intervention. EOs/BLs whose
# observed and predicted DI differ substantially are populations where more
# sampling would likely reveal additional alleles.
#
# Quadrants map directly to the breeding-strategies table:
#
#                  DI < 0.50 (mild)            DI >= 0.50 (severe)
#   P_compat >= 0.40   HEALTHY                   INFORMED BREEDING
#                      (random mating)           + ALLELE INJECTION
#                                                (preventive)
#   P_compat <  0.40   INFORMED BREEDING         INFORMED BREEDING
#                      (frequency skew)          + ALLELE INJECTION
#                                                (urgent)
#
# Quadrant colour code (urgency-graded; aligned with the breeding-strategies
# vocabulary):
#   HEALTHY                                  blue   (safe)
#   INFORMED BREEDING (frequency skew)       orange (caution)
#   INFORMED BREEDING + INJECTION preventive light red (warning)
#   INFORMED BREEDING + INJECTION urgent     dark red (emergency)
#
# Inputs:
#   Tables/SRK_EO_allele_richness.tsv
#
# Outputs (3 variants x 2 panels = 6 figures):
#   figures/SRK_depletion_ranking_{observed,predicted}_{blank,EOs,all}.{png,pdf}
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) message("Note: install ggrepel for non-overlapping labels.")

source("srk_bl_constants.R")

# Thresholds
DI_THRESHOLD       <- 0.50   # mild vs severe depletion
PCOMP_THRESHOLD    <- 0.40   # random-mating viability

K_SPECIES <- 59              # MM consensus — for subtitle reference

# =============================================================================
# 1. LOAD
# =============================================================================
cat("Loading SRK_EO_allele_richness_with_nulls.tsv ...\n")
dat <- read_tsv("Tables/Phase4/step17b_EO_allele_richness_with_nulls.tsv",
                show_col_types = FALSE) %>%
  mutate(
    level = factor(level, levels = c("EO", "BL")),
    BL    = factor(BL,    levels = BL_ORDER_NUMERIC),
    label = ifelse(level == "EO", paste0("EO", group), group),
    DI_observed  = as.numeric(DI_observed),
    DI_predicted = as.numeric(DI_predicted)
  )

# =============================================================================
# 2. PLOT BUILDER
# =============================================================================
build_panel <- function(di_col, title_tag, subtitle_extra, show = "all") {
  stopifnot(show %in% c("none", "EO", "all"))
  pdat <- dat %>%
    filter(!is.na(.data[[di_col]])) %>%
    mutate(DI = .data[[di_col]],
           P_compat = P_compat_L0,
           P_lo = P_compat_L0_lo,
           P_hi = P_compat_L0_hi)
  # Subset of data to draw points/whiskers/labels for
  pdat_show <- switch(show,
    "none" = pdat[0, ],
    "EO"   = pdat %>% filter(level == "EO"),
    "all"  = pdat
  )

  cat(sprintf("\nStrategy tally for %s:\n", title_tag))
  pdat %>% mutate(
    strategy = case_when(
      DI <  DI_THRESHOLD & P_compat >= PCOMP_THRESHOLD ~ "HEALTHY (random mating)",
      DI >= DI_THRESHOLD & P_compat >= PCOMP_THRESHOLD ~ "Informed breeding + injection (preventive)",
      DI <  DI_THRESHOLD & P_compat <  PCOMP_THRESHOLD ~ "Informed breeding (frequency skew)",
      DI >= DI_THRESHOLD & P_compat <  PCOMP_THRESHOLD ~ "Informed breeding + injection (urgent)",
    )
  ) %>% { print(table(.$level, .$strategy)) }

  ggplot(pdat, aes(x = DI, y = P_compat)) +
    annotate("rect", xmin = 0, xmax = DI_THRESHOLD,
             ymin = PCOMP_THRESHOLD, ymax = 1,
             fill = "#4575b4", alpha = 0.08) +
    annotate("text", x = 0.015, y = 0.985,
             label = "HEALTHY\n(random mating)",
             colour = "#4575b4", fontface = "bold", size = 4,
             hjust = 0, vjust = 1, lineheight = 0.95) +
    annotate("rect", xmin = DI_THRESHOLD, xmax = 1,
             ymin = PCOMP_THRESHOLD, ymax = 1,
             fill = "#fc9272", alpha = 0.14) +
    annotate("text", x = 0.985, y = 0.985,
             label = "INFORMED BREEDING\n+ ALLELE INJECTION\n(preventive)",
             colour = "#a50f15", fontface = "bold", size = 3.8,
             hjust = 1, vjust = 1, lineheight = 0.95) +
    annotate("rect", xmin = 0, xmax = DI_THRESHOLD,
             ymin = 0, ymax = PCOMP_THRESHOLD,
             fill = "#f16913", alpha = 0.20) +
    annotate("text", x = 0.015, y = 0.015,
             label = "INFORMED BREEDING\n(frequency skew)",
             colour = "#8c2d04", fontface = "bold", size = 3.8,
             hjust = 0, vjust = 0, lineheight = 0.95) +
    annotate("rect", xmin = DI_THRESHOLD, xmax = 1,
             ymin = 0, ymax = PCOMP_THRESHOLD,
             fill = "#cb181d", alpha = 0.20) +
    annotate("text", x = 0.985, y = PCOMP_THRESHOLD - 0.015,
             label = "INFORMED BREEDING\n+ ALLELE INJECTION\n(urgent)",
             colour = "#67000d", fontface = "bold", size = 3.8,
             hjust = 1, vjust = 1, lineheight = 0.95) +
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
    scale_fill_manual(values = BL_COLORS, name = "Bottleneck lineage",
                      breaks = BL_ORDER_NUMERIC, drop = FALSE) +
    scale_shape_manual(values = c(EO = 21, BL = 24), name = "Level",
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
    guides(fill  = guide_legend(override.aes = list(shape = 21, size = 4.5,
                                                    alpha = 1),
                                order = 1),
           shape = guide_legend(override.aes = list(size = 4.5,
                                                    fill = "grey80"),
                                order = 2)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    labs(
      title = paste0("Conservation ranking — Compatibility x Depletion (",
                     title_tag, ")"),
      subtitle = paste0(
        "x = Depletion Index (1 - k / k_species; k_species = ", K_SPECIES,
        " from MM consensus).\n",
        "y = Compatible-pair fraction under strict tetraploid SI (L = 0); ",
        "vertical whisker = bootstrap 95% CI.\n",
        subtitle_extra
      ),
      x = "Depletion Index (1 = no alleles, 0 = at species equilibrium)",
      y = "Compatible-pair fraction (fraction of random pairs that produce seed)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", legend.box = "vertical")
}

# =============================================================================
# 3. RENDER
# =============================================================================
dir.create("figures/Phase5", recursive = TRUE, showWarnings = FALSE)

render_set <- function(di_col, title_tag, subtitle_extra, base_stem) {
  variants <- list(blank = "none", EOs = "EO", all = "all")
  for (label in names(variants)) {
    p <- build_panel(di_col, title_tag, subtitle_extra,
                      show = variants[[label]])
    out_stem <- paste0("figures/Phase4/", base_stem, "_", label)
    ggsave(paste0(out_stem, ".png"), plot = p,
           width = 11, height = 8, dpi = 200)
    ggsave(paste0(out_stem, ".pdf"), plot = p,
           width = 11, height = 8)
    cat(sprintf("  Written: %s.{png,pdf}\n", out_stem))
  }
}

render_set(
  "DI_observed", "observed",
  "Observed Depletion Index uses k_rarefied30 (sample-size-corrected).",
  "step17b_depletion_ranking_observed_with_nulls"
)
render_set(
  "DI_predicted", "predicted",
  paste0("Predicted Depletion Index uses MM-extrapolated k per group ",
         "(Step 15 asymptote)."),
  "step17b_depletion_ranking_predicted_with_nulls"
)

cat("\nStep 17 (depletion-ranking) complete.\n")
cat("Variants per panel: _blank (zones only), _EOs (focal EOs only), ",
    "_all (BL + EO).\n", sep = "")
