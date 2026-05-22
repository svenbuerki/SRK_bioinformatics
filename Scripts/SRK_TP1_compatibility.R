#!/usr/bin/env Rscript
# =============================================================================
# SRK_TP1_compatibility.R — revised Step 17 of Canu_amplicon pipeline
# Tipping Point 1 (TP1) — frequency evenness x mating compatibility
# =============================================================================
#
# Replaces the legacy SRK_TP1_tipping_point.R (prop_optimum x evenness).
# Depletion against the species optimum is already shown by the Step 16
# S-allele erosion barplots; TP1 now answers a complementary question:
#
#     Given that a population is depleted, is it still demographically
#     functioning, and if not, what is the conservation intervention?
#
# Axes:
#   x = J (Shannon evenness, H / ln(k))      — NFDS distance
#   y = P_compat (tetraploid sporophytic SI) — mating-pool functionality
#
# Point size encodes rarefied S-allele richness (k at N = 30) so depletion
# is visible as small dots even though the axis is no longer k itself.
#
# Quadrant readings (conservation actions):
#   top-right    Monitor + augment for long-term sustainability
#   top-left     Augment to restore evenness
#   bottom-right Augment urgently — too few compatible mates
#   bottom-left  Biobank + restoration source — collapsed mating pool
#
# Four figures are produced (EO vs BL on separate panels, strict vs leaky SI):
#   EO_strict   EO scatter at L = 0
#   EO_leaky    EO scatter at L = 0.25
#   BL_strict   BL scatter at L = 0
#   BL_leaky    BL scatter at L = 0.25
#
# Inputs:
#   Tables/SRK_EO_allele_richness.tsv  (from srk_eo_allele_richness.py)
#
# Outputs (in figures/):
#   SRK_TP1_compatibility_EO_strict.png    EO scatter, strict SI (L = 0)
#   SRK_TP1_compatibility_EO_leaky.png     EO scatter, leaky SI (L = 0.25)
#   SRK_TP1_compatibility_BL_strict.png    BL scatter, strict SI (L = 0)
#   SRK_TP1_compatibility_BL_leaky.png     BL scatter, leaky SI (L = 0.25)
#   …_blank.png variants for presentations
# Also writes PDF versions at repo root.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) message("Note: install ggrepel for non-overlapping labels.")

source("srk_bl_constants.R")
BL_PALETTE <- BL_COLORS

# Thresholds
TP1_J        <- 0.80   # evenness floor — NFDS-equilibrium expectation
TP1_P_COMPAT <- 0.40   # compatibility floor — functional mating pool

# Species optimum (subtitle reference only; not used in computation).
SPECIES_K <- 58

# =============================================================================
# 1. LOAD
# =============================================================================
cat("Loading SRK_EO_allele_richness.tsv ...\n")

dat <- read_tsv("Tables/SRK_EO_allele_richness.tsv", show_col_types = FALSE) %>%
  mutate(
    level = factor(level, levels = c("EO", "BL")),
    BL    = factor(BL,    levels = BL_ORDER_NUMERIC),
    label = ifelse(level == "EO", paste0("EO", group), group)
  )

cat(sprintf("  %d rows loaded (%d BL + %d EO)\n",
            nrow(dat), sum(dat$level == "BL"), sum(dat$level == "EO")))

# =============================================================================
# 2. PLOT BUILDER
# =============================================================================
build_tp1 <- function(level_filter, p_col, title_suffix, subtitle_extra) {
  pdat <- dat %>%
    filter(level == level_filter) %>%
    mutate(P_compat = .data[[p_col]],
           status = case_when(
             evenness_J <  TP1_J & P_compat <  TP1_P_COMPAT ~ "BIOBANK",
             evenness_J >= TP1_J & P_compat <  TP1_P_COMPAT ~ "AUGMENT_RICHNESS",
             evenness_J <  TP1_J & P_compat >= TP1_P_COMPAT ~ "AUGMENT_EVENNESS",
             TRUE                                            ~ "MONITOR"
           ))

  cat(sprintf("\n%s · %s status tally:\n", level_filter, title_suffix))
  print(table(pdat$status))

  # Zone labels — placed in corners so they never collide with points.
  p_base <- ggplot(pdat, aes(x = evenness_J, y = P_compat)) +
    annotate("rect", xmin = TP1_J, xmax = 1,
             ymin = TP1_P_COMPAT, ymax = 1,
             fill = "#4575b4", alpha = 0.08) +
    annotate("text", x = 0.985, y = 0.985,
             label = "MONITOR\n(augment for sustainability)",
             colour = "#4575b4", fontface = "bold", size = 3.4,
             lineheight = 0.95, hjust = 1, vjust = 1) +
    annotate("rect", xmin = 0, xmax = TP1_J,
             ymin = TP1_P_COMPAT, ymax = 1,
             fill = "#fc8d59", alpha = 0.08) +
    annotate("text", x = 0.015, y = 0.985,
             label = "AUGMENT\n(restore evenness)",
             colour = "#b25a3a", fontface = "bold", size = 3.4,
             lineheight = 0.95, hjust = 0, vjust = 1) +
    annotate("rect", xmin = TP1_J, xmax = 1,
             ymin = 0, ymax = TP1_P_COMPAT,
             fill = "#fc8d59", alpha = 0.08) +
    annotate("text", x = 0.985, y = 0.015,
             label = "AUGMENT URGENTLY\n(too few compatible mates)",
             colour = "#b25a3a", fontface = "bold", size = 3.4,
             lineheight = 0.95, hjust = 1, vjust = 0) +
    annotate("rect", xmin = 0, xmax = TP1_J,
             ymin = 0, ymax = TP1_P_COMPAT,
             fill = "#d73027", alpha = 0.08) +
    annotate("text", x = 0.015, y = 0.015,
             label = "BIOBANK + RESTORE\n(mating pool collapsed)",
             colour = "#d73027", fontface = "bold", size = 3.4,
             lineheight = 0.95, hjust = 0, vjust = 0) +
    geom_vline(xintercept = TP1_J,        linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    geom_hline(yintercept = TP1_P_COMPAT, linetype = "dashed",
               colour = "grey50", linewidth = 0.6) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    labs(
      title = paste0("TP1 (", level_filter, " level) — Evenness x compatibility (",
                     title_suffix, ")"),
      subtitle = paste0(
        if (level_filter == "EO")
          "Circles = EO (N >= 15)."
        else
          "Triangles = BL aggregates (all five lineages).",
        " Colours = parent BL. Point size = rarefied S-allele richness ",
        "(k at N = 30; species optimum ~", SPECIES_K, ").",
        subtitle_extra
      ),
      x = "Frequency evenness  J  (Shannon H / ln(k))",
      y = "Compatible-pair fraction  (tetraploid sporophytic SI)"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", legend.box = "vertical")

  list(base = p_base, dat = pdat)
}

# =============================================================================
# 3. RENDER
# =============================================================================
add_points_and_labels <- function(p_base, pdat, shape_val) {
  p <- p_base +
    geom_point(data = pdat,
               aes(fill = BL, size = k_rarefied30_mean),
               shape = shape_val, stroke = 0.9, colour = "grey20",
               alpha = 0.95) +
    scale_fill_manual(values = BL_PALETTE, name = "Bottleneck lineage",
                      breaks = BL_ORDER_NUMERIC, drop = FALSE) +
    scale_size_continuous(name = "Rarefied k\n(at N = 30)",
                          range = c(3, 10),
                          breaks = c(5, 10, 15, 20)) +
    guides(
      fill = guide_legend(override.aes = list(shape = shape_val, size = 4.5),
                          order = 1),
      size = guide_legend(order = 2)
    )

  if (use_repel) {
    p <- p + ggrepel::geom_text_repel(
      data = pdat,
      aes(label = label),
      fontface = ifelse(pdat$level == "BL", "bold.italic", "bold"),
      size = 4.2, show.legend = FALSE,
      box.padding = 0.55, point.padding = 0.4,
      seed = 42, max.overlaps = 20
    )
  } else {
    p <- p + geom_text(data = pdat, aes(label = label), vjust = -1.2,
                       size = 4, fontface = "bold", show.legend = FALSE)
  }
  p
}

render_panel <- function(level_filter, p_col, title_suffix, subtitle_extra,
                          shape_val, tag) {
  built <- build_tp1(level_filter, p_col, title_suffix, subtitle_extra)
  p <- add_points_and_labels(built$base, built$dat, shape_val)
  base_path <- sprintf("figures/SRK_TP1_compatibility_%s_%s", level_filter, tag)
  ggsave(paste0(base_path, "_blank.png"), plot = built$base,
         width = 10, height = 8, dpi = 200)
  ggsave(paste0(base_path, ".png"),       plot = p,
         width = 10, height = 8, dpi = 200)
  ggsave(sprintf("SRK_TP1_compatibility_%s_%s.pdf", level_filter, tag),
         plot = p, width = 10, height = 8)
  cat(sprintf("  Written: %s{,_blank}.png\n", base_path))
  built$dat
}

# EO panels
eo_strict_dat <- render_panel("EO", "P_compat_L0", "strict SI",
  subtitle_extra = "  Strict SI (L = 0).", shape_val = 21, tag = "strict")
eo_leaky_dat  <- render_panel("EO", "P_compat_L0.25", "with SI leakage L = 0.25",
  subtitle_extra = paste0("  SI leakage L = 0.25 ",
                          "(close to empirical L_hat ~ 0.18)."),
  shape_val = 21, tag = "leaky")

# BL panels
bl_strict_dat <- render_panel("BL", "P_compat_L0", "strict SI",
  subtitle_extra = "  Strict SI (L = 0).", shape_val = 24, tag = "strict")
bl_leaky_dat  <- render_panel("BL", "P_compat_L0.25", "with SI leakage L = 0.25",
  subtitle_extra = paste0("  SI leakage L = 0.25 ",
                          "(close to empirical L_hat ~ 0.18)."),
  shape_val = 24, tag = "leaky")

# =============================================================================
# 4. CONSOLE SUMMARY (each panel)
# =============================================================================
print_panel_summary <- function(d, tag) {
  cat(sprintf("\n--- %s: non-MONITOR groups ---\n", tag))
  d %>%
    filter(status != "MONITOR") %>%
    select(level, group, BL, N, k_rarefied30_mean, evenness_J,
           P_compat, L_hat_from_AAAA, status) %>%
    arrange(status, P_compat) %>%
    { cat(capture.output(print(., n = Inf)), sep = "\n"); . }
}

print_panel_summary(eo_strict_dat, "EO · strict SI (L = 0)")
print_panel_summary(eo_leaky_dat,  "EO · leaky SI (L = 0.25)")
print_panel_summary(bl_strict_dat, "BL · strict SI (L = 0)")
print_panel_summary(bl_leaky_dat,  "BL · leaky SI (L = 0.25)")

cat("\nStep 17 (revised, four-panel) complete.\n")
