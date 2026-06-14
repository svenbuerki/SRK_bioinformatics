# =============================================================================
# SRK_SI_status_figures.R — Step 25 figures (7-tier pSI-severity breakdown)
# -----------------------------------------------------------------------------
# Reads Tables/SRK_individual_SI_status.tsv (Step 25 aggregator) and renders
# stacked-bar figures answering "What is the status of the SI system?" at
# three nested scales (species / BL / EO) in TWO variants:
#
#   *_full.png    All 401 ingroup individuals — includes pSI · low confidence
#                 and Insufficient_data tiers. Used for QC / re-do triage.
#   *_robust.png  254 individuals with defensible SI/pSI/SC calls — drops
#                 pSI_low (23) and Insufficient_data (124). Used for the
#                 biological-process figures and the binomial-fit model.
#
# Seven-tier scheme on a green → red gradient:
#   SI                  4 functional copies                #1b7837 (dark green)
#   pSI · 1 NF copy     mild leakiness, 3 functional        #fee0b6 (pale amber)
#   pSI · 2 NF copies   modal pSI, half broken              #e08214 (mid amber)
#   pSI · 3 NF copies   severe, one copy left               #a63603 (deep amber)
#   pSI · low conf      frac_NF < 0.25, ambiguous            #d9d9d9 (light grey)
#   SC                  4 non-functional, robust            #b2182b (dark red)
#   Insufficient_data   < 4 raw haplotypes, redo            #737373 (dark grey)
#
# In-bar text labels are intentionally omitted — they cluttered the stacked
# segments and the side panel (species) + `n=` totals atop (BL, EO) carry the
# numbers cleanly.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(readr)
  library(patchwork)
})

.this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
HERE <- if (!is.null(.this_file)) dirname(normalizePath(.this_file)) else normalizePath(".")
setwd(HERE)

source("srk_bl_constants.R")

IN_TSV   <- "Tables/Phase5/step25b_individual_SI_status.tsv"
FIG_DIR  <- "figures/Phase5"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

STATUS_LEVELS_FULL <- c("SI",
                        "pSI_1nf", "pSI_2nf", "pSI_3nf",
                        "pSI_low",
                        "SC",
                        "Insufficient_data")
STATUS_LEVELS_ROBUST <- c("SI", "pSI_1nf", "pSI_2nf", "pSI_3nf", "SC")

STATUS_COLORS <- c(
  SI                = "#1b7837",
  pSI_1nf           = "#fee0b6",
  pSI_2nf           = "#e08214",
  pSI_3nf           = "#a63603",
  pSI_low           = "#d9d9d9",
  SC                = "#b2182b",
  Insufficient_data = "#737373"
)
STATUS_LABELS <- c(
  SI                = "SI — 4 functional copies",
  pSI_1nf           = "pSI — 1 of 4 copies non-functional",
  pSI_2nf           = "pSI — 2 of 4 copies non-functional",
  pSI_3nf           = "pSI — 3 of 4 copies non-functional",
  pSI_low           = "pSI — low confidence (frac_NF < 0.25)",
  SC                = "SC — 4 non-functional",
  Insufficient_data = "Insufficient data (<4 haps)"
)

si_fill_scale <- function(levels) {
  scale_fill_manual(
    values = STATUS_COLORS[levels],
    labels = STATUS_LABELS[levels],
    name   = NULL, drop = FALSE,
    guide  = guide_legend(nrow = if (length(levels) > 5) 4 else 2,
                          byrow = TRUE,
                          keyheight = unit(0.6, "lines"),
                          keywidth  = unit(1.2, "lines"))
  )
}

# ─── Load + tidy ──────────────────────────────────────────────────────────────
si_raw <- read_tsv(IN_TSV, show_col_types = FALSE) |>
  mutate(
    cp_nf_int = suppressWarnings(as.integer(copies_nonfunctional)),
    SI_status_detailed = case_when(
      SI_status == "pSI" & pSI_confidence == "high" & cp_nf_int == 1 ~ "pSI_1nf",
      SI_status == "pSI" & pSI_confidence == "high" & cp_nf_int == 2 ~ "pSI_2nf",
      SI_status == "pSI" & pSI_confidence == "high" & cp_nf_int == 3 ~ "pSI_3nf",
      SI_status == "pSI" & pSI_confidence == "low"                   ~ "pSI_low",
      TRUE                                                            ~ as.character(SI_status)
    )
  )

cat(sprintf("Loaded %d individuals from %s\n", nrow(si_raw), IN_TSV))

# =============================================================================
# Render block — reused for the FULL and ROBUST variants
# =============================================================================
render_variant <- function(suffix, status_levels, subtitle_tag) {
  si <- si_raw |>
    filter(SI_status_detailed %in% status_levels) |>
    mutate(SI_status = factor(SI_status_detailed, levels = status_levels))

  cat(sprintf("\n--- Rendering %s variant (n = %d) ---\n", suffix, nrow(si)))

  out_sp <- file.path(FIG_DIR, sprintf("step25b_SI_status_species_%s.png", suffix))
  out_bl <- file.path(FIG_DIR, sprintf("step25b_SI_status_by_BL_%s.png",  suffix))
  out_eo <- file.path(FIG_DIR, sprintf("step25b_SI_status_by_EO_%s.png",  suffix))

  fill_sc <- si_fill_scale(status_levels)

  # (1) Species-level — bar + companion side table
  species <- si |>
    count(SI_status, .drop = FALSE) |>
    mutate(pct = 100 * n / sum(n))

  p_sp_bar <- ggplot(species, aes(x = 1, y = pct, fill = SI_status)) +
    geom_col(width = 0.55, colour = "white") +
    scale_fill_manual(values = STATUS_COLORS[status_levels],
                      guide = "none", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      title    = sprintf("SI system status — Lepidium papilliferum (n = %d)", nrow(si)),
      subtitle = subtitle_tag,
      x = NULL, y = "Proportion of individuals (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin        = margin(10, 14, 10, 14)
    )

  side_tbl <- species |>
    mutate(
      label = sprintf("%s\n%d ind. (%.1f%%)",
                      STATUS_LABELS[as.character(SI_status)], n, pct),
      y     = rev(seq_len(n()))
    )
  p_sp_side <- ggplot(side_tbl, aes(x = 1, y = y, fill = SI_status)) +
    geom_tile(width = 0.35, height = 0.85, colour = "white") +
    geom_text(aes(x = 1.4, label = label),
              hjust = 0, vjust = 0.5, size = 4.2, lineheight = 0.95) +
    scale_fill_manual(values = STATUS_COLORS[status_levels], guide = "none") +
    scale_x_continuous(limits = c(0.7, 4.5)) +
    scale_y_continuous(limits = c(0.4, nrow(side_tbl) + 0.6)) +
    theme_void() +
    theme(plot.margin = margin(10, 8, 10, 0))

  p_sp <- p_sp_bar + p_sp_side + plot_layout(widths = c(1, 1.4))
  ggsave(out_sp, p_sp, width = 14, height = 6, dpi = 300)
  cat(sprintf("Wrote %s\n", out_sp))

  # (2) Per-BL — no in-bar counts; `n=` totals atop each bar only
  bl_df <- si |>
    filter(BL_inferred %in% BL_ORDER) |>
    count(BL_inferred, SI_status, .drop = FALSE) |>
    group_by(BL_inferred) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(BL_inferred = factor(BL_inferred, levels = BL_ORDER))

  bl_totals <- bl_df |>
    group_by(BL_inferred) |>
    summarise(total = sum(n), .groups = "drop") |>
    mutate(label = sprintf("n = %d", total))

  p_bl <- ggplot(bl_df, aes(x = BL_inferred, y = pct, fill = SI_status)) +
    geom_col(width = 0.7, colour = "white") +
    geom_text(data = bl_totals,
              aes(x = BL_inferred, y = 106, label = label),
              inherit.aes = FALSE, size = 4.2, fontface = "bold") +
    fill_sc +
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "SI status by bottleneck lineage (BL)",
      subtitle = paste0("BLs ordered by total habitat area (Ne proxy). ",
                        subtitle_tag),
      x = NULL, y = "Proportion of individuals (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x        = element_text(size = 12, face = "bold"),
      legend.position    = "bottom",
      legend.text        = element_text(size = 11),
      plot.margin        = margin(10, 18, 10, 14)
    )

  ggsave(out_bl, p_bl, width = 11, height = 7, dpi = 300)
  cat(sprintf("Wrote %s\n", out_bl))

  # (3) Per-EO — facet by BL, `n=` atop each bar only
  eo_codes <- si |>
    filter(BL_inferred %in% BL_ORDER) |>
    pull(EO_normalised) |>
    unique()
  eo_codes <- eo_codes[nzchar(eo_codes)]
  eo_order <- get_eo_order_within_bl(eo_codes)

  eo_df <- si |>
    filter(EO_normalised %in% eo_order) |>
    count(BL_inferred, EO_normalised, SI_status, .drop = FALSE) |>
    group_by(EO_normalised) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(
      EO_normalised = factor(EO_normalised, levels = eo_order),
      BL_inferred   = factor(BL_inferred,   levels = BL_ORDER)
    )

  eo_totals <- eo_df |>
    group_by(EO_normalised, BL_inferred) |>
    summarise(total = sum(n), .groups = "drop") |>
    mutate(label = sprintf("n=%d", total))

  p_eo <- ggplot(eo_df, aes(x = EO_normalised, y = pct, fill = SI_status)) +
    geom_col(width = 0.78, colour = "white") +
    geom_text(data = eo_totals,
              aes(x = EO_normalised, y = 107, label = label),
              inherit.aes = FALSE, size = 3.4, fontface = "bold") +
    facet_grid(. ~ BL_inferred, scales = "free_x", space = "free_x", switch = "x") +
    fill_sc +
    scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "SI status by population (EO), grouped by bottleneck lineage",
      subtitle = paste0("EOs ordered within each BL by ascending drift index. ",
                        subtitle_tag),
      x = NULL, y = "Proportion of individuals (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
      legend.position    = "bottom",
      legend.text        = element_text(size = 11),
      strip.text         = element_text(face = "bold", size = 12),
      strip.placement    = "outside",
      panel.spacing      = unit(0.6, "lines"),
      plot.margin        = margin(10, 18, 10, 14)
    )

  ggsave(out_eo, p_eo, width = 16, height = 8, dpi = 300)
  cat(sprintf("Wrote %s\n", out_eo))
}

render_variant("full",
               STATUS_LEVELS_FULL,
               "All ingroup individuals (includes pSI low-confidence + Insufficient_data)")
render_variant("robust",
               STATUS_LEVELS_ROBUST,
               "Robust calls only — drops pSI low-confidence + Insufficient_data")

cat("\nDone.\n")
