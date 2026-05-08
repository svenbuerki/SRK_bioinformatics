#!/usr/bin/env Rscript
# =============================================================================
# SRK_GFS_reproductive_effort.R   (Step 21 - BL-stratified)
# Proportion of individuals supporting reproductive effort by GFS, at the
# Element Occurrence (EO) and Bottleneck Lineage (BL) levels.
# =============================================================================
#
# In a self-incompatible (SI) plant, an individual's value as a breeding
# partner scales with the number of distinct allele types it carries:
#   - AAAA (GFS = 0.000): all four genome copies share one allele - produces
#     only one gamete type, contributing no allelic diversity to crosses.
#   - AAAB-ABCD (GFS > 0): at least two distinct alleles present - gametes
#     carry heterozygous allele combinations that can drive compatible crosses.
#
# This script visualises, for each EO and each BL, the proportion of
# individuals at each GFS tier, highlighting the fraction that "supports
# reproductive effort" (GFS > 0) vs. those that cannot contribute allelic
# diversity to managed crosses. EOs are the central plotted units; BLs are
# the grouping / test framework. Y-axis labels are coloured by parent BL
# (Set1 palette, matching LEPA_EO_spatial_clustering).
#
# =============================================================================
# INPUTS
#   SRK_individual_GFS.tsv        -- per-individual GFS (Step 19/20 output,
#                                    now carries EO, BL, Drift_index)
#   SRK_individual_zygosity.tsv   -- allele compositions for AAAA panel
#
# OUTPUTS (root)
#   SRK_GFS_reproductive_effort.pdf       -- 2 pages: EO panel + BL panel
#   SRK_GFS_AAAA_allele_composition.pdf   -- 2 pages: EO panel + BL panel
#   SRK_BL_reproductive_effort_summary.tsv (NEW)
#
# OUTPUTS (figures/)
#   SRK_GFS_reproductive_effort_EO.png
#   SRK_GFS_reproductive_effort_BL.png
#   SRK_GFS_AAAA_allele_composition_EO.png
#   SRK_GFS_AAAA_allele_composition_BL.png
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# =============================================================================
# SETTINGS
# =============================================================================

GFS_FILE   <- "SRK_individual_GFS.tsv"
ZYGO_FILE  <- "SRK_individual_zygosity.tsv"
EO_MIN_N   <- 5     # matches Step 17 / Step 19 convention

TP2_PROP_AAAA <- 0.30   # TP2 threshold: >30% AAAA = CRITICAL

# Locked Set1 BL palette (matches Steps 14-20 and LEPA_EO_spatial_clustering)
BL_PALETTE <- c(
  BL1 = "#E41A1C", BL2 = "#377EB8", BL3 = "#4DAF4A",
  BL4 = "#984EA3", BL5 = "#FF7F00"
)

GFS_LEVELS  <- c("AAAA (0.000)", "AAAB (0.500)", "AABB (0.667)",
                 "AABC (0.833)", "ABCD (1.000)")
GFS_COLOURS <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(GFS_COLOURS) <- GFS_LEVELS

GFS_LABELS <- c(
  "AAAA (0.000)" = "AAAA  GFS = 0.000\n(no allelic diversity)",
  "AAAB (0.500)" = "AAAB  GFS = 0.500",
  "AABB (0.667)" = "AABB  GFS = 0.667",
  "AABC (0.833)" = "AABC  GFS = 0.833",
  "ABCD (1.000)" = "ABCD  GFS = 1.000\n(maximum diversity)"
)

dir.create("figures", showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading GFS data...\n")
gfs <- read_tsv(GFS_FILE, show_col_types = FALSE) %>%
  mutate(GFS_class = factor(GFS_class, levels = GFS_LEVELS))

cat(sprintf("  Total individuals: %d\n", nrow(gfs)))
cat(sprintf("  BL-assigned: %d  Unassigned: %d\n",
            sum(!is.na(gfs$BL)), sum(is.na(gfs$BL))))

# =============================================================================
# 2. EO-LEVEL SUBSET (N >= EO_MIN_N) AND BL-LEVEL SUBSET
# =============================================================================

eo_ok <- gfs %>%
  filter(!is.na(EO), !is.na(BL)) %>%
  group_by(EO) %>%
  summarise(n = n(), BL = first(BL), .groups = "drop") %>%
  filter(n >= EO_MIN_N) %>%
  arrange(BL, EO)

EO_FOCUS <- eo_ok$EO
cat(sprintf("\n  EOs with N >= %d: %s\n",
            EO_MIN_N, paste(EO_FOCUS, collapse = ", ")))

gfs_eo <- gfs %>% filter(EO %in% EO_FOCUS)
gfs_bl <- gfs %>% filter(!is.na(BL))

cat(sprintf("  Individuals in focus EOs: %d\n", nrow(gfs_eo)))
cat(sprintf("  Individuals with BL assignment: %d\n", nrow(gfs_bl)))

# =============================================================================
# 3. PER-EO AND PER-BL SUMMARIES
# =============================================================================

eo_ann <- gfs_eo %>%
  group_by(EO) %>%
  summarise(
    BL              = first(BL),
    n               = n(),
    n_supporting    = sum(GFS > 0),
    prop_supporting = mean(GFS > 0),
    prop_AAAA       = mean(GFS == 0),
    mean_GFS        = mean(GFS),
    .groups = "drop"
  ) %>%
  # Order EOs by parent BL, then by mean GFS (worst at top of figure)
  arrange(BL, mean_GFS) %>%
  mutate(
    EO      = factor(EO, levels = rev(EO)),       # worst at top
    BL      = factor(BL, levels = names(BL_PALETTE)),
    label_x = (prop_AAAA + 1) / 2
  )

bl_ann <- gfs_bl %>%
  group_by(BL) %>%
  summarise(
    n               = n(),
    n_supporting    = sum(GFS > 0),
    prop_supporting = mean(GFS > 0),
    prop_AAAA       = mean(GFS == 0),
    mean_GFS        = mean(GFS),
    .groups = "drop"
  ) %>%
  arrange(mean_GFS) %>%
  mutate(
    BL      = factor(BL, levels = rev(BL)),       # worst at top
    label_x = (prop_AAAA + 1) / 2
  )

# Persist BL-level summary (NEW output)
write_tsv(
  bl_ann %>%
    mutate(BL = as.character(BL)) %>%
    select(BL, n, n_supporting, prop_supporting, prop_AAAA, mean_GFS),
  "SRK_BL_reproductive_effort_summary.tsv"
)
cat("  Written: SRK_BL_reproductive_effort_summary.tsv\n")

cat("\n-- EO reproductive support summary --\n")
print(eo_ann %>% select(EO, BL, n, n_supporting, prop_supporting, mean_GFS),
      n = Inf)

cat("\n-- BL reproductive support summary --\n")
print(bl_ann %>% select(BL, n, n_supporting, prop_supporting, mean_GFS),
      n = Inf)

# Apply factor ordering to per-individual data
gfs_eo <- gfs_eo %>%
  mutate(EO = factor(EO, levels = levels(eo_ann$EO)))
gfs_bl <- gfs_bl %>%
  mutate(BL = factor(BL, levels = levels(bl_ann$BL)))

# Pre-compute proportions for stacked bars
eo_prop <- gfs_eo %>%
  count(EO, GFS_class, .drop = FALSE) %>%
  group_by(EO) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(EO, GFS_class)

bl_prop <- gfs_bl %>%
  count(BL, GFS_class, .drop = FALSE) %>%
  group_by(BL) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(BL, GFS_class)

# Y-axis label colours: parent BL for EO panel, BL itself for BL panel
eo_label_colours <- BL_PALETTE[as.character(eo_ann$BL[match(levels(eo_ann$EO),
                                                            as.character(eo_ann$EO))])]
bl_label_colours <- BL_PALETTE[as.character(levels(bl_ann$BL))]

# =============================================================================
# 4. REPRODUCTIVE EFFORT FIGURES (2-page PDF)
# =============================================================================
cat("\nGenerating reproductive effort figures...\n")

ann_x <- 1.03

build_repro_plot <- function(prop_df, ann_df, group_var, label_colours,
                             title, subtitle) {
  ann_df <- ann_df %>%
    mutate(
      .group = .data[[group_var]],
      label  = sprintf("%.0f%% support\n(n = %d / %d; mean GFS = %.3f)",
                       prop_supporting * 100, n_supporting, n, mean_GFS)
    )
  prop_df$.group <- prop_df[[group_var]]

  ggplot(prop_df, aes(y = .group, x = prop, fill = GFS_class)) +
    geom_col(position = position_stack(reverse = TRUE),
             colour = "white", linewidth = 0.35) +
    geom_vline(xintercept = TP2_PROP_AAAA, linetype = "dashed",
               colour = "grey35", linewidth = 0.7) +
    annotate("text",
             x = TP2_PROP_AAAA - 0.01,
             y = nlevels(prop_df$.group) + 0.5,
             label = sprintf("TP2 threshold\n(%d%% AAAA)",
                             round(TP2_PROP_AAAA * 100)),
             hjust = 1, vjust = 1, size = 3.1, colour = "grey35") +
    geom_text(data = ann_df,
              aes(y = .group, x = ann_x, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 3.2, colour = "grey20", lineheight = 0.95) +
    scale_fill_manual(values = GFS_COLOURS, labels = GFS_LABELS,
                      name = "Genotype class") +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 1.40),
                       expand = expansion(mult = c(0, 0)),
                       breaks = seq(0, 1, 0.25)) +
    labs(title = title, subtitle = subtitle,
         x = "Proportion of individuals", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position    = "bottom",
      legend.key.size    = unit(0.55, "cm"),
      legend.text        = element_text(size = 9),
      legend.title       = element_text(size = 10, face = "bold"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, colour = "grey30",
                                        lineheight = 1.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(face = "bold", size = 12,
                                        colour = label_colours)
    ) +
    guides(fill = guide_legend(nrow = 1, label.position = "bottom"))
}

p_eo <- build_repro_plot(
  prop_df       = eo_prop,
  ann_df        = eo_ann,
  group_var     = "EO",
  label_colours = eo_label_colours,
  title         = "Reproductive effort support per Element Occurrence",
  subtitle      = paste0(
    "Individuals with GFS > 0 carry 2+ distinct SRK alleles and can ",
    "contribute allelic diversity to compatible crosses.\n",
    "EOs sorted by parent BL then by mean GFS (worst at top). ",
    "Y-axis labels coloured by parent BL (Set1 palette).\n",
    "Dashed line = TP2 threshold (<= 30% AAAA required to avoid CRITICAL)."
  )
)

p_bl <- build_repro_plot(
  prop_df       = bl_prop,
  ann_df        = bl_ann,
  group_var     = "BL",
  label_colours = bl_label_colours,
  title         = "Reproductive effort support per Bottleneck Lineage",
  subtitle      = paste0(
    "BL aggregates all BL-assigned individuals (n = ", nrow(gfs_bl), "). ",
    "BLs sorted by mean GFS (worst at top).\n",
    "Y-axis labels coloured by BL (Set1 palette, matches ",
    "LEPA_EO_spatial_clustering).\n",
    "Dashed line = TP2 threshold (<= 30% AAAA required to avoid CRITICAL)."
  )
)

pdf("SRK_GFS_reproductive_effort.pdf", width = 11, height = 6)
print(p_eo)
print(p_bl)
invisible(dev.off())
cat("  Written: SRK_GFS_reproductive_effort.pdf\n")

ggsave("figures/SRK_GFS_reproductive_effort_EO.png", plot = p_eo,
       width = 11, height = 6, dpi = 200)
cat("  Written: figures/SRK_GFS_reproductive_effort_EO.png\n")
ggsave("figures/SRK_GFS_reproductive_effort_BL.png", plot = p_bl,
       width = 11, height = 6, dpi = 200)
cat("  Written: figures/SRK_GFS_reproductive_effort_BL.png\n")

# =============================================================================
# 5. AAAA ALLELE IDENTITY FIGURES (2-page PDF)
# =============================================================================
# Unpacks the AAAA bar in the reproductive effort figure: for each EO/BL,
# shows how many AAAA individuals carry each allele. Pan-group alleles are
# present in AAAA individuals across (nearly) every focus group. Allele_050
# and Allele_057 (W-group 1, HV-identical, likely the same SI specificity)
# appear as pan-BL across all 5 BLs and as pan-EO across the 5 large EOs;
# EO18 (smallest, only 4 AAAA individuals) carries a different drift
# signature so a >=80% threshold is used for pan-EO classification.
# =============================================================================

PAN_EO_FRAC <- 0.80   # allele must be in AAAA individuals in >=80% of focus EOs

cat("\nGenerating AAAA allele identity figures...\n")

zygo <- read_tsv(ZYGO_FILE, show_col_types = FALSE)

# Build AAAA-allele table for both EO and BL panels
aaaa_all <- gfs %>%
  filter(Genotype_Pattern == "AAAA") %>%
  select(Individual, EO, BL) %>%
  left_join(zygo %>% select(Individual, Allele_composition),
            by = "Individual") %>%
  mutate(Allele = str_extract(Allele_composition, "Allele_\\d+"))

# ---- EO panel data ----
aaaa_eo <- aaaa_all %>%
  filter(EO %in% EO_FOCUS) %>%
  mutate(EO = factor(EO, levels = levels(eo_ann$EO)))

pan_eo_min   <- ceiling(length(EO_FOCUS) * PAN_EO_FRAC)
pan_eo_table <- aaaa_eo %>%
  group_by(Allele) %>%
  summarise(n_eos = n_distinct(EO), .groups = "drop") %>%
  filter(n_eos >= pan_eo_min) %>%
  arrange(desc(n_eos), Allele)

pan_eo_alleles <- pan_eo_table$Allele
pan_eo_counts  <- setNames(pan_eo_table$n_eos, pan_eo_table$Allele)

cat(sprintf("  Pan-EO AAAA alleles (in >= %d of %d focus EOs, n=%d): %s\n",
            pan_eo_min, length(EO_FOCUS), length(pan_eo_alleles),
            paste(pan_eo_alleles, collapse = ", ")))

# ---- BL panel data ----
aaaa_bl <- aaaa_all %>%
  filter(!is.na(BL)) %>%
  mutate(BL = factor(BL, levels = levels(bl_ann$BL)))

n_bls_total  <- nlevels(aaaa_bl$BL)
pan_bl_table <- aaaa_bl %>%
  group_by(Allele) %>%
  summarise(n_bls = n_distinct(BL), .groups = "drop") %>%
  filter(n_bls == n_bls_total) %>%
  arrange(Allele)

pan_bl_alleles <- pan_bl_table$Allele
pan_bl_counts  <- setNames(pan_bl_table$n_bls, pan_bl_table$Allele)

cat(sprintf("  Pan-BL AAAA alleles (in every BL, n=%d): %s\n",
            length(pan_bl_alleles),
            paste(pan_bl_alleles, collapse = ", ")))

build_aaaa_plot <- function(aaaa_df, group_var, pan_alleles, pan_count_map,
                            n_groups, ann_df, label_colours,
                            title, subtitle_template) {
  if (length(pan_alleles) > 0) {
    pan_labels <- setNames(
      paste0(pan_alleles, "  (",
             pan_count_map[pan_alleles], "/", n_groups, " ", group_var,
             "s, W-group 1)"),
      pan_alleles
    )
  } else {
    pan_labels <- character(0)
  }
  group_levels <- c(pan_labels, "Other alleles")

  aaaa_df <- aaaa_df %>%
    mutate(
      Allele_group = if_else(Allele %in% pan_alleles,
                             pan_labels[Allele],
                             "Other alleles"),
      Allele_group = factor(Allele_group, levels = group_levels)
    )
  aaaa_df$.group <- aaaa_df[[group_var]]

  counts <- aaaa_df %>% count(.group, Allele_group, .drop = FALSE)

  totals <- aaaa_df %>%
    group_by(.group) %>%
    summarise(
      n_total = n(),
      n_pan   = sum(Allele %in% pan_alleles),
      pct_pan = round(100 * mean(Allele %in% pan_alleles)),
      .groups = "drop"
    )

  pan_colours <- c("#e6550d", "#3182bd",
                   "#31a354", "#756bb1")[seq_along(pan_alleles)]
  names(pan_colours) <- pan_labels
  allele_colours <- c(pan_colours, "Other alleles" = "#bdbdbd")

  ann_x <- max(totals$n_total) * 1.03

  ggplot(counts, aes(y = .group, x = n, fill = Allele_group)) +
    geom_col(position = position_stack(reverse = TRUE),
             colour = "white", linewidth = 0.35) +
    geom_text(data = totals,
              aes(y = .group, x = ann_x,
                  label = sprintf("n = %d  (%d%% W-group 1)",
                                  n_total, pct_pan)),
              inherit.aes = FALSE, hjust = 0, size = 3.2, colour = "grey20") +
    scale_fill_manual(values = allele_colours, name = "Allele",
                      drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30)),
                       breaks = breaks_pretty()) +
    labs(title = title, subtitle = subtitle_template,
         x = "Number of AAAA individuals", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position    = "bottom",
      legend.key.size    = unit(0.55, "cm"),
      legend.text        = element_text(size = 9),
      legend.title       = element_text(size = 10, face = "bold"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 9, colour = "grey30",
                                        lineheight = 1.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(face = "bold", size = 12,
                                        colour = label_colours)
    ) +
    guides(fill = guide_legend(nrow = 1))
}

p_aaaa_eo <- build_aaaa_plot(
  aaaa_df        = aaaa_eo,
  group_var      = "EO",
  pan_alleles    = pan_eo_alleles,
  pan_count_map  = pan_eo_counts,
  n_groups       = length(EO_FOCUS),
  label_colours  = eo_label_colours,
  title          = "Allele identity of AAAA individuals per Element Occurrence",
  subtitle_template = paste0(
    "Each AAAA individual carries a single SRK allele across all four ",
    "genome copies. W-group 1 alleles (",
    paste(pan_eo_alleles, collapse = ", "),
    ")\nappear in AAAA individuals in >= ", round(PAN_EO_FRAC * 100),
    "% of focus EOs (HV-identical - likely the same SI specificity ",
    "hammered by drift across multiple EOs).\n",
    "EOs sorted by parent BL then by mean GFS. Y-axis labels coloured by ",
    "parent BL (Set1 palette)."
  )
)

p_aaaa_bl <- build_aaaa_plot(
  aaaa_df        = aaaa_bl,
  group_var      = "BL",
  pan_alleles    = pan_bl_alleles,
  pan_count_map  = pan_bl_counts,
  n_groups       = n_bls_total,
  label_colours  = bl_label_colours,
  title          = "Allele identity of AAAA individuals per Bottleneck Lineage",
  subtitle_template = paste0(
    "AAAA individuals aggregated across all BL-assigned EOs (n = ",
    nrow(aaaa_bl), " AAAA individuals).\n",
    "Pan-BL W-group 1 alleles (",
    if (length(pan_bl_alleles) > 0)
      paste(pan_bl_alleles, collapse = ", ") else "none",
    ") are present in AAAA individuals across every BL - confirms ",
    "shared W-group 1 fixation despite independent bottlenecks.\n",
    "BLs sorted by mean GFS. Y-axis labels coloured by BL ",
    "(Set1 palette, matches LEPA_EO_spatial_clustering)."
  )
)

pdf("SRK_GFS_AAAA_allele_composition.pdf", width = 11, height = 6)
print(p_aaaa_eo)
print(p_aaaa_bl)
invisible(dev.off())
cat("  Written: SRK_GFS_AAAA_allele_composition.pdf\n")

ggsave("figures/SRK_GFS_AAAA_allele_composition_EO.png", plot = p_aaaa_eo,
       width = 11, height = 6, dpi = 200)
cat("  Written: figures/SRK_GFS_AAAA_allele_composition_EO.png\n")
ggsave("figures/SRK_GFS_AAAA_allele_composition_BL.png", plot = p_aaaa_bl,
       width = 11, height = 6, dpi = 200)
cat("  Written: figures/SRK_GFS_AAAA_allele_composition_BL.png\n")

cat("\nDone.\n")
