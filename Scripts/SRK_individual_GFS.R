#!/usr/bin/env Rscript
# =============================================================================
# SRK_individual_GFS.R - Steps 19 + 20 of Canu_amplicon pipeline
# Individual Genotypic Fitness Score (GFS) and Tipping Point 2 (TP2)
# =============================================================================
#
# Computes the Heterozygous Gamete Proportion for each individual - the
# fraction of diploid gametes (sampling 2 of 4 allele copies) that carry
# two distinct alleles. This is a direct measure of individual reproductive
# potential under the self-incompatibility (SI) system.
#
#   GFS_i = 1 - sum_k [ n_k * (n_k - 1) ] / 12
#
# Genotype class mapping:
#   ABCD (1,1,1,1)  -> GFS = 1.000
#   AABC (2,1,1)    -> GFS = 0.833
#   AABB (2,2)      -> GFS = 0.667
#   AAAB (3,1)      -> GFS = 0.500
#   AAAA (4)        -> GFS = 0.000
#
# Two stratification levels (added 2026-05-07):
#   - Element Occurrence (EO; sorted by parent BL on bar plots; circles
#     on the TP2 scatter colored by parent BL)
#   - Bottleneck Lineage (BL1-BL5; aggregate stats; triangles on TP2
#     scatter colored by BL)
#
# Tipping Point 2 (TP2) thresholds applied to per-group summaries:
#   mean_GFS   < 0.667   group average below the AABB level
#   prop_AAAA  > 0.30    more than 30 % of individuals are SI dead-ends
#
# =============================================================================
# INPUTS
#   SRK_individual_zygosity.tsv          - genotype pattern per individual (Step 12)
#   SRK_individual_BL_assignments.tsv    - EO/BL labels per individual (Step 13)
#
# OUTPUTS
#   SRK_individual_GFS.tsv               - per-individual GFS + EO + BL
#   SRK_EO_GFS_summary.tsv               - EO-level summary with TP2 flags + BL column
#   SRK_BL_GFS_summary.tsv               - BL-level summary with TP2 flags (NEW)
#   SRK_GFS_plots.pdf                    - 4-page PDF at root
#   figures/SRK_GFS_plots_p1_composition_proportional.png
#   figures/SRK_GFS_plots_p2_individual_jitter.png
#   figures/SRK_GFS_plots_p3_TP2_tipping_point.png
#   figures/SRK_GFS_plots_p3_TP2_tipping_point_blank.png
#   figures/SRK_GFS_plots_p4_composition_counts.png
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) {
  message("Note: install ggrepel for non-overlapping labels.")
}

# =============================================================================
# USER SETTINGS
# =============================================================================

ZYG_FILE <- "SRK_individual_zygosity.tsv"
BL_FILE  <- "SRK_individual_BL_assignments.tsv"

# Locked Set1 BL palette (matches Steps 14-18 and LEPA_EO_spatial_clustering)
BL_PALETTE <- c(
  BL1 = "#E41A1C", BL2 = "#377EB8", BL3 = "#4DAF4A",
  BL4 = "#984EA3", BL5 = "#FF7F00"
)

# Minimum sample size for an EO to appear in plots (BLs are always plotted)
EO_MIN_N <- 5

# Tipping Point 2 thresholds
TP2_MEAN_GFS  <- 0.667
TP2_PROP_AAAA <- 0.30

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data (Steps 19 + 20)...\n")

zyg     <- read_tsv(ZYG_FILE, show_col_types = FALSE)
bl_data <- read_tsv(BL_FILE,  show_col_types = FALSE)

cat(sprintf("  Genotype file: %d individuals\n", nrow(zyg)))
cat(sprintf("  BL file: %d individuals\n", nrow(bl_data)))

# Restrict to BL-assigned individuals (Assigned + Inferred). The 10
# unresolved germplasm sub-codes are dropped because they cannot be
# placed into an EO/BL.
bl_ok <- bl_data %>%
  filter(BL_status %in% c("Assigned", "Inferred")) %>%
  select(Individual, Pop, EO, BL, Drift_index)

geno <- zyg %>%
  rename(Genotype_Pattern = Genotype) %>%
  mutate(Method = "zygosity_model") %>%
  inner_join(bl_ok, by = "Individual")

cat(sprintf("  BL-assigned individuals retained: %d\n", nrow(geno)))

# =============================================================================
# 2. COMPUTE GFS
# =============================================================================
cat("Computing GFS per individual...\n")

gfs_from_pattern <- function(pattern) {
  dplyr::case_when(
    pattern == "ABCD" ~ 1.000,
    pattern == "AABC" ~ 5 / 6,
    pattern == "AABB" ~ 4 / 6,
    pattern == "AAAB" ~ 3 / 6,
    pattern == "AAAA" ~ 0.000,
    TRUE              ~ NA_real_
  )
}

gfs_class_label <- function(gfs) {
  dplyr::case_when(
    gfs > 0.917 ~ "ABCD (1.000)",
    gfs > 0.750 ~ "AABC (0.833)",
    gfs > 0.583 ~ "AABB (0.667)",
    gfs > 0.250 ~ "AAAB (0.500)",
    TRUE         ~ "AAAA (0.000)"
  )
}

geno <- geno %>%
  mutate(
    GFS       = gfs_from_pattern(Genotype_Pattern),
    GFS_class = gfs_class_label(GFS)
  )

# =============================================================================
# 3. PER-INDIVIDUAL OUTPUT
# =============================================================================
ind_out <- geno %>%
  select(Individual, Pop, EO, BL, Drift_index,
         Genotype_Pattern, GFS, GFS_class, Method)

write_tsv(ind_out, "SRK_individual_GFS.tsv")
cat("  Written: SRK_individual_GFS.tsv\n")

# =============================================================================
# 4. PER-GROUP SUMMARIES (EO and BL)
# =============================================================================
summarise_group <- function(df, key_col) {
  key_sym <- rlang::sym(key_col)
  df %>%
    group_by(!!key_sym) %>%
    summarise(
      n_individuals = n(),
      mean_GFS      = round(mean(GFS), 3),
      sd_GFS        = round(sd(GFS), 3),
      se_GFS        = round(sd(GFS) / sqrt(n()), 3),
      prop_ABCD     = round(mean(GFS > 0.917), 3),
      prop_AABC     = round(mean(GFS > 0.750 & GFS <= 0.917), 3),
      prop_AABB     = round(mean(GFS > 0.583 & GFS <= 0.750), 3),
      prop_AAAB     = round(mean(GFS > 0.250 & GFS <= 0.583), 3),
      prop_AAAA     = round(mean(GFS <= 0.250), 3),
      n_ABCD        = sum(GFS > 0.917),
      n_AABC        = sum(GFS > 0.750 & GFS <= 0.917),
      n_AABB        = sum(GFS > 0.583 & GFS <= 0.750),
      n_AAAB        = sum(GFS > 0.250 & GFS <= 0.583),
      n_AAAA        = sum(GFS <= 0.250),
      .groups = "drop"
    ) %>%
    mutate(
      TP2_mean_breach = mean_GFS  < TP2_MEAN_GFS,
      TP2_AAAA_breach = prop_AAAA > TP2_PROP_AAAA,
      TP2_status = case_when(
        TP2_mean_breach & TP2_AAAA_breach ~ "CRITICAL",
        TP2_mean_breach | TP2_AAAA_breach ~ "AT RISK",
        TRUE                               ~ "OK"
      )
    )
}

# ---- EO summary ----
eo_summary <- summarise_group(geno, "EO")

# Attach parent BL and arrange by BL then mean_GFS
eo_to_bl <- geno %>% distinct(EO, BL)
eo_summary <- eo_summary %>%
  left_join(eo_to_bl, by = "EO") %>%
  arrange(BL, mean_GFS)

write_tsv(eo_summary, "SRK_EO_GFS_summary.tsv")
cat("  Written: SRK_EO_GFS_summary.tsv\n")

cat("\n-- EO-level GFS / TP2 Summary (sorted by BL) --\n")
print(eo_summary %>%
        select(BL, EO, n_individuals, mean_GFS, sd_GFS,
               prop_AAAA, prop_AAAB, prop_AABB, prop_AABC,
               TP2_status),
      n = Inf)

# ---- BL summary ----
bl_summary <- summarise_group(geno, "BL") %>%
  mutate(BL = factor(BL, levels = names(BL_PALETTE))) %>%
  arrange(BL)

write_tsv(bl_summary, "SRK_BL_GFS_summary.tsv")
cat("  Written: SRK_BL_GFS_summary.tsv\n")

cat("\n-- BL-level GFS / TP2 Summary --\n")
print(bl_summary %>%
        select(BL, n_individuals, mean_GFS, sd_GFS,
               prop_AAAA, prop_AAAB, prop_AABB, prop_AABC,
               TP2_status),
      n = Inf)

# =============================================================================
# 5. PLOT SETUP
# =============================================================================
GFS_LEVELS  <- c("AAAA (0.000)", "AAAB (0.500)", "AABB (0.667)",
                 "AABC (0.833)", "ABCD (1.000)")
GFS_COLOURS <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(GFS_COLOURS) <- GFS_LEVELS

# Filter to plotted EOs (N >= EO_MIN_N), sorted by parent BL then EO
eo_levels <- eo_summary %>%
  filter(n_individuals >= EO_MIN_N) %>%
  pull(EO)

geno_plot <- geno %>%
  filter(EO %in% eo_levels) %>%
  mutate(
    EO        = factor(EO, levels = eo_levels),
    GFS_class = factor(GFS_class, levels = GFS_LEVELS),
    BL        = factor(BL, levels = names(BL_PALETTE))
  )

eo_plot <- eo_summary %>%
  filter(EO %in% eo_levels) %>%
  mutate(EO = factor(EO, levels = eo_levels),
         BL = factor(BL, levels = names(BL_PALETTE)))

bl_plot <- bl_summary

# Helper: colored x-axis labels by parent BL
eo_label_colors <- BL_PALETTE[as.character(eo_plot$BL)]

# =============================================================================
# 6. PLOTS
# =============================================================================
cat("\nGenerating plots...\n")

pdf("SRK_GFS_plots.pdf", width = 12, height = 9)

# ---- Plot 1: proportional stacked bar (EO sorted by BL) ----
p1 <- ggplot(geno_plot, aes(x = EO, fill = GFS_class)) +
  geom_bar(position = position_fill(reverse = TRUE),
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = GFS_COLOURS, name = "Genotype class") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  geom_hline(yintercept = 0.30, linetype = "dashed", colour = "grey30",
             linewidth = 0.6) +
  annotate("text", x = length(eo_levels) + 0.4, y = 0.31,
           label = "30% AAAA\n(TP2 threshold)",
           hjust = 1, vjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "Genotype class composition per EO (sorted by parent BL)",
    subtitle = paste0("GFS = Genotypic Fitness Score (heterozygous gamete ",
                      "proportion). Dashed = TP2 AAAA threshold (30%). ",
                      "X-axis labels colored by parent BL (Set1 palette)."),
    x = "Element Occurrence", y = "Proportion of individuals"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        axis.text.x = element_text(colour = eo_label_colors,
                                   face = "bold"))
print(p1)
ggsave("figures/SRK_GFS_plots_p1_composition_proportional.png", p1,
       width = 12, height = 7, dpi = 200)
cat("  Written: figures/SRK_GFS_plots_p1_composition_proportional.png\n")

# ---- Plot 2: individual GFS jitter + mean crossbar (EO sorted by BL) ----
p2 <- ggplot(geno_plot, aes(x = EO, y = GFS, colour = GFS_class)) +
  geom_jitter(width = 0.18, height = 0, size = 2.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55,
               colour = "black", linewidth = 0.8, fatten = 2) +
  geom_hline(yintercept = 0.667, linetype = "dashed",
             colour = "#e41a1c", linewidth = 0.7) +
  geom_hline(yintercept = 0.500, linetype = "dotted",
             colour = "#ff7f00", linewidth = 0.7) +
  scale_colour_manual(values = GFS_COLOURS, name = "Genotype class") +
  scale_y_continuous(
    breaks = c(0, 0.500, 0.667, 0.833, 1.000),
    labels = c("0.000\nAAAA", "0.500\nAAAB", "0.667\nAABB",
               "0.833\nAABC", "1.000\nABCD"),
    limits = c(-0.05, 1.05)
  ) +
  labs(
    title    = "Individual GFS per EO (sorted by parent BL)",
    subtitle = paste0("Bar = EO mean. Red dashed = AABB threshold (0.667); ",
                      "orange dotted = AAAB threshold (0.500). ",
                      "X-axis labels colored by parent BL."),
    x = "Element Occurrence", y = "Genotypic Fitness Score (GFS)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        axis.text.x = element_text(colour = eo_label_colors,
                                   face = "bold"))
print(p2)
ggsave("figures/SRK_GFS_plots_p2_individual_jitter.png", p2,
       width = 12, height = 7, dpi = 200)
cat("  Written: figures/SRK_GFS_plots_p2_individual_jitter.png\n")

# ---- Plot 3: TP2 tipping point scatter (EO + BL on same plot) ----
# Same layout as TP1: zone polygons + threshold lines, EOs as circles +
# BLs as triangles, both colored by BL using the locked Set1 palette.

tp2_eo_pl <- eo_plot %>%
  mutate(label = as.character(EO), level = "EO")
tp2_bl_pl <- bl_plot %>%
  mutate(label = as.character(BL), level = "BL", EO = as.character(BL))
tp2_combined <- bind_rows(
  tp2_eo_pl  %>% mutate(EO = as.character(EO)),
  tp2_bl_pl
) %>%
  mutate(BL = factor(BL, levels = names(BL_PALETTE)),
         level = factor(level, levels = c("EO", "BL")))

p3_base <- ggplot(tp2_combined,
                  aes(x = prop_AAAA, y = mean_GFS)) +
  # OK zone: top-left
  annotate("rect",
           xmin = 0, xmax = TP2_PROP_AAAA,
           ymin = TP2_MEAN_GFS, ymax = 1,
           fill = "#4575b4", alpha = 0.08) +
  annotate("text",
           x = TP2_PROP_AAAA / 2,
           y = (TP2_MEAN_GFS + 1) / 2,
           label = "OK",
           colour = "#4575b4", fontface = "bold", size = 4.5) +
  # AT RISK zone A: top-right
  annotate("rect",
           xmin = TP2_PROP_AAAA, xmax = 1,
           ymin = TP2_MEAN_GFS,  ymax = 1,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = (TP2_PROP_AAAA + 1) / 2,
           y = (TP2_MEAN_GFS + 1) / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # AT RISK zone B: bottom-left
  annotate("rect",
           xmin = 0, xmax = TP2_PROP_AAAA,
           ymin = 0, ymax = TP2_MEAN_GFS,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = TP2_PROP_AAAA / 2,
           y = TP2_MEAN_GFS / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # CRITICAL zone: bottom-right (label nudged toward corner)
  annotate("rect",
           xmin = TP2_PROP_AAAA, xmax = 1,
           ymin = 0, ymax = TP2_MEAN_GFS,
           fill = "#d73027", alpha = 0.08) +
  annotate("text",
           x = TP2_PROP_AAAA + (1 - TP2_PROP_AAAA) * 0.85,
           y = TP2_MEAN_GFS * 0.18,
           label = "CRITICAL",
           colour = "#d73027", fontface = "bold", size = 4) +
  geom_vline(xintercept = TP2_PROP_AAAA, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = TP2_MEAN_GFS, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = 0.05)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.05)) +
  labs(
    title    = "Tipping Point 2 - Genotypic fitness status",
    subtitle = paste0(
      "Circles = EO (N >= ", EO_MIN_N, "), triangles = BL aggregates. ",
      "Colors = parent BL (Set1, matching LEPA_EO_spatial_clustering)."
    ),
    x = "Proportion AAAA individuals",
    y = "Mean GFS"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        legend.box = "vertical")

ggsave("figures/SRK_GFS_plots_p3_TP2_tipping_point_blank.png", p3_base,
       width = 10, height = 8, dpi = 200)
cat("  Written: figures/SRK_GFS_plots_p3_TP2_tipping_point_blank.png\n")

p3 <- p3_base +
  geom_point(data = tp2_combined,
             aes(fill = BL, shape = level),
             size = 5.5, stroke = 0.9, colour = "grey20", alpha = 0.95) +
  scale_fill_manual(values = BL_PALETTE,
                    name = "Bottleneck lineage", drop = FALSE) +
  scale_shape_manual(values = c(EO = 21, BL = 24), name = "Level") +
  guides(fill  = guide_legend(override.aes = list(shape = 21, size = 4.5)),
         shape = guide_legend(override.aes = list(size = 4.5,
                                                  fill = "grey80")))

if (use_repel) {
  p3 <- p3 + ggrepel::geom_text_repel(
    data = tp2_combined,
    aes(label = label,
        fontface = ifelse(level == "BL", "bold.italic", "bold")),
    size = 4.2, show.legend = FALSE,
    box.padding = 0.55, point.padding = 0.4,
    seed = 42, max.overlaps = 20
  )
} else {
  p3 <- p3 + geom_text(data = tp2_combined, aes(label = label),
                       vjust = -1.2, size = 4, fontface = "bold",
                       show.legend = FALSE)
}

print(p3)
ggsave("figures/SRK_GFS_plots_p3_TP2_tipping_point.png", p3,
       width = 11, height = 8, dpi = 200)
cat("  Written: figures/SRK_GFS_plots_p3_TP2_tipping_point.png\n")

# ---- Plot 4: absolute count stacked bar ----
p4 <- ggplot(geno_plot, aes(x = EO, fill = GFS_class)) +
  geom_bar(position = position_stack(reverse = TRUE),
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = GFS_COLOURS, name = "Genotype class") +
  labs(
    title = "Genotype class counts per EO (sorted by parent BL)",
    x = "Element Occurrence",
    y = "Number of individuals"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right",
        axis.text.x = element_text(colour = eo_label_colors,
                                   face = "bold"))
print(p4)
ggsave("figures/SRK_GFS_plots_p4_composition_counts.png", p4,
       width = 12, height = 7, dpi = 200)
cat("  Written: figures/SRK_GFS_plots_p4_composition_counts.png\n")

dev.off()
cat("  Written: SRK_GFS_plots.pdf\n")

# =============================================================================
# 7. CONSOLE SUMMARY
# =============================================================================
cat("\n-- Species-wide GFS summary (BL-assigned subset) --\n")
cat(sprintf("  Total individuals:  %d\n", nrow(geno)))
cat(sprintf("  Overall mean GFS:   %.3f (SD = %.3f)\n",
            mean(geno$GFS), sd(geno$GFS)))
cat("\n  Genotype class breakdown:\n")
geno %>%
  count(GFS_class) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  arrange(GFS_class) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\n-- Groups breaching TP2 thresholds --\n")
cat("\nEO level:\n")
eo_summary %>%
  filter(EO %in% eo_levels, TP2_status != "OK") %>%
  select(BL, EO, mean_GFS, prop_AAAA, TP2_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nBL level:\n")
bl_summary %>%
  filter(TP2_status != "OK") %>%
  select(BL, mean_GFS, prop_AAAA, TP2_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\n-- Seed production priorities (top GFS individuals per EO) --\n")
geno %>%
  filter(EO %in% eo_levels) %>%
  arrange(EO, desc(GFS)) %>%
  group_by(EO) %>%
  slice_head(n = 3) %>%
  select(BL, EO, Individual, Genotype_Pattern, GFS) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nSteps 19 + 20 complete.\n")
