#!/usr/bin/env Rscript
# =============================================================================
# SRK_GFS_reproductive_effort.R
# Proportion of individuals supporting reproductive effort by GFS at EO level
# =============================================================================
#
# In a self-incompatible (SI) plant, an individual's value as a breeding
# partner scales with the number of distinct allele types it carries:
#   - AAAA (GFS = 0.000): all four genome copies share one allele — produces
#     only one gamete type, contributing no allelic diversity to crosses.
#   - AAAB–ABCD (GFS > 0): at least two distinct alleles present — gametes
#     carry heterozygous allele combinations that can drive compatible crosses.
#
# This script visualises, for each Element Occurrence (EO), the proportion of
# individuals at each GFS tier, highlighting the fraction that "supports
# reproductive effort" (GFS > 0) vs. those that cannot contribute allelic
# diversity to managed crosses.
#
# =============================================================================
# INPUT
#   SRK_individual_GFS.tsv  — per-individual GFS scores (output of Step 18)
#
# OUTPUTS
#   SRK_GFS_reproductive_effort.pdf
#   SRK_GFS_reproductive_effort.png
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# =============================================================================
# SETTINGS
# =============================================================================

GFS_FILE   <- "SRK_individual_GFS.tsv"
EO_FOCUS   <- c("EO25", "EO27", "EO67", "EO70", "EO76")

TP2_PROP_AAAA <- 0.30   # TP2 threshold: >30% AAAA = CRITICAL

GFS_LEVELS  <- c("AAAA (0.000)", "AAAB (0.500)", "AABB (0.667)",
                 "AABC (0.833)", "ABCD (1.000)")
GFS_COLOURS <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(GFS_COLOURS) <- GFS_LEVELS

# Labels for the legend (shorter form for horizontal layout)
GFS_LABELS <- c(
  "AAAA (0.000)" = "AAAA  GFS = 0.000\n(no allelic diversity)",
  "AAAB (0.500)" = "AAAB  GFS = 0.500",
  "AABB (0.667)" = "AABB  GFS = 0.667",
  "AABC (0.833)" = "AABC  GFS = 0.833",
  "ABCD (1.000)" = "ABCD  GFS = 1.000\n(maximum diversity)"
)

# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================
cat("Loading GFS data...\n")
gfs <- read_tsv(GFS_FILE, show_col_types = FALSE)

gfs_eo <- gfs %>%
  filter(EO %in% EO_FOCUS) %>%
  mutate(GFS_class = factor(GFS_class, levels = GFS_LEVELS))

cat(sprintf("  Individuals in focus EOs: %d\n", nrow(gfs_eo)))

# =============================================================================
# 2. EO-LEVEL SUMMARIES FOR ANNOTATION
# =============================================================================
eo_ann <- gfs_eo %>%
  group_by(EO) %>%
  summarise(
    n               = n(),
    n_supporting    = sum(GFS > 0),
    prop_supporting = mean(GFS > 0),
    prop_AAAA       = mean(GFS == 0),
    mean_GFS        = mean(GFS),
    .groups = "drop"
  ) %>%
  # Order worst → best by mean GFS, then reverse so ggplot2 places worst at top
  # (ggplot2 y-axis: first factor level = bottom, last = top)
  arrange(mean_GFS) %>%
  mutate(
    EO      = factor(EO, levels = rev(EO)),   # EO67 first (bottom), EO27 last (top)
    label_x = (prop_AAAA + 1) / 2             # mid-point of the "supporting" segment
  )

# Apply the same EO factor ordering to the individual-level data
gfs_eo <- gfs_eo %>%
  mutate(EO = factor(EO, levels = levels(eo_ann$EO)))

# Pre-compute proportions for reliable stacking order (avoids position_fill
# data-order artefacts with horizontal geom_col)
gfs_prop <- gfs_eo %>%
  count(EO, GFS_class, .drop = FALSE) %>%
  group_by(EO) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(EO, GFS_class)   # explicit factor-level order within each EO

cat("\nEO reproductive support summary:\n")
print(
  eo_ann %>% select(EO, n, n_supporting, prop_supporting, mean_GFS),
  n = Inf
)

# =============================================================================
# 3. FIGURE
# =============================================================================
cat("\nGenerating figure...\n")

# Annotation text placed just beyond x = 1
ann_x <- 1.03

p <- ggplot(gfs_prop, aes(y = EO, x = prop, fill = GFS_class)) +

  # Proportional stacked bars (AAAA at left, ABCD at right)
  geom_col(
    position  = position_stack(reverse = TRUE),
    colour    = "white",
    linewidth = 0.35
  ) +

  # TP2 AAAA threshold line (marks maximum acceptable AAAA proportion)
  geom_vline(
    xintercept = TP2_PROP_AAAA,
    linetype   = "dashed",
    colour     = "grey35",
    linewidth  = 0.7
  ) +
  annotate(
    "text",
    x      = TP2_PROP_AAAA - 0.01,
    y      = 5.6,
    label  = sprintf("TP2 threshold\n(%d%% AAAA)", round(TP2_PROP_AAAA * 100)),
    hjust  = 1, vjust = 1, size = 3.1, colour = "grey35"
  ) +

  # Right-side annotation: proportion, count, and mean GFS
  geom_text(
    data        = eo_ann,
    aes(
      y     = EO,
      x     = ann_x,
      label = sprintf("%.0f%% support\n(n = %d / %d; mean GFS = %.3f)",
                      prop_supporting * 100, n_supporting, n, mean_GFS)
    ),
    inherit.aes = FALSE,
    hjust       = 0,
    size        = 3.2,
    colour      = "grey20",
    lineheight  = 0.95
  ) +

  scale_fill_manual(
    values = GFS_COLOURS,
    labels = GFS_LABELS,
    name   = "Genotype class"
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.40),
    expand = expansion(mult = c(0, 0)),
    breaks = seq(0, 1, 0.25)
  ) +

  labs(
    title    = "Proportion of individuals supporting reproductive effort per Element Occurrence",
    subtitle = paste0(
      "Individuals with GFS > 0 carry 2+ distinct SRK alleles and can contribute allelic diversity to compatible crosses.\n",
      "AAAA individuals (GFS = 0) carry a single allele across all four genome copies and cannot diversify SI gametes.\n",
      "EOs ordered by mean GFS (ascending). Dashed line = TP2 threshold: <=30% AAAA required to avoid CRITICAL status."
    ),
    x = "Proportion of individuals",
    y = "Element Occurrence"
  ) +

  theme_bw(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.key.size    = unit(0.55, "cm"),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 10, face = "bold"),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey30", lineheight = 1.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "bold", size = 12)
  ) +
  guides(
    fill = guide_legend(nrow = 1, label.position = "bottom")
  )

# =============================================================================
# 4. SAVE
# =============================================================================
ggsave("SRK_GFS_reproductive_effort.pdf", plot = p,
       width = 11, height = 6)
cat("  Written: SRK_GFS_reproductive_effort.pdf\n")

ggsave("SRK_GFS_reproductive_effort.png", plot = p,
       width = 11, height = 6, dpi = 200)
cat("  Written: SRK_GFS_reproductive_effort.png\n")

cat("\nDone.\n")
