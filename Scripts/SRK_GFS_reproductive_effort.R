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

# =============================================================================
# 5. AAAA ALLELE IDENTITY FIGURE
# =============================================================================
# Unpacks the AAAA bar in the reproductive effort figure: for each EO, shows
# how many AAAA individuals carry each allele. Two alleles — Allele_050 and
# Allele_057 — are present in AAAA individuals across all five EOs and are
# both members of W-group 1 (likely the same SI specificity). All other alleles
# are grouped as "Other alleles" to keep the figure readable.
# =============================================================================
cat("\nGenerating AAAA allele identity figure...\n")

zygo <- read_tsv(ZYGO_FILE, show_col_types = FALSE)

# Extract each AAAA individual's allele from Allele_composition
# (AAAA individuals carry one allele — Allele_composition is e.g. "Allele_050(2)")
aaaa_raw <- gfs %>%
  filter(Genotype_Pattern == "AAAA", EO %in% EO_FOCUS) %>%
  select(Individual, EO) %>%
  left_join(zygo %>% select(Individual, Allele_composition), by = "Individual") %>%
  mutate(
    Allele = str_extract(Allele_composition, "Allele_\\d+"),
    EO     = factor(EO, levels = levels(eo_ann$EO))   # match existing figure order
  )

# Pan-EO alleles: present in AAAA individuals in every focus EO
pan_eo <- aaaa_raw %>%
  group_by(Allele) %>%
  summarise(n_eos = n_distinct(EO), .groups = "drop") %>%
  filter(n_eos == length(EO_FOCUS)) %>%
  pull(Allele) %>%
  sort()

cat(sprintf("  Pan-EO AAAA alleles (present in all %d EOs): %s\n",
            length(EO_FOCUS), paste(pan_eo, collapse = ", ")))

# Assign display group: pan-EO alleles labelled individually, rest grouped
pan_labels <- setNames(
  paste0(pan_eo, "  (pan-EO, W-group 1)"),
  pan_eo
)
group_levels <- c(pan_labels, "Other alleles")

aaaa_raw <- aaaa_raw %>%
  mutate(
    Allele_group = if_else(
      Allele %in% pan_eo,
      pan_labels[Allele],
      "Other alleles"
    ),
    Allele_group = factor(Allele_group, levels = group_levels)
  )

# Counts per EO × allele group
aaaa_counts <- aaaa_raw %>%
  count(EO, Allele_group, .drop = FALSE)

# Per-EO totals and pan-EO proportion for annotation
aaaa_ann <- aaaa_raw %>%
  group_by(EO) %>%
  summarise(
    n_total  = n(),
    n_pan    = sum(Allele %in% pan_eo),
    pct_pan  = round(100 * mean(Allele %in% pan_eo)),
    .groups  = "drop"
  )

cat("\nAAA allele summary per EO:\n")
print(aaaa_ann, n = Inf)

# Colours: one per pan-EO allele + grey for Other
pan_colours <- c("#e6550d", "#3182bd")[seq_along(pan_eo)]
names(pan_colours) <- pan_labels
allele_colours <- c(pan_colours, "Other alleles" = "#bdbdbd")

ann_x <- max(aaaa_ann$n_total) * 1.03

p_alleles <- ggplot(aaaa_counts, aes(y = EO, x = n, fill = Allele_group)) +
  geom_col(
    position  = position_stack(reverse = TRUE),
    colour    = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data        = aaaa_ann,
    aes(
      y     = EO,
      x     = ann_x,
      label = sprintf("n = %d  (%d%% pan-EO)", n_total, pct_pan)
    ),
    inherit.aes = FALSE,
    hjust       = 0,
    size        = 3.2,
    colour      = "grey20"
  ) +
  scale_fill_manual(values = allele_colours, name = "Allele") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.30)),
    breaks = breaks_pretty()
  ) +
  labs(
    title    = "Allele identity of AAAA individuals per Element Occurrence",
    subtitle = paste0(
      "Each AAAA individual carries a single SRK allele across all four genome copies. ",
      paste(pan_eo, collapse = " and "),
      " are present in AAAA\n",
      "individuals across all five EOs and are both members of W-group 1 ",
      "(HV-identical - likely the same SI specificity).\n",
      "Annotation: total AAAA count and percentage carrying a pan-EO allele. ",
      "EOs ordered by mean GFS (ascending), matching Figure 20a."
    ),
    x = "Number of AAAA individuals",
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
  guides(fill = guide_legend(nrow = 1))

ggsave("SRK_GFS_AAAA_allele_composition.pdf", plot = p_alleles,
       width = 11, height = 6)
cat("  Written: SRK_GFS_AAAA_allele_composition.pdf\n")

ggsave("SRK_GFS_AAAA_allele_composition.png", plot = p_alleles,
       width = 11, height = 6, dpi = 200)
cat("  Written: SRK_GFS_AAAA_allele_composition.png\n")

cat("\nDone.\n")
