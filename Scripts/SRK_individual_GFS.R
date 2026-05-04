#!/usr/bin/env Rscript
# =============================================================================
# SRK_individual_GFS.R
# Individual Genotypic Fitness Score (GFS) for tetraploid SRK genotypes
# =============================================================================
#
# Computes the Heterozygous Gamete Proportion for each individual — the
# fraction of diploid gametes (sampling 2 of 4 allele copies) that carry
# two distinct alleles.  This is a direct measure of individual reproductive
# potential under the self-incompatibility (SI) system and negative
# frequency-dependent selection (NFDS).
#
# Formula:
#   GFS_i = 1 - [ sum_k  n_k * (n_k - 1) ] / 12
#
#   n_k  = copy number of allele k in individual i (sum over all alleles)
#   12   = 2 * C(4,2) — normalisation over the tetraploid gamete space
#
# Genotype class mapping:
#   ABCD  (1,1,1,1)  →  GFS = 1.000   all 6 gametes are heterozygous
#   AABC  (2,1,1)    →  GFS = 0.833   5/6 gametes are heterozygous
#   AABB  (2,2)      →  GFS = 0.667   4/6 gametes are heterozygous
#   AAAB  (3,1)      →  GFS = 0.500   3/6 gametes are heterozygous
#   AAAA  (4)        →  GFS = 0.000   no heterozygous gametes
#
# Key distinction from zygosity analysis: AABB (0.667) > AAAB (0.500)
# even though both have 2 unique alleles — because AABB produces more
# diverse pollen/ovules.
#
# Tipping Point 2 (TP2) thresholds applied to EO-level summaries:
#   mean_GFS < 0.667   EO average below AABB level
#   prop_AAAA > 0.30   more than 30 % of individuals are SI dead-ends
#
# =============================================================================
# INPUTS
#   SRK_individual_zygosity.tsv                   — genotype pattern per individual from Step 12
#   sampling_metadata.csv                          — population metadata (Pop, Ingroup)
#
# OUTPUTS (written to current working directory)
#   SRK_individual_GFS.tsv    — per-individual GFS with genotype class
#   SRK_EO_GFS_summary.tsv   — EO-level summary with TP2 status flags
#   SRK_GFS_plots.pdf         — four diagnostic plots
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ── optional: ggrepel for labelled scatter plot ───────────────────────────────
use_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!use_repel) {
  message("Note: install ggrepel for labelled scatter plot (Plot 3).")
}

# =============================================================================
# USER SETTINGS
# =============================================================================

# Paths (relative to the Canu_amplicon project root)
GENO_FILE <- "SRK_individual_zygosity.tsv"
META_FILE <- "sampling_metadata.csv"

# Element Occurrence → Population mapping.
# Sub-populations of the same EO share the same EO label.
# Extend this vector if additional sub-populations are identified.
EO_MAP <- c(
  "25"    = "EO25",
  "27"    = "EO27", "27pv" = "EO27", "27rt" = "EO27",
  "67"    = "EO67",
  "70"    = "EO70",
  "76"    = "EO76"
)

# Focus EOs for main plots
EO_FOCUS <- c("EO25", "EO27", "EO67", "EO70", "EO76")

# Tipping Point 2 thresholds
TP2_MEAN_GFS  <- 0.667   # mean GFS below AABB level
TP2_PROP_AAAA <- 0.30    # proportion of AAAA individuals

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("Loading data...\n")

zyg  <- read_tsv(GENO_FILE, show_col_types = FALSE)
meta <- read_csv(META_FILE, show_col_types = FALSE)

cat(sprintf("  Genotype file: %d individuals\n", nrow(zyg)))
cat(sprintf("  Metadata file: %d rows\n", nrow(meta)))

# Join Pop and Ingroup from metadata
meta_join <- meta %>%
  select(SampleID, Pop, Ingroup) %>%
  rename(Individual = SampleID) %>%
  mutate(Pop = as.character(Pop))

geno <- zyg %>%
  rename(Genotype_Pattern = Genotype) %>%
  mutate(Method = "zygosity_model") %>%
  left_join(meta_join, by = "Individual")

# =============================================================================
# 2. COMPUTE GFS
# =============================================================================
cat("Computing GFS per individual...\n")

# GFS is computed directly from the called genotype pattern.
# This avoids dependence on partial haplotype recovery (N_total_proteins < 4):
# the zygosity model assigns the pattern from the number of distinct alleles,
# and GFS reflects the expected heterozygous gamete proportion for that pattern.
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
  # Use midpoint thresholds to avoid floating-point exact-equality issues.
  # True GFS values: AAAA=0, AAAB=1/2, AABB=2/3, AABC=5/6, ABCD=1
  # Midpoints:       0.25,  0.583,  0.75,  0.917
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
# 3. ASSIGN EO & FILTER INGROUP
# =============================================================================
geno <- geno %>%
  mutate(
    EO = EO_MAP[Pop],
    EO = if_else(is.na(EO), paste0("Pop_", Pop), EO)
  ) %>%
  filter(Ingroup == 1)

cat(sprintf("  Ingroup individuals retained: %d\n", nrow(geno)))

# =============================================================================
# 4. INDIVIDUAL-LEVEL OUTPUT
# =============================================================================
ind_out <- geno %>%
  select(Individual, Pop, EO, Genotype_Pattern, GFS, GFS_class, Method)

write_tsv(ind_out, "SRK_individual_GFS.tsv")
cat("  Written: SRK_individual_GFS.tsv\n")

# =============================================================================
# 5. EO-LEVEL SUMMARY
# =============================================================================
cat("Computing EO-level summaries...\n")

eo_summary <- geno %>%
  group_by(EO) %>%
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
  ) %>%
  arrange(mean_GFS)

write_tsv(eo_summary, "SRK_EO_GFS_summary.tsv")
cat("  Written: SRK_EO_GFS_summary.tsv\n")

cat("\n── EO-level GFS Summary ────────────────────────────────────────────────\n")
print(
  eo_summary %>%
    select(EO, n_individuals, mean_GFS, sd_GFS,
           prop_AAAA, prop_AAAB, prop_AABB, prop_AABC, prop_ABCD,
           TP2_status),
  n = Inf
)

# =============================================================================
# 6. PLOTS
# =============================================================================
cat("\nGenerating plots...\n")

GFS_LEVELS  <- c("AAAA (0.000)", "AAAB (0.500)", "AABB (0.667)",
                 "AABC (0.833)", "ABCD (1.000)")
GFS_COLOURS <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(GFS_COLOURS) <- GFS_LEVELS

geno_plot <- geno %>%
  filter(EO %in% EO_FOCUS) %>%
  mutate(
    EO        = factor(EO, levels = EO_FOCUS),
    GFS_class = factor(GFS_class, levels = GFS_LEVELS)
  )

eo_plot <- eo_summary %>%
  filter(EO %in% EO_FOCUS) %>%
  mutate(EO = factor(EO, levels = EO_FOCUS))

pdf("SRK_GFS_plots.pdf", width = 12, height = 9)

# ── Plot 1: Proportional stacked bar ─────────────────────────────────────────
p1 <- ggplot(geno_plot, aes(x = EO, fill = GFS_class)) +
  geom_bar(position = position_fill(reverse = TRUE), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = GFS_COLOURS, name = "Genotype class") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  geom_hline(yintercept = 0.30, linetype = "dashed", colour = "grey30",
             linewidth = 0.6) +
  annotate("text", x = 5.4, y = 0.31, label = "30% AAAA\n(TP2 threshold)",
           hjust = 1, vjust = 0, size = 3.2, colour = "grey30") +
  labs(
    title    = "Genotype class composition per Element Occurrence",
    subtitle = paste0("GFS = Genotypic Fitness Score (heterozygous gamete proportion). ",
                      "Dashed = TP2 AAAA threshold (30%)."),
    x = "Element Occurrence",
    y = "Proportion of individuals"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

print(p1)

# ── Plot 2: Individual GFS jitter + mean crossbar ────────────────────────────
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
    title    = "Individual GFS per Element Occurrence",
    subtitle = paste0("Bar = EO mean. ",
                      "Red dashed = AABB threshold (0.667); ",
                      "orange dotted = AAAB threshold (0.500)."),
    x = "Element Occurrence",
    y = "Genotypic Fitness Score (GFS)"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

print(p2)

# ── Plot 3: TP2 tipping point map (prop_AAAA × mean_GFS) ─────────────────────
# Base plot: zone polygons + threshold lines only (no data points)
p3_base <- ggplot(eo_plot,
                  aes(x = prop_AAAA, y = mean_GFS,
                      colour = TP2_status, label = EO)) +
  # --- background zone polygons ---
  # OK zone: top-left (neither threshold breached)
  annotate("rect",
           xmin = 0,            xmax = TP2_PROP_AAAA,
           ymin = TP2_MEAN_GFS, ymax = 1,
           fill = "#4575b4", alpha = 0.08) +
  annotate("text",
           x = TP2_PROP_AAAA / 2,
           y = (TP2_MEAN_GFS + 1) / 2,
           label = "OK",
           colour = "#4575b4", fontface = "bold", size = 4.5) +
  # AT RISK zone A: top-right (AAAA breach only)
  annotate("rect",
           xmin = TP2_PROP_AAAA, xmax = 1,
           ymin = TP2_MEAN_GFS,  ymax = 1,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = (TP2_PROP_AAAA + 1) / 2,
           y = (TP2_MEAN_GFS + 1) / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # AT RISK zone B: bottom-left (mean GFS breach only)
  annotate("rect",
           xmin = 0,            xmax = TP2_PROP_AAAA,
           ymin = 0,            ymax = TP2_MEAN_GFS,
           fill = "#fc8d59", alpha = 0.08) +
  annotate("text",
           x = TP2_PROP_AAAA / 2,
           y = TP2_MEAN_GFS / 2,
           label = "AT RISK",
           colour = "#fc8d59", fontface = "bold", size = 4) +
  # CRITICAL zone: bottom-right (both thresholds breached)
  annotate("rect",
           xmin = TP2_PROP_AAAA, xmax = 1,
           ymin = 0,             ymax = TP2_MEAN_GFS,
           fill = "#d73027", alpha = 0.08) +
  annotate("text",
           x = (TP2_PROP_AAAA + 1) / 2,
           y = TP2_MEAN_GFS / 2,
           label = "CRITICAL",
           colour = "#d73027", fontface = "bold", size = 4) +
  # --- threshold lines and formatting ---
  geom_vline(xintercept = TP2_PROP_AAAA, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = TP2_MEAN_GFS, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  scale_colour_manual(
    values = c("CRITICAL" = "#d73027", "AT RISK" = "#fc8d59", "OK" = "#4575b4"),
    name = "TP2 status"
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = 0.05)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.05)) +
  labs(
    title    = "Tipping Point 2 — EO genotypic fitness status",
    subtitle = paste0("X axis: proportion AAAA individuals (TP2 threshold = ",
                      percent(TP2_PROP_AAAA, accuracy = 1), "). ",
                      "Y axis: mean GFS (TP2 threshold = ", TP2_MEAN_GFS, ")."),
    x = "Proportion AAAA individuals",
    y = "Mean GFS"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("SRK_GFS_plots_p3_TP2_tipping_point_blank.png", plot = p3_base,
       width = 10, height = 8, dpi = 200)
cat("  Written: SRK_GFS_plots_p3_TP2_tipping_point_blank.png\n")

# Full plot: add data points and labels
p3 <- p3_base + geom_point(size = 5.5, alpha = 0.9)

if (use_repel) {
  p3 <- p3 + ggrepel::geom_text_repel(
    size = 4.5, fontface = "bold", show.legend = FALSE,
    box.padding = 0.5, seed = 42
  )
} else {
  p3 <- p3 + geom_text(vjust = -1.2, size = 4, fontface = "bold",
                        show.legend = FALSE)
}

print(p3)

ggsave("SRK_GFS_plots_p3_TP2_tipping_point.png", plot = p3,
       width = 10, height = 8, dpi = 200)
cat("  Written: SRK_GFS_plots_p3_TP2_tipping_point.png\n")

# ── Plot 4: Absolute count stacked bar ───────────────────────────────────────
p4 <- ggplot(geno_plot, aes(x = EO, fill = GFS_class)) +
  geom_bar(position = position_stack(reverse = TRUE), colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = GFS_COLOURS, name = "Genotype class") +
  labs(
    title = "Genotype class counts per Element Occurrence",
    x     = "Element Occurrence",
    y     = "Number of individuals"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

print(p4)

dev.off()
cat("  Written: SRK_GFS_plots.pdf\n")

# =============================================================================
# 7. CONSOLE SUMMARY
# =============================================================================
cat("\n── Species-wide GFS summary ────────────────────────────────────────────\n")
cat(sprintf("  Total ingroup individuals:  %d\n", nrow(geno)))
cat(sprintf("  Overall mean GFS:           %.3f (SD = %.3f)\n",
            mean(geno$GFS), sd(geno$GFS)))
cat("\n  Genotype class breakdown (all ingroup individuals):\n")
geno %>%
  count(GFS_class) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  arrange(GFS_class) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\n── EOs breaching TP2 thresholds ────────────────────────────────────────\n")
eo_summary %>%
  filter(EO %in% EO_FOCUS, TP2_status != "OK") %>%
  select(EO, mean_GFS, prop_AAAA, TP2_status) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\n── Seed production priorities (top individuals within each EO) ─────────\n")
geno %>%
  filter(EO %in% EO_FOCUS) %>%
  arrange(EO, desc(GFS)) %>%
  group_by(EO) %>%
  slice_head(n = 3) %>%
  select(EO, Individual, Pop, Genotype_Pattern, GFS) %>%
  { cat(capture.output(print(., n = Inf)), sep = "\n"); . }

cat("\nDone.\n")
