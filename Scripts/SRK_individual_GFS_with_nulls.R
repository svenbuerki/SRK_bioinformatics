# =============================================================================
# SRK_individual_GFS_with_nulls.R — Step 19b + 20b (null-aware variant)
#
# Recomputes the Genotypic Fitness Score (GFS) and TP2 tipping point per
# individual using the null-aware genotype matrix from Step 26. Null copies
# carry zero functional contribution and the score is rescaled so that
# loss of functional copies (pSI / SC) reduces fitness monotonically.
#
# Null-aware formula (Sven, 2026-06-11)
# -------------------------------------
# Let n_k = count of functional allele k in the 4-slot tetraploid genotype,
# n_func = sum_k n_k = total functional copies (0..4), and Allele_NULL = 4 -
# n_func.
#
#   GFS_func          = 1 - sum_k[ n_k(n_k - 1) ] / (n_func * (n_func - 1))
#                       (functional-allele evenness within the functional pool;
#                       0 when n_func <= 1)
#   GFS_null_aware    = GFS_func * (n_func / 4)
#                       (scaling by the fraction of functional copies — the
#                       broken portion of the tetraploid contributes nothing)
#
# Resulting ordering (lowest to highest):
#   0.000  0000  A000  AA00  AAA0  AAAA
#   0.500  AAAB  AB00  AAB0
#   0.667  AABB
#   0.750  ABC0
#   0.833  AABC
#   1.000  ABCD
#
# Compare to canonical GFS (Step 19): AAAA / AAAB / AABB / AABC / ABCD →
# 0 / 0.5 / 0.667 / 0.833 / 1.0. Null-aware adds a band of intermediate
# genotypes (AB00 = 0.5, ABC0 = 0.75) and reclassifies SC individuals plus
# severely null-loaded genotypes (A000, AA00, AAA0) as zero-fitness.
#
# Outputs
#   SRK_individual_GFS_with_nulls.tsv
#   SRK_EO_GFS_summary_with_nulls.tsv
#   SRK_BL_GFS_summary_with_nulls.tsv
#   figures/SRK_GFS_with_nulls_composition.png
#   figures/SRK_GFS_with_nulls_TP2_scatter.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source("srk_bl_constants.R")

ZYG <- "Tables/Phase4/step23_individual_zygosity_with_nulls.tsv"
GEN <- "Tables/Phase4/step23_individual_allele_genotypes_with_nulls.tsv"

TP2_MEAN_GFS  <- 0.667
TP2_PROP_AAAA <- 0.30   # in null-aware: GFS <= 0.0 (AAAA-equivalent + worse)

cat("Loading null-aware genotype tables\n")
zyg <- read_tsv(ZYG, show_col_types = FALSE)
gen <- read_tsv(GEN, show_col_types = FALSE)

allele_cols <- setdiff(
  grep("^Allele_", names(gen), value = TRUE),
  "Allele_NULL"
)

# Restrict to BL-assigned individuals (BL_inferred ∈ BL_ORDER)
zyg <- zyg %>% filter(BL_inferred %in% BL_ORDER)
gen <- gen %>% filter(BL_inferred %in% BL_ORDER)
cat(sprintf("Retained %d BL-assigned individuals\n", nrow(zyg)))

# Drift index from EO_group_BL_summary
eo_di <- read.csv("Tables/EO_group_BL_summary.csv", stringsAsFactors = FALSE)
parts <- strsplit(eo_di$EO, "[;,]\\s*")
eo_di <- eo_di[rep(seq_len(nrow(eo_di)), lengths(parts)), ]
eo_di$EO <- trimws(unlist(parts))
eo_drift_vec <- tapply(eo_di$Drift_index, eo_di$EO,
                       function(x) round(mean(x, na.rm = TRUE), 3))

# Null-aware GFS per individual (vectorised over the functional allele matrix)
mat_func <- as.matrix(gen[, allele_cols, drop = FALSE])
storage.mode(mat_func) <- "numeric"

gfs_null_aware_row <- function(counts) {
  func_counts <- counts[counts > 0]
  n_func <- sum(func_counts)
  if (n_func <= 1) return(0)
  homozygosity_sum <- sum(func_counts * (func_counts - 1))
  max_homo <- n_func * (n_func - 1)
  gfs_func <- 1 - homozygosity_sum / max_homo
  gfs_func * (n_func / 4)
}

gen$GFS_null_aware <- apply(mat_func, 1, gfs_null_aware_row)

# Attach zygosity / genotype info
ind <- gen %>%
  select(Individual, EO_normalised, BL_inferred, SI_status, pSI_confidence,
         Genotype_class_flag, GFS_null_aware) %>%
  left_join(zyg %>% select(Individual, Genotype, N_null_copies,
                            N_distinct_functional_alleles),
            by = "Individual") %>%
  mutate(Drift_index = as.numeric(eo_drift_vec[EO_normalised]))

# Genotype tier label (for figures)
ind$GFS_tier <- cut(
  ind$GFS_null_aware,
  breaks = c(-Inf, 0.001, 0.501, 0.668, 0.751, 0.834, Inf),
  labels = c("SC / homozygous_null (0.000)", "AAAB / AB00 / AAB0 (0.500)",
             "AABB (0.667)", "ABC0 (0.750)", "AABC (0.833)", "ABCD (1.000)"),
  right = FALSE
)

write_tsv(ind, "Tables/Phase4/step19b_individual_GFS_with_nulls.tsv")
cat("Written SRK_individual_GFS_with_nulls.tsv\n")

# Group summaries
summarise_group <- function(df, key_col) {
  key <- rlang::sym(key_col)
  df %>%
    group_by(!!key) %>%
    summarise(
      n_individuals = n(),
      mean_GFS_null = round(mean(GFS_null_aware), 3),
      sd_GFS_null   = round(sd(GFS_null_aware), 3),
      prop_zero     = round(mean(GFS_null_aware <= 0.001), 3),
      prop_SI       = round(mean(SI_status == "SI" |
                                  (SI_status == "pSI" & pSI_confidence == "low")), 3),
      prop_pSI      = round(mean(SI_status == "pSI" & pSI_confidence == "high"), 3),
      prop_SC       = round(mean(SI_status == "SC"), 3),
      n_SC          = sum(SI_status == "SC"),
      .groups = "drop"
    ) %>%
    mutate(
      TP2_mean_breach = mean_GFS_null  <  TP2_MEAN_GFS,
      TP2_zero_breach = prop_zero      >  TP2_PROP_AAAA,
      TP2_status = case_when(
        TP2_mean_breach & TP2_zero_breach ~ "CRITICAL",
        TP2_mean_breach | TP2_zero_breach ~ "AT RISK",
        TRUE                              ~ "OK"
      )
    )
}

eo_summary <- summarise_group(ind, "EO_normalised") %>%
  left_join(ind %>% distinct(EO_normalised, BL_inferred), by = "EO_normalised") %>%
  select(BL = BL_inferred, EO = EO_normalised, everything()) %>%
  mutate(
    BL = factor(BL, levels = BL_ORDER),
    EO = factor(EO, levels = get_eo_order_within_bl(EO))
  ) %>%
  arrange(BL, EO)

bl_summary <- summarise_group(ind, "BL_inferred") %>%
  rename(BL = BL_inferred) %>%
  mutate(BL = factor(BL, levels = BL_ORDER)) %>%
  arrange(BL)

write_tsv(eo_summary, "Tables/Phase4/step19b_EO_GFS_summary_with_nulls.tsv")
write_tsv(bl_summary, "Tables/Phase4/step19b_BL_GFS_summary_with_nulls.tsv")
cat("Written EO + BL GFS summaries (null-aware)\n\n")

cat("=== BL-LEVEL GFS SUMMARY (null-aware) ===\n")
print(as.data.frame(bl_summary %>%
  select(BL, n_individuals, mean_GFS_null, prop_zero,
         prop_SI, prop_pSI, prop_SC, n_SC, TP2_status)),
  row.names = FALSE)

cat("\n=== EO-LEVEL GFS SUMMARY (null-aware) — sorted by BL ===\n")
print(as.data.frame(eo_summary %>%
  select(BL, EO, n_individuals, mean_GFS_null, prop_zero,
         prop_pSI, prop_SC, n_SC, TP2_status)),
  row.names = FALSE)

# Figures
dir.create("figures/Phase5", recursive = TRUE, showWarnings = FALSE)

# (1) GFS composition by BL
TIER_COLORS <- c(
  "SC / homozygous_null (0.000)" = "#b2182b",
  "AAAB / AB00 / AAB0 (0.500)"   = "#fc8d59",
  "AABB (0.667)"                  = "#fee090",
  "ABC0 (0.750)"                  = "#e0f3f8",
  "AABC (0.833)"                  = "#91bfdb",
  "ABCD (1.000)"                  = "#1b7837"
)

comp_df <- ind %>%
  count(BL_inferred, GFS_tier, .drop = FALSE) %>%
  group_by(BL_inferred) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(BL_inferred = factor(BL_inferred, levels = BL_ORDER))

p_comp <- ggplot(comp_df,
                 aes(x = BL_inferred, y = pct, fill = GFS_tier)) +
  geom_col(width = 0.75, colour = "white") +
  scale_fill_manual(values = TIER_COLORS, name = "GFS tier",
                    guide = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "Null-aware genotype tier composition by BL",
    subtitle = "Including pSI and SC individuals; null copies penalised by formula",
    x = NULL, y = "Proportion of individuals (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        plot.margin        = margin(10, 18, 10, 14))

ggsave("figures/Phase4/step19b_GFS_with_nulls_composition.png",
       p_comp, width = 11, height = 7, dpi = 300)

# (2) TP2 scatter: mean GFS vs prop_zero
scatter_df <- bind_rows(
  eo_summary %>% mutate(level = "EO", label = as.character(EO)) %>%
    select(level, label, BL, n_individuals, mean_GFS_null, prop_zero, TP2_status),
  bl_summary %>% mutate(level = "BL", label = as.character(BL), BL = BL) %>%
    select(level, label, BL, n_individuals, mean_GFS_null, prop_zero, TP2_status)
)

p_tp2 <- ggplot(scatter_df,
                aes(x = prop_zero, y = mean_GFS_null,
                    fill = BL, shape = level, size = n_individuals)) +
  geom_hline(yintercept = TP2_MEAN_GFS, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = TP2_PROP_AAAA, linetype = "dashed", colour = "grey50") +
  geom_point(colour = "black", alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = label), size = 3.2, show.legend = FALSE,
                           max.overlaps = Inf) +
  scale_fill_manual(values = BL_COLORS, name = "BL",
                    breaks = BL_ORDER_NUMERIC) +
  scale_shape_manual(values = c("EO" = 21, "BL" = 24),
                     name = "Level") +
  scale_size_continuous(range = c(2, 7), name = "N individuals") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "TP2 (null-aware) — group-level mean GFS vs proportion zero-GFS",
    subtitle = "Dashed lines: TP2 thresholds (mean GFS < 0.667 OR prop_zero > 30 % → AT RISK)",
    x = "Proportion of individuals with GFS = 0 (AAAA-equivalent or worse)",
    y = "Mean GFS (null-aware)"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

ggsave("figures/Phase4/step20b_TP2_with_nulls_scatter.png",
       p_tp2, width = 10, height = 7, dpi = 300)

cat("\nFigures written to figures/SRK_GFS_with_nulls_{composition,TP2_scatter}.png\n")
cat("\nStep 19b + 20b complete.\n")
