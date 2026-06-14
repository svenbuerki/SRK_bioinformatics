#!/usr/bin/env Rscript
# =============================================================================
# SRK_donor_ranking_figure.R — Step 19 companion
# Allele-richness recovery ladder: best within-BL + best cross-BL donor
# =============================================================================
#
# EOs are the management units we are restoring; BLs are NOT recipients
# here — they enter only as candidate source pools (local-adaptation
# context). Recipients are the 5 focal EOs that the Step 17/18 ranking
# flagged for injection (EO27 is excluded; it sits on the mild side of the
# Depletion Index threshold).
#
# Per recipient (sorted worst-first by current k), a horizontal stacked bar:
#   1) current k(R)                         — solid recipient-BL colour
#   2) + best within-BL donor EO (novel)    — lighter shade of same colour
#   3) + best cross-BL source (additional)  — grey
#
# "Additional" for cross-BL means alleles the cross-BL source brings that the
# within-BL donor did NOT already bring (set difference on novel_allele_ids).
# This makes the bar tell a true sequential story rather than double-counting.
#
# Recipients without a within-BL donor (EO70 alone in BL2; EO76 alone in
# BL3) show only segments 1 and 3.
#
# Reference lines:
#   k_pool    = alleles observed across all groups (= 49)
#   k_species = MM consensus (= 59)
#
# Input:  Tables/SRK_injection_donor_ranking.tsv
# Output: figures/SRK_donor_recovery_ladder.{png,pdf}
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

source("srk_bl_constants.R")

K_SPECIES <- 59
K_POOL    <- 49

dat <- read_tsv("Tables/Phase5/step28_injection_donor_ranking.tsv",
                show_col_types = FALSE)

parse_set <- function(s) {
  if (is.na(s) || s == "-" || s == "") return(character(0))
  strsplit(s, ",", fixed = TRUE)[[1]]
}

label_group <- function(level, group) {
  ifelse(level == "EO", paste0("EO", group), group)
}

tops <- dat %>%
  filter(rank_in_tier == 1) %>%
  mutate(donor_label = label_group(donor_level, donor_group),
         recipient_label = label_group(recipient_level, recipient_group))

recipients <- tops %>%
  filter(recipient_level == "EO") %>%        # EOs only — BLs are source pools
  distinct(recipient_level, recipient_group, recipient_label,
           recipient_BL, recipient_N, k_recipient) %>%
  arrange(k_recipient)

build_row <- function(i) {
  r <- recipients[i, ]
  w <- tops %>% filter(recipient_label == r$recipient_label,
                        tier == "within-BL")
  c <- tops %>% filter(recipient_label == r$recipient_label,
                        tier == "cross-BL")
  w_set <- if (nrow(w) > 0) parse_set(w$novel_allele_ids[1]) else character(0)
  c_set <- if (nrow(c) > 0) parse_set(c$novel_allele_ids[1]) else character(0)
  c_addn <- setdiff(c_set, w_set)
  tibble(
    recipient_label  = r$recipient_label,
    recipient_BL     = r$recipient_BL,
    recipient_N      = r$recipient_N,
    k_current        = r$k_recipient,
    within_donor     = if (nrow(w) > 0) w$donor_label[1] else NA_character_,
    within_novel     = length(w_set),
    cross_donor      = if (nrow(c) > 0) c$donor_label[1] else NA_character_,
    cross_additional = length(c_addn),
    k_after_within   = r$k_recipient + length(w_set),
    k_after_cross    = r$k_recipient + length(w_set) + length(c_addn)
  )
}

rec_df <- bind_rows(lapply(seq_len(nrow(recipients)), build_row))

# Sort: worst (lowest k_current) at TOP of plot. ggplot puts first factor level
# at the bottom, so reverse the factor levels.
rec_df <- rec_df %>%
  mutate(recipient_label = factor(recipient_label,
                                   levels = rev(recipient_label)))

# Long format for stacked bar
long_df <- rec_df %>%
  select(recipient_label, recipient_BL,
         within_donor, cross_donor,
         k_current, within_novel, cross_additional) %>%
  pivot_longer(cols = c(k_current, within_novel, cross_additional),
               names_to = "segment", values_to = "alleles") %>%
  filter(alleles > 0) %>%
  mutate(
    segment = factor(segment,
                     levels = c("k_current", "within_novel", "cross_additional")),
    bl_hex = BL_COLORS[as.character(recipient_BL)],
    fill_color = case_when(
      segment == "k_current"        ~ bl_hex,
      segment == "within_novel"     ~ alpha(bl_hex, 0.45),
      segment == "cross_additional" ~ "grey80"
    ),
    seg_label = case_when(
      segment == "k_current"        ~ paste0("k=", alleles),
      segment == "within_novel"     ~ paste0("+", alleles, " (", within_donor, ")"),
      segment == "cross_additional" ~ paste0("+", alleles, " (", cross_donor, ")")
    )
  )

# Hide labels for very narrow segments so text doesn't overflow
long_df <- long_df %>%
  mutate(seg_label = ifelse(alleles < 3, "", seg_label))

# Legend dummy data: 3 fixed swatches (recipient k, within-BL gain, cross-BL gain)
legend_df <- tibble(
  segment_label = factor(c("Current k(R)",
                            "+ within-BL donor (best)",
                            "+ cross-BL donor (additional, best)"),
                          levels = c("Current k(R)",
                                     "+ within-BL donor (best)",
                                     "+ cross-BL donor (additional, best)")),
  swatch = c("grey40", alpha("grey40", 0.45), "grey80")
)

n_rec <- nrow(recipients)

# Legend block, placed in the empty area between the bars (max ~36 alleles)
# and the k_pool / k_species reference lines.
legend_x       <- 41
legend_y_top   <- n_rec - 0.4
legend_step    <- 0.55
swatch_w       <- 2.2
swatch_h       <- 0.20

p <- ggplot(long_df,
             aes(x = alleles, y = recipient_label, fill = fill_color,
                 group = segment)) +
  geom_col(colour = "grey25", linewidth = 0.35, width = 0.72,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = seg_label),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3.4, colour = "grey15", fontface = "bold") +
  scale_fill_identity() +
  geom_vline(xintercept = K_POOL, linetype = "dashed",
             colour = "grey35", linewidth = 0.5) +
  geom_vline(xintercept = K_SPECIES, linetype = "solid",
             colour = "grey15", linewidth = 0.6) +
  annotate("text", x = K_POOL - 0.4, y = n_rec + 0.55,
           label = paste0("k_pool = ", K_POOL),
           hjust = 1, vjust = 0, size = 3.4, colour = "grey25",
           fontface = "italic") +
  annotate("text", x = K_SPECIES - 0.4, y = n_rec + 0.55,
           label = paste0("k_species = ", K_SPECIES),
           hjust = 1, vjust = 0, size = 3.4, colour = "grey15",
           fontface = "italic") +
  # ---- Legend block ----
  annotate("rect",
           xmin = legend_x - 0.5,
           xmax = legend_x + 17,
           ymin = legend_y_top - 3 * legend_step - 0.25,
           ymax = legend_y_top + 0.25,
           fill = "white", colour = "grey60", linewidth = 0.3) +
  annotate("text", x = legend_x, y = legend_y_top,
           label = "Segments", hjust = 0, vjust = 0.5,
           fontface = "bold", size = 3.6) +
  annotate("rect",
           xmin = legend_x, xmax = legend_x + swatch_w,
           ymin = legend_y_top - 1 * legend_step - swatch_h,
           ymax = legend_y_top - 1 * legend_step + swatch_h,
           fill = "grey40", colour = "grey25") +
  annotate("text", x = legend_x + swatch_w + 0.4,
           y = legend_y_top - 1 * legend_step,
           label = "Current k(R) (BL colour)",
           hjust = 0, vjust = 0.5, size = 3.3) +
  annotate("rect",
           xmin = legend_x, xmax = legend_x + swatch_w,
           ymin = legend_y_top - 2 * legend_step - swatch_h,
           ymax = legend_y_top - 2 * legend_step + swatch_h,
           fill = alpha("grey40", 0.45), colour = "grey25") +
  annotate("text", x = legend_x + swatch_w + 0.4,
           y = legend_y_top - 2 * legend_step,
           label = "+ within-BL donor (best)",
           hjust = 0, vjust = 0.5, size = 3.3) +
  annotate("rect",
           xmin = legend_x, xmax = legend_x + swatch_w,
           ymin = legend_y_top - 3 * legend_step - swatch_h,
           ymax = legend_y_top - 3 * legend_step + swatch_h,
           fill = "grey80", colour = "grey25") +
  annotate("text", x = legend_x + swatch_w + 0.4,
           y = legend_y_top - 3 * legend_step,
           label = "+ cross-BL donor (additional)",
           hjust = 0, vjust = 0.5, size = 3.3) +
  scale_x_continuous(limits = c(0, K_SPECIES + 2),
                     breaks = seq(0, 60, 10),
                     expand = expansion(mult = c(0.005, 0.02))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  labs(
    title = "Allele-richness recovery ladder - best donor per tier",
    subtitle = paste0(
      "Stacked: current k(R)  ->  + best within-BL donor EO (novel)  ->  ",
      "+ best cross-BL donor source (additional, beyond within-BL).\n",
      "Within-BL bar missing = no other focal EO in the same BL ",
      "(EO70 in BL2; EO76 in BL3). Cross-BL donor is the best BL/EO ",
      "source pool outside the recipient's BL."),
    x = "Distinct S-alleles",
    y = NULL,
    caption = paste0(
      "Dashed line: alleles observed across the entire dataset (k_pool). ",
      "Solid line: species-level MM consensus (k_species).")
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.5, lineheight = 1.05,
                                  margin = margin(b = 8)),
    plot.caption = element_text(size = 9, colour = "grey30",
                                 hjust = 0)
  )

dir.create("figures/Phase5", recursive = TRUE, showWarnings = FALSE)
ggsave("figures/Phase5/step28_donor_recovery_ladder.png", p,
       width = 11, height = 7, dpi = 200)
ggsave("figures/Phase5/step28_donor_recovery_ladder.pdf", p,
       width = 11, height = 7)
cat("Written: figures/SRK_donor_recovery_ladder.{png,pdf}\n")
