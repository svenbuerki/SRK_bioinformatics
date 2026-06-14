# =============================================================================
# SRK_inheritance_figures.R — Step 27 figures
#
# Reads Tables/SRK_inheritance_trajectories.tsv + ..._time_to_sc.tsv produced
# by SRK_inheritance_simulator.py and renders three figures:
#
#   figures/SRK_inheritance_pNULL_trajectories.png
#       Per-BL p_NULL trajectories — one line per replicate, scenario in colour.
#       The headline projection figure.
#   figures/SRK_inheritance_SC_progression.png
#       Per-BL SC-frequency trajectories with ribbon (median + IQR across reps).
#   figures/SRK_inheritance_time_to_sc.png
#       Median first-passage time to SC_THRESHOLD per (BL × scenario) as bars.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
})

source("srk_bl_constants.R")

TRAJ <- "Tables/Phase5/step27_inheritance_trajectories.tsv"
TTSC <- "Tables/Phase5/step27_inheritance_time_to_sc.tsv"

stopifnot(file.exists(TRAJ), file.exists(TTSC))

traj <- read_tsv(TRAJ, show_col_types = FALSE) |>
  mutate(
    p_NULL = (n_SI*0 + n_pSI1*1 + n_pSI2*2 + n_pSI3*3 + n_SC*4) / (4 * n_individuals),
    SC_frac = n_SC / n_individuals,
    BL = factor(BL, levels = BL_ORDER),
    scenario = factor(scenario,
                      levels = c("baseline","rescue_low","rescue_high","high_drift"))
  )

ttsc <- read_tsv(TTSC, show_col_types = FALSE) |>
  mutate(BL = factor(BL, levels = BL_ORDER),
         scenario = factor(scenario,
                           levels = c("baseline","rescue_low","rescue_high","high_drift")))

SCEN_COLOURS <- c(
  baseline    = "#666666",
  rescue_low  = "#3182bd",
  rescue_high = "#08519c",
  high_drift  = "#b2182b"
)

dir.create("figures/Phase5", recursive = TRUE, showWarnings = FALSE)

# (1) p_NULL trajectories — one line per replicate per scenario, facet by BL
p1 <- ggplot(traj,
             aes(x = generation, y = p_NULL,
                 colour = scenario, group = interaction(scenario, replicate))) +
  geom_line(alpha = 0.25, linewidth = 0.4) +
  stat_summary(aes(group = scenario), fun = median,
               geom = "line", linewidth = 1.0) +
  facet_wrap(~ BL, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = SCEN_COLOURS) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "Per-copy broken-allele frequency over time, per BL",
    subtitle = paste("Each thin line = 1 replicate; bold line = median across replicates.",
                     "Initial state = empirical Step 26 per-BL genotypes."),
    x = "Generations", y = "Mean p_NULL (broken-copy frequency)",
    colour = "Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("figures/Phase5/step27_inheritance_pNULL_trajectories.png",
       p1, width = 15, height = 5.5, dpi = 200)

# (2) SC fraction progression with ribbon
sc_summary <- traj |>
  group_by(scenario, BL, generation) |>
  summarise(
    median_SC = median(SC_frac),
    lo = quantile(SC_frac, 0.25),
    hi = quantile(SC_frac, 0.75),
    .groups = "drop"
  )

p2 <- ggplot(sc_summary,
             aes(x = generation, y = median_SC,
                 colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
  annotate("text", x = max(sc_summary$generation), y = 0.52,
           label = "SC threshold (50 %)", hjust = 1, size = 3.2,
           colour = "grey30") +
  facet_wrap(~ BL, ncol = 5) +
  scale_colour_manual(values = SCEN_COLOURS) +
  scale_fill_manual(values = SCEN_COLOURS) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "Projected SC-frequency progression per BL",
    subtitle = "Bands = inter-quartile range across replicates; line = median",
    x = "Generations", y = "Proportion of individuals SC",
    colour = "Scenario", fill = "Scenario"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("figures/Phase5/step27_inheritance_SC_progression.png",
       p2, width = 15, height = 5.5, dpi = 200)

# (3) Time-to-SC bars (median + IQR; right-censored shown as ▲ at top)
max_gen <- max(traj$generation)
ttsc_summary <- ttsc |>
  group_by(scenario, BL) |>
  summarise(
    n_rep        = n(),
    n_reached_SC = sum(!is.na(time_to_sc)),
    median_t     = if (any(!is.na(time_to_sc))) median(time_to_sc, na.rm = TRUE)
                   else NA_real_,
    lo           = if (any(!is.na(time_to_sc))) quantile(time_to_sc, 0.25, na.rm = TRUE)
                   else NA_real_,
    hi           = if (any(!is.na(time_to_sc))) quantile(time_to_sc, 0.75, na.rm = TRUE)
                   else NA_real_,
    .groups = "drop"
  )

p3 <- ggplot(ttsc_summary,
             aes(x = BL, y = median_t, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = "white") +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.8), width = 0.2,
                colour = "grey25", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d/%d", n_reached_SC, n_rep),
                y = hi + 5),
            position = position_dodge(width = 0.8), size = 3.2) +
  scale_fill_manual(values = SCEN_COLOURS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Median time to 50 % SC frequency, per BL × scenario",
    subtitle = paste0("Bars = median across replicates that reached SC; ",
                      "labels = #reps reaching SC / #total. Whisker = IQR."),
    x = NULL, y = "Generations to 50 % SC", fill = "Scenario"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "bottom",
        panel.grid.major.x = element_blank())

ggsave("figures/Phase5/step27_inheritance_time_to_sc.png",
       p3, width = 12, height = 6, dpi = 200)

cat("Wrote three figures to figures/SRK_inheritance_*.png\n")
