#!/usr/bin/env Rscript
# =============================================================================
# SRK_pipeline_overview_figure.R
# Produces a one-page horizontal flow diagram summarising the 5 phases of the
# SRK Bioinformatics Pipeline. Output: figures/SRK_pipeline_overview.{pdf,png}.
# Embedded at the top of README.md and index.Rmd as the at-a-glance visual.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ---- Phase definitions (the data shown in the figure) ----------------------

phases <- data.frame(
  phase   = 1:5,
  label   = c("Phase 1", "Phase 2", "Phase 3", "Phase 4", "Phase 5"),
  title   = c("Assembly\n& Phasing",
              "Functional Proteins,\nAllele Definition,\nGenotyping",
              "Population Genetics\n& Conservation\nDiagnostics",
              "Per-individual\nSelf-Incompatibility Status\n& Forward Simulation",
              "Hypothesis Testing\n& Crossing Design"),
  steps   = c("Steps 1-8",
              "Steps 9-12c",
              "Steps 13-21",
              "Steps 22-25",
              "Steps 26-27"),
  bullets = c(
    "- Multi-replicate\n  assembly (Canu)\n- Assembly polishing\n  (RACON)\n- Polyploid phasing\n  (WhatsHap)\n- Phased haplotypes\n  per individual",
    "- Translate &\n  abundance-filter\n  functional proteins\n- Distance-based\n  S-allele clustering\n- Tetraploid\n  genotyping\n  (AAAA -> ABCD)\n- Data-quality gate\n  -> lab redo CSV",
    "- Bottleneck Lineage\n  integration\n- Species-level\n  allele pool\n  estimation\n- Random-mating\n  viability +\n  allele depletion\n- Fitness scores +\n  reproductive\n  collapse risk\n- Reproductive effort\n  per population\n  / lineage",
    "- Per-individual\n  classification\n  (full / partial /\n   lost SI)\n- Broken-allele\n  aware genotype\n  rebuild\n- Re-runs of\n  population genetics,\n  mating viability,\n  fitness diagnostics\n- Forward-time\n  inheritance\n  simulator\n- Donor ranking for\n  allele injection",
    "- Cross-species\n  variability scan\n  (Brassica +\n   Arabidopsis)\n- Hierarchical\n  clustering into\n  Class I / II\n- Synonymy networks\n  (functionally\n   redundant alleles)\n- Hypothesis-tested\n  cross plan\n  (uses Phase 4\n   SI parents)"
  ),
  out_key = c(
    "Phased haplotypes\nper individual",
    "Functional proteins,\nallele bins,\nindividual genotypes",
    "Lineage stratification,\nmating + fitness metrics,\nlab re-sequencing CSV",
    "Self-incompatibility\nstatus per individual,\ndonor ranking",
    "Hypothesis-tested\ncross plan,\nfunctional allele groups"
  ),
  stringsAsFactors = FALSE
)

# ---- Geometry --------------------------------------------------------------

n_phases <- nrow(phases)
box_w    <- 2.3
box_h    <- 8.5
gap      <- 0.80
total_w  <- n_phases * box_w + (n_phases - 1) * gap

# x-centres for each phase box
xc <- seq(box_w/2, total_w - box_w/2, by = box_w + gap)

# colour palette — distinct per phase, picked to read on screen + print
phase_cols <- c("#377EB8",  # blue   Phase 1
                "#4DAF4A",  # green  Phase 2
                "#FF7F00",  # orange Phase 3
                "#E41A1C",  # red    Phase 4
                "#984EA3")  # purple Phase 5

# ---- Geom layers -----------------------------------------------------------

draw_phase_box <- function(i) {
  x  <- xc[i]
  col <- phase_cols[i]

  list(
    # Header band (coloured)
    annotate("rect",
             xmin = x - box_w/2, xmax = x + box_w/2,
             ymin = box_h - 0.85, ymax = box_h,
             fill = col, colour = NA),

    # Header text — phase label + step range, in white
    annotate("text",
             x = x, y = box_h - 0.30,
             label = phases$label[i],
             colour = "white", fontface = "bold", size = 4.8),
    annotate("text",
             x = x, y = box_h - 0.62,
             label = phases$steps[i],
             colour = "white", size = 3.0),

    # White body
    annotate("rect",
             xmin = x - box_w/2, xmax = x + box_w/2,
             ymin = 0, ymax = box_h - 0.85,
             fill = "white", colour = col, linewidth = 0.7),

    # Title (bold, below header band)
    annotate("text",
             x = x, y = box_h - 1.50,
             label = phases$title[i],
             fontface = "bold", size = 3.7, lineheight = 0.95),

    # Bullet body
    annotate("text",
             x = x - box_w/2 + 0.12,
             y = box_h - 2.80,
             label = phases$bullets[i],
             hjust = 0, vjust = 1, size = 3.05,
             lineheight = 1.18,
             family = "mono"),

    # Output key callout at base
    annotate("rect",
             xmin = x - box_w/2 + 0.10, xmax = x + box_w/2 - 0.10,
             ymin = 0.15, ymax = 1.55,
             fill = "grey92", colour = col, linewidth = 0.4),
    annotate("text",
             x = x, y = 1.32,
             label = "Output",
             fontface = "italic", size = 2.8, colour = "grey30"),
    annotate("text",
             x = x, y = 0.65,
             label = phases$out_key[i],
             size = 2.9, lineheight = 1.15, fontface = "bold")
  )
}

draw_arrow_between <- function(i) {
  x0 <- xc[i]   + box_w / 2
  x1 <- xc[i+1] - box_w / 2
  ymid <- box_h / 2

  annotate("segment",
           x = x0 + 0.05, xend = x1 - 0.05,
           y = ymid, yend = ymid,
           arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
           colour = "grey40", linewidth = 0.7)
}

# ---- Compose plot ----------------------------------------------------------

p <- ggplot() +
  # Title above the boxes
  annotate("text",
           x = total_w / 2, y = box_h + 0.55,
           label = "SRK Bioinformatics Pipeline — 5-Phase Overview",
           fontface = "bold", size = 5.5) +
  annotate("text",
           x = total_w / 2, y = box_h + 0.20,
           label = "Nanopore amplicons -> functional S-alleles -> conservation diagnostics. 28 steps, 5 phases.",
           size = 3.0, colour = "grey30", fontface = "italic")

# Boxes
for (i in seq_len(n_phases)) {
  p <- Reduce(`+`, draw_phase_box(i), accumulate = FALSE, init = p)
}
# Arrows between boxes
for (i in seq_len(n_phases - 1)) {
  p <- p + draw_arrow_between(i)
}

# ─── Chimera filter highlight overlay on Phase 1 (key innovation callout) ───
x1 <- xc[1]
p <- p +
  annotate("rect",
           xmin = x1 - box_w/2 + 0.14, xmax = x1 + box_w/2 - 0.14,
           ymin = 1.75, ymax = 3.25,
           fill = "#FFF0F0", colour = "#E41A1C", linewidth = 0.9) +
  annotate("text",
           x = x1, y = 2.95,
           label = "KEY INNOVATION",
           fontface = "bold", size = 2.4, colour = "#E41A1C") +
  annotate("text",
           x = x1, y = 2.55,
           label = "Detect & discard",
           fontface = "bold", size = 3.1, colour = "#E41A1C") +
  annotate("text",
           x = x1, y = 2.20,
           label = "chimeric sequences",
           fontface = "bold", size = 3.1, colour = "#E41A1C") +
  annotate("text",
           x = x1, y = 1.90,
           label = "(coverage-based filter)",
           fontface = "italic", size = 2.7, colour = "#E41A1C")

# Bottom caption with key external resources
p <- p +
  annotate("text",
           x = 0, y = -0.55,
           label = "Outputs: tables/Phase{2,3,4,5}/stepXX_*.{tsv,csv,fasta}  +  figures/Phase{2,3,4,5}/stepXX_*.{pdf,png}",
           hjust = 0, size = 2.5, colour = "grey25", family = "mono") +
  annotate("text",
           x = 0, y = -0.95,
           label = "External reference FASTAs (Brassica + Arabidopsis SRK) live in FASTA/.  Sample registry at tables/sampling_metadata.csv.",
           hjust = 0, size = 2.5, colour = "grey25") +
  coord_fixed(clip = "off") +
  xlim(-0.2, total_w + 0.2) +
  ylim(-1.2, box_h + 1.0) +
  theme_void() +
  theme(plot.margin = margin(8, 8, 8, 8))

# ---- Save (figures/ root — site-level overview, not Phase-specific) --------

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/SRK_pipeline_overview.pdf", plot = p,
       width = 20, height = 11, units = "in")
ggsave("figures/SRK_pipeline_overview.png", plot = p,
       width = 20, height = 11, units = "in", dpi = 300, bg = "white")

cat("Written:\n",
    "  figures/SRK_pipeline_overview.pdf\n",
    "  figures/SRK_pipeline_overview.png\n", sep = "")
