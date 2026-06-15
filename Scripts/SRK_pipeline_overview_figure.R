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
              "Hypothesis Testing\n& Crossing Design",
              "Per-individual SI Status\n& Forward Simulation"),
  steps   = c("Steps 1-8",
              "Steps 9-12c",
              "Steps 13-21",
              "Steps 22-23",
              "Steps 25-28"),
  bullets = c(
    "- Multi-CANU assembly\n- Coverage chimera filter\n- RACON polishing\n- WhatsHap polyphase\n- 4 phased haplotypes\n  per tetraploid",
    "- Translate & abundance\n  filter functional proteins\n- Distance-based S-allele\n  clustering (N = 55)\n- Tetraploid genotyping\n  (AAAA - ABCD)\n- Step 12c QC gate\n  -> lab redo CSV",
    "- BL integration\n  (5 lineages BL1-BL5)\n- Accumulation curves\n  (MM = 59, Chao1 = 54)\n- TP1 P_compat x DI\n- GFS + TP2\n- Reproductive effort\n  per EO / BL",
    "- Cross-Brassicaceae HV\n  variability scan\n- UPGMA Class I / II\n- Synonymy networks\n- Cross plan H0/H1a/H1b/\n  H2/H3 (819 attempts)",
    "- Per-individual SI / pSI / SC\n  classification\n- Null-aware genotype\n  rebuild (Step 26)\n- Cascade re-runs\n  (14b, 17b, 19b, 20b)\n- Forward-time inheritance\n  simulator (Step 27)\n- Donor ranking (Step 28)"
  ),
  out_key = c(
    "phased haplotypes\nper barcode",
    "349 functional proteins\n49 allele bins\n367 ingroup genotypes",
    "BL stratification\nP_compat, GFS, TP1/TP2\nlab redo CSV (257 samples)",
    "Cross plan TSVs\n(H0 - H3)\nSynonymy groups",
    "247 SI / 15 pSI / 1 SC /\n234 Insufficient_data\nDonor ranking"
  ),
  stringsAsFactors = FALSE
)

# ---- Geometry --------------------------------------------------------------

n_phases <- nrow(phases)
box_w    <- 1.7
box_h    <- 4.6
gap      <- 0.65
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
             x = x, y = box_h - 1.30,
             label = phases$title[i],
             fontface = "bold", size = 3.2, lineheight = 0.95),

    # Bullet body
    annotate("text",
             x = x - box_w/2 + 0.07,
             y = box_h - 2.45,
             label = phases$bullets[i],
             hjust = 0, vjust = 1, size = 2.55,
             lineheight = 1.08,
             family = "mono"),

    # Output key callout at base
    annotate("rect",
             xmin = x - box_w/2 + 0.10, xmax = x + box_w/2 - 0.10,
             ymin = 0.08, ymax = 0.95,
             fill = "grey92", colour = col, linewidth = 0.4),
    annotate("text",
             x = x, y = 0.78,
             label = "Output",
             fontface = "italic", size = 2.4, colour = "grey30"),
    annotate("text",
             x = x, y = 0.40,
             label = phases$out_key[i],
             size = 2.5, lineheight = 1.05, fontface = "bold")
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
           label = paste0("Nanopore amplicons -> functional S-alleles -> conservation diagnostics. ",
                          "28 steps, 5 phases. Current dataset: 367 ingroup individuals / 49 alleles."),
           size = 3.0, colour = "grey30", fontface = "italic")

# Boxes
for (i in seq_len(n_phases)) {
  p <- Reduce(`+`, draw_phase_box(i), accumulate = FALSE, init = p)
}
# Arrows between boxes
for (i in seq_len(n_phases - 1)) {
  p <- p + draw_arrow_between(i)
}

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
       width = 14, height = 5.5, units = "in")
ggsave("figures/SRK_pipeline_overview.png", plot = p,
       width = 14, height = 5.5, units = "in", dpi = 200, bg = "white")

cat("Written:\n",
    "  figures/SRK_pipeline_overview.pdf\n",
    "  figures/SRK_pipeline_overview.png\n", sep = "")
