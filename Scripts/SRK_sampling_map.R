# =============================================================================
# SRK_sampling_map.R — Sampling overview map for the SRK pipeline
# -----------------------------------------------------------------------------
# Adapts EO_BL_geographic_context_map.png (sibling LEPA_EO_spatial_clustering
# project) to highlight, against the full spatial catalogue of 39 populations
# in 19 EOs, which populations have been successfully sampled by the SRK
# pipeline. Purpose: show colleagues what we have done AND what remains.
#
# Visual conventions
#   - EVERY catalogued location is plotted (39 points), at its true
#     lat/lon coords. Fill = parent BL colour (Set1 palette). Size scales
#     with population habitat area (area_ha proxy).
#   - SAMPLED status — for each sampled EO (`BL_status ∈ Assigned, Inferred`
#     in the SRK BL_assignments), we cannot in general resolve a compound
#     Pop code ("27-3", "84-2") to a single location, so we conventionally
#     mark the LARGEST location within each sampled EO as "sampled" and
#     ring it with a dashed black overlay. Labels follow the same rule:
#     the EO code is placed on the largest sampled location.
#   - Focus EOs (EO18, EO25, EO27, EO67, EO70, EO76) — the 6 EOs with N ≥ 5
#     that drive most SRK analyses — receive bold white-fill labels.
#   - JAR region — successful samples whose Pop code has no entry in the
#     spatial catalogue — is summarised as a corner annotation, not plotted.
#
# Inputs (canonical sources):
#   Tables/SRK_individual_allele_genotypes.tsv  — 335 successful samples
#   SRK_individual_BL_assignments.tsv           — Pop → EO → BL resolution
#   EO_location_groups.csv                      — 39 location coords + area
#   Tables/EO_group_BL_summary.csv              — EO → BL (for unsampled EOs)
#
# Outputs:
#   figures/SRK_sampling_map.{png,pdf}  + _presentation variants
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(elevatr)
  library(tidyterra)
  library(rnaturalearth)
  library(terra)
})

source("srk_bl_constants.R")  # BL_COLORS, BL_ORDER

# -----------------------------------------------------------------------------
# 1. Load canonical sources
# -----------------------------------------------------------------------------
geno <- read.csv("Tables/SRK_individual_allele_genotypes.tsv", sep = "\t",
                 fileEncoding = "UTF-8-BOM", check.names = FALSE,
                 stringsAsFactors = FALSE)
successful <- geno[[names(geno)[1]]]
n_total <- length(successful)
cat(sprintf("[geno]  %d successful samples in allele matrix\n", n_total))

bla <- read.csv("SRK_individual_BL_assignments.tsv", sep = "\t",
                fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
bla <- bla[bla$Individual %in% successful, ]
df_in  <- bla[bla$BL_status %in% c("Assigned", "Inferred"), ]
df_jar <- bla[bla$BL_status == "Unassigned", ]
cat(sprintf("[BL]    %d BL-resolved + %d JAR (Unassigned) = %d\n",
            nrow(df_in), nrow(df_jar), nrow(df_in) + nrow(df_jar)))

# -----------------------------------------------------------------------------
# 2. Load location catalogue (lat/lon for every population)
# -----------------------------------------------------------------------------
LOC_PATHS <- c(
  "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility/Brainstorm_NSF/Data/EO_location_groups.csv",
  "Tables/EO_location_groups.csv"
)
loc_path <- LOC_PATHS[file.exists(LOC_PATHS)][1]
if (is.na(loc_path))
  stop("EO_location_groups.csv not found in any of:\n  ",
       paste(LOC_PATHS, collapse = "\n  "))
loc <- read.csv(loc_path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
cat(sprintf("[loc]   %d catalogued populations across %d EOs (%s)\n",
            nrow(loc), length(unique(loc$EOCode)), loc_path))

# -----------------------------------------------------------------------------
# 3. Build EO -> BL map (covers EVERY catalogued EO, even unsampled ones,
#    so the locations of unsampled EOs still get their parent BL colour).
#    Two-pass lookup: EOCode first, then Group as fallback for locations whose
#    EO isn't in the BL summary (e.g., composite entries already split).
# -----------------------------------------------------------------------------
eo_bl_csv <- read.csv("Tables/EO_group_BL_summary.csv",
                      fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
# Composite EO entries like "EO118; EO76" mean one geographic group spanning
# both EOs — split so each contributing EO inherits the group's BL.
split_eos <- function(s) trimws(unlist(strsplit(s, "[;,]")))
eo_bl_long <- do.call(rbind, lapply(seq_len(nrow(eo_bl_csv)), function(i) {
  data.frame(EOCode = split_eos(eo_bl_csv$EO[i]),
             Group  = eo_bl_csv$Group[i],
             BL     = eo_bl_csv$BL[i], stringsAsFactors = FALSE)
}))
eo_bl_long <- unique(eo_bl_long)

loc$BL <- eo_bl_long$BL[match(loc$EOCode, eo_bl_long$EOCode)]
# Fallback: match by Group (each location has a group column matching
# EO_group_BL_summary.csv's Group column)
need <- is.na(loc$BL)
if (any(need)) {
  loc$BL[need] <- eo_bl_long$BL[match(loc$group[need], eo_bl_long$Group)]
}
n_unmapped <- sum(is.na(loc$BL))
if (n_unmapped > 0)
  cat(sprintf("[loc]   %d locations still have no BL after Group fallback\n",
              n_unmapped))

# -----------------------------------------------------------------------------
# 4. Mark sampled locations + per-EO sample counts
#    Convention: for each sampled EO, the location with the largest area_ha
#    is the "marked" sampled location (where the dashed ring + label go).
# -----------------------------------------------------------------------------
loc$is_sampled  <- FALSE
loc$is_focus_eo <- FALSE
FOCUS_EOS <- c("EO18", "EO25", "EO27", "EO67", "EO70", "EO76")

sampled_eos <- unique(df_in$EO)
cat(sprintf("[mark]  %d sampled EOs in spatial catalogue: %s\n",
            length(sampled_eos), paste(sort(sampled_eos), collapse = ", ")))

eo_counts <- as.data.frame(table(df_in$EO), stringsAsFactors = FALSE)
names(eo_counts) <- c("EOCode", "n_samples")

for (eo in sampled_eos) {
  eo_locs <- which(loc$EOCode == eo)
  if (length(eo_locs) == 0) next
  largest <- eo_locs[which.max(loc$area_ha[eo_locs])]
  loc$is_sampled[largest] <- TRUE
  loc$is_focus_eo[largest] <- eo %in% FOCUS_EOS
}

loc$n_samples <- eo_counts$n_samples[match(loc$EOCode, eo_counts$EOCode)]
loc$n_samples[is.na(loc$n_samples)] <- 0L

cat(sprintf("[mark]  %d locations marked sampled (one per sampled EO; the largest by area)\n",
            sum(loc$is_sampled)))

# -----------------------------------------------------------------------------
# 5. Project to UTM 32611 for plotting
# -----------------------------------------------------------------------------
loc_sf <- st_transform(st_as_sf(loc, coords = c("lon", "lat"), crs = 4326), 32611)
xy <- as.data.frame(st_coordinates(loc_sf))
loc$X <- xy$X
loc$Y <- xy$Y

# -----------------------------------------------------------------------------
# 6. Topo backdrop + Snake River + reference cities (light-touch copy)
# -----------------------------------------------------------------------------
city_coords <- data.frame(
  name = c("Boise", "Mountain Home", "Glenns Ferry", "New Plymouth"),
  lat  = c(43.6150, 43.1322, 42.9550, 43.9685),
  lon  = c(-116.2023, -115.6912, -115.3026, -116.8200),
  # Per-city label nudge: Mountain Home and Glenns Ferry nudged south to keep
  # their labels away from BL1 populations (e.g., EO29) sitting just north of
  # the cities. Boise and New Plymouth still nudge north.
  nudge_y_m = c(4000, -5500, -5500, 4000)
)
cities_sf <- st_transform(st_as_sf(city_coords, coords = c("lon", "lat"), crs = 4326), 32611)
city_xy   <- as.data.frame(st_coordinates(cities_sf))
city_coords$X <- city_xy$X
city_coords$Y <- city_xy$Y

bbox_pts <- rbind(data.frame(X = loc$X, Y = loc$Y), city_xy)
xrange   <- range(bbox_pts$X) + c(-15000, 15000)
yrange   <- range(bbox_pts$Y) + c(-15000, 15000)
bbox_utm <- st_as_sfc(st_bbox(c(xmin = xrange[1], xmax = xrange[2],
                                ymin = yrange[1], ymax = yrange[2]), crs = 32611))
bbox_wgs <- st_transform(bbox_utm, 4326)

Sys.setenv(PROJ_LIB = system.file("proj", package = "sf"))

cat("[topo]  fetching DEM (zoom 9)...\n")
dem <- get_elev_raster(locations = st_as_sf(bbox_wgs), z = 9,
                       clip = "locations", verbose = FALSE)
dem_utm <- terra::project(terra::rast(dem), "EPSG:32611")

cat("[topo]  fetching Snake River centreline...\n")
rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines",
                      category = "physical", returnclass = "sf",
                      destdir = tempdir())
snake <- rivers[!is.na(rivers$name) & rivers$name == "Snake", ]
snake_utm <- suppressWarnings(st_crop(st_transform(snake, 32611),
                                      st_bbox(bbox_utm)))

# -----------------------------------------------------------------------------
# 7. JAR annotation
# -----------------------------------------------------------------------------
jar_text <- sprintf(
  "JAR region: %d samples / %d germplasm sub-codes\n(not in spatial catalogue)",
  nrow(df_jar), length(unique(df_jar$Pop))
)

# -----------------------------------------------------------------------------
# 8. Build the figure
# -----------------------------------------------------------------------------
# Lock BL factor order so the legend matches BL_ORDER (area-then-connectivity)
loc$BL_fct <- factor(loc$BL, levels = BL_ORDER)

# Label data: ONLY sampled locations get labels (one per sampled EO, on the
# largest location per the convention above). Use distinct subsets so focus
# EOs draw with the bold white-fill style.
label_df <- loc[loc$is_sampled, c("X", "Y", "EOCode", "n_samples",
                                  "is_focus_eo", "BL")]
label_df$label_text <- sprintf("%s\n(n = %d)", label_df$EOCode, label_df$n_samples)

# Per-EO nudge overrides for labels that ggrepel struggles to place cleanly.
# EO38 sits ~12 km NW of Boise so its label tends to overlap the Boise marker;
# push it north-west to anchor it in the empty quadrant above the river.
label_df$nudge_x <- 0
label_df$nudge_y <- 0
label_df$nudge_x[label_df$EOCode == "EO38"] <-  8000
label_df$nudge_y[label_df$EOCode == "EO38"] <-  9000

p <- ggplot() +
  geom_spatraster(data = dem_utm, maxcell = 5e5) +
  scale_fill_gradientn(
    colours  = c("#fff8e7", "#f0d9a8", "#d4a574", "#a87545", "#6b4423", "#3d2817"),
    name     = "Elevation (m)",
    na.value = NA,
    guide    = guide_colorbar(barheight = unit(3.2, "cm"),
                              barwidth  = unit(0.4, "cm"),
                              order     = 4)
  ) +
  geom_sf(data = snake_utm, aes(color = "Snake River"),
          linewidth = 0.8, alpha = 0.85, show.legend = "line") +
  scale_color_manual(values = c("Snake River" = "#1f77b4"),
                     name   = "Map features",
                     guide  = guide_legend(order = 5,
                                           override.aes = list(linewidth = 1.2))) +
  ggnewscale::new_scale_fill() +

  # --- All 39 catalogued locations ---
  # Unsampled: filled CIRCLE   (shape 21)
  # Sampled:   filled TRIANGLE (shape 24)
  # Both filled by parent BL colour, sized by habitat area, with dark border.
  geom_point(
    data = loc,
    aes(x = X, y = Y, fill = BL_fct, shape = is_sampled, size = area_ha),
    color = "grey15", stroke = 0.55, alpha = 0.95
  ) +
  scale_fill_manual(
    values = BL_COLORS, name = "Bottleneck lineage", drop = FALSE,
    na.value = "grey80",
    guide = guide_legend(
      order = 1,
      override.aes = list(shape = 21, size = 6, color = "grey15",
                          stroke = 0.55, alpha = 0.95)
    )
  ) +
  scale_shape_manual(
    values = c("FALSE" = 21, "TRUE" = 24),
    name   = "SRK status",
    labels = c("FALSE" = "Not yet sampled", "TRUE" = "Sampled"),
    guide  = guide_legend(
      order = 2,
      override.aes = list(fill = "grey80", color = "grey15",
                          stroke = 0.55, size = 6, alpha = 0.95)
    )
  ) +
  scale_size_area(
    name     = "Habitat area (ha)",
    max_size = 11,
    breaks   = c(0.05, 0.5, 1, 5, 10, 20),
    guide    = guide_legend(order = 3,
                            override.aes = list(shape = 21, fill = "grey85",
                                                color = "grey15"))
  ) +

  # --- Reference cities ---
  geom_sf(data = cities_sf, shape = 22, fill = "black",
          color = "white", size = 3.2, stroke = 0.5) +
  ggrepel::geom_text_repel(
    data = city_coords, aes(x = X, y = Y, label = name),
    fontface = "bold", size = 4.8, color = "grey15",
    bg.color = "white", bg.r = 0.15,
    point.padding = unit(0.55, "lines"),
    box.padding = unit(0.75, "lines"),
    nudge_y = city_coords$nudge_y_m,
    min.segment.length = 0.2,
    segment.color = "grey35", segment.size = 0.3,
    max.overlaps = Inf, seed = 42
  ) +

  # --- Focus EO labels (bold white-fill) ---
  ggrepel::geom_label_repel(
    data = label_df[label_df$is_focus_eo, ],
    aes(x = X, y = Y, label = label_text),
    size = 4.0, fontface = "bold", fill = "white", color = "grey10",
    label.size = 0.4, alpha = 0.95,
    box.padding = 0.9, point.padding = 0.7,
    min.segment.length = 0.1, segment.color = "grey35", segment.size = 0.4,
    max.overlaps = Inf, lineheight = 0.9, seed = 7
  ) +

  # --- Other sampled-EO labels (smaller plain text) ---
  # Higher box.padding pushes labels away from city labels (Mountain Home
  # specifically conflicts with EO29 in BL1 north of Mountain Home).
  # Per-row nudge_x/y let ggrepel honour custom offsets (e.g., EO38 -> Boise).
  ggrepel::geom_text_repel(
    data = label_df[!label_df$is_focus_eo, ],
    aes(x = X, y = Y, label = label_text),
    size = 3.0, color = "grey20",
    bg.color = "white", bg.r = 0.12,
    box.padding = 1.1, point.padding = 0.4,
    nudge_x = label_df$nudge_x[!label_df$is_focus_eo],
    nudge_y = label_df$nudge_y[!label_df$is_focus_eo],
    force_pull = 0.15,
    min.segment.length = 0.1, segment.color = "grey55", segment.size = 0.3,
    max.overlaps = Inf, lineheight = 0.9, max.iter = 5000, seed = 13
  ) +

  labs(
    title = sprintf(
      "SRK sampling overview: %d successful samples across %d sampled EOs of %d catalogued",
      n_total, length(sampled_eos), length(unique(loc$EOCode))
    ),
    subtitle = paste0(
      "All ", nrow(loc), " catalogued populations shown, coloured by parent BL and sized by habitat area. ",
      "Sampled populations are drawn as triangles; not-yet-sampled populations are circles. ",
      "Focus EOs used in most analyses (EO18, EO25, EO27, EO67, EO70, EO76) labelled in white-fill boxes; ",
      "other sampled EOs in smaller text."
    ),
    caption = paste0(
      "Population coordinates are not accurately represented: Lepidium papilliferum ",
      "is federally threatened and exact coordinates are confidential. ",
      "Compound germplasm sub-codes can match multiple sub-populations within an EO, ",
      "so the triangle conventionally marks the largest of those candidates."
    ),
    x = NULL, y = NULL
  ) +
  coord_sf(crs = 32611, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13, hjust = 0.5,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 10, color = "grey30", hjust = 0.5,
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = 8, color = "grey45", hjust = 0.5,
                                    margin = margin(t = 10)),
    legend.position  = "right",
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    axis.text        = element_text(size = 7, color = "grey50"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(14, 14, 14, 14)
  )

# -----------------------------------------------------------------------------
# 9. Save (regular + presentation variant)
# -----------------------------------------------------------------------------
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
cairo_pdf("SRK_sampling_map.pdf", width = 13, height = 10); print(p); dev.off()
png("figures/SRK_sampling_map.png", width = 13, height = 10,
    units = "in", res = 300, type = "cairo"); print(p); dev.off()

p_pres <- p + labs(title = NULL, subtitle = NULL, caption = NULL)
cairo_pdf("SRK_sampling_map_presentation.pdf", width = 13, height = 9); print(p_pres); dev.off()
png("figures/SRK_sampling_map_presentation.png", width = 13, height = 9,
    units = "in", res = 300, type = "cairo"); print(p_pres); dev.off()

cat("\nDone.\n")
cat("  SRK_sampling_map.pdf / figures/SRK_sampling_map.png\n")
cat("  SRK_sampling_map_presentation.pdf / figures/SRK_sampling_map_presentation.png\n")
