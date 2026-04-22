# ============================================================
# EO Spatial Connectivity & Clustering — Stage 1 C3 Hypothesis
# Goal: Identify geographic groups of locations based on 500m
#       pollinator-dispersal connectivity to proxy shared vs.
#       independent demographic history (ancestral bottleneck)
# Input:  Peggy_EOs_Germplasm_w_lat_long_from_Events_22Apr2026.csv
# Output: EO_location_groups.csv, EO_connectivity_summary.csv,
#         EO_connectivity_map.pdf/.png, EO_clustering_dendrogram.pdf
# ============================================================

library(sf)
library(igraph)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(maps)
library(cluster)   # silhouette
library(ggnewscale) # dual fill scales in ggplot2

set.seed(42)

out_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) ""
)
if (!nzchar(out_dir) || out_dir == ".") {
  out_dir <- "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility/Brainstorm_NSF/Data"
}

# ============================================================
# 1. LOAD DATA
# ============================================================
dat <- read.csv(
  file.path(out_dir, "Peggy_EOs_Germplasm_w_lat_long_from_Events_22Apr2026.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# Remove rows with missing coordinates
dat <- dat[!is.na(dat$eventDecimalLatitude) & !is.na(dat$eventDecimalLongitude), ]
dat$lat <- as.numeric(dat$eventDecimalLatitude)
dat$lon <- as.numeric(dat$eventDecimalLongitude)
dat <- dat[!is.na(dat$lat) & !is.na(dat$lon), ]

cat(sprintf("Records after removing missing coords: %d\n", nrow(dat)))
cat(sprintf("Unique EOs (EOCode):     %d\n", length(unique(dat$EOCode))))
cat(sprintf("Unique locations (locationID): %d\n\n", length(unique(dat$locationID))))

# ============================================================
# 2. BUILD SF OBJECTS — event points, then convex hulls per location
# ============================================================

# Event points as sf (WGS84)
pts_sf <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)

# For each locationID: collect all event points → convex hull
# st_convex_hull on a MULTIPOINT gives a POLYGON (or POINT/LINESTRING for 1-2 pts)
loc_hulls_wgs <- lapply(split(pts_sf, pts_sf$locationID), function(sub) {
  eo   <- unique(sub$EOCode)[1]
  loc  <- unique(sub$locationID)[1]
  n    <- nrow(sub)
  geom <- st_union(sub)              # MULTIPOINT
  hull <- st_convex_hull(geom)       # POLYGON / LINESTRING / POINT
  st_sf(
    locationID = loc,
    EOCode     = eo,
    n_events   = n,
    geometry   = hull
  )
})
loc_hulls_wgs <- do.call(rbind, loc_hulls_wgs)
rownames(loc_hulls_wgs) <- NULL

cat("=== Location hull geometry types ===\n")
print(table(st_geometry_type(loc_hulls_wgs)))
cat("\n")

# Project to UTM Zone 11N (EPSG:32611) for meter-accurate distances
loc_hulls <- st_transform(loc_hulls_wgs, crs = 32611)

# Also keep event centroids per location for labelling
loc_centroids <- aggregate(
  cbind(lat = dat$lat, lon = dat$lon) ~ locationID + EOCode,
  data = dat,
  FUN  = mean
)
event_counts <- table(dat$locationID)
loc_centroids$n_events <- as.integer(
  event_counts[as.character(loc_centroids$locationID)]
)
loc_centroids$n_events[is.na(loc_centroids$n_events)] <- 1L

# ============================================================
# 3. PAIRWISE MINIMUM DISTANCES BETWEEN LOCATION HULLS (metres)
# ============================================================
n_loc    <- nrow(loc_hulls)
loc_ids  <- loc_hulls$locationID

dist_mat <- st_distance(loc_hulls)           # units: metres (UTM)
dist_mat <- matrix(as.numeric(dist_mat),
                   nrow = n_loc, ncol = n_loc,
                   dimnames = list(loc_ids, loc_ids))

cat("=== Pairwise hull-to-hull distance range (m) ===\n")
off_diag <- dist_mat[dist_mat > 0]
cat(sprintf("Min (off-diagonal): %.0f m\n", min(off_diag)))
cat(sprintf("Max:                %.0f m\n\n", max(dist_mat)))

# ============================================================
# 4. 500m CONNECTIVITY → ADJACENCY MATRIX → CONNECTED COMPONENTS
# ============================================================
THRESHOLD_M <- 500   # pollinator dispersal limit

# Exclude diagonal (self-loops) but include dist=0 pairs (overlapping hulls = connected)
adj <- dist_mat <= THRESHOLD_M
diag(adj) <- FALSE

g      <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
V(g)$name <- loc_ids
comps  <- components(g)

loc_hulls$group <- factor(comps$membership[match(loc_hulls$locationID,
                                                  names(comps$membership))])
loc_centroids$group <- factor(
  comps$membership[match(loc_centroids$locationID, names(comps$membership))]
)

cat(sprintf("=== Connectivity (threshold = %d m) ===\n", THRESHOLD_M))
cat(sprintf("Number of geographic groups (connected components): %d\n", comps$no))
cat(sprintf("Largest group size: %d locations\n\n", max(comps$csize)))

# Summary per group
grp_summary <- data.frame(
  group      = seq_len(comps$no),
  n_locations = comps$csize,
  EOs        = sapply(seq_len(comps$no), function(g_id) {
    locs <- names(comps$membership)[comps$membership == g_id]
    paste(sort(unique(loc_hulls$EOCode[loc_hulls$locationID %in% locs])),
          collapse = ", ")
  }),
  stringsAsFactors = FALSE
)
cat("=== Group summary ===\n")
print(grp_summary)
cat("\n")

# ============================================================
# 4b. HULL AREA & GENETIC DRIFT PROXY PER GROUP
#     Union of member location hulls (UTM) gives total connected
#     habitat area. Smaller area → smaller Ne → stronger drift.
#     Drift index: 0 = weakest drift (largest group), 1 = strongest.
# ============================================================
loc_hulls$hull_area_m2 <- as.numeric(st_area(loc_hulls))

grp_area <- do.call(rbind, lapply(
  sort(unique(as.integer(as.character(loc_hulls$group)))), function(gid) {
    sub_h      <- loc_hulls[as.integer(as.character(loc_hulls$group)) == gid, ]
    union_area <- as.numeric(st_area(st_union(sub_h)))
    data.frame(
      group       = gid,
      n_locations = nrow(sub_h),
      EOs         = paste(sort(unique(sub_h$EOCode)), collapse = ", "),
      area_m2     = round(union_area, 1),
      area_ha     = round(union_area / 10000, 4),
      stringsAsFactors = FALSE
    )
  }
))

a_min <- min(grp_area$area_ha)
a_max <- max(grp_area$area_ha)
grp_area$drift_index <- round(
  if (a_max > a_min) 1 - (grp_area$area_ha - a_min) / (a_max - a_min)
  else rep(0.5, nrow(grp_area)),
  3
)

cat("=== Group areas and drift index (highest drift first) ===\n")
print(grp_area[order(grp_area$drift_index, decreasing = TRUE), ], row.names = FALSE)
cat("\n")

write.csv(grp_area, file.path(out_dir, "EO_group_areas.csv"), row.names = FALSE)

# Attach area and drift index to loc_centroids
loc_centroids <- merge(
  loc_centroids,
  grp_area[, c("group", "area_ha", "drift_index")],
  by.x = "group", by.y = "group", all.x = TRUE
)

# ============================================================
# 5. CONNECTIVITY EDGE LIST (for mapping)
# ============================================================
edge_list <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)
colnames(edge_list) <- c("from", "to")

# Attach coordinates for each endpoint
edge_list <- merge(edge_list,
                   loc_centroids[, c("locationID", "lon", "lat")],
                   by.x = "from", by.y = "locationID")
edge_list <- merge(edge_list,
                   loc_centroids[, c("locationID", "lon", "lat")],
                   by.x = "to", by.y = "locationID",
                   suffixes = c("_from", "_to"))

# ============================================================
# 6. OUTPUT TABLES
# ============================================================
loc_out <- loc_centroids[order(as.integer(as.character(loc_centroids$group)),
                               loc_centroids$EOCode, loc_centroids$locationID), ]
loc_out <- merge(loc_out, grp_summary[, c("group", "EOs")],
                 by = "group", suffixes = c("", "_all_EOs"))
write.csv(loc_out, file.path(out_dir, "EO_location_groups.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 6a. FULL PAIRWISE LOCATION CONNECTIVITY TABLE
#     All n*(n-1)/2 location pairs with distance and connectivity status
# ------------------------------------------------------------------
loc_ids_vec <- loc_hulls$locationID
pair_rows <- do.call(rbind, lapply(seq_len(n_loc - 1), function(i) {
  do.call(rbind, lapply(seq(i + 1, n_loc), function(j) {
    li <- loc_ids_vec[i]
    lj <- loc_ids_vec[j]
    gi <- as.integer(as.character(loc_hulls$group[i]))
    gj <- as.integer(as.character(loc_hulls$group[j]))
    eoi <- loc_hulls$EOCode[i]
    eoj <- loc_hulls$EOCode[j]
    d   <- dist_mat[i, j]
    data.frame(
      locationID_A   = li,
      EOCode_A       = eoi,
      group_A        = gi,
      locationID_B   = lj,
      EOCode_B       = eoj,
      group_B        = gj,
      distance_m     = round(d, 1),
      connected_500m = d <= THRESHOLD_M,
      same_EO        = eoi == eoj,
      same_group     = gi == gj,
      link_type      = ifelse(d > THRESHOLD_M, "not connected",
                       ifelse(eoi == eoj, "within-EO", "between-EO")),
      stringsAsFactors = FALSE
    )
  }))
}))
pair_rows <- pair_rows[order(pair_rows$distance_m), ]

write.csv(pair_rows, file.path(out_dir, "EO_pairwise_connectivity.csv"),
          row.names = FALSE)

cat("=== Pairwise connectivity table (connected pairs only) ===\n")
print(pair_rows[pair_rows$connected_500m, ], row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------
# 6b. GROUP-TO-GROUP MINIMUM DISTANCE TABLE
#     Minimum hull-to-hull distance between every pair of groups
# ------------------------------------------------------------------
grp_ids <- sort(unique(as.integer(as.character(loc_hulls$group))))
grp_pair_rows <- do.call(rbind, lapply(seq_len(length(grp_ids) - 1), function(gi) {
  do.call(rbind, lapply(seq(gi + 1, length(grp_ids)), function(gj) {
    g1 <- grp_ids[gi]; g2 <- grp_ids[gj]
    locs1 <- which(as.integer(as.character(loc_hulls$group)) == g1)
    locs2 <- which(as.integer(as.character(loc_hulls$group)) == g2)
    min_d  <- min(dist_mat[locs1, locs2, drop = FALSE])
    eos1   <- paste(sort(unique(loc_hulls$EOCode[locs1])), collapse = ", ")
    eos2   <- paste(sort(unique(loc_hulls$EOCode[locs2])), collapse = ", ")
    data.frame(
      group_A      = g1,
      EOs_A        = eos1,
      n_locs_A     = length(locs1),
      group_B      = g2,
      EOs_B        = eos2,
      n_locs_B     = length(locs2),
      min_dist_m   = round(min_d, 1),
      connected    = min_d <= THRESHOLD_M,
      stringsAsFactors = FALSE
    )
  }))
}))
grp_pair_rows <- grp_pair_rows[order(grp_pair_rows$min_dist_m), ]

write.csv(grp_pair_rows, file.path(out_dir, "EO_group_distances.csv"),
          row.names = FALSE)

cat("=== Closest group pairs (top 15) ===\n")
print(head(grp_pair_rows, 15), row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------
# 6c. CONNECTIVITY SUMMARY
# ------------------------------------------------------------------
conn_summary <- data.frame(
  metric = c(
    "Total location pairs evaluated",
    "Connected pairs (<=500m)",
    "  Within-EO connections",
    "  Between-EO connections",
    "Isolated locations (no neighbour within 500m)",
    "Number of geographic groups (connected components)",
    "Groups with >1 location",
    "Closest unconnected group pair (m)",
    "Farthest connected pair (m)"
  ),
  value = c(
    nrow(pair_rows),
    sum(pair_rows$connected_500m),
    sum(pair_rows$connected_500m & pair_rows$same_EO),
    sum(pair_rows$connected_500m & !pair_rows$same_EO),
    sum(degree(g) == 0),
    comps$no,
    sum(comps$csize > 1),
    round(min(grp_pair_rows$min_dist_m[!grp_pair_rows$connected]), 1),
    round(max(pair_rows$distance_m[pair_rows$connected_500m]), 1)
  )
)
write.csv(conn_summary, file.path(out_dir, "EO_connectivity_summary.csv"),
          row.names = FALSE)
cat("=== Connectivity summary ===\n")
print(conn_summary, row.names = FALSE)
cat("\n")

# Colours and map extent — used in network, heatmap, and map figures
eos        <- sort(unique(loc_hulls$EOCode))
n_eos      <- length(eos)
eo_colours <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_eos),
  eos
)
lat_range <- range(loc_centroids$lat) + c(-0.3, 0.3)
lon_range <- range(loc_centroids$lon) + c(-0.3, 0.3)

# ============================================================
# 7. CONNECTIVITY NETWORK FIGURE  (graph layout — no basemap)
# ============================================================
NEAR_MISS_M <- 2000   # show near-miss edges up to this distance

near_miss_pairs <- pair_rows[!pair_rows$connected_500m &
                               pair_rows$distance_m <= NEAR_MISS_M, ]

# Build layout graph: use connected + near-miss edges with different weights
# so that nearby locations are pulled together in the layout
g_layout <- g
if (nrow(near_miss_pairs) > 0) {
  g_layout <- add_edges(
    g_layout,
    as.vector(rbind(as.character(near_miss_pairs$locationID_A),
                    as.character(near_miss_pairs$locationID_B)))
  )
}
E(g_layout)$weight <- c(rep(2.0, ecount(g)),
                         rep(0.4, nrow(near_miss_pairs)))

# Fruchterman-Reingold layout
set.seed(42)
lay <- layout_with_fr(g_layout, weights = E(g_layout)$weight)
rownames(lay) <- V(g_layout)$name

# Node data frame
node_df <- data.frame(
  locationID = as.integer(rownames(lay)),
  x          = lay[, 1],
  y          = lay[, 2],
  stringsAsFactors = FALSE
)
node_df <- merge(node_df,
                 loc_centroids[, c("locationID", "EOCode", "group", "n_events",
                                   "area_ha", "drift_index")],
                 by = "locationID")
node_df$label    <- paste0(node_df$EOCode, "\n(", node_df$locationID, ")")
node_df$isolated <- degree(g)[as.character(node_df$locationID)] == 0

# Attach layout positions to edge tables
add_layout <- function(pairs, lay) {
  if (nrow(pairs) == 0) return(pairs)
  pairs$x_A <- lay[as.character(pairs$locationID_A), 1]
  pairs$y_A <- lay[as.character(pairs$locationID_A), 2]
  pairs$x_B <- lay[as.character(pairs$locationID_B), 1]
  pairs$y_B <- lay[as.character(pairs$locationID_B), 2]
  pairs$mid_x <- (pairs$x_A + pairs$x_B) / 2
  pairs$mid_y <- (pairs$y_A + pairs$y_B) / 2
  pairs
}
conn_edges_net <- add_layout(pair_rows[pair_rows$connected_500m, ],    lay)
near_edges_net <- add_layout(near_miss_pairs,                           lay)

# Convex hull enclosures for multi-location groups (with drift_index)
grp_hull_net <- do.call(rbind, Filter(Negate(is.null), lapply(
  sort(unique(as.integer(as.character(node_df$group)))), function(gid) {
    sub <- node_df[as.integer(as.character(node_df$group)) == gid, ]
    if (nrow(sub) < 2) return(NULL)
    pts <- sub[, c("x", "y")]
    if (nrow(pts) == 2) {
      pts <- rbind(pts, data.frame(x = mean(pts$x) + 0.3,
                                   y = mean(pts$y) + 0.3))
    }
    idx <- chull(pts$x, pts$y); idx <- c(idx, idx[1])
    di  <- grp_area$drift_index[grp_area$group == gid]
    data.frame(x = pts$x[idx], y = pts$y[idx],
               group = factor(gid), drift_index = di)
  }
)))

# Group centroid labels for area / drift annotations (all groups)
grp_label_net <- do.call(rbind, lapply(
  sort(unique(as.integer(as.character(node_df$group)))), function(gid) {
    sub <- node_df[as.integer(as.character(node_df$group)) == gid, ]
    ai  <- grp_area$area_ha[grp_area$group == gid]
    di  <- grp_area$drift_index[grp_area$group == gid]
    data.frame(x = mean(sub$x), y = mean(sub$y),
               label       = sprintf("%.2f ha\nDI=%.2f", ai, di),
               drift_index = di,
               stringsAsFactors = FALSE)
  }
))

p_net <- ggplot() +

  # ── FIRST FILL SCALE: drift index ──────────────────────────
  # Group enclosures (multi-location groups) — fill = drift index
  {if (!is.null(grp_hull_net) && nrow(grp_hull_net) > 0)
    geom_polygon(data = grp_hull_net,
                 aes(x = x, y = y, group = group, fill = drift_index),
                 color = "grey40", linewidth = 0.55, alpha = 0.60)
  } +

  # Isolated node halos (single-location groups) — fill = drift index
  geom_point(data = node_df[node_df$isolated, ],
             aes(x = x, y = y, fill = drift_index),
             shape = 21, size = 13, color = "grey40",
             stroke = 0.5, alpha = 0.55) +

  scale_fill_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#D6604D"),
    name    = "Drift index\n(area proxy)",
    limits  = c(0, 1),
    breaks  = c(0, 0.5, 1),
    labels  = c("0\n(large\nhabitat)", "0.5", "1\n(small\nhabitat)")
  ) +
  new_scale_fill() +

  # ── Edges ──────────────────────────────────────────────────
  # Near-miss edges
  {if (nrow(near_edges_net) > 0)
    geom_segment(data = transform(near_edges_net, link_type = "near-miss"),
                 aes(x = x_A, y = y_A, xend = x_B, yend = y_B,
                     color = link_type),
                 linewidth = 0.5, linetype = "dashed", alpha = 0.75)
  } +

  # Connected edges — coloured by link type
  {if (nrow(conn_edges_net) > 0)
    geom_segment(data = conn_edges_net,
                 aes(x = x_A, y = y_A, xend = x_B, yend = y_B,
                     color = link_type),
                 linewidth = 1.6, alpha = 0.92)
  } +

  # Distance labels on connected edges
  {if (nrow(conn_edges_net) > 0)
    geom_label(data = conn_edges_net,
               aes(x = mid_x, y = mid_y,
                   label = paste0(round(distance_m), " m")),
               size = 2.5, fill = "white", color = "grey20",
               linewidth = 0.2, alpha = 0.90,
               label.padding = unit(0.10, "lines"))
  } +

  # Area / drift index labels at group centroids
  geom_text(data = grp_label_net,
            aes(x = x, y = y, label = label),
            size = 2.1, color = "grey20", lineheight = 0.85,
            fontface = "italic", vjust = 2.2) +

  # ── SECOND FILL SCALE: EO identity ─────────────────────────
  # Nodes — shape reflects isolation status
  geom_point(data = node_df,
             aes(x = x, y = y, fill = EOCode, size = n_events,
                 shape = isolated),
             color = "white", stroke = 0.7, alpha = 0.95) +

  # Node labels
  geom_label_repel(data = node_df,
                   aes(x = x, y = y, label = label, fill = EOCode),
                   color = "white", fontface = "bold", size = 2.6,
                   lineheight = 0.80, box.padding = 0.35,
                   point.padding = 0.30, label.size = NA,
                   alpha = 0.88, max.overlaps = 60) +

  scale_shape_manual(values  = c("TRUE" = 21, "FALSE" = 23),
                     name    = "Status",
                     labels  = c("TRUE" = "Isolated", "FALSE" = "Connected")) +
  scale_fill_manual(values   = eo_colours, name = "EO", guide = "none") +
  scale_color_manual(
    values = c("within-EO"  = "#1A5276", "between-EO" = "#B03A2E",
               "near-miss"  = "grey60"),
    limits = c("within-EO", "between-EO", "near-miss"),
    name   = "Connection",
    labels = c("within-EO"  = paste0("Within-EO (<=", THRESHOLD_M, " m)"),
               "between-EO" = paste0("Between-EO (<=", THRESHOLD_M, " m)"),
               "near-miss"  = paste0("Near-miss (", THRESHOLD_M, " m\u20132 km,\nlayout only)")),
    drop   = FALSE
  ) +
  scale_size_continuous(name = "N collection\nevents", range = c(3, 11)) +
  guides(
    shape = guide_legend(override.aes = list(fill = "grey40", size = 5)),
    size  = guide_legend(override.aes = list(fill = "grey40", shape = 21)),
    color = guide_legend(
      override.aes = list(
        linetype  = c("within-EO" = "solid", "between-EO" = "solid",
                      "near-miss" = "dashed"),
        linewidth = 1.2
      )
    )
  ) +

  labs(
    title    = "Habitat fragmentation severs pollinator connectivity and concentrates genetic drift across Lepidium papilliferum populations",
    subtitle = sprintf(
      "C3 Hypothesis \u2014 Stages 1 & 4  |  %d geographic groups  |  %d / %d locations isolated (%.0f%%)  |  %d connected pairs  |  dashed = near-miss (500 m\u20132 km)",
      comps$no, sum(degree(g) == 0), n_loc,
      100 * sum(degree(g) == 0) / n_loc,
      sum(pair_rows$connected_500m)
    ),
    caption  = paste0(
      "Layout: Fruchterman-Reingold algorithm. Near-miss edges used for layout only.\n",
      "Enclosure/halo colour = drift index (DI): blue = large habitat area (weak drift), red = small area (strong drift). ",
      "DI = 1 \u2212 (area \u2212 min) / (max \u2212 min).  Labels show group union area (ha) and DI.\n",
      "Diamonds = connected locations; circles = isolated.  Node size = N collection events.  C3 Hypothesis Stage 1."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5,
                                   margin = margin(b = 4)),
    plot.subtitle   = element_text(size = 9, color = "grey35", hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.caption    = element_text(size = 7.5, color = "grey50", hjust = 0.5,
                                   margin = margin(t = 6)),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(12, 12, 12, 12)
  )

ggsave(file.path(out_dir, "EO_connectivity_network.pdf"),
       p_net, width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "EO_connectivity_network.png"),
       p_net, width = 14, height = 10, dpi = 300)

# ============================================================
# 8. PAIRWISE DISTANCE HEATMAP
# ============================================================

# Order locations: by group, then EO, then locationID
loc_order <- loc_centroids[
  order(as.integer(as.character(loc_centroids$group)),
        loc_centroids$EOCode,
        loc_centroids$locationID),
  "locationID"
]
# Axis label: EO + locationID
loc_labels <- setNames(
  paste0(loc_centroids$EOCode[match(loc_order, loc_centroids$locationID)],
         "\n(", loc_order, ")"),
  as.character(loc_order)
)

# Long-format distance matrix
heat_df <- do.call(rbind, lapply(as.character(loc_order), function(la) {
  data.frame(
    loc_A      = la,
    loc_B      = as.character(loc_order),
    distance_m = dist_mat[la, as.character(loc_order)],
    stringsAsFactors = FALSE
  )
}))
heat_df$loc_A <- factor(heat_df$loc_A, levels = as.character(loc_order))
heat_df$loc_B <- factor(heat_df$loc_B, levels = as.character(loc_order))

# Connectivity category
heat_df$category <- cut(
  heat_df$distance_m,
  breaks = c(-Inf, 0, 500, 2000, 10000, Inf),
  labels = c("Same location", "Connected (<=500 m)",
             "Near-miss (500 m-2 km)", "Close (2-10 km)", "Distant (>10 km)")
)

cat_colours <- c(
  "Same location"          = "#1A1A1A",
  "Connected (<=500 m)"    = "#1E8449",
  "Near-miss (500 m-2 km)" = "#F39C12",
  "Close (2-10 km)"        = "#CB4335",
  "Distant (>10 km)"       = "#7B241C"
)

# EO colour strip (top annotation bar)
strip_df <- data.frame(
  loc   = factor(as.character(loc_order), levels = as.character(loc_order)),
  EOCode = loc_centroids$EOCode[match(loc_order, loc_centroids$locationID)],
  y     = -0.6,
  stringsAsFactors = FALSE
)

# Group strip
strip_grp <- data.frame(
  loc   = factor(as.character(loc_order), levels = as.character(loc_order)),
  group = as.integer(as.character(
    loc_centroids$group[match(loc_order, loc_centroids$locationID)]
  )),
  y = -1.2,
  stringsAsFactors = FALSE
)

p_heat <- ggplot(heat_df, aes(x = loc_B, y = loc_A, fill = category)) +
  geom_tile(color = "white", linewidth = 0.15) +

  # EO colour strip on both axes (shown as a separate tile row/col)
  geom_tile(data = strip_df,
            aes(x = loc, y = factor("EO"), fill = EOCode),
            color = "white", linewidth = 0.1,
            inherit.aes = FALSE) +

  scale_fill_manual(
    values = c(cat_colours, eo_colours),
    breaks = names(cat_colours),   # only show distance categories in legend
    name   = "Pairwise distance",
    guide  = guide_legend(order = 1)
  ) +

  scale_x_discrete(labels = loc_labels) +
  scale_y_discrete(labels = c(loc_labels, "EO" = "EO")) +

  labs(
    title    = "Pairwise distance matrix — Lepidium papilliferum locations",
    subtitle = sprintf(
      "Locations ordered by geographic group  |  %d groups  |  %d connected pairs (<=500 m)  |  %d isolated locations",
      comps$no, sum(pair_rows$connected_500m), sum(degree(g) == 0)
    ),
    caption  = "Colour indicates pairwise hull-to-hull distance. Top strip = EO identity. C3 Hypothesis Stage 1.",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y       = element_text(size = 7),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 9, color = "grey35"),
    plot.caption      = element_text(size = 8, color = "grey50", hjust = 0),
    legend.position   = "right",
    legend.key.size   = unit(0.45, "cm"),
    panel.grid        = element_blank()
  )

ggsave(file.path(out_dir, "EO_distance_heatmap.pdf"),
       p_heat, width = 11, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "EO_distance_heatmap.png"),
       p_heat, width = 11, height = 10, dpi = 300)

# ============================================================
# 9. HIERARCHICAL CLUSTERING OF GROUPS (group centroid distances)
#    Only if ≥ 3 groups exist
# ============================================================
if (comps$no >= 3) {

  # Group centroids (mean of member location centroids)
  grp_coords <- do.call(rbind, lapply(seq_len(comps$no), function(g_id) {
    locs <- loc_centroids[loc_centroids$group == g_id, ]
    data.frame(group = g_id,
               lat   = mean(locs$lat),
               lon   = mean(locs$lon))
  }))

  # Haversine-like distance via st_distance on group centroids
  grp_sf    <- st_as_sf(grp_coords, coords = c("lon", "lat"), crs = 4326)
  grp_sf_utm <- st_transform(grp_sf, crs = 32611)
  grp_dist  <- matrix(as.numeric(st_distance(grp_sf_utm)),
                      nrow = comps$no,
                      dimnames = list(paste0("G", seq_len(comps$no)),
                                      paste0("G", seq_len(comps$no)))) / 1000  # km

  hc <- hclust(as.dist(grp_dist), method = "ward.D2")

  k_max <- min(6, comps$no - 1)
  if (k_max >= 2) {
    sil_scores <- sapply(2:k_max, function(k) {
      cl  <- cutree(hc, k = k)
      sil <- silhouette(cl, as.dist(grp_dist))
      mean(sil[, "sil_width"])
    })
    names(sil_scores) <- paste0("k=", 2:k_max)
    cat("=== Silhouette scores (group-level clustering) ===\n")
    print(round(sil_scores, 3))
    best_k <- which.max(sil_scores) + 1
    cat(sprintf("\nOptimal k (group clusters): %d\n\n", best_k))
    grp_coords$macro_cluster <- factor(cutree(hc, k = best_k))
    loc_centroids$macro_cluster <- grp_coords$macro_cluster[
      match(as.integer(as.character(loc_centroids$group)), grp_coords$group)
    ]
  } else {
    best_k <- 1
    grp_coords$macro_cluster <- factor(1)
    loc_centroids$macro_cluster <- factor(1)
    cat("Too few groups for silhouette analysis; single macro-cluster assigned.\n\n")
  }

  # Dendrogram — relabel leaves: Group | EO(s) | location IDs
  hc_plot <- hc
  hc_plot$labels <- sapply(seq_len(comps$no), function(i) {
    eos  <- grp_summary$EOs[grp_summary$group == i]
    locs <- sort(loc_centroids$locationID[
      as.integer(as.character(loc_centroids$group)) == i
    ])
    paste0("G", i, " | ", eos, " | pop.", paste(locs, collapse = ","))
  })

  # Height at which to draw the cut line for best_k clusters
  cut_h <- if (k_max >= 2 && best_k <= comps$no - 1) {
    mean(hc$height[c(comps$no - best_k, comps$no - best_k + 1)])
  } else NA

  clust_cols <- RColorBrewer::brewer.pal(max(best_k, 3), "Set1")[seq_len(best_k)]

  dend_plot <- function() {
    par(mar = c(10, 4.5, 5, 1))
    plot(hc_plot,
         main  = "",
         xlab  = "",
         ylab  = "Ward's D2 linkage criterion",
         sub   = "",
         cex   = 0.72,
         hang  = -1,
         axes  = TRUE)
    if (k_max >= 2) {
      rect.hclust(hc_plot, k = best_k, border = clust_cols)
    }
    if (!is.na(cut_h)) {
      abline(h = cut_h, lty = 2, col = "grey40", lwd = 1.2)
      text(x = 0.5, y = cut_h * 1.04,
           labels = sprintf("k = %d independent bottleneck lineages", best_k),
           adj = c(0, 0), cex = 0.78, col = "grey20", font = 3)
    }
    # BL labels — colors must match rect.hclust, which colors boxes left-to-right.
    # Map each cutree cluster to its left-to-right rank in the dendrogram.
    cl_assign <- cutree(hc_plot, k = best_k)
    cl_order  <- hc_plot$order
    first_pos <- tapply(match(seq_len(comps$no), cl_order), cl_assign, min)
    lr_rank   <- rank(first_pos, ties.method = "first")  # lr_rank[k] = LR position of cluster k
    for (k_i in seq_len(best_k)) {
      leaves_in_k <- which(cl_assign == k_i)
      positions   <- which(cl_order %in% leaves_in_k)
      mid_pos     <- mean(range(positions))
      col_i       <- lr_rank[k_i]
      mtext(sprintf("BL%d", col_i), side = 1, line = -0.5,
            at = mid_pos, col = clust_cols[col_i], cex = 0.80, font = 2)
    }
    # Title block
    mtext(
      "Independent ancestral bottleneck lineages of Lepidium papilliferum",
      side = 3, line = 3.2, cex = 1.00, font = 2, adj = 0.5
    )
    mtext(
      sprintf(
        "Ward's D2 hierarchical clustering of %d geographic groups (centroid distances)  |  Silhouette-optimal k = %d  |  C3 Hypothesis \u2014 Stage 1",
        comps$no, best_k
      ),
      side = 3, line = 1.9, cex = 0.75, col = "grey30", adj = 0.5
    )
    mtext(
      sprintf(
        "Each bottleneck lineage (BL1\u2013BL%d) defines an independent evolutionary unit for S-allele sampling.\nS-allele surveys must span all %d lineages to capture landscape-wide diversity.",
        best_k, best_k
      ),
      side = 1, line = 5.5, cex = 0.70, col = "grey20", adj = 0.5
    )
  }

  pdf(file.path(out_dir, "EO_clustering_dendrogram.pdf"), width = 14, height = 7)
  dend_plot()
  dev.off()

  png(file.path(out_dir, "EO_clustering_dendrogram.png"),
      width = 14, height = 7, units = "in", res = 300)
  dend_plot()
  dev.off()

} else {
  cat("Fewer than 3 groups — skipping hierarchical clustering of groups.\n\n")
  loc_centroids$macro_cluster <- factor(1)
}

# ============================================================
# 8. MAP
# ============================================================

# Colour by EO
id_map <- map_data("state", region = "idaho")

# Convert hull polygons back to WGS84 for ggplot
loc_hulls_wgs84 <- st_transform(loc_hulls, crs = 4326)
hull_df <- do.call(rbind, lapply(seq_len(nrow(loc_hulls_wgs84)), function(i) {
  geom_type <- st_geometry_type(loc_hulls_wgs84[i, ])
  if (geom_type %in% c("POLYGON")) {
    coords <- st_coordinates(loc_hulls_wgs84[i, ])
    data.frame(
      lon        = coords[, "X"],
      lat        = coords[, "Y"],
      locationID = loc_hulls_wgs84$locationID[i],
      EOCode     = loc_hulls_wgs84$EOCode[i],
      group      = loc_hulls_wgs84$group[i],
      stringsAsFactors = FALSE
    )
  } else {
    # POINT or LINESTRING (1-2 events): use centroid with small pseudo-polygon
    crd <- st_coordinates(loc_hulls_wgs84[i, ])
    cx  <- mean(crd[, "X"])
    cy  <- mean(crd[, "Y"])
    d   <- 0.002   # ~220m visual buffer for single-point locations
    data.frame(
      lon        = cx + d * cos(seq(0, 2*pi, length.out = 12)),
      lat        = cy + d * sin(seq(0, 2*pi, length.out = 12)),
      locationID = loc_hulls_wgs84$locationID[i],
      EOCode     = loc_hulls_wgs84$EOCode[i],
      group      = loc_hulls_wgs84$group[i],
      stringsAsFactors = FALSE
    )
  }
}))

p_map <- ggplot() +
  # Idaho basemap
  geom_polygon(data = id_map,
               aes(x = long, y = lat, group = group),
               fill = "#F5F5F0", color = "grey60", linewidth = 0.4) +
  coord_fixed(ratio = 1.3, xlim = lon_range, ylim = lat_range) +

  # Convex hull polygons per location (filled, semi-transparent)
  geom_polygon(data = hull_df,
               aes(x = lon, y = lat, group = locationID, fill = EOCode),
               alpha = 0.30, color = NA) +
  geom_polygon(data = hull_df,
               aes(x = lon, y = lat, group = locationID, color = EOCode),
               fill = NA, linewidth = 0.55) +

  # Connectivity edges (500m links)
  {if (nrow(edge_list) > 0)
    geom_segment(data = edge_list,
                 aes(x = lon_from, y = lat_from,
                     xend = lon_to,   yend = lat_to),
                 color = "grey30", linewidth = 0.55, linetype = "dashed",
                 alpha = 0.70)
  } +

  # Location centroids sized by n_events
  geom_point(data = loc_centroids,
             aes(x = lon, y = lat, fill = EOCode,
                 size = n_events),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.95) +

  # Location labels
  geom_label_repel(data = loc_centroids,
                   aes(x = lon, y = lat,
                       label = paste0(EOCode, "\n", locationID),
                       fill  = EOCode),
                   color = "white", fontface = "bold", size = 2.5,
                   lineheight = 0.80,
                   box.padding = 0.30, point.padding = 0.25,
                   label.size = NA, alpha = 0.88,
                   max.overlaps = 30) +

  # Group label at group centroid
  {
    gl <- do.call(rbind, lapply(sort(unique(loc_centroids$group)), function(g_id) {
      sub <- loc_centroids[loc_centroids$group == g_id, ]
      data.frame(group = g_id,
                 lon   = mean(sub$lon),
                 lat   = mean(sub$lat),
                 label = paste0("Group ", g_id,
                                "\n(", nrow(sub), " loc.)"),
                 stringsAsFactors = FALSE)
    }))
    geom_label(data = gl,
               aes(x = lon, y = lat, label = label),
               fill = "white", color = "grey20", size = 2.8,
               fontface = "bold", alpha = 0.75, label.size = 0.3)
  } +

  scale_fill_manual(values  = eo_colours, name = "EO") +
  scale_color_manual(values = eo_colours, name = "EO") +
  scale_size_continuous(name  = "N collection events", range = c(2, 9)) +

  labs(
    title    = "Spatial connectivity of Lepidium papilliferum collection locations",
    subtitle = paste0("Convex hull polygons per location  |  500 m connectivity threshold  |  ",
                      comps$no, " geographic group(s) (connected components)"),
    caption  = paste0("Dashed lines indicate location pairs within 500 m (pollinator dispersal limit).\n",
                      "Point size proportional to number of seed collection events. ",
                      "C3 Hypothesis — Stage 1."),
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    plot.caption  = element_text(size = 8,  color = "grey50", hjust = 0),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "EO_connectivity_map.pdf"),
       p_map, width = 11, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, "EO_connectivity_map.png"),
       p_map, width = 11, height = 8, dpi = 300)

# ============================================================
# 9. FINAL CONSOLE SUMMARY
# ============================================================
cat("=== Location → Group assignments ===\n")
print(loc_out[, c("EOCode", "locationID", "lat", "lon",
                  "n_events", "group")])

message("\nDone. Outputs written to: ", out_dir)
message("  EO_location_groups.csv")
message("  EO_group_areas.csv")
message("  EO_connectivity_summary.csv")
message("  EO_connectivity_map.pdf / .png")
if (comps$no >= 3) message("  EO_clustering_dendrogram.pdf / .png")
