# =============================================================================
# srk_bl_constants.R — single source of truth for BL ordering and colours
# -----------------------------------------------------------------------------
# Loaded by every SRK pipeline script that orders or colours bottleneck
# lineages so that ordering and palette stay consistent across the project.
#
# Both orderings are DERIVED FROM CSVs mirrored from the
# LEPA_EO_spatial_clustering project, so the entire pipeline reorders
# automatically when sampling grows — only the CSVs in Tables/ need updating.
#
#   BL_ORDER (built from Tables/EO_BL_summary.csv):
#     BLs sorted high-to-low by total habitat area, with within-BL
#     connectivity as tie-breaker:
#       total_area_ha                              (primary, desc)
#       connectivity = n_locations - n_groups      (secondary, desc)
#       BL name                                    (deterministic tie-break)
#     Area is the primary Ne proxy (carrying capacity -> drift floor);
#     connectivity is the secondary stratification for BLs of similar size.
#     Same sort key drives the panel order in the BL drift figure produced
#     by the spatial-clustering project.
#
#   BL_COLORS:
#     Set1 palette mapped to BL by cluster-index, matching the rendered
#     figures in LEPA_EO_spatial_clustering (dendrogram, BL geographic
#     context map, drift panel, fragmentation network). Verified against
#     live bl_strip_cols on 2026-05-18.
#
#   get_eo_order_within_bl(eo_codes = NULL):
#     EO codes ordered by (BL_ORDER, ascending mean Drift_index within BL).
#     Reads Tables/EO_group_BL_summary.csv. EO name is the tie-break so the
#     order is deterministic and matches the Python helper.
# =============================================================================

.BL_SUMMARY_FILE       <- "Tables/EO_BL_summary.csv"
.EO_GROUP_BL_SUMMARY   <- "Tables/EO_group_BL_summary.csv"

# BL -> hex colour (cluster-index Set1 mapping)
BL_COLORS <- c(
  BL1 = "#984EA3",  # purple
  BL2 = "#377EB8",  # blue
  BL3 = "#E41A1C",  # red
  BL4 = "#FF7F00",  # orange
  BL5 = "#4DAF4A"   # green
)

# Build BL_ORDER from the per-BL summary CSV at load time.
# Sort key: total habitat area (ha) DESC, then within-BL connectivity DESC,
# then BL name ASC. Area is the primary Ne proxy; connectivity is the
# secondary stratification for BLs of similar size.
.derive_bl_order <- function(path = .BL_SUMMARY_FILE) {
  if (!file.exists(path))
    stop("BL summary file not found: ", path,
         " — mirror it from LEPA_EO_spatial_clustering/data/EO_BL_summary.csv")
  s <- read.csv(path, stringsAsFactors = FALSE)
  s$connectivity <- s$n_locations - s$n_groups
  s <- s[order(-s$total_area_ha, -s$connectivity, s$BL), ]
  s$BL
}
BL_ORDER <- .derive_bl_order()

# Numerical (alphanumeric) BL order — used for SCATTER-PLOT LEGENDS only
# (TP1, TP2, GFS scatters, allele-accumulation curve legend, etc.) where the
# x/y axes are not BL and the legend just needs an easy-to-read order.
# Pass to ggplot scales as `breaks = BL_ORDER_NUMERIC`. Factor levels on the
# data should still be BL_ORDER (area-then-connectivity) so axis / facet /
# table order remains driven by the Ne-proxy ranking.
BL_ORDER_NUMERIC <- sort(BL_ORDER)

# Order EOs by (BL_ORDER, ascending mean Drift_index, EO name).
get_eo_order_within_bl <- function(eo_codes = NULL,
                                   bl_summary_path = .EO_GROUP_BL_SUMMARY) {
  if (!file.exists(bl_summary_path))
    stop("EO group BL summary file not found: ", bl_summary_path)
  s <- read.csv(bl_summary_path, stringsAsFactors = FALSE)
  # Composite EO entries like "EO118; EO76" represent one geographic group
  # spanning two EOs — split them so downstream code sees both codes (each
  # inherits the group's mean Drift_index).
  parts <- strsplit(s$EO, "[;,]\\s*")
  s <- s[rep(seq_len(nrow(s)), lengths(parts)), ]
  s$EO <- trimws(unlist(parts))
  eo_di <- aggregate(Drift_index ~ EO + BL, data = s, FUN = mean)
  eo_di$BL <- factor(eo_di$BL, levels = BL_ORDER)
  # Secondary sort on EO name keeps ties deterministic and identical to Python.
  eo_di <- eo_di[order(eo_di$BL, eo_di$Drift_index, eo_di$EO), ]
  ordered_eos <- eo_di$EO

  if (!is.null(eo_codes)) {
    # Preserve order for EOs present in the summary; append any extras lexically.
    extras <- setdiff(eo_codes, ordered_eos)
    ordered_eos <- c(intersect(ordered_eos, eo_codes), sort(extras))
  }
  ordered_eos
}
