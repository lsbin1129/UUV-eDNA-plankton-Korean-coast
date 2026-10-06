# =============================================================================
# Figure 9. Geospatial distribution of the top 20 dominant plankton genera
#           across five Korean coastal sea regions
#
#   Layout reproduces the original Figure 9:
#     - 4 x 5 grid of portrait map panels, tagged (a)-(t) at the top-left
#     - genus name (italic) and read summary inside each panel, top-right
#     - grey land, blue rivers, light graticule
#     - glow bubbles per station (soft light halo -> saturated core),
#       coloured by sea region and scaled to log reads
#     - white dot = genus detected at the station; grey x = not detected
#   Organism drawings are added manually afterwards (space is left at the
#   top-right of each panel, under the genus label).
# =============================================================================
# install.packages(c("readxl", "dplyr", "ggplot2", "sf", "patchwork",
#                    "rnaturalearth", "rnaturalearthdata"))

library(readxl)
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

sf_use_s2(FALSE)

# =============================================================================
# 0. User settings
# =============================================================================
# River layer. If you have the river shapefile used for the original
# Figures 1 and 9, give its path here (any line shapefile / GeoPackage in
# WGS84 or any CRS). If NULL, Natural Earth 10m rivers are downloaded
# (only the Han, Namhan and Nakdong rivers are included at that scale).
RIVER_FILE <- NULL            # e.g. "rivers/korea_rivers.shp"

# Coastline detail: "medium" (Natural Earth 50m) or "large" (10m; needs the
# rnaturalearthhires package). LAND_FILE overrides both with a local file.
MAP_SCALE <- "medium"
LAND_FILE <- NULL

# Read summary shown under each genus name: "mean" (mean reads per profile,
# matching Results 3.2) or "total" (sum over all profiles, as in the
# original figure).
SUBTITLE_MODE <- "mean"

# Region colours. "unified" = palette shared with Figures 1 and 8
# (response to Reviewer 1, comment 13); "original" = colours of the
# original Figure 9.
PALETTE <- "unified"

# Output size (cm); the original figure is ~ 0.69 width:height
FIG_WIDTH  <- 22
FIG_HEIGHT <- 32

# =============================================================================
# 1. Load and prepare data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx",
                  col_names = FALSE, .name_repair = "minimal")

locations  <- as.character(unlist(raw[2, 3:ncol(raw)]))
lats       <- as.numeric(unlist(raw[3, 3:ncol(raw)]))
lons       <- as.numeric(unlist(raw[4, 3:ncol(raw)]))
genus_vec  <- as.character(unlist(raw[6:nrow(raw), 1]))

counts_mat <- apply(raw[6:nrow(raw), 3:ncol(raw)], 2,
                    function(x) as.numeric(as.character(x)))
counts_mat[is.na(counts_mat)] <- 0   # empty cells in the spreadsheet = 0 reads
rownames(counts_mat) <- genus_vec

# =============================================================================
# 2. Top 20 genera by mean reads per profile (all 125 profiles,
#    UUV and Net pooled; identical ranking to the total)
# =============================================================================
mean_reads   <- rowMeans(counts_mat)
top20_idx    <- order(mean_reads, decreasing = TRUE)[1:20]
top20_genera <- genus_vec[top20_idx]
global_max   <- max(counts_mat[top20_idx, ])   # common bubble scale

print(data.frame(Panel = paste0("(", letters[1:20], ")"),
                 Genus = top20_genera,
                 MeanReadsPerProfile = round(unname(mean_reads[top20_idx])),
                 TotalReads = unname(rowSums(counts_mat)[top20_idx])))

# =============================================================================
# 3. Plot settings
# =============================================================================
sea_regions <- c("Yellow Sea", "South Sea", "East Sea",
                 "Jeju Sea", "Seogwipo Sea")

region_core <- switch(PALETTE,
  unified  = c("Yellow Sea" = "#D4AC0D", "South Sea" = "#8E44AD",
               "East Sea"   = "#16A085", "Jeju Sea"  = "#A0522D",
               "Seogwipo Sea" = "#C2185B"),
  original = c("Yellow Sea" = "#f3a01b", "South Sea" = "#e85344",
               "East Sea"   = "#389adb", "Jeju Sea"  = "#33cd74",
               "Seogwipo Sea" = "#9b59b6"))

# Glow: N_RINGS concentric discs, drawn largest (pale, transparent) first
# and smallest (saturated, opaque) last, so each bubble has a soft halo
# around a solid core.
N_RINGS    <- 12
MAX_SIZE   <- 15     # halo diameter for the maximum read count (pt-size units)
MIN_SIZE   <- 2.5    # halo diameter for the smallest non-zero count
CORE_FRAC  <- 0.45   # core diameter as a fraction of the halo diameter
ALPHA_OUT  <- 0.10   # opacity of the outermost ring
ALPHA_IN   <- 0.85   # opacity of the innermost ring

RIVER_COL  <- "#5B8FD6"
LAND_FILL  <- "grey92"
LAND_LINE  <- "grey50"

lon_range <- c(124.0, 132.0)
lat_range <- c(32.2,  39.6)

# =============================================================================
# 4. Base map layers
# =============================================================================
bbox_sf <- st_as_sfc(st_bbox(c(xmin = lon_range[1] - 1, xmax = lon_range[2] + 1,
                               ymin = lat_range[1] - 1, ymax = lat_range[2] + 1),
                             crs = 4326))

land <- if (!is.null(LAND_FILE)) {
  st_read(LAND_FILE, quiet = TRUE)
} else {
  rnaturalearth::ne_countries(scale = MAP_SCALE, returnclass = "sf")
}
land <- suppressWarnings(st_crop(st_make_valid(st_transform(land, 4326)), bbox_sf))

rivers <- tryCatch({
  if (!is.null(RIVER_FILE)) {
    st_read(RIVER_FILE, quiet = TRUE)
  } else {
    rnaturalearth::ne_download(scale = 10, type = "rivers_lake_centerlines",
                               category = "physical", returnclass = "sf")
  }
}, error = function(e) {
  message("River layer not available (", conditionMessage(e), ") - drawn without rivers.")
  NULL
})
if (!is.null(rivers)) {
  rivers <- suppressWarnings(st_crop(st_transform(rivers, 4326), bbox_sf))
  # keep rivers on land only (avoids lines running into the sea)
  rivers <- suppressWarnings(st_intersection(rivers, st_union(land)))
}

# =============================================================================
# 5. Function: glow bubble map for one genus
# =============================================================================
make_glow_map <- function(genus_name, genus_index) {

  reads_vec <- as.numeric(counts_mat[genus_index, ])
  ratio_vec <- log1p(pmax(reads_vec, 0)) / log1p(global_max)   # 0-1

  summary_txt <- if (SUBTITLE_MODE == "mean") {
    paste0("Mean: ", format(round(mean(reads_vec)), big.mark = ","), " reads/profile")
  } else {
    paste0("Total: ", format(sum(reads_vec), big.mark = ","))
  }

  p <- ggplot() +
    geom_sf(data = land, fill = LAND_FILL, colour = LAND_LINE, linewidth = 0.25)

  if (!is.null(rivers) && nrow(rivers) > 0) {
    p <- p + geom_sf(data = rivers, colour = RIVER_COL, linewidth = 0.35,
                     lineend = "round")
  }

  # ---- glow rings (largest/palest first, smallest/most saturated last) ----
  ring_df <- list()
  for (sea in sea_regions) {
    idx <- which(locations == sea & reads_vec > 0)
    if (length(idx) == 0) next
    ramp <- colorRampPalette(c("white", region_core[[sea]]))(N_RINGS + 2)[-(1:2)]
    halo <- MIN_SIZE + ratio_vec[idx] * (MAX_SIZE - MIN_SIZE)
    for (k in seq_len(N_RINGS)) {           # k = 1 outermost ... N_RINGS core
      f <- (k - 1) / (N_RINGS - 1)          # 0 -> 1
      ring_df[[length(ring_df) + 1]] <- data.frame(
        lon   = lons[idx],
        lat   = lats[idx],
        size  = halo * (1 - f * (1 - CORE_FRAC)),
        alpha = ALPHA_OUT + f^1.5 * (ALPHA_IN - ALPHA_OUT),
        col   = ramp[k],
        k     = k
      )
    }
  }
  if (length(ring_df) > 0) {
    ring_df <- do.call(rbind, ring_df)
    ring_df <- ring_df[order(ring_df$k, -ring_df$size), ]   # draw order
    p <- p +
      geom_point(data = ring_df,
                 aes(x = lon, y = lat, size = size, alpha = alpha, colour = col),
                 shape = 16, stroke = 0) +
      scale_size_identity() +
      scale_alpha_identity() +
      scale_colour_identity()
  }

  # ---- station markers ----
  detected_df <- data.frame(lon = lons[reads_vec > 0],  lat = lats[reads_vec > 0])
  absent_df   <- data.frame(lon = lons[reads_vec == 0], lat = lats[reads_vec == 0])
  if (nrow(absent_df) > 0) {
    p <- p + geom_point(data = absent_df, aes(x = lon, y = lat),
                        shape = 4, size = 0.9, stroke = 0.35, colour = "grey55")
  }
  if (nrow(detected_df) > 0) {
    p <- p + geom_point(data = detected_df, aes(x = lon, y = lat),
                        shape = 21, size = 0.75, stroke = 0.25,
                        fill = "white", colour = "grey20")
  }

  # ---- labels inside the panel (top-right) ----
  p +
    annotate("text", x = lon_range[2] - 0.18, y = lat_range[2] - 0.22,
             label = genus_name, hjust = 1, vjust = 1,
             fontface = "italic", size = 2.9, colour = "black") +
    annotate("text", x = lon_range[2] - 0.18, y = lat_range[2] - 0.62,
             label = summary_txt, hjust = 1, vjust = 1,
             size = 2.0, colour = "grey35") +
    coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE,
             label_axes = list(bottom = "E", left = "N")) +
    scale_x_continuous(breaks = seq(124, 132, 1)) +
    scale_y_continuous(breaks = seq(33, 39, 1)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.2),
      panel.border     = element_rect(colour = "grey45", linewidth = 0.4),
      axis.text        = element_text(size = 3.6, colour = "grey30"),
      axis.ticks       = element_line(linewidth = 0.2, colour = "grey50"),
      axis.ticks.length = unit(0.6, "mm"),
      legend.position  = "none",
      plot.margin      = margin(12, 3, 2, 3)   # top space for (a)-(t) tags
    )
}

# =============================================================================
# 6. Panels and composite (4 columns x 5 rows, tagged (a)-(t))
# =============================================================================
message("Generating bubble maps for top 20 genera...")
map_list <- lapply(seq_along(top20_genera), function(i) {
  message(sprintf("  [%2d / 20]  %s", i, top20_genera[i]))
  make_glow_map(top20_genera[i], top20_idx[i])
})

figure9 <- wrap_plots(map_list, ncol = 4, nrow = 5) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag          = element_text(family = "serif", size = 13,
                                         hjust = 0, vjust = 1),
        plot.tag.position = "topleft",
        plot.background   = element_rect(fill = "white", colour = NA))

# =============================================================================
# 7. Save
# =============================================================================
ggsave("Figure_09.tiff", plot = figure9,
       width = FIG_WIDTH, height = FIG_HEIGHT, units = "cm",
       dpi = 300, bg = "white", compression = "lzw")
ggsave("Figure_09.pdf", plot = figure9,
       width = FIG_WIDTH, height = FIG_HEIGHT, units = "cm", bg = "white")

message("Figure 9 saved: Figure_09.tiff / Figure_09.pdf")
