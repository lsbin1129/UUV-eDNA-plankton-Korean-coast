# =============================================================================
# Figure 9. Geospatial distribution of the top 20 dominant plankton genera
#           across five Korean coastal sea regions
#
# Visualization approach: glow-style bubble maps
# =============================================================================
# Required packages:
#   install.packages(c("readxl", "dplyr", "ggplot2",
#                      "rnaturalearth", "rnaturalearthdata",
#                      "sf", "patchwork"))

library(readxl)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# =============================================================================
# 1. Load and prepare data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx",
                  col_names = FALSE, .name_repair = "minimal")

locations  <- as.character(unlist(raw[2, 3:ncol(raw)]))
lats       <- as.numeric(unlist(raw[3, 3:ncol(raw)]))
lons       <- as.numeric(unlist(raw[4, 3:ncol(raw)]))
genus_vec  <- as.character(unlist(raw[6:nrow(raw), 1]))
tax_vec    <- as.character(unlist(raw[6:nrow(raw), 2]))

counts_raw <- raw[6:nrow(raw), 3:ncol(raw)]
counts_mat <- matrix(
  as.numeric(unlist(counts_raw)),
  nrow = nrow(counts_raw),
  ncol = ncol(counts_raw)
)
rownames(counts_mat) <- genus_vec

# ── Reclassify Noctiluca as Zooplankton ──────────────────────────────────────
tax_vec[genus_vec == "Noctiluca"] <- "동물플랑크톤"

# =============================================================================
# 2. Select top 20 genera by total read abundance
# =============================================================================
total_reads  <- rowSums(counts_mat)
top20_idx    <- order(total_reads, decreasing = TRUE)[1:20]
top20_genera <- genus_vec[top20_idx]
global_max   <- max(counts_mat[top20_idx, ])   # used for log-normalisation

# =============================================================================
# 3. Plot settings
# =============================================================================
# Sea regions
sea_regions <- c("Yellow Sea", "South Sea", "East Sea",
                 "Jeju Sea", "Seogwipo Sea")

# Colour palettes per region (light outer ring → saturated core)
sea_palettes <- list(
  "Yellow Sea"   = c("#fef9ec", "#fde3a7", "#f3a01b"),
  "South Sea"    = c("#fef0ef", "#f4a49e", "#e85344"),
  "East Sea"     = c("#eef6fd", "#95c8ee", "#389adb"),
  "Jeju Sea"     = c("#eafaf2", "#8ce8b8", "#33cd74"),
  "Seogwipo Sea" = c("#f5eefa", "#c99fdc", "#9b59b6")
)

N_RINGS   <- 10     # number of concentric glow rings per bubble
MAX_SIZE  <- 18     # maximum bubble size (ggplot2 size units)
MIN_SIZE  <- 1.5    # minimum bubble size

# Map coordinate extent
lon_range <- c(124.0, 132.0)
lat_range <- c(32.2,  39.6)

# =============================================================================
# 4. Load base map
# =============================================================================
world <- ne_countries(scale = "medium", returnclass = "sf")
korea_region <- world[world$name %in%
                        c("South Korea", "North Korea", "China", "Japan"), ]

# =============================================================================
# 5. Function: generate glow bubble map for one genus
# =============================================================================
make_glow_map <- function(genus_name, genus_index) {

  reads_vec <- as.numeric(counts_mat[genus_index, ])

  # Log-normalised abundance ratio (0–1 scale)
  ratio_vec <- log1p(pmax(reads_vec, 0)) / log1p(global_max)

  # Initialise base map
  p <- ggplot() +
    geom_sf(data  = korea_region,
            fill  = "grey93",
            colour = "grey55",
            linewidth = 0.3) +
    coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE)

  # Add glow rings — outermost ring drawn first (k = N_RINGS),
  # innermost last so it renders on top (k = 1)
  for (k in N_RINGS:1) {
    frac    <- k / N_RINGS
    pt_size <- MIN_SIZE + (1 - frac) * (MAX_SIZE - MIN_SIZE)
    alpha_k <- (1 - frac)^1.3 * 0.82

    for (sea in sea_regions) {
      idx      <- which(locations == sea)
      pos_mask <- reads_vec[idx] > 0
      if (sum(pos_mask) == 0) next

      idx_pos   <- idx[pos_mask]
      ratio_pos <- ratio_vec[idx_pos]

      # Interpolated ring colour for this sea and ring index
      pal_func <- colorRampPalette(sea_palettes[[sea]])
      ring_cols <- pal_func(N_RINGS)
      col_k     <- ring_cols[N_RINGS - k + 1]

      df_ring <- data.frame(
        lon      = lons[idx_pos],
        lat      = lats[idx_pos],
        pt_size  = pt_size * ratio_pos,
        alpha_pt = pmin(alpha_k * ratio_pos + 0.01, 1)
      )

      p <- p +
        geom_point(
          data    = df_ring,
          mapping = aes(x = lon, y = lat,
                        size = pt_size, alpha = alpha_pt),
          shape   = 16,
          colour  = col_k
        ) +
        scale_size_identity() +
        scale_alpha_identity()
    }
  }

  # Detected stations: open white circle with grey border
  detected_df <- data.frame(
    lon = lons[reads_vec > 0],
    lat = lats[reads_vec > 0]
  )
  if (nrow(detected_df) > 0) {
    p <- p +
      geom_point(data   = detected_df,
                 aes(x = lon, y = lat),
                 shape  = 21, size = 1.4, stroke = 0.4,
                 fill   = "white", colour = "grey30", alpha = 0.9)
  }

  # Absent stations: grey cross symbol
  absent_df <- data.frame(
    lon = lons[reads_vec == 0],
    lat = lats[reads_vec == 0]
  )
  if (nrow(absent_df) > 0) {
    p <- p +
      geom_point(data   = absent_df,
                 aes(x = lon, y = lat),
                 shape  = 4, size = 1.5, stroke = 0.55,
                 colour = "grey65", alpha = 0.7)
  }

  # Panel labels and theme
  p <- p +
    labs(
      title    = bquote(italic(.(genus_name))),
      subtitle = paste0("Total: ",
                        format(as.integer(sum(reads_vec)), big.mark = ",")),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title       = element_text(size = 10, face = "bold.italic",
                                      hjust = 0.5, margin = margin(b = 1)),
      plot.subtitle    = element_text(size = 7.5, colour = "grey40",
                                      hjust = 0.5, margin = margin(b = 2)),
      axis.text        = element_text(size = 6),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
      panel.border     = element_rect(colour = "grey50", linewidth = 0.5),
      legend.position  = "none",
      plot.margin      = margin(4, 4, 4, 4)
    )

  return(p)
}

# =============================================================================
# 6. Generate individual panels for all top 20 genera
# =============================================================================
message("Generating bubble maps for top 20 genera...")

map_list <- vector("list", length(top20_genera))
for (i in seq_along(top20_genera)) {
  message(sprintf("  [%2d / 20]  %s", i, top20_genera[i]))
  map_list[[i]] <- make_glow_map(top20_genera[i], top20_idx[i])
}

# =============================================================================
# 7. Assemble composite figure (4 columns × 5 rows)
# =============================================================================
message("\nAssembling composite figure...")

figure9 <- wrap_plots(map_list, ncol = 4, nrow = 5)

# =============================================================================
# 8. Save output
# =============================================================================
ggsave("Figure_09.tiff", plot = figure9,
       width = 22, height = 28, units = "cm",
       dpi = 300, bg = "white", compression = "lzw")

ggsave("Figure_09.pdf", plot = figure9,
       width = 22, height = 28, units = "cm",
       bg = "white")

message("Figure 9 saved: Figure_09.tiff / Figure_09.pdf")
