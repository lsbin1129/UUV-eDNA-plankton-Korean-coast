# =============================================================================
# Figure 1. Spatial distribution of sampling stations across five Korean
#           coastal sea regions — 6 separate output files
# =============================================================================
# Required packages:
#   install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata",
#                      "ggspatial", "readxl", "dplyr"))

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(readxl)
library(dplyr)

# -----------------------------------------------------------------------------
# 0. Load and prepare sampling data
# -----------------------------------------------------------------------------
raw <- read_excel("MetaK90_location.xlsx", col_names = FALSE)

methods    <- as.character(raw[1, 3:ncol(raw)])
locations  <- as.character(raw[2, 3:ncol(raw)])
lats       <- as.numeric(raw[3, 3:ncol(raw)])
lons       <- as.numeric(raw[4, 3:ncol(raw)])

stations_raw <- data.frame(
  Method    = methods,
  Location  = locations,
  Latitude  = lats,
  Longitude = lons,
  stringsAsFactors = FALSE
) %>% filter(!is.na(Latitude), !is.na(Longitude))

# ── Collapse to one marker per station, classified by sampling pairing ──────
# UUV and net hauls were deployed within a 100-m radius of the same station,
# so the two methods are visually indistinguishable at map scale. Each UUV
# profile is paired with the nearest Net profile in the same region within
# PAIR_MAX_KM (greedy matching in order of increasing distance, each profile
# used at most once). Paired stations are drawn once (at the UUV position);
# unpaired profiles are drawn as "UUV only" / "Net only". The same pairing
# rule is used in Statistics_paired.R.
PAIR_MAX_KM <- 0.2

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad  <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

uuv_i <- which(stations_raw$Method == "Drone")
net_i <- which(stations_raw$Method == "Net")
cand  <- expand.grid(u = uuv_i, n = net_i)
cand  <- cand[stations_raw$Location[cand$u] == stations_raw$Location[cand$n], ]
cand$d <- haversine_km(stations_raw$Latitude[cand$u], stations_raw$Longitude[cand$u],
                       stations_raw$Latitude[cand$n], stations_raw$Longitude[cand$n])
cand  <- cand[cand$d <= PAIR_MAX_KM, ]
cand  <- cand[order(cand$d), ]
pairs <- cand[0, ]
for (k in seq_len(nrow(cand))) {
  if (!(cand$u[k] %in% pairs$u) && !(cand$n[k] %in% pairs$n))
    pairs <- rbind(pairs, cand[k, ])
}
message(sprintf("Paired stations: %d (max UUV-Net distance %.0f m)",
                nrow(pairs), 1000 * max(c(0, pairs$d))))

uuv_unpaired <- setdiff(uuv_i, pairs$u)
net_unpaired <- setdiff(net_i, pairs$n)

stations <- rbind(
  data.frame(stations_raw[pairs$u, c("Location", "Latitude", "Longitude")],
             Pairing = "Paired (UUV + Net)"),
  data.frame(stations_raw[uuv_unpaired, c("Location", "Latitude", "Longitude")],
             Pairing = rep("UUV only", length(uuv_unpaired))),
  data.frame(stations_raw[net_unpaired, c("Location", "Latitude", "Longitude")],
             Pairing = rep("Net only", length(net_unpaired)))
)
print(table(stations$Location, stations$Pairing))

# Consistent region colours across all panels
region_colours <- c(
  "East Sea"     = "#16A085",  # teal
  "Yellow Sea"   = "#D4AC0D",  # mustard/gold
  "South Sea"    = "#8E44AD",  # purple
  "Jeju Sea"     = "#A0522D",  # sienna/brown
  "Seogwipo Sea" = "#C2185B"   # magenta
)

# Factor levels for legend order
stations$Location <- factor(stations$Location,
                            levels = c("East Sea", "Yellow Sea",
                                       "South Sea", "Jeju Sea",
                                       "Seogwipo Sea"))
stations$Pairing <- factor(stations$Pairing,
                          levels = c("Paired (UUV + Net)", "UUV only", "Net only"))

# Convert to sf
stations_sf <- st_as_sf(stations, coords = c("Longitude", "Latitude"),
                        crs = 4326)

# Base map layers
world <- ne_countries(scale = "medium", returnclass = "sf")

# Common theme
map_theme <- theme_bw() +
  theme(
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    panel.border     = element_rect(colour = "black", linewidth = 0.8),
    legend.position  = "none"
  )

# Shared scale/shape settings for all panels
shared_colour <- scale_colour_manual(values = region_colours)
shared_shape  <- scale_shape_manual(
  values = c("Paired (UUV + Net)" = 16, "UUV only" = 24, "Net only" = 17),
  drop = FALSE
)

# -----------------------------------------------------------------------------
# Helper: save one panel
# -----------------------------------------------------------------------------
save_map <- function(plot_obj, filename, width = 12, height = 14) {
  ggsave(paste0(filename, ".tiff"), plot = plot_obj,
         width = width, height = height, units = "cm",
         dpi = 300, compression = "lzw")
  ggsave(paste0(filename, ".pdf"), plot = plot_obj,
         width = width, height = height, units = "cm")
  message("Saved: ", filename)
}

# =============================================================================
# 1. Korea Overview Map
# =============================================================================
fig_overview <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = stations_sf,
          aes(colour = Location, shape = Pairing),
          size = 1.8, stroke = 0.4, alpha = 0.9,
          key_glyph = "point") +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(124.5, 130.5), ylim = c(32.8, 39.0), expand = FALSE, clip = "off") +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tl",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.1, "cm"), width = unit(1.1, "cm")
  ) +
  annotate("text", x = 127.6, y = 36.5, label = "Republic\nof Korea",
           size = 4, fontface = "bold", colour = "grey30") +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  # Legend placed OUTSIDE the map panel (to the right) so it never
  # overlaps sampling points, regardless of point density in any region
  theme_bw() +
  theme(
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 8),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    panel.border      = element_rect(colour = "black", linewidth = 0.8),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.background = element_rect(fill = "white", colour = "grey60",
                                     linewidth = 0.4),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.1, "cm")
  ) +
  guides(
    shape  = guide_legend(title = "Sampling Coverage", order = 1,
                          override.aes = list(size = 2.5, colour = "black",
                                              fill = NA)),
    colour = guide_legend(title = "Location",        order = 2,
                          override.aes = list(shape = 16, size = 2.5))
  )

save_map(fig_overview, "Figure_01_overview", width = 17, height = 14)

# =============================================================================
# 2. South Sea
# =============================================================================
ss <- filter(stations_sf, Location == "South Sea")

fig_south <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = ss, aes(colour = Location, shape = Pairing),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.8, 129.5), ylim = c(34.0, 35.5), expand = FALSE, clip = "off") +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.7) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)",
       title = "South Sea") +
  map_theme +
  theme(
    plot.title = element_text(size = 13, face = "bold",
                              colour = region_colours["South Sea"])
  )

save_map(fig_south, "Figure_01_SouthSea", width = 13, height = 10)

# =============================================================================
# 3. East Sea
# =============================================================================
es <- filter(stations_sf, Location == "East Sea")

fig_east <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = es, aes(colour = Location, shape = Pairing),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(128.2, 130.2), ylim = c(34.9, 39.0), expand = FALSE, clip = "off") +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tl",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.0, "cm"), width = unit(1.0, "cm")
  ) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)",
       title = "East Sea") +
  map_theme +
  theme(
    plot.title = element_text(size = 13, face = "bold",
                              colour = region_colours["East Sea"])
  )

save_map(fig_east, "Figure_01_EastSea", width = 9, height = 14)

# =============================================================================
# 4. Yellow Sea
# =============================================================================
ys <- filter(stations_sf, Location == "Yellow Sea")

fig_yellow <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = ys, aes(colour = Location, shape = Pairing),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.6, 127.8), ylim = c(34.0, 38.2), expand = FALSE, clip = "off") +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tl",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.0, "cm"), width = unit(1.0, "cm")
  ) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)",
       title = "Yellow Sea") +
  map_theme +
  theme(
    plot.title = element_text(size = 13, face = "bold",
                              colour = region_colours["Yellow Sea"])
  )

save_map(fig_yellow, "Figure_01_YellowSea", width = 10, height = 14)

# =============================================================================
# 5. Jeju Sea & Seogwipo Sea
# =============================================================================
jeju <- filter(stations_sf, Location %in% c("Jeju Sea", "Seogwipo Sea"))

fig_jeju <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = jeju, aes(colour = Location, shape = Pairing),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.9, 127.3), ylim = c(32.95, 34.0), expand = FALSE, clip = "off") +
  annotation_scale(location = "bl", width_hint = 0.35, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tl",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.0, "cm"), width = unit(1.0, "cm")
  ) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)",
       title = "Jeju Sea & Seogwipo Sea") +
  map_theme +
  theme(
    plot.title = element_text(size = 13, face = "bold",
                              colour = region_colours["Jeju Sea"])
  )

save_map(fig_jeju, "Figure_01_JejuSeogwipo", width = 13, height = 9)

message("\n All Figure 1 panels saved successfully.")
