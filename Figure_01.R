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

stations <- data.frame(
  Method    = methods,
  Location  = locations,
  Latitude  = lats,
  Longitude = lons,
  stringsAsFactors = FALSE
) %>% filter(!is.na(Latitude), !is.na(Longitude))

# Consistent region colours across all panels
region_colours <- c(
  "East Sea"     = "#4472C4",
  "Yellow Sea"   = "#ED7D31",
  "South Sea"    = "#E00000",
  "Jeju Sea"     = "#70AD47",
  "Seogwipo Sea" = "#7030A0"
)

# Factor levels for legend order
stations$Location <- factor(stations$Location,
                            levels = c("East Sea", "Yellow Sea",
                                       "South Sea", "Jeju Sea",
                                       "Seogwipo Sea"))
stations$Method <- factor(stations$Method,
                          levels = c("Drone", "Net"),
                          labels = c("Uncrewed Underwater Vehicle (UUV)", "Net"))

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
  values = c("Uncrewed Underwater Vehicle (UUV)" = 16, "Net" = 17)
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
          aes(colour = Location, shape = Method),
          size = 1.8, stroke = 0.4, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(124.5, 130.5), ylim = c(32.8, 39.0), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
  annotation_north_arrow(
    location = "tl",
    style    = north_arrow_fancy_orienteering(),
    height   = unit(1.1, "cm"), width = unit(1.1, "cm")
  ) +
  annotate("text", x = 127.6, y = 36.5, label = "Republic\nof Korea",
           size = 4, fontface = "bold", colour = "grey30") +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  # Legend inside the overview panel
  theme_bw() +
  theme(
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 8),
    panel.grid.major  = element_line(colour = "grey90", linewidth = 0.3),
    panel.border      = element_rect(colour = "black", linewidth = 0.8),
    legend.position   = c(0.82, 0.25),
    legend.background = element_rect(fill = "white", colour = "grey60",
                                     linewidth = 0.4),
    legend.title      = element_text(size = 8, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.1, "cm")
  ) +
  guides(
    shape  = guide_legend(title = "Sampling Method", order = 1,
                          override.aes = list(size = 2.5, colour = "black")),
    colour = guide_legend(title = "Location",        order = 2,
                          override.aes = list(shape = 16, size = 2.5))
  )

save_map(fig_overview, "Figure_01_overview", width = 13, height = 14)

# =============================================================================
# 2. South Sea
# =============================================================================
ss <- filter(stations_sf, Location == "South Sea")

fig_south <- ggplot() +
  geom_sf(data = world, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_sf(data = ss, aes(colour = Location, shape = Method),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.8, 129.5), ylim = c(34.0, 35.5), expand = FALSE) +
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
  geom_sf(data = es, aes(colour = Location, shape = Method),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(128.8, 130.2), ylim = c(34.9, 39.0), expand = FALSE) +
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
  geom_sf(data = ys, aes(colour = Location, shape = Method),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.6, 127.8), ylim = c(34.0, 38.2), expand = FALSE) +
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
  geom_sf(data = jeju, aes(colour = Location, shape = Method),
          size = 2, stroke = 0.5, alpha = 0.9) +
  shared_colour +
  shared_shape +
  coord_sf(xlim = c(125.9, 127.3), ylim = c(32.95, 34.0), expand = FALSE) +
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
