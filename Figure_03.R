# =============================================================================
# Figure 3. Heatmap of read abundance for dominant plankton genera
#           by sampling method (UUV vs Net)
# =============================================================================
# install.packages(c("readxl", "dplyr", "ggplot2", "scales", "tidyr"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =============================================================================
# 0. Load data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx", col_names = FALSE)

methods   <- as.character(raw[1, 3:ncol(raw)])
genus_vec <- as.character(raw[-(1:5), 1][[1]])
tax_vec   <- as.character(raw[-(1:5), 2][[1]])

asv_matrix <- raw[-(1:5), -(1:2)]
asv_matrix <- apply(asv_matrix, 2,
                    function(x) as.numeric(as.character(x)))
asv_matrix[is.na(asv_matrix)] <- 0

# ── Reclassify Noctiluca as Zooplankton ──────────────────────────────────────
tax_vec[genus_vec == "Noctiluca"] <- "동물플랑크톤"

uuv_idx <- which(methods == "Drone")
net_idx <- which(methods == "Net")

# Per-method MEAN reads per sample (matches original code)
uuv_mean <- rowMeans(asv_matrix[, uuv_idx, drop = FALSE])
net_mean  <- rowMeans(asv_matrix[, net_idx, drop = FALSE])

# Log10 transform
uuv_log <- log10(uuv_mean + 1)
net_log  <- log10(net_mean  + 1)

# =============================================================================
# 1. Manual genus order (top → bottom): Phytoplankton first, Zooplankton after
#    Within each group: sorted by descending total ASV (from data)
# =============================================================================
genus_order_top_to_bottom <- c(
  # ── Phytoplankton (green) ──────────────────────────────────────────────────
  "Chaetoceros", "Thalassiosira", "Stephanodiscus", "Skeletonema",
  "Ulva", "Hemiaulus", "Rhizosolenia", "Urospora", "Lithodesmium",
  "Ebria", "Gyrodiniellum", "Nitzschia", "Navicula", "Scytosiphon",
  "Pseudo-nitzschia", "Conticribra", "Gymnodinium",
  # ── Zooplankton (orange) ──────────────────────────────────────────────────
  "Noctiluca", "Neoturris", "Centropages", "Pseudocalanus", "Pontocythere",
  "Heterocapsa", "Beroe", "Polykrikos", "Strombidium",
  "Tisbe", "Tigriopus", "Corycaeus", "Parastrombidinopsis"
)

# =============================================================================
# 2. Build long data frame
# =============================================================================
heatmap_df <- data.frame(
  Genus   = genus_vec,
  UUV     = uuv_log,
  Net     = net_log,
  Tax     = tax_vec
) %>%
  filter(Genus %in% genus_order_top_to_bottom) %>%
  pivot_longer(cols = c("UUV", "Net"),
               names_to  = "Method",
               values_to = "log_ASV") %>%
  mutate(
    Group  = ifelse(Tax == "식물플랑크톤", "Phytoplankton", "Zooplankton"),
    # Factor: bottom → top for ggplot2 y-axis (drawn bottom-up)
    Genus  = factor(Genus, levels = rev(genus_order_top_to_bottom)),
    Method = factor(Method, levels = c("UUV", "Net"))
  )

# =============================================================================
# 3. Y-axis label colours
# =============================================================================
col_phyto <- "#2E8B57"   # green
col_zoo   <- "#E07B00"   # orange

level_order <- levels(heatmap_df$Genus)   # bottom → top

label_colours <- sapply(level_order, function(g) {
  tax <- unique(heatmap_df$Tax[heatmap_df$Genus == g])
  ifelse(length(tax) > 0 && tax[1] == "식물플랑크톤", col_phyto, col_zoo)
})

# =============================================================================
# 4. Plot
# =============================================================================
figure3 <- ggplot(heatmap_df,
                  aes(x = Method, y = Genus, fill = log_ASV)) +
  geom_tile(colour = "grey65", linewidth = 0.2) +
  # Color gradient: cream → yellow → orange → red → dark red
  scale_fill_gradientn(
    colours = c("#FFFACD", "#FFD700", "#FF8C00", "#E03000", "#8B0000"),
    values  = rescale(c(0, 0.8, 1.6, 2.4, 3.2)),
    limits  = c(0, 3.6),
    oob     = squish,
    name    = expression(log[10](Reads + 1)),
    guide   = guide_colorbar(
      barheight    = unit(14, "cm"),
      barwidth     = unit(0.6, "cm"),
      ticks.colour = "grey40",
      frame.colour = "grey40",
      label.theme  = element_text(size = 9)
    )
  ) +
  # Vertical divider between UUV and Net columns
  geom_vline(xintercept = 1.5, colour = "grey40", linewidth = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Sampling Method", y = "Genus") +
  theme_classic(base_size = 11) +
  theme(
    # Axis titles
    axis.title.x  = element_text(size = 13, face = "bold",
                                 margin = margin(t = 10)),
    axis.title.y  = element_text(size = 13, face = "bold",
                                 margin = margin(r = 10)),
    # X-axis text (UUV, Net)
    axis.text.x   = element_text(size = 12, colour = "black"),
    # Y-axis text: italic + group colour
    axis.text.y   = element_text(size = 9.5, face = "italic",
                                 colour = label_colours, hjust = 1),
    # Remove axis lines/ticks on tile sides
    axis.line     = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.ticks.x  = element_line(colour = "grey50", linewidth = 0.3),
    # Legend
    legend.title  = element_text(size = 10, angle = 90,
                                 hjust = 0.5, vjust = 0.5),
    legend.text   = element_text(size = 9),
    legend.position = "right",
    plot.margin   = margin(15, 10, 10, 10)
  )

# =============================================================================
# 5. Save
# =============================================================================
ggsave("Figure_03.tiff", plot = figure3,
       width = 18, height = 24, units = "cm",
       dpi = 300, compression = "lzw")

ggsave("Figure_03.pdf", plot = figure3,
       width = 18, height = 24, units = "cm")

message("Figure 3 saved: Figure_03.tiff / Figure_03.pdf")
