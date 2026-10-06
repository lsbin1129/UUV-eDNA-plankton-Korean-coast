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

reads_matrix <- raw[-(1:5), -(1:2)]
reads_matrix <- apply(reads_matrix, 2,
                    function(x) as.numeric(as.character(x)))
reads_matrix[is.na(reads_matrix)] <- 0
rownames(reads_matrix) <- genus_vec

# ── Guild labels (Phytoplankton / Zooplankton) are read directly from the
#    corrected "Taxonomic Group" column of Supplementary Data S1
#    (MetaK90_location.xlsx); no genus-level overrides are applied in code.

uuv_idx <- which(methods == "Drone")
net_idx <- which(methods == "Net")

# Per-method MEAN reads per profile
uuv_mean <- rowMeans(reads_matrix[, uuv_idx, drop = FALSE])
net_mean  <- rowMeans(reads_matrix[, net_idx, drop = FALSE])

# Log10 transform
uuv_log <- log10(uuv_mean + 1)
net_log  <- log10(net_mean  + 1)

# =============================================================================
# 1. Genus selection and order (top -> bottom)
#    Top 30 genera by mean reads per profile across all 125 profiles
#    (UUV and Net pooled), displayed as Phytoplankton first, then Zooplankton,
#    each block in descending order of mean reads. Guild membership is taken
#    from the "Taxonomic Group" column, so no genus is placed by hand.
# =============================================================================
N_TOP       <- 30
overall_mean <- setNames(rowMeans(reads_matrix), genus_vec)
top_genera   <- names(sort(overall_mean, decreasing = TRUE))[1:N_TOP]
top_tax      <- tax_vec[match(top_genera, genus_vec)]

genus_order_top_to_bottom <- c(top_genera[top_tax == "Phytoplankton"],
                               top_genera[top_tax == "Zooplankton"])

# Report differences from the genus set used in the previous figure version
previous_set <- c("Chaetoceros", "Thalassiosira", "Stephanodiscus", "Skeletonema",
  "Ulva", "Hemiaulus", "Rhizosolenia", "Urospora", "Lithodesmium", "Nitzschia",
  "Navicula", "Scytosiphon", "Pseudo-nitzschia", "Conticribra", "Gymnodinium",
  "Heterocapsa", "Noctiluca", "Ebria", "Polykrikos", "Neoturris", "Centropages",
  "Pseudocalanus", "Pontocythere", "Beroe", "Gyrodiniellum", "Strombidium",
  "Tisbe", "Tigriopus", "Corycaeus", "Parastrombidinopsis")
message("Top-30 genera (data-driven): ",
        paste(genus_order_top_to_bottom, collapse = ", "))
message("  Added vs previous version:   ",
        paste(setdiff(genus_order_top_to_bottom, previous_set), collapse = ", "))
message("  Dropped vs previous version: ",
        paste(setdiff(previous_set, genus_order_top_to_bottom), collapse = ", "))

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
               values_to = "log_reads") %>%
  mutate(
    Group  = ifelse(Tax == "Phytoplankton", "Phytoplankton", "Zooplankton"),
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
  ifelse(length(tax) > 0 && tax[1] == "Phytoplankton", col_phyto, col_zoo)
})

# =============================================================================
# 4. Plot
# =============================================================================
figure3 <- ggplot(heatmap_df,
                  aes(x = Method, y = Genus, fill = log_reads)) +
  geom_tile(colour = "grey65", linewidth = 0.2) +
  # Color gradient: cream → yellow → orange → red → dark red
  scale_fill_gradientn(
    colours = c("#FFFACD", "#FFD700", "#FF8C00", "#E03000", "#8B0000"),
    values  = rescale(c(0, 0.8, 1.6, 2.4, 3.2)),
    limits  = c(0, 3.6),
    oob     = squish,
    name    = expression(log[10]("Mean reads per profile" + 1)),
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
