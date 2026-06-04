# =============================================================================
# Figure 7. Relative read abundance of plankton genera across five Korean
#           coastal sea regions
#
# Visualization: stacked bar chart (relative abundance, %)
#
# Output: Figure_07.tiff and Figure_07.pdf
# =============================================================================
# Required packages:
#   install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "scales"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =============================================================================
# 1. Load and prepare data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx",
                  col_names = FALSE, .name_repair = "minimal")

# Extract metadata and count matrix
locations  <- as.character(unlist(raw[2, 3:ncol(raw)]))
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
# 2. Aggregate read counts by sea region
# =============================================================================
sea_regions <- c("Yellow Sea", "South Sea", "East Sea",
                 "Jeju Sea", "Seogwipo Sea")

# Sum reads per genus within each sea region
sea_totals <- sapply(sea_regions, function(sea) {
  idx <- which(locations == sea)
  rowSums(counts_mat[, idx, drop = FALSE])
})

sea_df        <- as.data.frame(sea_totals)
sea_df$Genus  <- genus_vec
sea_df$TaxGroup <- ifelse(tax_vec == "식물플랑크톤",
                          "Phytoplankton", "Zooplankton")

# =============================================================================
# 3. Calculate relative abundance (%) per sea region
# =============================================================================
rel_df <- sea_df
for (sea in sea_regions) {
  total       <- sum(rel_df[[sea]])
  rel_df[[sea]] <- rel_df[[sea]] / total * 100
}

# Mean relative abundance across all five sea regions
rel_df$MeanRelAbund <- rowMeans(rel_df[, sea_regions])

# Sort genera by mean relative abundance (descending)
rel_df <- rel_df[order(-rel_df$MeanRelAbund), ]
rownames(rel_df) <- NULL

# =============================================================================
# 4. Classify genera: top 30 vs. remainder
# =============================================================================
N_TOP        <- 30
top_genera   <- rel_df$Genus[1:N_TOP]
rest_genera  <- rel_df$Genus[(N_TOP + 1):nrow(rel_df)]
n_rest       <- length(rest_genera)

# =============================================================================
# 5. Assign colour palettes
# =============================================================================
# Top 30: distinct qualitative colours
top_colours <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363"
)
names(top_colours) <- top_genera

# Remainder: muted, desaturated tones (not shown in legend)
set.seed(42)
rest_colours <- colorRampPalette(
  c("#d9d9d9", "#bdbdbd", "#c6dbef", "#dadaeb",
    "#d9f0a3", "#fdd0a2", "#fde0dd")
)(n_rest)
rest_colours <- rest_colours[sample(n_rest)]
names(rest_colours) <- rest_genera

all_colours <- c(top_colours, rest_colours)

# =============================================================================
# 6. Reshape data to long format for ggplot2
# =============================================================================
long_df <- rel_df %>%
  pivot_longer(cols      = all_of(sea_regions),
               names_to  = "Sea",
               values_to = "RelAbund") %>%
  mutate(
    Sea   = factor(Sea, levels = sea_regions),
    Genus = factor(Genus, levels = rev(rel_df$Genus))  # abundant stacked at bottom
  )

# Cumulative sum of top-30 per sea region (for dashed boundary line)
top30_boundary <- rel_df %>%
  filter(Genus %in% top_genera) %>%
  summarise(across(all_of(sea_regions), sum)) %>%
  pivot_longer(everything(),
               names_to  = "Sea",
               values_to = "Top30Sum") %>%
  mutate(Sea = factor(Sea, levels = sea_regions))

# =============================================================================
# 7. Build italic legend labels for top-30 genera
# =============================================================================
italic_labels <- setNames(
  paste0("italic('", top_genera, "')"),
  top_genera
)

# =============================================================================
# 8. Plot
# =============================================================================
figure7 <- ggplot(long_df, aes(x = Sea, y = RelAbund, fill = Genus)) +

  # Stacked bar
  geom_col(width = 0.70, colour = NA) +

  # Dashed boundary line between top-30 and remainder
  geom_segment(
    data        = top30_boundary,
    aes(x       = as.numeric(Sea) - 0.32,
        xend    = as.numeric(Sea) + 0.32,
        y       = Top30Sum,
        yend    = Top30Sum),
    inherit.aes = FALSE,
    colour      = "#333333",
    linewidth   = 0.65,
    linetype    = "dashed"
  ) +

  # Colour scale: only top-30 labelled in legend
  scale_fill_manual(
    values  = all_colours,
    breaks  = top_genera,
    labels  = parse(text = italic_labels[top_genera]),
    guide   = guide_legend(
      title         = paste0("Top ", N_TOP, " genera"),
      title.theme   = element_text(size = 9.5, face = "bold"),
      label.theme   = element_text(size = 8.5, face = "italic"),
      keywidth      = unit(0.8, "cm"),
      keyheight     = unit(0.5, "cm"),
      ncol          = 1,
      byrow         = TRUE,
      override.aes  = list(colour = "white", linewidth = 0.3)
    )
  ) +

  # Axes
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.01)),
    labels = function(x) paste0(x, "%")
  ) +
  scale_x_discrete(
    labels = c(
      "Yellow Sea"   = "Yellow\nSea",
      "South Sea"    = "South\nSea",
      "East Sea"     = "East\nSea",
      "Jeju Sea"     = "Jeju\nSea",
      "Seogwipo Sea" = "Seogwipo\nSea"
    )
  ) +

  labs(x = NULL, y = "Relative Abundance (%)") +

  theme_classic(base_size = 12) +
  theme(
    axis.text.x        = element_text(size = 11, face = "bold",
                                      colour = "black"),
    axis.text.y        = element_text(size = 10),
    axis.title.y       = element_text(size = 11, face = "bold",
                                      margin = margin(r = 8)),
    axis.line.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    legend.position    = "right",
    legend.margin      = margin(l = 10),
    legend.title       = element_text(size = 9.5, face = "bold"),
    legend.text        = element_text(size = 8.5, face = "italic"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    plot.margin        = margin(t = 10, r = 5, b = 10, l = 10)
  )

# =============================================================================
# 9. Save output
# =============================================================================
ggsave("Figure_07.tiff", plot = figure7,
       width = 36, height = 20, units = "cm",
       dpi = 300, bg = "white", compression = "lzw")

ggsave("Figure_07.pdf", plot = figure7,
       width = 36, height = 20, units = "cm",
       bg = "white")

message("Figure 7 saved: Figure_07.tiff / Figure_07.pdf")
