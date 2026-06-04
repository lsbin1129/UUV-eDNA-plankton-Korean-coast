# =============================================================================
# Figure 4. Macro-regional benchmarking of plankton generic richness
#
#   Output 1 — Figure_04a.tiff      : (a) bar chart of mean total richness
#   Output 2 — Figure_04bc.tiff     : (b) phytoplankton boxplot
#                                     (c) zooplankton boxplot
#   Output 3 — Figure_04dh.tiff     : (d)–(h) scatter plots per region
#
# =============================================================================
# install.packages(c("readxl","dplyr","ggplot2","cowplot","scales"))

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)

# =============================================================================
# 0. Load and prepare data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx", col_names = FALSE)

methods   <- as.character(raw[1, 3:ncol(raw)])
locations <- as.character(raw[2, 3:ncol(raw)])
genus_vec <- as.character(raw[-(1:5), 1][[1]])
tax_vec   <- as.character(raw[-(1:5), 2][[1]])

asv_matrix <- raw[-(1:5), -(1:2)]
asv_matrix <- apply(asv_matrix, 2,
                    function(x) as.numeric(as.character(x)))
asv_matrix[is.na(asv_matrix)] <- 0

# ── Reclassify Noctiluca as Zooplankton ──────────────────────────────────────
tax_vec[genus_vec == "Noctiluca"] <- "동물플랑크톤"

phyto_rows <- which(tax_vec == "식물플랑크톤")
zoo_rows   <- which(tax_vec == "동물플랑크톤")

# ── Per-sample richness ──────────────────────────────────────────────────────
n_samples <- ncol(asv_matrix)

sample_df <- data.frame(
  Method    = methods,
  Location  = locations,
  TotalRich = sapply(1:n_samples, function(j) sum(asv_matrix[, j] > 0)),
  PhytoRich = sapply(1:n_samples, function(j) sum(asv_matrix[phyto_rows, j] > 0)),
  ZooRich   = sapply(1:n_samples, function(j) sum(asv_matrix[zoo_rows,   j] > 0))
)

# Factor levels
loc_levels <- c("South Sea", "East Sea", "Yellow Sea", "Jeju Sea", "Seogwipo Sea")
loc_labels <- c("South\nSea", "East\nSea", "Yellow\nSea", "Jeju\nSea", "Seogwipo\nSea")

sample_df <- sample_df %>%
  mutate(
    Location = factor(Location, levels = loc_levels),
    Method   = factor(Method,   levels = c("Drone", "Net"),
                                labels = c("UUV", "Net"))
  )

# Colour palette
col_uuv <- "#4472C4"
col_net <- "#E05252"

# =============================================================================
# (a) Bar chart — mean total richness per Location × Method
# =============================================================================
bar_df <- sample_df %>%
  group_by(Location, Method) %>%
  summarise(MeanRich = mean(TotalRich), .groups = "drop") %>%
  mutate(LocLabel = factor(Location,
                           levels = loc_levels,
                           labels = loc_labels))

fig_a <- ggplot(bar_df, aes(x = LocLabel, y = MeanRich, fill = Method)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 130),
                     breaks = seq(0, 120, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Richness") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(size = 10, colour = "black"),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 12, face = "bold"),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    legend.position = c(0.15, 0.88),
    legend.background = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.4),
    legend.key.size = unit(0.5, "cm"),
    legend.text     = element_text(size = 10),
    plot.margin     = margin(10, 10, 5, 10)
  )

# Save (a)
ggsave("Figure_04a.tiff", plot = fig_a,
       width = 14, height = 10, units = "cm",
       dpi = 300, compression = "lzw")
ggsave("Figure_04a.pdf",  plot = fig_a,
       width = 14, height = 10, units = "cm")
message("Saved: Figure_04a")

# =============================================================================
# Helper: boxplot function for (b) and (c)
# =============================================================================
make_boxplot <- function(data, yvar, ytitle, ylim_max) {
  data$yval <- data[[yvar]]
  data <- data %>%
    mutate(LocLabel = factor(Location,
                             levels = loc_levels,
                             labels = loc_labels))

  ggplot(data, aes(x = LocLabel, y = yval, fill = Method)) +
    geom_boxplot(
      position     = position_dodge(width = 0.75),
      width        = 0.65,
      outlier.shape = 1,
      outlier.size  = 1.5,
      outlier.stroke = 0.5,
      linewidth    = 0.45,
      colour       = "black"
    ) +
    scale_fill_manual(
      values = c("UUV" = col_uuv, "Net" = col_net),
      name   = NULL,
      labels = c("UUV" = "Uncrewed underwater vehicle (UUV)", "Net" = "Net")
    ) +
    scale_y_continuous(
      limits = c(0, ylim_max),
      breaks = seq(0, ylim_max, 20),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(x = NULL, y = ytitle) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x       = element_text(size = 9.5, colour = "black"),
      axis.text.y       = element_text(size = 10),
      axis.title.y      = element_text(size = 11, face = "bold"),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      legend.position   = c(0.28, 0.10),
      legend.background = element_rect(fill = "white", colour = "grey70",
                                       linewidth = 0.4),
      legend.key.size   = unit(0.5, "cm"),
      legend.text       = element_text(size = 9),
      plot.margin       = margin(10, 10, 5, 10)
    )
}

# =============================================================================
# (b) Phytoplankton richness boxplot
# (c) Zooplankton richness boxplot
# =============================================================================
fig_b <- make_boxplot(sample_df, "PhytoRich",
                      "Phytoplankton Genus Richness", 100) +
  ggtitle("Phytoplankton Richness by Location & Method") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

fig_c <- make_boxplot(sample_df, "ZooRich",
                      "Zooplankton Genus Richness", 60) +
  ggtitle("Zooplankton Richness by Location & Method") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

fig_bc <- plot_grid(fig_b, fig_c,
                    ncol   = 2,
                    labels = c("(b)", "(c)"),
                    label_size = 13, label_fontface = "bold")

ggsave("Figure_04bc.tiff", plot = fig_bc,
       width = 28, height = 12, units = "cm",
       dpi = 300, compression = "lzw")
ggsave("Figure_04bc.pdf",  plot = fig_bc,
       width = 28, height = 12, units = "cm")
message("Saved: Figure_04bc")

# =============================================================================
# (d)–(h) Scatter plots: Phytoplankton vs Zooplankton richness per region
# =============================================================================
make_scatter <- function(data, region, xlim = c(0, 100), ylim_max = 70) {

  sub <- data %>% filter(Location == region)

  # Global axis limits based on all regions for consistency
  ggplot(sub, aes(x = PhytoRich, y = ZooRich, colour = Method)) +
    # Diagonal reference line (y = x)
    geom_abline(slope = 1, intercept = 0,
                colour = "grey60", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_colour_manual(values = c("UUV" = col_uuv, "Net" = col_net),
                        name = NULL) +
    scale_x_continuous(limits = xlim,
                       breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_y_continuous(limits = c(0, ylim_max),
                       breaks = seq(0, ylim_max, 10),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "Phytoplankton Richness",
         y = "Zooplankton Richness",
         title = region) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.text       = element_text(size = 8.5, colour = "black"),
      axis.title      = element_text(size = 9.5),
      axis.line       = element_line(colour = "black", linewidth = 0.4),
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", colour = "grey70",
                                       linewidth = 0.3),
      legend.key.size = unit(0.35, "cm"),
      legend.text     = element_text(size = 8),
      plot.margin     = margin(8, 8, 5, 8)
    )
}

fig_d <- make_scatter(sample_df, "South Sea",  ylim_max = 70)
fig_e <- make_scatter(sample_df, "East Sea",   ylim_max = 60)
fig_f <- make_scatter(sample_df, "Yellow Sea", ylim_max = 60)
fig_g <- make_scatter(sample_df, "Jeju Sea",   ylim_max = 60)
fig_h <- make_scatter(sample_df, "Seogwipo Sea", ylim_max = 50)

fig_dh <- plot_grid(
  fig_d, fig_e, fig_f, fig_g, fig_h,
  ncol        = 5,
  labels      = c("(d)", "(e)", "(f)", "(g)", "(h)"),
  label_size  = 13, label_fontface = "bold"
)

ggsave("Figure_04dh.tiff", plot = fig_dh,
       width = 40, height = 10, units = "cm",
       dpi = 300, compression = "lzw")
ggsave("Figure_04dh.pdf",  plot = fig_dh,
       width = 40, height = 10, units = "cm")
message("Saved: Figure_04dh")

message("\nAll Figure 4 panels saved successfully.")
