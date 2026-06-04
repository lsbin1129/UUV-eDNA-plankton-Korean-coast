# =============================================================================
# Figure 5. Multidimensional benchmarking of plankton assay efficiencies
#
#   (a) Total ASV boxplot
#   (b) Genus richness boxplot
#   (c) Phytoplankton / Zooplankton read ratio stacked bar
#   (d) Shannon diversity boxplot
#   (e) Shannon diversity bar by group (Phyto / Zoo)
#   (f) Method-exclusive genus counts bar
#
# =============================================================================
# install.packages(c("readxl","dplyr","ggplot2","scales","tidyr"))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =============================================================================
# 0. Load and prepare data
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

phyto_rows <- which(tax_vec == "식물플랑크톤")
zoo_rows   <- which(tax_vec == "동물플랑크톤")
uuv_cols   <- which(methods == "Drone")
net_cols   <- which(methods == "Net")

# Colours
col_uuv   <- "#4472C4"
col_net   <- "#E05252"
col_phyto <- "#2E8B57"
col_zoo   <- "#E07B00"

# Shannon diversity function
shannon_H <- function(v) {
  v <- v[v > 0]
  if (length(v) == 0) return(0)
  p <- v / sum(v)
  -sum(p * log(p))
}

# =============================================================================
# Helper: save figure
# =============================================================================
save_fig <- function(plot_obj, name, width, height) {
  ggsave(paste0(name, ".tiff"), plot = plot_obj,
         width = width, height = height, units = "cm",
         dpi = 300, compression = "lzw")
  ggsave(paste0(name, ".pdf"),  plot = plot_obj,
         width = width, height = height, units = "cm")
  message("Saved: ", name)
}

# Common boxplot theme
box_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 11, colour = "black"),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.line    = element_line(colour = "black", linewidth = 0.5),
    plot.margin  = margin(10, 15, 10, 10)
  )

# =============================================================================
# (a) Total ASV per sample — boxplot
# =============================================================================
asv_df <- data.frame(
  Method = factor(c(rep("UUV", length(uuv_cols)),
                    rep("Net", length(net_cols))),
                  levels = c("UUV", "Net")),
  ASV    = c(colSums(asv_matrix[, uuv_cols]),
             colSums(asv_matrix[, net_cols]))
)

mean_asv <- asv_df %>%
  group_by(Method) %>%
  summarise(m = round(mean(ASV)), .groups = "drop")

fig_a <- ggplot(asv_df, aes(x = Method, y = ASV, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.8,
               linewidth = 0.5, width = 0.55) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 60000),
                     breaks = seq(0, 60000, 10000),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Total ASV") +
  box_theme +
  theme(legend.position = "none")

save_fig(fig_a, "Figure_05a", 10, 12)

# =============================================================================
# (b) Total genus richness per sample — boxplot
# =============================================================================
rich_df <- data.frame(
  Method  = factor(c(rep("UUV", length(uuv_cols)),
                     rep("Net", length(net_cols))),
                   levels = c("UUV", "Net")),
  Richness = c(colSums(asv_matrix[, uuv_cols] > 0),
               colSums(asv_matrix[, net_cols] > 0))
)

mean_rich <- rich_df %>%
  group_by(Method) %>%
  summarise(m = round(mean(Richness), 1), .groups = "drop")

fig_b <- ggplot(rich_df, aes(x = Method, y = Richness, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.8,
               linewidth = 0.5, width = 0.55) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(limits = c(30, 145),
                     breaks = seq(40, 140, 20),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Number of Detected Genus") +
  box_theme +
  theme(legend.position = "none")

save_fig(fig_b, "Figure_05b", 10, 12)

# =============================================================================
# (c) Phytoplankton / Zooplankton read ratio — stacked bar
# =============================================================================
uuv_phyto <- sum(asv_matrix[phyto_rows, uuv_cols])
uuv_zoo   <- sum(asv_matrix[zoo_rows,   uuv_cols])
net_phyto <- sum(asv_matrix[phyto_rows, net_cols])
net_zoo   <- sum(asv_matrix[zoo_rows,   net_cols])

ratio_df <- data.frame(
  Method = factor(c("UUV", "UUV", "Net", "Net"),
                  levels = c("UUV", "Net")),
  Group  = factor(c("Phytoplankton", "Zooplankton",
                    "Phytoplankton", "Zooplankton"),
                  levels = c("Zooplankton", "Phytoplankton")),
  Pct = c(
    uuv_phyto / (uuv_phyto + uuv_zoo) * 100,
    uuv_zoo   / (uuv_phyto + uuv_zoo) * 100,
    net_phyto / (net_phyto + net_zoo)  * 100,
    net_zoo   / (net_phyto + net_zoo)  * 100
  )
)

fig_c <- ggplot(ratio_df, aes(x = Method, y = Pct, fill = Group)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  # Boundary line between Phyto and Zoo
  geom_hline(yintercept = c(
    uuv_phyto / (uuv_phyto + uuv_zoo) * 100,
    net_phyto / (net_phyto + net_zoo)  * 100
  ), colour = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Phytoplankton" = col_phyto, "Zooplankton" = col_zoo),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Ratio (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x       = element_text(size = 11, colour = "black"),
    axis.text.y       = element_text(size = 10),
    axis.title.y      = element_text(size = 11, face = "bold"),
    axis.line         = element_line(colour = "black", linewidth = 0.5),
    legend.position   = c(0.82, 0.50),
    legend.background = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.4),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 9),
    plot.margin       = margin(10, 15, 10, 10)
  )

save_fig(fig_c, "Figure_05c", 10, 12)

# =============================================================================
# (d) Shannon diversity — boxplot
# =============================================================================
shan_df <- data.frame(
  Method  = factor(c(rep("UUV", length(uuv_cols)),
                     rep("Net", length(net_cols))),
                   levels = c("UUV", "Net")),
  Shannon = c(
    sapply(uuv_cols, function(j) shannon_H(asv_matrix[, j])),
    sapply(net_cols,  function(j) shannon_H(asv_matrix[, j]))
  )
)

mean_shan <- shan_df %>%
  group_by(Method) %>%
  summarise(m = round(mean(Shannon), 3), .groups = "drop")

fig_d <- ggplot(shan_df, aes(x = Method, y = Shannon, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.8,
               linewidth = 0.5, width = 0.55) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(limits = c(0.8, 4.2),
                     breaks = seq(1.0, 4.0, 0.5),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Shannon Diversity Index") +
  box_theme +
  theme(legend.position = "none")

save_fig(fig_d, "Figure_05d", 10, 12)

# =============================================================================
# (e) Shannon diversity by group (Phyto / Zoo) — grouped bar
# =============================================================================
shan_grp <- data.frame(
  Group  = factor(rep(c("Phytoplankton", "Zooplankton"), each = 2),
                  levels = c("Phytoplankton", "Zooplankton")),
  Method = factor(rep(c("UUV", "Net"), 2),
                  levels = c("UUV", "Net")),
  Shannon = c(
    mean(sapply(uuv_cols, function(j) shannon_H(asv_matrix[phyto_rows, j]))),
    mean(sapply(net_cols,  function(j) shannon_H(asv_matrix[phyto_rows, j]))),
    mean(sapply(uuv_cols, function(j) shannon_H(asv_matrix[zoo_rows, j]))),
    mean(sapply(net_cols,  function(j) shannon_H(asv_matrix[zoo_rows, j])))
  )
)

fig_e <- ggplot(shan_grp, aes(x = Group, y = Shannon, fill = Method)) +
  geom_col(position = position_dodge(width = 0.65),
           width = 0.6, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("UUV" = col_uuv, "Net" = col_net),
    name   = NULL,
    guide  = guide_legend(keywidth = unit(0.5, "cm"),
                          keyheight = unit(0.45, "cm"))
  ) +
  scale_y_continuous(limits = c(0, 2.6),
                     breaks = seq(0, 2.5, 0.5),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Shannon Diversity Index") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x       = element_text(size = 11, colour = "black"),
    axis.text.y       = element_text(size = 10),
    axis.title.y      = element_text(size = 11, face = "bold"),
    axis.line         = element_line(colour = "black", linewidth = 0.5),
    legend.position   = c(0.88, 0.85),
    legend.background = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.4),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 10),
    plot.margin       = margin(10, 15, 10, 10)
  )

save_fig(fig_e, "Figure_05e", 12, 12)

# =============================================================================
# (f) Method-exclusive genera counts — bar
# =============================================================================
uuv_det <- genus_vec[rowSums(asv_matrix[, uuv_cols] > 0) > 0]
net_det  <- genus_vec[rowSums(asv_matrix[, net_cols] > 0) > 0]
uuv_only <- setdiff(uuv_det, net_det)
net_only  <- setdiff(net_det, uuv_det)

tax_map <- setNames(tax_vec, genus_vec)

excl_df <- data.frame(
  Category = factor(
    c("Phytoplankton\nUUV only", "Phytoplankton\nNet only",
      "Zooplankton\nUUV only",   "Zooplankton\nNet only"),
    levels = c("Phytoplankton\nUUV only", "Phytoplankton\nNet only",
               "Zooplankton\nUUV only",   "Zooplankton\nNet only")
  ),
  Count = c(
    sum(tax_map[uuv_only] == "식물플랑크톤"),  # Phyto UUV only
    sum(tax_map[net_only]  == "식물플랑크톤"),  # Phyto Net only
    sum(tax_map[uuv_only] == "동물플랑크톤"),  # Zoo UUV only
    sum(tax_map[net_only]  == "동물플랑크톤")   # Zoo Net only
  ),
  Method = factor(c("UUV", "Net", "UUV", "Net"),
                  levels = c("UUV", "Net"))
)

fig_f <- ggplot(excl_df, aes(x = Category, y = Count, fill = Method)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(limits = c(0, 58),
                     breaks = seq(0, 50, 10),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "Number of Method-Specific Genus") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x     = element_text(size = 9.5, colour = "black"),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 10.5, face = "bold"),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(10, 15, 10, 10)
  )

save_fig(fig_f, "Figure_05f", 14, 12)

message("\nAll Figure 5 panels (a–f) saved successfully.")
