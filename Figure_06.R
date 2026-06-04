# =============================================================================
# Figure 6. Phytoplankton and Zooplankton richness & top-20 ASV comparison
#
#   (a) Phytoplankton genus richness boxplot  → Figure_06a.tiff
#   (b) Top-20 phytoplankton total ASV bar   → Figure_06b.tiff
#   (c) Zooplankton genus richness boxplot   → Figure_06c.tiff
#   (d) Top-20 zooplankton total ASV bar     → Figure_06d.tiff
# =============================================================================
# install.packages(c("readxl","dplyr","ggplot2","scales","ggtext"))

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(ggtext)   # for element_markdown (bold italic axis labels)

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

phyto_rows <- which(tax_vec == "식물플랑크톤")
zoo_rows   <- which(tax_vec == "동물플랑크톤")
uuv_cols   <- which(methods == "Drone")
net_cols   <- which(methods == "Net")

# Colours
col_uuv   <- "#4472C4"
col_net   <- "#E05252"
col_phyto <- "#2E8B57"
col_zoo   <- "#E07B00"

# Helper: save
save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".tiff"), plot = p,
         width = w, height = h, units = "cm",
         dpi = 300, compression = "lzw")
  ggsave(paste0(name, ".pdf"),  plot = p,
         width = w, height = h, units = "cm")
  message("Saved: ", name)
}

# =============================================================================
# (a) Phytoplankton genus richness boxplot
# =============================================================================
phyto_rich <- data.frame(
  Method  = factor(c(rep("UUV", length(uuv_cols)),
                     rep("Net", length(net_cols))),
                   levels = c("UUV", "Net")),
  Richness = c(
    colSums(asv_matrix[phyto_rows, uuv_cols] > 0),
    colSums(asv_matrix[phyto_rows, net_cols]  > 0)
  )
)

fig_a <- ggplot(phyto_rich, aes(x = Method, y = Richness, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.8,
               linewidth = 0.5, width = 0.55) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(limits = c(10, 100),
                     breaks = seq(20, 100, 20),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Number of Phytoplankton Genus") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(size = 12, colour = "black"),
    axis.text.y     = element_text(size = 11),
    axis.title.y    = element_text(size = 12, face = "bold"),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(10, 12, 10, 10)
  )

save_fig(fig_a, "Figure_06a", 10, 13)

# =============================================================================
# (b) Top-20 phytoplankton total ASV — grouped horizontal bar
#     Order: descending by total (UUV + Net) ASV → top genus at top
# =============================================================================
phyto_names  <- genus_vec[phyto_rows]
phyto_uuv    <- rowSums(asv_matrix[phyto_rows, uuv_cols, drop = FALSE])
phyto_net    <- rowSums(asv_matrix[phyto_rows, net_cols,  drop = FALSE])
phyto_total  <- phyto_uuv + phyto_net

top20p_idx   <- order(phyto_total, decreasing = TRUE)[1:20]
top20p_order <- rev(phyto_names[top20p_idx])   # bottom→top for coord_flip

# Build data frame
top20p_df <- data.frame(
  Genus  = rep(phyto_names[top20p_idx], 2),
  Method = factor(rep(c("UUV", "Net"), each = 20),
                  levels = c("UUV", "Net")),
  ASV    = c(phyto_uuv[top20p_idx], phyto_net[top20p_idx])
) %>%
  mutate(
    # Bold italic markdown labels for ggtext
    GenusMD = paste0("***", Genus, "***"),
    GenusMD = factor(GenusMD,
                     levels = paste0("***", top20p_order, "***"))
  )

fig_b <- ggplot(top20p_df, aes(x = GenusMD, y = ASV, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c("UUV" = col_uuv, "Net" = col_net),
    name   = NULL,
    guide  = guide_legend(reverse = TRUE)   # Net on top in legend
  ) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 250000),
                     breaks = seq(0, 200000, 50000),
                     expand = expansion(mult = c(0, 0.03))) +
  coord_flip() +
  labs(x = NULL, y = "Total ASV") +
  theme_classic(base_size = 11) +
  theme(
    # Bold italic green labels via ggtext
    axis.text.y       = element_markdown(size = 9.5, colour = col_phyto,
                                         hjust = 1),
    axis.text.x       = element_text(size = 9, colour = "black"),
    axis.title.x      = element_text(size = 10.5, face = "bold"),
    axis.line         = element_line(linewidth = 0.4),
    legend.position   = c(0.80, 0.12),
    legend.background = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.4),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 9.5),
    plot.margin       = margin(10, 12, 10, 5)
  )

save_fig(fig_b, "Figure_06b", 18, 18)

# =============================================================================
# (c) Zooplankton genus richness boxplot
# =============================================================================
zoo_rich <- data.frame(
  Method  = factor(c(rep("UUV", length(uuv_cols)),
                     rep("Net", length(net_cols))),
                   levels = c("UUV", "Net")),
  Richness = c(
    colSums(asv_matrix[zoo_rows, uuv_cols] > 0),
    colSums(asv_matrix[zoo_rows, net_cols]  > 0)
  )
)

fig_c <- ggplot(zoo_rich, aes(x = Method, y = Richness, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1.8,
               linewidth = 0.5, width = 0.55) +
  scale_fill_manual(values = c("UUV" = col_uuv, "Net" = col_net)) +
  scale_y_continuous(limits = c(10, 65),
                     breaks = seq(10, 60, 10),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = NULL, y = "Number of Zooplankton Genus") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x     = element_text(size = 12, colour = "black"),
    axis.text.y     = element_text(size = 11),
    axis.title.y    = element_text(size = 12, face = "bold"),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(10, 12, 10, 10)
  )

save_fig(fig_c, "Figure_06c", 10, 13)

# =============================================================================
# (d) Top-20 zooplankton total ASV — grouped horizontal bar
#     Noctiluca now included as Zooplankton (rank 1 in Net)
# =============================================================================
zoo_names  <- genus_vec[zoo_rows]
zoo_uuv    <- rowSums(asv_matrix[zoo_rows, uuv_cols, drop = FALSE])
zoo_net    <- rowSums(asv_matrix[zoo_rows, net_cols,  drop = FALSE])
zoo_total  <- zoo_uuv + zoo_net

top20z_idx   <- order(zoo_total, decreasing = TRUE)[1:20]
top20z_order <- rev(zoo_names[top20z_idx])   # bottom→top

top20z_df <- data.frame(
  Genus  = rep(zoo_names[top20z_idx], 2),
  Method = factor(rep(c("UUV", "Net"), each = 20),
                  levels = c("UUV", "Net")),
  ASV    = c(zoo_uuv[top20z_idx], zoo_net[top20z_idx])
) %>%
  mutate(
    GenusMD = paste0("***", Genus, "***"),
    GenusMD = factor(GenusMD,
                     levels = paste0("***", top20z_order, "***"))
  )

fig_d <- ggplot(top20z_df, aes(x = GenusMD, y = ASV, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7, colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c("UUV" = col_uuv, "Net" = col_net),
    name   = NULL,
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_y_continuous(labels = comma,
                     limits = c(0, 175000),
                     breaks = seq(0, 160000, 40000),
                     expand = expansion(mult = c(0, 0.03))) +
  coord_flip() +
  labs(x = NULL, y = "Total ASV") +
  theme_classic(base_size = 11) +
  theme(
    # Bold italic orange labels via ggtext
    axis.text.y       = element_markdown(size = 9.5, colour = col_zoo,
                                         hjust = 1),
    axis.text.x       = element_text(size = 9, colour = "black"),
    axis.title.x      = element_text(size = 10.5, face = "bold"),
    axis.line         = element_line(linewidth = 0.4),
    legend.position   = c(0.80, 0.12),
    legend.background = element_rect(fill = "white", colour = "grey70",
                                     linewidth = 0.4),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 9.5),
    plot.margin       = margin(10, 12, 10, 5)
  )

save_fig(fig_d, "Figure_06d", 18, 18)

message("\nAll Figure 6 panels (a–d) saved successfully.")
