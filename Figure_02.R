# =============================================================================
# Figure 2. Comparison of plankton genus detection between UUV and Net
#   (a) Venn diagram   (b) Bar chart   (c) Top-20 UUV   (d) Top-20 Net
#
# =============================================================================
# install.packages(c("readxl","dplyr","ggplot2","cowplot","scales"))

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
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

uuv_idx    <- which(methods == "Drone")
net_idx    <- which(methods == "Net")
uuv_totals <- setNames(rowSums(asv_matrix[, uuv_idx, drop = FALSE]), genus_vec)
net_totals <- setNames(rowSums(asv_matrix[, net_idx, drop = FALSE]), genus_vec)

uuv_det <- genus_vec[uuv_totals > 0]
net_det <- genus_vec[net_totals > 0]

n_uuv  <- length(setdiff(uuv_det, net_det))   # 82
n_net  <- length(setdiff(net_det, uuv_det))   # 78
n_both <- length(intersect(uuv_det, net_det))  # 436
n_total <- n_uuv + n_net + n_both

pct_uuv  <- round(n_uuv  / n_total * 100, 1)
pct_net  <- round(n_net  / n_total * 100, 1)
pct_both <- round(n_both / n_total * 100, 1)

# Colour palette
col_uuv   <- "#4472C4"
col_net   <- "#E05252"
col_both  <- "#9B59B6"
col_phyto <- "#2E8B57"
col_zoo   <- "#E07B00"

# =============================================================================
# (a) Venn diagram — ellipses built from parametric equations
# =============================================================================

# Helper: generate ellipse outline points
make_ellipse <- function(x0, y0, a, b, n = 300) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + a * cos(t),
             y = y0 + b * sin(t))
}

# Left ellipse (UUV): centre (-0.5, 0), a=1.4, b=0.85
# Right ellipse (Net): centre ( 0.5, 0), a=1.4, b=0.85
ell_uuv <- make_ellipse(-0.5, 0, 1.4, 0.85)
ell_net <- make_ellipse( 0.5, 0, 1.4, 0.85)

fig_a <- ggplot() +
  # Filled ellipses (drawn as polygons)
  geom_polygon(data = ell_uuv,
               aes(x = x, y = y),
               fill = "#AEC6E8", colour = col_uuv,
               linewidth = 1.0, alpha = 0.6) +
  geom_polygon(data = ell_net,
               aes(x = x, y = y),
               fill = "#F4AAAA", colour = col_net,
               linewidth = 1.0, alpha = 0.6) +
  # Numbers
  annotate("text", x = -1.10, y =  0,
           label = n_uuv,
           size = 10, fontface = "bold", colour = col_uuv) +
  annotate("text", x =  0.00, y =  0,
           label = n_both,
           size = 10, fontface = "bold", colour = col_both) +
  annotate("text", x =  1.10, y =  0,
           label = n_net,
           size = 10, fontface = "bold", colour = col_net) +
  # Set labels
  annotate("text", x = -0.90, y =  1.05,
           label = "UUV Only",
           size = 5.5, fontface = "bold", colour = col_uuv) +
  annotate("text", x =  0.90, y =  1.05,
           label = "Net Only",
           size = 5.5, fontface = "bold", colour = col_net) +
  annotate("text", x =  0.00, y = -1.10,
           label = "Both Methods",
           size = 5.5, fontface = "bold", colour = col_both) +
  coord_fixed(xlim = c(-2.1, 2.1), ylim = c(-1.4, 1.4)) +
  theme_void() +
  theme(plot.margin = margin(10, 5, 15, 5))

# =============================================================================
# (b) Bar chart
# =============================================================================
bar_df <- data.frame(
  Category = factor(c("UUV Only", "Both Methods", "Net Only"),
                    levels = c("UUV Only", "Both Methods", "Net Only")),
  Count    = c(n_uuv, n_both, n_net),
  Pct      = c(pct_uuv, pct_both, pct_net),
  Fill     = c(col_uuv, col_both, col_net)
)

fig_b <- ggplot(bar_df, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = Count),
            vjust = -0.55, fontface = "bold",
            size = 4.5, colour = "black") +
  geom_text(aes(label = paste0(Pct, "%"), y = Count / 2),
            fontface = "bold", size = 4.2, colour = "white") +
  scale_fill_manual(values = setNames(as.character(bar_df$Fill),
                                      as.character(bar_df$Category))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = seq(0, 500, 100),
                     limits = c(0, 510)) +
  labs(y = "Number of Genus", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(size = 10, face = "bold",
                                   colour = c(col_uuv, col_both, col_net)),
    axis.text.y     = element_text(size = 9),
    axis.title.y    = element_text(size = 10),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(5, 15, 5, 5)
  )

# =============================================================================
# Helper: horizontal bar plot for top-20 genera
# =============================================================================
make_top20_bar <- function(totals, tax_vec, genus_vec,
                           col_phyto, col_zoo) {
  top20 <- data.frame(
    Genus = names(sort(totals, decreasing = TRUE)[1:20]),
    ASV   = as.numeric(sort(totals, decreasing = TRUE)[1:20])
  ) %>%
    mutate(
      Tax   = tax_vec[match(Genus, genus_vec)],
      Group = ifelse(Tax == "식물플랑크톤", "Phytoplankton", "Zooplankton"),
      Rank  = 1:20,
      Label = paste0(Rank, ". ", Genus)
    ) %>%
    arrange(ASV) %>%
    mutate(Label = factor(Label, levels = Label))

  # Axis label colours (top of factor = rank 1 at top of flipped plot)
  label_cols <- ifelse(
    top20$Group[match(levels(top20$Label), top20$Label)] == "Phytoplankton",
    col_phyto, col_zoo
  )

  ggplot(top20, aes(x = Label, y = ASV, fill = Group)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
    scale_fill_manual(
      values = c("Phytoplankton" = col_phyto, "Zooplankton" = col_zoo),
      name   = NULL,
      guide  = guide_legend(
        keywidth  = unit(0.5, "cm"),
        keyheight = unit(0.4, "cm")
      )
    ) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.03))) +
    coord_flip() +
    labs(x = NULL, y = "Total ASV") +
    theme_classic() +
    theme(
      axis.text.y       = element_text(size = 9, face = "italic",
                                       colour = label_cols, hjust = 1),
      axis.text.x       = element_text(size = 8),
      axis.title.x      = element_text(size = 9.5),
      axis.line         = element_line(linewidth = 0.4),
      legend.position   = c(0.80, 0.10),
      legend.background = element_rect(fill     = "white",
                                       colour   = "grey75",
                                       linewidth = 0.4),
      legend.text       = element_text(size = 9),
      plot.margin       = margin(5, 12, 5, 5)
    )
}

# =============================================================================
# (c) Top-20 UUV   (d) Top-20 Net
# =============================================================================
fig_c <- make_top20_bar(uuv_totals, tax_vec, genus_vec, col_phyto, col_zoo)
fig_d <- make_top20_bar(net_totals, tax_vec, genus_vec, col_phyto, col_zoo)

# =============================================================================
# Assemble
# =============================================================================
top_row <- plot_grid(
  fig_a, fig_b,
  ncol        = 2,
  rel_widths  = c(1.1, 1),
  labels      = c("(a)", "(b)"),
  label_size  = 14, label_fontface = "bold"
)

bottom_row <- plot_grid(
  fig_c, fig_d,
  ncol        = 2,
  labels      = c("(c)", "(d)"),
  label_size  = 14, label_fontface = "bold",
  align       = "h", axis = "tb"
)

figure2 <- plot_grid(
  top_row, bottom_row,
  nrow        = 2,
  rel_heights = c(1, 1.8)
)

# =============================================================================
# Save
# =============================================================================
ggsave("Figure_02.tiff", plot = figure2,
       width = 36, height = 30, units = "cm",
       dpi = 300, compression = "lzw")

ggsave("Figure_02.pdf", plot = figure2,
       width = 36, height = 30, units = "cm")

message("Figure 2 saved: Figure_02.tiff / Figure_02.pdf")
