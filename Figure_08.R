# =============================================================================
# Figure 8. Top 10 phytoplankton and zooplankton genera by total read counts
#           across five Korean coastal sea regions
#
# Panel layout (5 output files, PNG + PDF):
#   Figure_08ab — (a) South Sea phytoplankton + (b) South Sea zooplankton
#   Figure_08cd — (c) East Sea phytoplankton  + (d) East Sea zooplankton
#   Figure_08ef — (e) Yellow Sea phytoplankton + (f) Yellow Sea zooplankton
#   Figure_08gh — (g) Jeju Sea phytoplankton  + (h) Jeju Sea zooplankton
#   Figure_08ij — (i) Seogwipo Sea phytoplankton + (j) Seogwipo Sea zooplankton
# =============================================================================
# Required packages:
#   install.packages(c("readxl", "dplyr", "ggplot2", "scales", "cowplot"))

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

# =============================================================================
# 1. Load and prepare data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx",
                  col_names = FALSE, .name_repair = "minimal")

locations <- as.character(unlist(raw[2, 3:ncol(raw)]))
genus_vec <- as.character(unlist(raw[6:nrow(raw), 1]))
tax_vec   <- as.character(unlist(raw[6:nrow(raw), 2]))

counts_raw <- raw[6:nrow(raw), 3:ncol(raw)]
counts_mat <- matrix(
  as.numeric(unlist(counts_raw)),
  nrow = nrow(counts_raw),
  ncol = ncol(counts_raw)
)
rownames(counts_mat) <- genus_vec

# ── Reclassify Noctiluca as Zooplankton ──────────────────────────────────────
tax_vec[genus_vec == "Noctiluca"] <- "동물플랑크톤"

phyto_rows <- which(tax_vec == "식물플랑크톤")
zoo_rows   <- which(tax_vec == "동물플랑크톤")

# =============================================================================
# 2. Regional colour scheme (consistent across all figures)
# =============================================================================
region_colours <- c(
  "South Sea"    = "#E05252",
  "East Sea"     = "#4472C4",
  "Yellow Sea"   = "#ED7D31",
  "Jeju Sea"     = "#70AD47",
  "Seogwipo Sea" = "#9B59B6"
)

# =============================================================================
# 3. Function: build top-10 horizontal bar chart
# =============================================================================
make_bar <- function(region, guild) {

  # Select row indices for the taxonomic guild
  row_idx <- if (guild == "Phytoplankton") phyto_rows else zoo_rows

  # Select column indices for the sea region
  col_idx <- which(locations == region)

  # Compute total reads per genus within the region
  totals      <- rowSums(counts_mat[row_idx, col_idx, drop = FALSE])
  genus_names <- genus_vec[row_idx]

  # Keep top 10 genera
  n_top   <- min(10, sum(totals > 0))
  top_ord <- order(totals, decreasing = TRUE)[1:n_top]

  top_df <- data.frame(
    Genus = genus_names[top_ord],
    Reads = totals[top_ord],
    stringsAsFactors = FALSE
  ) %>%
    arrange(Reads) %>%
    mutate(Genus = factor(Genus, levels = Genus))

  fill_col <- region_colours[[region]]

  ggplot(top_df, aes(x = Genus, y = Reads)) +
    geom_col(fill      = fill_col,
             width     = 0.72,
             colour    = "white",
             linewidth = 0.15) +
    scale_y_continuous(labels = comma,
                       expand = expansion(mult = c(0, 0.05))) +
    coord_flip() +
    labs(
      title = paste0("Top 10 ", guild, "\n", region),
      x     = NULL,
      y     = "Total Reads"
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title   = element_text(size       = 10.5,
                                  face       = "bold",
                                  hjust      = 0.5,
                                  lineheight = 1.2,
                                  margin     = margin(b = 4)),
      axis.text.y  = element_text(size   = 9,
                                  face   = "bold.italic",
                                  hjust  = 1,
                                  colour = "black"),
      axis.text.x  = element_text(size = 8.5),
      axis.title.x = element_text(size = 9.5, face = "bold"),
      axis.line    = element_line(linewidth = 0.4),
      plot.margin  = margin(8, 12, 8, 8)
    )
}

# =============================================================================
# 4. Function: assemble and save one paired panel
# =============================================================================
save_pair <- function(region, label_left, label_right, out_name) {

  fig_phyto <- make_bar(region, "Phytoplankton")
  fig_zoo   <- make_bar(region, "Zooplankton")

  combined <- plot_grid(
    fig_phyto, fig_zoo,
    ncol           = 2,
    labels         = c(label_left, label_right),
    label_size     = 13,
    label_fontface = "bold"
  )

  # ── PNG output ──────────────────────────────────────────────────────────────
  png(filename = paste0(out_name, ".png"),
      width    = 28,
      height   = 14,
      units    = "cm",
      res      = 300,
      bg       = "white")
  print(combined)
  dev.off()

  # ── PDF output ──────────────────────────────────────────────────────────────
  pdf(file   = paste0(out_name, ".pdf"),
      width  = 28 / 2.54,   # convert cm to inches
      height = 14 / 2.54)
  print(combined)
  dev.off()

  message("Saved: ", out_name, ".png  |  ", out_name, ".pdf")
}

# =============================================================================
# 5. Generate and save all five paired figures
# =============================================================================

# (a) + (b): South Sea
save_pair("South Sea",
          label_left  = "(a)",
          label_right = "(b)",
          out_name    = "Figure_08ab")

# (c) + (d): East Sea
save_pair("East Sea",
          label_left  = "(c)",
          label_right = "(d)",
          out_name    = "Figure_08cd")

# (e) + (f): Yellow Sea
save_pair("Yellow Sea",
          label_left  = "(e)",
          label_right = "(f)",
          out_name    = "Figure_08ef")

# (g) + (h): Jeju Sea
save_pair("Jeju Sea",
          label_left  = "(g)",
          label_right = "(h)",
          out_name    = "Figure_08gh")

# (i) + (j): Seogwipo Sea
save_pair("Seogwipo Sea",
          label_left  = "(i)",
          label_right = "(j)",
          out_name    = "Figure_08ij")

message("\nAll Figure 8 panels saved successfully.")
message("Output location: ", getwd())
