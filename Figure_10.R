# =============================================================================
# Figure 10. Phytoplankton–Zooplankton co-occurrence patterns
#
#   Figure_10a.png  — Spearman correlation heatmap (top15 phyto × top15 zoo)
#   Figure_10b.png  — South Sea co-occurrence network
#   Figure_10c.png  — East Sea co-occurrence network
#   Figure_10d.png  — Yellow Sea co-occurrence network
#   Figure_10e.png  — Jeju Sea co-occurrence network
#   Figure_10f.png  — Seogwipo Sea co-occurrence network
# =============================================================================
# Required packages:
#   install.packages(c("readxl","dplyr","ggplot2","scales",
#                      "igraph","ggraph","tidygraph","reshape2"))

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(igraph)
library(ggraph)
library(tidygraph)
library(reshape2)

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

# ── Global top 15 genera for emphasis labelling ───────────────────────────────
phyto_totals  <- rowSums(counts_mat[phyto_rows, ])
zoo_totals    <- rowSums(counts_mat[zoo_rows,   ])
top15p_global <- genus_vec[phyto_rows[order(phyto_totals, decreasing = TRUE)[1:15]]]
top15z_global <- genus_vec[zoo_rows  [order(zoo_totals,   decreasing = TRUE)[1:15]]]
label_top15   <- c(top15p_global, top15z_global)

# =============================================================================
# 2. Colour palette
# =============================================================================
col_phyto <- "#33A145"
col_zoo   <- "#E07B00"
col_pos   <- "#E05252"
col_neg   <- "#5B8DD9"

# =============================================================================
# 3. Figure 10a — Spearman correlation heatmap (top15 phyto × top15 zoo)
# =============================================================================
top15p_rows <- phyto_rows[order(phyto_totals, decreasing = TRUE)[1:15]]
top15z_rows <- zoo_rows  [order(zoo_totals,   decreasing = TRUE)[1:15]]

# Compute pairwise Spearman correlations across ALL samples
rmat <- matrix(NA, nrow = 15, ncol = 15,
               dimnames = list(genus_vec[top15p_rows],
                               genus_vec[top15z_rows]))

for (i in 1:15) {
  for (j in 1:15) {
    ct        <- cor.test(counts_mat[top15p_rows[i], ],
                          counts_mat[top15z_rows[j], ],
                          method = "spearman", exact = FALSE)
    rmat[i, j] <- ct$estimate
  }
}

# Reshape to long format
heat_df <- melt(rmat, varnames = c("Phyto", "Zoo"), value.name = "r")
heat_df$label <- as.character(round(heat_df$r, 2))
heat_df$Phyto <- factor(heat_df$Phyto,
                        levels = rev(genus_vec[top15p_rows]))
heat_df$Zoo   <- factor(heat_df$Zoo,
                        levels = genus_vec[top15z_rows])

fig_a <- ggplot(heat_df, aes(x = Zoo, y = Phyto, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 2.5, colour = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-0.8, 0.8),
    name     = "Spearman\nCorrelation",
    guide    = guide_colorbar(barheight = unit(8, "cm"),
                              barwidth  = unit(0.5, "cm"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold.italic",
                                colour = col_zoo,   size = 8.5),
    axis.text.y  = element_text(face  = "bold.italic",
                                colour = col_phyto, size = 8.5),
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text  = element_text(size = 7.5),
    plot.margin  = margin(10, 10, 10, 10)
  )

# Save Figure 10a
png("Figure_10a.png",
    width = 22, height = 16, units = "cm",
    res = 300, bg = "white")
print(fig_a)
dev.off()

pdf("Figure_10a.pdf",
    width = 22 / 2.54, height = 16 / 2.54)
print(fig_a)
dev.off()

message("Saved: Figure_10a.png / Figure_10a.pdf")

# =============================================================================
# 4. Function: build co-occurrence network (Option B — edge nodes only)
# =============================================================================
make_network <- function(sea, n_pos_keep, n_neg_keep) {

  # ── Subset data for this sea region ────────────────────────────────────────
  sea_cols  <- which(locations == sea)
  n_samp    <- length(sea_cols)
  sea_mat   <- counts_mat[, sea_cols]

  det_rows  <- which(rowSums(sea_mat) > 0)
  det_phyto <- intersect(det_rows, phyto_rows)
  det_zoo   <- intersect(det_rows, zoo_rows)

  if (length(det_phyto) == 0 || length(det_zoo) == 0) return(NULL)

  # ── Compute all phyto × zoo Spearman correlations ─────────────────────────
  n_pairs  <- length(det_phyto) * length(det_zoo)
  rs       <- numeric(n_pairs)
  ps       <- numeric(n_pairs)
  from_vec <- character(n_pairs)
  to_vec   <- character(n_pairs)

  k <- 1
  for (pi in det_phyto) {
    for (zi in det_zoo) {
      ct          <- cor.test(sea_mat[pi, ], sea_mat[zi, ],
                               method = "spearman", exact = FALSE)
      rs[k]       <- ct$estimate
      ps[k]       <- ct$p.value
      from_vec[k] <- genus_vec[pi]
      to_vec[k]   <- genus_vec[zi]
      k <- k + 1
    }
  }

  # ── Select edges: top N positive + top N negative (p < 0.05) ──────────────
  all_df <- data.frame(from = from_vec, to = to_vec,
                       r = rs, p = ps,
                       stringsAsFactors = FALSE) %>%
    filter(p < 0.05)

  pos_df <- all_df %>% filter(r > 0) %>% arrange(desc(r)) %>% head(n_pos_keep)
  neg_df <- all_df %>% filter(r < 0) %>% arrange(r)       %>% head(n_neg_keep)
  sig_df <- bind_rows(pos_df, neg_df)

  if (nrow(sig_df) == 0) return(NULL)

  n_pos  <- nrow(pos_df)
  n_neg  <- nrow(neg_df)

  # ── Option B: only nodes involved in edges ─────────────────────────────────
  edge_nodes <- unique(c(sig_df$from, sig_df$to))

  node_df <- data.frame(
    name       = edge_nodes,
    group      = ifelse(edge_nodes %in% genus_vec[phyto_rows],
                        "Phytoplankton", "Zooplankton"),
    total      = rowSums(sea_mat[match(edge_nodes, genus_vec), , drop = FALSE]),
    show_top15 = edge_nodes %in% label_top15,
    stringsAsFactors = FALSE
  )

  # Log-scaled node size
  node_df$size <- log1p(node_df$total)
  node_df$size <- rescale(node_df$size, to = c(2, 14))

  n_nodes <- nrow(node_df)
  n_edges <- nrow(sig_df)

  subtitle <- paste0(
    sea, "\n(",
    n_samp, " samples, ", n_nodes, " nodes, ", n_edges, " edges)\n",
    "Positive: ", n_pos, ", Negative: ", n_neg
  )

  # ── Build tidygraph object ─────────────────────────────────────────────────
  g <- tbl_graph(
    nodes    = node_df,
    edges    = sig_df %>%
      filter(from %in% node_df$name & to %in% node_df$name),
    directed = FALSE
  )

  # ── Plot ──────────────────────────────────────────────────────────────────
  set.seed(42)
  p <- ggraph(g, layout = "fr") +

    # Edges
    geom_edge_link(
      aes(colour   = ifelse(r > 0, "Positive", "Negative"),
          linetype = ifelse(r > 0, "solid", "dashed")),
      linewidth = 0.7, alpha = 0.80
    ) +
    scale_edge_colour_manual(
      values = c("Positive" = col_pos, "Negative" = col_neg),
      name   = NULL,
      guide  = guide_legend(override.aes = list(linewidth = 1.2))
    ) +
    scale_edge_linetype_identity() +

    # Nodes
    geom_node_point(
      aes(size = size, fill = group),
      shape  = 21, colour = "white", stroke = 0.5
    ) +
    scale_fill_manual(
      values = c("Phytoplankton" = col_phyto, "Zooplankton" = col_zoo),
      name   = NULL
    ) +
    scale_size_identity() +

    # All nodes: bold italic black with white outline
    # point.padding: distance between label and node circle edge
    # box.padding: minimum distance between label boxes
    geom_node_text(
      data           = function(x) x,
      aes(label      = name),
      size           = 2.8,
      fontface       = "bold.italic",
      colour         = "black",
      repel          = TRUE,
      max.overlaps   = Inf,
      point.padding  = unit(1.2, "lines"),
      box.padding    = unit(0.6, "lines"),
      segment.colour = "grey50",
      segment.size   = 0.3,
      segment.alpha  = 0.6,
      bg.colour      = "white",
      bg.r           = 0.15
    ) +

    # Subtitle
    annotate("text", x = -Inf, y = -Inf,
             label  = subtitle,
             hjust  = 0, vjust = -0.3,
             size   = 3.2, fontface = "bold",
             colour = "grey20") +

    labs(x = NULL, y = NULL) +
    theme_void(base_size = 10, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.text     = element_text(size = 8.5, family = "sans"),
      legend.key.size = unit(0.45, "cm"),
      plot.margin     = margin(10, 10, 25, 10)
    ) +
    guides(
      fill        = guide_legend(override.aes = list(size = 5),
                                 title = NULL, order = 2),
      edge_colour = guide_legend(title = NULL, order = 1)
    )

  return(p)
}

# =============================================================================
# 5. Generate and save network figures (b–f)
# =============================================================================
sea_panel <- list(
  list(sea = "South Sea",    label = "(b)", file = "Figure_10b",
       n_pos = 28, n_neg = 12),
  list(sea = "East Sea",     label = "(c)", file = "Figure_10c",
       n_pos = 18, n_neg = 22),
  list(sea = "Yellow Sea",   label = "(d)", file = "Figure_10d",
       n_pos = 29, n_neg = 11),
  list(sea = "Jeju Sea",     label = "(e)", file = "Figure_10e",
       n_pos = 22, n_neg = 18),
  list(sea = "Seogwipo Sea", label = "(f)", file = "Figure_10f",
       n_pos = 12, n_neg = 28)
)

for (sp in sea_panel) {
  message("Building: ", sp$sea, " ...")

  fig <- make_network(sea      = sp$sea,
                      n_pos_keep = sp$n_pos,
                      n_neg_keep = sp$n_neg)

  if (!is.null(fig)) {
    png(paste0(sp$file, ".png"),
        width = 20, height = 18, units = "cm",
        res = 300, bg = "white")
    print(fig)
    dev.off()

    pdf(paste0(sp$file, ".pdf"),
        width = 20 / 2.54, height = 18 / 2.54)
    print(fig)
    dev.off()

    message("  Saved: ", sp$file, ".png / ", sp$file, ".pdf")
  }
}

message("\nAll Figure 10 panels (a–f) saved successfully.")
message("Output location: ", getwd())
