# =============================================================================
# Figure 10. Phytoplankton-Zooplankton co-occurrence patterns
#
#   Reconstructed to match the criteria stated in Methods 2.5:
#     - Spearman correlations computed on within-sample RELATIVE ABUNDANCE
#       (not raw read counts), to avoid spurious positive correlation driven
#       by sample-to-sample sequencing-depth variation.
#     - Only "dominant" genera are tested: mean relative abundance > 1% AND
#       detected in >= 10% of profiles (within the relevant sample set).
#     - Edges retained require |rho| > 0.6 AND Benjamini-Hochberg
#       FDR-adjusted p < 0.01. No fixed cap on the number of edges is
#       applied (previous version hardcoded n_pos_keep/n_neg_keep per
#       region, which pre-determined the reported positive/negative
#       proportions rather than letting them emerge from the data).
#     - UUV and Net profiles are pooled within each region (Methods 2.5);
#       a method-sensitivity analysis is written to
#       Figure_10_method_sensitivity.csv (section 8).
#     - Jeju Sea and Seogwipo Sea are pooled into a single "Jeju + Seogwipo
#       Sea" region (see Methods 2.5 and Discussion 4.3 for rationale:
#       these are the northern/southern halves of the same island shelf
#       system, and pooling improves statistical power for the network
#       analysis given the small per-region sample sizes, n = 10 and n = 8).
#
#   Figure_10a.png  - Spearman correlation heatmap (dominant phyto x dominant zoo,
#                     nationwide; * = edge meeting |rho| > 0.6 & FDR p < 0.01)
#   Figure_10_edges.csv              - every retained edge (all panels)
#   Figure_10_method_sensitivity.csv - pooled edges re-tested by method
#   Figure_10b.png  - South Sea co-occurrence network
#   Figure_10c.png  - East Sea co-occurrence network
#   Figure_10d.png  - Yellow Sea co-occurrence network
#   Figure_10e.png  - Jeju + Seogwipo Sea co-occurrence network
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

methods   <- as.character(unlist(raw[1, 3:ncol(raw)]))
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
counts_mat[is.na(counts_mat)] <- 0   # empty cells in the spreadsheet = 0 reads

# ── Guild labels (Phytoplankton / Zooplankton) are read directly from the
#    corrected "Taxonomic Group" column of Supplementary Data S1
#    (MetaK90_location.xlsx); no genus-level overrides are applied in code.

phyto_rows <- which(tax_vec == "Phytoplankton")
zoo_rows   <- which(tax_vec == "Zooplankton")

# ── Pool Jeju Sea and Seogwipo Sea into a single region ─────────────────────
# Rationale (Methods 2.5 / Discussion 4.3): Jeju and Seogwipo are the
# northern and southern coastal shelves of the same volcanic island system,
# not independent basins. Pooling substantially improves the statistical
# power of the per-region correlation network (n = 10 and n = 8 individually
# vs. n = 18 pooled) without merging genuinely distinct oceanographic regimes.
locations_grouped <- ifelse(locations %in% c("Jeju Sea", "Seogwipo Sea"),
                            "Jeju + Seogwipo Sea", locations)

# ── Convert to within-sample RELATIVE ABUNDANCE ──────────────────────────────
# Raw read counts are not directly comparable across samples because total
# sequencing depth varies (see Results 3.1). Correlating raw counts across
# samples inflates spurious positive correlations, since most taxa's counts
# rise and fall together with a sample's total depth. Relative abundance
# (each genus's reads divided by that sample's total reads) removes this
# depth confound and is used throughout this script.
sample_totals <- colSums(counts_mat)
rel_mat        <- sweep(counts_mat, 2, sample_totals, FUN = "/")
rownames(rel_mat) <- genus_vec

# =============================================================================
# 2. Colour palette
# =============================================================================
col_phyto <- "#2E8B57"
col_zoo   <- "#E07B00"
col_pos   <- "#C0392B"
col_neg   <- "#2874A6"

# =============================================================================
# 3. Helper: select "dominant" genera within a given sample set
#    (Methods 2.5: mean relative abundance > 1% AND detected in >= 10% of
#    profiles). Returns row indices into rel_mat/counts_mat/genus_vec.
# =============================================================================
select_dominant <- function(rel_sub_mat, candidate_rows,
                            min_mean_rel = 0.01, min_prevalence = 0.10) {
  sub    <- rel_sub_mat[candidate_rows, , drop = FALSE]
  m_rel  <- rowMeans(sub)
  preval <- rowMeans(sub > 0)
  candidate_rows[m_rel > min_mean_rel & preval >= min_prevalence]
}

# =============================================================================
# 4. Helper: Spearman correlations between two sets of genera, with BH-FDR
#    correction and the |rho| > 0.6 & p_adj < 0.01 edge-selection rule.
#    No cap is placed on the number of edges retained per sign.
# =============================================================================
build_edges <- function(rel_sub_mat, phyto_rows_sub, zoo_rows_sub) {
  combos <- expand.grid(p = phyto_rows_sub, z = zoo_rows_sub)
  rs <- numeric(nrow(combos))
  ps <- numeric(nrow(combos))
  for (k in seq_len(nrow(combos))) {
    ct     <- suppressWarnings(cor.test(rel_sub_mat[combos$p[k], ],
                                        rel_sub_mat[combos$z[k], ],
                                        method = "spearman", exact = FALSE))
    rs[k] <- ct$estimate
    ps[k] <- ct$p.value
  }
  padj <- p.adjust(ps, method = "BH")

  data.frame(
    from = genus_vec[combos$p],
    to   = genus_vec[combos$z],
    r    = rs,
    p    = ps,
    padj = padj,
    stringsAsFactors = FALSE
  ) %>%
    filter(abs(r) > 0.6, padj < 0.01)
}

# Helper: label an edge table with its panel (works for zero-row tables)
tag_edges <- function(df, panel, n_samples) {
  if (nrow(df) == 0) return(NULL)
  data.frame(panel = panel, n_samples = n_samples, df, stringsAsFactors = FALSE)
}

# =============================================================================
# 5. Figure 10a - nationwide correlation heatmap, dominant genera only
# =============================================================================
dom_phyto_nat <- select_dominant(rel_mat, phyto_rows)
dom_zoo_nat   <- select_dominant(rel_mat, zoo_rows)

message(sprintf("Nationwide dominant genera: %d phytoplankton, %d zooplankton",
               length(dom_phyto_nat), length(dom_zoo_nat)))

rmat <- matrix(NA, nrow = length(dom_phyto_nat), ncol = length(dom_zoo_nat),
              dimnames = list(genus_vec[dom_phyto_nat], genus_vec[dom_zoo_nat]))

for (i in seq_along(dom_phyto_nat)) {
  for (j in seq_along(dom_zoo_nat)) {
    ct <- suppressWarnings(cor.test(rel_mat[dom_phyto_nat[i], ],
                                    rel_mat[dom_zoo_nat[j], ],
                                    method = "spearman", exact = FALSE))
    rmat[i, j] <- ct$estimate
  }
}

heat_df <- melt(rmat, varnames = c("Phyto", "Zoo"), value.name = "r")
heat_df$Phyto <- as.character(heat_df$Phyto)
heat_df$Zoo   <- as.character(heat_df$Zoo)

# Nationwide edges under the same rule as the regional networks
# (|rho| > 0.6 and BH-FDR-adjusted p < 0.01, FDR applied across all
#  dominant-phyto x dominant-zoo tests nationwide)
nat_edges <- build_edges(rel_mat, dom_phyto_nat, dom_zoo_nat)
message(sprintf("Nationwide edges retained: %d (positive %d, negative %d)",
                nrow(nat_edges), sum(nat_edges$r > 0), sum(nat_edges$r < 0)))
print(nat_edges)

heat_df$sig   <- paste(heat_df$Phyto, heat_df$Zoo) %in%
                 paste(nat_edges$from, nat_edges$to)
heat_df$label <- ifelse(heat_df$sig,
                        paste0(sprintf("%.2f", heat_df$r), "*"),
                        sprintf("%.2f", heat_df$r))
# order axes by descending total relative abundance for readability
phyto_order <- genus_vec[dom_phyto_nat[order(rowMeans(rel_mat[dom_phyto_nat, ]), decreasing = TRUE)]]
zoo_order   <- genus_vec[dom_zoo_nat[order(rowMeans(rel_mat[dom_zoo_nat, ]), decreasing = TRUE)]]
heat_df$Phyto <- factor(heat_df$Phyto, levels = rev(phyto_order))
heat_df$Zoo   <- factor(heat_df$Zoo,   levels = zoo_order)

fig_a <- ggplot(heat_df, aes(x = Zoo, y = Phyto, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = label, fontface = ifelse(sig, "bold", "plain")),
            size = 2.8, colour = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Spearman\nCorrelation\n(rel. abund.)",
    guide    = guide_colorbar(barheight = unit(8, "cm"),
                              barwidth  = unit(0.5, "cm"))
  ) +
  labs(x = NULL, y = NULL,
       caption = paste0("Dominant genera only (mean relative abundance > 1% in >= 10% of profiles, nationwide). ",
                        "* |rho| > 0.6 and BH-FDR-adjusted p < 0.01")) +
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
    plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
    plot.margin  = margin(10, 10, 10, 10)
  )

png("Figure_10a.png", width = 20, height = 14, units = "cm", res = 300, bg = "white")
print(fig_a)
dev.off()

pdf("Figure_10a.pdf", width = 20 / 2.54, height = 14 / 2.54)
print(fig_a)
dev.off()

message("Saved: Figure_10a.png / Figure_10a.pdf")

# =============================================================================
# 6. Function: build co-occurrence network for one region
# =============================================================================
make_network <- function(sea) {

  sea_cols <- which(locations_grouped == sea)
  n_samp   <- length(sea_cols)
  sea_rel  <- rel_mat[, sea_cols, drop = FALSE]

  dom_phyto <- select_dominant(sea_rel, phyto_rows)
  dom_zoo   <- select_dominant(sea_rel, zoo_rows)

  if (length(dom_phyto) == 0 || length(dom_zoo) == 0) {
    message(sprintf("%s: no dominant genera meeting the prevalence filter - skipped", sea))
    return(NULL)
  }

  sig_df <- build_edges(sea_rel, dom_phyto, dom_zoo)

  if (nrow(sig_df) == 0) {
    message(sprintf(
      "%s: %d dominant phyto x %d dominant zoo tested, but no edge met |rho| > 0.6 & FDR p < 0.01 - skipped",
      sea, length(dom_phyto), length(dom_zoo)))
    return(NULL)
  }

  n_pos <- sum(sig_df$r > 0)
  n_neg <- sum(sig_df$r < 0)

  # keep a copy of the edge table for export / sensitivity analysis
  all_edges[[sea]] <<- tag_edges(sig_df, sea, n_samp)

  edge_nodes <- unique(c(sig_df$from, sig_df$to))
  node_df <- data.frame(
    name  = edge_nodes,
    group = ifelse(edge_nodes %in% genus_vec[phyto_rows], "Phytoplankton", "Zooplankton"),
    total = rowMeans(sea_rel[match(edge_nodes, genus_vec), , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  node_df$size <- rescale(log1p(node_df$total), to = c(4, 14))

  n_nodes <- nrow(node_df)
  n_edges <- nrow(sig_df)

  caption_txt <- paste0(
    sea, " (", n_samp, " samples, ", length(dom_phyto), " dominant phyto, ",
    length(dom_zoo), " dominant zoo tested)\n",
    n_nodes, " nodes, ", n_edges, " edges - Positive: ", n_pos, ", Negative: ", n_neg
  )

  g <- tbl_graph(
    nodes    = node_df,
    edges    = sig_df %>% filter(from %in% node_df$name & to %in% node_df$name),
    directed = FALSE
  )

  set.seed(42)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(colour   = ifelse(r > 0, "Positive", "Negative"),
          linetype = ifelse(r > 0, "solid", "dashed")),
      linewidth = 0.8, alpha = 0.85
    ) +
    scale_edge_colour_manual(
      values = c("Positive" = col_pos, "Negative" = col_neg),
      name   = NULL,
      guide  = guide_legend(override.aes = list(linewidth = 1.2))
    ) +
    scale_edge_linetype_identity() +
    geom_node_point(
      aes(size = size, fill = group),
      shape = 21, colour = "white", stroke = 0.5
    ) +
    scale_fill_manual(
      values = c("Phytoplankton" = col_phyto, "Zooplankton" = col_zoo),
      name   = NULL
    ) +
    scale_size_identity() +
    geom_node_text(
      aes(label = name),
      size           = 3.0,
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
    # Expand the panel so force-directed nodes never sit flush against the
    # edge, leaving clean room for the caption drawn OUTSIDE the panel below.
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = NULL, y = NULL, caption = caption_txt) +
    theme_void(base_size = 10, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.text     = element_text(size = 8.5, family = "sans"),
      legend.key.size = unit(0.45, "cm"),
      plot.caption    = element_text(size = 9, face = "bold", colour = "grey20",
                                     hjust = 0, margin = margin(t = 12)),
      plot.margin     = margin(10, 10, 10, 10)
    ) +
    guides(
      fill        = guide_legend(override.aes = list(size = 5), title = NULL, order = 2),
      edge_colour = guide_legend(title = NULL, order = 1)
    )

  return(p)
}

# =============================================================================
# 7. Generate and save network figures (b-e) - 4 regions
#    (Jeju Sea and Seogwipo Sea pooled; see section 1 above)
# =============================================================================
sea_panel <- list(
  list(sea = "South Sea",           label = "(b)", file = "Figure_10b"),
  list(sea = "East Sea",            label = "(c)", file = "Figure_10c"),
  list(sea = "Yellow Sea",          label = "(d)", file = "Figure_10d"),
  list(sea = "Jeju + Seogwipo Sea", label = "(e)", file = "Figure_10e")
)

all_edges <- list("Nationwide" = tag_edges(nat_edges, "Nationwide", ncol(rel_mat)))

for (sp in sea_panel) {
  message("Building: ", sp$sea, " ...")

  fig <- make_network(sea = sp$sea)

  if (!is.null(fig)) {
    png(paste0(sp$file, ".png"), width = 20, height = 19, units = "cm",
        res = 300, bg = "white")
    print(fig)
    dev.off()

    pdf(paste0(sp$file, ".pdf"), width = 20 / 2.54, height = 19 / 2.54)
    print(fig)
    dev.off()

    message("Saved: ", sp$file, ".png / ", sp$file, ".pdf")
  } else {
    message("No network plotted for ", sp$sea,
           " (no edges met the |rho| > 0.6 & FDR-adjusted p < 0.01 criteria).")
  }
}

# =============================================================================
# 8. Export edge tables and method-sensitivity analysis
#
#    Because UUV and Net profiles are pooled within each region, a systematic
#    method difference (e.g., a genus enriched in net hauls vs. a genus
#    enriched in UUV profiles) could by itself generate a correlation. Each
#    pooled edge is therefore re-evaluated (i) within UUV profiles only,
#    (ii) within Net profiles only, and (iii) with a method-stratified
#    Spearman correlation (ranks computed separately within each method,
#    which removes between-method shifts in relative abundance).
#    Edges whose sign is not reproduced within methods should be interpreted
#    with caution.
# =============================================================================
edges_out <- do.call(rbind, all_edges)
if (is.null(edges_out) || nrow(edges_out) == 0) stop("No edges retained in any panel.")
rownames(edges_out) <- NULL
edges_out$sign <- ifelse(edges_out$r > 0, "positive", "negative")
write.csv(edges_out, "Figure_10_edges.csv", row.names = FALSE)
message("Saved: Figure_10_edges.csv")
print(table(edges_out$panel, edges_out$sign))

panel_cols <- function(panel) {
  if (panel == "Nationwide") seq_len(ncol(rel_mat))
  else which(locations_grouped == panel)
}

spearman_safe <- function(x, y) {
  if (length(x) < 4 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

stratified_rho <- function(x, y, grp) {
  rx <- ave(x, grp, FUN = function(v) rank(v) / length(v))
  ry <- ave(y, grp, FUN = function(v) rank(v) / length(v))
  if (sd(rx) == 0 || sd(ry) == 0) return(NA_real_)
  cor(rx, ry)
}

sens <- do.call(rbind, lapply(seq_len(nrow(edges_out)), function(k) {
  e    <- edges_out[k, ]
  cols <- panel_cols(e$panel)
  xi   <- match(e$from, genus_vec)
  yi   <- match(e$to,   genus_vec)
  cu   <- cols[methods[cols] == "Drone"]
  cn   <- cols[methods[cols] == "Net"]
  data.frame(
    panel        = e$panel,
    from         = e$from,
    to           = e$to,
    rho_pooled   = round(e$r, 3),
    rho_UUV_only = round(spearman_safe(rel_mat[xi, cu], rel_mat[yi, cu]), 3),
    n_UUV        = length(cu),
    rho_Net_only = round(spearman_safe(rel_mat[xi, cn], rel_mat[yi, cn]), 3),
    n_Net        = length(cn),
    rho_method_stratified = round(stratified_rho(rel_mat[xi, cols],
                                                 rel_mat[yi, cols],
                                                 methods[cols]), 3)
  )
}))
sens$sign_consistent_within_methods <-
  with(sens, sign(rho_UUV_only) == sign(rho_pooled) &
             sign(rho_Net_only) == sign(rho_pooled))
write.csv(sens, "Figure_10_method_sensitivity.csv", row.names = FALSE)
message("Saved: Figure_10_method_sensitivity.csv")
message(sprintf("Edges with the same sign in both UUV-only and Net-only subsets: %d of %d",
                sum(sens$sign_consistent_within_methods, na.rm = TRUE), nrow(sens)))
