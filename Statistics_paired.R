# =============================================================================
# Statistics_paired.R
#   Paired UUV vs Net statistics reported in Results 3.1-3.2 (Methods 2.5)
#
#   1. Station-level pairing of UUV and Net profiles (62 pairs expected)
#   2. Per-profile metrics: reads, genus richness (total / phyto / zoo),
#      Shannon diversity (total / phyto / zoo), phytoplankton read proportion
#   3. Descriptive statistics by method (mean, SD, median, IQR)
#   4. Paired Wilcoxon signed-rank tests on the matched stations
#   5. Library size vs. richness (Spearman)
#   6. Rarefaction to the minimum library size (Hurlbert 1971; vegan::rarefy)
#      + rarefied Shannon diversity (mean of 100 random subsamplings)
#   7. Regional summaries (Figure 4 values; Results 3.2 phyto read proportion)
#   8. Method-exclusive genera by guild (Figure 5f) and pooled read fractions
#   9. Beta diversity: PERMANOVA and NMDS (Bray-Curtis on relative abundance)
#
#   All outputs are written to ./stats_output/
#   Shannon diversity is computed at the genus level (natural log),
#   matching Figure 5.
# =============================================================================
# install.packages(c("readxl", "dplyr", "tidyr", "vegan", "ggplot2", "permute"))

library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(permute)
library(ggplot2)

dir.create("stats_output", showWarnings = FALSE)
set.seed(2024)

# =============================================================================
# 0. Load data
# =============================================================================
raw <- read_excel("MetaK90_location.xlsx", col_names = FALSE,
                  .name_repair = "minimal")

methods   <- as.character(unlist(raw[1, 3:ncol(raw)]))
locations <- as.character(unlist(raw[2, 3:ncol(raw)]))
lats      <- as.numeric(unlist(raw[3, 3:ncol(raw)]))
lons      <- as.numeric(unlist(raw[4, 3:ncol(raw)]))
genus_vec <- as.character(unlist(raw[6:nrow(raw), 1]))
tax_vec   <- as.character(unlist(raw[6:nrow(raw), 2]))   # corrected guild labels

reads_matrix <- apply(raw[6:nrow(raw), 3:ncol(raw)], 2,
                      function(x) as.numeric(as.character(x)))
reads_matrix[is.na(reads_matrix)] <- 0
rownames(reads_matrix) <- genus_vec

phyto_rows <- which(tax_vec == "Phytoplankton")
zoo_rows   <- which(tax_vec == "Zooplankton")
method_lab <- ifelse(methods == "Drone", "UUV", "Net")

# =============================================================================
# 1. Station-level pairing
#    Each UUV profile is paired with the nearest Net profile in the same
#    region within PAIR_MAX_KM (greedy matching by increasing distance; each
#    profile used once). Same rule as Figure_01.R.
#    >> Check that 62 pairs are found and that the unpaired profile is the
#       South Sea UUV profile whose net sample failed QC.
# =============================================================================
PAIR_MAX_KM <- 0.2

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad  <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

uuv_i <- which(methods == "Drone")
net_i <- which(methods == "Net")
cand  <- expand.grid(u = uuv_i, n = net_i)
cand  <- cand[locations[cand$u] == locations[cand$n], ]
cand$d <- haversine_km(lats[cand$u], lons[cand$u], lats[cand$n], lons[cand$n])
cand  <- cand[cand$d <= PAIR_MAX_KM, ]
cand  <- cand[order(cand$d), ]
pairs <- cand[0, ]
for (k in seq_len(nrow(cand))) {
  if (!(cand$u[k] %in% pairs$u) && !(cand$n[k] %in% pairs$n))
    pairs <- rbind(pairs, cand[k, ])
}
pairs$pair_id  <- paste0("P", sprintf("%02d", seq_len(nrow(pairs))))
pairs$Location <- locations[pairs$u]

cat(sprintf("Paired stations: %d | max UUV-Net distance: %.0f m\n",
            nrow(pairs), 1000 * max(c(0, pairs$d))))
cat("Unpaired UUV profiles (column index in sheet, region): ",
    paste0(setdiff(uuv_i, pairs$u) + 2, " (", locations[setdiff(uuv_i, pairs$u)], ")",
           collapse = "; "), "\n")
cat("Unpaired Net profiles: ",
    paste0(setdiff(net_i, pairs$n) + 2, " (", locations[setdiff(net_i, pairs$n)], ")",
           collapse = "; "), "\n")
print(table(pairs$Location))
write.csv(data.frame(pair_id = pairs$pair_id, Location = pairs$Location,
                     UUV_sheet_column = pairs$u + 2, Net_sheet_column = pairs$n + 2,
                     distance_m = round(1000 * pairs$d)),
          "stats_output/station_pairs.csv", row.names = FALSE)

# =============================================================================
# 2. Per-profile metrics
# =============================================================================
shannon_H <- function(v) {
  v <- v[v > 0]
  if (length(v) == 0) return(0)
  p <- v / sum(v)
  -sum(p * log(p))
}

sample_df <- data.frame(
  col         = seq_len(ncol(reads_matrix)),
  Method      = factor(method_lab, levels = c("UUV", "Net")),
  Location    = locations,
  Reads       = colSums(reads_matrix),
  Richness    = colSums(reads_matrix > 0),
  PhytoRich   = colSums(reads_matrix[phyto_rows, ] > 0),
  ZooRich     = colSums(reads_matrix[zoo_rows, ] > 0),
  Shannon     = apply(reads_matrix, 2, shannon_H),
  ShannonPhyto = apply(reads_matrix[phyto_rows, ], 2, shannon_H),
  ShannonZoo  = apply(reads_matrix[zoo_rows, ], 2, shannon_H),
  PhytoProp   = colSums(reads_matrix[phyto_rows, ]) / colSums(reads_matrix)
)
sample_df$pair_id <- NA_character_
sample_df$pair_id[pairs$u] <- pairs$pair_id
sample_df$pair_id[pairs$n] <- pairs$pair_id

cat(sprintf("\nMinimum library size (rarefaction depth): %d reads\n",
            min(sample_df$Reads)))
cat(sprintf("Genus-assigned reads: total %s, mean %.0f per profile, range %d-%d\n",
            format(sum(sample_df$Reads), big.mark = ","), mean(sample_df$Reads),
            min(sample_df$Reads), max(sample_df$Reads)))

# =============================================================================
# 6a. Rarefaction (computed here so it enters the same tables)
#     Expected richness at the minimum depth: Hurlbert (1971) analytical
#     estimator (vegan::rarefy). Rarefied Shannon: mean over 100 random
#     subsamplings without replacement (vegan::rrarefy).
# =============================================================================
depth     <- min(sample_df$Reads)
count_int <- t(round(reads_matrix))            # samples x genera, integers
sample_df$RichnessRarefied <- as.numeric(rarefy(count_int, sample = depth))

N_ITER <- 100
shan_rar <- replicate(N_ITER, diversity(rrarefy(count_int, sample = depth),
                                        index = "shannon"))
sample_df$ShannonRarefied <- rowMeans(shan_rar)

write.csv(sample_df, "stats_output/per_profile_metrics.csv", row.names = FALSE)

# =============================================================================
# 3. Descriptive statistics by method (all profiles: UUV n = 63, Net n = 62)
# =============================================================================
metrics <- c("Reads", "Richness", "PhytoRich", "ZooRich", "Shannon",
             "ShannonPhyto", "ShannonZoo", "PhytoProp",
             "RichnessRarefied", "ShannonRarefied")

desc <- sample_df %>%
  select(Method, all_of(metrics)) %>%
  pivot_longer(-Method, names_to = "Metric", values_to = "Value") %>%
  group_by(Metric, Method) %>%
  summarise(n = n(), Mean = mean(Value), SD = sd(Value),
            Median = median(Value),
            Q1 = quantile(Value, 0.25), Q3 = quantile(Value, 0.75),
            Min = min(Value), Max = max(Value), .groups = "drop") %>%
  mutate(Metric = factor(Metric, levels = metrics)) %>%
  arrange(Metric, Method)

cat("\n===== Descriptive statistics (all profiles) =====\n")
print(as.data.frame(desc), digits = 4)
write.csv(desc, "stats_output/descriptive_by_method.csv", row.names = FALSE)

# =============================================================================
# 4. Paired Wilcoxon signed-rank tests (62 matched stations)
# =============================================================================
paired_wide <- sample_df %>%
  filter(!is.na(pair_id)) %>%
  select(pair_id, Location, Method, all_of(metrics)) %>%
  pivot_wider(names_from = Method, values_from = all_of(metrics))

paired_test <- function(metric) {
  x <- paired_wide[[paste0(metric, "_UUV")]]
  y <- paired_wide[[paste0(metric, "_Net")]]
  wt <- wilcox.test(x, y, paired = TRUE, exact = FALSE)
  data.frame(Metric = metric, n_pairs = length(x),
             Mean_UUV = mean(x), Mean_Net = mean(y),
             Median_UUV = median(x), Median_Net = median(y),
             Median_diff = median(x - y),
             V = unname(wt$statistic), p = wt$p.value)
}
paired_res <- do.call(rbind, lapply(metrics, paired_test))
paired_res$p_formatted <- ifelse(paired_res$p < 0.0001, "< 0.0001",
                                 sprintf("%.4f", paired_res$p))

cat("\n===== Paired Wilcoxon signed-rank tests (station-matched) =====\n")
print(paired_res, digits = 4)
write.csv(paired_res, "stats_output/paired_wilcoxon.csv", row.names = FALSE)

# =============================================================================
# 5. Library size vs. genus richness (Spearman)
# =============================================================================
lib_cor <- function(d, label) {
  ct <- cor.test(d$Reads, d$Richness, method = "spearman", exact = FALSE)
  data.frame(Subset = label, n = nrow(d), rho = unname(ct$estimate),
             p = ct$p.value)
}
lib_res <- rbind(lib_cor(sample_df, "All"),
                 lib_cor(filter(sample_df, Method == "UUV"), "UUV"),
                 lib_cor(filter(sample_df, Method == "Net"), "Net"))
cat("\n===== Library size vs richness (Spearman) =====\n")
print(lib_res, digits = 3)
write.csv(lib_res, "stats_output/library_size_vs_richness.csv", row.names = FALSE)

# =============================================================================
# 7. Regional summaries
# =============================================================================
loc_levels <- c("South Sea", "East Sea", "Yellow Sea", "Jeju Sea", "Seogwipo Sea")

# (a) Figure 4 values: mean total richness and median guild richness
regional_rich <- sample_df %>%
  mutate(Location = factor(Location, levels = loc_levels)) %>%
  group_by(Location, Method) %>%
  summarise(n = n(),
            MeanRichness   = mean(Richness),
            MedianPhytoRich = median(PhytoRich),
            MedianZooRich  = median(ZooRich),
            .groups = "drop")
cat("\n===== Regional richness by method (Figure 4) =====\n")
print(as.data.frame(regional_rich), digits = 4)
write.csv(regional_rich, "stats_output/regional_richness.csv", row.names = FALSE)

# (b) Results 3.2, first paragraph: per-profile phytoplankton read proportion
#     by region (UUV and Net pooled, and by method)
prop_summary <- function(d) {
  d %>% summarise(n = n(),
                  Median = 100 * median(PhytoProp),
                  Q1 = 100 * quantile(PhytoProp, 0.25),
                  Q3 = 100 * quantile(PhytoProp, 0.75),
                  Min = 100 * min(PhytoProp),
                  Max = 100 * max(PhytoProp), .groups = "drop")
}
regional_prop_pooled <- sample_df %>%
  mutate(Location = factor(Location, levels = loc_levels)) %>%
  group_by(Location) %>% prop_summary()
regional_prop_method <- sample_df %>%
  mutate(Location = factor(Location, levels = loc_levels)) %>%
  group_by(Location, Method) %>% prop_summary()
cat("\n===== Phytoplankton read proportion (%) by region, methods pooled =====\n")
print(as.data.frame(regional_prop_pooled), digits = 3)
cat("\n===== Phytoplankton read proportion (%) by region and method =====\n")
print(as.data.frame(regional_prop_method), digits = 3)
write.csv(regional_prop_pooled, "stats_output/regional_phyto_proportion_pooled.csv",
          row.names = FALSE)
write.csv(regional_prop_method, "stats_output/regional_phyto_proportion_by_method.csv",
          row.names = FALSE)

# =============================================================================
# 8. Method-exclusive genera (Figure 5f) and pooled read fractions (Figure 5c)
# =============================================================================
uuv_cols <- which(method_lab == "UUV")
net_cols <- which(method_lab == "Net")
uuv_det  <- genus_vec[rowSums(reads_matrix[, uuv_cols] > 0) > 0]
net_det  <- genus_vec[rowSums(reads_matrix[, net_cols] > 0) > 0]
uuv_only <- setdiff(uuv_det, net_det)
net_only <- setdiff(net_det, uuv_det)
tax_map  <- setNames(tax_vec, genus_vec)

excl <- data.frame(
  Category = c("Shared", "UUV only", "Net only",
               "Phyto UUV only", "Phyto Net only", "Zoo UUV only", "Zoo Net only"),
  Count = c(length(intersect(uuv_det, net_det)), length(uuv_only), length(net_only),
            sum(tax_map[uuv_only] == "Phytoplankton"),
            sum(tax_map[net_only] == "Phytoplankton"),
            sum(tax_map[uuv_only] == "Zooplankton"),
            sum(tax_map[net_only] == "Zooplankton"))
)
pooled_frac <- data.frame(
  Method = c("UUV", "Net"),
  PhytoReadPct = 100 * c(sum(reads_matrix[phyto_rows, uuv_cols]) / sum(reads_matrix[, uuv_cols]),
                         sum(reads_matrix[phyto_rows, net_cols]) / sum(reads_matrix[, net_cols]))
)
pooled_frac$ZooReadPct <- 100 - pooled_frac$PhytoReadPct

cat("\n===== Genus inventory (Figure 2a / 5f) =====\n");  print(excl)
cat("\n===== Pooled read fractions (Figure 5c) =====\n"); print(pooled_frac, digits = 3)
write.csv(excl, "stats_output/method_exclusive_genera.csv", row.names = FALSE)
write.csv(pooled_frac, "stats_output/pooled_read_fractions.csv", row.names = FALSE)

# =============================================================================
# 9. Beta diversity: PERMANOVA and NMDS (Bray-Curtis, relative abundance)
#    (a) Region x Method on all 125 profiles (free permutations)
#    (b) Method effect within stations: 62 matched pairs, permutations
#        restricted within station (blocks), so the test respects pairing
# =============================================================================
rel_t <- t(sweep(reads_matrix, 2, colSums(reads_matrix), "/"))   # samples x genera
bc    <- vegdist(rel_t, method = "bray")

meta <- data.frame(Location = factor(locations, levels = loc_levels),
                   Method   = factor(method_lab, levels = c("UUV", "Net")))

perm_a <- adonis2(bc ~ Location * Method, data = meta,
                  permutations = 9999, by = "terms")
cat("\n===== PERMANOVA (a): Region x Method, all profiles =====\n")
print(perm_a)

keep    <- !is.na(sample_df$pair_id)
bc_pair <- vegdist(rel_t[keep, ], method = "bray")
meta_p  <- data.frame(Method = meta$Method[keep],
                      pair   = factor(sample_df$pair_id[keep]))
ctrl    <- how(blocks = meta_p$pair, nperm = 9999)
perm_b  <- adonis2(bc_pair ~ Method, data = meta_p, permutations = ctrl)
cat("\n===== PERMANOVA (b): Method within stations (62 pairs, blocked) =====\n")
print(perm_b)

disp <- permutest(betadisper(bc, meta$Method), permutations = 9999)
cat("\n===== Homogeneity of multivariate dispersion (Method) =====\n")
print(disp)

capture.output(perm_a, perm_b, disp, file = "stats_output/permanova.txt")

nmds <- metaMDS(rel_t, distance = "bray", k = 2, trymax = 100,
                autotransform = FALSE, trace = 0)
cat(sprintf("\nNMDS stress: %.3f\n", nmds$stress))

nmds_df <- data.frame(scores(nmds, display = "sites"), meta)
region_colours <- c("East Sea" = "#16A085", "Yellow Sea" = "#D4AC0D",
                    "South Sea" = "#8E44AD", "Jeju Sea" = "#A0522D",
                    "Seogwipo Sea" = "#C2185B")
fig_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, colour = Location, shape = Method)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = region_colours, name = NULL) +
  scale_shape_manual(values = c("UUV" = 16, "Net" = 17), name = NULL) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, size = 3.2,
           label = sprintf("Stress = %.3f", nmds$stress)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")
ggsave("stats_output/Figure_S_NMDS.tiff", fig_nmds, width = 16, height = 12,
       units = "cm", dpi = 300, compression = "lzw")
ggsave("stats_output/Figure_S_NMDS.pdf", fig_nmds, width = 16, height = 12,
       units = "cm")

cat("\nAll statistics written to ./stats_output/\n")
