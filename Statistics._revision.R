# =============================================================================
# Statistics_revision.R
#   Additional analyses added during revision (Methods 2.4-2.5)
#
#   1. Station-level pairing of UUV and Net profiles (same rule as
#      Statistics_paired.R; 62 pairs expected)
#   2. Macroalgae sensitivity analysis (Methods 2.4, Results 3.1):
#      phytoplankton comparisons repeated after excluding 97 benthic
#      macroalgal genera (red, brown, and green macroalgae and
#      filamentous algae), paired Wilcoxon signed-rank tests
#   3. Method-exclusive phytoplankton genera split into microalgae and
#      macroalgae (Results 3.1)
#   4. Region-specific paired comparisons (Supplementary Table S2):
#      total, phytoplankton, and zooplankton richness, Shannon diversity,
#      and phytoplankton read proportion, tested within each region
#
#   All outputs are written to ./stats_output/
#   Input: MetaK90_location.xlsx (see README.md)
# =============================================================================
# install.packages(c("readxl", "dplyr"))

library(readxl)
library(dplyr)

dir.create("stats_output", showWarnings = FALSE)

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
tax_vec   <- as.character(unlist(raw[6:nrow(raw), 2]))

reads_matrix <- apply(raw[6:nrow(raw), 3:ncol(raw)], 2,
                      function(x) as.numeric(as.character(x)))
reads_matrix[is.na(reads_matrix)] <- 0
rownames(reads_matrix) <- genus_vec

method_lab <- ifelse(methods == "Drone", "UUV", "Net")
is_phyto   <- tax_vec == "Phytoplankton"
is_zoo     <- tax_vec == "Zooplankton"

# =============================================================================
# 1. Station-level pairing (identical to Statistics_paired.R)
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
pairs$Location <- locations[pairs$u]
cat(sprintf("Paired stations: %d\n", nrow(pairs)))
print(table(pairs$Location))

# Two-sided paired Wilcoxon signed-rank test. Zero differences are dropped
# (Wilcoxon's method). The p-value is
#   - exact, when there are no zero and no tied differences (n <= 50);
#   - from complete enumeration of all 2^n sign permutations, when there are
#     ties or zeros and n <= 13 (Jeju Sea, Seogwipo Sea);
#   - otherwise from the normal approximation with tie correction and
#     without continuity correction.
signed_rank_p <- function(d) {
  n  <- length(d)
  nz <- d[d != 0]
  m  <- length(nz)
  if (m == 0) return(NA_real_)
  has_ties <- any(duplicated(abs(nz)))
  if (n <= 50 && !has_ties && m == n) {
    return(wilcox.test(nz, mu = 0, exact = TRUE)$p.value)
  }
  if (n <= 13) {
    r     <- rank(abs(nz))
    obs   <- sum(r[nz > 0])
    signs <- as.matrix(expand.grid(rep(list(c(0, 1)), m)))
    null  <- as.vector(signs %*% r)
    eps   <- 1e-9
    return(min(1, 2 * min(mean(null >= obs - eps), mean(null <= obs + eps))))
  }
  suppressWarnings(wilcox.test(nz, mu = 0, exact = FALSE, correct = FALSE)$p.value)
}

# Returns medians, IQRs, median difference, and p-value for the paired
# profiles in `pr`, plus medians and means across all profiles by method
paired_test <- function(v, pr) {
  x <- v[pr$u]; y <- v[pr$n]
  p <- signed_rank_p(x - y)
  data.frame(n_pairs = nrow(pr),
             UUV_median = median(x),
             UUV_Q1 = unname(quantile(x, 0.25)), UUV_Q3 = unname(quantile(x, 0.75)),
             Net_median = median(y),
             Net_Q1 = unname(quantile(y, 0.25)), Net_Q3 = unname(quantile(y, 0.75)),
             Median_difference = median(x - y),
             UUV_median_all = median(v[uuv_i]), Net_median_all = median(v[net_i]),
             UUV_mean_all = mean(v[uuv_i]), Net_mean_all = mean(v[net_i]),
             p_value = p)
}

shannon_H <- function(v) {
  v <- v[v > 0]
  if (length(v) == 0) return(0)
  p <- v / sum(v)
  -sum(p * log(p))
}

# =============================================================================
# 2. Macroalgae sensitivity analysis
#    Benthic macroalgal genera: Florideophyceae, Bangiophyceae,
#    Compsopogonophyceae, Phaeophyceae, Ulvophyceae, and the filamentous
#    algae Klebsormidium, Spirogyra, Stigeoclonium, and Xanthonema.
# =============================================================================
macroalgal_genera <- c(
  "Acrothrix", "Aglaothamnion", "Ahnfeltiopsis", "Amphiroa", "Anotrichium",
  "Antithamnion", "Antithamnionella", "Bangia", "Batrachospermum",
  "Blidingia", "Bolbocoleon", "Callophyllis", "Capsosiphon", "Caulacanthus",
  "Ceramium", "Ceratodictyon", "Champia", "Chondria", "Chorda", "Cladophora",
  "Colaconema", "Collinsiella", "Corallina", "Cryptopleura", "Dasya",
  "Desmarestia", "Dichotomaria", "Dictyota", "Distromium", "Dumontia",
  "Ectocarpus", "Ellisolandia", "Erythrocladia", "Gelidiophycus", "Gelidium",
  "Gloiopeltis", "Gloiosiphonia", "Gracilaria", "Grateloupia", "Halothrix",
  "Helminthocladia", "Herposiphonia", "Hildenbrandia", "Hypnea", "Jania",
  "Klebsormidium", "Kuckuckia", "Lithophyllum", "Lithothamnion",
  "Lomentaria", "Mazzaella", "Membranoptera", "Mesophyllum", "Monostroma",
  "Myrionema", "Neorhodomela", "Neosiphonia", "Osmundea", "Padina",
  "Palmaria", "Peyssonnelia", "Phaeophila", "Pikea", "Plocamium",
  "Pneophyllum", "Polysiphonia", "Porolithon", "Porphyra", "Protomonostroma",
  "Pseudendoclonium", "Pterocladiella", "Pterosiphonia", "Pyropia",
  "Ralfsia", "Rhizoclonium", "Rhodochorton", "Rhodophysema", "Rugulopteryx",
  "Saccharina", "Scytosiphon", "Sonderophycus", "Spatoglossum",
  "Sphacelaria", "Spirogyra", "Spyridia", "Stenogramma", "Stigeoclonium",
  "Symphyocladia", "Tichocarpus", "Ulothrix", "Ulva", "Ulvaria", "Umbraulva",
  "Urospora", "Wrangelia", "Xanthonema", "Zonaria"
)


is_macro <- is_phyto & genus_vec %in% macroalgal_genera
is_micro <- is_phyto & !is_macro
cat(sprintf("\nMacroalgal genera in dataset: %d of %d phytoplankton genera\n",
            sum(is_macro), sum(is_phyto)))
missing <- setdiff(macroalgal_genera, genus_vec[is_phyto])
if (length(missing) > 0)
  cat("  Not found as Phytoplankton in data:", paste(missing, collapse = ", "), "\n")

pres <- reads_matrix > 0
tot  <- colSums(reads_matrix)
sens_metrics <- list(
  "Phytoplankton richness (all genera)"         = colSums(pres[is_phyto, ]),
  "Phytoplankton richness (excl. macroalgae)"   = colSums(pres[is_micro, ]),
  "Macroalgal richness"                         = colSums(pres[is_macro, ]),
  "Phytoplankton read proportion (%, all)"      = 100 * colSums(reads_matrix[is_phyto, ]) / tot,
  "Microalgal read proportion (%, excl. macroalgal reads)" =
    100 * colSums(reads_matrix[is_micro, ]) / (tot - colSums(reads_matrix[is_macro, ])),
  "Macroalgal read proportion (%)"              = 100 * colSums(reads_matrix[is_macro, ]) / tot,
  "Total richness (excl. macroalgae)"           = colSums(pres[!is_macro, ])
)
# Paired medians (62 pairs) and medians/means across all profiles
# (63 UUV, 62 Net; the latter are the values quoted in Results 3.1)
sens <- bind_rows(lapply(names(sens_metrics), function(k)
  cbind(Metric = k, paired_test(sens_metrics[[k]], pairs))))
print(sens, digits = 3)
write.csv(sens, "stats_output/macroalgae_sensitivity.csv", row.names = FALSE)

# =============================================================================
# 3. Method-exclusive genera by guild, with phytoplankton split into
#    microalgae and macroalgae
# =============================================================================
in_uuv <- rowSums(pres[, uuv_i, drop = FALSE]) > 0
in_net <- rowSums(pres[, net_i, drop = FALSE]) > 0
excl <- bind_rows(lapply(list(
  "Phytoplankton (all)" = is_phyto,
  "Microalgae"          = is_micro,
  "Macroalgae"          = is_macro,
  "Zooplankton"         = is_zoo), function(m)
    data.frame(UUV_only = sum(m & in_uuv & !in_net),
               Net_only = sum(m & in_net & !in_uuv),
               Both     = sum(m & in_uuv & in_net))), .id = "Group")
print(excl)
write.csv(excl, "stats_output/macroalgae_exclusive_genera.csv", row.names = FALSE)
write.csv(data.frame(Genus = genus_vec[is_macro],
                     UUV_only = (in_uuv & !in_net)[is_macro],
                     Net_only = (in_net & !in_uuv)[is_macro]),
          "stats_output/macroalgal_genera_list.csv", row.names = FALSE)

# =============================================================================
# 4. Region-specific paired comparisons (Supplementary Table S2)
#    Jeju Sea and Seogwipo Sea have too few pairs for a meaningful test;
#    their p-values are reported for completeness only.
# =============================================================================
reg_metrics <- list(
  "Total richness"                    = colSums(pres),
  "Phytoplankton richness"            = colSums(pres[is_phyto, ]),
  "Zooplankton richness"              = colSums(pres[is_zoo, ]),
  "Shannon diversity"                 = apply(reads_matrix, 2, shannon_H),
  "Phytoplankton read proportion (%)" = 100 * colSums(reads_matrix[is_phyto, ]) / tot
)
regions <- c("South Sea", "East Sea", "Yellow Sea", "Jeju Sea", "Seogwipo Sea",
             "All regions")
reg <- bind_rows(lapply(regions, function(R) {
  pr <- if (R == "All regions") pairs else pairs[pairs$Location == R, ]
  bind_rows(lapply(names(reg_metrics), function(k) {
    out <- paired_test(reg_metrics[[k]], pr)
    cbind(Region = R, Metric = k, out[, !grepl("_all$", names(out))])
  }))
}))
print(reg, digits = 3)
write.csv(reg, "stats_output/regional_paired_wilcoxon.csv", row.names = FALSE)
