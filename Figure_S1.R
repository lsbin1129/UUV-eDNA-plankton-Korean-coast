# =============================================================================
# Figure_S1.R
#   Figure S1: NMDS ordination of plankton community composition
#   (Bray-Curtis dissimilarity on genus-level relative read abundance)
#
#   Input : MetaK90_location.xlsx (same file used by Statistics_paired.R)
#           row 1 = sampling method (Drone / Net), row 2 = region,
#           row 5 = sample IDs, rows 6- = genera x samples read counts
#   Output: Figure_S1_NMDS.tiff (300 dpi) and Figure_S1_NMDS.pdf
#
#   PERMANOVA / PERMDISP statistics reported with this figure (Results 3.2)
#   are produced by Statistics_paired.R (section 9).
# =============================================================================
# install.packages(c("readxl", "vegan", "ggplot2"))

library(readxl)
library(vegan)
library(ggplot2)

# Set the working directory to the folder that contains MetaK90_location.xlsx
# (use forward slashes "/" in Windows paths)
# setwd("C:/path/to/folder")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
raw <- read_excel("MetaK90_location.xlsx", col_names = FALSE,
                  .name_repair = "minimal")

methods   <- as.character(unlist(raw[1, 3:ncol(raw)]))
locations <- as.character(unlist(raw[2, 3:ncol(raw)]))
samples   <- as.character(unlist(raw[5, 3:ncol(raw)]))
genus_vec <- as.character(unlist(raw[6:nrow(raw), 1]))

reads_matrix <- apply(raw[6:nrow(raw), 3:ncol(raw)], 2,
                      function(x) as.numeric(as.character(x)))
reads_matrix[is.na(reads_matrix)] <- 0
rownames(reads_matrix) <- genus_vec
colnames(reads_matrix) <- samples

method_lab <- ifelse(methods == "Drone", "UUV", "Net")
loc_levels <- c("South Sea", "East Sea", "Yellow Sea", "Jeju Sea", "Seogwipo Sea")

cat("Genera :", nrow(reads_matrix), "\n")
cat("Profiles:", ncol(reads_matrix),
    "(UUV", sum(method_lab == "UUV"), "/ Net", sum(method_lab == "Net"), ")\n")

# -----------------------------------------------------------------------------
# 2. Relative abundance and NMDS (Bray-Curtis)
# -----------------------------------------------------------------------------
rel_t <- t(sweep(reads_matrix, 2, colSums(reads_matrix), "/"))   # samples x genera

set.seed(2024)
nmds <- metaMDS(rel_t, distance = "bray", k = 2, trymax = 100,
                autotransform = FALSE, trace = 0)
cat(sprintf("NMDS stress: %.3f\n", nmds$stress))

nmds_df <- data.frame(scores(nmds, display = "sites"),
                      Location = factor(locations, levels = loc_levels),
                      Method   = factor(method_lab, levels = c("UUV", "Net")))
write.csv(nmds_df, "Figure_S1_NMDS_scores.csv")

# -----------------------------------------------------------------------------
# 3. Plot
# -----------------------------------------------------------------------------
region_colours <- c("South Sea"    = "#8E44AD",
                    "East Sea"     = "#16A085",
                    "Yellow Sea"   = "#D4AC0D",
                    "Jeju Sea"     = "#A0522D",
                    "Seogwipo Sea" = "#C2185B")

fig_s1 <- ggplot(nmds_df, aes(NMDS1, NMDS2, colour = Location, shape = Method)) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = region_colours, name = NULL) +
  scale_shape_manual(values = c("UUV" = 16, "Net" = 17), name = NULL) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6, size = 3.2,
           label = sprintf("Stress = %.3f", nmds$stress)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

ggsave("Figure_S1_NMDS.tiff", fig_s1, width = 16, height = 12,
       units = "cm", dpi = 300, compression = "lzw")
ggsave("Figure_S1_NMDS.pdf", fig_s1, width = 16, height = 12, units = "cm")

cat("Saved: Figure_S1_NMDS.tiff, Figure_S1_NMDS.pdf, Figure_S1_NMDS_scores.csv\n")
