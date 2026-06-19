R Script Files (Figure_01.R to Figure_10.R)

* These script files contain the complete downstream statistical and visualization pipelines scripted in R.

* By executing these code files sequentially with the prepared dataset, all primary data analysis results and comparative publication-quality graphs (corresponding to Figures 1 through 10) will be automatically generated and outputted.

Files and variables
File: Figure_01.R
Description: Generates the spatial sampling design map displaying the distribution of collection stations across five macro-regions within Korean coastal waters (Yellow Sea, East Sea, South Sea, Jeju Sea, and Seogwipo Sea), distinguishing between UUV-mediated molecular collection and conventional net hauls.

File: Figure_02.R
Description: Generates a Venn diagram and bar plot illustrating the shared and method-specific plankton genus richness between UUV and Net methods, along with the ranking of total ASV yields for the top 20 dominant genera.

File: Figure_03.R
Description: Generates a comparative heatmap presenting the log-transformed sequencing read abundances (log 10 + 1) of major phytoplankton and zooplankton genera detected by UUV versus Net sampling.

File: Figure_04.R
Description: Generates a panel of bar plots, box plots, and scatter plots that compare alpha diversity (genus richness) for both phytoplankton and zooplankton partitioned by sampling method and geographic location.

File: Figure_05.R
Description: Generates box plots and bar charts comparing overall data structures between the two methods, including total ASV counts, detected genus richness, taxonomic ratios (Phytoplankton vs. Zooplankton), Shannon diversity indices, and the number of method-specific unique genera.

File: Figure_06.R
Description: Generates box plots for genus richness and horizontal bar charts of total ASV yields, explicitly separated into independent profiles for Phytoplankton (a, b) and Zooplankton (c, d) under both sampling platforms.

File: Figure_07.R
Description: Generates a stacked bar plot illustrating the relative abundance (%) of major phytoplankton and zooplankton genera pooled across the five macro-regional marine basins.

File: Figure_08.R
Description: Generates horizontal bar charts ranking the top 10 most abundant phytoplankton and zooplankton genera based on total sequence reads (ASVs) within each of the five individual macro-regions.

File: Figure_09.R
Description: Generates macro-regional distribution maps for the top 20 dominant plankton genera (a to t) across the Korean Peninsula, mapping total read abundances at each sampling station.

File: Figure_10.R
Description: Generates an inter-generic Spearman rank correlation matrix heatmap (a) and macro-regional co-occurrence network architectures (b to f) showing the ecological interactions between dominant phytoplankton and zooplankton genera across the five sea areas.

File: MetaK90_location.csv
Description: The primary processed dataset containing spatial coordinates (latitude and longitude), regional classifications, sampling methods (UUV vs. Net), and the absolute sequencing read abundances for each individual plankton genus across all 125 pelagic profiles. This consolidated matrix serves as the direct data source for generating the spatial distribution maps, abundance plots, and ecological network analyses in the R scripts.
