# UUV-eDNA-plankton-Korean-coast

R scripts for the statistical analyses and figures of a comparative study of plankton diversity detection based on environmental DNA (eDNA) collected by an uncrewed underwater vehicle (UUV) and by conventional net sampling across Korean coastal waters.

## Overview

The study compares two methods at paired stations across five coastal regions bordering the Korean Peninsula: the South Sea, East Sea, Yellow Sea, Jeju Sea, and Seogwipo Sea. Samples were collected in March 2024 (early spring).

- **UUV water sampling:** a tethered, remotely operated UUV (FIFISH V6 EXPERT M100) was launched from shore. It collected 1 L of unfiltered water at 1.0 m depth with an onboard peristaltic pump, and the water was filtered after retrieval through 0.45-µm mixed cellulose ester membranes.
- **Net sampling:** within a 100-m radius, a 20-µm mesh conical plankton net was cast from the same shoreline site and retrieved obliquely by hand. A 100-mL subsample of the concentrated retentate was filtered in the same way.

Both sample types were processed through an identical DNA extraction, 18S rRNA V9 amplification (primers 1380F/1510R), Illumina MiSeq sequencing (2 × 150 bp), and bioinformatic pipeline.

The dataset comprises 125 profiles: 63 UUV samples and 62 net hauls from 63 stations. Sixty-two stations are fully paired; at one South Sea station, the net sample was excluded following sequencing quality-control failure.

## Data

The scripts use one processed dataset, `MetaK90_location.csv`, deposited in Dryad (see the Data Availability Statement of the manuscript).

- **Contents:** a genus × profile read-count matrix for all 125 profiles and 596 plankton genera.
  - Rows 1–4: sampling method, region, latitude, and longitude of each profile.
  - Row 5: profile IDs.
  - Remaining rows: the number of reads assigned to each genus in each profile.
  - First two columns: genus name and functional guild.
- **How it was generated:**
  - Reads were processed in QIIME2 (v2024.10) with Cutadapt (v4.9) and DADA2, yielding 5,538,609 reads and 11,324 ASVs.
  - ASVs were assigned with BLASTN against the NCBI nt database restricted to Eukaryota (≥80% identity, ≥90% query coverage, e-value ≤ 1×10⁻²⁰).
  - Hits to uncultured or environmental sequences, and hits to genera not recorded in the National Species List of Korea (National Institute of Biological Resources), were discarded. Each ASV was assigned to the genus of its highest-identity remaining hit.
  - Only assignments with ≥90% identity to genera classified as phytoplankton or zooplankton were retained (2,616,666 reads).

### Variables

| Variable | Values |
|---|---|
| Sampling method | `Drone` (uncrewed underwater vehicle; UUV), `Net` |
| Region | South Sea, East Sea, Yellow Sea, Jeju Sea, Seogwipo Sea |
| Guild (Taxonomic Group) | `Phytoplankton` (photosynthetic genera, including mixotrophic dinoflagellates); `Zooplankton` (heterotrophic genera, including heterotrophic protists and metazoans) |
| Read counts | Number of sequence reads assigned to each genus in each profile |

## Scripts

| File | Output |
|---|---|
| `Statistics_paired.R` | All statistics reported in the manuscript, written to `./stats_output/`. These cover station-level pairing of UUV and net profiles (62 pairs); per-profile metrics (reads, generic richness, Shannon diversity, phytoplankton read proportion); descriptive statistics; paired Wilcoxon signed-rank tests; library-size correlations; rarefaction to the minimum library size; regional summaries; method-exclusive genera; pooled read fractions; and PERMANOVA, PERMDISP, and NMDS. |
| `Figure_01.R` | Map of sampling stations across the five regions. Each station is classified by pairing status (paired, or net sample excluded following QC failure). |
| `Figure_02.R` | Venn diagram and bar plot of genera detected by both methods or by only one method, and mean per-sample read counts of the 20 most abundant genera in UUV and net profiles. |
| `Figure_03.R` | Heatmap of log10(mean reads per profile + 1) for the 30 most abundant genera, grouped by guild. |
| `Figure_04.R` | Generic richness (total, phytoplankton, and zooplankton) by sampling method and region. |
| `Figure_05.R` | Reads per profile, generic richness, phytoplankton/zooplankton read proportions, Shannon diversity (total and by guild), and the number of genera detected by only one method. |
| `Figure_06.R` | Guild-level comparisons: phytoplankton and zooplankton generic richness (a, c) and mean per-sample read counts of the 20 most abundant phytoplankton and zooplankton genera (b, d). |
| `Figure_07.R` | Relative read abundance (%) of genera in each region (UUV and net profiles pooled). |
| `Figure_08.R` | The 10 most abundant phytoplankton and zooplankton genera, by mean reads per profile, in each region. |
| `Figure_09.R` | Distribution maps of the 20 most abundant genera (a–t), scaled by log-transformed read counts. |
| `Figure_10.R` | Phytoplankton–zooplankton co-occurrence analysis based on within-sample relative abundance. Panel (a) is the nationwide Spearman correlation matrix. Panels (b–e) are the regional networks for the South Sea, East Sea, Yellow Sea, and pooled Jeju + Seogwipo Sea (\|ρ\| > 0.6, Benjamini–Hochberg FDR-adjusted p < 0.01). The script also exports the retained edges and the method-sensitivity results, in which each edge is recalculated within UUV and net profiles separately and by a method-stratified rank correlation. |
| `Figure_S1.R` | Supplementary Figure S1: NMDS ordination of community composition based on Bray–Curtis dissimilarities of genus-level relative read abundance. |

## How to run

1. **Prepare the data:** download `MetaK90_location.csv`, open it in Microsoft Excel or LibreOffice Calc, and save it as `MetaK90_location.xlsx` in your R working directory. The scripts read the Excel format.
2. **Statistics:** run `Statistics_paired.R`. Results are written to `./stats_output/`.
3. **Figures:** run `Figure_01.R` to `Figure_10.R` and `Figure_S1.R`. Each script reads `MetaK90_location.xlsx` and generates the corresponding figure.

On Windows, use forward slashes (`/`) in any file path passed to `setwd()`.

## Requirements

- R 4.4.2 was used; R ≥ 4.3.0 is recommended. RStudio is optional.
- R packages:
  - Data handling: `readxl`, `dplyr`, `tidyr`, `reshape2`
  - Graphics and cartography: `ggplot2` (v3.4.2), `cowplot` (v1.1.1), `patchwork`, `scales`, `ggtext`, `sf` (v1.0-13), `ggspatial`, `rnaturalearth`, `rnaturalearthdata`
  - Community ecology and statistics: `vegan` (v2.7-2) and `permute`, used for rarefaction, Shannon diversity, PERMANOVA, PERMDISP, and NMDS. Paired Wilcoxon signed-rank tests and Spearman correlations use base R (`stats`).
  - Network analysis: `igraph` (v1.4.3), `tidygraph`, `ggraph` (v2.1.0)

```r
install.packages(c("readxl", "dplyr", "tidyr", "reshape2", "ggplot2", "cowplot",
                   "patchwork", "scales", "ggtext", "sf", "ggspatial",
                   "rnaturalearth", "rnaturalearthdata", "vegan", "permute",
                   "igraph", "tidygraph", "ggraph"))
```

## Data source

The data are primary field data collected and analyzed by the authors in March 2024; they were not derived from external or published sources. All field samples were collected under the regulatory guidelines of the relevant Korean coastal management authorities. No specific commercial licensing or proprietary permits were required for the water and plankton sampling.
