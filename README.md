# 2026 BSH Substrate Screening

This repository contains the source data and R scripts used to reproduce **Figure 3** and **Supplementary Figure S7** from the BSH/T enzyme substrate-screening study. 
All analysis scripts are provided in the `Code/` folder and the corresponding input/source data are provided in the `Source_data/` folder. Scripts should be run from the root directory of the repository.

## Figure 3 – Enzyme-substrate landscape of bile acid-amine conjugates
Figure 3 summarizes the substrate preference landscape of the screened BSH/T enzymes across bile acid-amine conjugates. Replicate-level peak areas from the enzyme screen were summarized across triplicate measurements. For each enzyme-product pair, the mean peak area was calculated and log10 transformed for visualization.
Rows represent individual BSH/T enzymes and columns represent detected bile amidate products. Enzymes are annotated by bacterial phylum and predicted signal peptide length. Hierarchical clustering was performed using Euclidean distance and complete linkage.
Bile amidate nomenclature uses **Mono-, Di-, or Tri-** to indicate the bile acid core, followed by the conjugated amine. Numerical identifiers indicate cases in which multiple isomeric products were detected but could not be structurally resolved further.
Black outlines identify robust enzyme-substrate pairs meeting both criteria:
`CV < 0.4` and `log10 mean peak area > 1`
White cells indicate no detectable signal.

### Figure 3 legend
**Figure 3. Enzyme-substrate landscape of BBAs across bacterial phyla.** Heatmap of log10 mean peak areas for bile acid-amine conjugates across BSH/T enzymes in triplicate. Rows are enzymes colored by phylum and annotated with signal peptide length; columns are bile amidate products. Nomenclature for products are Mono-, Di-, or Tri-, representing the bile acid core followed by the amine added, and the numerical ID represents that the BBA has more than one isomer detected in the assays that cannot be structurally resolved further. Black outlines mark robust enzyme-substrate pairs (CV < 0.4, log10 mean > 1). White indicates no detectable signal.

## Supplementary Figure S7 – Frequency of commonly detected BBAAs
Supplementary Figure S7 summarizes how frequently individual BBAAs were detected across the enzyme screen. A BBAA was considered detected for an enzyme when it was observed in at least one replicate. Detection frequency was calculated as the percentage of screened enzymes for which each BBAA was detected. The ten most frequently detected products were then ranked and visualized as a bar plot.

### Supplementary Figure S7 legend
**Supplementary Figure S7. Frequency of the most commonly detected BBAs across the enzyme screen.** Bar plot showing the top 10 BBAAs ranked by detection frequency. Detection frequency was calculated as the percentage of screened enzymes for which a given BBAA was detected in at least one replicate. Percentages are indicated at the end of each bar.
