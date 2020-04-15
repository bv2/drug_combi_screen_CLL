This repository contains the code for the analysis of the drug combination screen on CLL from Lukas, Velten et al. *Survey of ex vivo drug combination effects in chronic lymphocytic leukemia reveals synergistic drug effects and genetic dependencies*.

The raw data files for this analysis are stored in the EMBL-EBI Biostudies repository (accession number XXX).

The processed data objects can be reproduced by running the preprocessing scripts and are used for the visualization in the [Shiny App](http://mozi.embl.de/public/combiScreen/).

The analysis contains the following parts:

*A. Pre-processing*
- `DrugCombi_DataImport.Rmd`: This script takes the raw data files as input and generates data object containing the normalized viability values for each well. Furthermore, it generates a data frame combining the single effect of each drug and the effects of each drug-drug combination.
- `DrugCombi_QC.Rmd`: Some quality control on the data (plate plots, reproducibility between replicates, etc.). Replicates are collapsed to their average values and a filtering threshold is applied on the normalized viability values to remove outliers.
- `DrugCombi_AddOmicsData.Rmd`: Check available omics data for the patients included in the combination screen and show an overview of genetic heterogeneity.

*B. Visualization of base drug effects*
- `DrugCombi_BaseAlone.Rmd`: Analyse effect of base compounds alone and their associations to most frequent mutations.

*C. Analysis of drug combination effect*
- `DrugCombi_Combi.Rmd`: Analyse effect of drug combinations, associations to most frequent mutations and heatmaps.
- `DrugCombi_CombiSynegy.Rmd`: Investigate synergistic effects
- `DrugCombi_10x10.Rmd`: Analysis of 10x10 screens

The script `plot_utils.R` contains some helper functions for recurrent plots. The script `runAll.R` runs the whole analysis starting from the raw data by rendering all the .Rmd files above.
The script `Afatinib_targets.Rmd` investigates the expression values of potential tagets of Afatinib.

To reproduce the analysis you need to clone the repository and download the data from the EMBL-EBI Biostudies repository into the data directory. Afterwards, the script `runAll.R` can be used to reproduce the full analysis and generate the processed data objects, figures and tables used in the manuscript. The file `sessionInfo.txt` gives details about the package versions that were used to reproduce the analysis. 
