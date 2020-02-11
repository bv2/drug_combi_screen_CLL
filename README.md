This repostitory contains the code for the analysis of the drug combination screen on CLL from Lukas, Velten et al.

The raw data files for this analysis are stored in `/g/huber/projects/nct/cll/RawData/DrugScreens/Marina_CombiScreen`. 

The processed data objects after running the preprocessing scripts can be found in `/g/huber/projects/nct/cll/ProcessedData/DrugScreens/Marina_CombiScreen`.

The analysis contains the following parts:

*A. Pre-processing*
- `DrugCombi_DataImport.Rmd`: This script takes the raw data files as input and generates data object containing the normalized viability values for each well as well as a data frame combining the single effect of each drug and the effects of each drug-drug combination
- `DrugCombi_QC.Rmd`: Some quality control on the data (plate plots, reproducibility between replicates ...). Replicates are collapsed to their average values and a filtering threshold is applied on the noramlized viability values
- `DrugCombi_AddOmicsData.Rmd`: Check available omics data for the patients included in the combi screen and show an overview of genetic heterogeneity

*B. Visualization of base drug effects*
- `DrugCombi_BaseAlone.Rmd`: Analyse effect of base compounds alone, associations to most frequent mutations and heatmaps

*C. Analysis of drug combination effect*
- `DrugCombi_Combi.Rmd`: Analyse effect of drug combinations, associations to most frequent mutations and heatmaps
- `DrugCombi_CombiSynegy.Rmd`: Look for synergisitc effect
- `DrugCombi_10x10.Rmd`: analysis of 10x10 screens

The script `plot_utils.R` contains some functions for plottings. `runAll.R` runs the whole analysis starting from the raw data by rendering all the .Rmd files above.
