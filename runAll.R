# This file reproduces the whole analysis starting from the raw data by rendering the Rmds.
# start from a new session
# before running the script please download the raw data from EMBL-EBI Biostudies into the data directory of this repository
library(knitr)
rm(list=ls())
today <- "200415"
rmarkdown::render("DrugCombi_DataImport.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_QC.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_AddOmicsData.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_BaseAlone.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_Combi.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_10x10.Rmd", params = list(today=today))
rmarkdown::render("DrugCombi_CombiSynergy.Rmd", params = list(today=today))
rmarkdown::render("Afatinib_targets.Rmd", params = list(today=today))

# store sessionInfo
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
