# start from a new session, requires connection to the server
library(knitr)
rm(list=ls())
today <- "180903"
rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_DataImport.Rmd", params = list(today=today))
rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_QC.Rmd", params = list(today=today))
rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_AddOmicsData.Rmd", params = list(today=today))
rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_BaseAlone.Rmd", params = list(today=today))
rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_Combi.Rmd", params = list(today=today))
# rmarkdown::render("~/Documents/cll/MarinaDrugComb/DrugCombi_Code4Ms/DrugCombi_CombiSynergy.Rmd", params = list(today=today))
