load("out/df10x10_180914.Rdata")
for(drtmp in c("Afatinib", "Navitoclax", "Idelalisib", "Venetoclax")){
  for(pattmp in unique(filter(df10x10, BaseDrugName == drtmp)$PatientID)){
    if(!(drtmp == "Idelalisib" & pattmp == "P0050")){ # fits do not converge in this case - check!
    rmarkdown::render("Exploration_Loewe_ZIP.Rmd",
                      params = list(dr=drtmp, pat=pattmp),
                      output_file = paste0(paste("out_explorations_LoeweZIP/",drtmp,pattmp,sep="_"),".html"))
    }
  }
}