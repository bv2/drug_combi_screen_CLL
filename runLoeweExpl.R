load("out/df10x10_180914.Rdata")
for(drtmp in c("Afatinib", "Navitoclax", "Idelalisib", "Venetoclax")){
  for(pattmp in unique(filter(df10x10, BaseDrugName == dr)$PatientID)){
    rmarkdown::render("Exploration_Loewe_ZIP.Rmd",
                      params = list(dr=drtmp, pat=pattmp),
                      output_file = paste0(paste(drtmp,pattmp,sep="_"),".html"))
  }
}