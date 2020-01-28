library(tidyverse)
vaf <- read.csv("~/Google Drive/CombiScreen/revision/add_info/TP53_VAF.csv", sep = ";", dec = ",")
minvaf <- min(vaf$VAF, na.rm = T)

pdf("~/Google Drive/CombiScreen/revision/add_info/VAF.pdf", height = 4)
vaf %>% mutate(VAF = as.numeric(VAF)) %>% ggplot(aes(x=VAF)) +
    geom_histogram(fill = "navy", bin =100) + theme_bw() + geom_vline(xintercept = minvaf, col = "red", lty = 2)
dev.off()