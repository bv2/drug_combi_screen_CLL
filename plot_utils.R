# function to extract legend from a ggplot object
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

myHeatmap <- function(test, minV=NULL, maxV = NULL, paletteLength=100,
                      colors=c("black", "blue", "white", "orange","red"),  ...){
  
  if(is.null(minV)) minV <- min(test)
  if(is.null(maxV)) maxV <- max(test)
  
  myColor <- colorRampPalette(colors)(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(minV, 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(maxV/paletteLength,maxV, length.out=floor(paletteLength/2)))
  
  # Plot the heatmap
  pheatmap(test, color=myColor, breaks=myBreaks, ...)
}

##############################
# Drug Combination Viability Comparisons
#############################

# Function to plot comparison of two combination compunds across all base drugs
plotComparsionCDrugs <- function(df, CDrugAbrv.x, CDrugAbrv.y, range = c(0,1.4), type=c("scatter_joint", "scatter_factBDrug", "boxplot_joint")){

  df.x <- filter(df, CDrugAbrv==CDrugAbrv.x)
  df.y <- filter(df, CDrugAbrv==CDrugAbrv.y)

  # only takes common patient samples of the two df
  df <- merge(df.x, df.y, by=c("PatientID","BDrugID", "BDrugConcId",
                             "BDrugConc", "BDrugName"))
  print(paste0("n_pats=", length(unique(df$PatientID))))
  print(paste0("n_points=", nrow(df)))
  
    if(type == "scatter_joint"){
    gg <- ggplot(df, aes(x=effectBC.x, y=effectBC.y, col=BDrugName)) +
      geom_point() +
      xlab(paste("Combination effect with", CDrugAbrv.x, "\n (viability relative to DMSO control)")) +
      ylab(paste("Combination effect with", CDrugAbrv.y, "\n (viability relative to DMSO control)")) +    
      annotate("text", x=range[1]+0.1, y=range[2]-0.1, label=paste("cor ==", (round(cor(df$effectBC.x,df$effectBC.y),2))), parse=T) +
      geom_abline(slope=1, intercept=0, lty="dashed") + guides(col=guide_legend(title="Base compound"))  +
      theme_bw(base_size=16)

} else if(type == "scatter_factBDrug"){
  gg <- ggplot(df, aes(x=effectBC.x, y=effectBC.y, col=BDrugName)) +
    geom_point() +
    xlab(paste("Combination effect with",CDrugAbrv.x,  "\n (viability relative to DMSO control)")) +
    ylab(paste("Combination effect with",CDrugAbrv.y,  "\n (viability relative to DMSO control)")) +
    facet_wrap(~BDrugName) +
    annotate("text", x=range[1]+0.3, y=range[2]-0.1,
             label=paste("cor ==", (round(summarize(group_by(df, BDrugName), cor=cor(effectBC.x, effectBC.y))$cor,2))), parse=T, size=3) +
    guides(col=guide_legend(ncol=1)) + geom_abline(slope=1, intercept=0)
  
} else if( type=="boxplot_joint"){
  df_rbind <- data.frame(effect=c(df$effectBC.x, df$effectBC.y),
                         Cdrug=c(df$CDrugAbrv.x, df$CDrugAbrv.y),
                         BDrugName = c(df$BDrugName, df$BDrugName))
  t.out <- t.test(effect ~ Cdrug, df_rbind, var.equal=TRUE)
  print(paste("Median in drug",CDrugAbrv.x, round(median(df$effectBC.x),3)))
  print(paste("Median in drug",CDrugAbrv.y, round(median(df$effectBC.y),3)))
  
  gg <- ggplot(df_rbind, aes(x=Cdrug, y=effect)) + 
    ggbeeswarm::geom_beeswarm(alpha=0.7, col="gray", cex=0.7) + 
    geom_boxplot(alpha=0.4, width=0.2, outlier.shape=NA) +
    ggpubr::stat_compare_means(method = "t.test", aes(label =  paste0("p = ",..p.format..), group = Cdrug)) +
    ylab(paste("Combination effect \n (viability relative to DMSO control)")) +
    xlab("Combination compound  \n") +
    theme_bw(base_size=16)
}
  
return(gg)
}


##############################
# Drug Combination Viability Heatmaps
#############################

# Function to plot heatmaps based on combination viability values
# CDrugAbrv - Combination drug to plot
# type patient by patient or patient by drug heatmap?
# useAverage - take avaerage across concentrations specified in conc4average?
# nOcc - number of occurences for genetic aberrations to be annotated 
plotHeatmap <- function(df4ana, dfMuts, CDrugAbrv4plot, type=c("PatPat", "PatDrug", "DrugDrug"), useAverage=FALSE,
                        conc4average = paste0("c",1:5), nOcc = 10, DrugDrugbyIGHV=FALSE, returnMat = FALSE){
  
  dfBPlusC <- filter(df4ana, CDrugAbrv == CDrugAbrv4plot)
  dfBPlusC %<>% mutate(BDrugNameConc = paste(BDrugName, BDrugConcId, sep="_"))
  
  if(!useAverage){
  effectBPlusC_mat <- dfBPlusC %>% 
    select(BDrugNameConc, PatientID, effectBC)  %>%
    spread(key=BDrugNameConc, value=effectBC) %>% 
    column_to_rownames("PatientID") %>%
    as.matrix()
  } else {
    effectBPlusC_mat <- dfBPlusC %>% 
      select(BDrugConcId,BDrugName, PatientID, effectBC)  %>%
      filter(BDrugConcId %in% conc4average) %>%
      group_by(BDrugName,PatientID) %>%
      summarise(effectBC = mean(effectBC, na.rm=TRUE)) %>%
      ungroup() %>%
      spread(key=BDrugName, value=effectBC) %>% 
      column_to_rownames("PatientID") %>%
      as.matrix()
  }
  
  # add genetic background
  patsMCLL <- filter(dfMuts, IGHV == 1)$PatientID
  patsMCLL <- intersect(patsMCLL, rownames(effectBPlusC_mat))
  patsUCLL <- filter(dfMuts, IGHV == 0)$PatientID
  patsUCLL <- intersect(patsUCLL, rownames(effectBPlusC_mat))
  
  dfanno <- dfMuts %>% column_to_rownames("PatientID") %>% select(which(colSums(.) >= nOcc))
  dfanno <- as.data.frame(ifelse(dfanno==1, "mut", "wt"))
  
  # set colors for mutation status
  cols_mut        <- c("white", "black")
  names(cols_mut) <- c("wt", "mut")
  anno_colors <- lapply(colnames(dfanno), function(x) cols_mut)
  names(anno_colors) <- colnames(dfanno)
  
  if(type=="PatPat"){
    #only annotate by IGHV
    dfanno_IGHV <- select(dfMuts, IGHV, PatientID ) %>%
      mutate(IGHV = ifelse(IGHV ==0, "U-CLL", "M-CLL")) %>%
      column_to_rownames("PatientID")
    
    cols_IGHV      <- c("black", "white")
    names(cols_IGHV) <- c("M-CLL", "U-CLL")
    corPat <- cor(t(effectBPlusC_mat), use="complete.obs")
    pheatmap(corPat, annotation_row = dfanno_IGHV, annotation_colors =list(IGHV=cols_IGHV),
             show_rownames = FALSE, show_colnames = FALSE, annotation_legend = FALSE,
             treeheight_col = 12, treeheight_row = 12)
    
  } else if(type == "PatDrug"){
    
    if(!useAverage){
      # outlying values cut off at 1.4
      effectBPlusC_mat[is.na(effectBPlusC_mat)] <- 1.4
    }
      pheatmap(effectBPlusC_mat, na_col="gray", clustering_distance_rows="correlation",
               clustering_distance_cols = "correlation",
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(140),
               breaks=seq(0,filter_th,0.01), show_rownames = FALSE, show_colnames = TRUE,
               treeheight_row = 15, treeheight_col = 15, annotation_row = dfanno, annotation_colors =anno_colors,
               annotation_legend = FALSE, fontsize_col=6)

    
  } else if(type == "DrugDrug"){
      if(DrugDrugbyIGHV) {
        mat <- effectBPlusC_mat[patsMCLL, ]
        corDrug <- cor(mat, use="complete.obs")
        hmMCLL <- pheatmap(corDrug, na_col="gray", clustering_distance_rows="correlation",
               clustering_distance_cols = "correlation",
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
               treeheight_row = 15, treeheight_col = 15, main = "M-CLL",breaks=seq(0,1,0.01))
        mat <- effectBPlusC_mat[patsUCLL, ]
        corDrug <- cor(mat, use="complete.obs")
        hmUCLL <- pheatmap(corDrug, na_col="gray", clustering_distance_rows="correlation",
                 clustering_distance_cols = "correlation",
                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
                 treeheight_row = 15, treeheight_col = 15, main = "U-CLL",
                 breaks=seq(0,1,0.01))
        grid.arrange(hmMCLL$gtable, hmUCLL$gtable, ncol =2)
        
      } else {
        corDrug <- cor(effectBPlusC_mat, use="complete.obs")
        pheatmap(corDrug, na_col="gray", clustering_distance_rows="correlation",
                 clustering_distance_cols = "correlation",
                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(140),
                 treeheight_row = 15, treeheight_col = 15)
      }
      
  } 
  if(returnMat) return(effectBPlusC_mat)
}

##############################
# VOLCANOE
#############################

#function of Gosia for volcanoe plot with some minor edits
# Ycut should be a dataframe containing 
# X - x-calues on volcanoe
# Y - y-values on volcanoe
# optional: Grey - boolean: points to make grey and not label
# Label - labe foe each point

ggvolc = function(df, title ="", Ycut, color=c("deeppink", "navy"), xlab = "") {
  
  # by default label all
  if(!"Grey" %in%  colnames(df)) df$Grey <- FALSE
    
  # check if dataq.frame have required columns
  stopifnot(all(c("X","Y","Label","Grey") %in% colnames(df)))
  
  # check if colors are for the palette # if not use the default
  col.default = c("deeppink", "navy")
  color = ifelse(color %in% colors(), color,  col.default)
  
  # function for colors
  col2hex = function(cols, alpha=1, names=NA) {
    tmp = col2rgb(cols) 
    max=255
    tmp = apply(tmp,2, function(t) rgb(red=t[1], green=t[2], blue=t[3], maxColorValue=max, alpha=alpha*max))
    if(all(!is.na(names)) & length(names)==length(tmp)) tmp = setNames(tmp,nm=names)
    tmp
  }
  
  # xlim
  minX = min(df$X)
  maxX = max(df$X)
  
  # y axis labels
  maxY = max(ceiling(max(df$Y)), 4)  # minimum ylim is set to 4
  axisMark =if(maxY<5) 1 else if(maxY<10) 2 else if(maxY<15) 4 else 5
  
  # direction of the effect
  df$Direction = 0
  df$Direction[df$Y>=Ycut & df$X<0] = -1
  df$Direction[df$Y>=Ycut & df$X>0] = 1
  df$Direction[df$Y>=Ycut & df$Grey] = 2
  df$Direction = factor(df$Direction, levels=c(-1,1,0,2))
  
  # add column saying if association is significant
  df$IsSignificant = df$Y>=Ycut & !df$Grey
  
  ## ONLY IF THERE IS ANY SIGNIFICANT RESULT
  if(sum(df$IsSignificant)>0) {
    
    # hanging ends for labels
    hend = max(abs(c(minX, maxX)))*0.7 #0.42
    
    # shift of labels and lines
    shiftX = max(abs(c(minX, maxX)))*0.05
    
    # positions for labels
    calcY = function(y.org) {
      rng = range(y.org)
      inc = max((rng[2]-rng[1])/length(y.org), 0.1)
      newY = rng[1]+inc*(1:length(y.org))
      y.org[order(y.org)] = newY
      y.org
    }
    df$labY[df$IsSignificant] = calcY(y.org=df$Y[df$IsSignificant])
    df$labX = with(df, ifelse(IsSignificant, ifelse(Direction==1, maxX+shiftX, minX-shiftX), NA))
    df$hjust = ifelse(sign(df$labX)==1, 0, 1)
    
  } else {
    # hanging ends for labels
    hend = 0
  }
  
  gg = ggplot() + geom_point(data=df, aes(x=X, y=Y, fill=Direction, colour=Direction), size=ifelse(df$IsSignificant,3,2), shape=21) + scale_colour_manual(values=c(`-1`=color[1], `1`=color[2], `0`="black", `2`="darkgrey"), guide=FALSE) + scale_fill_manual(values=c(col2hex(c(color[1], color[2]), names=c("-1","1"), alpha=0.7), "0"="black", `2`="grey"), guide=FALSE) + theme_bw() + geom_hline(yintercept=Ycut, colour="navy", linetype="dashed") + geom_vline(xintercept=0, colour="navy", size=0.1) + ylab("") + xlab(xlab) + ggtitle(title) + scale_y_continuous("P-value", breaks=seq(1,maxY,axisMark), labels=10^(-seq(1,maxY,axisMark)), limits=c(0,maxY)) + xlim(minX-hend, maxX+hend)
  
  if(sum(df$IsSignificant)>0) {
    df = df[df$IsSignificant,]
    gg = gg + geom_segment(data=df, aes(x=X, y=Y, xend=labX, yend=labY), colour="darkgrey", alpha=0.7, linetype="dotted") + geom_text(data=df, mapping=aes(x=labX, y=labY, label=Label, colour=Direction, hjust=hjust), size=2.5)
  }
  gg
}


#############################
# Drug Combination Synergy
#############################

# plot the single drug responses and the response to the combination as well as the additive effect as curves (mean + SE) for a pairs of drugs drB and drC.
plotResponseCurves <- function(df, drC , drB, th = filter_th, sep_by_IGHV =FALSE, sep_by_TP53=FALSE){
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName == drB) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c",5:1))) %>%
    select(CDrugAbrv, BDrugName, BDrugConcId, starts_with("viab"), PatientID) %>%
    gather(key="type", value = "viability", starts_with("viab")) %>%
    mutate(type = ifelse(type == "viabB", paste(drB, "(A)"),
                         ifelse(type == "viabC", paste(drC, "(B)"),
                                ifelse(type == "viabBC", paste(drB, "+", drC, "(AB)"),
                                       "additive effect (A*B)"))))
  df4plot <- left_join(df4plot, dfMuts4testing, by ="PatientID")
  df4plot %<>% mutate(IGHV = ifelse(IGHV == 0, "U-CLL", "M-CLL"))
  df4plot %<>% mutate(TP53 = ifelse(TP53 == 0, "TP53-wt", "TP53-mut"))
  
 
  gg <- ggplot(data=df4plot, aes(x=BDrugConcId, y=viability, col =type, group=type)) +
    stat_summary(fun.data = "mean_se") +
    stat_summary(fun.y = "mean", geom="line") + 
    theme_bw(base_size = 20) + xlab(paste0("Concentration of ", drB)) +
    guides(col = guide_legend(title="")) + ylim(c(0,th))
  
  if(sep_by_IGHV & !sep_by_TP53){
    gg <- gg + facet_wrap(~IGHV)
  } else if(!sep_by_IGHV & sep_by_TP53){
    gg <- gg + facet_wrap(~TP53)
  } else if(sep_by_IGHV & sep_by_TP53){
    gg <- gg + facet_wrap(TP53~IGHV)
  }
  
  return(gg)
}

# plot the single drug responses and the response to the combination as well as the additive effect as curves (mean + SE) for a pairs of drugs drB and drC.
plotMultipleBResponseCurves <- function(df, drC , drsB, th = filter_th){
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName %in% drsB) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c",5:1))) %>%
    select(CDrugAbrv, BDrugName, BDrugConcId, starts_with("viab"), PatientID) %>%
    mutate(BDrugName = factor(BDrugName, levels =drsB)) %>% #keep order as specified
    gather(key="type", value = "viability", starts_with("viab")) %>%
    mutate(type = ifelse(type == "viabB", paste("Base compound", "(A)"),
                         ifelse(type == "viabC", paste(drC, "(B)"),
                                ifelse(type == "viabBC", paste("Base compound", "+", drC, "(AB)"),
                                       "additive effect (A*B)"))))
  df4plot <- left_join(df4plot, dfMuts4testing, by ="PatientID")
  df4plot %<>% mutate(IGHV = ifelse(IGHV == 0, "U-CLL", "M-CLL"))
  
  
  gg <- ggplot(data=df4plot, aes(x=BDrugConcId, y=viability, col =type, group=type)) +
    stat_summary(fun.data = "mean_se") +
    stat_summary(fun.y = "mean", geom="line") + 
    theme_bw(base_size = 20) + xlab(paste0("Concentration")) +
    guides(col = guide_legend(title="")) + ylim(c(0,th)) +
    facet_wrap(~ BDrugName) +theme(legend.position = "top") +
    guides(col=guide_legend(nrow=2,byrow=TRUE, title = ""))
  
  
  return(gg)
}

# plot Boxplot of CI for drug - drug combination
plotBoxplotCI <- function(df, drC , drB, CI_type = c("Bliss", "hsa", "SI")){
  CI_type <- match.arg(CI_type)
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName == drB) %>% 
    mutate(BDrugConcId =factor(BDrugConcId, levels = paste0("c",5:1)))
  if(CI_type == "Bliss"){
    df4plot$CI <- df4plot$BlissCI
  } else if(CI_type == "hsa") {
    df4plot$CI <- df4plot$hsaCI
  } else {
    df4plot$CI <- df4plot$addModelSI
  }

  ggplot(df4plot, aes(x=BDrugConcId, y=CI)) +
    ggbeeswarm::geom_beeswarm(col = "gray", alpha =0.7) + 
    geom_boxplot(outlier.shape = NA, alpha=0.2, width = 0.3) +
    geom_hline(yintercept = ifelse(CI_type =="SI",0,1), lty = "dashed") +
    theme_bw()  + xlab(paste0("Concentration of ", drB)) + ylab(paste0("Combination index (", CI_type, ")"))
}

# plot bar plot of individual combination indices for each sample
plotWaterfallCI <-  function(df, drC , drB, CI_type = c("Bliss", "hsa", "SI"), annotate=NULL){
  
  CI_type <- match.arg(CI_type)
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName == drB) 
    
  if(CI_type == "Bliss"){
    df4plot$CI <- df4plot$BlissCImean
    df4plot$CIse <- df4plot$BlissCIse
  } else if(CI_type == "hsa") {
    df4plot$CI <- df4plot$hsaCImean
    df4plot$CIse <- df4plot$hsaCIse
  } else {
    df4plot$CI <- df4plot$addModelSImean
    df4plot$CIse <- df4plot$addModelSIse
  }
  
  # annotate by genetic features
  if(!is.null(annotate)){
    df4plot %<>% left_join(dfMuts4testing, by = "PatientID") 
  }
  
  # order patients based on the CI value
  df4plot %<>% arrange(CI)
  df4plot$PatientID = factor(df4plot$PatientID, levels = df4plot$PatientID)
  

  gg <- ggplot(df4plot, aes(x=PatientID, y=CI, fill=PatientID)) +
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymin = CI - CIse, ymax = CI + CIse), width = 0.3)+
    geom_hline(yintercept = ifelse(CI_type =="SI",0,1), lty = "dashed") +
    theme_bw()  + xlab(paste0("Patient sample")) +
    ylab(paste0("Combination index (", CI_type, ")")) + 
    scale_fill_manual(values = patcol) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(fill=FALSE)
  
  # annotate by genetic features
  if(!is.null(annotate)){
    gg <- gg + ylim(c(min(df4plot$CI - df4plot$CIse), max(df4plot$CI + df4plot$CIse) +0.04))
    if("IGHV" %in% annotate){
    gg <- gg + geom_text(aes(y = max(df4plot$CI + df4plot$CIse) +0.01, label = ifelse(IGHV==1, "*", "")), size=5)
    gg <- gg + geom_text(x=nrow(df4plot)-10, y= min(df4plot$CI), label = "* M-CLL", col="black", size=5)
    }
    if("TP53" %in% annotate){
    gg <- gg + geom_text(aes(y = max(df4plot$CI + df4plot$CIse) +0.03, label = ifelse(TP53==1, "*", "")), col="red", size=5)
    gg <- gg + geom_text(x=nrow(df4plot)-3, y= min(df4plot$CI), label = "* TP53 mut", col="red", size=5)
    }
  }
  
  return(gg)
}

# dotplots add therotical effect vs measured combination effect
plotScattter <- function(df, drB, drC, th = filter_th){
  
  df4plot <- filter(df, BDrugName==drB)
  range = c(0, th)

  gg <- ggplot(df4plot, aes(x=viabBC_add, y=viabBC, color=PatientID))+
    geom_point(alpha=0.7)+
    geom_hline(aes(yintercept=1), colour="grey", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="grey", linetype="dashed") +
    scale_color_manual(values=patcol)+facet_wrap(~BDrugConcId  , ncol=5) +
    geom_abline(intercept = 0, slope = 1, colour="black", linetype="solid") +
    # ggtitle(paste("Base drug (A):", drB, "\n Combination drug (B):", drC)) +
    coord_fixed() + scale_x_continuous(limits=range) + 
    scale_y_continuous( limits=range) +
    theme_bw(base_size = 14) +
    ylab(paste0("viability (combination)")) +
    xlab(paste0("viability (", drB,") * viability (", drC, ")")) +
    annotate("text", x=0.8*th, y=0.3, label= "synergy", size=4, alpha=0.4) + 
    guides(col=FALSE)
  
  return(gg)
}

# dotplots base compound effect vs measured combination effect
plotScattterVsCombi <- function(df, drB, drC, th = filter_th){
  
  df4plot <- filter(df, BDrugName==drB)
  range = c(0, th)
  
  gg <- ggplot(df4plot, aes(y=viabC, x=viabBC, color=PatientID))+
    geom_point(alpha=0.7)+
    geom_hline(aes(yintercept=1), colour="grey", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="grey", linetype="dashed") +
    scale_color_manual(values=patcol) +
    facet_wrap(~BDrugConcId  , ncol=5) +
    geom_abline(intercept = 0, slope = 1, colour="black", linetype="solid") +
    # ggtitle(paste("Base drug (A):", drB, "Combination drug (B):", drC)) + 
    coord_fixed() + scale_x_continuous(limits=range) + 
    scale_y_continuous( limits=range) +
    theme_bw(base_size = 14) +
    ylab(paste0("viability (", drC, ")")) +
    xlab(paste0("viability (", drB,"+", drC, ")")) +
    annotate("text", y=0.8*th, x=0.3, label= "additivity", size=4, alpha=0.4) + 
    guides(col=FALSE)
  
  return(gg)
}

#############################
# 10x10 screens
#############################

plotTiles10x10 <- function(df, drB, pat, colviab = colorRampPalette(c("darkgreen", "red"))(13)){

    df4plot <- filter(df, PatientID==pat, BaseDrug == drB, CombiDrug == "Ibrutinib") 
    df4plot %<>% mutate(base = paste(BaseDrug, concBaseDrug, sep="_"))
    df4plot %<>% mutate(combi = paste(CombiDrug, concCombiDrug, sep="_"))
    df4plot %<>%  select(combi, base, normalizedValue)
    levelBase <- unique(df4plot$base)
    levelCombi <- unique(df4plot$combi)
    gg <- ggplot(df4plot, aes(x=factor(base,levels=levelBase), y=factor(combi, levels=levelCombi), fill=normalizedValue)) +
        geom_tile() + 
        ggtitle(paste(drB, pat)) +
        scale_fill_gradient(limits=c(0,1.4), breaks=seq(0,1.2,0.2))
      
    return(gg)
    }  
