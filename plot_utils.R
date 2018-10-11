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


plotBaseResponseCurves <- function(dfBmuts, drugs2plot, gene = "TP53", round2=3){
  df4plot <- dfBmuts %>% filter(mutation == gene, BDrugName %in% drugs2plot) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c", 5:1))) %>%
    mutate(BDrugConc = factor(round(BDrugConc, round2))) %>%
    mutate(BDrugName = factor(BDrugName, levels = drugs2plot)) # order as specified
  
  if(gene == "TP53") df4plot %<>% mutate(status = ifelse(status ==0, "wt","mut"))
  else if(gene == "IGHV") df4plot %<>% mutate(status = ifelse(status ==0, "U-CLL","M-CLL"))
    
  gg <- ggplot(df4plot, aes(x=BDrugConc, y=effectB, group=status)) +
    stat_summary(fun.data = "mean_se", aes(col=factor(status)), geom="line", , fun.args = list(mult = 2)) + facet_wrap(~BDrugName, ncol=3, scales = "free_x") +
    stat_summary(fun.data = "mean_se",aes(col=factor(status)), geom="errorbar", width=0.2, fun.args = list(mult =2)) +
    ggpubr::stat_compare_means(method = "t.test", aes(group=status, x=BDrugConc, label =  ..p.signif..), hide.ns = TRUE) +
    guides(col = guide_legend(title=gene))+ xlab("Concentration (µM)") +ylab("viability")+theme_bw(base_size = 14) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, vjust=1, hjust =1))
  # if(length(unique(drugs2plot)) <4) gg <- gg + theme(legend.position = "top")
  if(gene == "IGHV") gg  <- gg +scale_color_manual(values= c("M-CLL" = "red", "U-CLL" = "blue"))
  else if(gene == "TP53") gg  <- gg +scale_color_manual(values= c("wt" = "black", "mut" = "orange"))
  
  gg
}

plotBaseBoxplots <- function(dfBmuts, drugs2plot, gene = "TP53", round2=3){
  df4plot <- dfBmuts %>% filter(mutation == gene, BDrugName %in% drugs2plot) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c", 5:1))) %>%
    mutate(BDrugConc = factor(round(BDrugConc, round2))) %>%
    mutate(BDrugName = factor(BDrugName, levels = drugs2plot)) # order as specified
  
    if(gene == "TP53") df4plot %<>% mutate(status = ifelse(status ==0, "wt","mut"))
    else if(gene == "IGHV") df4plot %<>% mutate(status = ifelse(status ==0, "U-CLL","M-CLL"))
  
    gg <- ggplot(df4plot, aes(x=BDrugConc, y=effectB)) + geom_boxplot(aes(fill=factor(status))) + facet_wrap(~BDrugName, ncol=3 , scales = "free_x") +
    ggpubr::stat_compare_means(method = "t.test", aes(group=status, label =  ..p.signif..), hide.ns = TRUE) +
    guides(fill = guide_legend(title=gene)) + xlab("Concentration (µM)") +ylab("viability") +theme_bw(base_size = 15) +
      theme(strip.background = element_blank(),
            axis.text.x = element_text(angle=90, vjust=1, hjust =1))
    if(gene == "IGHV") gg  <- gg +scale_fill_manual(values= c("M-CLL" = "red", "U-CLL" = "blue"))
    else if(gene == "TP53") gg  <- gg +scale_fill_manual(values= c("wt" = "gray", "mut" = "orange"))
    
    gg
}

plotCombiBoxplots <- function(dfBCmuts, drugs2plot, gene = "TP53", round2=3){
  df4plot <- dfBCmuts %>% filter(mutation == gene, combi %in% drugs2plot) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c", 5:1))) %>%
    mutate(BDrugConc = factor(round(BDrugConc, round2)))
  
  if(gene == "TP53") df4plot %<>% mutate(status = ifelse(status ==0, "wt","mut"))
  else if(gene == "IGHV") df4plot %<>% mutate(status = ifelse(status ==0, "U-CLL","M-CLL"))
  
  gg <- ggplot(df4plot, aes(x=BDrugConc, y=effectBC)) + geom_boxplot(aes(fill=factor(status))) + facet_wrap(~combi, ncol=3 , scales = "free_x") +
    ggpubr::stat_compare_means(method = "t.test", aes(group=status, label =  ..p.signif..), hide.ns = TRUE) +
    guides(fill = guide_legend(title=gene)) + xlab("Concentration (µM)") +ylab("viability") +theme_bw(base_size = 15) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, vjust=1, hjust =1))
  if(gene == "IGHV") gg  <- gg +scale_fill_manual(values= c("M-CLL" = "red", "U-CLL" = "blue"))
  else if(gene == "TP53") gg  <- gg +scale_fill_manual(values= c("wt" = "gray", "mut" = "orange"))
  
  gg
}

plotCombiResponseCurves <- function(dfBCmuts, drugs2plot, gene = "TP53", round2=3){
  df4plot <- dfBCmuts %>% filter(mutation == gene, combi %in% drugs2plot) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c", 5:1))) %>%
    mutate(BDrugConc = factor(round(BDrugConc, round2)))
  
  if(gene == "TP53") df4plot %<>% mutate(status = ifelse(status ==0, "wt","mut"))
  else if(gene == "IGHV") df4plot %<>% mutate(status = ifelse(status ==0, "U-CLL","M-CLL"))
  
  gg <- ggplot(df4plot, aes(x=BDrugConc, y=effectBC, group=status)) +
    stat_summary(fun.data = "mean_se", aes(col=factor(status)), geom="line") + facet_wrap(~combi, ncol=3, scales = "free_x") +
    stat_summary(fun.data = "mean_se",aes(col=factor(status)), geom="errorbar", width=0.2, fun.args = list(mult =2)) +
    ggpubr::stat_compare_means(method = "t.test", aes(group=status, x=BDrugConc, label =  ..p.signif..), hide.ns = TRUE) +
    guides(col = guide_legend(title=gene))+ xlab("Concentration (µM)") +ylab("viability")+theme_bw(base_size = 14) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, vjust=1, hjust =1))
  # if(length(unique(drugs2plot)) <4) gg <- gg + theme(legend.position = "top")
  if(gene == "IGHV") gg  <- gg +scale_color_manual(values= c("M-CLL" = "red", "U-CLL" = "blue"))
  else if(gene == "TP53") gg  <- gg +scale_color_manual(values= c("wt" = "black", "mut" = "orange"))
  
  gg
}

##############################
# Drug Combination Viability Comparisons
#############################

# Function to plot comparison of two combination compunds across all base drugs
plotComparsionCDrugs <- function(df, CDrugAbrv.x, CDrugAbrv.y, range = c(0,1.4), type=c("scatter_joint", "scatter_factBDrug", "boxplot_joint")){

  type <- match.arg(type)
  stopifnot(type %in% c("scatter_joint", "scatter_factBDrug", "boxplot_joint"))
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
      annotate("text", x=range[1]+0.1, y=range[2]-0.1,
               label=paste("cor ==", (round(cor(df$effectBC.x,df$effectBC.y),2))), parse=T) +
      # annotate("text", x=range[2]-0.1, y=range[2]-0.1,
      #          label=paste("R2 ==", round(1 - sum((df$effectBC.x - df$effectBC.y)^2)/sum(df$effectBC.x - mean(df$effectBC.x)^2),2),
      #                      round( 1- sum((df$effectBC.x - df$effectBC.y)^2)/sum(df$effectBC.y - mean(df$effectBC.y)^2),2), parse=T)) +
      geom_abline(slope=1, intercept=0, lty="dashed") + guides(col=guide_legend(title="Base compound"))  +
      theme_bw(base_size=16)

} else if(type == "scatter_factBDrug"){
  gg <- ggplot(df, aes(x=effectBC.x, y=effectBC.y, col=BDrugName)) +
    geom_point() +
    xlab(paste("Combination effect with",CDrugAbrv.x,  "\n (viability relative to DMSO control)")) +
    ylab(paste("Combination effect with",CDrugAbrv.y,  "\n (viability relative to DMSO control)")) +
    facet_wrap(~BDrugName, ncol=5) +
    # annotate("text", x=range[1]+0.3, y=range[2]-0.1,
    #          label=paste("cor ==", (round(summarize(group_by(df, BDrugName), cor=cor(effectBC.x, effectBC.y))$cor,2))), parse=T, size=3) +
    annotate("text", x=range[1]+0.3, y=range[2]-0.1,
             label=paste("R^2 ==", pmin(round(summarize(group_by(df, BDrugName), r2= 1 - sum((effectBC.x - effectBC.y)^2)/sum(effectBC.x - mean(effectBC.x)^2))$r2,2),
                         round(summarize(group_by(df, BDrugName), r2= 1 - sum((effectBC.x - effectBC.y)^2)/sum(effectBC.y - mean(effectBC.y)^2))$r2,2))), parse=T, size=3) +
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
                        conc4average = paste0("c",1:5), nOcc = 10, DrugDrugbyIGHV=FALSE, returnMat = FALSE,
                        dist2usecols = "euclidean", dist2userows = "euclidean"){
  
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
    #rotate tree on patients by IGHV status
    if(!useAverage){
      # outlying values cut off at 1.4
      effectBPlusC_mat[is.na(effectBPlusC_mat)] <- 1.4
      fsz <- 5
    } else fsz <- 12
    
    callbackIGHV = function(hc, mat){
      ighv <- dfMuts$IGHV
      names(ighv) <- dfMuts$PatientID
      dend = reorder(as.dendrogram(hc), wts = ighv[rownames(mat)])
      as.hclust(dend)
    }
      pheatmap(effectBPlusC_mat, na_col="gray", clustering_distance_rows=dist2userows,
               clustering_distance_cols = dist2usecols,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(140),
               breaks=seq(0,filter_th,0.01), show_rownames = FALSE, show_colnames = TRUE,
               treeheight_row = 15, treeheight_col = 15, annotation_row = dfanno, annotation_colors =anno_colors,
               annotation_legend = FALSE, fontsize_col=fsz,
               clustering_callback = callbackIGHV)

    
  } else if(type == "DrugDrug"){
      if(DrugDrugbyIGHV) {
        mat <- effectBPlusC_mat[patsMCLL, ]
        corDrug <- cor(mat, use="complete.obs")
        hmMCLL <- pheatmap(corDrug, na_col="gray", clustering_distance_rows=dist2userows,
               clustering_distance_cols = dist2usecols,
               color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(120),
               treeheight_row = 0, treeheight_col = 0, #main = "M-CLL",
               breaks=seq(-0.2,1,0.01),
               legend=FALSE, show_colnames = FALSE, show_rownames = FALSE,  fontsize = 16)
        mat <- effectBPlusC_mat[patsUCLL, ]
        corDrug <- cor(mat, use="complete.obs")
        hmUCLL <- pheatmap(corDrug, na_col="gray", clustering_distance_rows=dist2userows,
                 clustering_distance_cols = dist2usecols,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(120),
                 treeheight_row = 0, treeheight_col = 0, #main = "U-CLL",
                 breaks=seq(-0.2,1,0.01), show_colnames = FALSE, show_rownames = FALSE,
                 fontsize = 16, legend=FALSE)
        grid.arrange(hmMCLL$gtable, hmUCLL$gtable, ncol =2)
        
      } else {
        corDrug <- cor(effectBPlusC_mat, use="complete.obs")
        pheatmap(corDrug, na_col="gray", clustering_distance_rows=dist2userows,
                 clustering_distance_cols = dist2usecols,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(140),
                 treeheight_row = 15, treeheight_col = 15, fontsize = 12)
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

# plot the single drug responses and the response to the combination as well as the independent effect as curves (mean + SE) for a pairs of drugs drB and drC.
plotResponseCurves <- function(df, drC , drB, th = filter_th, CItype = "SI",
                               sep_by_IGHV =FALSE, sep_by_TP53=FALSE, annoSI=FALSE, annoP = TRUE){

  stopifnot(CItype %in% c("SI", "HSA"))
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName == drB) %>% 
    select(CDrugAbrv, BDrugName, BDrugConcId,BDrugConc, starts_with("viab"), PatientID) %>%
    gather(key="type", value = "viability", starts_with("viab")) %>%
    mutate(type = ifelse(type == "viabB", paste(drB, "(A)"),
                         ifelse(type == "viabC", paste(drC, "(B)"),
                                ifelse(type == "viabBC", paste(drB, "+", drC, "(AB)"),
                                       "Expected effect (A*B)")))) %>%
    mutate(type = factor(type, levels = c(paste(drB, "(A)"), paste(drC, "(B)"),  paste(drB, "+", drC, "(AB)"), "Expected effect (A*B)")))
  
  df4plot <- left_join(df4plot, dfMuts4testing, by ="PatientID")
  df4plot %<>% mutate(IGHV = ifelse(IGHV == 0, "U-CLL", "M-CLL"))
  df4plot %<>% mutate(TP53 = ifelse(TP53 == 0, "TP53-wt", "TP53-mut"))
  
  # annotate by SI
  dfanno <- dfsynSummaryPat %>% filter(CDrugAbrv == drC, BDrugName == drB) %>%
    select(BDrugConcId, addModelSImed, hsaCImed)
  df4plot %<>% left_join(dfanno, by = "BDrugConcId")
  
  # annotate by p-values
  if(CItype == "SI") dfsig <- dfsigSI else dfsig <- dfsigHSA
  dfannop <- dfsig %>% filter(CDrugAbrv == drC, BDrugName == drB) %>%
    select(BDrugConcId, pval)
  df4plot %<>% left_join(dfannop, by = "BDrugConcId")
  
  df4plot %<>% mutate(BDrugConc = factor(round(BDrugConc*1000,1)))
  
  gg <- ggplot(data=df4plot, aes(x=BDrugConc, y=viability, col = type, group=type, linetype = type)) +
    stat_summary(fun.data = "mean_se", fun.args = list(mult = 2), geom="errorbar", width=0.05) +
    stat_summary(fun.y = "mean", geom="line", fun.args = list(mult = 2)) + 
    theme_bw(base_size = 20) + xlab(paste0("Concentration of ", drB, " (nM)")) +
    theme(legend.position = "top", legend.title = element_blank()) + 
    guides(col = guide_legend(ncol=1), linetype = guide_legend(ncol=1)) +
    #guides(col = guide_legend(title="", ncol=1)) +
    ylim(c(0,th)) +
    scale_color_manual(values =  brewer.pal(name ="Set1",9)[c(2,3,5,4)]) + 
    scale_linetype_manual(values = c(rep("solid", 3), "dashed"))

  if(sep_by_IGHV & !sep_by_TP53){
    gg <- gg + facet_wrap(~IGHV)
  } else if(!sep_by_IGHV & sep_by_TP53){
    gg <- gg + facet_wrap(~TP53)
  } else if(sep_by_IGHV & sep_by_TP53){
    gg <- gg + facet_wrap(TP53~IGHV)
  }
  
  if(annoSI){
    if(sep_by_IGHV | sep_by_IGHV) stop("SI annotation for separate IGHV or TP53 not implemented yet")
    if(CItype == "SI") {
      gg <- gg + geom_text(aes(x = BDrugConc, label = round(addModelSImed,2)), y=1.3, col="black", size=5)
    } else{
      gg <- gg + geom_text(aes(x = BDrugConc, label = round(hsaCImed,2)), y=1.3, col="black", size=5)
    }
    
  }
  
  if(annoP){
  if(sep_by_IGHV | sep_by_IGHV) stop("SI p-value annotation for separate IGHV or TP53 not implemented yet")
  gg <- gg + geom_text(aes(x = BDrugConc, label = paste0("p=",format(pval, digits=2))), y=0, col="gray", size=4)
  }
  
  return(gg)
}

# plot the single drug responses and the response to the combination as well as the independent effect as curves (mean + SE) for a pairs of drugs drB and drC.
plotMultipleBResponseCurves <- function(df, drC , drsB, th = filter_th, annoSI = FALSE){
  df4plot <- df %>% filter(CDrugAbrv == drC, BDrugName %in% drsB) %>% 
    mutate(BDrugConcId = factor(BDrugConcId, levels = paste0("c",5:1))) %>%
    select(CDrugAbrv, BDrugName, BDrugConcId,BDrugConc, starts_with("viab"), PatientID) %>%
    gather(key="type", value = "viability", starts_with("viab")) %>%
    mutate(type = ifelse(type == "viabB", paste("Base compound", "(A)"),
                         ifelse(type == "viabC", paste(drC, "(B)"),
                                ifelse(type == "viabBC", paste("Base compound", "+", drC, "(AB)"),
                                       "Expected effect (A*B)")))) %>%
    mutate(type = factor(type, levels = c(paste("Base compound", "(A)"), paste(drC, "(B)"),  paste("Base compound", "+", drC, "(AB)"), "Expected effect (A*B)")))
  
  df4plot <- left_join(df4plot, dfMuts4testing, by ="PatientID")
  df4plot %<>% mutate(IGHV = ifelse(IGHV == 0, "U-CLL", "M-CLL"))
  
  # annotate by SI
  dfanno <- dfsynSummaryPat %>% filter(CDrugAbrv == drC, BDrugName %in% drsB) %>%
    select(BDrugConcId, addModelSImed, BDrugName)
  df4plot %<>% left_join(dfanno, by = c("BDrugConcId", "BDrugName"))
  df4plot %<>% mutate(BDrugName = factor(BDrugName, levels =drsB)) #keep order as specified
  df4plot %<>% mutate(BDrugConc = factor(round(BDrugConc*1000,1)))
  
  gg <- ggplot(data=df4plot, aes(x=BDrugConc, y=viability, col =type, group=type, linetype =type)) +
    stat_summary(fun.data = "mean_se", fun.args = list(mult = 2), geom="errorbar", width=0.15) +
    stat_summary(fun.y = "mean", geom="line") + 
    theme_bw(base_size = 20) + xlab(paste0("Concentration (nM)")) +
    guides(col = guide_legend(title="")) + ylim(c(0,th)) +
    facet_wrap(~ BDrugName, ncol = 4) +theme(legend.position = "top", legend.title = element_blank(), legend.spacing.x = unit(0.4, 'cm')) +
    theme(strip.background = element_blank()) +
    guides(col = guide_legend(ncol=2), linetype = guide_legend(ncol=2)) +
    #guides(col = guide_legend(title="", ncol=1)) +
    scale_color_manual(values =  brewer.pal(name ="Set1",9)[c(2,3,5,4)]) + 
    scale_linetype_manual(values = c(rep("solid", 3), "dashed"))
  
  if(annoSI){
    gg <- gg + geom_text(aes(x = BDrugConc, label = round(addModelSImed,2)), y=1.3, col="black", size=5)
  }
  
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
  
    df4plot %<>% mutate(BDrugConc = factor(round(BDrugConc*1000,1)))
  
    ggplot(df4plot, aes(x=BDrugConc, y=CI)) +
    ggbeeswarm::geom_beeswarm(col = "gray", alpha =0.7) + 
    geom_boxplot(outlier.shape = NA, alpha=0.2, width = 0.3) +
    geom_hline(yintercept = ifelse(CI_type =="Bliss",1,0), lty = "dashed") +
    theme_bw(base_size = 15)  + xlab(paste0("Concentration of ", drB, " (nM)")) + ylab(paste0("Combination index (", CI_type, ")"))
}

# plot bar plot of individual combination indices for each sample
plotWaterfallCI <-  function(df, drC , drB, CI_type = c("Bliss", "hsa", "SI"),
                             annotate=NULL, pats2label = NULL, y_nudge = 0.1, label_size=5){
  
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
    # geom_errorbar(aes(ymin = CI - CIse, ymax = CI + CIse), width = 0.3)+
    geom_hline(yintercept = ifelse(CI_type =="SI",0,1), lty = "dashed") +
    theme_bw(base_size = 15)  + xlab(paste0("Patient sample")) +
    ylab(paste0("Combination index (", CI_type, ")")) + 
    scale_fill_manual(values = patcol) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(fill=FALSE)
  
  # annotate by genetic features
  if(!is.null(annotate)){
    gg <- gg + ylim(c(min(df4plot$CI), max(df4plot$CI) +0.04))
    if("IGHV" %in% annotate){
    gg <- gg + geom_text(aes(y = max(df4plot$CI) +0.01, label = ifelse(IGHV==1, "*", "")),col="red", size=5)
    gg <- gg + geom_text(x=nrow(df4plot)-10, y= min(df4plot$CI), label = "* M-CLL", col="red", size=5)
    }
    if("TP53" %in% annotate){
    gg <- gg + geom_text(aes(y = max(df4plot$CI) +0.03, label = ifelse(TP53==1, "*", "")), col="black", size=5)
    gg <- gg + geom_text(x=nrow(df4plot)-3, y= min(df4plot$CI), label = "* TP53 mut", col="black", size=5)
    }
  }
  
  gg  <- gg + geom_text(y= max(df4plot$CI)- y_nudge,
                        aes(label = ifelse(PatientID %in% pats2label, as.character(PatientID), "")),
                        angle = 90, size = label_size)
  return(gg)
}

# dotplots add therotical effect vs measured combination effect
plotScattter <- function(df, drB, drC, th = filter_th){
  
  df4plot <- filter(df, BDrugName==drB)
  range = c(0, th)
  df4plot %<>% mutate(label = factor(paste0(round(BDrugConc * 1000,1), " (nM)"),
                                     levels = paste0(sort(unique(round(BDrugConc * 1000,1))), " (nM)")))
  gg <- ggplot(df4plot, aes(x=viabBC_add, y=viabBC, color=PatientID))+
    geom_point(alpha=0.7)+
    geom_hline(aes(yintercept=1), colour="grey", linetype="dashed") +
    geom_vline(aes(xintercept=1), colour="grey", linetype="dashed") +
    scale_color_manual(values=patcol)+facet_wrap(~ label  , ncol=5) +
    geom_abline(intercept = 0, slope = 1, colour="black", linetype="solid") +
    # ggtitle(paste("Base drug (A):", drB, "\n Combination drug (B):", drC)) +
    coord_fixed() + scale_x_continuous(limits=range) + 
    scale_y_continuous( limits=range) +
    theme_bw(base_size = 12) +
    ylab(paste0("Measured viability")) +
    xlab(paste0("Expected viability of ", drB," with ", drC)) +
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

plotBoxplotSIMuts <- function(drB, drC, gene = "IGHV") {
  gg <- left_join(dfsyn,dfMuts4testing, by = "PatientID") %>%
    filter(BDrugName == drB, CDrugAbrv == drC) %>%
    mutate(IGHV = ifelse(IGHV == 0, "U-CLL", "M-CLL")) %>%
    mutate(TP53 = ifelse(TP53 == 0, "wt", "mut")) %>%
    mutate(BDrugConc = factor(round(BDrugConc*1000,0))) %>%
    ggplot(aes_string(x="BDrugConc", y="addModelSI", fill= gene)) +geom_boxplot(width=0.4) + #facet_wrap(~BDrugConc, nrow=1) + 
    # ggtitle(paste(drC, drB, sep="+")) +
    ylab("SI") + 
    ggpubr::stat_compare_means(aes(label =  paste0("p=",..p.format..)), method = "t.test") +
    theme_bw(base_size = 14) +
    theme(strip.background = element_blank(),
          legend.position = "top",        
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,-7,0)) + # move legend closer to panel 
    xlab("Concentration (nM)")
  
  if(gene == "IGHV") gg  <- gg + scale_fill_manual(values= c("M-CLL" = "red", "U-CLL" = "blue"))
  else if(gene == "TP53") gg  <- gg +scale_fill_manual(values= c("wt" = "gray", "mut" = "orange"))
  gg
}

#############################
# 10x10 screens
#############################

plotTiles10x10 <- function(df, drB, pat, type = c("tile", "contour")){
    type <- match.arg(type)
    df4plot <- filter(df, PatientID==pat, BaseDrugName == drB, CombiDrug == "Ibrutinib") 
    df4plot %<>% rename(viability = normalizedValue)
    df4plot %<>%  select(concCvalue, concBvalue, viability)
    if(type == "tile"){
    gg <- ggplot(df4plot, aes(x=factor(round(concBvalue*1000,1)), y=factor(round(concCvalue*1000,1)), fill = viability)) +
        geom_tile() + 
        ggtitle(pat) +
        scale_fill_gradient(low = "white",high = "navy", limits=c(0,1.4)) +
      ylab("Ibrutinib (nM)") + xlab(paste(drB, "(nM)")) +
      theme_bw(base_size = 16) + coord_fixed() +
      theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1))
    gg
    } else if(type == "contour"){
    gg <- ggplot(df4plot, aes(x=(concBvalue), y=(concCvalue), z=viability)) + geom_contour(aes())
    gg
    }
    return(gg)
    }  

plotCITiles <- function(df, CItype, cutoff = 100){
  dfres <- data.frame() # replace by lapply, bind_rows
  for(dr in unique(df$BaseDrugName)){
    if(dr!="DM"){
      print(dr)
      data <- filter(df, BaseDrugName == dr, CombiDrug == "Ibrutinib")
      for(pat in unique(data$PatientID)){
        print(pat)
        dfpat <- filter(data, PatientID == pat) %>%
                mutate(base_conc = factor(base_conc, levels = paste0("c",1:10)),
                       combi_conc = factor(combi_conc, levels = paste0("c",1:10))) %>%
          select(Replicate = PatientID, DrugRow = BaseDrugName, DrugCol = CombiDrug,
                 ConcRow = concBvalue, ConcCol = concCvalue, Response = normalizedValue,
                 Row = base_conc, Col = combi_conc) %>%
          mutate(BlockID = 1, ConcRowUnit = "μM", ConcColUnit = "μM")   %>% # only one drug-drug combination
          mutate(Row = sub("c", "", Row),
                 Col = sub("c", "", Col),
                 Response = pmin(cutoff,100 * Response)) # need percentage
        dfpat <- ReshapeData(dfpat, data.type = "viability") # does not work with multiple replicates
        
        #re-order by concentrations (in ReshapeData ordered as charatcers 1, 10,2,...)
        dfpat$dose.response.mats[[1]] <- dfpat$dose.response.mats[[1]][order(as.numeric(rownames(dfpat$dose.response.mats[[1]]))), order(as.numeric(colnames(dfpat$dose.response.mats[[1]])))]
        # print(PlotDoseResponse(dfpat))
        
        if(CItype != "myLoewe"){
        synergy.score <- CalculateSynergy(dfpat,method = CItype, correction = TRUE)
        df_score <- melt(synergy.score$scores[[1]], varnames = c("concB", "concC"), value.name = "score") %>%
          filter(concB !=0 & concC!=0)
        } else {
        synergy.score <- myLoewe(dfpat$dose.response.mats[[1]])
        df_score <- melt(synergy.score$scores, varnames = c("concB", "concC"), value.name = "score") %>%
          filter(concB !=0 & concC!=0)
        }
        # dfsynCon <- filter(df_score, score>0) %>% select(concB,concC)
        # print(paste(dr,":", paste0(unique(dfsynCon$concB), collapse = ", ")))
        # print(paste("Ibrutnib :", paste0(unique(dfsynCon$concC), collapse = ", ")))
        
        gg <- ggplot(df_score, aes(x = factor(round(concB*1000,1)), y=factor(round(concC*1000,1)), fill = score)) +
          geom_tile() +
          ylab("Ibrutinib (nM)") + xlab(paste(dr, "(nM)")) +
          # ggtitle(paste(pat, ": average score (", CItype, ") ", round(mean(df_score$score),3), sep="")) +
          ggtitle(pat) +
          scale_fill_gradient2(low= "blue", high="red", mid="white", midpoint = 0) +
          theme_bw(base_size = 16) + coord_fixed() +
          theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1),
                plot.title = element_text(colour =  "black"))
        print(gg)
        # PlotSynergy(synergy.score, type = "2D", save.file = TRUE)
        dfres <- rbind(dfres,cbind(df_score, PatientID = pat, BDrugName = dr, CDrugName = "Ibrutinib"))
      }
    }
  }
  return(dfres)
}


plotSummaryLoewe <- function(dfLoewe, dr, summarize_by = "mean", type = "col") {
  stopifnot(type %in% c("row", "col"))
  
  dfLoewe_dr <- filter(dfLoewe, BDrugName == dr) 
  
  orderConc <- as.numeric(factor(rank(dfLoewe_dr$concB)))
  dfLoewe_dr$concBrank <- as.numeric(orderConc)
  orderConcC <- as.numeric(factor(rank(dfLoewe_dr$concC)))
  dfLoewe_dr$concCrank <- as.numeric(orderConcC)
  # smooth along B and C concentrations
  dfLoewe_dr$sumScore <- sapply(1:nrow(dfLoewe_dr), function(i){
    df <- filter(dfLoewe_dr, abs(concBrank-dfLoewe_dr$concBrank[i])<= 1,
                 abs(concC - dfLoewe_dr$concC[i])<= 1) # take mean across two neighboring distr in both direction
    mean(df$score, na.rm=TRUE)
  })
    
  
  # if(length(unique(dfLoewe_dr$concB)) < 15){
  #   dfLoewe_dr %<>% group_by(concB, concC) %>%
  #     summarize(sumScore = switch(summarize_by,
  #                                mean  = mean(score),
  #                                median  = median(score)))
  # } else {
  #   orderConc <- as.numeric(factor(rank(dfLoewe_dr$concB)))
  #   dfLoewe_dr$concBrank <- as.numeric(orderConc)
  #   orderConcC <- as.numeric(factor(rank(dfLoewe_dr$concC)))
  #   dfLoewe_dr$concCrank <- as.numeric(orderConcC)
  #   # smooth along B concentration to reconcile different conc from pats
  #   dfLoewe_dr$sumScore <- sapply(1:nrow(dfLoewe_dr), function(i){
  #     df <- filter(dfLoewe_dr, abs(concBrank-dfLoewe_dr$concBrank[i])<= 2,
  #                  concC == dfLoewe_dr$concC[i]) # take mean across two neighboring distr in both direction
  #     mean(df$score)
  #   })
  # }
  
  gg_medScore <- ggplot(dfLoewe_dr, aes(x = factor(round(concB*1000,1)), y=factor(round(concC*1000,1)), fill = sumScore)) +
    geom_tile() +
    ylab("Ibrutinib (nM)") + xlab(paste(dr, "(nM)")) +
    scale_fill_gradient2(low= "blue", high="red", mid="white", midpoint = 0, breaks = seq(-100,100,20)) +
    theme_bw(base_size = 16) + coord_fixed() +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1),
          plot.title = element_text(colour =  "black")) +
    guides(fill = guide_colorbar(title = paste0("Loewe \n score")))
  
  if(type == "row"){
    gg_medScore <- gg_medScore + theme(axis.text.y = element_blank(),
                                       axis.title.y = element_blank(),
                                       axis.ticks.y = element_blank(),
                                       legend.position = "top",
                                       plot.margin = unit(c(0,-1,0,0.5), "cm"),
                                       legend.margin=margin(0,0,0,0),
                                       legend.box.margin=margin(0,0,-7,0))  + 
      guides(fill = guide_colorbar(title.position="left", title.hjust = 1, title.vjust = 1,
                                   title = paste0("Loewe  score")))
  }
  

  df10x10_dr <- filter(df10x10, BaseDrugName == dr, CombiDrug== "Ibrutinib")
  
  orderConc <- as.numeric(factor(rank(df10x10_dr$concBvalue)))
  df10x10_dr$concBrank <- as.numeric(orderConc)
  orderConcC <- as.numeric(factor(rank(df10x10_dr$concCvalue)))
  df10x10_dr$concCrank <- as.numeric(orderConcC)
  # smooth along B concentration to reconcile different conc from pats
  df10x10_dr$sumViab <- sapply(1:nrow(df10x10_dr), function(i){
    df <- filter(df10x10_dr, abs(concBrank-df10x10_dr$concBrank[i])<= 1,
                 abs(concCvalue - df10x10_dr$concCvalue[i]) <=1) # take mean across two neighboring distr
    mean(df$normalizedValue)
  })
    
  # if(length(unique(df10x10_dr$concBvalue)) < 15){
  #   df10x10_dr %<>% group_by(concBvalue, concCvalue) %>%
  #   summarize(sumViab = switch(summarize_by,
  #                              mean  = mean(normalizedValue),
  #                              median  = median(normalizedValue)))
  # } else {
  # orderConc <- as.numeric(factor(rank(df10x10_dr$concBvalue)))
  # df10x10_dr$concBrank <- as.numeric(orderConc)
  # orderConcC <- as.numeric(factor(rank(df10x10_dr$concCvalue)))
  # df10x10_dr$concCrank <- as.numeric(orderConcC)
  # # smooth along B concentration to reconcile different conc from pats
  # df10x10_dr$sumViab <- sapply(1:nrow(df10x10_dr), function(i){
  #   df <- filter(df10x10_dr, abs(concBrank-df10x10_dr$concBrank[i])<= 2,
  #                concCvalue == df10x10_dr$concCvalue[i]) # take mean across two neighboring distr
  #   mean(df$normalizedValue)
  # })
# }
  
  gg_medViab <- ggplot(df10x10_dr, aes(x = factor(round(concBvalue*1000,1)), y=factor(round(concCvalue*1000,1)), fill = sumViab)) +
    geom_tile() +
    ylab("Ibrutinib (nM)") + xlab(paste(dr, "(nM)")) +
    scale_fill_gradient(low = "white",high = "navy", limits=c(0,1.4)) +
    theme_bw(base_size = 16) + coord_fixed() +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1),
          plot.title = element_text(colour =  "black"))  + 
    guides(fill = guide_colorbar(title = paste0("viability")))
  
  if(type == "row"){
    gg_medViab <- gg_medViab + theme(plot.margin = unit(c(0,0.5,0,-1), "cm"),
                                     legend.position = "top",
                                     legend.margin=margin(0,0,0,0),
                                     legend.box.margin=margin(0,0,-7,0))+
      guides(fill = guide_colorbar(title.position="left", title.hjust = 1, title.vjust = 1, title = paste0("viability")))
  }
  cowplot::plot_grid(gg_medViab, gg_medScore, nrow=ifelse(type == "row",1,2), align = "hv", axis ="lb")
}
