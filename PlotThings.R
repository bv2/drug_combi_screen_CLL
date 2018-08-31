get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#plotting function
plotThings<-function(what=c("ABC", "DiffPerConc", "Curves","CombiBenefit"), IGHV=NULL, CombiDrugs,
                     BaseDrugs, df4anaAvreplicates, outputdir){
  what=match.arg(what)
  filename=paste(what,".pdf", sep="")
  if(!is.null(IGHV)) filename<-paste(what,"_IGHVann.pdf", sep="")
  height<-50
  width<-30
  if(what=="DiffPerConc") {height<-20; width=15} #200 
  if(what=="CombiBenefit") height<-30
  pdf(file.path(outputdir, filename), width=width, height=height)
  plot_list<- list()
  CombiBenefit<-list()
  i=1
  par(mfrow=c(16,2))
  
  for(drC in CombiDrugs) {
    for(drB in BaseDrugs) {
      
      idx = which(df4anaAvreplicates$BDrugName ==drB & df4anaAvreplicates$CDrugAbrv==drC)
      subDF<-df4anaAvreplicates[idx,]
      subDF$PatientID<-as.character(subDF$PatientID)
      
      if(what=="Curves"){
        gg<-ggplot(data=subDF) +
          geom_line(aes(x=BDrugConc, y=effectBC, group = PatientID, colour = PatientID,alpha=0.1)) +
          geom_point(aes(x=0, y=effectC, group = PatientID ,colour = PatientID)) +
          xlab(paste("conc of", drB))+ylab("viability")+ggtitle(paste(drB, drC, sep="+"))+
          geom_line(aes(x=BDrugConc, y=effectB*effectC, group = PatientID, colour = PatientID, alpha=0.1),linetype = "dashed")+
          scale_color_manual(values=patcol)+theme_bw()
        
        
        #add median curve for U-CLL and M-CLL seperately
        if(!is.null(IGHV)) {
          #remove other curves again to make plot less complex
          subDF$IGHV<-IGHV[subDF$PatientID,]
          dfUCLL<-filter(subDF, IGHV==0)
          dfMCLL<-filter(subDF, IGHV==1)
          DFbyConc<-aggregate(select(dfUCLL, effectB, effectC, effectBC), by=list(conc=dfUCLL$BDrugConc), FUN=median)
          gg<-ggplot(data=DFbyConc,aes(x=conc, y=effectBC))+geom_line(colour="darkgreen", linetype="solid", size=1)+
            geom_line(data=DFbyConc,aes(x=conc, y=effectB*effectC), colour="darkgreen", linetype="dashed", size=1)+
            xlab(paste("conc of", drB))+ylab("viability")+ggtitle(paste(drB, drC, sep="+"))+theme_bw()
          DFbyConc<-aggregate(select(dfMCLL, effectB, effectC,effectBC), by=list(conc=dfMCLL$BDrugConc), FUN=median)
          gg<-gg+geom_line(data=DFbyConc,aes(x=conc, y=effectB*effectC), colour="red", linetype="dashed", size=1)+
            geom_line(data=DFbyConc,aes(x=conc, y=effectBC), colour="red", linetype="solid", size=1)
        }
        
        #add median curve
        DFbyConc<-aggregate(select(subDF, effectB, effectC, effectBC), by=list(conc=subDF$BDrugConc), FUN=median)
        gg<-gg+geom_line(data=DFbyConc,aes(x=conc, y=effectBC), colour="black", linetype="solid", size=1.5)+
          geom_line(data=DFbyConc,aes(x=conc, y=effectB*effectC), colour="black", linetype="dashed", size=1.5)
        
        
                        
        if(i==1&is.null(IGHV)) legend <- get_legend(gg)
        gg<-gg+  theme(legend.position="none")
        plot_list[[drC]][[drB]]<-gg
        i<-2
      }
      
      if(what=="ABC"|what=="CombiBenefit"){
        #calculate area between curves
        by_AB<--by(subDF, subDF$PatientID, function(x) trapz(x$BDrugConc, x$effectBC), simplify = TRUE)
        AUC_AB<-as.vector(by_AB); names(AUC_AB)<-names(by_AB)
        
        by_addAB<--by(subDF, subDF$PatientID, function(x) trapz(x$BDrugConc, x$effectB*x$effectC), simplify = TRUE)
        AUC_addAB<-as.vector(by_addAB); names(AUC_addAB)<-names(by_addAB)
        
        ABC<- - (-AUC_addAB+AUC_AB)/as.numeric(filter(DrugBase,drB==Substance)[4])
        orderABC<-order(ABC)
        ABC_ord<-ABC[orderABC]
        if(what=="ABC"){    
          bplt<-barplot(ABC_ord, las=2, col=patcol[substr(names(ABC_ord),1,5)], 
                        main = paste("lacking area under viability curves of", drC, "plus", 
                                     drB, "compared to an additive effect model"),
                        cex.main=1.5, cex.axis = 1, names.arg = substr(names(ABC_ord),1,5),
                        ylab= "<-- less than additive  more than additive -->") 
          if(!is.null(IGHV)) text(x=bplt[IGHV[names(ABC_ord),]==1], y=median(ABC), labels="*", cex=5)
          abline(h=median(ABC), col="red")
        } else {
          #normalize by highest concentration used
          CombiBenefit[[drC]][[drB]]<-median(ABC)
        }
      }
      
      
      if(what=="DiffPerConc"){
        if(!is.null(IGHV)) {
          subDF$IGHV<-IGHV[subDF$PatientID,]
          subDF$IGHV[ subDF$IGHV!=1]<-""
          subDF$IGHV[ subDF$IGHV==1]<-"*"
        }
        #calculate differences at each concentration
        subDF$diff<- -(-subDF$effectC*subDF$effectB+subDF$effectBC)
        CombiBenefit[[drC]][[drB]]<-mean(by(subDF$diff, subDF$PatientID, median))

        #if average over patients wanted, comment print(gg) below
        # par(mfrow=c(1,1))
        # barplot((by(subDF$diff, subDF$BDrugConcId, mean)), las=2, main=paste(drB, drC))
        
        gg<-ggplot(subDF, aes(x=paste(PatientID, BDrugConcId), y=diff, fill=PatientID, color="black") )+
          geom_bar(position = "dodge", stat="identity")+
          scale_fill_manual(values=patcol)+scale_color_manual(values="black")+coord_flip()+
          ggtitle(paste("Lacking viability after combination treatment of", drC, "and", drB, "compared to an additive model")) +
          ylab("<-- less than additive  more than additive -->")
        
        if(!is.null(IGHV)) gg<-gg +geom_text( aes(label=IGHV), size=10)
        print(gg)
        
        #     subDFgr<-arrange(group_by(subDF, PatientID))
        #     plot_list[[drC]][[drB]]<-barchart(diff~PatientID,data=subDFgr,groups=BDrugConcId, origin=0, 
        #                main = paste("Differnce of additive effect and combination effect of", drC, "and", drB),
        #                auto.key=list(columns=5), scales=list(x=list(rot=90,cex=0.8)),
        #                col)
      }
    }}
  
  #make plots
  if(what=="Curves"){
    for(i in 1:length(plot_list)) {
      do.call(grid.arrange, c(plot_list[[i]], ncol=4))
    }
    if(is.null(IGHV)) grid.arrange(legend)
  }
  
  if(what=="CombiBenefit"){
    par(mfrow=c(7,2), mar=c(10,4,4,2))
    for(i in 1:length(CombiBenefit)) {
      barplot(CombiBenefit[[i]][order(CombiBenefit[[i]])], las=2, 
              main=paste("Additional benefit from combination with", 
                         names(CombiBenefit)[i]),
              ylab= "<-- less than additive  more than additive -->")   
    }
  }
  
  #   if(what=="DiffPerConc"){
  #     for(i in 1:length(plot_list)) do.call(grid.arrange, c(plot_list[[i]], ncol=1))
  #   }
  
  dev.off()
  if(what=="CombiBenefit"|what=="DiffPerConc") return(CombiBenefit)
}
