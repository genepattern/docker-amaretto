AMARETTO_VisualizeModule <- 
  function(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,ModuleNr) {
    common_samples <- Reduce(intersect, list(colnames(CNV_matrix),colnames(MET_matrix))); 
    CNV_matrix <- CNV_matrix[,common_samples]; MET_matrix <- MET_matrix[,common_samples]
    # getting the data
    if (ModuleNr>AMARETTOresults$NrModules){
      cat('\tCannot plot Module',ModuleNr,'since the total number of modules is',AMARETTOresults$N,'.\n')
    }
    else {
      ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
      currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
      ModuleGenes=rownames(ModuleData)
      cat('Module',ModuleNr,'has',length(rownames(ModuleData)),'genes and',length(currentRegulators),'regulators for',length(colnames(ModuleData)),'samples.\n')
      
      # Clustering the module itself
      SampleClustering=hclust(dist(t(ModuleData)), method = "complete", members = NULL)    
      GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
      ClustRegulatorData <- RegulatorData[,SampleClustering$order]
      ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]
      ClustCombinedData <- rbind(ClustModuleData,ClustRegulatorData)
      
      # create annotations 
      Alterations <- rep(0,nrow(ClustCombinedData))
      Alterations[(nrow(ModuleData)+1):nrow(ClustCombinedData)] <- rowSums(cbind(10*AMARETTOinit$RegulatorAlterations$Summary[currentRegulators,1],AMARETTOinit$RegulatorAlterations$Summary[currentRegulators,2]))
      if (length(which(Alterations==11))>0){
        if (length(which(Alterations==1))>0){
          if (length(which(Alterations==10))>0){
            case = 1 # everything
          } else {
            case = 2 # both and MET
          }
        } else {
          if (length(which(Alterations==10))>0){
            case = 3 #both and CNV
          } else {
            case =4 #Both only
          }
        }
      } else {
        if (length(which(Alterations==1))>0){
          if (length(which(Alterations==10))>0){
            case = 5 # CNV and MET
          } else {
            case = 6 # MET only
          }
        } else {
          if (length(which(Alterations==10))>0){
            case = 7 # CNV only
          }
        }
      }
      Alterations[which(Alterations==1)] <- rep("Methylation aberrations",length(which(Alterations==1)))
      Alterations[which(Alterations=="10")] <- rep("Copy number alterations",length(which(Alterations=="10")))
      Alterations[which(Alterations=="11")] <- rep("Both methylation and copy number alterations",length(which(Alterations=="11")))
      Alterations[which(Alterations=="0")] <- rep(" ",length(which(Alterations=="0")))
      Alterations <- data.frame(Alterations)
      
      if (case==1){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue"," "= "white","Both methylation and copy number alterations"="bisque4")),
                               which = "row", width = unit(1, "cm"),name="")
      } 
      if (case==2){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Both methylation and copy number alterations"="bisque4"," "="white","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      if (case==3){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Both methylation and copy number alterations"="bisque4"," "="white","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      if (case==4){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Both methylation and copy number alterations"="bisque4"," "="white","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      if (case==5){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      if (case==6){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue"," "="white")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      if (case==7){
        ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white")),
                               which = "row", width = unit(1, "cm"),name="")
      }
      CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV),currentRegulators)
      METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET),currentRegulators)
      if (length(CNVreg)>0){
        CNVData <- matrix(0,nrow=length(CNVreg),ncol=ncol(ModuleData))
        colnames(CNVData) <- colnames(ModuleData)[SampleClustering$order]
        rownames(CNVData) <- CNVreg
        CNVData[,colnames(CNV_matrix)] <- CNV_matrix[rownames(CNVData),]  
        CNVData[which(CNVData>0)] <- "Amplified"  # amplification
        CNVData[which(CNVData<0)] <- "Deleted"  # deletion
        CNVData[which(CNVData==0)] <- " " # nothing      
        
        if (length(METreg)>0){  
          METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
          colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
          rownames(METData) <- METreg
          METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
          METData[which(METData>0)] <- "Hyper-methylated"  # hyper
          METData[which(METData<0)] <- "Hypo-methylated"  # hypo
          METData[which(METData==0)] <- " " # nothing      
          Genes <- data.frame(t(METData),No=rep(" ",ncol(ModuleData)),t(CNVData))
        } else {
          Genes <- data.frame(No=rep(" ",ncol(ModuleData)),t(CNVData))
        }
      } else {
        METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
        colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
        rownames(METData) <- METreg
        METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
        METData[which(METData>0)] <- "Hyper-methylated"  # hyper
        METData[which(METData<0)] <- "Hypo-methylated"  # hypo
        METData[which(METData==0)] <- " " # nothing      
        Genes <- data.frame(t(METData),No=rep(" ",ncol(ModuleData)))
      }
      ColAnnotation <- c("Hyper-methylated" = "yellow", "Hypo-methylated" = "cornflowerblue"," "="white","Amplified"="sienna1","Deleted"="cyan")
      ColAnnotation <- rep(list(ColAnnotation),ncol(Genes))
      names(ColAnnotation) <- colnames(Genes)
      haRow = HeatmapAnnotation(df=Genes ,name="test",
                                col = ColAnnotation,which="column",show_legend=FALSE)
      
      # plotting
      heatmap <- Heatmap(ClustCombinedData, name = "Gene expression", column_title = paste('Module',ModuleNr), cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ModuleData)),rep("black",nrow(RegulatorData))),fontsize=10),
                         column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Module Genes",nrow(ModuleData)),rep(" Regulators",nrow(RegulatorData))),gap = unit(5, "mm"),
                         #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                         col=colorRamp2(c(-max(abs(ClustCombinedData)), 0, max(abs(ClustCombinedData))), c("green", "black", "red")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=haRow)
      ComplexHeatmap::draw(heatmap+ha)
      
      #    nf <- layout(matrix(c(1,2),2,1, byrow=T), widths=c(6),heights=c(2,3), respect=T)
      #    layout.show(nf)
      #    image(1:length(colnames(ModuleData)),1:length(currentRegulators),t(RegulatorData[,SampleClustering$order]),col=rev(brewer.pal(11,"RdBu")),xlab='',ylab='',main=paste('Module',ModuleNr),axes=FALSE)
      #    axis(LEFT<-2, at=1:length(currentRegulators), labels=currentRegulators, las= 2,cex.axis=0.7)        
      #    image(1:length(colnames(ModuleData)),1:length(ModuleGenes),t(ModuleData[GeneClustering$order,SampleClustering$order]),col=rev(brewer.pal(11,"RdBu")),xlab='Samples',ylab='ModuleGenes',main='',axes=FALSE)
      #axis(LEFT<-2,at=1:length(ModuleGenes), labels=ModuleGenes, las= 2,cex.axis=0.7)    
      if (length(CNVreg)>0){
        if (length(METreg)>0){
          MeanMET <- floor(nrow(METData)/2)+1
          MeanCNV <- floor(nrow(CNVData)/2)+1+nrow(METData)+1
          for(an in colnames(Genes)) {
            decorate_annotation(an, {
              if (an=="No"){
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
              } else {
                # annotation names on the right
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
              }
              if (an==colnames(Genes)[MeanMET]){
                grid.text("MET", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
              }
              if (an==colnames(Genes)[MeanCNV]){
                grid.text("CNV", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
              }
              # annotation names on the left
              # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
            })
          }
        } else {
          MeanCNV <- floor(nrow(CNVData)/2)+1+1
          for(an in colnames(Genes)) {
            decorate_annotation(an, {
              if (an=="No"){
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
              } else { 
                # annotation names on the right
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
              }
              if (an==colnames(Genes)[MeanCNV]){
                grid.text("CNV", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
              }
              # annotation names on the left
              # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
            })
          }
        }
      } else {
        if (length(METreg)>0){
          MeanMET <- floor(nrow(METData)/2)+1+1
          for(an in colnames(Genes)) {
            decorate_annotation(an, {         
              if (an=="No"){
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
              } else {
                # annotation names on the right
                grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
              }
              if (an==colnames(Genes)[MeanMET]){
                grid.text("MET", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
              }
              # annotation names on the left
              # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
            })
          }
        }
      }  
    }
  }
