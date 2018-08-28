rileen <- function(Rdata_address, AMARETTOinit, AMARETTOresults){

# file_wd=dirname(rstudioapi::getSourceEditorContext()$path)
# setwd(file_wd)
# file_wd

# Rdata_address
readRDS(Rdata_address)

all_genes<-unique(c(AMARETTOresults$AllGenes,AMARETTOresults$AllRegulators))
annotation<-data.frame(Gene=all_genes)
annotation$Regulator<-annotation$Gene%in%AMARETTOresults$AllRegulators

all_data<-unique(rbind(AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData))


rp<-t(AMARETTOresults$RegulatoryPrograms)
ma<-AMARETTOresults$ModuleMembership

module_weights<-matrix(0,nrow=length(all_genes),ncol=dim(rp)[2])
rownames(module_weights)<-all_genes
colnames(module_weights)<-colnames(rp)
for(i in 1:dim(module_weights)[1]){
  g<-rownames(module_weights)[i]
  w<-rep(0,dim(rp)[2])
  if(g %in% rownames(rp))
    w<-w+rp[g,]
  if(g %in% rownames(ma))
    w[ma[g,]]<-1+w[ma[g,]]
  module_weights[g,]<-w
}

annotation<-cbind(annotation,module_weights)

in_module<-module_weights!=0
colnames(in_module)<-unlist(lapply(colnames(in_module),paste0,"_B"))
annotation<-cbind(annotation,in_module)


gmt_file="hyper_geo_test/TCGA_modules_target_only.gmt"
for(i in 1:dim(module_weights)[2]){
  genes<-rownames(in_module)[which(module_weights[,i]==1)]
  if(i==1){
    write(c(colnames(module_weights)[i],colnames(module_weights)[i],genes),file=gmt_file,sep = "\t",ncolumns=length(genes)+2)
  }else{
    write(c(colnames(module_weights)[i],colnames(module_weights)[i],genes),file=gmt_file,sep = "\t",ncolumns=length(genes)+2,append = T)
  }
}

aa=c(1)
return(aa)

}
