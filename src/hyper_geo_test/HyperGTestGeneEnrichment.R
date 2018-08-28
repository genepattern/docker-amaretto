HyperGTestGeneEnrichment<-function(gmtfile,testgmtfile,outfile,show.overlapping.genes=FALSE,filter.genes=TRUE,show.unrecognized=FALSE )
  
{
library(GSEABase)
require(plyr)
source('/usr/local/bin/amaretto/hyper_geo_test/readGMT.R')
  


# 
# ## Help section
# if("--help" %in% args) {
#   cat("
#       The R Script HyperGTestGeneEnrichment is used to perform the hypergeometric test. 
#  
#       Arguments:
#       1) the signature collection
#       2) the genelist
#       3) the output file
#       4) Default FALSE: output the overlapping genes
#       5) Default TRUE: only count the genes in all_genes.txt (useful for e.g. mouse signatures) 
#       6) Default TRUE: output genes that are not recognized
#       --help              - print this text
#  
#       Example:
#       Rscript HyperGTestGeneEnrichment.R h.all.v5.0.symbols.gmt my_genesignatures.gmt output.txt\n\n")
#   
#   q(save="no")
# }


#MSigDB 45956
ref.num<-45956

all_genes<-scan("/usr/local/bin/amaretto/hyper_geo_test/all_genes.txt", what="", sep="\n")
#ref.num<-length(all_genes)

# gmtfile<-args[1]     #signature collection
# testgmtfile<-args[2] #genelist
# outfile<-args[3]     #output file
# show.overlapping.genes<-args[4]
if(show.overlapping.genes){
  
}
# filter.genes<-as.logical(args[5])
# show.unrecognized<-as.logical(args[6])

test.gmt<-readGMT(testgmtfile)
gmt.path<-readGMT(gmtfile)
out<-c()
out.genes<-c()


if(filter.genes && show.unrecognized){
  genes.in.test.gmt<-unlist(test.gmt$genesets)
  if(sum(!genes.in.test.gmt %in% all_genes)>0) {
    cat("The following genes are not recognized: ",genes.in.test.gmt[!genes.in.test.gmt %in% all_genes],"\n")
  }
}

for(i in 1:length(gmt.path$genesets)){
	for(j in 1:length(test.gmt$genesets)){
	  
	  
	  
		set.num<-length(gmt.path$genesets[[i]])
		k<-sum(gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]])
		l<-set.num
		m<-ref.num
    if(filter.genes){
		  n<-sum(test.gmt$genesets[[j]] %in% all_genes)
    }else{
      n<-length(test.gmt$genesets[[j]])
    }
		p1<-phyper(k-1,l,m-l,n,lower.tail=FALSE)
		r<-c(gmt.path$geneset.names[i],gmt.path$geneset.descriptions[[i]],test.gmt$geneset.names[[j]],ref.num,set.num,n,k,p1)
		out<-rbind(out,r)
    if(show.overlapping.genes){
      overlapping.genes<-gmt.path$genesets[[i]][gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]]]
      
      
      overlapping.genes<-gsub('\t',',',as.character(overlapping.genes))
      overlapping.genes<-gsub('  ',',',overlapping.genes)
      overlapping.genes<-gsub(' ',',',overlapping.genes)
      overlapping.genes<-paste(overlapping.genes,collapse = ', ')
      #if(!identical(overlapping.genes, character(0))){print(overlapping.genes)}

      rr<-c(gmt.path$geneset.names[i],gmt.path$geneset.descriptions[[i]],test.gmt$geneset.names[[j]],p1,k,overlapping.genes)
      out.genes<-rbind(out.genes,rr)
      # write(out.genes,file=outfile.genes,append=T,sep='\t',ncol=length(out.genes))
      # print(c(i,j))
      # print('sham')
    }
		
}
}
pp<-out[,8]
pp.adj<-p.adjust(pp,method='BH')
out<-cbind(out,pp.adj)





col.names<-c('GENESETNAME','GENSETDESCRIPTION','TESTSETNAME','ALLGENE(N)','GENESET(K)','TESTSET(n)','OVERLAP(k)','p-value','q-value(FDR)')
colnames(out)<-col.names
write.table(out,file=outfile,sep='\t',quote=F,col.names=T,row.names=F)

outfile.genes<-paste(tools::file_path_sans_ext(outfile),'.genes.',tools::file_ext(outfile),sep="")
out.genes<-cbind(out.genes,pp.adj)

col.names2<-c('Geneset','Description','Testset','p-value','n-Overlapping','Overlapping genes','q-value')
colnames(out.genes)<-col.names2
write.table(out.genes,file=outfile.genes,sep='\t',quote=F,col.names=T,row.names=F)

return(c(0))
}


