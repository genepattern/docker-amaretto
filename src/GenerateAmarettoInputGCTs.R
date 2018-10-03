#
# #######  Generate all the gct files
#
source("/Users/liefeld/GenePattern/gp_dev/docker/docker-amaretto/src/common.R")
library(AMARETTO)
TargetDirectory=getwd()


#codes = c("BLCA","BRCA","CESC","CHOL", "COAD","ESCA","GBM","HNSC","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","READ","SARC","STAD","THCA","THYM","UCEC")

#codes = c("ACC","FPPP","KICH","DLBC","MESO","PRAD","SKCM","STAD","TGCT","UCS","UVM")
codes = c("KIRP", "HSNC", "COADREAD")

data(MethylStates)

for (CancerCode in codes){
	# keep going if one code fails and retry later
	try ({
	
		DataSetDirectories=Download_CancerSite(CancerCode, TargetDirectory)
		
		ProcessedData = Preprocess_CancerSite(CancerCode,DataSetDirectories)

		data(ProcessedData)

		# now write a GCT for each example matrix
		# will need another module to do the downloads possibly
		dirName = "/Users/liefeld/GenePattern/gp_dev/docker/docker-amaretto/test"
		output.file = paste("TCGA_", CancerCode)
		
		gctExp <-list(data=ProcessedData$MA_TCGA)
		write.gct(gctExp, file.path(dirName, paste(output.file,"_Expression.gct", sep="")))	

		gctCnv <-list(data=ProcessedData$CNV_TCGA)
		write.gct(gctCnv, file.path(dirName, paste(output.file,"_CNV.gct", sep="")))
	
		gctMet <-list(data=ProcessedData$MET_TCGA)
		write.gct(gctMet, file.path(dirName, paste(output.file,"_Methylation.gct", sep="")))
	})
}

#
# #######  end of generating the gct files
#


