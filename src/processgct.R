source("/Users/liefeld/GenePattern/gp_dev/docker/docker-amaretto/src/common.R")


file.names <- dir(".", pattern ="*DM.gct")
for(i in 1:length(file.names)){
	print(file.names[i])
    code = strsplit(file.names[i], '_')[[1]][1] 
	ov = read.table(file = file.names[i], sep = '\t')
	
	cn2 <- gsub(x = colnames(ov), pattern="\\.", replacement="-")
	colnames(ov) <- cn2
	colnames(ov) <- paste(colnames(ov),"-01", sep = "")
	
	
	gctExp = list(data=as.matrix(ov))
	write.gct(gctExp, file.path("/Users/liefeld/GenePattern/gp_dev/docker/docker-amaretto/test/from_Jay/PanCancer_DMValues", paste("TCGA_",code,"_Methylation.gct", sep="")))
}