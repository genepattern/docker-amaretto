rm(list=ls(all=TRUE))
library("AMARETTO")
source('./mohsen_report_function.R')
source('./hyper_geo_test/read_gct.R')
########################################################
# Demo data
########################################################

data(ProcessedDataLAML)
MA_matrix<-data.matrix(ProcessedDataLAML$MA_TCGA)
CNV_matrix<-data.matrix(ProcessedDataLAML$CNV_TCGA)
MET_matrix<-data.matrix(ProcessedDataLAML$MET_TCGA)



dim(MA_matrix)

NrModules <- 10
VarPercentage <- 20

########################################################
# LIHC Data
########################################################
MA_matrix <- data.matrix(read_gct('TCGA_LIHC_Expression.gct'))
CNV_matrix <-data.matrix(read_gct('TCGA_LIHC_CNV.gct'))
MET_matrix <-data.matrix(read_gct('TCGA_LIHC_Methylation.gct'))


NrModules <- 150
VarPercentage <- 75
########################################################
# Running AMARETTO
########################################################

AMARETTOinit<-AMARETTO_Initialize(MA_matrix=MA_matrix,
                                  CNV_matrix=CNV_matrix,MET_matrix=MET_matrix, 
                                  NrModules=NrModules,VarPercentage=VarPercentage)

AMARETTOresults<-AMARETTO_Run(AMARETTOinit)

save(AMARETTOinit,file="hyper_geo_test/AMARETTOinit.Rda")
save(AMARETTOresults,file="hyper_geo_test/AMARETTOresults.Rda")

save(AMARETTOinit,file="LIHC_AMARETTO_RESULTS/AMARETTOinit.Rda")
save(AMARETTOresults,file="LIHC_AMARETTO_RESULTS/AMARETTOresults.Rda")


################################################################################################################
################################################################################################################


dim(MA_matrix)
dim(CNV_matrix)
dim(MET_matrix)


load(file="hyper_geo_test/AMARETTOinit.Rda")
load(file="hyper_geo_test/AMARETTOresults.Rda")

########################## LIHC
load(file="LIHC_AMARETTO_RESULTS/AMARETTOinit.Rda")
load(file="LIHC_AMARETTO_RESULTS/AMARETTOresults.Rda")
#################################


################################################################################################################
################################################################################################################

res<-amaretto_html_report(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,TRUE)






