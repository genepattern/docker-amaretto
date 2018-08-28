read_gct<-function(file_address){
  data_fr=read.table(file_address, skip = 2,sep = '\t',header=TRUE,row.names = 1)
  data_fr <- subset(data_fr, select = -c(Description))
}

