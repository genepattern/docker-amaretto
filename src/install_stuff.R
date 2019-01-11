## try http:// if https:// URLs are not supported
# module specific packages first 
source("http://bioconductor.org/biocLite.R")
install.packages(c("optparse"))

#biocLite("ComplexHeatmap")

#install.packages(c("RCurl", "limma", "foreach", "parallel", "doParallel", "glmnet", "matrixStats", "RColorBrewer", "impute", "Matrix", "BiocStyle", "stringr", "circlize", "R.utils" , "svglite", "callr", "Rcpp", "DT", "htmltools", "reshape2"))

#install.packages(c(   "GSEABase", "rstudioapi", "R2HTML", "plyr" , "tm", "SnowballC", "wordcloud", "snow", "V8", "randomcoloR", "tidyverse"))

biocLite(c("ComplexHeatmap", "curatedTCGAData", "TCGAutils"))

#remove.packages("callr")
#install.packages("callr")
#source(file.path("/usr/local/bin/amaretto/", "callr.R"))

print("Installing devtools")
install.packages('devtools')




