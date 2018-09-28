## try http:// if https:// URLs are not supported
# module specific packages first 

source("http://bioconductor.org/biocLite.R")
install.packages(c("optparse"))

biocLite("ComplexHeatmap")

install.packages(c("RCurl", "limma", "foreach", "doParallel", "glmnet", "matrixStats", "RColorBrewer", "impute", "BiocStyle", "circlize", "R.utils" , "svglite"))

# install.packages("/usr/local/bin/AMARETTO_0.99.1.tar.gz", repos = NULL, type="source")
install.packages(c(   "GSEABase", "rstudioapi", "R2HTML", "plyr" , "tm", "SnowballC", "wordcloud", "snow"))
install.packages("/source/AMARETTO", repos = NULL, type="s:qource")



