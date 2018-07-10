## try http:// if https:// URLs are not supported
# module specific packages first 

source("http://bioconductor.org/biocLite.R")
install.packages(c("optparse"))

biocLite("ComplexHeatmap")

install.packages(c("RCurl", "limma", "foreach", "doParallel", "glmnet", "matrixStats", "RColorBrewer", "impute", "BiocStyle", "circlize", "R.utils" ))

# install dependencies hosted on CRAN
#cran_packages <- c("Rcpp", "RColorBrewer", "gplots", "cluster", "shiny", 
#    "doParallel", "foreach", "ggplot2", "reshape", "testthat", "lintr", "knitr",
#    "rmarkdown", "BH")

#install.packages(cran_packages, repos='http://cran.us.r-project.org')

# install dependencies hosted on BiocConductor
#bioc_packages <- c("BiocStyle")
#biocLite()
#biocLite(bioc_packages)

install.packages("/usr/local/bin/AMARETTO_0.99.1.tar.gz", repos = NULL, type="source")

