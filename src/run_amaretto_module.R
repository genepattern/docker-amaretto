## The Regents of the University of California and The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2018) by the
## Regents of the University of California abd the 
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# Load any packages used to in our code to interface with GenePattern.
# Note the use of suppressMessages and suppressWarnings here.  The package
# loading process is often noisy on stderr, which will (by default) cause
# GenePattern to flag the job as failing even when nothing went wrong. 
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(AMARETTO)))

# Print the sessionInfo so that there is a listing of loaded packages, 
# the current version of R, and other environmental information in our
# stdout file.  This can be useful for reproducibility, troubleshooting
# and comparing between runs.
sessionInfo()

# Get the command line arguments.  We'll process these with optparse.
# https://cran.r-project.org/web/packages/optparse/index.html
arguments <- commandArgs(trailingOnly=TRUE)

print(packageVersion("AMARETTO"))
# Declare an option list for optparse to use in parsing the command line.
option_list <- list(
  # Note: it's not necessary for the names to match here, it's just a convention
  # to keep things consistent.
  make_option("--expression.file", dest="expression.file"),
  make_option("--copy.number.file", dest="copy.number.file"),
  make_option("--methylation.file", dest="methylation.file"),
  make_option("--output.file", dest="output.file"),
  make_option("--percent.genes", type="integer", dest="percent.genes"),
  make_option("--number.of.modules", type="integer", dest="number.of.modules")
  )

# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=arguments)
print(opt)
opts <- opt$options
# Load some common GP utility code for handling GCT files and so on.  This is included
# with the module and so it will be found in the same location as this script (libdir).
source(file.path("/usr/local/bin/amaretto/", "common.R"))
source(file.path("/usr/local/bin/amaretto/","AMARETTO_VisualizeModulePatch_0.99.2.R"))

# Optparse will validate increment.value and convert it to a numeric value or give it the
# default value of 10 if missing.  We must check for NA however (and NULL, to be safe) as
# it will use that for non-numeric values.
if (is.null(opts$percent.genes) || is.na(opts$percent.genes)) {
   stop("Parameter percent.genes must be numeric")
} else {
   percent.genes <- strtoi(opts$percent.genes)
}
if (is.null(opts$number.of.modules) || is.na(opts$number.of.modules)) {
   stop("Parameter number.of.modules must be numeric")
} else {
   number.of.modules <- strtoi(opts$number.of.modules)
}


patternRange <- seq(1,number.of.modules)


# Load the GCT input file2.
print("Loading gct files now")
if (!file.exists(opts$expression.file)){
     print("Expression file does not exist")
}
gct_exp <- read.gct(opts$expression.file)

if (!file.exists(opts$copy.number.file)){
     print("Copy number file does not exist")
}
gct_cn <- read.gct(opts$copy.number.file)

if (!file.exists(opts$methylation.file)){
     print("Methylation file does not exist")
}
gct_meth <- read.gct(opts$methylation.file)


AMARETTOinit2 = AMARETTO_Initialize(gct_exp$data,gct_cn$data, gct_meth$data,number.of.modules,percent.genes)
AMARETTOresults2 = AMARETTO_Run(AMARETTOinit2)

dirName="."

for (modNum in patternRange)
{
	pdf( paste(opts$output.file,"_module_", modNum, ".pdf"))
	print(AMARETTO_VisualizeModule(AMARETTOinit2,AMARETTOresults2,gct_cn$data,gct_meth$data,ModuleNr=modNum))
	dev.off()
}


tsvFiles = c("NrModules","AllRegulators","AllGenes")
for(res_file in tsvFiles)
{
  resdata <- AMARETTOresults2[[res_file]]
  write.table(resdata,file.path(getwd(),paste(res_file,"_amaretto.tsv",sep = "")),row.names=T,sep="\t",quote=F)#,col.names=F
}

gctFiles = c("ModuleMembership","ModuleData","RegulatoryProgramData","RegulatoryPrograms")
for(res_file in gctFiles)
{
    gct <-list(data=AMARETTOresults2[[res_file]])
    write.gct(gct, file.path(getwd(),paste(res_file,"_amaretto.gct",sep = "")))
}
