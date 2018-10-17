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
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(plyr)))

# Print the sessionInfo so that there is a listing of loaded packages, 
# the current version of R, and other environmental information in our
# stdout file.  This can be useful for reproducibility, troubleshooting
# and comparing between runs.
sessionInfo()
suppressMessages(suppressWarnings(library('svglite')))
suppressMessages(suppressWarnings(library(R2HTML)))
suppressMessages(suppressWarnings(library("GSEABase")))
suppressMessages(suppressWarnings(library("rstudioapi")))
suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(doParallel)))
suppressMessages(suppressWarnings(require("tm")))
suppressMessages(suppressWarnings(require("SnowballC")))
suppressMessages(suppressWarnings(require("wordcloud")))
suppressMessages(suppressWarnings(require("RColorBrewer")))
suppressMessages(suppressWarnings(require(plyr)))

##################
#
# XXX Temp workaround - we need to copy /usr/local/bin/amaretto/hyper_geo_test/all_genes.txt
# to a hyper_geo_test directory under the $PWD.  Why this is not checked into the AMARETTO
# github until it is not needed anymore (some day they say) is a mystery to me.
#
#################
dir.create('./hyper_geo_test')
file.copy('/usr/local/bin/amaretto/hyper_geo_test/all_genes.txt', './hyper_geo_test/all_genes.txt')
file.copy('/usr/local/bin/amaretto/hyper_geo_test/H.C2CP.genesets.gmt', './hyper_geo_test/H.C2CP.genesets.gmt')
################
#
# XXX end of Temp Workaround 
#
###############

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
  make_option("--number.of.modules", type="integer", dest="number.of.modules"),
  make_option("--driver.gene.list", dest="driver.gene.list"),
  make_option("--driver.gene.list.file", dest="driver.gene.list.file"),
  make_option("--driver.gene.list.selection.mode", dest="driver.gene.list.selection.mode"),
  make_option("--num.cpu", type="integer", dest="num.cpu")
  )

# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=arguments)
print(opt)
opts <- opt$options
# Load some common GP utility code for handling GCT files and so on.  This is included
# with the module and so it will be found in the same location as this script (libdir).
source(file.path("/usr/local/bin/amaretto/", "common.R"))
#source(file.path("/usr/local/bin/amaretto/","AMARETTO_VisualizeModulePatch_0.99.2.R"))

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

#
# if driver gene list selection mode != computed, we need to specify the list and/or load a file
#
intersect = FALSE
geneList = NULL
gene_list_combination_method = NULL

if (opts$driver.gene.list.selection.mode == "computed") {
   # just compute and ignore any passed in or predefined gene lists
   geneList = NULL
   gene_list_combination_method = NULL
} else {
    # for predefined, union and intersect we need a list
    # if a file was provided, thats the predefined list.  If not, use the driver.gene.list dropdown
    # to pick one of the preset lists 
    
    if ((!is.null(opts$driver.gene.list.file))){
        if (file.exists(opts$driver.gene.list.file)){
 			geneList = readLines(opts$driver.gene.list.file)
 		}
	}	
		
    
    if (is.null(geneList)) {
     
		# get the preformed gene lists.  They use names that match what is passed in except for "Su-In Lee" which 
        # has problematic spaces in the name
		data(Driver_Genes)
		x = names(Driver_Genes)
		x[1] = 'su-in-lee'
		names(Driver_Genes) <- x
	    geneList = Driver_Genes[opts$driver.gene.list][]
    }
    
    if (opts$driver.gene.list.selection.mode == "predefined"){
         gene_list_combination_method = NULL
    } else if (opts$driver.gene.list.selection.mode == "intersect"){
         gene_list_combination_method =  "intersect"     
    } else if (opts$driver.gene.list.selection.mode == "union"){
         gene_list_combination_method = "union"
    }
}

patternRange <- seq(1,number.of.modules)


# Load the GCT input file2.
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

# num.cpu is the full number allocated to the container.  If there are a few, we take one off the top for the OS/Container itself to avoid thrashing
NrCores = as.integer( 1 * opts$num.cpu)

#
# Driver list - if provided neglects MET and/or CNV data
#    Su-In Lee or MSigDB or allow user to provide a file
# If CNV/MET and a list provided, can use generated, provided, or intersection of both
#
print("===== gene list following ==== ")
print(geneList)
print("===== gene combo method following ==== ")
print(gene_list_combination_method)

AMARETTOinit = AMARETTO_Initialize(gct_exp$data,gct_cn$data, gct_meth$data,number.of.modules,percent.genes, Driver_list = geneList, NrCores = NrCores, method = gene_list_combination_method)
AMARETTOresults = AMARETTO_Run(AMARETTOinit)

dirName="."


#
# Below is old code for exporting results.  This now seems to be redundant with the data under the report_html
#

#for (modNum in patternRange)
#{
#	pdf( paste(opts$output.file,"_module_", modNum, ".pdf", sep = " "))
#	print(AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults,gct_cn$data,gct_meth$data,ModuleNr=modNum))
#	dev.off()
#}
#tsvFiles = c("NrModules","AllRegulators","AllGenes")
#for(res_file in tsvFiles)
#{
#  resdata <- AMARETTOresults[[res_file]]
#  write.table(resdata,file.path(getwd(),paste(res_file,"_amaretto.tsv",sep = "")),row.names=T,sep="\t",quote=F)#,col.names=F
#}
#
#gctFiles = c("ModuleMembership","ModuleData","RegulatoryProgramData","RegulatoryPrograms")
#for(res_file in gctFiles)
#{
#    gct <-list(data=AMARETTOresults[[res_file]])
#    write.gct(gct, file.path(getwd(),paste(res_file,"_amaretto.gct",sep = "")))
#}

# New output file for Community-amaretto follow on module
print(paste("About to save results ", file.path(getwd(),paste(opts$output.file,".RData",sep = ""))))
save(AMARETTOresults, file=file.path(getwd(),paste(opts$output.file,".RData",sep = "")))

dir.create("hyper_geo_test")

res<-amaretto_html_report(AMARETTOinit,AMARETTOresults,gct_cn$data, gct_meth$data,percent.genes,hyper_geo_test_bool=TRUE,n_cluster=NrCores,wordcloud_bool=FALSE)


# copy the tsv files to the top for easy acccess in the notebook
# report_html/htmls/tables/module_hyper_geo_test/
#flist2 <- list.files("report_html/htmls/tables/module_hyper_geo_test/",  full.names = TRUE)
#file.copy(flist2, ".")

#unlink("report_html", recursive = TRUE)





