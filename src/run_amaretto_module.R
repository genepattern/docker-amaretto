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
  make_option("--number.of.modules", type="integer", dest="number.of.modules"),
  make_option("--driver.gene.list", dest="driver.gene.list"),
  make_option("--driver.gene.list.file", dest="driver.gene.list.file"),
  make_option("--driver.gene.list.selection.mode", dest="driver.gene.list.selection.mode"),
  make_option("--num.cpu", type="integer", dest="num.cpu"),
  make_option("--hyper.geo.ref.file", dest="hyper.geo.ref.file"),
  make_option("--from.msigdb", type="logical", dest="from.msigdb"),
  make_option("--gmt.url.present", type="logical", dest="gmt.url.present")
  )


# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=arguments)
print(opt)
opts <- opt$options
# Load some common GP utility code for handling GCT files and so on.  This is included
# with the module and so it will be found in the same location as this script (libdir).
source(file.path("/usr/local/bin/amaretto/", "common.R"))

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
gene_list_combination_method = "union"

if (opts$driver.gene.list.selection.mode == "computed") {
   # just compute and ignore any passed in or predefined gene lists
   geneList = NULL
   gene_list_combination_method = "union"
} else {
    # for predefined, union and intersect we need a list
    # if a file was provided, thats the predefined list.  If not, use the driver.gene.list dropdown
    # to pick one of the preset lists 
    
    if ((!is.null(opts$driver.gene.list.file))){
        if (file.exists(opts$driver.gene.list.file)){
 			geneList = as.character(read.delim(opts$driver.gene.list.file)$V1)
 		}
	}	
		
    
    if (is.null(geneList)) {
     
	# get the preformed gene lists.  They use names that match what is passed in except for "Su-In Lee" which 
        # has problematic spaces in the name
	data(Driver_Genes)
	x = names(Driver_Genes)
	x[1] = 'su-in-lee'
	names(Driver_Genes) <- x
        geneList = Driver_Genes[[opts$driver.gene.list]]
    }
    
    if (opts$driver.gene.list.selection.mode == "predefined"){
         gene_list_combination_method = "union"
    } else if (opts$driver.gene.list.selection.mode == "intersect"){
         gene_list_combination_method =  "intersect"     
    } else if (opts$driver.gene.list.selection.mode == "union"){
         gene_list_combination_method = "union"
    }
}

hyper.geo.ref = NULL
from.msigdb = TRUE
gmt.url.present = FALSE

if ((!is.null(opts$hyper.geo.ref.file))){
        if (file.exists(opts$hyper.geo.ref.file)){
             hyper.geo.ref = as.character(read.delim(opts$hyper.geo.ref.file)$V1)
             from.msigdb=opts$from.msigdb
             gmt.url.present=opts$gmt.url.present
        } else {
             print("Hyper geo ref file does not exist")
             hyper.geo.ref="/source/AMARETTO/inst/templates/H.C2CP.genesets.gmt"
             from.msigdb=TRUE
             gmt.url.present=FALSE
        }
} else {
	hyper.geo.ref="/source/AMARETTO/inst/templates/H.C2CP.genesets.gmt"
        from.msigdb=TRUE
        gmt.url.present=FALSE 
}     

patternRange <- seq(1,number.of.modules)


# Load the GCT input file2.
if (!file.exists(opts$expression.file)){
     print("Expression file does not exist")
}
gct_exp <- read.gct(opts$expression.file)

# when the selection mode is "predefined list" then the met and CNV files are optional.  For all other modes
# they are required.  Check now and pop an error if they are needed but not here
if (opts$driver.gene.list.selection.mode == "predefined"){
    if (!is.null(opts$copy.number.file)){
        print("Using predefined gene list but a copy number file was provided.  It will be ignored.")
    } 
    if (!is.null(opts$methylation.file)){
        print("Using predefined gene list but a methylation file was provided.  It will be ignored.")
    }

} else {
    allFilesPresent = TRUE
    if (is.null(opts$copy.number.file) || !file.exists(opts$copy.number.file)){
         print(paste("Required copy number file: '", opts$copy.number.file , "'  does not exist"))
         allFilesPresent = FALSE
    }

    if (is.null(opts$methylation.file) || !file.exists(opts$methylation.file)){
	print(paste("Required methylation file: '", opts$methylation.file , "'  does not exist"))
        allFilesPresent = FALSE
    }

    if (!allFilesPresent){
         # missing files, quit with a non-zero exit code
         quit(status=999)
    }

    gct_cn <- read.gct(opts$copy.number.file)
    gct_meth <- read.gct(opts$methylation.file)

}


AMARETTOinit = AMARETTO_Initialize(gct_exp$data,gct_cn$data, gct_meth$data,number.of.modules,VarPercentage = percent.genes, Driver_list = geneList,  method = gene_list_combination_method)
AMARETTOresults = AMARETTO_Run(AMARETTOinit)


# New output file for Community-amaretto follow on module
# you can now use this function
AMARETTO_ExportResults(AMARETTOinit, AMARETTOresults, ".", Heatmaps = FALSE)


# gmt_file<-"/source/AMARETTO/inst/templates/H.C2CP.genesets.gmt"

AMARETTO_HTMLreport(AMARETTOinit,AMARETTOresults,CNV_matrix=gct_cn$data, MET_matrix=gct_meth$data, VarPercentage=percent.genes, hyper_geo_test_bool=TRUE, hyper_geo_reference=hyper.geo.ref, MSIGDB=from.msigdb, GMTURL=gmt.url.present, output_address='./')






