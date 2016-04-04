#!/usr/bin/env Rscript

## Script to install dependencies
install.packages("argparse",repos="http://cran.us.r-project.org")
install.packages("data.table",repos="http://cran.us.r-project.org")
install.packages("depmixS4",repos="http://cran.us.r-project.org")
install.packages("doParallel",repos="http://cran.us.r-project.org")
source("http://bioconductor.org/biocLite.R")
biocLite("rtracklayer")

