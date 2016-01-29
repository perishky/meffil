#' Only a few steps are needed to install meffil in R:
#' 1. install.packages("devtools")
#' 2. library(devtools)
#' 3. source("http://bioconductor.org/biocLite.R")
#' 4. devtools::install_github("perishky/meffil")

#' The steps below are needed only to regenerate
#' the data objects included with the package. 
#' These have already been generated and are included
#' with the package.

#' install.packages("devtools")
#' devtools::install_github("klutometis/roxygen")
library(devtools)
library(roxygen2)

document("meffil")

system("R CMD INSTALL meffil")
reload(inst("meffil"))

source("meffil/data-raw/globals.r",chdir=T)
system("R CMD INSTALL meffil") 

source("meffil/tests/run-all.r", chdir=T)
