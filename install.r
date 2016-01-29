#' Only two steps are needed to install meffil in R:
#' 1. library(devtools)
#' 2. devtools::install_github("perishky/meffil")

#' The following steps are needed only to regenerate
#' the included with the package. This data has already 
#' present in the package.

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
