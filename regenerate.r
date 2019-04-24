#' The steps below are needed to regenerate
#' the data objects and documentation files
#' included with the package and then
#' run all tests.

#' install.packages("devtools")
#' devtools::install_github("klutometis/roxygen")
library(devtools)
library(roxygen2)

document("meffil")

system("R CMD INSTALL meffil")
reload(inst("meffil"))

system("R CMD Rd2pdf meffil")
system("mv meffil.pdf meffil/docs")

source("meffil/data-raw/globals.r",chdir=T)
system("R CMD INSTALL meffil") 

