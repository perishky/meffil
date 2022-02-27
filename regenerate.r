#' The steps below are needed to regenerate
#' the data objects and documentation files
#' included with the package and then
#' run all tests.

#' install.packages("devtools")
#' devtools::install_github("klutometis/roxygen")
library(devtools)
library(roxygen2)
library(pkgload)

document("meffil")

source("meffil/data-raw/globals.r",chdir=T)

install("meffil")

reload(inst("meffil"))

system("R CMD Rd2pdf meffil")
system("mv meffil.pdf meffil/docs")
