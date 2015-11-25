#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

library(devtools)
library(roxygen2)

document("meffil")

system("R CMD INSTALL meffil")
reload(inst("meffil"))

source("meffil/data-raw/globals.r",chdir=T)
system("R CMD INSTALL meffil") 

## source("meffil/tests/run-all.r", chdir=T)

