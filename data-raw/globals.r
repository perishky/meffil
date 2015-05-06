library(devtools)
library(meffil)

options(mc.cores=10)

probe.info <- meffil:::collate.probe.info(verbose=T)
save(probe.info, file="../inst/probe-info.rda")

meffil.set.probe.info(probe.info)

references <- meffil:::create.gse35069.references(verbose=T)
save(list=names(references),
     file="../inst/gse35069-references.rda",
     envir=list2env(references))


