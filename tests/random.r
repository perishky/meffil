# Adjusting for random effects directly

library(meffil)

source("dataset-450k-demo.r")
samplesheet <- meffil.create.samplesheet(path)

options(mc.cores=5)

qc.file <- "random/qc-report.html"
author <- "Prickett, et al."
study <- "Silver-Russell syndrome patients (GEO:GSE55491)"
number.pcs <- 2
norm.file <- "random/normalization-report.html"
cell.type.reference <- "blood gse35069"

data <- meffil.normalize.dataset(samplesheet,
                                 qc.file="random/qc-report.html",
                                 author=author,
                                 study=study,
                                 number.pcs=number.pcs,
                                 random.effects="Slide",
                                 norm.file="random/normalization-report.html",
                                 cell.type.reference=cell.type.reference,
                                 verbose=T)

