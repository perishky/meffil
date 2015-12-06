library(meffil)

source("dataset-450k-demo.r")
path <- download.450k.demo.dataset()
samplesheet <- meffil.create.samplesheet(path)

options(mc.cores=3)

qc.file <- "ewas/qc-report.html"
author <- "Prickett, et al."
study <- "Silver-Russell syndrome patients (GEO:GSE55491)"
number.pcs <- 2
norm.file <- "ewas/normalization-report.html"
cell.type.reference <- "blood gse35069"


data <- meffil.normalize.dataset(samplesheet,
                                 qc.file=qc.file,
                                 author=author,
                                 study=study,
                                 number.pcs=number.pcs,
                                 norm.file=norm.file,
                                 cell.type.reference=cell.type.reference,
                                 verbose=T)
samples <- read.csv(file.path(path, "samples.csv"), row.names=1)
samples <- samples[match(colnames(data$beta), rownames(samples)),]

table(samples$sex,
      with(data$qc.summary$sex.summary,
           tab[match(rownames(samples),tab$sample.name), "predicted.sex"]))

counts <- meffil.cell.count.estimates(data$norm.objects)
counts <- t(counts[,match(colnames(data$beta), colnames(counts))])

variable <- samples$sex
covariates <- data.frame(group=samples$group, counts)
ewas.ret <- meffil.ewas(data$beta, variable=variable, covariates=covariates, most.variable=50000, verbose=T)

ewas.parameters <- meffil.ewas.parameters(sig.threshold=1e-20, max.plots=5)
candidate.sites <- c("cg04946709","cg06710937","cg12177922","cg15817705","cg20299935","cg21784396")
ewas.summary <- meffil.ewas.summary(ewas.ret,
                                    data$beta,
                                    selected.cpg.sites=candidate.sites,
                                    parameters=ewas.parameters, verbose=T)
meffil.ewas.report(ewas.summary,
                   output.file="ewas/ewas-report.html",
                   author="Me!",
                   study="Sex differences (GEO:GSE55491)")

