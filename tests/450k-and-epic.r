library(meffil)

source("dataset-epic-demo.r")
path.epic <- download.epic.demo.dataset()
source("dataset-450k-demo.r")
path.450k <- download.450k.demo.dataset()

options(mc.cores=10)

samplesheet.epic <- meffil.read.samplesheet(path.epic, pattern="Demo_SampleSheet.csv")
samplesheet.450k <-  meffil.create.samplesheet(path.450k)[1:10,]

common.columns <- intersect(colnames(samplesheet.epic), colnames(samplesheet.450k))
samplesheet <- rbind(samplesheet.epic[,common.columns],
                     samplesheet.450k[,common.columns])

qc.file <- "450k-and-epic/qc-report.html"
author <- "Illumina, et al. and Prickett, et al."
study <- "Normalizing EPIC and 450K microarrays together"
number.pcs <- 5
norm.file <- "450k-and-epic/normalization-report.html"
cell.type.reference <- "blood gse35069"
featureset <- "common" ## select a featureset compatible with both 450k and epic

qc.objects <- meffil.qc(samplesheet, cell.type.reference=cell.type.reference,
                        featureset=featureset, verbose=T)

qc.summary <- meffil.qc.summary(qc.objects, verbose=T)
meffil.qc.report(qc.summary,
                 output.file=qc.file,
                 author=author,
                 study=study)

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=number.pcs, verbose=T)

norm.meffil <- meffil.normalize.samples(norm.objects,
                                        just.beta=F, 
                                        cpglist.remove=qc.summary$bad.cpgs$name,
                                        verbose=T)

beta.meffil <- meffil.get.beta(norm.meffil$M, norm.meffil$U)
parameters <- meffil.normalization.parameters(norm.objects)
parameters$batch.threshold <- 0.01
parameters$probe.pcs <- 1:min(10,length(norm.objects))
parameters$control.pcs <- 1:min(10,length(norm.objects))
norm.summary <- meffil.normalization.summary(beta.meffil,
                                             norm.objects=norm.objects,
                                             parameters=parameters, verbose=T)

meffil.normalization.report(norm.summary,
                            output.file=norm.file,
                            author=author,
                            study=study)
