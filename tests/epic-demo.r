## Information about the new EPIC microarray can be obtained here:
## http://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html

library(meffil)

source("dataset-epic-demo.r")
path <- download.epic.demo.dataset()

qc.file <- "epic-demo/qc-report.html"
author <- "Illumina, et al."
study <- "EPIC demo dataset"
number.pcs <- 2
norm.file <- "epic-demo/normalization-report.html"

samplesheet <- meffil.read.samplesheet(base=path, pattern="Demo_SampleSheet.csv")

qc.objects <- meffil.qc(samplesheet, cell.type.reference=NULL, verbose=T)

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
