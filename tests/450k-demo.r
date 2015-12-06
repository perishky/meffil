library(meffil)

source("dataset-450k-demo.r")
path <- download.450k.demo.dataset()

samplesheet <- meffil.create.samplesheet(path)

options(mc.cores=3)

meffil.list.cell.type.references()

qc.file <- "450k-demo/qc-report.html"
author <- "Prickett, et al."
study <- "Silver-Russell syndrome patients (GEO:GSE55491)"
number.pcs <- 2
norm.file <- "450k-demo/normalization-report.html"
cell.type.reference <- "blood gse35069"

qc.objects <- meffil.qc(samplesheet, cell.type.reference=cell.type.reference, verbose=T)

qc.summary <- meffil.qc.summary(qc.objects, verbose=T)
meffil.qc.report(qc.summary,
                 output.file=qc.file,
                 author=author,
                 study=study)

if (nrow(qc.summary$bad.samples) > 0)
    qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)


print(meffil.plot.pc.fit(qc.objects, n.cross=3)$data)
##   n        M        U
## 1 1 95142.16 757128.5
## 2 2 28414.64 109172.8
## 3 3 29346.49 109954.2
## 4 4 29714.49 111386.4

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

counts.meffil <- t(meffil.cell.count.estimates(norm.objects))

counts.beta <- meffil.estimate.cell.counts.from.betas(beta.meffil, cell.type.reference)

for (cell.type in colnames(counts.meffil)) {
    cat(cell.type, cor(counts.meffil[,cell.type], counts.beta[,cell.type]), "\n")
}
## Bcell 0.9814654 
## CD4T 0.998171 
## CD8T 0.9966193 
## Gran 0.9993144 
## Mono 0.9818829 
## NK 0.9810345 
