library(meffil)
options(mc.cores=24)


path <- "~/data/test_meffil/"
B <- meffil.normalize.dataset(path=path, number.pcs=2)


basenames <- meffil.basenames(path)
samplesheet <- data.frame(Sample_Name = paste("id",1:24, sep=""), sex=sample(c("M", "F"), 24, replace=T))

norm.objects <- mclapply(basenames, meffil.compute.normalization.object, detection.threshold = 0.05)
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=10)


B.long <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object))
}))


a <- sapply(norm.objects, function(x) x$xy.diff)
b <- sapply(norm.objects, function(x) x$predicted.sex)

boxplot(a ~ b)

