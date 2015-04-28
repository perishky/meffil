library(meffil)
options(mc.cores=24)


dir.create(path <- "~/data/test_meffil", recursive=TRUE)

if (length(list.files(path, "*.idat$")) == 0) {
  filename <-  file.path(path, "gse55491.tar")
  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
  cat(date(), "Extracting files from GEO archive.\n")
  system(paste("cd", path, ";", "tar xvf", basename(filename)))
  unlink(filename)
  cat(date(), "Unzipping IDAT files.\n")
  system(paste("cd", path, ";", "gunzip *.idat.gz"))
}

# B <- meffil.normalize.dataset(path=path, number.pcs=2)

basenames <- meffil.basenames(path)

samplesheet <- data.frame(do.call(rbind, strsplit(basename(basenames), split="_")), basenames)
names(samplesheet) <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position", "Basename")


samplesheet <- data.frame(Sample_Name = paste("id",1:24, sep=""), sex=sample(c("M", "F"), 24, replace=T))



norm.objects <- mclapply(basenames[1:5], meffil.compute.normalization.object)
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=10)


res <- meffil.pre.processing(
	samplesheet,
	norm.objects,
	colour.code = NULL,
	control.categories=names(norm.objects[[1]]$controls),
	sex.outlier.sd = 3,
	meth.unmeth.outlier.sd = 3,
	control.means.outlier.sd = 5,
	detectionp.samples.threshold = 0.05,
	beadnum.samples.threshold = 0.05,
	detectionp.cpgs.threshold = 0.05,
	beadnum.cpgs.threshold = 0.05
)


B.long <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object))
}))

design.matrix <- meffil.design.matrix(norm.objects)





# Potential alternative workflow:


# Get sample sheet
# - Sample_Name
# - Sex
# - Batch variables (optional)
samplesheet <- meffil.read.samplesheet(path)
# or
samplesheet <- meffil.create.samplesheet(path)

# Do background and dye correction
qc.objects <- meffil.qc(samplesheet)

# Find individuals and probes that should be removed prior to performing normalisation
# Generate html doc with some graphs
qc.report <- meffil.qc.report(qc.objects)

# Remove individuals (and probes at this stage? Maybe not because bad probes shouldn't have a big effect on normalisation I don't think)
qc.objects <- meffil.remove.ids(qc.objects, qc.report$bad.ids)

# Normalise
norm.objects <- meffil.normalize.objects(qc.objects)

# Generate betas
B <- meffil.normalize.samples(norm.objects)

# Generate html doc with summary of normalisation
post.normalization.report <- meffil.post.normalization(samplesheet, B, norm.objects)



library(meffil)
options(mc.cores=16)

samplesheet <- meffil.create.samplesheet("~/data/test_meffil")
# samplesheet$Sex <- 
# samplesheet$Sample_Name <- 
qc.objects <- meffil.qc(samplesheet, verbose=TRUE)
qc.summary <- meffil.qc.summary(qc.objects)
meffil.qc.report(qc.summary)
qc.objects <- meffil.remove.ids(qc.objects, qc.summary$bad.ids)
qc.quantiles <- meffil.normalize.quantiles(qc.objects, number.pcs=2)
normalized.beta <- meffil.normalize.samples(qc.quantiles, cpglist.remove=qc.report$bad.cpgs)
normalization.summary <- meffil.normalization.summary(normalized.beta, qc.quantiles)
meffil.normalization.report(normalization.summary)
