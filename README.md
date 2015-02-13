# meffil
Efficient algorithms for analyzing DNA methylation data.

Code for applying the `minfi` functional normalization algorithm
to a dataset (source directory unspecified).
```r
library(minfi)
data.dir <- ....
sheet <- read.450k.sheet(data.dir, pattern="samplesheet")
example <- read.450k.exp(base ="/", targets = sheet)
example.norm <- preprocessFunnorm(example, nPCs=2, sex=NULL,
                                  bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE)
```

Load the `meffil` code and probe annotation.
```r
source("normalize-450k.r")
probes <- meffil.probe.info()
```

Extract the control information for each sample.
```r
control.matrix <- sapply(sheet$Basename, function(basename) {
    meffil.extract.controls(basename, probes) 
})
```

Background and dye correct each sample and the compute probe quantiles
for each sample.
```r
norm.objects <- lapply(sheet$Basename, function(basename) {
    meffil.compute.normalization.object(basename, control.matrix, probes=probes)
})
```

Normalize the resulting quantiles together. 
```r
norm.objects <- meffil.normalize.objects(norm.objects, control.matrix,
                                         number.pcs=2, probes=probes)
```

Transform the methylation data to fit the normalized quantiles
and return resulting beta-values.
```r
B <- sapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object, probes)) 
})
```


