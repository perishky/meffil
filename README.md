# meffil
Efficient algorithms for analyzing DNA methylation data.

## Functional normalization
The `minfi` version consists of a single function call
that requires access to all data.
The `meffil` version splits the method up into several functions
in order to allow parallelization
and to reduce the amount of data that needs to be loaded.

### `minfi` 
Code for applying the `minfi` functional normalization algorithm
to a dataset (source directory unspecified).
```r
library(minfi)
data.dir <- ....
example <- read.450k.exp(data.dir)
example.norm <- preprocessFunnorm(example, nPCs=2, sex=NULL,
                                  bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE)
```
### `meffil`
Load the `meffil` code, probe annotation and
raw data file information.
```r
source("R/functional-normalization.r")
probes <- meffil.probe.info()
basenames <- meffil.basenames(data.dir)
```
Extract the control information for each sample.
```r
control.matrix <- sapply(basenames, function(basename) {
    meffil.extract.controls(basename, probes) 
})
```

Background and dye correct each sample and the compute probe quantiles
for each sample.
```r
norm.objects <- lapply(basenames, function(basename) {
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


