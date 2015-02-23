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
### `meffil` (short version)
Here is the simplest way to apply the `meffil` version.
```r
library(meffil)
B <- meffil.normalize.dataset(data.dir=data.dir, number.pcs=2)
```
In this implementation, data for each sample is handled one at a time
until the final step when normalized beta values are merged into
a single dataset-wide matrix `B`.
Consequently memory use is minimized.

### `meffil` (long version)
The long `meffil` version is available in order to
improve performance using parallel computing.
The example below uses `mclapply` to spread computation across multiple processors.

Load the `meffil` code, probe annotation and
raw data file information.
```r
library(meffil)
basenames <- meffil.basenames(data.dir)
```

Background and dye correct each sample and the compute probe controls and quantiles
for each sample.
```r
norm.objects <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename)
})
```

Normalize the resulting quantiles together. 
```r
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=2)
```

Transform the methylation data to fit the normalized quantiles
and return resulting beta-values.
```r
B.long <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object)) 
}))
```


