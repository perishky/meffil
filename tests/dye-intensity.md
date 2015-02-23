# Comparison of functional normalization implementations

## Load example data set 

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55491).


```r
dir.create(data.dir <- "data")
```

```
## Warning in dir.create(data.dir <- "data"): 'data' already exists
```

```r
if (length(list.files(data.dir, "*.idat$")) == 0) {
  filename <-  file.path(data.dir, "gse55491.tar")
  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
  cat(date(), "Extracting files from GEO archive.\n")
  system(paste("cd", data.dir, ";", "tar xvf", basename(filename)))
  unlink(filename)
  cat(date(), "Unzipping IDAT files.\n")
  system(paste("cd", data.dir, ";", "gunzip *.idat.gz"))
}
```

## Load and normalize using `meffil`
Load the code, probe annotation and sample filename information.

```r
source("../R/functional-normalization.r")
probes <- meffil.probe.info()
```

```
## [probe.characteristics] Mon Feb 23 01:05:05 2015 extracting I-Red 
## [probe.characteristics] Mon Feb 23 01:05:07 2015 extracting I-Green 
## [probe.characteristics] Mon Feb 23 01:05:07 2015 extracting II 
## [probe.characteristics] Mon Feb 23 01:05:07 2015 extracting Control 
## [meffil.probe.info] Mon Feb 23 01:05:07 2015 reorganizing type information 
## [probe.locations] Mon Feb 23 01:05:35 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```
## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```r
basenames <- meffil.basenames(data.dir)
```

Next we collect controls, perform background and dye bias correction and then
compute quantiles for each sample
using two different dye intensity targets, one very low and the other very high
to see if the resulting beta values are affected.

```r
norm.objects.lo <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename, probes=probes, dye.intensity=500) 
})
norm.objects.hi <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename, probes=probes, dye.intensity=5000) 
})
```

Normalize quantiles from each sample using the control
matrix to identify batch effects.

```r
norm.objects.lo <- meffil.normalize.objects(norm.objects.lo, number.pcs=2)
```

```
## [meffil.normalize.objects] Mon Feb 23 01:11:44 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 23 01:11:44 2015 selecting dye correction reference 
## [meffil.normalize.objects] Mon Feb 23 01:11:44 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 23 01:11:44 2015 normalizing quantiles 
```

```r
norm.objects.hi <- meffil.normalize.objects(norm.objects.hi, number.pcs=2)
```

```
## [meffil.normalize.objects] Mon Feb 23 01:11:45 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 23 01:11:45 2015 selecting dye correction reference 
## [meffil.normalize.objects] Mon Feb 23 01:11:45 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 23 01:11:45 2015 normalizing quantiles 
```

Apply quantile normalization to methylation levels of each sample.

```r
B.lo <- do.call(cbind, mclapply(norm.objects.lo, function(object) {
    meffil.get.beta(meffil.normalize.sample(object, probes)) 
}))
B.hi <- do.call(cbind, mclapply(norm.objects.hi, function(object) {
    meffil.get.beta(meffil.normalize.sample(object, probes)) 
})) 
```


```r
quantile(B.lo - B.hi)
```

```
##            0%           25%           50%           75%          100% 
## -3.330669e-16  0.000000e+00  0.000000e+00  0.000000e+00  3.330669e-16
```
