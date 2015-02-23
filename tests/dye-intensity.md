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
library(meffil)
```

```
## Loading required package: illuminaio
## Loading required package: MASS
## Loading required package: IlluminaHumanMethylation450kmanifest
## Loading required package: minfi
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: lattice
## Loading required package: GenomicRanges
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: Biostrings
## Loading required package: XVector
## Loading required package: bumphunter
## Loading required package: foreach
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
## Loading required package: iterators
## Loading required package: locfit
## locfit 1.5-9.1 	 2013-03-22
## Setting options('download.file.method.GEOquery'='curl')
## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```r
probes <- meffil.probe.info()
```

```
## [probe.characteristics] Mon Feb 23 12:40:32 2015 extracting I-Red 
## [probe.characteristics] Mon Feb 23 12:40:35 2015 extracting I-Green 
## [probe.characteristics] Mon Feb 23 12:40:35 2015 extracting II 
## [probe.characteristics] Mon Feb 23 12:40:35 2015 extracting Control 
## [meffil.probe.info] Mon Feb 23 12:40:35 2015 reorganizing type information 
## [probe.locations] Mon Feb 23 12:41:02 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
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
## [meffil.normalize.objects] Mon Feb 23 12:48:51 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 23 12:48:51 2015 selecting dye correction reference 
## [meffil.normalize.objects] Mon Feb 23 12:48:51 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 23 12:48:51 2015 normalizing quantiles 
## [FUN] Mon Feb 23 12:48:51 2015 genomic.iG M 
## [FUN] Mon Feb 23 12:48:52 2015 genomic.iG U 
## [FUN] Mon Feb 23 12:48:52 2015 genomic.iR M 
## [FUN] Mon Feb 23 12:48:52 2015 genomic.iR U 
## [FUN] Mon Feb 23 12:48:52 2015 genomic.ii M 
## [FUN] Mon Feb 23 12:48:52 2015 genomic.ii U 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.iG M 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.iG U 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.iR M 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.iR U 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.ii M 
## [FUN] Mon Feb 23 12:48:52 2015 autosomal.ii U 
## [FUN] Mon Feb 23 12:48:52 2015 not.y.iG M 
## [FUN] Mon Feb 23 12:48:52 2015 not.y.iG U 
## [FUN] Mon Feb 23 12:48:52 2015 not.y.iR M 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.iR U 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.ii M 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.ii U 
## [FUN] Mon Feb 23 12:48:53 2015 sex M 
## [FUN] Mon Feb 23 12:48:53 2015 sex U 
## [FUN] Mon Feb 23 12:48:53 2015 chry M 
## [FUN] Mon Feb 23 12:48:53 2015 chry U 
## [FUN] Mon Feb 23 12:48:53 2015 chrx M 
## [FUN] Mon Feb 23 12:48:53 2015 chrx U
```

```r
norm.objects.hi <- meffil.normalize.objects(norm.objects.hi, number.pcs=2)
```

```
## [meffil.normalize.objects] Mon Feb 23 12:48:53 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 23 12:48:53 2015 selecting dye correction reference 
## [meffil.normalize.objects] Mon Feb 23 12:48:53 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 23 12:48:53 2015 normalizing quantiles 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.iG M 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.iG U 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.iR M 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.iR U 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.ii M 
## [FUN] Mon Feb 23 12:48:53 2015 genomic.ii U 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.iG M 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.iG U 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.iR M 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.iR U 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.ii M 
## [FUN] Mon Feb 23 12:48:53 2015 autosomal.ii U 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.iG M 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.iG U 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.iR M 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.iR U 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.ii M 
## [FUN] Mon Feb 23 12:48:53 2015 not.y.ii U 
## [FUN] Mon Feb 23 12:48:53 2015 sex M 
## [FUN] Mon Feb 23 12:48:53 2015 sex U 
## [FUN] Mon Feb 23 12:48:53 2015 chry M 
## [FUN] Mon Feb 23 12:48:53 2015 chry U 
## [FUN] Mon Feb 23 12:48:53 2015 chrx M 
## [FUN] Mon Feb 23 12:48:53 2015 chrx U
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
