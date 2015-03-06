# Comparison of functional normalization implementations

## Download example data set 


```r
source("get-dataset.r")
```

```
## Warning in dir.create(path <- "data"): 'data' already exists
```

## Load and normalize using `minfi`

The current version `minfi::preprocessFunnorm()` contains a bug (Feb 25, 2015).
To make a proper comparison with `meffil`, this bug should be fixed
and `minfi` reinstalled.
This can be done with the included patch file: `fix-minfi-control-matrix.patch`.

Load the data.


```r
library(minfi, quietly=TRUE)
```

```
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
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: XVector
## Loading required package: foreach
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
## Loading required package: iterators
## Loading required package: locfit
## locfit 1.5-9.1 	 2013-03-22
## Setting options('download.file.method.GEOquery'='curl')
```

```r
raw.minfi <- read.450k.exp(base = path)
```

Normalize using the minfi package.

```r
system.time(norm.minfi <- preprocessFunnorm(raw.minfi, nPCs=2,
                                            sex=NULL, bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE))
```

```
## [preprocessFunnorm] Background and dye bias correction with noob
```

```
## Loading required package: IlluminaHumanMethylation450kmanifest
## Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```
## Using sample number 4 as reference level...
## [preprocessFunnorm] Mapping to genome
## [preprocessFunnorm] Quantile extraction
## [preprocessFunnorm] Normalization
```

```
##    user  system elapsed 
## 113.376   1.871 115.530
```

## Load and normalize using `meffil`


```r
options(mc.cores=10)
```

Load the code, probe annotation and sample filename information.

```r
library(meffil)
```

```
## Loading required package: illuminaio
## Loading required package: MASS
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
probes <- meffil.probe.info()
basenames <- meffil.basenames(path)
```

Collect controls, perform background and dye bias correction and then
compute quantiles for each sample, one at a time.

```r
norm.objects <- mclapply(basenames, meffil.compute.normalization.object)
```
 
Normalize quantiles from each sample using the control
matrix to identify batch effects.

```r
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=2)
```

```
## [meffil.normalize.objects] Fri Mar  6 02:12:37 2015 selecting dye correction reference 
## [meffil.normalize.objects] Fri Mar  6 02:12:37 2015 predicting sex 
## [meffil.normalize.objects] Fri Mar  6 02:12:37 2015 creating control matrix 
## [meffil.normalize.objects] Fri Mar  6 02:12:37 2015 normalizing quantiles 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.iG M 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.iG U 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.iR M 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.iR U 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.ii M 
## [FUN] Fri Mar  6 02:12:37 2015 genomic.ii U 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.iG M 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.iG U 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.iR M 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.iR U 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.ii M 
## [FUN] Fri Mar  6 02:12:37 2015 autosomal.ii U 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.iG M 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.iG U 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.iR M 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.iR U 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.ii M 
## [FUN] Fri Mar  6 02:12:37 2015 not.y.ii U 
## [FUN] Fri Mar  6 02:12:37 2015 sex M 
## [FUN] Fri Mar  6 02:12:37 2015 sex U 
## [FUN] Fri Mar  6 02:12:37 2015 chry M 
## [FUN] Fri Mar  6 02:12:37 2015 chry U 
## [FUN] Fri Mar  6 02:12:37 2015 chrx M 
## [FUN] Fri Mar  6 02:12:37 2015 chrx U
```

Apply quantile normalization to methylation levels of each sample.

```r
B.meffil <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object)) 
})) 
```

Make sure both methods for normalizing samples produce the same output.

```r
B.meffil.2 <- meffil.normalize.samples(norm.objects)
B.meffil.3 <- meffil.normalize.dataset(path, number.pcs=2)
```

```
## [meffil.normalize.objects] Fri Mar  6 02:15:02 2015 selecting dye correction reference 
## [meffil.normalize.objects] Fri Mar  6 02:15:02 2015 predicting sex 
## [meffil.normalize.objects] Fri Mar  6 02:15:02 2015 creating control matrix 
## [meffil.normalize.objects] Fri Mar  6 02:15:02 2015 normalizing quantiles 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.iG M 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.iG U 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.iR M 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.iR U 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.ii M 
## [FUN] Fri Mar  6 02:15:02 2015 genomic.ii U 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.iG M 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.iG U 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.iR M 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.iR U 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.ii M 
## [FUN] Fri Mar  6 02:15:02 2015 autosomal.ii U 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.iG M 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.iG U 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.iR M 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.iR U 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.ii M 
## [FUN] Fri Mar  6 02:15:02 2015 not.y.ii U 
## [FUN] Fri Mar  6 02:15:02 2015 sex M 
## [FUN] Fri Mar  6 02:15:02 2015 sex U 
## [FUN] Fri Mar  6 02:15:02 2015 chry M 
## [FUN] Fri Mar  6 02:15:02 2015 chry U 
## [FUN] Fri Mar  6 02:15:02 2015 chrx M 
## [FUN] Fri Mar  6 02:15:02 2015 chrx U
```

```r
quantile(B.meffil - B.meffil.2)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

```r
quantile(B.meffil - B.meffil.3)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

## Compare normalizations

First make sure that the raw data is the same.

```r
raw.meffil <- meffil.read.rg(basenames[1])
```

```
## [read.idat] Fri Mar  6 02:15:43 2015 Reading data/GSM1338100_6057825094_R01C01_Grn.idat 
## [read.idat] Fri Mar  6 02:15:44 2015 Reading data/GSM1338100_6057825094_R01C01_Red.idat
```

```r
all(raw.meffil$R == getRed(raw.minfi)[names(raw.meffil$R),1])
```

```
## [1] TRUE
```

### Background correction

Apply `meffil` background correction.

```r
bc.meffil <- meffil.background.correct(raw.meffil)
M.meffil <- meffil.rg.to.mu(bc.meffil)$M
```

Apply `minfi` background correction.

```r
bc.minfi <- preprocessNoob(raw.minfi, dyeCorr = FALSE)
M.minfi <- getMeth(bc.minfi)[names(M.meffil),1]
```

Compare results, should be identical.

```r
quantile(M.meffil - M.minfi)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

### Dye bias correction
Apply `meffil` correction.

```r
dye.meffil <- meffil.dye.bias.correct(bc.meffil,
                                      intensity=norm.objects[[1]]$reference.intensity)
M.meffil <- meffil.rg.to.mu(dye.meffil)$M
U.meffil <- meffil.rg.to.mu(dye.meffil)$U
```

Apply `minfi` background and dye bias correction.

```r
dye.minfi <- preprocessNoob(raw.minfi, dyeCorr = TRUE)
M.minfi <- getMeth(dye.minfi)[names(M.meffil),1]
U.minfi <- getUnmeth(dye.minfi)[names(U.meffil),1]
```

Compare results, should be identical.

```r
quantile(M.meffil - M.minfi)
```

```
##            0%           25%           50%           75%          100% 
## -7.275958e-12  0.000000e+00  0.000000e+00  0.000000e+00  7.275958e-12
```

```r
quantile(U.meffil - U.minfi)
```

```
##            0%           25%           50%           75%          100% 
## -7.275958e-12  0.000000e+00  0.000000e+00  0.000000e+00  7.275958e-12
```

```r
quantile(M.meffil/(M.meffil+U.meffil+100) - M.minfi/(M.minfi+U.minfi+100))
```

```
##            0%           25%           50%           75%          100% 
## -2.220446e-16  0.000000e+00  0.000000e+00  0.000000e+00  2.220446e-16
```

### Control matrix extracted

Extract `minfi` control matrix.

```r
library(matrixStats, quietly=TRUE)
```

```
## 
## Attaching package: 'matrixStats'
## 
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```r
control.matrix.minfi <- minfi:::.buildControlMatrix450k(minfi:::.extractFromRGSet450k(raw.minfi))
```

Control matrix for `meffil` was extracted earlier.
Preprocess the matrix using the same procedure used by `minfi`.

```r
control.matrix.meffil <- t(meffil.control.matrix(norm.objects))
```

The control variable order (columns) may be different between `minfi` and `meffil`
so they are reordered.

```r
i <- apply(control.matrix.minfi,2,function(v) {
    which.min(apply(control.matrix.meffil,2,function(w) {
        sum(abs(v-w))
    }))
})
control.matrix.meffil <- control.matrix.meffil[,i]
```

Compare resulting matrices, should be identical.

```r
quantile(control.matrix.meffil - control.matrix.minfi)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

### Final beta values


```r
B.minfi <- getBeta(norm.minfi)
quantile(B.meffil - B.minfi[rownames(B.meffil),])
```

```
##            0%           25%           50%           75%          100% 
## -3.193965e-01 -1.579708e-08 -2.670409e-10  0.000000e+00  2.225038e-01
```

On what chromosomes are the biggest differences appearing?

```r
probe.chromosome <- probes$chr[match(rownames(B.meffil), probes$name)]
sex <- sapply(norm.objects, function(object) object$sex)
is.diff <- abs(B.meffil - B.minfi[rownames(B.meffil),]) > 1e-4
table(probe.chromosome[which(is.diff, arr.ind=T)[,"row"]])
```

```
## 
##   chrX   chrY 
## 166640   9837
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="M")], arr.ind=T)[,"row"]])
```

```
## 
##  chrX  chrY 
## 61671  5308
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="F")], arr.ind=T)[,"row"]])
```

```
## 
##   chrX   chrY 
## 104969   4529
```
These results not surprising
because chromosome Y not handled like `minfi` in either males or females, 
and chromosome X is handled just like `minfi` in females but not in males.


```r
autosomal.cgs <- unique(probes$name[which(probes$chr %in% paste("chr", 1:22, sep=""))])
quantile(B.meffil[autosomal.cgs,] - B.minfi[autosomal.cgs,])
```

```
##            0%           25%           50%           75%          100% 
## -7.580066e-08 -1.537428e-08 -2.670804e-10  0.000000e+00  4.986327e-09
```

In spite of the differences, CG correlations between methods are pretty close to 1.

```r
male.idx <- which(rowSums(is.diff[,sex=="M"]) >= 5)
male.cg.r <- unlist(mclapply(male.idx, function(idx) {
    cor(B.meffil[idx, sex=="M"], B.minfi[rownames(B.meffil)[idx], sex=="M"])
}))
quantile(male.cg.r, probs=c(0.05,0.1,0.25,0.5))
```

```
##        5%       10%       25%       50% 
## 0.9830653 0.9997521 0.9998964 0.9999567
```

```r
female.idx <- which(rowSums(is.diff[,sex=="F"]) > 0)
female.cg.r <- sapply(female.idx[1:200], function(idx) {
    cor(B.meffil[idx, sex=="F"], B.minfi[rownames(B.meffil)[idx], sex=="F"])
})
quantile(female.cg.r, probs=c(0.05,0.1,0.25, 0.5))
```

```
##        5%       10%       25%       50% 
## 0.8622201 0.8948659 0.9509070 0.9993283
```


This file was generated from R markdown
using `knitr::knit("functional-normalization.rmd")`.
