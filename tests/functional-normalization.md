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

## Load and normalize using `minfi`
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
raw.minfi <- read.450k.exp(base = data.dir)
```

Normalize using the minfi package.

```r
norm.minfi <- preprocessFunnorm(raw.minfi, nPCs=2, sex=NULL, bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE)
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

## Load and normalize using `meffil`
Load the code, probe annotation and sample filename information.

```r
library(meffil)
```

```
## Loading required package: illuminaio
## Loading required package: MASS
```

```r
probes <- meffil.probe.info()
```

```
## [probe.characteristics] Mon Feb 23 13:29:06 2015 extracting I-Red 
## [probe.characteristics] Mon Feb 23 13:29:06 2015 extracting I-Green 
## [probe.characteristics] Mon Feb 23 13:29:06 2015 extracting II 
## [probe.characteristics] Mon Feb 23 13:29:06 2015 extracting Control 
## [meffil.probe.info] Mon Feb 23 13:29:06 2015 reorganizing type information 
## [probe.locations] Mon Feb 23 13:29:32 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```r
basenames <- meffil.basenames(data.dir)
```

Collect controls, perform background and dye bias correction and then
compute quantiles for each sample, one at a time.

```r
norm.objects <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename, probes=probes) 
})
```

Normalize quantiles from each sample using the control
matrix to identify batch effects.

```r
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=2)
```

```
## [meffil.normalize.objects] Mon Feb 23 13:33:50 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 23 13:33:50 2015 selecting dye correction reference 
## [meffil.normalize.objects] Mon Feb 23 13:33:50 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 23 13:33:50 2015 normalizing quantiles 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.iG M 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.iG U 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.iR M 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.iR U 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.ii M 
## [FUN] Mon Feb 23 13:33:50 2015 genomic.ii U 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.iG M 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.iG U 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.iR M 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.iR U 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.ii M 
## [FUN] Mon Feb 23 13:33:50 2015 autosomal.ii U 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.iG M 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.iG U 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.iR M 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.iR U 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.ii M 
## [FUN] Mon Feb 23 13:33:50 2015 not.y.ii U 
## [FUN] Mon Feb 23 13:33:50 2015 sex M 
## [FUN] Mon Feb 23 13:33:50 2015 sex U 
## [FUN] Mon Feb 23 13:33:50 2015 chry M 
## [FUN] Mon Feb 23 13:33:50 2015 chry U 
## [FUN] Mon Feb 23 13:33:50 2015 chrx M 
## [FUN] Mon Feb 23 13:33:50 2015 chrx U
```

Apply quantile normalization to methylation levels of each sample.

```r
B.meffil <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object, probes)) 
})) 
```

## Compare normalizations

First make sure that the raw data is the same.

```r
raw.meffil <- meffil.read.rg(basenames[1], probes)
```

```
## [read.idat] Mon Feb 23 13:37:08 2015 Reading data/GSM1338100_6057825094_R01C01_Grn.idat 
## [read.idat] Mon Feb 23 13:37:09 2015 Reading data/GSM1338100_6057825094_R01C01_Red.idat
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
bc.meffil <- meffil.background.correct(raw.meffil, probes)
M.meffil <- meffil.rg.to.mu(bc.meffil, probes)$M
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
                                      intensity=norm.objects[[1]]$reference.intensity,
                                      probes=probes)
M.meffil <- meffil.rg.to.mu(dye.meffil, probes)$M
U.meffil <- meffil.rg.to.mu(dye.meffil, probes)$U
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
control.matrix.meffil <- sapply(norm.objects, function(object) object$controls)
control.matrix.meffil <- t(meffil.preprocess.control.matrix(control.matrix.meffil))
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
## -3.193965e-01 -1.581674e-08 -2.892948e-10  0.000000e+00  2.227712e-01
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
##  chrX  chrY 
## 42489  9839
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="M")], arr.ind=T)[,"row"]])
```

```
## 
##  chrX  chrY 
## 42488  5307
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="F")], arr.ind=T)[,"row"]])
```

```
## 
## chrX chrY 
##    1 4532
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
## 0.9587299 0.9989199 0.9999511 0.9999821
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
## 0.7515260 0.8373563 0.9188234 0.9677237
```


This file was generated from R markdown
using `knitr::knit("functional-normalization.rmd")`.
