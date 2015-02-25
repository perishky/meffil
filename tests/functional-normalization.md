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
raw.minfi <- read.450k.exp(base = path)
```

Normalize using the minfi package.

```r
norm.minfi <- preprocessFunnorm(raw.minfi, nPCs=2, sex=NULL, bgCorr=TRUE, dyeCorr=TRUE, verbose=TRUE)
```

```
## [preprocessFunnorm] Background and dye bias correction with noob 
## Using sample number 4 as reference level...
## [preprocessFunnorm] Mapping to genome
## [preprocessFunnorm] Quantile extraction
## [preprocessFunnorm] Normalization
```

## Load and normalize using `meffil`

Load the code, probe annotation and sample filename information.

```r
library(meffil)
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
## [meffil.normalize.objects] Wed Feb 25 01:56:59 2015 selecting dye correction reference 
## [meffil.normalize.objects] Wed Feb 25 01:56:59 2015 predicting sex 
## [meffil.normalize.objects] Wed Feb 25 01:56:59 2015 creating control matrix 
## [meffil.normalize.objects] Wed Feb 25 01:56:59 2015 normalizing quantiles 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.iG M 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.iG U 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.iR M 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.iR U 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.ii M 
## [FUN] Wed Feb 25 01:56:59 2015 genomic.ii U 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.iG M 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.iG U 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.iR M 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.iR U 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.ii M 
## [FUN] Wed Feb 25 01:56:59 2015 autosomal.ii U 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.iG M 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.iG U 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.iR M 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.iR U 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.ii M 
## [FUN] Wed Feb 25 01:56:59 2015 not.y.ii U 
## [FUN] Wed Feb 25 01:56:59 2015 sex M 
## [FUN] Wed Feb 25 01:56:59 2015 sex U 
## [FUN] Wed Feb 25 01:56:59 2015 chry M 
## [FUN] Wed Feb 25 01:56:59 2015 chry U 
## [FUN] Wed Feb 25 01:56:59 2015 chrx M 
## [FUN] Wed Feb 25 01:56:59 2015 chrx U
```

Apply quantile normalization to methylation levels of each sample.

```r
B.meffil <- do.call(cbind, mclapply(norm.objects, function(object) {
    meffil.get.beta(meffil.normalize.sample(object)) 
})) 
```

## Compare normalizations

First make sure that the raw data is the same.

```r
raw.meffil <- meffil.read.rg(basenames[1])
```

```
## [read.idat] Wed Feb 25 01:59:13 2015 Reading data/GSM1338100_6057825094_R01C01_Grn.idat 
## [read.idat] Wed Feb 25 01:59:14 2015 Reading data/GSM1338100_6057825094_R01C01_Red.idat
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
