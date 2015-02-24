# Comparison of functional normalization implementations

## Load example data set 

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55491).


```r
dir.create(path <- "data")
```

```
## Warning in dir.create(path <- "data"): 'data' already exists
```

```r
if (length(list.files(path, "*.idat$")) == 0) {
  filename <-  file.path(path, "gse55491.tar")
  download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
  cat(date(), "Extracting files from GEO archive.\n")
  system(paste("cd", path, ";", "tar xvf", basename(filename)))
  unlink(filename)
  cat(date(), "Unzipping IDAT files.\n")
  system(paste("cd", path, ";", "gunzip *.idat.gz"))
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
raw.minfi <- read.450k.exp(base = path)
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
## [probe.characteristics] Tue Feb 24 16:29:51 2015 extracting I-Red 
## [probe.characteristics] Tue Feb 24 16:29:52 2015 extracting I-Green 
## [probe.characteristics] Tue Feb 24 16:29:52 2015 extracting II 
## [probe.characteristics] Tue Feb 24 16:29:52 2015 extracting Control 
## [meffil.probe.info] Tue Feb 24 16:29:52 2015 reorganizing type information 
## [probe.locations] Tue Feb 24 16:30:11 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```
## Error in meffil.probe.info(): cannot change value of locked binding for 'saved.probe.info'
```

```r
basenames <- meffil.basenames(path)
```

Collect controls, perform background and dye bias correction and then
compute quantiles for each sample, one at a time.

```r
norm.objects <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename) 
})
```

```
## Warning in mclapply(basenames, function(basename) {: all scheduled cores
## encountered errors in user code
```
 
Normalize quantiles from each sample using the control
matrix to identify batch effects.

```r
norm.objects <- meffil.normalize.objects(norm.objects, number.pcs=2)
```

```
## Error: all(sapply(objects, is.normalization.object)) is not TRUE
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
## [read.idat] Tue Feb 24 16:30:54 2015 Reading data/GSM1338100_6057825094_R01C01_Grn.idat 
## [read.idat] Tue Feb 24 16:30:54 2015 Reading data/GSM1338100_6057825094_R01C01_Red.idat 
## [probe.characteristics] Tue Feb 24 16:30:55 2015 extracting I-Red 
## [probe.characteristics] Tue Feb 24 16:30:55 2015 extracting I-Green 
## [probe.characteristics] Tue Feb 24 16:30:55 2015 extracting II 
## [probe.characteristics] Tue Feb 24 16:30:55 2015 extracting Control 
## [meffil.probe.info] Tue Feb 24 16:30:55 2015 reorganizing type information 
## [probe.locations] Tue Feb 24 16:31:15 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```
## Error in meffil.probe.info(): cannot change value of locked binding for 'saved.probe.info'
```

```r
all(raw.meffil$R == getRed(raw.minfi)[names(raw.meffil$R),1])
```

```
## Error in eval(expr, envir, enclos): object 'raw.meffil' not found
```

### Background correction

Apply `meffil` background correction.

```r
bc.meffil <- meffil.background.correct(raw.meffil)
```

```
## Error in match(x, table, nomatch = 0L): object 'raw.meffil' not found
```

```r
M.meffil <- meffil.rg.to.mu(bc.meffil)$M
```

```
## Error in match(x, table, nomatch = 0L): object 'bc.meffil' not found
```

Apply `minfi` background correction.

```r
bc.minfi <- preprocessNoob(raw.minfi, dyeCorr = FALSE)
M.minfi <- getMeth(bc.minfi)[names(M.meffil),1]
```

```
## Error in eval(expr, envir, enclos): object 'M.meffil' not found
```

Compare results, should be identical.

```r
quantile(M.meffil - M.minfi)
```

```
## Error in quantile(M.meffil - M.minfi): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'M.meffil' not found
```

### Dye bias correction
Apply `meffil` correction.

```r
dye.meffil <- meffil.dye.bias.correct(bc.meffil,
                                      intensity=norm.objects[[1]]$reference.intensity)
```

```
## Error in match(x, table, nomatch = 0L): object 'bc.meffil' not found
```

```r
M.meffil <- meffil.rg.to.mu(dye.meffil)$M
```

```
## Error in match(x, table, nomatch = 0L): object 'dye.meffil' not found
```

```r
U.meffil <- meffil.rg.to.mu(dye.meffil)$U
```

```
## Error in match(x, table, nomatch = 0L): object 'dye.meffil' not found
```

Apply `minfi` background and dye bias correction.

```r
dye.minfi <- preprocessNoob(raw.minfi, dyeCorr = TRUE)
M.minfi <- getMeth(dye.minfi)[names(M.meffil),1]
```

```
## Error in eval(expr, envir, enclos): object 'M.meffil' not found
```

```r
U.minfi <- getUnmeth(dye.minfi)[names(U.meffil),1]
```

```
## Error in eval(expr, envir, enclos): object 'U.meffil' not found
```

Compare results, should be identical.

```r
quantile(M.meffil - M.minfi)
```

```
## Error in quantile(M.meffil - M.minfi): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'M.meffil' not found
```

```r
quantile(U.meffil - U.minfi)
```

```
## Error in quantile(U.meffil - U.minfi): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'U.meffil' not found
```

```r
quantile(M.meffil/(M.meffil+U.meffil+100) - M.minfi/(M.minfi+U.minfi+100))
```

```
## Error in quantile(M.meffil/(M.meffil + U.meffil + 100) - M.minfi/(M.minfi + : error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'M.meffil' not found
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
```

```
## Error in object$controls: $ operator is invalid for atomic vectors
```

```r
control.matrix.meffil <- t(meffil.preprocess.control.matrix(control.matrix.meffil))
```

```
## Error in t(meffil.preprocess.control.matrix(control.matrix.meffil)): error in evaluating the argument 'x' in selecting a method for function 't': Error: could not find function "meffil.preprocess.control.matrix"
```

The control variable order (columns) may be different between `minfi` and `meffil`
so they are reordered.

```r
i <- apply(control.matrix.minfi,2,function(v) {
    which.min(apply(control.matrix.meffil,2,function(w) {
        sum(abs(v-w))
    }))
})
```

```
## Error in which.min(apply(control.matrix.meffil, 2, function(w) {: error in evaluating the argument 'x' in selecting a method for function 'which.min': Error in apply(control.matrix.meffil, 2, function(w) { (from <text>#2) : 
##   object 'control.matrix.meffil' not found
```

```r
control.matrix.meffil <- control.matrix.meffil[,i]
```

```
## Error in eval(expr, envir, enclos): object 'control.matrix.meffil' not found
```

Compare resulting matrices, should be identical.

```r
quantile(control.matrix.meffil - control.matrix.minfi)
```

```
## Error in quantile(control.matrix.meffil - control.matrix.minfi): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'control.matrix.meffil' not found
```

### Final beta values


```r
B.minfi <- getBeta(norm.minfi)
quantile(B.meffil - B.minfi[rownames(B.meffil),])
```

```
## Error in quantile(B.meffil - B.minfi[rownames(B.meffil), ]): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error in B.meffil - B.minfi[rownames(B.meffil), ] : 
##   non-numeric argument to binary operator
```

On what chromosomes are the biggest differences appearing?

```r
probe.chromosome <- probes$chr[match(rownames(B.meffil), probes$name)]
```

```
## Error in eval(expr, envir, enclos): object 'probes' not found
```

```r
sex <- sapply(norm.objects, function(object) object$sex)
```

```
## Error in object$sex: $ operator is invalid for atomic vectors
```

```r
is.diff <- abs(B.meffil - B.minfi[rownames(B.meffil),]) > 1e-4
```

```
## Error in B.meffil - B.minfi[rownames(B.meffil), ]: non-numeric argument to binary operator
```

```r
table(probe.chromosome[which(is.diff, arr.ind=T)[,"row"]])
```

```
## Error in eval(expr, envir, enclos): object 'probe.chromosome' not found
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="M")], arr.ind=T)[,"row"]])
```

```
## Error in eval(expr, envir, enclos): object 'probe.chromosome' not found
```

```r
table(probe.chromosome[which(is.diff[,which(sex=="F")], arr.ind=T)[,"row"]])
```

```
## Error in eval(expr, envir, enclos): object 'probe.chromosome' not found
```
These results not surprising
because chromosome Y not handled like `minfi` in either males or females, 
and chromosome X is handled just like `minfi` in females but not in males.


```r
autosomal.cgs <- unique(probes$name[which(probes$chr %in% paste("chr", 1:22, sep=""))])
```

```
## Error in unique(probes$name[which(probes$chr %in% paste("chr", 1:22, sep = ""))]): error in evaluating the argument 'x' in selecting a method for function 'unique': Error: object 'probes' not found
```

```r
quantile(B.meffil[autosomal.cgs,] - B.minfi[autosomal.cgs,])
```

```
## Error in quantile(B.meffil[autosomal.cgs, ] - B.minfi[autosomal.cgs, ]): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'autosomal.cgs' not found
```

In spite of the differences, CG correlations between methods are pretty close to 1.

```r
male.idx <- which(rowSums(is.diff[,sex=="M"]) >= 5)
```

```
## Error in which(rowSums(is.diff[, sex == "M"]) >= 5): error in evaluating the argument 'x' in selecting a method for function 'which': Error in is.data.frame(x) : object 'is.diff' not found
```

```r
male.cg.r <- unlist(mclapply(male.idx, function(idx) {
    cor(B.meffil[idx, sex=="M"], B.minfi[rownames(B.meffil)[idx], sex=="M"])
}))
```

```
## Error in unlist(mclapply(male.idx, function(idx) {: error in evaluating the argument 'x' in selecting a method for function 'unlist': Error in mclapply(male.idx, function(idx) { : object 'male.idx' not found
```

```r
quantile(male.cg.r, probs=c(0.05,0.1,0.25,0.5))
```

```
## Error in quantile(male.cg.r, probs = c(0.05, 0.1, 0.25, 0.5)): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'male.cg.r' not found
```

```r
female.idx <- which(rowSums(is.diff[,sex=="F"]) > 0)
```

```
## Error in which(rowSums(is.diff[, sex == "F"]) > 0): error in evaluating the argument 'x' in selecting a method for function 'which': Error in is.data.frame(x) : object 'is.diff' not found
```

```r
female.cg.r <- sapply(female.idx[1:200], function(idx) {
    cor(B.meffil[idx, sex=="F"], B.minfi[rownames(B.meffil)[idx], sex=="F"])
})
```

```
## Error in sapply(female.idx[1:200], function(idx) {: error in evaluating the argument 'X' in selecting a method for function 'sapply': Error: object 'female.idx' not found
```

```r
quantile(female.cg.r, probs=c(0.05,0.1,0.25, 0.5))
```

```
## Error in quantile(female.cg.r, probs = c(0.05, 0.1, 0.25, 0.5)): error in evaluating the argument 'x' in selecting a method for function 'quantile': Error: object 'female.cg.r' not found
```


This file was generated from R markdown
using `knitr::knit("functional-normalization.rmd")`.
