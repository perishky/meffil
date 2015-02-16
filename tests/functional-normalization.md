# Comparison of functional normalization implementations

## Load example data set 

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55491).


```r
dir.create(data.dir <- "data")
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

```
## Mon Feb 16 14:26:30 2015 Extracting files from GEO archive.
## Mon Feb 16 14:27:34 2015 Unzipping IDAT files.
```

## Load and normalize using `minfi`
Load the data.


```r
library(minfi)
```

```
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
source("../R/functional-normalization.r")
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
```

```
## [probe.characteristics] Mon Feb 16 14:31:52 2015 extracting I-Red 
## [probe.characteristics] Mon Feb 16 14:31:52 2015 extracting I-Green 
## [probe.characteristics] Mon Feb 16 14:31:52 2015 extracting II 
## [probe.characteristics] Mon Feb 16 14:31:52 2015 extracting Control 
## [meffil.probe.info] Mon Feb 16 14:31:52 2015 reorganizing type information 
## [probe.locations] Mon Feb 16 14:32:13 2015 loading probe genomic location annotation IlluminaHumanMethylation450kanno.ilmn12.hg19
```

```r
basenames <- meffil.basenames(data.dir)
```

Extract the controls for each sample.

```r
control.matrix <- do.call(cbind, mclapply(basenames, function(basename) {
    meffil.extract.controls(basename, probes) 
}))
colnames(control.matrix) <- basenames
```

Perform background and dye bias correction and then
compute quantiles for each sample.

```r
norm.objects <- mclapply(basenames, function(basename) {
    meffil.compute.normalization.object(basename, control.matrix, probes=probes) 
})
```

Normalize quantiles from each sample using the control
matrix to identify batch effects.

```r
norm.objects <- meffil.normalize.objects(norm.objects, control.matrix, number.pcs=2, probes=probes)
```

```
## [meffil.normalize.objects] Mon Feb 16 14:36:12 2015 preprocessing the control matrix 
## [meffil.normalize.objects] Mon Feb 16 14:36:12 2015 predicting sex 
## [meffil.normalize.objects] Mon Feb 16 14:36:12 2015 normalizing quantiles
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
raw.meffil <- meffil.read.rg(basenames[1])
```

```
## [read.idat] Mon Feb 16 14:38:47 2015 Reading data/GSM1338100_6057825094_R01C01_Grn.idat 
## [read.idat] Mon Feb 16 14:38:47 2015 Reading data/GSM1338100_6057825094_R01C01_Red.idat
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
                                      norm.objects[[1]]$dye.bias.factors$R,
                                      norm.objects[[1]]$dye.bias.factors$G)
M.meffil <- meffil.rg.to.mu(dye.meffil)$M
```

Apply `minfi` background and dye bias correction.

```r
dye.minfi <- preprocessNoob(raw.minfi, dyeCorr = TRUE)
M.minfi <- getMeth(dye.minfi)[names(M.meffil),1]
```

Compare results, should be identical.

```r
quantile(M.meffil - M.minfi)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

### Control matrix extracted

Extract `minfi` control matrix.

```r
library(matrixStats)                                      
```

```
## matrixStats v0.13.1 (2015-01-21) successfully loaded. See ?matrixStats for help.
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
control.matrix.meffil <- t(meffil.preprocess.control.matrix(control.matrix))
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
Since sex chromosomes are handled differently than minfi,
final beta values should be identical on the autosomes only.

```r
autosomal.sites <- unique(probes$name[which(probes$chr.type == "autosomal")])
B.minfi <- getBeta(norm.minfi)[rownames(B.meffil),]
quantile(B.meffil[autosomal.sites,] - B.minfi[autosomal.sites,])
```

```
##            0%           25%           50%           75%          100% 
## -3.360023e-06 -1.537325e-08 -2.677564e-10  0.000000e+00  2.393849e-06
```

This file was generated from R markdown
using `knitr::knit("functional-normalization.rmd")`.
