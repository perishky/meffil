

# Using `gdsfmt` for large datasets

## Download example data set 


> Sen A, Cingolani P, Senut MC, Land S et al. Lead exposure induces
> changes in 5-hydroxymethylcytosine clusters in CpG islands in human
> embryonic stem cells and umbilical cord blood. Epigenetics
> 2015;10(7):607-21. PMID: 26046694

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69633).





```r
download.450k.lead.dataset <- function() {
    dir.create(path <- "data-450k-lead")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        filename <-  file.path(path, "gse69633.tar")
        download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69633&format=file", filename)
        cat(date(), "Extracting files from GEO archive.\n")
        system(paste("cd", path, ";", "tar xvf", basename(filename)))
        unlink(filename)
        cat(date(), "Unzipping IDAT files.\n")
        system(paste("cd", path, ";", "gunzip *.idat.gz"))
        
        library(GEOquery)
        geo <- getGEO("GSE69633", GSEMatrix=F)
        ## this output from geo is really messed up.
        geo <- geo@gsms[[1]]@header
        characteristics <- matrix(geo$characteristics_ch1, nrow=length(geo$title), byrow=T)
        colnames(characteristics) <- sub("([^:]+):.*", "\\1", characteristics[1,])
        rownames(characteristics) <- geo$geo_accession
        characteristics <- apply(characteristics, 2, function(x) sub("[^:]+: (.*)", "\\1", x))
        characteristics <- as.data.frame(characteristics, stringsAsFactors=F)
        for (name in c("socioeconomic score", "gestational age", "birth weight", "pbconc (ng/dl)"))
            characteristics[[name]] <- as.numeric(characteristics[[name]])
        write.csv(characteristics, file=file.path(path, "samples.csv"))
    }
    
    path
}
```


```r
path <- download.450k.lead.dataset()
```

## Normalize dataset 

Create samplesheet

```r
library(meffil)
```

```
## Loading required package: illuminaio
```

```
## Loading required package: MASS
```

```
## Loading required package: limma
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```
## Loading required package: sva
```

```
## Loading required package: mgcv
```

```
## Loading required package: nlme
```

```
## This is mgcv 1.8-24. For overview type 'help("mgcv-package")'.
```

```
## Loading required package: genefilter
```

```
## 
## Attaching package: 'genefilter'
```

```
## The following object is masked from 'package:MASS':
## 
##     area
```

```
## Loading required package: BiocParallel
```

```
## Loading required package: ggplot2
```

```
## Loading required package: plyr
```

```
## Loading required package: reshape2
```

```
## Loading required package: gridExtra
```

```
## 
## Attaching package: 'gridExtra'
```

```
## The following object is masked from 'package:Biobase':
## 
##     combine
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     combine
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following object is masked from 'package:plyr':
## 
##     count
```

```
## The following objects are masked from 'package:genefilter':
## 
##     rowSds, rowVars
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: multcomp
```

```
## Loading required package: mvtnorm
```

```
## Loading required package: survival
```

```
## Loading required package: TH.data
```

```
## 
## Attaching package: 'TH.data'
```

```
## The following object is masked from 'package:MASS':
## 
##     geyser
```

```
## Loading required package: lme4
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'lme4'
```

```
## The following object is masked from 'package:nlme':
## 
##     lmList
```

```
## Loading required package: fastICA
```

```
## Loading required package: DNAcopy
```

```
## Loading required package: quadprog
```

```
## Loading required package: statmod
```

```
## Loading required package: gdsfmt
```

```r
options(mc.cores=10)
samplesheet <- meffil.create.samplesheet(path)

samples <- read.csv(file.path(path, "samples.csv"), check.names=F, row.names=1)
samplesheet <- data.frame(samplesheet,
                          samples[match(samplesheet$Sample_Name, rownames(samples)),],
                          stringsAsFactors=F, check.names=F)

samplesheet <- samplesheet[which(samplesheet[["sample type"]] == "HM450K"),]
```

Parameters.

```r
qc.file <- "gds/qc-report.html"
author <- "Sen, et al."
study <- "Cord blood DNA methylation and lead exposure (GSE69633)"
norm.file <- "gds/normalization-report.html"
cell.type.reference <- "gervin and lyle cord blood"
```

Generate QC objects for each sample and QC report.

```r
qc.objects <- meffil.qc(samplesheet, cell.type.reference=cell.type.reference, verbose=T)
```

```
## [read.idat] Thu Nov 29 08:07:02 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Grn.idat 
## [read.idat] Thu Nov 29 08:07:03 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Red.idat 
## [extract.detection.pvalues] Thu Nov 29 08:07:06 2018  
## [extract.beadnum] Thu Nov 29 08:07:12 2018  
## [extract.snp.betas] Thu Nov 29 08:07:15 2018  
## [extract.controls] Thu Nov 29 08:07:17 2018  
## [background.correct] Thu Nov 29 08:07:18 2018 background correction for dye = R 
## [background.correct] Thu Nov 29 08:07:21 2018 background correction for dye = G 
## [dye.bias.correct] Thu Nov 29 08:07:24 2018  
## [FUN] Thu Nov 29 08:07:31 2018 predicting sex
```

```r
qc.summary <- meffil.qc.summary(qc.objects, verbose=T)
```

```
## [meffil.qc.summary] Thu Nov 29 08:10:35 2018 Sex summary TRUE 
## [meffil.qc.summary] Thu Nov 29 08:10:35 2018 Meth vs unmeth summary 
## [meffil.qc.summary] Thu Nov 29 08:10:35 2018 Control means summary 
## [meffil.qc.summary] Thu Nov 29 08:10:35 2018 Sample detection summary 
## [meffil.qc.summary] Thu Nov 29 08:10:37 2018 CpG detection summary 
## [meffil.qc.summary] Thu Nov 29 08:10:38 2018 Sample bead numbers summary 
## [meffil.qc.summary] Thu Nov 29 08:10:39 2018 CpG bead numbers summary 
## [meffil.qc.summary] Thu Nov 29 08:10:41 2018 Cell count summary 
## [meffil.qc.summary] Thu Nov 29 08:10:41 2018 Genotype concordance
```

```r
meffil.qc.report(qc.summary,
                 output.file=qc.file,
                 author=author,
                 study=study)
```

```
## [meffil.qc.report] Thu Nov 29 08:10:41 2018 Writing report as html file to gds/qc-report.html
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

Remove any low quality samples.

```r
if (nrow(qc.summary$bad.samples) > 0)
    qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

samplesheet <- samplesheet[match(names(qc.objects), rownames(samplesheet)),]
```

Check how many principal components to include.

```r
print(meffil.plot.pc.fit(qc.objects, n.cross=3)$plot)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

Ten seems about right.

```r
number.pcs <- 10
```

Normalize dataset and generate normalization report.

```r
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=number.pcs, verbose=T)
```

```
## [meffil.normalize.quantiles] Thu Nov 29 08:11:41 2018 selecting dye correction reference 
## [meffil.normalize.quantiles] Thu Nov 29 08:11:41 2018 creating control matrix 
## [meffil.normalize.quantiles] Thu Nov 29 08:11:41 2018 normalizing quantiles 
## [FUN] Thu Nov 29 08:11:41 2018 genomic.iG M 
## [FUN] Thu Nov 29 08:11:41 2018 genomic.iG U 
## [FUN] Thu Nov 29 08:11:42 2018 genomic.iR M 
## [FUN] Thu Nov 29 08:11:42 2018 genomic.iR U 
## [FUN] Thu Nov 29 08:11:42 2018 genomic.ii M 
## [FUN] Thu Nov 29 08:11:42 2018 genomic.ii U 
## [FUN] Thu Nov 29 08:11:42 2018 autosomal.iG M 
## [FUN] Thu Nov 29 08:11:43 2018 autosomal.iG U 
## [FUN] Thu Nov 29 08:11:43 2018 autosomal.iR M 
## [FUN] Thu Nov 29 08:11:43 2018 autosomal.iR U 
## [FUN] Thu Nov 29 08:11:43 2018 autosomal.ii M 
## [FUN] Thu Nov 29 08:11:43 2018 autosomal.ii U 
## [FUN] Thu Nov 29 08:11:43 2018 not.y.iG M 
## [FUN] Thu Nov 29 08:11:43 2018 not.y.iG U 
## [FUN] Thu Nov 29 08:11:43 2018 not.y.iR M 
## [FUN] Thu Nov 29 08:11:44 2018 not.y.iR U 
## [FUN] Thu Nov 29 08:11:44 2018 not.y.ii M 
## [FUN] Thu Nov 29 08:11:44 2018 not.y.ii U 
## [FUN] Thu Nov 29 08:11:44 2018 sex M 
## [FUN] Thu Nov 29 08:11:44 2018 sex U 
## [FUN] Thu Nov 29 08:11:44 2018 chrx M 
## [FUN] Thu Nov 29 08:11:44 2018 chrx U 
## [FUN] Thu Nov 29 08:11:44 2018 chry M 
## [FUN] Thu Nov 29 08:11:45 2018 chry U
```

```r
beta <- meffil.normalize.samples(norm.objects,
                                 just.beta=T, 
                                 cpglist.remove=qc.summary$bad.cpgs$name,
                                 verbose=T)
```

```
## [read.idat] Thu Nov 29 08:11:45 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Grn.idat 
## [read.idat] Thu Nov 29 08:11:46 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Red.idat 
## [background.correct] Thu Nov 29 08:11:46 2018 background correction for dye = R 
## [background.correct] Thu Nov 29 08:11:48 2018 background correction for dye = G 
## [dye.bias.correct] Thu Nov 29 08:11:50 2018  
## [meffil.normalize.sample] Thu Nov 29 08:11:53 2018 Normalizing methylated and unmethylated signals.
```

Normalize while saving to a GDS file.

```r
gds.filename <- "gds/beta.gds"
if (!file.exists(gds.filename)) 
    meffil.normalize.samples(norm.objects,
                             just.beta=T, 
                             cpglist.remove=qc.summary$bad.cpgs$name,
                             gds.filename=gds.filename,
                             verbose=T)
```

```
## [read.idat] Thu Nov 29 08:13:39 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Grn.idat 
## [read.idat] Thu Nov 29 08:13:40 2018 Reading data-450k-lead/GSM1704996_5815073001_R04C01_Red.idat 
## [background.correct] Thu Nov 29 08:13:41 2018 background correction for dye = R 
## [background.correct] Thu Nov 29 08:13:43 2018 background correction for dye = G 
## [dye.bias.correct] Thu Nov 29 08:13:45 2018  
## [meffil.normalize.sample] Thu Nov 29 08:13:47 2018 Normalizing methylated and unmethylated signals.
```

```
## [1] "gds/beta.gds"
```

Load the matrix in the GDS file for comparison.

```r
gds.file <- openfn.gds(gds.filename)
matrix.node <- index.gdsn(gds.file, "matrix")
beta.gds <- read.gdsn(matrix.node)
rownames(beta.gds) <- read.gdsn(index.gdsn(gds.file, "row.names"))
colnames(beta.gds) <- read.gdsn(index.gdsn(gds.file, "col.names"))
closefn.gds(gds.file)
```

It should be the same as the one generated the standard way.

```r
all(beta.gds == beta)
```

```
## [1] TRUE
```

```r
identical(colnames(beta.gds), colnames(beta))
```

```
## [1] TRUE
```

```r
identical(rownames(beta.gds), rownames(beta))
```

```
## [1] TRUE
```

We've implemented an approach for generating principal components
that is equivalent to the standard approach.

```r
pcs <- meffil.methylation.pcs(beta)
pcs.gds <- meffil.methylation.pcs(gds.filename)
```


```r
all(pcs.gds == pcs)
```

```
## [1] TRUE
```

EWAS is not yet ready but here is the basic idea.

```r
cl <- makeCluster(getOption("mc.cores",1))
clusterExport(cl, c("samplesheet"))
system.time(p.gds <- clusterApply.gdsn(cl=cl,
                           gds.fn=gds.filename,
                           node.name="matrix",
                           margin=1,
                           as.is="double",
                           FUN=function(x) {
                               fit <- lm(x ~ lead + gender, data=data.frame(lead=samplesheet$pbconc, gender=samplesheet$gender))
                               coef(summary(fit))["lead","Pr(>|t|)"]  
                           }))
```

```
##    user  system elapsed 
##   0.957   0.253 141.898
```

We test associations in the beta matrix.

```r
system.time(p <- unlist(mclapply(1:nrow(beta), function(i) {
    fit <- lm(beta[i,] ~ lead + gender, data=data.frame(lead=samplesheet$pbconc, gender=samplesheet$gender))
    coef(summary(fit))["lead","Pr(>|t|)"]
})))
```

```
##     user   system  elapsed 
## 1847.544   88.447  219.088
```

We obtain identical results.

```r
table(p < 1e-3, p.gds < 1e-3)
```

```
##        
##          FALSE   TRUE
##   FALSE 485103      0
##   TRUE       0    272
```
