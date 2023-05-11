

# Normalize an 450K, EPIC and EPIC v2 dataset

## Download example data set 


> Prickett AR, Ishida M, Böhm S, Frost JM et al. Genome-wide methylation
> analysis in Silver-Russell syndrome patients. Hum Genet 2015
> Mar;134(3):317-32. PMID: 25563730

Retrieve the data from the Gene Expression Omnibus (GEO) website
(http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55491).





```r
download.450k.demo.dataset <- function() {
    dir.create(path <- "data-450k-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        filename <-  file.path(path, "gse55491.tar")
        download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE55491&format=file", filename)
        cat(date(), "Extracting files from GEO archive.\n")
        system(paste("cd", path, ";", "tar xvf", basename(filename)))
        unlink(filename)
        cat(date(), "Unzipping IDAT files.\n")
        system(paste("cd", path, ";", "gunzip *.idat.gz"))

        library(GEOquery)
        geo <- getGEO("GSE55491", GSEMatrix=F)
        geo <- lapply(geo@gsms, function(gsm) unlist(gsm@header))
        geo <- do.call(rbind, geo)
        geo <- as.data.frame(geo, stringAsFactors=F)
        geo$group <- geo$characteristics_ch13
        geo$sex <-   geo$characteristics_ch11
        write.csv(geo, file=file.path(path, "samples.csv"))
    }
    
    path
}
```

Information about the new EPIC microarray can be obtained here:

http://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html





```r
download.epic.demo.dataset <- function() {
    dir.create(path <- "data-epic-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        download.zip <- function(url, path) {
            filename <- file.path(path, "data.zip")
            download.file(url, filename)
            filenames <- unzip(filename, junkpaths=T, exdir=path)
            unlink(filename)
            invisible(filenames)
        }    
        
        ftp.url <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC"
                
        download.zip(file.path(ftp.url, "infinium-methylationepic-demo-dataset.zip"), path)
    }
    
    path
}
```



Information about the new EPIC v2 microarray can be obtained here:

https://www.illumina.com/products/by-type/microarray-kits/infinium-methylation-epic.html




```r
download.epic2.demo.dataset <- function() {
    dir.create(path <- "data-epic2-demo")
    
    if (length(list.files(path, "*.idat$")) == 0) {
        download.zip <- function(url, path) {
            filename <- file.path(path, "data.zip")
            download.file(url, filename)
            filenames <- unzip(filename, junkpaths=T, exdir=path)
            unlink(filename)
            invisible(filenames)
        }    

        url <- "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/DemoDataEPIC_v2.zip"
        download.zip(url, file.path(path, "data.zip"))
    }
    
    path
}
```




```r
path.450k <- download.450k.demo.dataset()
path.epic <- download.epic.demo.dataset()
path.epic2 <- download.epic2.demo.dataset()
```

## Normalize dataset 

Create samplesheet

```r
library(meffil)
samplesheet.450k <- meffil.create.samplesheet(path=path.450k)
samplesheet.epic <- meffil.read.samplesheet(base=path.epic, pattern="Demo_SampleSheet.csv")
```

```
## [read.450k.sheet] Found the following CSV files:
## [1] "data-epic-demo/Demo_SampleSheet.csv"
```

```r
samplesheet.epic2 <- meffil.read.samplesheet(base=path.epic2, pattern="SampleSheet")
```

```
## [read.450k.sheet] Found the following CSV files:
## [1] "data-epic2-demo/Demo_EPIC-8v2-0_A1_SampleSheet_16.csv"
```

```r
samplesheet.450k$chip <- "450k"
samplesheet.epic$chip <- "epic"
samplesheet.epic2$chip <- "epic2"
common.cols <- intersect(colnames(samplesheet.450k),colnames(samplesheet.epic))
common.cols <- intersect(common.cols, colnames(samplesheet.epic2))
samplesheet <- rbind(
    samplesheet.450k[1:10,common.cols],
    samplesheet.epic[,common.cols],
    samplesheet.epic2[,common.cols])
```

Parameters.

```r
qc.file <- "complex-demo/qc-report.html"
author <- "Illumina, et al."
study <- "450k/epic/epic2 demo dataset"
number.pcs <- 2
norm.file <- "complex-demo/normalization-report.html"
```

Generate QC objects.

```r
options(mc.core=20)
qc.objects <- meffil.qc(samplesheet, featureset="450k:epic:epic2", cell.type.reference="blood gse35069 complete", verbose=T)
```

QC report.

```r
qc.summary <- meffil.qc.summary(qc.objects, verbose=T)
```

```
## [meffil.qc.summary] Thu May 11 16:25:07 2023 Sex summary TRUE 
## [meffil.qc.summary] Thu May 11 16:25:07 2023 Meth vs unmeth summary 
## [meffil.qc.summary] Thu May 11 16:25:07 2023 Control means summary 
## [meffil.qc.summary] Thu May 11 16:25:07 2023 Sample detection summary 
## [meffil.qc.summary] Thu May 11 16:25:09 2023 CpG detection summary 
## [meffil.qc.summary] Thu May 11 16:25:10 2023 Sample bead numbers summary 
## [meffil.qc.summary] Thu May 11 16:25:11 2023 CpG bead numbers summary 
## [meffil.qc.summary] Thu May 11 16:25:11 2023 Cell count summary 
## [meffil.qc.summary] Thu May 11 16:25:11 2023 Genotype concordance
```

```r
meffil.qc.report(
    qc.summary,
    output.file=qc.file,
    author=author,
    study=study)
```

```
## [meffil.qc.report] Thu May 11 16:25:11 2023 Writing report as html file to complex-demo/qc-report.html
```

Normalization dataset.

```r
norm.objects <- meffil.normalize.quantiles(
    qc.objects,
    number.pcs=number.pcs,
    verbose=T)
```

```
## [meffil.normalize.quantiles] Thu May 11 16:25:21 2023 selecting dye correction reference 
## [meffil.normalize.quantiles] Thu May 11 16:25:21 2023 creating control matrix 
## [meffil.normalize.quantiles] Thu May 11 16:25:21 2023 normalizing quantiles 
## [FUN] Thu May 11 16:25:21 2023 genomic.iG M 
## [FUN] Thu May 11 16:25:21 2023 genomic.iG U 
## [FUN] Thu May 11 16:25:21 2023 genomic.iR M 
## [FUN] Thu May 11 16:25:21 2023 genomic.iR U 
## [FUN] Thu May 11 16:25:21 2023 genomic.ii M 
## [FUN] Thu May 11 16:25:21 2023 genomic.ii U 
## [FUN] Thu May 11 16:25:21 2023 autosomal.iG M 
## [FUN] Thu May 11 16:25:21 2023 autosomal.iG U 
## [FUN] Thu May 11 16:25:21 2023 autosomal.iR M 
## [FUN] Thu May 11 16:25:21 2023 autosomal.iR U 
## [FUN] Thu May 11 16:25:21 2023 autosomal.ii M 
## [FUN] Thu May 11 16:25:21 2023 autosomal.ii U 
## [FUN] Thu May 11 16:25:21 2023 not.y.iG M 
## [FUN] Thu May 11 16:25:22 2023 not.y.iG U 
## [FUN] Thu May 11 16:25:22 2023 not.y.iR M 
## [FUN] Thu May 11 16:25:22 2023 not.y.iR U 
## [FUN] Thu May 11 16:25:22 2023 not.y.ii M 
## [FUN] Thu May 11 16:25:22 2023 not.y.ii U 
## [FUN] Thu May 11 16:25:22 2023 sex M 
## [FUN] Thu May 11 16:25:22 2023 sex U 
## [FUN] Thu May 11 16:25:22 2023 chrx M 
## [FUN] Thu May 11 16:25:22 2023 chrx U 
## [FUN] Thu May 11 16:25:22 2023 chry M 
## [FUN] Thu May 11 16:25:22 2023 chry U
```

```r
beta.meffil <- meffil.normalize.samples(
    norm.objects,
    just.beta=T, 
    cpglist.remove=qc.summary$bad.cpgs$name,
    verbose=T)
```

Compute principal components of the normalized methylation matrix.

```r
pcs <- meffil.methylation.pcs(
    beta.meffil,
    sites=meffil.get.autosomal.sites("epic"),
    verbose=T)
```

```
## [meffil.methylation.pcs] Thu May 11 16:27:27 2023 Calculating CpG variance 
## [meffil.methylation.pcs] Thu May 11 16:27:29 2023 Calculating beta PCs
```

Normalization report.

```r
parameters <- meffil.normalization.parameters(norm.objects)
parameters$batch.threshold <- 0.01
norm.summary <- meffil.normalization.summary(
    norm.objects=norm.objects,
    pcs=pcs,
    parameters=parameters,
    verbose=T)
```

```
## [meffil.plot.control.batch] Thu May 11 16:27:29 2023 Extracting batch variables 
## [meffil.plot.control.batch] Thu May 11 16:27:29 2023 Testing associations 
## [meffil.plot.probe.batch] Thu May 11 16:27:33 2023 Extracting batch variables 
## [meffil.plot.probe.batch] Thu May 11 16:27:33 2023 Testing associations
```

```r
meffil.normalization.report(
    norm.summary,
    output.file=norm.file,
    author=author,
    study=study)
```

```
## [meffil.normalization.report] Thu May 11 16:27:36 2023 Writing report as html file to complex-demo/normalization-report.html
```

It looks like the normalized methylation retains
has a strong association with chip.

```r
chip <- samplesheet$chip
stats <- t(apply(pcs,2,function(pc) {
    fit <- coef(summary(lm(pc ~ chip)))
    fit[which.min(fit[-1,"Pr(>|t|)"])+1,]
}))
stats[stats[,"Pr(>|t|)"] < 0.05,]
```

```
##     Estimate Std. Error  t value     Pr(>|t|)
## PC1 90.42404   15.10665 5.985712 2.554382e-06
## PC2 31.68942   12.30709 2.574892 1.606936e-02
## PC3 58.88792   13.92749 4.228180 2.572938e-04
## PC5 27.00016    3.95589 6.825306 3.033529e-07
```
This isn't surprising given that the data for each
chip comes from a different study.
If there was reason to believe that the the chip
associations were entirely technical,
then we could remove some of this variation
by including 'chip' as a random effect
in the functional normalization.
