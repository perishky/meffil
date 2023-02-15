

# Normalize an EPIC v2 dataset

## Download example data set 


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
path <- download.epic2.demo.dataset()
```

## Normalize dataset 

Create samplesheet

```r
library(meffil)
samplesheet <- meffil.read.samplesheet(base=path, pattern="SampleSheet")
```

```
## [read.450k.sheet] Found the following CSV files:
## [1] "data-epic2-demo/Demo_EPIC-8v2-0_A1_SampleSheet_16.csv"
```

Parameters.

```r
qc.file <- "epic2-demo/qc-report.html"
author <- "Illumina, et al."
study <- "EPIC demo dataset"
number.pcs <- 2
norm.file <- "epic2-demo/normalization-report.html"
```

Generate QC objects.

```r
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose=T)
## If you get an error relating to a function from
## the 'preprocessCore' R package, 
## then you may need to reinstall it as follows: 
## BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
```

QC report.

```r
qc.summary <- meffil.qc.summary(qc.objects, verbose=T)
```

```
## [meffil.qc.summary] Wed Feb 15 10:33:13 2023 Sex summary TRUE 
## [meffil.qc.summary] Wed Feb 15 10:33:13 2023 Meth vs unmeth summary 
## [meffil.qc.summary] Wed Feb 15 10:33:13 2023 Control means summary 
## [meffil.qc.summary] Wed Feb 15 10:33:13 2023 Sample detection summary 
## [meffil.qc.summary] Wed Feb 15 10:33:17 2023 CpG detection summary 
## [meffil.qc.summary] Wed Feb 15 10:33:17 2023 Sample bead numbers summary 
## [meffil.qc.summary] Wed Feb 15 10:33:19 2023 CpG bead numbers summary 
## [meffil.qc.summary] Wed Feb 15 10:33:20 2023 Cell count summary 
## [meffil.qc.summary] Wed Feb 15 10:33:20 2023 Genotype concordance
```

```r
meffil.qc.report(qc.summary,
                 output.file=qc.file,
                 author=author,
                 study=study)
```

```
## [meffil.qc.report] Wed Feb 15 10:33:21 2023 Writing report as html file to epic2-demo/qc-report.html
```

Normalization dataset.

```r
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=number.pcs, verbose=T)
```

```
## [meffil.normalize.quantiles] Wed Feb 15 10:33:31 2023 selecting dye correction reference 
## [meffil.normalize.quantiles] Wed Feb 15 10:33:31 2023 creating control matrix 
## [meffil.normalize.quantiles] Wed Feb 15 10:33:31 2023 normalizing quantiles 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.iG M 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.iG U 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.iR M 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.iR U 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.ii M 
## [FUN] Wed Feb 15 10:33:31 2023 genomic.ii U 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.iG M 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.iG U 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.iR M 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.iR U 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.ii M 
## [FUN] Wed Feb 15 10:33:31 2023 autosomal.ii U 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.iG M 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.iG U 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.iR M 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.iR U 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.ii M 
## [FUN] Wed Feb 15 10:33:31 2023 not.y.ii U 
## [FUN] Wed Feb 15 10:33:31 2023 sex M 
## [FUN] Wed Feb 15 10:33:31 2023 sex U 
## [FUN] Wed Feb 15 10:33:32 2023 chrx M 
## [FUN] Wed Feb 15 10:33:32 2023 chrx U 
## [FUN] Wed Feb 15 10:33:32 2023 chry M 
## [FUN] Wed Feb 15 10:33:32 2023 chry U
```

```r
norm.dataset <- meffil.normalize.samples(norm.objects,
                                         just.beta=F, 
                                         cpglist.remove=qc.summary$bad.cpgs$name,
                                         verbose=T)
```


```r
beta.meffil <- meffil.get.beta(norm.dataset$M, norm.dataset$U)
```

Compute principal components of the normalized methylation matrix.

```r
pcs <- meffil.methylation.pcs(beta.meffil, sites=meffil.get.autosomal.sites("epic"), verbose=T)
```

```
## [meffil.methylation.pcs] Wed Feb 15 10:35:06 2023 Calculating CpG variance 
## [meffil.methylation.pcs] Wed Feb 15 10:35:07 2023 Calculating beta PCs
```

Normalization report.

```r
parameters <- meffil.normalization.parameters(norm.objects)
parameters$batch.threshold <- 0.01
norm.summary <- meffil.normalization.summary(norm.objects=norm.objects,
                                             pcs=pcs,
                                             parameters=parameters, verbose=T)
```

```
## [meffil.plot.control.batch] Wed Feb 15 10:35:07 2023 Extracting batch variables 
## [meffil.plot.control.batch] Wed Feb 15 10:35:07 2023 Testing associations 
## [meffil.plot.probe.batch] Wed Feb 15 10:35:13 2023 Extracting batch variables 
## [meffil.plot.probe.batch] Wed Feb 15 10:35:13 2023 Testing associations
```

```r
meffil.normalization.report(norm.summary,
                            output.file=norm.file,
                            author=author,
                            study=study)
```

```
## [meffil.normalization.report] Wed Feb 15 10:35:19 2023 Writing report as html file to epic2-demo/normalization-report.html
```


## Horvath's clock and the EPIC v2 microarray


```r
#require(RCurl)
#clock <- read.csv(textConnection(getURL("https://labs.genetics.ucla.edu/horvath/dnamage/AdditionalFile3.csv")), comment.char="#", stringsAsFactors=F)
clock <- read.csv("AdditionalFile3.csv", comment.char="#", stringsAsFactors=F)
length(setdiff(clock$CpGmarker, rownames(beta.meffil)))
```

```
## [1] 14
```

```r
nrow(clock)
```

```
## [1] 354
```

```r
length(intersect(clock$CpGmarker, rownames(beta.meffil)))/nrow(clock)
```

```
## [1] 0.960452
```

20 of the 'clock' sites are missing.

## Hannum predictor CpG sites and the EPIC microarray

71 CpG sites were used in the Hannum et al. age predictor:

> Hannum G, et al.
> Genome-wide methylation profiles reveal quantitative views of human aging rates.
> Mol Cell. 2013 Jan 24;49(2):359-67.

The list can be obtained from here:
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3780611/bin/NIHMS418935-supplement-02.xlsx


```r
hannum.sites <- readLines("hannum.txt")
length(hannum.sites)
```

```
## [1] 71
```

```r
length(setdiff(hannum.sites, rownames(beta.meffil)))
```

```
## [1] 7
```

```r
length(intersect(hannum.sites, rownames(beta.meffil)))/length(hannum.sites)
```

```
## [1] 0.9014085
```

6 of the 71 sites are missing.
