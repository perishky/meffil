

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
```

```
## [read.idat] Mon Apr 17 00:40:23 2023 Reading data-epic2-demo/206891110001_R01C01_Grn.idat 
## [read.idat] Mon Apr 17 00:40:23 2023 Reading data-epic2-demo/206891110001_R01C01_Red.idat 
## [extract.detection.pvalues] Mon Apr 17 00:40:27 2023  
## [extract.beadnum] Mon Apr 17 00:40:31 2023  
## [extract.snp.betas] Mon Apr 17 00:40:35 2023  
## [extract.controls] Mon Apr 17 00:40:35 2023  
## [background.correct] Mon Apr 17 00:40:37 2023 background correction for dye = R 
## [background.correct] Mon Apr 17 00:40:39 2023 background correction for dye = G 
## [dye.bias.correct] Mon Apr 17 00:40:41 2023  
## [FUN] Mon Apr 17 00:40:46 2023 predicting sex
```

```r
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
## [meffil.qc.summary] Mon Apr 17 00:49:13 2023 Sex summary TRUE 
## [meffil.qc.summary] Mon Apr 17 00:49:13 2023 Meth vs unmeth summary 
## [meffil.qc.summary] Mon Apr 17 00:49:13 2023 Control means summary 
## [meffil.qc.summary] Mon Apr 17 00:49:13 2023 Sample detection summary 
## [meffil.qc.summary] Mon Apr 17 00:49:15 2023 CpG detection summary 
## [meffil.qc.summary] Mon Apr 17 00:49:15 2023 Sample bead numbers summary 
## [meffil.qc.summary] Mon Apr 17 00:49:16 2023 CpG bead numbers summary 
## [meffil.qc.summary] Mon Apr 17 00:49:17 2023 Cell count summary 
## [meffil.qc.summary] Mon Apr 17 00:49:17 2023 Genotype concordance
```

```r
meffil.qc.report(
    qc.summary,
    output.file=qc.file,
    author=author,
    study=study)
```

```
## [meffil.qc.report] Mon Apr 17 00:49:18 2023 Writing report as html file to epic2-demo/qc-report.html
```

Normalization dataset.

```r
norm.objects <- meffil.normalize.quantiles(
    qc.objects,
    number.pcs=number.pcs,
    verbose=T)
```

```
## [meffil.normalize.quantiles] Mon Apr 17 00:41:33 2023 selecting dye correction reference 
## [meffil.normalize.quantiles] Mon Apr 17 00:41:33 2023 creating control matrix 
## [meffil.normalize.quantiles] Mon Apr 17 00:41:33 2023 normalizing quantiles 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.iG M 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.iG U 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.iR M 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.iR U 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.ii M 
## [FUN] Mon Apr 17 00:41:33 2023 genomic.ii U 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.iG M 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.iG U 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.iR M 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.iR U 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.ii M 
## [FUN] Mon Apr 17 00:41:33 2023 autosomal.ii U 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.iG M 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.iG U 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.iR M 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.iR U 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.ii M 
## [FUN] Mon Apr 17 00:41:33 2023 not.y.ii U 
## [FUN] Mon Apr 17 00:41:34 2023 sex M 
## [FUN] Mon Apr 17 00:41:34 2023 sex U 
## [FUN] Mon Apr 17 00:41:34 2023 chrx M 
## [FUN] Mon Apr 17 00:41:34 2023 chrx U 
## [FUN] Mon Apr 17 00:41:34 2023 chry M 
## [FUN] Mon Apr 17 00:41:34 2023 chry U
```

```r
beta.meffil <- meffil.normalize.samples(
    norm.objects,
    just.beta=T, 
    cpglist.remove=qc.summary$bad.cpgs$name,
    verbose=T)
```

```
## [read.idat] Mon Apr 17 00:41:34 2023 Reading data-epic2-demo/206891110001_R01C01_Grn.idat 
## [read.idat] Mon Apr 17 00:41:34 2023 Reading data-epic2-demo/206891110001_R01C01_Red.idat 
## [background.correct] Mon Apr 17 00:41:34 2023 background correction for dye = R 
## [background.correct] Mon Apr 17 00:41:37 2023 background correction for dye = G 
## [dye.bias.correct] Mon Apr 17 00:41:39 2023  
## [meffil.normalize.sample] Mon Apr 17 00:41:42 2023 Normalizing methylated and unmethylated signals.
```

Compute principal components of the normalized methylation matrix.

```r
pcs <- meffil.methylation.pcs(
    beta.meffil,
    sites=meffil.get.autosomal.sites("epic"),
    verbose=T)
```

```
## [meffil.methylation.pcs] Mon Apr 17 00:42:11 2023 Calculating CpG variance 
## [meffil.methylation.pcs] Mon Apr 17 00:42:12 2023 Calculating beta PCs
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
## [meffil.plot.control.batch] Mon Apr 17 00:49:28 2023 Extracting batch variables 
## [meffil.plot.control.batch] Mon Apr 17 00:49:28 2023 Testing associations 
## [meffil.plot.probe.batch] Mon Apr 17 00:49:33 2023 Extracting batch variables 
## [meffil.plot.probe.batch] Mon Apr 17 00:49:33 2023 Testing associations
```

```r
meffil.normalization.report(
    norm.summary,
    output.file=norm.file,
    author=author,
    study=study)
```

```
## [meffil.normalization.report] Mon Apr 17 00:49:39 2023 Writing report as html file to epic2-demo/normalization-report.html
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

Some of the 'clock' sites are missing.

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

Some but not all of the sites are missing.

## Duplicated probes

Several thousand probes appear in duplicate across the microarray.
By default, the `meffil.normalize.samples` function collapses
duplicated probes by taking their median.
It is possible to a different summarizing function by setting the
`dup.fun` argument.
It is also possible to keep the duplicate probes by
setting `dup.fun` to `NULL` as follows:


```r
beta.dup <- meffil.normalize.samples(
    norm.objects,
    just.beta=T, 
    cpglist.remove=qc.summary$bad.cpgs$name,
    dup.fun=NULL,
    verbose=T)
```

```
## [read.idat] Mon Apr 17 00:42:32 2023 Reading data-epic2-demo/206891110001_R01C01_Grn.idat 
## [read.idat] Mon Apr 17 00:42:32 2023 Reading data-epic2-demo/206891110001_R01C01_Red.idat 
## [background.correct] Mon Apr 17 00:42:32 2023 background correction for dye = R 
## [background.correct] Mon Apr 17 00:42:34 2023 background correction for dye = G 
## [dye.bias.correct] Mon Apr 17 00:42:36 2023  
## [meffil.normalize.sample] Mon Apr 17 00:42:40 2023 Normalizing methylated and unmethylated signals.
```

We can then collapse the duplicate probes as follows:


```r
beta.nodup <- meffil.collapse.dups(beta.dup)
```

This function also has a `dup.fun` argument allowing the user
to specify a different summarizing function.

For example, cg06373096 appears
10 times.


```r
quantile(beta.nodup["cg06373096",]-beta.meffil["cg06373096",],na.rm=T)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

```r
quantile(colMedians(beta.dup[grepl("cg06373096",rownames(beta.dup)),],na.rm=T)-beta.meffil["cg06373096",],na.rm=T)
```

```
##   0%  25%  50%  75% 100% 
##    0    0    0    0    0
```

Here are the correlations between these duplicates:


```r
cor(t(beta.dup[grep("cg06373096",rownames(beta.dup)),]))
```

```
##                  cg06373096 cg06373096_TC12 cg06373096_TC13 cg06373096_TC14
## cg06373096        1.0000000       0.8243699       0.7592742       0.7793531
## cg06373096_TC12   0.8243699       1.0000000       0.7897399       0.7397427
## cg06373096_TC13   0.7592742       0.7897399       1.0000000       0.7173957
## cg06373096_TC14   0.7793531       0.7397427       0.7173957       1.0000000
## cg06373096_TC15   0.8752993       0.7641724       0.8347072       0.8772968
## cg06373096_TC16   0.8825980       0.8154338       0.7713652       0.8267420
## cg06373096_TC17   0.7712993       0.7361161       0.6659681       0.6809560
## cg06373096_TC18   0.8289695       0.6955757       0.6189466       0.5692075
## cg06373096_TC19   0.7045560       0.7007004       0.8266329       0.7450534
## cg06373096_TC110  0.8118309       0.5930558       0.7121591       0.7932888
##                  cg06373096_TC15 cg06373096_TC16 cg06373096_TC17
## cg06373096             0.8752993       0.8825980       0.7712993
## cg06373096_TC12        0.7641724       0.8154338       0.7361161
## cg06373096_TC13        0.8347072       0.7713652       0.6659681
## cg06373096_TC14        0.8772968       0.8267420       0.6809560
## cg06373096_TC15        1.0000000       0.8176969       0.7212996
## cg06373096_TC16        0.8176969       1.0000000       0.6575714
## cg06373096_TC17        0.7212996       0.6575714       1.0000000
## cg06373096_TC18        0.6435677       0.8375359       0.5615424
## cg06373096_TC19        0.7625179       0.6663931       0.5967316
## cg06373096_TC110       0.8022731       0.7850622       0.6125044
##                  cg06373096_TC18 cg06373096_TC19 cg06373096_TC110
## cg06373096             0.8289695       0.7045560        0.8118309
## cg06373096_TC12        0.6955757       0.7007004        0.5930558
## cg06373096_TC13        0.6189466       0.8266329        0.7121591
## cg06373096_TC14        0.5692075       0.7450534        0.7932888
## cg06373096_TC15        0.6435677       0.7625179        0.8022731
## cg06373096_TC16        0.8375359       0.6663931        0.7850622
## cg06373096_TC17        0.5615424       0.5967316        0.6125044
## cg06373096_TC18        1.0000000       0.6419688        0.6092255
## cg06373096_TC19        0.6419688       1.0000000        0.6502503
## cg06373096_TC110       0.6092255       0.6502503        1.0000000
```

```r
colMeans(cor(t(beta.dup[grep("cg06373096",rownames(beta.dup)),])))
```

```
##       cg06373096  cg06373096_TC12  cg06373096_TC13  cg06373096_TC14 
##        0.8237550        0.7658907        0.7696189        0.7729036 
##  cg06373096_TC15  cg06373096_TC16  cg06373096_TC17  cg06373096_TC18 
##        0.8098831        0.8060399        0.7003989        0.7006540 
##  cg06373096_TC19 cg06373096_TC110 
##        0.7294804        0.7369650
```

How do matrices differ when we collapse
duplicates (`beta.meffil` and `beta.dup`)
or not (`beta.dup`)?


```r
identical(rownames(beta.nodup), rownames(beta.meffil))
```

```
## [1] TRUE
```

```r
all(rownames(beta.nodup) %in% sub("_.*", "", rownames(beta.dup)))
```

```
## [1] TRUE
```

```r
all(sub("_.*", "", rownames(beta.dup)) %in% rownames(beta.nodup))
```

```
## [1] TRUE
```


```r
dim(beta.meffil)
```

```
## [1] 917191     16
```

```r
dim(beta.nodup)
```

```
## [1] 917191     16
```

```r
dim(beta.dup)
```

```
## [1] 923547     16
```
