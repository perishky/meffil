library(meffil)

source("dataset-epic-demo.r")
path <- download.epic.demo.dataset()

samplesheet <- meffil.create.samplesheet(path)

options(mc.cores=1)

## Note that we are generating CNV profiles
## for EPIC microarrays using baselines
## derived from 450K microarrays.
cnv.list <- meffil.calculate.cnv(samplesheet, verbose=T)

cnv.matrix <- meffil.cnv.matrix(cnv.list)

## > quantile(cnv.matrix)
##      0%     25%     50%     75%    100% 
## -3.9896 -0.0601 -0.0177  0.0323  0.9058 
## > cor(cnv.matrix)
##              200144450018 200144450019 200144450021
## 200144450018    1.0000000    0.9441116    0.9727863
## 200144450019    0.9441116    1.0000000    0.9393877
## 200144450021    0.9727863    0.9393877    1.0000000

