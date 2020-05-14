library(devtools)
library(meffil)

library(openxlsx) ## load.450k.manifest() uses read.xlsx()
library(FlowSorted.CordBlood.450k) ## create.andrews.bakulski.reference() 

options(mc.cores=5)

## add 450k annotation
source("load-450k-manifest.r")
manifest.450k <- load.450k.manifest()
meffil.add.chip("450k", manifest.450k)

## add epic annotation
source("load-epic-manifest.r")
manifest.epic <- load.epic.manifest()
meffil.add.chip("epic", manifest.epic)

## create featureset common to both epic and 450k microarrays
featureset.450k <- meffil.featureset("450k")
featureset.epic <- meffil.featureset("epic")
featureset.both <- featureset.450k[which(with(featureset.450k, paste(type,target,name))
                                         %in%
                                         with(featureset.epic, paste(type,target,name))),]

## add this common featureset
meffil.add.featureset("common", featureset.both)

## create blood cell type references 
source("gse35069-references.r")
create.gse35069.references() ## ~15 minutes

source("gse68456-reference.r")
create.gse68456.reference()

source("gervin-lyle-reference.r")
create.gervin.lyle.reference("~/work/data/gervin-lyle-cord-blood-reference")

source("andrews-bakulski-reference.r")
create.andrews.bakulski.reference()

source("saliva-reference.r")
create.saliva.reference()

source("dlpfc-reference.r")
create.dlpfc.reference()

## save the global variables so they can be loaded by the package
meffil:::save.globals("../inst") ## see ../R/globals.r



