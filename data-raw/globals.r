library(devtools)
library(meffil)

options(mc.cores=12)

## add 450k annotation
source("load-450k-manifest.r")
manifest.450k <- load.450k.manifest()
meffil.add.manifest("450k", manifest.450k)

## add epic annotation
source("load-epic-manifest.r")
manifest.epic <- load.epic.manifest()
meffil.add.manifest("epic", manifest.epic)

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
reference.globals <- create.gse35069.references() ## ~15 minutes

## save the global variables so they can be loaded by the package
meffil:::save.globals("../inst") ## see ../R/globals.r



