library(meffil)
library(pkgload)

library(openxlsx) ## load.450k.manifest() uses read.xlsx()

options(mc.cores=5)

## add 450k annotation
if (!"450k" %in% meffil.list.chips()) {
    source("load-450k-manifest.r")
    manifest.450k <- load.450k.manifest()
    meffil.add.chip("450k", manifest.450k)
}

## add epic annotation
if (!"epic" %in% meffil.list.chips()) {
    source("load-epic-manifest.r")
    manifest.epic <- load.epic.manifest()
    meffil.add.chip("epic", manifest.epic)
}

if (!"epic2" %in% meffil.list.chips()) {
    source("load-epic2-manifest.r")
    manifest.epic2 <- load.epic2.manifest()
    meffil.add.chip("epic2",manifest.epic2)
}

if (!"common" %in% meffil.list.featuresets()) {
    ## create featureset common to both epic and 450k microarrays
    featureset.450k <- meffil.featureset("450k")
    featureset.epic <- meffil.featureset("epic")
    idx <- which(with(featureset.450k,
                      (paste(type,target,name))
                      %in%
                      with(featureset.epic, paste(type,target,name))))
    featureset.both <- featureset.450k[idx,]

    ## add this common featureset
    meffil.add.featureset("common", featureset.both)
}

if (!"mouse" %in% meffil.list.featuresets()) {
    source("load-mouse-manifest.r")
    manifest.mouse <- load.mouse.manifest()
    meffil.add.chip("mouse",manifest.mouse,intersections=F)
}

## create blood cell type references
if (!"blood gse35069" %in% meffil.list.cell.type.references()) {
    source("gse35069-references.r")
    create.gse35069.references() ## ~15 minutes
}

if (!"cord blood gse68456" %in% meffil.list.cell.type.references()) {
    source("gse68456-reference.r")
    create.gse68456.reference()
}

if (!"blood idoloptimized" %in% meffil.list.cell.type.references()) {
    source("idoloptimized-references.r")
    create.idoloptimized.references()
}

if (!"gervin and lyle cord blood" %in% meffil.list.cell.type.references()) {
    source("gervin-lyle-reference.r")
    create.gervin.lyle.reference("~/work/data/gervin-lyle-cord-blood-reference")
}

if (!"andrews and bakulski cord blood" %in% meffil.list.cell.type.references()) {
    library(FlowSorted.CordBlood.450k)
    source("andrews-bakulski-reference.r")
    create.andrews.bakulski.reference()
}

if (!"combined cord blood" %in% meffil.list.cell.type.references()) {
    source("combined-cord-reference.r")
    create.combined.cord.reference()
}

if (!"saliva gse48472" %in% meffil.list.cell.type.references()) {
    source("saliva-reference.r")
    create.saliva.reference()
}

if (!"guintivano dlpfc" %in% meffil.list.cell.type.references()) {
    source("dlpfc-reference.r")
    create.dlpfc.reference()
}

if (!"blood gse167998" %in% meffil.list.cell.type.references()) {
    source("gse167998-reference.r")
    create.gse167998.reference()
}

## save the global variables so they can be loaded by the package
meffil:::save.globals("../inst") ## see ../R/globals.r



