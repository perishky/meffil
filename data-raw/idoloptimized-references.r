## Salas LA, Koestler DC (2021). FlowSorted.Blood.EPIC: Illumina EPIC
## data on immunomagnetic sorted peripheral adult blood cells. R package
## version 1.12.1,

create.idoloptmized.references <- function(verbose=T) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!requireNamespace("minfi", quietly = TRUE))
        BiocManager::install("minfi")
    if (!requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE))
        BiocManager::install("FlowSorted.Blood.EPIC")
    if (!requireNamespace("ExperimentHub", quietly = TRUE))
        BiocManager::install("ExperimentHub")
    
    if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
        BiocManager::install("IlluminaHumanMethylationEPICmanifest")

    require("minfi")
    source("minfiEPIC.r")
    require("FlowSorted.Blood.EPIC")
    
    ## module add languages/r/4.1.0
    require("ExperimentHub")
    hub <- ExperimentHub()
    query(hub, "FlowSorted.Blood.EPIC")
    ds <- hub[["EH1136"]]

    ## code from minfi::preprocessFunnorm
    ## (could not just use that function because Meth/Unmeth
    ## signals not retained and we need them here to create the cell type
    ## reference)
    reference <- preprocessNoob(ds)
    reference <- mapToGenome(reference)

    extracted.data <- extractFromRGSetEPIC(ds)
    samplesheet <- colData(ds)
    sex <- getSex(reference, cutoff = -3)$predictedSex
    sex <- sign(sex=="F")+1
    reference <- normalizeFunnormEPIC(
        object=reference, extractedData=extracted.data,
        sex=sex, nPCs=10, verbose=F)
    M <- getMeth(reference)
    U <- getUnmeth(reference)

    cell.types <- c("CD8T", "CD4T", "NK",
                    "Bcell", "Mono", "Neu")
    stopifnot(all(cell.types %in% samplesheet$CellType))
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference(
        "blood idoloptimized epic",
        M[,selected], U[,selected],
        cell.types=samplesheet$CellType[selected],
        chip="epic",
        featureset="epic",
        specific.sites=IDOLOptimizedCpGs,
        description="Derived from FlowSorted.Blood.EPIC",
        verbose=verbose)
    meffil.add.cell.type.reference(
        "blood idoloptimized",
        M[,selected], U[,selected],
        cell.types=samplesheet$CellType[selected],
        chip="epic",
        featureset="common",
        specific.sites=IDOLOptimizedCpGs450klegacy,
        description="Derived from FlowSorted.Blood.450k",
        verbose=verbose)
}
