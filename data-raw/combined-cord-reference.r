## Gervin, K., Salas, L.A., Bakulski, K.M. et al. Systematic evaluation and validation of reference 
## and library selection methods for deconvolution of cord blood DNA methylation data. 
## Clin Epigenet 11, 125 (2019). https://doi.org/10.1186/s13148-019-0717-y

## recreate the FlowSorted.CordBloodCombined.450k reference in meffil
create.combined.cord.reference <- function(verbose=T) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!requireNamespace("minfi", quietly = TRUE))
        BiocManager::install("minfi")
    if (!requireNamespace("FlowSorted.CordBloodCombined.450k", quietly = TRUE))
        BiocManager::install("FlowSorted.CordBloodCombined.450k")
    if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE))
        BiocManager::install("IlluminaHumanMethylation450kmanifest")
    if (!requireNamespace("IDOLOptimizedCpGsCordBlood", quietly = TRUE))
        BiocManager::install("IDOLOptimizedCpGsCordBlood")
    require("minfi")
    require("FlowSorted.CordBloodCombined.450k")
    require("IDOLOptimizedCpGsCordBlood")
    ds <- FlowSorted.CordBloodCombined.450k()

    ## code from minfi::preprocessFunnorm
    ## (could not just use that function because Meth/Unmeth
    ## signals not retained and we need them here to create the cell type
    ## reference)
    reference <- preprocessNoob(ds)
    reference <- mapToGenome(reference)
    extracted.data <- minfi:::.extractFromRGSet450k(ds)
    samplesheet <- colData(ds) 
    sex <- sign(samplesheet$Sex == "F") + 1
    reference <- minfi:::.normalizeFunnorm450k(object=reference, extractedData=extracted.data,
                                               sex=sex, nPCs=10, verbose=F)
    M <- getMeth(reference)
    U <- getUnmeth(reference)
    
    data(IDOLOptimizedCpGsCordBlood)

    cell.types <- c("CD8T", "CD4T", "NK",
                    "Bcell", "Mono", "Gran", "nRBC")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference(
        "combined cord blood",
        M[,selected], U[,selected],
        cell.types=samplesheet$CellType[selected],
        chip="450k",
        featureset="common",
        specific.sites=IDOLOptimizedCpGsCordBlood,
        description="Derived from FlowSorted.CordBloodCombined.450k",
        verbose=verbose)
}
