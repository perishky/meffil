create.dlpfc.reference <- function(verbose=T) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    if (!requireNamespace("minfi", quietly = TRUE))
        BiocManager::install("minfi")
    if (!requireNamespace("FlowSorted.DLPFC.450k", quietly = TRUE))
        BiocManager::install("FlowSorted.DLPFC.450k")
    if (!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE))
        BiocManager::install("IlluminaHumanMethylation450kmanifest")
    require("minfi")
    require("FlowSorted.DLPFC.450k")
    data("FlowSorted.DLPFC.450k")

    ## code from minfi::preprocessFunnorm
    ## (could not just use that function because Meth/Unmeth
    ## signals not retained and we need them here to create the cell type
    ## reference)
    reference <- preprocessNoob(FlowSorted.DLPFC.450k)
    reference <- mapToGenome(reference)
    extracted.data <- minfi:::.extractFromRGSet450k(FlowSorted.DLPFC.450k)
    samplesheet <- colData(FlowSorted.DLPFC.450k) ## pData(phenoData(FlowSorted.DLPFC.450k))
    sex <- sign(samplesheet$sex == "Female") + 1
    reference <- minfi:::.normalizeFunnorm450k(object=reference, extractedData=extracted.data,
                                               sex=sex, nPCs=10, verbose=F)
    M <- getMeth(reference)
    U <- getUnmeth(reference)
    cell.types <- c("NeuN_neg","NeuN_pos")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference("guintivano dlpfc",
                                   M[,selected], U[,selected],
                                   cell.types=samplesheet$CellType[selected],
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)
}
