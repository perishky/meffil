create.andrews.bakulski.reference <- function(verbose=T) {
    source("http://bioconductor.org/biocLite.R")
    if (!require("FlowSorted.CordBlood.450k")) 
        biocLite("FlowSorted.CordBlood.450k")
    if (!require("minfi"))
        biocLite("minfi")
    require("minfi")
    require("FlowSorted.CordBlood.450k")
    data("FlowSorted.CordBlood.450k")

    ## code from minfi::preprocessFunnorm
    ## (could not just use that function because Meth/Unmeth
    ## signals not retained and we need them here to create the cell type
    ## reference)
    reference <- preprocessNoob(FlowSorted.CordBlood.450k)
    reference <- mapToGenome(reference)
    extracted.data <- minfi:::.extractFromRGSet450k(FlowSorted.CordBlood.450k)
    samplesheet <- pData(phenoData(FlowSorted.CordBlood.450k))
    sex <- sign(samplesheet$Sex == "F") + 1
    reference <- minfi:::.normalizeFunnorm450k(object=reference, extractedData=extracted.data,
                                               sex=sex, nPCs=10, verbose=F)
    M <- getMeth(reference)
    U <- getUnmeth(reference)
    cell.types <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
    selected <- samplesheet$CellType %in% cell.types
    meffil.add.cell.type.reference("andrews and bakulski cord blood",
                                   M[,selected], U[,selected],
                                   cell.types=samplesheet$CellType[selected],
                                   chip="450k",
                                   featureset="common",
                                   verbose=verbose)
}
