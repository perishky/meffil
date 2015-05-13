#' Defines cell type references "blood gse35069" and "blood gse35069 complete"
#' for estimating blood cell counts.
#' Both use methylation profiles from Reinius et al. 2012 (PMID 25424692)
#' for purified blood cell types.
#' The first is based on
#' six cell types: CD4T, CD8T, Mono, Bcell, NK, Gran.
#' The second is based on 
#' the same cell types but with Gran replaced by Neu and Eos.
create.gse35069.references <- function(number.pcs=5, temp.dir=NULL, verbose=F) {
    msg(verbose=verbose)

    if (is.null(temp.dir)) {
        dir.create(temp.dir <- tempfile(tmpdir="."))
        on.exit(unlink(temp.dir, recursive=TRUE))
    }
    
    samplesheet <- retrieve.gse35069(temp.dir, verbose)

    ds <- meffil.normalize.dataset(samplesheet,
                                   just.beta=F,
                                   qc.file="gse35069-qc-report.html",
                                   author="Reinius, et al.",
                                   study="Purified blood cell type methylation (GEO:GSE35069)",
                                   number.pcs=number.pcs,
                                   norm.file="gse35069-normalization-report.html",
                                   verbose=verbose)

    samplesheet <- samplesheet[match(colnames(ds$M), samplesheet$Sample_Name),]
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Neu","Eos")
    selected <- samplesheet$CellType %in% cell.types
    complete <- create.cell.type.reference(ds$M[,selected], ds$U[,selected],
                                           cell.types=samplesheet$CellType[selected],
                                           verbose=verbose)
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Gran")
    selected <- samplesheet$CellType %in% cell.types
    standard <- create.cell.type.reference(ds$M[,selected], ds$U[,selected],
                                           cell.types=samplesheet$CellType[selected],
                                           verbose=verbose)

    list("blood gse35069"=standard,
         "blood gse35069 complete"=complete)
}


## Purified blood cell type methylation profiles
## Reinius et al. 2012 (PMID 25424692)
## 
## GEO does not have the IDAT files
## (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069);
## however, the `FlowSorted.Blood.450k` R package
## contains the following code to download these
## files from another website as well as the sample
## information (see `getKereData.R` script).

retrieve.gse35069 <- function(temp.dir, verbose=F) {
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(temp.dir)

    msg("Downloading data ...", verbose=verbose)
    download.file("http://publications.scilifelab.se/f742bd8a08ad42c3921ccaaaf0e3997a/file/sample_sheet_IDAT.csv", destfile="samplesheet.csv", quiet=T)
    download.file("http://publications.scilifelab.se/f742bd8a08ad42c3921ccaaaf0e3997a/file/idat_files.zip", destfile="idat-files.zip",quiet=T)

    msg("Unzipping data ...", verbose=verbose)
    unzip("idat-files.zip")

    samples <- read.csv("samplesheet.csv", stringsAsFactors = TRUE)
    names(samples)[names(samples) == "Sample"] <- "Sample_Name"
    names(samples)[names(samples) == "Chip.ID"] <- "Slide"
    samples$Array <- paste0(samples$Chip.Row.Pos, samples$Chip.CO.position)
    samples$Chip.CO.position <- samples$Chip.Row.Pos <- NULL
    samples$Array <- gsub("O", "0", samples$Array)
    samples$Basename <- file.path(temp.dir, paste0(samples$Slide, "_", samples$Array))
    samples$CellType <- samples$Type
    samples$CellType <- sub("^CD14\\+ Monocytes", "Mono", samples$CellType)
    samples$CellType <- sub("^CD19\\+ B-cells", "Bcell", samples$CellType)
    samples$CellType <- sub("^CD4\\+ T-cells", "CD4T", samples$CellType)
    samples$CellType <- sub("^CD56\\+ NK-cells", "NK", samples$CellType)
    samples$CellType <- sub("^CD8\\+ T-cells", "CD8T", samples$CellType)
    samples$CellType <- sub("^Granulocytes", "Gran", samples$CellType)
    samples$CellType <- sub("^Whole blood", "WBC", samples$CellType)
    samples$CellType <- sub("^Neutrophils", "Neu", samples$CellType)
    samples$CellType <- sub("^Eosinophils", "Eos", samples$CellType)
    samples$Sex <- "M"

    check.samplesheet(samples)
    
    samples
}


