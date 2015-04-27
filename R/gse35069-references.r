#' Purified blood cell type methylation references
#'
#' Methylation data from Reinius et al. 2012 (PMID 25424692)
#' estimating cell blood cell counts in methylation profiles.
#'
#' @return A list of two objects,
#' each created by \code{meffil.create.cell.type.reference()}.
#' The first, named \code{standard}, is a reference based on
#' six cell types: CD4T, CD8T, Mono, Bcell, NK, Gran.
#' The second, named \code{complete}, is a reference based on 
#' the same cell types but with Gran replaced by Neu and Eos.
#' Each of these objects can be passed as the \code{reference}
#' parameter to the function \code{meffil.estimate.cell.counts()}.
#'
#' @export
meffil.gse35069.references <- function() {
    gse35069.references ## precomputed, see code in ../data-raw/
}


create.gse35069.references <- function(temp.dir=NULL, probes=meffil.probe.info(), verbose=F) {
    msg(verbose=verbose)

    if (is.null(temp.dir)) {
        dir.create(temp.dir <- tempfile(tmpdir="."))
        on.exit(unlink(temp.dir, recursive=TRUE))
    }
    gse35069 <- retrieve.gse35069(temp.dir, verbose) 
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Neu","Eos")
    selected <- gse35069$CellType %in% cell.types
    reference1 <- with(gse35069[selected,],
                       meffil.create.cell.type.reference(Basename, CellType,
                                                         temp.dir=temp.dir,
                                                         probes=probes,
                                                         verbose=verbose))
    
    cell.types <- c("CD4T","CD8T","Mono","Bcell","NK","Gran")
    selected <- gse35069$CellType %in% cell.types
    reference2 <- with(gse35069[selected,],
                       meffil.create.cell.type.reference(Basename, CellType,
                                                         temp.dir=temp.dir,
                                                         probes=probes,
                                                         verbose=verbose))

    list(complete=reference1, standard=reference2)
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
    samples$SampleID <- sub(".* ", "", samples$Sample_Name)
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

    samples
}


